"""SQLite 连接管理与事务。

并发模型：
  * 开启 WAL，允许多读单写；写事务一律 BEGIN IMMEDIATE 立即取写锁，
    使「样本上传」与「参考更新」严格串行化，网页读到的永远是
    某一版参考下完整、一致的汇总。
  * Python sqlite3 默认会在 DML 前自动 BEGIN，这里显式关闭
    （isolation_level=None），事务边界完全由本模块控制。
"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager

from flask import Flask, current_app, g

from .equivalence import equivalent_span
from .normalization import variant_kind


def get_db() -> sqlite3.Connection:
    if "db" not in g:
        conn = sqlite3.connect(
            current_app.config["DATABASE"],
            detect_types=sqlite3.PARSE_DECLTYPES,
        )
        conn.row_factory = sqlite3.Row
        conn.isolation_level = None  # 手动事务
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("PRAGMA busy_timeout=5000")
        g.db = conn
    return g.db


def close_db(_exc=None) -> None:
    conn: sqlite3.Connection | None = g.pop("db", None)
    if conn is not None:
        conn.close()


@contextmanager
def read_tx():
    """读事务：单条快照，保证一次请求内所有读取来自同一版本。"""
    conn = get_db()
    conn.execute("BEGIN")
    try:
        yield conn
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


@contextmanager
def write_tx():
    """写事务：立即取写锁，串行化所有上传/更新操作；失败整体回滚。"""
    conn = get_db()
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield conn
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


SCHEMA = """
CREATE TABLE IF NOT EXISTS reference_sequences (
    id            INTEGER PRIMARY KEY,
    sequence      TEXT NOT NULL,
    length        INTEGER NOT NULL,
    uploaded_at   TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS samples (
    id                INTEGER PRIMARY KEY,
    name              TEXT NOT NULL UNIQUE,
    reference_id      INTEGER REFERENCES reference_sequences(id),
    status            TEXT NOT NULL CHECK (status IN ('VALID', 'NEEDS_REUPLOAD')),
    error             TEXT,
    raw_tsv           TEXT NOT NULL,
    uploaded_at       TEXT NOT NULL DEFAULT (datetime('now')),
    revalidated_at    TEXT
);

CREATE TABLE IF NOT EXISTS sample_variants (
    id            INTEGER PRIMARY KEY,
    sample_id     INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
    pos           INTEGER NOT NULL,
    ref           TEXT NOT NULL,
    alt           TEXT NOT NULL,
    gt            TEXT NOT NULL,
    kind          TEXT NOT NULL,
    orig_pos      INTEGER NOT NULL,
    orig_ref      TEXT NOT NULL,
    orig_alt      TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_variant_canonical
    ON sample_variants(pos, ref, alt);
CREATE INDEX IF NOT EXISTS idx_variant_sample ON sample_variants(sample_id);

-- 每个 (参考, 规范变异) 的等价写法范围：闭区间 [lo, hi]，按变异类型口径取。
-- lo/hi 是派生数据，与检出记录同事务写入；旧库由 migrate() 补建。
CREATE TABLE IF NOT EXISTS variant_ranges (
    id            INTEGER PRIMARY KEY,
    reference_id  INTEGER NOT NULL REFERENCES reference_sequences(id) ON DELETE CASCADE,
    pos           INTEGER NOT NULL,
    ref           TEXT NOT NULL,
    alt           TEXT NOT NULL,
    kind          TEXT NOT NULL CHECK (kind IN ('SNP', 'INS', 'DEL')),
    lo            INTEGER NOT NULL,
    hi            INTEGER NOT NULL,
    UNIQUE (reference_id, pos, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_range_lookup ON variant_ranges(reference_id, lo, hi);
"""


def migration_tasks(conn: sqlite3.Connection) -> list[str]:
    """列出旧库里需要补齐的派生数据任务（只读检查，不修改任何数据）。"""
    latest = _latest_ref_id(conn)
    tasks: list[str] = []
    for ref in conn.execute(
        "SELECT id FROM reference_sequences ORDER BY id"
    ).fetchall():
        if _missing_range_count(conn, ref["id"], latest) > 0:
            tasks.append(f"build:{ref['id']}")
    if latest is not None:
        orphans = conn.execute(
            "SELECT 1 FROM variant_ranges r WHERE r.reference_id=? AND NOT EXISTS ("
            "  SELECT 1 FROM sample_variants v JOIN samples s ON s.id=v.sample_id "
            "  WHERE s.reference_id=? AND s.status='VALID' "
            "  AND v.pos=r.pos AND v.ref=r.ref AND v.alt=r.alt) LIMIT 1",
            (latest, latest),
        ).fetchone()
        if orphans is not None:
            # 当前参考上失去最后一份有效检出的范围不再参与扩展
            tasks.append("prune")
    return tasks


def _missing_range_count(
    conn: sqlite3.Connection, ref_id: int, latest: int | None
) -> int:
    if ref_id == latest:
        # 当前参考：有效检出必须有范围
        sql = (
            "SELECT COUNT(*) c FROM (SELECT DISTINCT v.pos, v.ref, v.alt "
            "FROM sample_variants v JOIN samples s ON s.id=v.sample_id "
            "WHERE s.reference_id=? AND s.status='VALID' "
            "AND NOT EXISTS (SELECT 1 FROM variant_ranges r WHERE r.reference_id=? "
            "  AND r.pos=v.pos AND r.ref=v.ref AND r.alt=v.alt))"
        )
    else:
        # 历史参考：其下所有检出都要有范围，换回该参考时可直接复用
        sql = (
            "SELECT COUNT(*) c FROM (SELECT DISTINCT v.pos, v.ref, v.alt "
            "FROM sample_variants v JOIN samples s ON s.id=v.sample_id "
            "WHERE s.reference_id=? "
            "AND NOT EXISTS (SELECT 1 FROM variant_ranges r WHERE r.reference_id=? "
            "  AND r.pos=v.pos AND r.ref=v.ref AND r.alt=v.alt))"
        )
    return conn.execute(sql, (ref_id, ref_id)).fetchone()["c"]


def build_reference_ranges(conn: sqlite3.Connection, ref_id: int, sequence: str) -> int:
    """为某个参考补齐全部检出对应的等价范围，返回新增行数。"""
    if ref_id == _latest_ref_id(conn):
        rows = conn.execute(
            "SELECT DISTINCT v.pos, v.ref, v.alt FROM sample_variants v "
            "JOIN samples s ON s.id=v.sample_id "
            "WHERE s.reference_id=? AND s.status='VALID'",
            (ref_id,),
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT DISTINCT v.pos, v.ref, v.alt FROM sample_variants v "
            "JOIN samples s ON s.id=v.sample_id WHERE s.reference_id=?",
            (ref_id,),
        ).fetchall()
    inserted = 0
    for row in rows:
        inserted += ensure_range(
            conn, ref_id, row["pos"], row["ref"], row["alt"], sequence
        )
    return inserted


def build_all_ranges(conn: sqlite3.Connection) -> int:
    """为全部参考补齐等价范围（迁移/测试用）。"""
    inserted = 0
    for ref in conn.execute("SELECT id, sequence FROM reference_sequences").fetchall():
        inserted += build_reference_ranges(conn, ref["id"], ref["sequence"])
    return inserted


def prune_reference_ranges(conn: sqlite3.Connection, ref_id: int) -> int:
    """删除当前参考上不再有任何有效检出的等价范围，返回删除行数。"""
    cur = conn.execute(
        "DELETE FROM variant_ranges WHERE reference_id=? AND NOT EXISTS ("
        "  SELECT 1 FROM sample_variants v JOIN samples s ON s.id=v.sample_id "
        "  WHERE s.reference_id=? AND s.status='VALID' "
        "  AND v.pos=variant_ranges.pos AND v.ref=variant_ranges.ref "
        "  AND v.alt=variant_ranges.alt)",
        (ref_id, ref_id),
    )
    return cur.rowcount


def ensure_range(
    conn: sqlite3.Connection,
    ref_id: int,
    pos: int,
    ref: str,
    alt: str,
    sequence: str,
) -> int:
    """确保 (参考, 规范变异) 的等价范围已入库（调用方负责在同一写事务内）。"""
    kind = variant_kind(ref, alt)
    if kind is None:
        return 0
    exists = conn.execute(
        "SELECT 1 FROM variant_ranges WHERE reference_id=? AND pos=? AND ref=? AND alt=?",
        (ref_id, pos, ref, alt),
    ).fetchone()
    if exists is not None:
        return 0
    lo, hi = equivalent_span(pos, ref, alt, sequence)
    cur = conn.execute(
        "INSERT OR IGNORE INTO variant_ranges "
        "(reference_id, pos, ref, alt, kind, lo, hi) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (ref_id, pos, ref, alt, kind, lo, hi),
    )
    return cur.rowcount


def ensure_ranges(
    conn: sqlite3.Connection, ref_id: int, sequence: str, variants
) -> int:
    """批量确保一份样本的等价范围入库；已存在的规范键不重复计算。"""
    existing = {
        (r["pos"], r["ref"], r["alt"])
        for r in conn.execute(
            "SELECT pos, ref, alt FROM variant_ranges WHERE reference_id=?",
            (ref_id,),
        ).fetchall()
    }
    inserted = 0
    for v in variants:
        key = (v["pos"], v["ref"], v["alt"])
        if key in existing:
            continue
        existing.add(key)
        inserted += ensure_range(
            conn, ref_id, v["pos"], v["ref"], v["alt"], sequence
        )
    return inserted


def _latest_ref_id(conn: sqlite3.Connection) -> int | None:
    row = conn.execute(
        "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
    ).fetchone()
    return None if row is None else row["id"]


def migrate(app: Flask) -> None:
    """启动时为旧库补建派生数据：等价范围索引。

    建表语句已在 init_db 中幂等执行；这里只补齐缺失的范围行并收敛当前参考。
    先只读检查，确实有缺口才取写锁；整个过程在一个写事务内完成，
    失败整体回滚，已有检出记录不受影响。
    """
    conn = sqlite3.connect(app.config["DATABASE"])
    conn.isolation_level = None
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    try:
        if not migration_tasks(conn):
            return
        conn.execute("BEGIN IMMEDIATE")
        try:
            # 取锁后重新判断，避免与启动瞬间的其它写入交错
            for task in migration_tasks(conn):
                if task.startswith("build:"):
                    ref_id = int(task.split(":", 1)[1])
                    ref = conn.execute(
                        "SELECT sequence FROM reference_sequences WHERE id=?",
                        (ref_id,),
                    ).fetchone()
                    if ref is not None:
                        build_reference_ranges(conn, ref_id, ref["sequence"])
                elif task == "prune":
                    latest = _latest_ref_id(conn)
                    if latest is not None:
                        prune_reference_ranges(conn, latest)
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise
    finally:
        conn.close()


def init_db(app: Flask) -> None:
    conn = sqlite3.connect(app.config["DATABASE"])
    try:
        conn.executescript(SCHEMA)
        conn.commit()
    finally:
        conn.close()


def init_app(app: Flask) -> None:
    app.teardown_appcontext(close_db)
    with app.app_context():
        init_db(app)
