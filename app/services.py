"""业务逻辑：参考/样本上传、重验、替换、汇总查询、区间核对。"""

from __future__ import annotations

import sqlite3

from .db import (
    build_reference_ranges,
    ensure_ranges,
    prune_reference_ranges,
    read_tx,
    write_tx,
)
from .normalization import (
    MAX_NAME_LEN,
    ValidationError,
    parse_reference,
    parse_tsv,
    validate_and_normalize,
)


# --------------------------------------------------------------------------- #
# 写操作
# --------------------------------------------------------------------------- #
def upload_reference(raw: bytes) -> dict:
    """上传新参考；若已存在样本，全部按新参考重新校验。

    成功后：
      * 通过重验的样本 -> VALID，变体替换为新参考下的归一化结果；
      * REF 不符 / 归一化后重复的样本 -> NEEDS_REUPLOAD，变体删除，记录原因，
        从汇总剔除；
      * 旧参考及旧变体在同一事务中整体替换，失败则全部回滚。
    """
    sequence = parse_reference(raw)
    with write_tx() as conn:
        cur = conn.execute("SELECT id FROM samples")
        old_samples = cur.fetchall()

        cur = conn.execute(
            "INSERT INTO reference_sequences (sequence, length) VALUES (?, ?)",
            (sequence, len(sequence)),
        )
        new_ref_id = cur.lastrowid

        revalidated, failed = [], []
        for row in old_samples:
            outcome = _revalidate_sample(conn, row["id"], new_ref_id, sequence)
            (revalidated if outcome else failed).append(
                conn.execute(
                    "SELECT name, status, error FROM samples WHERE id=?",
                    (row["id"],),
                ).fetchone()
            )
        # 新参考下的等价范围随重验结果一同建立（同一写事务）
        build_reference_ranges(conn, new_ref_id, sequence)
    return {
        "reference_id": new_ref_id,
        "length": len(sequence),
        "revalidated": revalidated,
        "failed": failed,
    }


def upload_sample(name: str, raw: bytes) -> dict:
    """上传（或同名重传）一份样本 TSV。

    重传时先完整校验新文件，成功才替换旧版本；旧变体在同一事务中删除、
    新变体写入，任何错误整体回滚，当前结果不受影响。
    """
    name = (name or "").strip()
    if not name:
        raise ValidationError("样本名称不能为空。")
    if len(name) > MAX_NAME_LEN:
        raise ValidationError(f"样本名称不能超过 {MAX_NAME_LEN} 个字符。")

    # 解析在事务外完成（不触碰数据库，不影响现有结果）
    rows = parse_tsv(raw)
    raw_text = raw.decode("utf-8")

    with write_tx() as conn:
        ref_row = conn.execute(
            "SELECT id, sequence FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if ref_row is None:
            raise ValidationError("请先上传参考序列，再上传样本。")
        ref_id, ref_seq = ref_row["id"], ref_row["sequence"]

        # 校验 + 左对齐 + 去冗余；失败抛 ValidationError，事务回滚
        variants = validate_and_normalize(rows, ref_seq)

        old = conn.execute("SELECT id FROM samples WHERE name=?", (name,)).fetchone()
        if old is not None:
            conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (old["id"],))
            conn.execute("DELETE FROM samples WHERE id=?", (old["id"],))

        cur = conn.execute(
            "INSERT INTO samples (name, reference_id, status, error, raw_tsv) "
            "VALUES (?, ?, 'VALID', NULL, ?)",
            (name, ref_id, raw_text),
        )
        sample_id = cur.lastrowid
        _insert_variants(conn, sample_id, variants)
        _index_variants(conn, ref_id, variants)
        # 重传可能删掉了旧版本独有的变异：失去最后一份有效检出的范围不再参与扩展
        prune_reference_ranges(conn, ref_id)
        replaced = old is not None

    return {"sample": name, "count": len(variants), "replaced": replaced}


def delete_sample(sample_id: int) -> None:
    with write_tx() as conn:
        row = conn.execute(
            "SELECT s.reference_id FROM samples s WHERE s.id=?", (sample_id,)
        ).fetchone()
        conn.execute("DELETE FROM samples WHERE id=?", (sample_id,))
        # 当前参考上已无有效检出的变异不再参与区间扩展
        if row is not None and row["reference_id"] is not None:
            latest = conn.execute(
                "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
            ).fetchone()
            if latest is not None and latest["id"] == row["reference_id"]:
                prune_reference_ranges(conn, latest["id"])

    return {"sample": name, "count": len(variants), "replaced": old is not None}


def delete_sample(sample_id: int) -> None:
    with write_tx() as conn:
        row = conn.execute(
            "SELECT s.reference_id FROM samples s WHERE s.id=?", (sample_id,)
        ).fetchone()
        conn.execute("DELETE FROM samples WHERE id=?", (sample_id,))
        # 当前参考上已无有效检出的变异不再参与区间扩展
        if row is not None and row["reference_id"] is not None:
            latest = conn.execute(
                "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
            ).fetchone()
            if latest is not None and latest["id"] == row["reference_id"]:
                prune_reference_ranges(conn, latest["id"])


def _revalidate_sample(
    conn: sqlite3.Connection, sample_id: int, ref_id: int, ref_seq: str
) -> bool:
    """按给定参考重新校验一个样本，返回是否通过。"""
    row = conn.execute(
        "SELECT raw_tsv FROM samples WHERE id=?", (sample_id,)
    ).fetchone()
    try:
        rows = parse_tsv(row["raw_tsv"])
        variants = validate_and_normalize(rows, ref_seq)
    except ValidationError as exc:
        conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
        conn.execute(
            "UPDATE samples SET reference_id=?, status='NEEDS_REUPLOAD', "
            "error=?, revalidated_at=datetime('now') WHERE id=?",
            (ref_id, str(exc), sample_id),
        )
        return False

    conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
    conn.execute(
        "UPDATE samples SET reference_id=?, status='VALID', error=NULL, "
        "revalidated_at=datetime('now') WHERE id=?",
        (ref_id, sample_id),
    )
    _insert_variants(conn, sample_id, variants)
    _index_variants(conn, ref_id, variants)
    return True


def _index_variants(
    conn: sqlite3.Connection, ref_id: int, variants: list[dict]
) -> None:
    """为一份样本的规范变异补写等价范围（与检出记录同一写事务）。"""
    row = conn.execute(
        "SELECT sequence FROM reference_sequences WHERE id=?", (ref_id,)
    ).fetchone()
    if row is None:
        return
    ensure_ranges(conn, ref_id, row["sequence"], variants)


def _insert_variants(
    conn: sqlite3.Connection, sample_id: int, variants: list[dict]
) -> None:
    conn.executemany(
        "INSERT INTO sample_variants "
        "(sample_id, pos, ref, alt, gt, kind, orig_pos, orig_ref, orig_alt) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                sample_id,
                v["pos"],
                v["ref"],
                v["alt"],
                v["gt"],
                v["kind"],
                v["orig_pos"],
                v["orig_ref"],
                v["orig_alt"],
            )
            for v in variants
        ],
    )


# --------------------------------------------------------------------------- #
# 读操作（同一读事务/快照内完成，页面只展示同一版参考的完整汇总）
# --------------------------------------------------------------------------- #
def get_dashboard() -> dict:
    with read_tx() as conn:
        ref = conn.execute(
            "SELECT * FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if ref is None:
            return {"reference": None, "stats": {}, "samples": [], "variants": []}

        ref_id = ref["id"]
        stats = {
            "valid": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id=? AND status='VALID'",
                (ref_id,),
            ).fetchone()["c"],
            "needs_reupload": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id=? "
                "AND status='NEEDS_REUPLOAD'",
                (ref_id,),
            ).fetchone()["c"],
            "stale": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id IS NOT NULL "
                "AND reference_id<>?",
                (ref_id,),
            ).fetchone()["c"],
        }
        samples = conn.execute(
            "SELECT s.*, (SELECT COUNT(*) FROM sample_variants v WHERE v.sample_id=s.id) "
            "AS variant_count FROM samples s ORDER BY s.name"
        ).fetchall()
        variants = conn.execute(
            "SELECT v.pos, v.ref, v.alt, v.kind, COUNT(DISTINCT v.sample_id) AS n_samples "
            "FROM sample_variants v "
            "JOIN samples s ON s.id = v.sample_id AND s.reference_id=? AND s.status='VALID' "
            "GROUP BY v.pos, v.ref, v.alt "
            "ORDER BY v.pos, v.ref, v.alt",
            (ref_id,),
        ).fetchall()
        return {
            "reference": ref,
            "stats": stats,
            "samples": samples,
            "variants": variants,
        }


def get_variant_detail(pos: int, ref: str, alt: str) -> dict | None:
    with read_tx() as conn:
        current = conn.execute(
            "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if current is None:
            return None
        rows = conn.execute(
            "SELECT s.name AS sample, v.gt, v.orig_pos, v.orig_ref, v.orig_alt, v.kind "
            "FROM sample_variants v "
            "JOIN samples s ON s.id = v.sample_id "
            "WHERE v.pos=? AND v.ref=? AND v.alt=? "
            "AND s.reference_id=? AND s.status='VALID' "
            "ORDER BY s.name",
            (pos, ref, alt, current["id"]),
        ).fetchall()
        if not rows:
            return None
        return {"pos": pos, "ref": ref, "alt": alt, "rows": rows}


# --------------------------------------------------------------------------- #
# 区间核对
# --------------------------------------------------------------------------- #
def get_interval_form() -> dict:
    """区间核对页面所需的基础数据：当前参考与表单默认值（只读快照）。"""
    with read_tx() as conn:
        ref = conn.execute(
            "SELECT id, sequence, length FROM reference_sequences "
            "ORDER BY id DESC LIMIT 1"
        ).fetchone()
    return {
        "reference": ref,
        "reference_id": ref["id"] if ref else None,
        "start": "",
        "end": "",
        "error": None,
        "checked": False,
        "result": None,
    }


def stale_reference_message(
    expected_reference_id: int | None, current_id: int | None
) -> str:
    """换版后用旧表单查询时的提示。"""
    return (
        f"查询表单对应参考 #{expected_reference_id}，当前参考已是 "
        f"#{current_id}，请重新选择参考后重试。"
    )


class StaleReferenceError(Exception):
    """查询表单携带的参考版本已不是当前参考。"""


def parse_interval(start: str, end: str, genome_len: int) -> tuple[int, int]:
    """解析从 1 起算的闭区间，返回 (start, end)。

    起点/终点必须是纯数字正整数；不做静默截断，非法输入直接报错。
    """
    def _to_int(value: str, name: str) -> int:
        text = str(value).strip()
        if not text.isdigit() or int(text) < 1:
            raise ValidationError(f"{name} 必须是从 1 起算的正整数。")
        return int(text)

    lo = _to_int(start, "区间起点")
    hi = _to_int(end, "区间终点")
    if lo > hi:
        raise ValidationError(f"区间起点（{lo}）不能大于终点（{hi}）。")
    if hi > genome_len:
        raise ValidationError(f"区间终点 {hi} 超出参考序列长度 {genome_len}。")
    return lo, hi


def check_interval(
    start: str, end: str, expected_reference_id: int | None = None
) -> dict:
    """核对区间：返回完整收录相交变异所需的最小扩展区间及新增收录的变异。

    * 只读同一快照中的等价范围索引与检出记录。
    * 扩展后新碰到的变异继续完整收录，直到区间边界稳定（迭代不动点）。
    * 表单携带的参考版本（``expected_reference_id``）不是当前参考时抛
      StaleReferenceError；未携带视为不指定版本。
    """
    with read_tx() as conn:
        ref = conn.execute(
            "SELECT id, sequence, length FROM reference_sequences "
            "ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if ref is None:
            raise ValidationError("请先上传参考序列，再使用区间核对。")
        if expected_reference_id is not None and expected_reference_id != ref["id"]:
            raise StaleReferenceError(
                stale_reference_message(expected_reference_id, ref["id"])
            )

        orig_lo, orig_hi = parse_interval(start, end, ref["length"])
        req_lo, req_hi = orig_lo, orig_hi

        # 扩展：收录所有与当前区间相交的变异；新收录的变异若越界则继续外扩，
        # 直到边界不再变化。范围行都来自当前参考的有效检出，不会跨版本混算。
        prev = None
        while prev != (req_lo, req_hi):
            prev = (req_lo, req_hi)
            row = conn.execute(
                "SELECT MIN(lo) AS lo, MAX(hi) AS hi FROM variant_ranges "
                "WHERE reference_id=? AND lo <= ? AND hi >= ?",
                (ref["id"], req_hi, req_lo),
            ).fetchone()
            if row["lo"] is not None:
                req_lo = min(req_lo, row["lo"])
                req_hi = max(req_hi, row["hi"])

        hits = conn.execute(
            "SELECT r.pos, r.ref, r.alt, r.kind, r.lo, r.hi, "
            "(SELECT COUNT(DISTINCT v.sample_id) FROM sample_variants v "
            " JOIN samples s ON s.id=v.sample_id "
            " WHERE s.reference_id=? AND s.status='VALID' "
            "   AND v.pos=r.pos AND v.ref=r.ref AND v.alt=r.alt) AS n_samples "
            "FROM variant_ranges r "
            "WHERE r.reference_id=? AND r.lo <= ? AND r.hi >= ? "
            "ORDER BY r.pos, r.ref, r.alt",
            (ref["id"], ref["id"], req_hi, req_lo),
        ).fetchall()

        return {
            "reference": ref,
            "requested": (orig_lo, orig_hi),
            "expanded": (req_lo, req_hi),
            "hits": hits,
            # 扩展后收录、但原区间未能完整覆盖的变异
            "uncovered": [
                h for h in hits if not (h["lo"] >= orig_lo and h["hi"] <= orig_hi)
            ],
        }
