"""区间核对与等价范围索引的测试。"""

import io
import os
import sqlite3
import tempfile
import unittest

from app import create_app, services
from app.equivalence import equivalent_span
from app.normalization import ValidationError, normalize_variant

REF_V1 = "GCTTTGCACGT"      # 11 bp：TTT 重复区
REF_V2 = "GCTATGCACGT"      # 第 4 位 T->A

# 位置:  1 2 3 4 5 6 7 8 9 10 11
#        G C T T T G C A C G  T


def tsv(text: str) -> bytes:
    return text.encode("utf-8")


class AppTestBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        os.environ["DATABASE_PATH"] = os.path.join(self.tmp.name, "test.db")
        self.app = create_app()
        self.app.config.update(TESTING=True, SECRET_KEY="test")
        self.client = self.app.test_client()

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    def upload_ref(self, seq: str):
        data = {"reference": (io.BytesIO(seq.encode()), "ref.txt")}
        return self.client.post(
            "/reference", data=data, content_type="multipart/form-data",
            follow_redirects=True,
        )

    def upload_sample(self, name: str, content: bytes):
        data = {"name": name, "sample": (io.BytesIO(content), f"{name}.tsv")}
        return self.client.post(
            "/samples", data=data, content_type="multipart/form-data",
            follow_redirects=True,
        )


class EquivalentSpanTest(unittest.TestCase):
    """等价范围口径：SNP 按位点，INS 按左锚，DEL 按实际删去的碱基。"""

    GENOME = "GCTTTGC"

    def test_snp_is_single_position(self):
        self.assertEqual(equivalent_span(3, "T", "A", self.GENOME), (3, 3))

    def test_insertion_uses_anchor_positions(self):
        # TTT 重复区插入一个 T，可锚在 2/3/4/5 位
        self.assertEqual(equivalent_span(2, "C", "CT", self.GENOME), (2, 5))

    def test_deletion_uses_deleted_bases_not_anchor(self):
        # (2,CT,C) 删掉位置 3 的 T；可等价删除位置 3/4/5
        self.assertEqual(equivalent_span(2, "CT", "C", self.GENOME), (3, 5))

    def test_insertion_of_multibase_repeat_rotates(self):
        # ATAT 重复区插入 "AT"：旋转后仍是同一个插入
        genome = "GATATATC"
        lo, hi = equivalent_span(1, "G", "GAT", genome)
        self.assertEqual((lo, hi), (1, hi))
        self.assertGreater(hi, 1)

    def test_span_covers_all_equivalent_writings(self):
        # 暴力枚举每个基因组上所有合法书写，按「单倍型 + 事件类型」分组，
        # 验证等价范围首尾正好覆盖该组所有写法的锚点/删除碱基。
        import itertools

        insertions = [
            "".join(t) for n in (1, 2) for t in itertools.product("ACGT", repeat=n)
        ]

        def writings(genome):
            out = []
            for p in range(1, len(genome) + 1):
                base = genome[p - 1]
                for b in "ACGT":
                    if b != base:
                        out.append((p, base, b, "SNP"))
                for ins in insertions:
                    out.append((p, base, base + ins, "INS"))
                for dlen in (1, 2):
                    if p + dlen <= len(genome):
                        out.append((p, genome[p - 1: p + dlen], base, "DEL"))
            return out

        def hap(genome, p, r, a):
            return genome[: p - 1] + a + genome[p - 1 + len(r):]

        for genome in ("GCTTTGC", "CAAA", "ATATAT", "CAACG", "AAAA"):
            groups: dict[tuple, list] = {}
            for p, r, a, kind in writings(genome):
                dlen = abs(len(r) - len(a))
                groups.setdefault((kind, dlen, hap(genome, p, r, a)), []).append(
                    (p, r, a)
                )
            for (kind, dlen, _), members in groups.items():
                p, r, a = normalize_variant(*members[0], genome)
                if (p, r, a) not in members:
                    continue  # 该组不含规范写法（例如起点无法重新锚定）
                lo, hi = equivalent_span(p, r, a, genome)
                if kind == "DEL":
                    positions = set()
                    for mp, mr, _ma in members:
                        positions.update(range(mp + 1, mp + len(mr)))
                else:
                    positions = {mp for mp, _mr, _ma in members}
                self.assertEqual(
                    (lo, hi), (min(positions), max(positions)),
                    f"{genome} {kind} canonical={(p, r, a)}",
                )


class IntervalCheckServiceTest(AppTestBase):
    def setUp(self):
        super().setUp()
        self.upload_ref(REF_V1)

    def check(self, start, end):
        with self.app.app_context():
            return services.check_interval(start, end)

    def test_requested_interval_without_expansion(self):
        # SNP 只占位点 8，区间 8-8 已完整覆盖
        self.upload_sample("S1", tsv("8\tA\tC\t0/1\n"))
        result = self.check(8, 8)
        self.assertEqual(result["requested"], (8, 8))
        self.assertEqual(result["expanded"], (8, 8))
        self.assertEqual(result["uncovered"], [])

    def test_insertion_written_outside_interval_forces_expansion(self):
        # 规范坐标为 (2,C,CT)，等价写法锚在 2..5；区间只含 5 位时
        # 规范坐标在区间外，但等价写法在区间内，必须扩展到 2
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        result = self.check(5, 5)
        self.assertEqual(result["expanded"], (2, 5))
        self.assertEqual([(h["pos"], h["ref"], h["alt"]) for h in result["hits"]],
                         [(2, "C", "CT")])
        self.assertEqual(len(result["uncovered"]), 1)

    def test_deletion_range_uses_deleted_bases(self):
        # 规范 (2,CT,C) 实际删除 3..5；区间 2-2 不应命中，2-3 应命中并扩到 5
        self.upload_sample("S1", tsv("3\tTT\tT\t0/1\n"))
        self.assertEqual(self.check(2, 2)["hits"], [])
        result = self.check(2, 3)
        self.assertEqual(result["expanded"], (2, 5))

    def test_expansion_iterates_until_stable(self):
        # 两个变异通过范围互相牵动：A 锚区 2..5，B 在 6..6
        # 区间 5-5 -> 收录 A 后扩到 2..5，不碰 B；
        # 区间 5-6 -> 收录 A 和 B，稳定在 2..6。
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n6\tG\tA\t1/1\n"))
        self.assertEqual(self.check(5, 6)["expanded"], (2, 6))
        self.assertEqual(len(self.check(5, 6)["hits"]), 2)

    def test_only_current_valid_detections_participate(self):
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        self.upload_sample("S2", tsv("5\tT\tTT\t1/1\n"))
        sid = self._sample_id("S1")
        self.client.post(f"/samples/{sid}/delete", follow_redirects=True)
        result = self.check(5, 5)
        # S2 仍有效，范围依然参与扩展
        self.assertEqual(result["expanded"], (2, 5))
        self.assertEqual(result["hits"][0]["n_samples"], 1)

        sid = self._sample_id("S2")
        self.client.post(f"/samples/{sid}/delete", follow_redirects=True)
        result = self.check(5, 5)
        # 失去最后一份有效检出后，不再参与扩展
        self.assertEqual(result["expanded"], (5, 5))
        self.assertEqual(result["hits"], [])

    def test_reupload_drops_obsolete_range(self):
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        self.assertEqual(self.check(5, 5)["expanded"], (2, 5))
        self.upload_sample("S1", tsv("8\tA\tC\t0/1\n"))
        self.assertEqual(self.check(5, 5)["expanded"], (5, 5))

    def test_invalid_intervals(self):
        for start, end in (("0", "5"), ("5", "0"), ("x", "5"), ("6", "5")):
            with self.assertRaises(ValidationError):
                self.check(start, end)
        with self.assertRaises(ValidationError):
            self.check(1, len(REF_V1) + 1)

    def test_stale_reference_form_rejected(self):
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        old_id = self._current_ref_id()
        self.upload_ref(REF_V2)
        with self.app.app_context():
            with self.assertRaises(services.StaleReferenceError):
                services.check_interval("5", "5", old_id)

    def test_server_side_check_without_expected_version_uses_current(self):
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        # 未携带 reference_id 时按当前参考查询（API/直连场景）
        with self.app.app_context():
            result = services.check_interval("5", "5")
        self.assertEqual(result["expanded"], (2, 5))

    def _sample_id(self, name):
        with self.app.app_context():
            from app.db import get_db
            return get_db().execute(
                "SELECT id FROM samples WHERE name=?", (name,)
            ).fetchone()["id"]

    def _current_ref_id(self):
        with self.app.app_context():
            from app.db import get_db
            return get_db().execute(
                "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
            ).fetchone()["id"]


class IntervalCheckNoReferenceTest(AppTestBase):
    def test_page_without_reference(self):
        page = self.client.get("/interval").get_data(as_text=True)
        self.assertIn("尚未上传参考序列", page)

    def test_check_without_reference_reports_error(self):
        page = self.client.get("/interval/check?start=1&end=2")
        self.assertEqual(page.status_code, 400)
        self.assertIn("请先上传参考序列", page.get_data(as_text=True))


class IntervalCheckPageTest(AppTestBase):
    def setUp(self):
        super().setUp()
        self.upload_ref(REF_V1)

    def test_page_lists_expansion(self):
        self.upload_sample("S1", tsv("5\tT\tTT\t0/1\n"))
        ref_id = self._current_ref_id()
        page = self.client.get(
            f"/interval/check?start=5&end=5&reference_id={ref_id}"
        ).get_data(as_text=True)
        self.assertIn("最小扩展区间", page)
        self.assertIn("2 – 5", page)
        self.assertIn("否", page)  # 原区间未完整覆盖

    def test_stale_reference_form_shows_prompt(self):
        old_id = self._current_ref_id()
        self.upload_ref(REF_V2)
        page = self.client.get(
            f"/interval/check?start=5&end=5&reference_id={old_id}"
        )
        self.assertEqual(page.status_code, 409)
        text = page.get_data(as_text=True)
        self.assertIn("请重新选择参考", text)

    def test_bad_interval_reports_error(self):
        page = self.client.get("/interval/check?start=9&end=2")
        self.assertEqual(page.status_code, 400)
        self.assertIn("不能大于终点", page.get_data(as_text=True))

    def test_index_links_to_interval_check(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("区间核对", page)
        self.assertIn("/interval", page)

    def _current_ref_id(self):
        with self.app.app_context():
            from app.db import get_db
            return get_db().execute(
                "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
            ).fetchone()["id"]


class LegacyDatabaseMigrationTest(unittest.TestCase):
    """旧库（没有 variant_ranges 表）启动时补建等价范围。"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.tmp.name, "legacy.db")

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    def _build_legacy_db(self):
        conn = sqlite3.connect(self.db_path)
        conn.executescript(
            """
            CREATE TABLE reference_sequences (
                id INTEGER PRIMARY KEY, sequence TEXT NOT NULL,
                length INTEGER NOT NULL,
                uploaded_at TEXT NOT NULL DEFAULT (datetime('now'))
            );
            CREATE TABLE samples (
                id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE,
                reference_id INTEGER REFERENCES reference_sequences(id),
                status TEXT NOT NULL, error TEXT, raw_tsv TEXT NOT NULL,
                uploaded_at TEXT NOT NULL DEFAULT (datetime('now')),
                revalidated_at TEXT
            );
            CREATE TABLE sample_variants (
                id INTEGER PRIMARY KEY,
                sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
                pos INTEGER NOT NULL, ref TEXT NOT NULL, alt TEXT NOT NULL,
                gt TEXT NOT NULL, kind TEXT NOT NULL,
                orig_pos INTEGER NOT NULL, orig_ref TEXT NOT NULL,
                orig_alt TEXT NOT NULL
            );
            """
        )
        conn.execute(
            "INSERT INTO reference_sequences (id, sequence, length) VALUES (1, ?, ?)",
            (REF_V1, len(REF_V1)),
        )
        conn.execute(
            "INSERT INTO samples (id, name, reference_id, status, raw_tsv) "
            "VALUES (1, 'OLD', 1, 'VALID', '5\tT\tTT\t0/1\n')"
        )
        conn.execute(
            "INSERT INTO sample_variants "
            "(sample_id, pos, ref, alt, gt, kind, orig_pos, orig_ref, orig_alt) "
            "VALUES (1, 2, 'C', 'CT', '0/1', 'INS', 5, 'T', 'TT')"
        )
        conn.commit()
        conn.close()

    def test_legacy_rows_get_ranges_on_startup(self):
        self._build_legacy_db()
        os.environ["DATABASE_PATH"] = self.db_path
        app = create_app()
        app.config.update(TESTING=True, SECRET_KEY="test")
        with app.app_context():
            result = services.check_interval("5", "5")
        self.assertEqual(result["expanded"], (2, 5))
        self.assertEqual(result["hits"][0]["pos"], 2)
        # 重验/新增写入时也保持范围存在
        client = app.test_client()
        client.post(
            "/samples",
            data={
                "name": "NEW",
                "sample": (io.BytesIO(tsv("6\tG\tA\t0/1\n")), "new.tsv"),
            },
            content_type="multipart/form-data",
            follow_redirects=True,
        )
        with app.app_context():
            from app.db import get_db
            rows = get_db().execute(
                "SELECT pos, ref, alt, lo, hi FROM variant_ranges "
                "WHERE reference_id=1 ORDER BY pos"
            ).fetchall()
        self.assertEqual(
            [(r["pos"], r["ref"], r["alt"], r["lo"], r["hi"]) for r in rows],
            [(2, "C", "CT", 2, 5), (6, "G", "A", 6, 6)],
        )

    def test_migration_keeps_detections_and_is_idempotent(self):
        self._build_legacy_db()
        os.environ["DATABASE_PATH"] = self.db_path
        app = create_app()
        app.config.update(TESTING=True, SECRET_KEY="test")
        with app.app_context():
            from app.db import get_db
            self.assertEqual(
                get_db().execute("SELECT COUNT(*) c FROM sample_variants")
                .fetchone()["c"],
                1,
            )
        # 第二次启动不应重复插入
        app2 = create_app()
        app2.config.update(TESTING=True, SECRET_KEY="test")
        with app2.app_context():
            from app.db import get_db
            count = get_db().execute(
                "SELECT COUNT(*) c FROM variant_ranges"
            ).fetchone()["c"]
        self.assertEqual(count, 1)


if __name__ == "__main__":
    unittest.main()
