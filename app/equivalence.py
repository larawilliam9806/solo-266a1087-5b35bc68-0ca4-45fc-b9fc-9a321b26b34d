"""规范变异在当前参考上的「等价写法范围」计算（纯函数，不依赖 Flask）。

同一次 SNP / 插入 / 缺失在重复序列中可能被写成多个合法坐标：
左对齐后的规范坐标可能落在下游区间之外，而某个等价写法却在区间之内。
为支持区间核对，这里求出变异覆盖的参考位置范围，口径为：

  * 单碱基替换：按位点（单个位置）。
  * 插入：按左锚位点（插入不消耗参考碱基，锚在哪个参考位点就记哪个）。
  * 缺失：按实际删去的碱基（锚碱基不计入范围）。

插入的等价写法是在重复区内左右滑动、插入序列随之旋转得到的；缺失的
等价写法是左右滑动删除窗口。两种滑动都是逐位可逆的，因此范围是一个
连续闭区间，用滑动边界即可求出，无需枚举全部写法。
"""

from __future__ import annotations


def equivalent_span(pos: int, ref: str, alt: str, genome: str) -> tuple[int, int]:
    """返回变异在 ``genome`` 上所有合法等价写法覆盖的最小/最大位置。

    ``pos`` 为 1-based 位点，``ref``/``alt`` 为等位基因。若传入的不是规范写法，
    会先按参考左对齐归一化，再按规范形状计算范围。返回 ``(lo, hi)``，闭区间。
    """
    from .normalization import variant_kind, normalize_variant

    kind = variant_kind(ref, alt)
    if kind is None:
        # 非平台支持的形状：保守地只覆盖自身
        return pos, pos + len(ref) - 1

    pos, ref, alt = normalize_variant(pos, ref, alt, genome)

    n = len(genome)
    if kind == "SNP":
        # 单碱基替换只有一个合法写法
        return pos, pos
    if kind == "INS":
        return _insertion_span(pos, alt[1:], genome, n)
    return _deletion_span(pos, len(ref) - 1, genome, n)


def _insertion_span(pos: int, inserted: str, genome: str, n: int) -> tuple[int, int]:
    """插入：枚举所有可等价书写的左锚位点。

    锚点为 ``a`` 的写法把 ``genome[:a] + inserted + genome[a:]`` 作为单倍型。
    向右滑一位合法当且仅当插入序列首碱基等于右侧参考碱基，插入序列循环
    左移一位；向左滑一位对称。两端即为所有等价锚点的最小/最大位点。
    """
    hi = pos
    a, seq = pos, inserted
    while a < n and seq and seq[0] == genome[a]:
        seq = seq[1:] + genome[a]
        a += 1
    hi = a

    lo = pos
    a, seq = pos, inserted
    while a > 1 and seq and seq[-1] == genome[a - 1]:
        seq = genome[a - 1] + seq[:-1]
        a -= 1
    lo = a
    return lo, hi


def _deletion_span(pos: int, deleted_len: int, genome: str, n: int) -> tuple[int, int]:
    """缺失：枚举所有可等价书写的左锚位点，再换算成实际删去的碱基范围。

    锚点 ``p`` 的写法删除 1-based 位置 ``p+1 .. p+deleted_len``。
    向右滑一位合法当且仅当离开删除窗口的碱基与进入窗口的碱基相同，
    即 ``genome[p] == genome[p + deleted_len]``（0-based）。
    """
    hi_anchor = pos
    p = pos
    while p + deleted_len < n and genome[p] == genome[p + deleted_len]:
        p += 1
    hi_anchor = p

    lo_anchor = pos
    p = pos
    while p > 1 and genome[p - 1] == genome[p + deleted_len - 1]:
        p -= 1
    lo_anchor = p

    # 锚点 p 删除的位置是 p+1 .. p+deleted_len，锚碱基本身不计入
    return lo_anchor + 1, hi_anchor + deleted_len
