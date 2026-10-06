"""BED12（转录本）解析器。

坐标按零基、左闭右开解释。与只取前 3 列的 :func:`app.bed.parse_bed`
不同，BED12 模式要求**严格 12 列**并逐列校验，任何一行格式错误都会
带物理行号收集，最终整次转换被拒绝：

1. 恰好 12 列，``chromStart < chromEnd``，``strand`` 为 ``+``/``-``；
2. ``score`` 为 0..1000 的整数；
3. ``thickStart``/``thickEnd`` 为区间内整数且 ``thickStart <= thickEnd``；
   编码记录（``thickStart < thickEnd``）的两个端点必须落在某个外显子
   （块）内；非编码记录保持零长 thick（端点同样必须落在外显子内）；
4. ``blockCount`` 为正整数，且与 ``blockSizes``/``blockStarts`` 的
   逗号分隔项数一致；size 为正整数、start 为非负整数；
5. 块按起点严格有序、互不重叠（允许相接），第一块起点为 0，
   最后一块的末端恰好顶到 ``chromEnd``，所有块都在区间范围内；
6. ``itemRgb`` 为 ``0`` 或 ``r,g,b``（各分量 0..255）。

``name`` / ``score`` / ``itemRgb`` 等列的原始文本保留在记录中，
迁移后原样写回（坐标类列重新计算）。
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Bed12Record:
    line_no: int           # 物理文件行号，从 1 开始
    chrom: str
    start: int
    end: int
    name: str
    score: int
    strand: str
    thick_start: int
    thick_end: int
    item_rgb: str          # 原始文本（"0" 或 "r,g,b"）
    block_count: int
    block_sizes: tuple[int, ...]
    block_starts: tuple[int, ...]
    raw: str

    def exons(self) -> list[tuple[int, int]]:
        """返回各外显子在染色体上的零基左闭右开绝对坐标（按源顺序）。"""
        return [
            (self.start + st, self.start + st + size)
            for st, size in zip(self.block_starts, self.block_sizes)
        ]


@dataclass
class Bed12ParseResult:
    records: list[Bed12Record]
    errors: list[str]


def _uint(text: str, line_no: int, field: str) -> int:
    t = text.strip()
    if not t or not t.isdigit():
        raise ValueError(f"BED12 第 {line_no} 行: {field} 必须是非负整数，得到 {text!r}")
    return int(t)


def _split_csv(text: str, line_no: int, field: str) -> list[str]:
    """解析 BED 逗号序列：允许且仅允许末尾一个空项（UCSC 风格尾逗号）。"""
    if text == "":
        raise ValueError(f"BED12 第 {line_no} 行: {field} 不能为空")
    parts = text.split(",")
    if parts[-1] == "":
        parts = parts[:-1]
    if any(p == "" for p in parts):
        raise ValueError(f"BED12 第 {line_no} 行: {field} 含空项，得到 {text!r}")
    return parts


def _check_thick_in_exons(
    line_no: int,
    thick_start: int,
    thick_end: int,
    exons: list[tuple[int, int]],
) -> list[str]:
    errs: list[str] = []

    def inside(pos: int) -> bool:
        # 零基左闭右开：端点必须落在某个外显子 [lo, hi) 内
        return any(lo <= pos < hi for lo, hi in exons)

    if not inside(thick_start):
        errs.append(
            f"BED12 第 {line_no} 行: thickStart ({thick_start}) 必须落在某个外显子内"
        )
    if not inside(thick_end):
        errs.append(
            f"BED12 第 {line_no} 行: thickEnd ({thick_end}) 必须落在某个外显子内"
        )
    return errs


def parse_bed12(text: str) -> Bed12ParseResult:
    records: list[Bed12Record] = []
    errors: list[str] = []

    for idx, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        # 与 BED3 解析保持一致：browser/track 指令行与 # 注释行允许出现
        if line.startswith(("#", "browser", "track")):
            continue

        f = line.split()
        if len(f) != 12:
            errors.append(
                f"BED12 第 {idx} 行: BED12 必须恰好是 12 列，实际 {len(f)} 列"
            )
            continue

        line_errors: list[str] = []

        def err(msg: str) -> None:
            line_errors.append(f"BED12 第 {idx} 行: {msg}")

        chrom = f[0]
        if not chrom:
            err("chrom 不能为空")

        start = end = None
        try:
            start = _uint(f[1], idx, "chromStart")
            end = _uint(f[2], idx, "chromEnd")
        except ValueError as exc:
            line_errors.append(str(exc))
        if start is not None and end is not None and start >= end:
            err(f"chromStart ({start}) 必须小于 chromEnd ({end})")

        name = f[3]
        score = 0
        score_ok = False
        try:
            score = _uint(f[4], idx, "score")
            score_ok = True
            if not 0 <= score <= 1000:
                err(f"score 必须在 0..1000 之间，得到 {score}")
        except ValueError as exc:
            line_errors.append(str(exc))

        strand = f[5]
        if strand not in ("+", "-"):
            err(f"strand 只能是 '+' 或 '-'，得到 {strand!r}")

        thick_start = thick_end = None
        try:
            thick_start = _uint(f[6], idx, "thickStart")
            thick_end = _uint(f[7], idx, "thickEnd")
        except ValueError as exc:
            line_errors.append(str(exc))
        if (
            start is not None and end is not None
            and thick_start is not None and thick_end is not None
        ):
            if not (start <= thick_start <= end and start <= thick_end <= end):
                err(
                    "thickStart/thickEnd 必须位于 "
                    f"[{start},{end}) 内（得到 {thick_start}/{thick_end}）"
                )
            if thick_start > thick_end:
                err(
                    f"thickStart ({thick_start}) 不能大于 thickEnd ({thick_end})"
                )

        item_rgb = f[8]
        if item_rgb != "0":
            rgb_parts = item_rgb.split(",")
            if (
                len(rgb_parts) != 3
                or not all(p.isdigit() and 0 <= int(p) <= 255 for p in rgb_parts)
            ):
                err(
                    "itemRgb 必须是 '0' 或 'r,g,b'"
                    f"（各分量 0..255），得到 {item_rgb!r}"
                )

        block_count = None
        try:
            block_count = _uint(f[9], idx, "blockCount")
            if block_count < 1:
                err(f"blockCount 必须是正整数，得到 {block_count}")
        except ValueError as exc:
            line_errors.append(str(exc))

        sizes: list[int] = []
        starts_rel: list[int] = []
        size_ok = start_ok = False
        size_tokens: list[str] = []
        start_tokens: list[str] = []
        try:
            size_tokens = _split_csv(f[10], idx, "blockSizes")
        except ValueError as exc:
            line_errors.append(str(exc))
        try:
            start_tokens = _split_csv(f[11], idx, "blockStarts")
        except ValueError as exc:
            line_errors.append(str(exc))
        if size_tokens and all(p.isdigit() for p in size_tokens):
            for p in size_tokens:
                v = int(p)
                if v <= 0:
                    err(f"blockSizes 每项必须为正整数，出现 {v}")
                sizes.append(v)
            size_ok = all(v > 0 for v in sizes)
        elif size_tokens:
            err(f"blockSizes 必须是正整数序列，得到 {f[10]!r}")
        if start_tokens and all(p.isdigit() for p in start_tokens):
            starts_rel = [int(p) for p in start_tokens]
            start_ok = True
        elif start_tokens:
            err(f"blockStarts 必须是非负整数序列，得到 {f[11]!r}")

        if block_count is not None:
            if size_tokens and len(size_tokens) != block_count:
                err(
                    f"blockCount={block_count} 与 blockSizes "
                    f"项数 {len(size_tokens)} 不一致"
                )
            if start_tokens and len(start_tokens) != block_count:
                err(
                    f"blockCount={block_count} 与 blockStarts "
                    f"项数 {len(start_tokens)} 不一致"
                )

        # 块几何校验（坐标、项数都齐全时才有意义）
        if (
            start is not None and end is not None
            and block_count is not None
            and len(sizes) == block_count
            and len(starts_rel) == block_count
            and size_ok and start_ok
        ):
            if starts_rel[0] != 0:
                err(f"第一个 blockStart 必须为 0，得到 {starts_rel[0]}")
            prev_end_rel = -1
            for bi, (rs, sz) in enumerate(zip(starts_rel, sizes)):
                if rs < prev_end_rel:
                    err(
                        f"第 {bi + 1} 个块 (start={rs}, size={sz}) "
                        "与前一块重叠或未按起点升序排列"
                    )
                    break
                prev_end_rel = rs + sz
            if prev_end_rel > end - start:
                err(
                    "块超出 chromStart/chromEnd 范围"
                    f"（块末端 {start + prev_end_rel} > chromEnd {end}）"
                )
            elif prev_end_rel != end - start:
                err(
                    f"最后一个块的末端（{start + prev_end_rel}）"
                    f"必须恰好顶到 chromEnd（{end}），外显子不得覆盖不到转录本末端"
                )

            exons_abs = [
                (start + rs, start + rs + sz)
                for rs, sz in zip(starts_rel, sizes)
            ]
            if thick_start is not None and thick_end is not None:
                line_errors += _check_thick_in_exons(
                    idx, thick_start, thick_end, exons_abs
                )

        if line_errors:
            errors.extend(line_errors)
            continue

        records.append(
            Bed12Record(
                line_no=idx,
                chrom=chrom,
                start=start,
                end=end,
                name=name,
                score=score,
                strand=strand,
                thick_start=thick_start,
                thick_end=thick_end,
                item_rgb=item_rgb,
                block_count=block_count,
                block_sizes=tuple(sizes),
                block_starts=tuple(starts_rel),
                raw=line,
            )
        )

    if not records and not errors:
        errors.append("BED12 文件为空（没有任何有效数据行）")

    return Bed12ParseResult(records=records, errors=errors)
