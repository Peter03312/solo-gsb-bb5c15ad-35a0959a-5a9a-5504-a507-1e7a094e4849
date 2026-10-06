"""BED12 转录本整体迁移引擎（单段 / 双段）。

与把区间切碎的 BED3 迁移不同，转录本只能**整体**迁移：一条 chain 路径
（双段时为“第一段 chain × 第二段 chain”组合）只有在**每一个外显子都
完整映射**（整个外显子落在单个对齐块内，不允许跨缺口）、全部目标外显子
落在**同一目标染色体**、相对方向**一致**、且目标块按源顺序**保持或
整体反转**、**互不重叠**时才成为候选。任何缺口或结构冲突只淘汰该条
路径，不产出“部分转录本”。

- 多条完整路径并存时全部保留，行状态标为 ``ambiguous``，不自动选取；
- 没有任何完整路径时行状态为 ``unmapped``，原因中列出失败发生的位置
  （双段另标 ``failed_stage``），**绝不**导出伪 BED12 行。

候选组装为输出 BED12 时：

- 按目标坐标升序重排块，重算 chromStart/chromEnd、blockStarts；
  blockSizes 长度在仿射映射中保持不变；
- thickStart/thickEnd 按“映射后外显子内的编码部分”重算；非编码记录
  （零长 thick）保持零长 thick（放在新的 chromStart）；
- 整体反向（各段相对方向之积为 -1）时翻转 strand，否则保留；
- name / score / itemRgb 原样保留。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .bed12 import Bed12Record
from .chain import Chain
from .engine import Fragment, LiftOverEngine


@dataclass
class ExonLift:
    """单个外显子的映射结果（坐标均为零基左闭右开的正向帧）。"""

    exon_index: int            # 源外显子序号（0 基，按源顺序）
    src_start: int
    src_end: int
    dest_start: int
    dest_end: int
    mid_start: int | None = None   # 双段时的中间坐标（第一段目标）
    mid_end: int | None = None
    # 源坐标 -> 目标坐标的块内前向映射参数：dest = fa * src + fb
    # （fa ∈ {+1,-1}，整体正向即 +1、反向即 -1）。thick 端点必须用它
    # 所在外显子的参数映射，不能复用第一个外显子的参数。
    fa: int = 1
    fb: int = 0


@dataclass
class Bed12Candidate:
    candidate_index: int
    dest_chrom: str
    dest_strand: str
    reversed: bool                  # 整体反向：输出时翻转源 strand
    exons: list[ExonLift]           # 按源顺序排列的外显子映射
    # ---- 组装好的输出 BED12 坐标列 ----
    out_start: int = 0
    out_end: int = 0
    out_strand: str = "+"
    out_thick_start: int = 0
    out_thick_end: int = 0
    out_block_sizes: tuple[int, ...] = ()
    out_block_starts: tuple[int, ...] = ()
    # ---- 路径溯源 ----
    chain_id: str = ""
    chain1_id: str = ""
    chain2_id: str = ""
    mid_chrom: str = ""
    mid_strand: str = ""


@dataclass
class Bed12Row:
    line_no: int
    src_chrom: str
    src_start: int
    src_end: int
    name: str
    status: str                     # mapped / ambiguous / unmapped
    reason: str = ""
    failed_stage: str = ""          # 仅双段无映射：segment1 / segment2
    candidates: list[Bed12Candidate] = field(default_factory=list)


@dataclass(frozen=True)
class _FullHit:
    """一个外显子在某条 chain 上的完整命中（整条外显子位于单个对齐块内）。"""

    chain_uid: int
    chain_id: str
    dest_chrom: str
    dest_strand: str
    d_lo: int
    d_hi: int
    frag: Fragment


@dataclass(frozen=True)
class _ExonLeg:
    """双段路径中一个外显子的两段命中与坐标。"""

    c1_uid: int
    c1_id: str
    c2_uid: int
    c2_id: str
    mid_chrom: str
    mid_strand: str
    dest_chrom: str
    dest_strand: str
    m_lo: int
    m_hi: int
    d_lo: int
    d_hi: int
    f1: Fragment
    f2: Fragment


def _full_hits(engine: LiftOverEngine, chrom: str, lo: int, hi: int) -> list[_FullHit]:
    """找出能**完整**包含外显子 [lo, hi) 的全部 chain 命中。

    跨缺口时引擎会把结果拆成多个片段或缺失端点，此类候选一律排除。
    """
    hits: list[_FullHit] = []
    for cand in engine.candidates_for(chrom, lo, hi):
        frags = cand.fragments
        if (
            len(frags) == 1
            and frags[0].src_start == lo
            and frags[0].src_end == hi
        ):
            hits.append(
                _FullHit(
                    chain_uid=cand.chain_uid,
                    chain_id=cand.chain_id,
                    dest_chrom=cand.dest_chrom,
                    dest_strand=cand.dest_strand,
                    d_lo=frags[0].dest_start,
                    d_hi=frags[0].dest_end,
                    frag=frags[0],
                )
            )
    return hits


def _frag_forward(frag: Fragment) -> tuple[int, int]:
    """从引擎逆映射参数（src = back_a*dest + back_b）推导前向映射
    dest = fa*src + fb 的参数。"""
    return frag.back_a, -frag.back_a * frag.back_b


def _summarize(notes: list[str], cap: int = 4) -> str:
    dedup = list(dict.fromkeys(notes))
    if len(dedup) > cap:
        return "；".join(dedup[:cap]) + f"；等 {len(dedup)} 处"
    return "；".join(dedup)


def _structure_check(exon_hits: list[tuple[int, int, int, int, int]]) -> str:
    """对一条路径的外显子目标区间做结构校验。

    入参为按源顺序的 ``(exon_index, d_lo, d_hi, src_lo, src_hi)``。
    返回非空字符串即淘汰原因；空串表示通过。校验：

    - 相邻外显子目标坐标严格保持源顺序或严格整体反转（不允许并列）；
    - 目标块互不重叠（允许相接）。
    """
    order = 0
    prev_lo = prev_hi = None
    for _ei, d_lo, d_hi, _s_lo, _s_hi in exon_hits:
        if prev_lo is not None:
            diff = d_lo - prev_lo
            if diff == 0:
                return f"外显子 {_ei} 与前一外显子映射到同一目标起点，源顺序无法保持"
            sign = 1 if diff > 0 else -1
            if order == 0:
                order = sign
            elif order != sign:
                return (
                    "外显子目标坐标既非按源顺序保持也非整体反转"
                    "（顺序在前/后之间发生变化）"
                )
            if sign > 0 and d_lo < prev_hi:
                return f"外显子 {_ei} 的目标区间与前一外显子重叠"
        prev_lo, prev_hi = d_lo, d_hi
    return ""


def _assemble(
    rec: Bed12Record,
    cand_index: int,
    dest_chrom: str,
    dest_strand: str,
    reversed_all: bool,
    exons: list[ExonLift],
    chain_id: str = "",
    chain1_id: str = "",
    chain2_id: str = "",
    mid_chrom: str = "",
    mid_strand: str = "",
) -> Bed12Candidate:
    """把通过结构校验的外显子映射组装成输出 BED12 候选。"""
    cand = Bed12Candidate(
        candidate_index=cand_index,
        dest_chrom=dest_chrom,
        dest_strand=dest_strand,
        reversed=reversed_all,
        exons=exons,
        chain_id=chain_id,
        chain1_id=chain1_id,
        chain2_id=chain2_id,
        mid_chrom=mid_chrom,
        mid_strand=mid_strand,
    )

    # 按目标坐标升序排列外显子（反向路径这里即整体反转后的顺序）
    ordered = sorted(exons, key=lambda e: e.dest_start)
    cand.out_start = ordered[0].dest_start
    cand.out_end = max(e.dest_end for e in exons)

    sizes: list[int] = []
    starts: list[int] = []
    cursor = cand.out_start
    for e in ordered:
        size = e.dest_end - e.dest_start   # 仿射映射保持长度
        sizes.append(size)
        starts.append(e.dest_start - cursor)
        cursor = e.dest_end
    cand.out_block_sizes = tuple(sizes)
    cand.out_block_starts = tuple(starts)

    cand.out_strand = ("-" if rec.strand == "+" else "+") if reversed_all else rec.strand

    if rec.thick_start == rec.thick_end:
        # 非编码记录：保持零长 thick
        cand.out_thick_start = cand.out_start
        cand.out_thick_end = cand.out_start
    else:
        # thick 端点必须使用其**所在外显子**的块内前向映射参数
        t1 = _map_thick_point(exons, rec.thick_start)
        t2 = _map_thick_point(exons, rec.thick_end)
        cand.out_thick_start = min(t1, t2)
        cand.out_thick_end = max(t1, t2)

    return cand


def _map_thick_point(exons: list[ExonLift], x: int) -> int:
    """按 thick 端点所在外显子的仿射参数做源→目标点映射。"""
    for e in exons:
        if e.src_start <= x < e.src_end:
            return e.fa * x + e.fb
    # 解析器已保证 thick 端点落在某个外显子内，正常不会走到这里
    raise ValueError(f"thick 端点 {x} 不落在任何已映射外显子内")


class Bed12Engine:
    """BED12 转录本单段整体迁移。"""

    def __init__(self, chains: list[Chain], side: str) -> None:
        self.engine = LiftOverEngine(chains, side)
        self.side = side

    def convert(self, rec: Bed12Record) -> Bed12Row:
        row = Bed12Row(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            name=rec.name,
            status="unmapped",
        )
        exons = rec.exons()

        if not self.engine.has_chrom(rec.chrom):
            row.reason = (
                f"染色体 {rec.chrom} 在 chain 文件的源侧不存在，"
                f"转录本的全部 {len(exons)} 个外显子均无法迁移"
            )
            return row

        # 每个外显子 -> 能完整映射它的 chain 命中
        per_exon: list[list[_FullHit]] = []
        notes: list[str] = []
        for ei, (lo, hi) in enumerate(exons, start=1):
            hits = _full_hits(self.engine, rec.chrom, lo, hi)
            per_exon.append(hits)
            if not hits:
                notes.append(
                    f"外显子 {ei} [{lo},{hi}) 无法完整落在任何对齐块内"
                    f"（跨链内缺口、落在未对齐区域或超出 chain 范围）"
                )

        # 以 chain_uid 归组路径：一条 chain 必须完整覆盖全部外显子
        chain_uids: dict[int, _FullHit] = {}
        for hits in per_exon:
            for h in hits:
                chain_uids[h.chain_uid] = h
        candidates: list[Bed12Candidate] = []
        reject_notes: list[str] = []
        ci = 0
        for uid in chain_uids:
            # 精确取该 chain 逐外显子的完整命中；任一外显子缺失即淘汰
            seq: list[_FullHit] = []
            for exon_hits in per_exon:
                h = next((x for x in exon_hits if x.chain_uid == uid), None)
                if h is None:
                    seq = []
                    break
                seq.append(h)
            if not seq:
                continue

            chroms = {h.dest_chrom for h in seq}
            if len(chroms) != 1:
                reject_notes.append(
                    f"chain {seq[0].chain_id}：各外显子目标染色体不一致"
                    f"（{', '.join(sorted(chroms))}）"
                )
                continue
            signs = {h.frag.back_a for h in seq}
            if len(signs) != 1:
                reject_notes.append(
                    f"chain {seq[0].chain_id}：各外显子映射方向不一致"
                )
                continue
            sign = next(iter(signs))
            struct = _structure_check(
                [(i, h.d_lo, h.d_hi, exons[i][0], exons[i][1])
                 for i, h in enumerate(seq)]
            )
            if struct:
                reject_notes.append(f"chain {seq[0].chain_id}：{struct}")
                continue

            ci += 1
            exon_lifts = [
                ExonLift(
                    exon_index=i,
                    src_start=exons[i][0],
                    src_end=exons[i][1],
                    dest_start=h.d_lo,
                    dest_end=h.d_hi,
                    fa=_frag_forward(h.frag)[0],
                    fb=_frag_forward(h.frag)[1],
                )
                for i, h in enumerate(seq)
            ]
            cand = _assemble(
                rec, ci, seq[0].dest_chrom, seq[0].dest_strand,
                sign == -1, exon_lifts,
                chain_id=seq[0].chain_id,
            )
            candidates.append(cand)

        if not candidates:
            bits = list(notes)
            if not bits:
                bits.append("没有任何一条 chain 能完整覆盖全部外显子")
            if reject_notes:
                bits.append("满足完整映射的路径经结构校验被淘汰："
                            + _summarize(reject_notes))
            row.reason = _summarize(bits, cap=6)
            return row

        row.candidates = candidates
        if len(candidates) > 1:
            row.status = "ambiguous"
            row.reason = (
                f"存在 {len(candidates)} 条完整路径（chain id: "
                f"{', '.join(c.chain_id for c in candidates)}），"
                "全部外显子均完整映射到不同目标位置；按要求不自动挑选，"
                "全部候选导出，需人工裁决"
            )
        else:
            row.status = "mapped"
        return row

    def convert_all(self, records: list[Bed12Record]) -> list[Bed12Row]:
        return [self.convert(r) for r in records]


class Bed12Engine2:
    """BED12 转录本双段整体迁移（chain 文件 1 × chain 文件 2）。"""

    def __init__(
        self,
        chains1: list[Chain],
        side1: str,
        chains2: list[Chain],
        side2: str,
    ) -> None:
        self.engine1 = LiftOverEngine(chains1, side1)
        self.engine2 = LiftOverEngine(chains2, side2)
        self.side1 = side1
        self.side2 = side2

    def convert(self, rec: Bed12Record) -> Bed12Row:
        row = Bed12Row(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            name=rec.name,
            status="unmapped",
        )
        exons = rec.exons()

        if not self.engine1.has_chrom(rec.chrom):
            row.failed_stage = "segment1"
            row.reason = (
                f"第一段无映射（chain 文件 1）：染色体 {rec.chrom} "
                f"在 chain 文件 1 的源侧不存在，转录本的全部 {len(exons)} 个外显子均无法迁移"
            )
            return row

        # 每个外显子枚举两条链上都完整命中的 (chain1, chain2) 组合
        per_exon_legs: list[list[_ExonLeg]] = []
        seg1_notes: list[str] = []
        seg2_notes: list[str] = []

        for ei, (lo, hi) in enumerate(exons, start=1):
            legs: list[_ExonLeg] = []
            hits1 = _full_hits(self.engine1, rec.chrom, lo, hi)
            if not hits1:
                seg1_notes.append(
                    f"外显子 {ei} [{lo},{hi}) 在第一段无法完整落在任何对齐块内"
                    f"（跨缺口、未对齐或超出 chain 文件 1 范围）"
                )
                per_exon_legs.append(legs)
                continue
            for h1 in hits1:
                hits2 = _full_hits(self.engine2, h1.dest_chrom, h1.d_lo, h1.d_hi)
                if not hits2:
                    if not self.engine2.has_chrom(h1.dest_chrom):
                        seg2_notes.append(
                            f"外显子 {ei} 的中间染色体 {h1.dest_chrom} "
                            f"（chain1 {h1.chain_id}）在 chain 文件 2 的源侧不存在"
                        )
                    else:
                        seg2_notes.append(
                            f"外显子 {ei} 经 chain1 {h1.chain_id} 映射到中间区间 "
                            f"{h1.dest_chrom}:[{h1.d_lo},{h1.d_hi}) 后，"
                            f"在第二段无法完整落在任何对齐块内"
                        )
                    continue
                for h2 in hits2:
                    legs.append(
                        _ExonLeg(
                            c1_uid=h1.chain_uid, c1_id=h1.chain_id,
                            c2_uid=h2.chain_uid, c2_id=h2.chain_id,
                            mid_chrom=h1.dest_chrom, mid_strand=h1.dest_strand,
                            dest_chrom=h2.dest_chrom, dest_strand=h2.dest_strand,
                            m_lo=h1.d_lo, m_hi=h1.d_hi,
                            d_lo=h2.d_lo, d_hi=h2.d_hi,
                            f1=h1.frag, f2=h2.frag,
                        )
                    )
            per_exon_legs.append(legs)

        # 路径以 (chain1_uid, 最终目标染色体) 为键：同一条第一段 chain
        # 把全部外显子送到同一个中间染色体；第二段允许每个外显子命中
        # 不同 chain，但最终必须落于同一目标染色体。
        route_keys: dict[tuple[int, str], list[_ExonLeg]] = {}
        for legs in per_exon_legs:
            for leg in legs:
                route_keys.setdefault((leg.c1_uid, leg.dest_chrom), []).append(leg)

        candidates: list[Bed12Candidate] = []
        reject_notes: list[str] = []
        ci = 0
        for (c1_uid, dest_chrom), _ in route_keys.items():
            seq: list[_ExonLeg] = []
            ok = True
            for exon_legs in per_exon_legs:
                leg = next(
                    (x for x in exon_legs
                     if x.c1_uid == c1_uid and x.dest_chrom == dest_chrom),
                    None,
                )
                if leg is None:
                    ok = False
                    break
                seq.append(leg)
            if not ok:
                continue

            mid_chroms = {leg.mid_chrom for leg in seq}
            if len(mid_chroms) != 1:
                reject_notes.append(
                    f"路径 chain1 {seq[0].c1_id} -> {dest_chrom}："
                    "各外显子的中间染色体不一致"
                )
                continue
            signs = {leg.f1.back_a * leg.f2.back_a for leg in seq}
            if len(signs) != 1:
                reject_notes.append(
                    f"路径 chain1 {seq[0].c1_id} -> {dest_chrom}："
                    "各外显子两段组合后的映射方向不一致"
                )
                continue
            sign = next(iter(signs))
            struct = _structure_check(
                [(i, leg.d_lo, leg.d_hi, exons[i][0], exons[i][1])
                 for i, leg in enumerate(seq)]
            )
            if struct:
                reject_notes.append(
                    f"路径 chain1 {seq[0].c1_id} -> {dest_chrom}：{struct}"
                )
                continue

            ci += 1
            exon_lifts = []
            for i, leg in enumerate(seq):
                fa1, fb1 = _frag_forward(leg.f1)
                fa2, fb2 = _frag_forward(leg.f2)
                exon_lifts.append(
                    ExonLift(
                        exon_index=i,
                        src_start=exons[i][0],
                        src_end=exons[i][1],
                        mid_start=leg.m_lo,
                        mid_end=leg.m_hi,
                        dest_start=leg.d_lo,
                        dest_end=leg.d_hi,
                        fa=fa2 * fa1,
                        fb=fa2 * fb1 + fb2,
                    )
                )
            cand = _assemble(
                rec, ci, dest_chrom, seq[0].dest_strand,
                sign == -1, exon_lifts,
                chain1_id=seq[0].c1_id,
                chain2_id="+".join(dict.fromkeys(leg.c2_id for leg in seq)),
                mid_chrom=seq[0].mid_chrom,
                mid_strand=seq[0].mid_strand,
            )
            candidates.append(cand)

        if not candidates:
            # 失败阶段：只要存在“第一段就无法完整命中”的外显子即记为
            # segment1（最先失败的一段）；否则所有外显子都通过第一段、
            # 仅在第二段失败/结构淘汰，记为 segment2。
            row.failed_stage = "segment1" if seg1_notes else "segment2"
            bits = seg1_notes + seg2_notes
            if reject_notes:
                bits.append("满足完整映射的路径经结构校验被淘汰："
                            + _summarize(reject_notes))
            if not bits:
                bits.append("没有任何 chain1 × chain2 组合能完整覆盖全部外显子")
            row.reason = _summarize(bits, cap=6)
            return row

        row.candidates = candidates
        if len(candidates) > 1:
            row.status = "ambiguous"
            row.reason = (
                f"存在 {len(candidates)} 条完整路径"
                f"（chain1 × chain2 组合），全部外显子均经两段完整映射；"
                "按要求不自动挑选，全部候选导出，需人工裁决"
            )
        else:
            row.status = "mapped"
        return row

    def convert_all(self, records: list[Bed12Record]) -> list[Bed12Row]:
        return [self.convert(r) for r in records]
