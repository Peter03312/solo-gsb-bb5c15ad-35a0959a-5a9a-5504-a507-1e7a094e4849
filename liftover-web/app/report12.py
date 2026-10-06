"""BED12 转录本迁移结果导出：BED12 文本与预览 JSON。

- :func:`build_bed12`：原始结果导出。mapped 行写一条 BED12；
  ambiguous 行的**全部候选**各写一条（不自动选取）；unmapped 行
  **不写任何行**（绝不导出伪 BED12）。
- :func:`build_bed12_adjudicated`：裁决版导出。已裁决的歧义行只写
  被选候选；未裁决的歧义行仍写全部候选；mapped 行照写；unmapped 行
  仍不写。

坐标列按引擎组装结果生成；name / score / itemRgb 原样保留。
下载文件为标准 BED12（无额外注释/裁决列），候选身份通过预览与文件名
区分；为方便下游核对，候选行在 name 为空时仍保留原名（不另加后缀）。
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping

from .bed12 import Bed12Record
from .engine12 import Bed12Candidate, Bed12Row


def _format_csv(values) -> str:
    return ",".join(str(v) for v in values) + ","


def candidate_to_bed12_line(rec: Bed12Record, cand: Bed12Candidate) -> str:
    fields = [
        cand.dest_chrom,
        str(cand.out_start),
        str(cand.out_end),
        rec.name,
        str(rec.score),
        cand.out_strand,
        str(cand.out_thick_start),
        str(cand.out_thick_end),
        rec.item_rgb,
        str(len(cand.out_block_sizes)),
        _format_csv(cand.out_block_sizes),
        _format_csv(cand.out_block_starts),
    ]
    return "\t".join(fields)


def _rows_with_records(
    rows: list[Bed12Row],
    records_by_line: dict[int, Bed12Record],
) -> list[tuple[Bed12Row, Bed12Record]]:
    return [(r, records_by_line[r.line_no]) for r in rows]


def build_bed12(
    rows: list[Bed12Row],
    records: list[Bed12Record],
) -> str:
    """原始 BED12 结果（全部歧义候选都导出；unmapped 行不导出）。"""
    by_line = {r.line_no: r for r in records}
    lines: list[str] = []
    for row, rec in _rows_with_records(rows, by_line):
        if not row.candidates:
            continue
        for cand in row.candidates:
            lines.append(candidate_to_bed12_line(rec, cand))
    return "\n".join(lines) + ("\n" if lines else "")


def build_bed12_adjudicated(
    rows: list[Bed12Row],
    records: list[Bed12Record],
    decisions: Mapping[int, int],
) -> str:
    """裁决版 BED12：已裁决歧义行只写被选候选，未裁决的写全部候选。"""
    by_line = {r.line_no: r for r in records}
    lines: list[str] = []
    for row, rec in _rows_with_records(rows, by_line):
        if not row.candidates:
            continue
        if row.status == "ambiguous":
            picked = decisions.get(row.line_no)
            if picked is not None and 1 <= picked <= len(row.candidates):
                lines.append(candidate_to_bed12_line(rec, row.candidates[picked - 1]))
                continue
        for cand in row.candidates:
            lines.append(candidate_to_bed12_line(rec, cand))
    return "\n".join(lines) + ("\n" if lines else "")


def _candidate_dict(cand: Bed12Candidate) -> dict[str, Any]:
    return asdict(cand)


def rows_to_bed12_preview(rows: list[Bed12Row]) -> dict[str, Any]:
    counts = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}
    out_rows = []
    exported = 0
    for r in rows:
        counts[r.status] += 1
        exported += len(r.candidates)
        out_rows.append(
            {
                "line_no": r.line_no,
                "src_chrom": r.src_chrom,
                "src_start": r.src_start,
                "src_end": r.src_end,
                "name": r.name,
                "status": r.status,
                "failed_stage": r.failed_stage,
                "reason": r.reason,
                "candidates": [_candidate_dict(c) for c in r.candidates],
            }
        )
    return {
        "total": len(rows),
        "counts": counts,
        "exported_bed12_lines": exported,
        "rows": out_rows,
    }
