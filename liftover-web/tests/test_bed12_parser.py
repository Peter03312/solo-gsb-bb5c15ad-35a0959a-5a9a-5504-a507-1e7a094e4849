"""BED12 解析器测试：12 列、块几何、thick 端点、错误汇总行号。"""

from app.bed12 import parse_bed12


GOOD = (
    "chr1\t0\t340\ttx\t900\t+\t60\t170\t30,120,200\t3\t16,120,40,\t0,56,300,\n"
)


def test_good_bed12_parses_exons_and_fields():
    result = parse_bed12(GOOD)
    assert result.errors == []
    rec = result.records[0]
    assert (rec.chrom, rec.start, rec.end) == ("chr1", 0, 340)
    assert (rec.name, rec.score, rec.strand) == ("tx", 900, "+")
    assert rec.block_count == 3
    assert rec.block_sizes == (16, 120, 40)
    assert rec.block_starts == (0, 56, 300)
    assert rec.exons() == [(0, 16), (56, 176), (300, 340)]
    assert rec.item_rgb == "30,120,200"


def test_bed12_requires_exactly_twelve_columns():
    text = "chr1\t0\t10\ttx\t0\t+\t0\t0\t0\t1\t10,\n"  # 11 列
    result = parse_bed12(text)
    assert result.records == []
    assert any("12 列" in e and "第 1 行" in e for e in result.errors)


def test_bed12_bad_columns_collected_by_line():
    text = (
        "chr1\t0\t10\ttx\t0\t+\t0\t0\t0\t1\t10,\n"           # 11 列
        "chr1\tx\t10\ttx\t0\t+\t0\t0\t0\t1\t10,\t0,\n"       # 非整数
        "chr1\t10\t0\ttx\t0\t+\t0\t0\t0\t1\t10,\t0,\n"       # start>=end
        "chr1\t0\t10\ttx\t.\t+\t0\t0\t0\t1\t10,\t0,\n"       # score 非整数
        "chr1\t0\t10\ttx\t1001\t+\t0\t0\t0\t1\t10,\t0,\n"    # score 越界
        "chr1\t0\t10\ttx\t0\t.\t0\t0\t0\t1\t10,\t0,\n"       # strand 非法
    )
    errors = parse_bed12(text).errors
    for line_no in range(1, 7):
        assert any(f"第 {line_no} 行" in e for e in errors), errors


def test_bed12_block_count_mismatch():
    text = "chr1\t0\t60\ttx\t0\t+\t0\t0\t0\t3\t10,20,30,\t0,10,\n"
    errors = parse_bed12(text).errors
    assert any("blockCount" in e and "blockStarts" in e and "不一致" in e for e in errors)


def test_bed12_first_block_start_must_be_zero():
    text = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t50,40,\t10,60,\n"
    errors = parse_bed12(text).errors
    assert any("第一个 blockStart 必须为 0" in e for e in errors)


def test_bed12_overlapping_blocks_rejected():
    # [0,60) 与 [50,100) 重叠
    text = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t60,50,\t0,50,\n"
    errors = parse_bed12(text).errors
    assert any("重叠" in e for e in errors)


def test_bed12_adjacent_blocks_allowed():
    text = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t50,50,\t0,50,\n"
    assert parse_bed12(text).errors == []


def test_bed12_last_block_must_reach_chrom_end():
    text = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t40,40,\t0,50,\n"
    errors = parse_bed12(text).errors
    assert any("顶到 chromEnd" in e for e in errors)


def test_bed12_zero_size_block_rejected():
    text = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t100,0,\t0,0,\n"
    errors = parse_bed12(text).errors
    assert any("正整数" in e for e in errors)


def test_bed12_thick_endpoints_must_lie_in_exons():
    # 外显子 [0,10) 与 [90,100)；thick 50/60 位于内含子
    text = "chr1\t0\t100\ttx\t0\t+\t50\t60\t0\t2\t10,10,\t0,90,\n"
    errors = parse_bed12(text).errors
    assert any("thickStart" in e for e in errors)
    assert any("thickEnd" in e for e in errors)


def test_bed12_thick_endpoint_at_exon_boundary_is_outside():
    # thickStart=10 是第一个外显子的末端（开区间端点），不算落在外显子内
    text = "chr1\t0\t100\ttx\t0\t+\t10\t95\t0\t2\t10,10,\t0,90,\n"
    errors = parse_bed12(text).errors
    assert any("thickStart" in e for e in errors)


def test_bed12_zero_length_thick_must_also_lie_in_exon():
    text = "chr1\t0\t100\ttx\t0\t+\t50\t50\t0\t2\t10,10,\t0,90,\n"
    errors = parse_bed12(text).errors
    assert any("thickStart" in e for e in errors) and any("thickEnd" in e for e in errors)


def test_bed12_item_rgb_validation():
    bad1 = "chr1\t0\t100\ttx\t0\t+\t0\t0\t256,0,0\t2\t50,50,\t0,50,\n"
    assert any("itemRgb" in e for e in parse_bed12(bad1).errors)
    bad2 = "chr1\t0\t100\ttx\t0\t+\t0\t0\tx,y,z\t2\t50,50,\t0,50,\n"
    assert any("itemRgb" in e for e in parse_bed12(bad2).errors)
    good = "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t50,50,\t0,50,\n"
    assert parse_bed12(good).errors == []


def test_bed12_comments_and_tracks_skipped_empty_is_error():
    text = "# comment\ntrack name=x\nbrowser position chr1:1\n"
    result = parse_bed12(text)
    assert result.records == []
    assert result.errors
