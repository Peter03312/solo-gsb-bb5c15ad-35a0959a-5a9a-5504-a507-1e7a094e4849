"""BED12 HTTP 接口测试：预览、下载 BED12、错误整次拒绝、裁决版本。"""

import os

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post_single(bed_name: str = "example_bed12.bed", side: str = "target",
                 chain_name: str = "demo.chain"):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, chain_name), "rb") as cf:
        return client.post(
            "/api/bed12/liftover",
            files={"bed": (bed_name, bf, "text/plain"),
                   "chain": (chain_name, cf, "text/plain")},
            data={"side": side},
        )


def _post_double(bed_name: str = "example_bed12_double.bed", side1: str = "target",
                 side2: str = "target"):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1, \
         open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb") as c2:
        return client.post(
            "/api/bed12/liftover2",
            files={
                "bed": (bed_name, bf, "text/plain"),
                "chain1": ("demo_seg1.chain", c1, "text/plain"),
                "chain2": ("demo_seg2.chain", c2, "text/plain"),
            },
            data={"side1": side1, "side2": side2},
        )


def test_bed12_single_preview_counts_and_kind():
    resp = _post_single()
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["kind"] == "bed12_single"
    preview = data["preview"]
    assert preview["total"] == 4
    assert preview["counts"] == {
        "mapped": 2, "partial": 0, "ambiguous": 1, "unmapped": 1,
    }
    # 3 条 mapped 行 + 2 条歧义候选 = 5 行 BED12；unmapped 行 0 行
    assert preview["exported_bed12_lines"] == 4


def test_bed12_download_is_valid_bed12_without_unmapped():
    token = _post_single().json()["download_token"]
    resp = client.get(f"/api/bed12/download/{token}")
    assert resp.status_code == 200
    assert resp.headers["content-type"].startswith("text/x-bed")
    body = resp.text
    lines = body.splitlines()
    assert len(lines) == 4
    for ln in lines:
        assert len(ln.split("\t")) == 12
    # 两个歧义候选都导出
    assert sorted(ln.split("\t")[0] for ln in lines if "tx-amb" in ln) == ["chrA", "chrB"]
    # 无映射的 tx-gap 不产生任何伪 BED12 行
    assert all("tx-gap" not in ln for ln in lines)


def test_bed12_bad_file_blocks_entire_conversion():
    bad = (
        "chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t2\t50,50,\t0,50,\n"
        "chr1\t0\t100\ttx\t0\t+\t50\t60\t0\t2\t10,10,\t0,90,\n"  # thick 不在外显子内
    )
    resp = client.post(
        "/api/bed12/liftover",
        files={
            "bed": ("bad.bed", bad.encode(), "text/plain"),
            "chain": ("demo.chain", open(os.path.join(EXAMPLES, "demo.chain"), "rb"),
                      "text/plain"),
        },
        data={"side": "target"},
    )
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("第 2 行" in e and "thickStart" in e for e in errors)


def test_bed12_chain_error_prefixes_double():
    resp = client.post(
        "/api/bed12/liftover2",
        files={
            "bed": ("b.bed",
                    b"chr1\t0\t100\ttx\t0\t+\t0\t0\t0\t1\t100,\t0,\n",
                    "text/plain"),
            "chain1": ("c1.chain", b"chain 1 chr1 100 + 0 10\n10\n", "text/plain"),
            "chain2": ("c2.chain", b"not a chain at all\n", "text/plain"),
        },
        data={"side1": "target", "side2": "target"},
    )
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any(e.startswith("chain 文件 1") for e in errors)
    assert any(e.startswith("chain 文件 2") for e in errors)


def test_bed12_bad_side_rejected():
    resp = client.post(
        "/api/bed12/liftover",
        files={
            "bed": ("b.bed", b"chr1\t0\t10\tx\t0\t+\t0\t0\t0\t1\t10,\t0,\n",
                    "text/plain"),
            "chain": ("c.chain", b"", "text/plain"),
        },
        data={"side": "upstream"},
    )
    assert resp.status_code == 400


def test_bed12_decisions_flow_and_adjudicated_download():
    token = _post_single().json()["download_token"]

    listing = client.get(f"/api/bed12/result/{token}/decisions").json()
    assert listing["version"] == 0
    assert listing["decisions"] == [
        {"source_line": 18, "candidate_index": None, "candidate_count": 2}
    ]

    resp = client.post(
        f"/api/bed12/result/{token}/decisions",
        json={"version": 0, "source_line": 18, "candidate_index": 2},
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["version"] == 1

    # 旧版本提交 -> 409
    stale = client.post(
        f"/api/bed12/result/{token}/decisions",
        json={"version": 0, "source_line": 18, "candidate_index": 1},
    )
    assert stale.status_code == 409
    assert stale.json()["current_version"] == 1

    # 非歧义行禁止裁决
    bad_row = client.post(
        f"/api/bed12/result/{token}/decisions",
        json={"version": 1, "source_line": 17, "candidate_index": 1},
    )
    assert bad_row.status_code == 400

    # 候选越界
    out_of_range = client.post(
        f"/api/bed12/result/{token}/decisions",
        json={"version": 1, "source_line": 18, "candidate_index": 9},
    )
    assert out_of_range.status_code == 400

    # 裁决版只保留被选候选（chrB）
    adj = client.get(f"/api/bed12/download/{token}/adjudicated")
    assert adj.status_code == 200
    lines = adj.text.splitlines()
    assert len(lines) == 3   # 2 mapped + 1 picked ambiguous
    assert not any("tx-amb" in ln and ln.startswith("chrA") for ln in lines)
    assert any("tx-amb" in ln and ln.startswith("chrB") for ln in lines)
    assert ".adjudicated.v1.bed" in adj.headers["content-disposition"]

    # 原始下载字节不随裁决改变
    raw = client.get(f"/api/bed12/download/{token}")
    assert len(raw.text.splitlines()) == 4


def test_bed12_token_separate_from_bed3_namespace():
    # BED12 token 在 BED3 下载接口上 404，反之互不影响
    token = _post_single().json()["download_token"]
    assert client.get(f"/api/download/{token}").status_code == 404
    assert client.get(f"/api/bed12/download/no-such-token").status_code == 404


def test_bed12_double_preview_and_failed_stage():
    data = _post_double().json()
    assert data["kind"] == "bed12_double"
    rows = data["preview"]["rows"]
    statuses = {r["name"]: (r["status"], r.get("failed_stage")) for r in rows}
    assert statuses["g1"] == ("mapped", "")
    assert statuses["g2"] == ("ambiguous", "")
    assert statuses["g3"] == ("unmapped", "segment1")
    assert statuses["g4"] == ("unmapped", "segment2")
    # 双段歧义候选带 chain1/chain2 id
    g2 = next(r for r in rows if r["name"] == "g2")
    assert len(g2["candidates"]) == 2
    assert {(c["chain1_id"], c["chain2_id"]) for c in g2["candidates"]} == {
        ("2", "11"), ("2", "12")
    }


def test_bed12_double_filter_and_audit_endpoints_are_bed3_only():
    token = _post_double().json()["download_token"]
    # 筛选/审计仅适用于 BED3 双段结果；BED12 token 在这些接口上 404
    assert client.post(f"/api/result/{token}/filter",
                       json={"min_coverage": 0.5,
                             "keep_statuses": ["mapped"]}).status_code == 404
    assert client.post(f"/api/result/{token}/audit",
                       json={}).status_code == 404
