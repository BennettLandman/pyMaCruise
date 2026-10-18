from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from pymacruise.asset_inventory import verify_bundled_assets
from pymacruise.atlas_resources import (
    AtlasResourceError, get_bundled_atlas_asset, read_bundled_atlas_source,
)


def test_native_plus_assets_include_every_original_descriptor_payload():
    manifest = json.loads(get_bundled_atlas_asset("ibsr11objV3/SOURCE_MANIFEST.json").read_text())
    descriptor = read_bundled_atlas_source("ibsr11objV3", manifest["descriptor"]).decode()
    references = [line.partition(":")[2].strip() for line in descriptor.splitlines()
                  if line.startswith("File:")]
    assert references == manifest["referenced_payloads"]
    assert len(references) == 13
    assert "Structures:\t12" in descriptor
    for reference in references:
        name = reference.split("/")[-1]
        payload = read_bundled_atlas_source("ibsr11objV3", name)
        assert hashlib.sha256(payload).hexdigest() == manifest["files"][name]["source_sha256"]
        dtype = "i1" if "topology" in name else "<f4"
        values = np.frombuffer(payload, dtype=dtype)
        assert values.size == 172 * 216 * 164
        assert np.isfinite(values).all()
        if dtype == "<f4":
            assert values.min() >= 0 and values.max() <= 1


@pytest.mark.parametrize("atlas,name", [("../escape", "file.raw"), ("ibsr11objV3", "../escape"),
                                        ("ibsr11objV3", "missing.raw"), ("..", "file.raw")])
def test_source_accessor_rejects_unsafe_or_missing_assets(atlas, name):
    with pytest.raises(AtlasResourceError):
        read_bundled_atlas_source(atlas, name)


def test_complete_installed_asset_manifest_and_decoding():
    receipt = verify_bundled_assets()
    assert receipt["runtime_tables_checked"] == 2
    assert receipt["ibsr11objV3_referenced_payloads_checked"] == 13
    assert receipt["ras_nifti_files_checked"] == 35


def test_source_accessor_detects_corrupted_packaged_payload(tmp_path, monkeypatch):
    from pymacruise import atlas_resources
    atlas = tmp_path / "ibsr11objV3"
    atlas.mkdir()
    (atlas / "bad.gz").write_bytes(b"corrupt stored data")
    (atlas / "SOURCE_MANIFEST.json").write_text(json.dumps({"files": {"shape.raw": {
        "stored_path": "bad.gz", "stored_sha256": "0" * 64, "encoding": "gzip",
        "source_sha256": "0" * 64, "source_bytes": 4,
    }}}))
    monkeypatch.setattr(atlas_resources, "get_bundled_atlas_root", lambda: tmp_path)
    with pytest.raises(AtlasResourceError, match="Stored asset checksum mismatch"):
        read_bundled_atlas_source("ibsr11objV3", "shape.raw")
