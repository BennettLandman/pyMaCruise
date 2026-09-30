from __future__ import annotations

import nibabel as nib

from pymacruise.atlas_resources import (
    AtlasResourceError,
    get_bundled_atlas_descriptor,
    get_bundled_atlas_root,
)
from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_topology,
)
from pymacruise.toads_atlas import load_toads_atlas
from pymacruise.toads_atlas_type_membership import read_toads_atlas_class_types


def test_source_checkout_resolver_loads_bundled_default_and_2012_descriptors() -> None:
    root = get_bundled_atlas_root()
    assert root.is_dir()
    assert root.name == "atlases"

    default = load_toads_atlas(
        get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt")
    )
    assert len(default.structures) == 11
    assert default.topology.path.suffixes[-2:] == [".nii", ".gz"]
    assert all(image.path.is_file() for image in default.shapes.values())

    toads2012_descriptor = get_bundled_atlas_descriptor(
        "lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt"
    )
    assert len(read_toads_atlas_class_types(toads2012_descriptor)) == 13
    toads2012_topology = load_lesion_atlas_topology(toads2012_descriptor)
    assert toads2012_topology.image.shape == (172, 216, 164)
    _assert_descriptor_assets_exist(toads2012_descriptor)


def test_source_checkout_resolver_loads_lesiontoads_numeric_and_topology_assets() -> None:
    descriptor = get_bundled_atlas_descriptor(
        "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
    )
    priors = load_lesion_atlas_priors(descriptor, ["T1_MPRAGE", "FLAIR"])
    assert len(priors.structures) == 13
    assert priors.modalities == ("T1_MPRAGE", "FLAIR")
    _assert_descriptor_assets_exist(descriptor)

    topology = load_lesion_atlas_topology(descriptor)
    assert topology.path.is_file()
    assert topology.path.suffixes[-2:] == [".nii", ".gz"]
    assert nib.aff2axcodes(topology.image.affine) == ("R", "A", "S")


def test_descriptor_resolver_rejects_paths_outside_bundle() -> None:
    for value in ("../outside.txt", "/tmp/outside.txt", "not-a-descriptor.json"):
        try:
            get_bundled_atlas_descriptor(value)
        except AtlasResourceError:
            continue
        raise AssertionError(f"Expected unsafe/non-descriptor path to fail: {value}")


def _assert_descriptor_assets_exist(descriptor) -> None:
    for line in descriptor.read_text(encoding="utf-8").splitlines():
        if line.startswith("File:"):
            relative = line.partition(":")[2].strip()
            assert relative.endswith(".nii.gz")
            assert (descriptor.parent / relative).is_file(), relative
