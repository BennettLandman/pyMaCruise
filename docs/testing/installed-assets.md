# Asset availability in Git and installed packages

All supplied atlas assets live under `data/atlases/` as regular Git files.
`MANIFEST.in` includes the entire tree in source distributions. The build hook
in `setup.py` stages the same tree under `pymacruise/data/atlases/` in wheels.
No runtime download, Git LFS pull, or adjacent MIPAV/VandyCRUISE checkout is
required to access the bundled data. The critical-point LUT and JIST mesh
case table live in `src/pymacruise/data/` and are included as package data.

The available data comprise TOADS/CRUISE 2008, LesionTOADS 2012, the historical
SPECTRE assets, and the owner-selected `ibsr11objV3` source package. Atlas and
runtime-table `SHA256SUMS` manifests detect missing or altered files.

After installing a wheel or source distribution, run:

```sh
python -m pymacruise.asset_inventory
```

This validates every stored atlas checksum, all NIfTI headers and descriptor
references, both runtime lookup tables, and every decompressed `ibsr11objV3`
source file. It reports asset availability, not workflow equivalence.

Resolve installed atlas resources through the package API:

```python
from pymacruise.atlas_resources import (
    get_bundled_atlas_descriptor,
    get_bundled_atlas_asset,
    read_bundled_atlas_source,
)

base = get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt")
native_inventory = get_bundled_atlas_asset("ibsr11objV3/SOURCE_MANIFEST.json")
descriptor_bytes = read_bundled_atlas_source(
    "ibsr11objV3", "LesionToads-11-obj-atlasV3.txt"
)
shape_bytes = read_bundled_atlas_source("ibsr11objV3", "atlas_CerebralGM.raw")
```

The `ibsr11objV3` source package includes all 29 supplied files: the descriptor,
14 raw payloads (including the alternate small topology), and 14 XML headers.
All original files are stored with deterministic gzip compression;
the accessor decompresses them and verifies source and stored hashes. The
CVS revision metadata is recorded in the source manifest; CVS working-copy
directories are omitted.

The native source package preserves original indexing and metadata. It is
available from a fresh installation, but it is not yet a converted AtlasEM
runtime descriptor: path/layout adaptation and checked RAS conversion remain
necessary. Packaging the files does not close that separate implementation gate.

The current Python workflows take T1, MaSeg, brain mask, and PLUS FLAIR as user
inputs. ANTsPy is a separately installed optional dependency, not a bundled
data asset. Historical MNI305 reference/mask files are not supplied and are
not needed by the implemented workflows; complete legacy MNI exports remain
unavailable until those exact assets are obtained.

For a distribution-level audit from a checkout, run:

```sh
python scripts/verify_atlas_distribution.py
```

It builds a source distribution, builds a wheel from that archive, verifies
the complete atlas file set and hashes in both archives, installs the wheel
in a separate target, and runs the installed asset check. Dependencies must
already be installed in the invoking Python environment. The receipt is
saved in the reported build directory as `verification.json`.
