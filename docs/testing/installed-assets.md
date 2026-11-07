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
source file. It reports asset availability, not workflow equivalence. Its `ibsr11objV3_runtime_status` field is copied verbatim from the preserved source-package manifest and can retain the historical "conversion remains open" label; current conversion support is described in the PLUS walkthrough below. That historical metadata is not a current implementation audit.

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

The native source package preserves original indexing and metadata. Run `python scripts/prepare_plus_testing.py --input-dir source-reference-run/inputs/100260 --output-dir build/plus-e2e` from a checkout to unpack verified native bytes, prepare a Java-compatible descriptor and convert the same atlas to a Python RAS descriptor plus NIfTIs. Conversion has round-trip/grid tests; complete original-runtime AtlasEM agreement remains unverified. See [the PLUS walkthrough](full-plus-walkthrough.md) for assumptions and both pipeline commands.

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


## Locate the bundled atlas for BASE

For a T1-only run, `pymacruise base-run` is the correct command. After wheel installation, atlas descriptors and NIfTIs live inside the installed Python package, not the extracted bundle's root `data/atlases/` folder (which holds provenance documents). Resolve the installed descriptor in the same activated environment:

```bash
python -m pymacruise.asset_inventory
ATLAS_DESCRIPTOR="$(python -c 'from pymacruise.atlas_resources import get_bundled_atlas_descriptor; print(get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt"))')"
printf '%s\n' "$ATLAS_DESCRIPTOR"
```

Use `--atlas-descriptor "$ATLAS_DESCRIPTOR"` in the existing `pymacruise base-run` command, retaining the other chosen inputs/settings. No atlas download or manual extraction is needed. BASE uses T1, MaSeg and a corresponding brain mask; for T1 plus FLAIR use the PLUS workflow above.

From version 0.1.1, a bare `base-run` uses full selected-source numerical
settings (TOADS 50, ACE diffusion 25, TGDM 6/4/7 and full registration/Dura
settings). Final outputs return to the original T1 grid and physical frame.
Older 0.1.0 wheels use short smoke-test limits; upgrade before a new human
comparison. See the [BASE comparison recipe](python-base-checks.md). Start a
new output directory after changing settings or upgrading; check `converged`
as well as `status` in `run.json`. These defaults support equivalence testing;
they do not establish full workflow agreement.
