# Bundled Atlas Assets

These atlas files are stored in this repository so PyMaCRUISE does not need the adjacent legacy source checkouts at runtime. The source tree and observed provenance for each dataset are recorded in [../../docs/assets-and-dependencies.md](../../docs/assets-and-dependencies.md).

Included datasets:

- `cruise-atlas-10obj/` and `cruise-atlas-10obj-toads2008.txt`: TOADS/CRUISE atlas used by the V3 base workflow. The root descriptor has been updated to resolve image payloads from this repository as `.nii.gz` files.
- `spectre/`: SPECTRE atlas assets used by the selected V3 skull-strip wrapper.
- `lesionToads-atlas-2012/`: lesion-TOADS assets used by the selected MaCRUISEplus workflow.

All assets are committed as regular Git files. Image payloads are compressed NIfTI (`.nii.gz`) with RAS voxel axes. `SHA256SUMS` lists SHA-256 hashes for each bundled file, excluding the checksum file itself. CVS working-copy metadata was omitted. Non-image MIPAV XML metadata and text descriptors remain in their native formats; descriptors point to the NIfTI payloads.

The user-provided NIfTI conversions were moved into their corresponding atlas folders and normalized by axis permutation/flip only. `CONVERSION_AUDIT.json` records source payload hashes, input/output geometry, exact voxel-value checks, output hashes, and metadata caveats. All 35 converted volumes match source voxel values exactly after RAS reorientation. The source Analyze examples do not establish a reliable world origin; some MIPAV XML entries report unknown orientation, and three rotated OAS1 raw files have an endian declaration that disagrees with their exact matching values. Resolve these source-coordinate details before claiming full atlas physical-space equivalence.

The selected base TOADS topology map has an additional known conflict: its text descriptor declares 0.9375 mm, while both source MIPAV XML and the converted NIfTI affine declare 1.0 mm. The Java atlas loader reads topology dimensions/resolution from the descriptor but reads voxel bytes without applying image-header geometry. `src/pymacruise/toads_atlas.py` preserves both spacings and reports a warning; it does not resample. This must be resolved or accepted explicitly as part of P18 before claiming physical atlas equivalence.

The MNI 305 reference image and mask are required for full MNI-space output, but were not present in the supplied source trees; they are not fabricated or substituted here. The associated output remains an explicit asset gap.

The project owner directs that required supplied assets remain bundled in this repository. Their provenance is recorded for the later tech-transfer review; this does not affect implementation progress.
