# Bundled Atlas Assets

These atlas files are stored in this repository so PyMaCRUISE does not need the adjacent legacy source checkouts at runtime. The source tree and observed provenance for each dataset are recorded in [../../docs/assets-and-dependencies.md](../../docs/assets-and-dependencies.md).

Included datasets:

- `cruise-atlas-10obj/`: TOADS/CRUISE atlas used by the V3 base workflow.
- `spectre/`: SPECTRE atlas assets used by the selected V3 skull-strip wrapper.
- `lesionToads-atlas-2012/`: lesion-TOADS assets used by the selected MaCRUISEplus workflow.

All assets are committed as regular Git files. `SHA256SUMS` lists SHA-256 hashes for each bundled file, excluding the checksum file itself. CVS working-copy metadata was omitted.

The MNI 305 reference image and mask are required for full MNI-space output, but were not present in the supplied source trees; they are not fabricated or substituted here. The associated output remains an explicit asset gap.

The owner's instruction authorizes bundling these supplied assets in this private project. This does not establish permission for a public release; source/data terms should be resolved before changing repository visibility or redistributing a release.
