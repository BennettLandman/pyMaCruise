# Owner-selected MaCRUISE v3.2.0 reference — 2026-10-16

The owner selected the v3.2.0 container as the external reference. This
supersedes the initial v3.1.0 walkthrough. The [current complete BASE
walkthrough](../testing/full-base-walkthrough.md) retains venv setup, both
independent full scientific runs, comparison, and return-file instructions.

## Immutable image and inspected artifacts

The public Docker registry returned these properties for
`masidocker/harmonizeapp:macruise`:

- Manifest digest:
  `sha256:96c6c7a033edb86f2be7b3a60c3e8cc658691de3f670b35a4147feaf287a28a6`
- Platform: `linux/amd64`.
- Entrypoint: `/usr/local/bin/run_harmonize.sh`; default command: `--help`.
- Compressed layer total: approximately 12.450 GB.
- Embedded runtime: `/singularity/macruise_classern_v3.2.0.simg`.
- Complete scientific executable: `/extra/MaCRUISE_v3_2_0_classern`.

Only registry manifest/configuration metadata and the final three small
launch-script layers (2,452, 1,572 and 1,566 compressed bytes) were retrieved.
Each compressed layer's SHA-256 matched its registry digest. No full image or
embedded runtime was downloaded or executed.

Inspected script contents have these SHA-256 values:

| Script | SHA-256 |
| --- | --- |
| `/scripts/run_macruise.sh` | `cc0c7e867c080ae870138b89d23e0a90fbc2812879101f88a969d956f55282a9` |
| `/usr/local/bin/run_harmonize.sh` | `8e451c7aa6b835b8f492c3f1d059becdd633fe1c7da0ec5fdd3ee6d75676113c` |

## Actual wrapper behavior

The top-level entrypoint accepts a BIDS directory, output directory,
`participant`, and session/subject filters, then invokes the MaCRUISE wrapper.
That wrapper discovers T1 from a SLANT-TICV `pre/FinalResult` directory and
MaSeg from its `post/FinalResult` directory, stages them under `/INPUTS` as
`T1.nii.gz` and `orig_target_seg.nii.gz`, and launches the complete embedded
Singularity application. It retains `--nv`, binds `/INPUTS` and `/OUTPUTS`,
uses Xvfb screen 3, and passes the scientific executable path twice. The
repeated path is present in the actual image and is preserved in the recipe.

After launching the application, the wrapper copies outputs to the named
subject directory, removes its temporary inputs/outputs, and prints completion.
It does not test or preserve the scientific process's exit status before these
operations; successful copying/cleanup can conceal an application failure.

## Chosen test invocation and scope

The revised test invokes the very same full scientific command directly via
the outer Docker image's Singularity executable. It retains the embedded
runtime, scientific binary, arguments, Xvfb setup and runtime flags, with host
directories mounted as `/INPUTS` and `/OUTPUTS`. This captures the actual
scientific process's exit status and keeps all scientific outputs.

This bypasses outer BIDS discovery, staging, output copying and cleanup. It is
not an unchanged invocation of the BIDS convenience wrapper. It runs the
complete external scientific application from raw T1/MaSeg without replacing
scientific stages or feeding it Python intermediates. Native skull stripping
and other preprocessing remain external application behavior. Python still
uses its supplied-mask interface; that boundary difference remains explicit.

The published container usage requires privileged Docker for nested
Singularity, and the inspected application command includes `--nv`. These
runtime requirements and any GPU-library warnings need to be recorded by the
tester. The exact runtime/export layout, native mesh frame, complete scientific
execution, and end-to-end numeric agreement remain unverified on this host.
