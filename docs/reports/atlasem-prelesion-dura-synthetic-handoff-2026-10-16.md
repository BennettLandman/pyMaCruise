# AtlasEM → pre-LesionTOADS Dura synthetic handoff — 2026-10-16

`run_atlas_em_prelesion_dura` connects the selected MATLAB handoff represented by `run_EM.m` and `run_RemoveDura.m`:

1. Accept source-grid original T1, class-last AtlasEM membership NIfTI, and its matching brain mask.
2. Reorient all images to RAS without interpolation and reject grid mismatches.
3. Derive the AtlasEM hard tissue-index map, then generate its GM membership, WM membership, and WM mask channels with the source `generateDuraRemovalOutputs` behavior.
4. Pass the three products with the original T1 to the existing selected Dura implementation using source wrapper settings: target/membership level 0.3, WM/GM topology regularization 0.05/0.02, curvature 0.5, alignment force/alpha 1.0/0.5, pressure 0.5, 5 passes by default, mask translation 2.0, and level-set translation 0.0.

The deterministic synthetic integration test builds a 17³ RAS spherical tissue case, runs the handoff and one Dura pass, and checks finite products, class/index maps, completed pass count, and RAS output affines. Focused command/result including the AtlasEM source-slice group: **19 passed in 5.47 s**.

This validates Python component composition on synthetic inputs. It does not validate the selected `ibsr8obj` resource bundle, atlas registration, complete source AtlasEM run, or real human PLUS equivalence. The actual source-provided atlas descriptor and payload remain unavailable, so the full PLUS route cannot yet use this adapter to generate source-faithful AtlasEM inputs.
