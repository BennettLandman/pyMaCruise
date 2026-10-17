# Selected AtlasEM registration defaults — 2026-10-16

The selected MIPAV plugin constructs `AlgorithmAtlasEMWrapper` with
`register=true`, `levels=4`, `initRegIter=50`, and `mainRegIter=1`; it passes
demons smoothing `1.0` and demons scale `2.0`. The selected MATLAB
`run_EM_jist.m` explicitly supplies smoothing `10`, maximum membership change
`0.01`, maximum EM iterations `99`, Atlas prior mode, rigid registration, and
inhomogeneity correction disabled. Its run command does not override the
registration settings, so the plugin constructor values apply.

Source files and SHA-256:

- `VandyCRUISE/src/edu/jhu/ece/iacl/plugins/classification/MedicAlgorithmAtlasEMSegmentation.java`:
  `96ad55161d772995135604792054c2f00d159b01159f82159a9fd06a33c29c10`
- `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmAtlasEMSegmentation.java`:
  `923b3decdb4826351e3285eaa437648763d93af75e34e61fd1565eda744e54d2`
- `MaCRUISE/pipeline/run_EM_jist.m`:
  `badc38dc6a867cfa403f459e3ce65ee07f6243820a741583a8af508b7cd2e446`

PyMaCRUISE now exports these registration defaults from
`pymacruise.atlas_em_registration`; the descriptor-driven two-stage function
uses 4 levels, 50 initial registration iterations, and 1 post-recrop iteration
when callers do not override them. Tests retain one-iteration caps for compact
synthetic fixtures. The standalone plugin-default assertions and descriptor
registration tests pass **4 tests in 0.56 s** on test commit
`61d884c7f65b7570fb02fa1975e3a442d931fae4` (defaults implemented in
`fec74662569f059e92d75a95acdf9ca4661c0e7e`).

This locks the selected call settings, not the optimizer output. The selected
`ibsr8obj` descriptor and image payload remain unavailable; no human AtlasEM
run or direct Java optimizer comparison is claimed.
