# Registration backend: NiftyReg to ANTsPy

This guide documents the selected legacy registration contract and the
ANTsPy operations currently used at the Python boundary. The papers describe
the methods/frameworks in context; the selected MATLAB/VandyCRUISE call sites
define the legacy behavior, and the ANTsPy API documentation defines the
Python library interface. Neither citations nor synthetic API tests establish
that ANTsPy reproduces NiftyReg Aladin numerically.

## References

1. Modat M, Cash DM, Daga P, Winston GP, Duncan JS, Ourselin S. Global image
   registration using a symmetric block-matching approach. *Journal of
   Medical Imaging*. 2014;1(2):024003.
   [doi:10.1117/1.JMI.1.2.024003](https://doi.org/10.1117/1.JMI.1.2.024003).
   [UCL author manuscript and record](https://discovery.ucl.ac.uk/id/eprint/1463769/).
2. Avants BB, Tustison NJ, Song G, Cook PA, Klein A, Gee JC. A reproducible
   evaluation of ANTs similarity metric performance in brain image
   registration. *NeuroImage*. 2011;54(3):2033–2044.
   [doi:10.1016/j.neuroimage.2010.09.025](https://doi.org/10.1016/j.neuroimage.2010.09.025).
   [PubMed record](https://pubmed.ncbi.nlm.nih.gov/20851191/).
3. Avants BB, Tustison NJ, Stauffer M, Song G, Wu B, Gee JC. The Insight
   ToolKit image registration framework. *Frontiers in Neuroinformatics*.
   2014;8:44.
   [doi:10.3389/fninf.2014.00044](https://doi.org/10.3389/fninf.2014.00044).
   [Full text](https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2014.00044/full).
4. Avants BB, Epstein CL, Grossman M, Gee JC. Symmetric diffeomorphic image
   registration with cross-correlation: evaluating automated labeling of
   elderly and neurodegenerative brain. *Medical Image Analysis*.
   2008;12(1):26–41.
   [doi:10.1016/j.media.2007.06.004](https://doi.org/10.1016/j.media.2007.06.004).
   [PubMed record](https://pubmed.ncbi.nlm.nih.gov/17659998/).
5. ANTsPy developers. `ants.registration` API reference.
   [Official API documentation](https://antspy.readthedocs.io/en/latest/api/ants.registration.html).
6. ANTsPy developers. Registration package reference, including
   `ants.apply_transforms`.
   [Official documentation](https://antspyx.readthedocs.io/en/stable/registration.html).
7. NiftyReg developers. `reg_resample` source and command-line defaults.
   [Official repository source](https://github.com/SuperElastix/niftyreg/blob/master/reg-apps/reg_resample.cpp).

## Paper-to-code map

| Code location | What the Python code does | What the reference explains | Scope boundary |
| --- | --- | --- | --- |
| [`registration.py`](../../src/pymacruise/registration.py): `RegistrationConfig`, `register_to_reference`, `register_flair_to_t1` | Exposes ANTsPy rigid/affine registration; maps legacy reference/floating inputs to ANTs fixed/moving; passes a fixed-space mask to all stages when present; saves transform files under a caller-supplied prefix. The selected defaults are Mattes metric, 32 bins, 0.2 random sampling, iterations `(1000, 500, 250, 0)`, shrink factors `(6, 4, 2, 1)`, smoothing sigmas `(3, 2, 1, 0)`, deterministic seed 123, and single precision. | The ANTs metric evaluation paper provides context for metric behavior and evaluation. The ITK framework paper describes registration as the composition of transform, metric, and optimizer choices and covers affine mappings and mutual-information metrics. The official ANTsPy API documents the actual arguments and their meanings. | These settings are explicit PyMaCRUISE choices/API settings; they are not prescribed by those papers and are not asserted to reproduce Aladin's block matching, objective, optimizer, initialization, or convergence. NiftyReg's symmetric block-matching paper describes the replaced backend, not the ANTsPy implementation. |
| [`ants_utils.py`](../../src/pymacruise/ants_utils.py): `nifti_to_ants` | Preserves voxel ordering and spatial geometry while converting values to `float32` by default, then represents origin and direction in ANTs/ITK LPS physical coordinates. | The ITK framework paper provides general context for physical-space transformations; it does not specify this RAS/LPS conversion. | The RAS↔LPS conversion and default dtype conversion are project adapters, not methods claimed from the papers. |
| [`registration.py`](../../src/pymacruise/registration.py): `apply_image_transform` | Calls ANTsPy `apply_transforms` on the fixed grid. Uses `nearestNeighbor` for integer labels and `bSpline` for continuous intensity and level-set images; then writes the result using the fixed RAS affine/header. | The official ANTsPy API documents fixed-grid resampling, transform-list ordering, and available interpolation names. NiftyReg's `reg_resample` source specifies nearest-neighbor option 0 and cubic spline as its default (option 3). | The wrapper's label interpolation intent is aligned, while the spline operation is a backend substitution. Similar interpolation names do not prove identical spline kernels, boundary behavior, precision, or voxel values. Outputs are constructed on the fixed RAS grid by project code. |
| [`registration.py`](../../src/pymacruise/registration.py): `apply_transforms_to_points_ras`, `apply_transforms_to_surface_ras` | Converts RAS points to LPS, applies the inverse of the image-resampling affine chain in reversed stack order, converts back to RAS, and reuses transformed point coordinates for surface vertices. The current point path accepts affine `.mat` transforms only. | ANTsPy documentation explains that transform lists are ordered and map a moving image into the fixed domain through resampling; the ITK framework paper describes composite mappings and inverse/forward concepts. | The point/surface direction handling and LPS/RAS conversion are implementation-specific API adapters. Nonlinear mesh warping is not supported here. This reference mapping does not establish a surface-registration equivalence result. |

## Deliberately excluded paper claim

Avants et al. (2008) describes symmetric diffeomorphic (SyN) registration with
cross-correlation. PyMaCRUISE's `RegistrationConfig` restricts this adapter to
`Rigid` or `Affine`; the SyN paper therefore does not explain the implemented
optimizer. It is listed for method-family context only, not as the method
implemented by these functions.

## Validation boundary

The selected wrappers use `reg_aladin` with reference/floating image roles;
the MNI path is rigid-only, and wrapper label resampling requests nearest
neighbor. The ANTsPy adapter is a documented software substitution. Its
synthetic known-transform tests cover coordinate direction, grid handling,
interpolation selection, and API behavior; they do not compare ANTsPy with
Aladin on shared inputs and do not establish numerical equivalence. See
[`registration-port.md`](../registration-port.md) for the current adapter
contract and its MNI asset limitation, and the P33/P38 rows in
[`plan.md`](../../plan.md) for their stated validation scope.
