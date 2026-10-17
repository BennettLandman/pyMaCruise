# AtlasEM descriptor-driven initial memberships

## Scope

Commit `f6a7962e7831abc462259abe2483aad669772965` adds `initialize_atlas_em_memberships_from_descriptor`. Given normalized, cropped image channels and source-order modality names, it reads the descriptor's intensity profiles, applies modality aliases, builds class centroids using source robust bounds, and computes inverse-distance memberships with the existing source-compared initialization routine. The returned class names preserve descriptor order.

This composes previously validated operations; it does not load or normalize a subject image, select the crop, consume the brain mask, register atlas shape priors, or run later EM iterations. It has not been matched end-to-end against the selected Java class.

## Validation

The focused descriptor initialization and preprocessing tests passed **4 tests in 0.93 s**. The broader AtlasEM regression passed **36 tests, 3 expected orientation warnings, in 5.96 s** on tested commit `f6a7962e7831abc462259abe2483aad669772965`. A synthetic descriptor confirms profile aliasing and exact output equality to the underlying centroid and membership functions. No full project suite was run.
