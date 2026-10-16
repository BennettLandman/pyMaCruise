"""AtlasEM rigid-registration preparation backed by shared TOADS kernels.

AtlasEM's ``DemonToadDeformableAtlas.registerShapesPyramid`` builds its
membership and per-class shape pyramids with ``ImageFunctions.subsample``.
The array operation is shared with the already source-checked TOADS pyramid
builder; the optimizer and transform update remain separate AtlasEM work.
"""

from __future__ import annotations

from numpy.typing import ArrayLike

from pymacruise.toads_registration_sampling import (
    ToadsRegistrationPyramidLevel,
    build_toads_registration_pyramid,
)


def build_atlas_em_rigid_registration_pyramid(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    levels: int,
) -> tuple[ToadsRegistrationPyramidLevel, ...]:
    """Build AtlasEM's source-ordered coarse membership/shape levels.

    For ``levels=L``, the selected Java wrapper visits scales
    ``2**(L-1), ..., 2`` from coarse to fine and does not run a scale-1
    registration stage. Inputs are native-grid ``(x,y,z,class)`` float32
    arrays; each level is block-mean subsampled directly from those arrays.
    This only prepares registration inputs. It does not estimate a rigid
    transform or claim that the TOADS optimizer equals AtlasEM registration.
    """
    return build_toads_registration_pyramid(memberships, atlas_shapes, levels=levels)
