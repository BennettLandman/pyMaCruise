"""LesionTOADS post-constructor atlas-registration setup handoff.

The active wrapper calls ``atlas.updateShapeRegistration`` with the
segmentation object's membership buffer after centroid setup. This bounded
adapter captures the selected default rigid setup arguments and membership
array reference; it does not perform transform-type conversion, coefficient
calculation, or registration optimization.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class LesionRegistrationSetupHandoff:
    """Wrapper values handed to the later atlas registration setup call."""

    memberships: NDArray[np.float32]
    transform_type: str
    iterations: int
    levels: int
    subsample: int
    coupling_factor: np.float32
    call_order: tuple[str, ...]


def lesion_registration_setup_handoff(
    memberships: NDArray[np.float32],
    *,
    transform_type: str = "rigid",
    iterations: int = 1,
    levels: int = 1,
    coupling_factor: float = 1.0,
) -> LesionRegistrationSetupHandoff:
    """Capture the selected single-rigid call arguments without copying maps.

    The selected default wrapper uses ``registrationMode='rigid'`` and
    ``alignIter=1``. ``updateShapeRegistration`` hard-codes subsample=2. This
    handoff represents only those source-visible arguments and the membership
    reference; the atlas method's transform-model reconfiguration remains a
    separate registration implementation step.
    """
    values = np.asarray(memberships)
    if values.ndim != 4 or any(size < 1 for size in values.shape):
        raise ValueError("memberships must have shape (x, y, z, classes)")
    if values.dtype != np.float32 or not np.isfinite(values).all():
        raise ValueError("memberships must be finite float32 values")
    if transform_type != "rigid":
        raise ValueError("this selected wrapper handoff covers rigid mode only")
    if isinstance(iterations, (bool, np.bool_)) or not isinstance(iterations, (int, np.integer)) or iterations < 0:
        raise ValueError("iterations must be a nonnegative integer")
    if isinstance(levels, (bool, np.bool_)) or not isinstance(levels, (int, np.integer)) or levels < 1:
        raise ValueError("levels must be a positive integer")
    factor = np.float32(coupling_factor)
    if not np.isfinite(factor):
        raise ValueError("coupling_factor must be finite")
    return LesionRegistrationSetupHandoff(
        memberships=values,
        transform_type=transform_type,
        iterations=int(iterations),
        levels=int(levels),
        subsample=2,
        coupling_factor=factor,
        call_order=("segmentation.getMemberships", "atlas.updateShapeRegistration"),
    )
