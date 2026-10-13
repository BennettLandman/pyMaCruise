from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.toads_shapes import ToadShapeBounds


def test_numba_registered_shape_sampling_matches_python_float32_oracle() -> None:
    rng = np.random.default_rng(20261012)
    atlas = tuple(rng.random((13, 14, 15), dtype=np.float32) for _ in range(3))
    bounds = tuple(
        ToadShapeBounds((0, 0, 0), (9, 10, 11), 1, None) for _ in atlas
    )
    registration = SimpleNamespace(
        transform=np.asarray([0.03, -0.02, 0.01, 0.2, -0.1, 0.15], dtype=np.float32),
        shape_bounds=bounds,
    )
    kwargs = {
        "image_dimensions_xyz": (9, 10, 11),
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
    }
    python = transform_registered_lesion_shape_priors(
        registration, atlas, backend="python", **kwargs
    )
    numba = transform_registered_lesion_shape_priors(
        registration, atlas, backend="numba", **kwargs
    )

    assert np.array_equal(numba.view(np.uint32), python.view(np.uint32))


def test_registered_shape_sampling_rejects_unknown_backend() -> None:
    registration = SimpleNamespace(
        transform=np.zeros(6, dtype=np.float32),
        shape_bounds=(ToadShapeBounds((0, 0, 0), (3, 3, 3), 1, None),),
    )
    atlas = (np.ones((4, 4, 4), dtype=np.float32),)
    try:
        transform_registered_lesion_shape_priors(
            registration,
            atlas,
            image_dimensions_xyz=(3, 3, 3),
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            backend="unknown",  # type: ignore[arg-type]
        )
    except ValueError as error:
        assert "backend" in str(error)
    else:
        raise AssertionError("unknown registered-shape backend was accepted")
