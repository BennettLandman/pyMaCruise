"""Historical ANTs N4 automatic Otsu estimation mask, without an ITK runtime.

ANTs 2.1.0 ``Examples/N4BiasFieldCorrection.cxx`` creates an Otsu mask
with 200 histogram bins when no mask argument is supplied. ITK 4.9's
ImageToHistogramFilter extends the floating-point maximum by one hundredth
of a bin, and its Otsu calculator returns the selected bin midpoint. Modern
ANTs instead uses the entire image. This ports that historical mask boundary;
N4 bias estimation remains delegated to ANTsPy.

Source: https://github.com/ANTsX/ANTs/blob/v2.1.0/Examples/N4BiasFieldCorrection.cxx
ITK source copyright Insight Software Consortium, Apache License 2.0.
"""

from __future__ import annotations

import numpy as np


def legacy_n4_otsu_mask(values: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the historical 200-bin automatic mask and float32 threshold."""
    image = np.asarray(values, dtype=np.float32)
    if image.ndim != 3 or not image.size or not np.isfinite(image).all():
        raise ValueError("Legacy N4 Otsu input must be a nonempty finite 3D array")
    minimum, maximum = float(image.min()), float(image.max())
    if maximum <= minimum:
        raise ValueError("Legacy N4 Otsu estimation needs a nonconstant image")
    bins = 200
    upper = maximum + (maximum - minimum) / bins / 100.0
    # ITK Histogram::Initialize stores a float interval even when its
    # measurement type is double. Preserve the float multiply at each edge.
    interval = np.float32(np.float32(np.float32(upper) - np.float32(minimum)) / np.float32(bins))
    edges = minimum + (np.arange(bins + 1, dtype=np.float32) * interval).astype(np.float64)
    edges[-1] = upper
    counts, _ = np.histogram(image, bins=edges)
    centres = (edges[:-1] + edges[1:]) / 2.0
    frequencies = counts.astype(np.float64)
    total = float(image.size)
    global_sum = float(sum(float(c) * float(x) for c, x in zip(frequencies, centres, strict=True)))
    lower_count = 0.0
    lower_sum = 0.0
    best_score = -np.inf
    best_index = 0
    for index in range(bins - 1):
        lower_count += frequencies[index]
        lower_sum += frequencies[index] * centres[index]
        upper_count = total - lower_count
        if lower_count <= 0 or upper_count <= 0:
            continue
        lower_mean = lower_sum / lower_count
        upper_mean = (global_sum - lower_sum) / upper_count
        score = (lower_count * lower_mean * lower_mean + upper_count * upper_mean * upper_mean) / total
        # ITK 4.9 updates a winning threshold only beyond a one-ULP tie.
        if score > best_score and score > np.nextafter(best_score, np.inf):
            best_score, best_index = float(score), index
    threshold = float(np.float32(centres[best_index]))
    mask = image > threshold
    if not np.any(mask):
        raise ValueError("Legacy N4 Otsu estimation mask is empty")
    return mask, threshold
