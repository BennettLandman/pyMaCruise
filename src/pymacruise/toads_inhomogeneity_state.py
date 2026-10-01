"""Selected TOADS inhomogeneity-field attachment state.

This models only ``MultipleToadSegmentation.addInhomogeneityCorrection``:
the supplied per-channel field object is stored by reference and enables
``useField``. Field fitting and downstream corrected-image arithmetic are
separate operations.
"""

from __future__ import annotations

import operator
from typing import Any


class ToadsInhomogeneityFieldState:
    """Hold per-channel field references with selected Java setter semantics."""

    def __init__(self, channel_count: int) -> None:
        if isinstance(channel_count, bool) or not isinstance(channel_count, int) or channel_count < 1:
            raise ValueError("channel_count must be a positive integer")
        self.field_by_channel: list[Any | None] = [None] * channel_count
        self.use_field = False

    def add_inhomogeneity_correction(self, field: Any | None, channel: int) -> None:
        """Assign ``field`` to ``channel`` without copying, then enable useField.

        The Java method performs direct array assignment and sets the flag even
        when the supplied field is ``null``. This method preserves that order
        and behavior; invalid channel indices naturally raise ``IndexError``.
        """
        if isinstance(channel, bool):
            raise TypeError("channel must be an integer index")
        index = operator.index(channel)
        if index < 0 or index >= len(self.field_by_channel):
            raise IndexError("channel index is outside the configured channel range")
        self.field_by_channel[index] = field
        self.use_field = True
