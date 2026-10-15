from __future__ import annotations

import importlib.util
from pathlib import Path

_SCRIPT = Path(__file__).parents[1] / "scripts" / "compare_p41z_registered_thinning_java.py"
_SPEC = importlib.util.spec_from_file_location("compare_p41z_registered_thinning_java", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def test_terminal_checkpoint_comparison_reports_unstored_loop_arrays():
    available, unavailable = _MODULE._available_checkpoint_arrays(
        {"summary", "segmentation"},
        ("segmentation", "skeleton", "availability", "ordering", "class_volume"),
    )

    assert available == ["segmentation"]
    assert unavailable == ["skeleton", "availability", "ordering", "class_volume"]


def test_resumable_checkpoint_comparison_reports_all_requested_arrays():
    names = ("segmentation", "skeleton", "availability", "ordering", "class_volume")
    available, unavailable = _MODULE._available_checkpoint_arrays(
        set(names), names
    )

    assert available == list(names)
    assert unavailable == []


def test_nonnegative_thinning_override_is_applied_without_changing_checkpoint():
    settings = {"thinning_limit": 1000.0}
    selected = _MODULE._settings_with_thinning_override(settings, 0.0)

    assert selected["thinning_limit"] == 0.0
    assert settings["thinning_limit"] == 1000.0


def test_negative_thinning_override_is_rejected():
    import pytest

    with pytest.raises(ValueError, match="finite and nonnegative"):
        _MODULE._settings_with_thinning_override({"thinning_limit": 3.0}, -1.0)
