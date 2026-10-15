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
