from __future__ import annotations

from pymacruise.base_workflow import BaseRunOptions


def test_base_tgdm_iteration_options_keep_legacy_fallback_and_allow_independent_counts() -> None:
    legacy = BaseRunOptions(tgdm_iterations=3)
    assert legacy.effective_tgdm_iterations() == (3, 3, 3)

    selected = BaseRunOptions(
        tgdm_iterations=2,
        tgdm_initial_iterations=6,
        tgdm_inner_iterations=4,
        tgdm_central_iterations=7,
    )
    assert selected.effective_tgdm_iterations() == (6, 4, 7)

    partial = BaseRunOptions(tgdm_iterations=2, tgdm_inner_iterations=5)
    assert partial.effective_tgdm_iterations() == (2, 5, 2)
