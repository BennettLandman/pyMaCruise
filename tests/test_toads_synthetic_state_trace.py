from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts/trace_synthetic_toads_state.py"
_SPEC = importlib.util.spec_from_file_location("trace_synthetic_toads_state", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
_array_summary = _MODULE._array_summary
_float32_bits = _MODULE._float32_bits
build_trace = _MODULE.build_trace


def test_python_only_toads_state_trace_schema_and_deterministic_summaries() -> None:
    values = np.asarray([1.0, -0.0, np.float32(0.125)], dtype=np.float32)
    assert _array_summary(values) == _array_summary(values.copy())
    assert _float32_bits(np.float32(1.0)) == "0x3f800000"

    trace = build_trace(maximum_iterations=2)
    assert trace["evidence_kind"] == "python_only_synthetic_diagnostic"
    assert trace["source_equivalence_claim"] is False
    assert trace["real_data_validation"] is False
    assert trace["do_growing"] is True
    assert len(trace["iterations"]) == 2
    assert trace["padding"] == 2
    assert trace["registration_memberships"]["pre_registration"]["sha256"]
    assert trace["registration_memberships"]["post_registration"]["sha256"]
    first = trace["iterations"][0]
    assert first["membership"]["objective_float32_bits"].startswith("0x")
    assert first["membership"]["memberships"]["sha256"]
    assert first["registration"]["transform"]["sha256"]
    assert first["post_registration"]["growing"] is not None
    assert first["post_registration"]["relations"]["closest_labels"]["sha256"]
    assert first["post_registration"]["next_thinning"] is not None
    terminal = trace["iterations"][-1]
    assert terminal["registration"]["stop"] is True
    assert terminal["post_registration"]["growing"] is not None
    assert terminal["post_registration"]["next_thinning"] is None
    assert trace["terminal"]["iteration_limit_reached"] is True
