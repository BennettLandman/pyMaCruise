"""Source-pinned reachability/order evidence for rigid atlas resampling."""

import json
import hashlib
from pathlib import Path


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_transformed_template"
EXPECTED_WRAPPER_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
EXPECTED_CALL_ORDER = [
    "if (alignShapes)",
    "preprocess.initialMemberships",
    "atlas.setImageInfo before registration",
    "atlas.initShapeRegistration rigid",
    "atlas.setTransform preprocessing transform",
    "atlas.updateShapeRegistration rigid",
    "atlas.registerShapesPyramid",
    "preprocess.uncropImages",
    "preprocess.findCroppingBoundaries",
    "preprocess.cropImages",
    "atlas.setImageInfo after recrop",
    "atlas.computeTransformedTemplate",
    "preprocess.updateTransformedTemplate(atlas.getTemplate())",
]


def test_rigid_transformed_template_is_reached_after_recrop_in_selected_wrapper():
    reference = json.loads((FIXTURE / "wrapper_call_order.json").read_text())
    assert reference["wrapper_source"] == "AlgorithmLesionToads.java"
    assert reference["wrapper_source_sha256"] == EXPECTED_WRAPPER_SHA256
    default_excerpt = (FIXTURE / "align_shapes_default_source_line.txt").read_bytes()
    order_excerpt = (FIXTURE / "wrapper_call_order_source_excerpt.java.txt").read_bytes()
    assert hashlib.sha256(default_excerpt).hexdigest() == reference[
        "align_shapes_default_excerpt_sha256"
    ]
    assert hashlib.sha256(order_excerpt).hexdigest() == reference[
        "wrapper_call_order_excerpt_sha256"
    ]
    assert reference["align_shapes_default"] is True
    assert reference["wrapper_call_order"] == EXPECTED_CALL_ORDER
    order = reference["wrapper_call_order"]
    excerpt_text = order_excerpt.decode("utf-8")
    for call in (
        "if (alignShapes)",
        "preprocess.initialMemberships",
        "atlas.registerShapesPyramid",
        "preprocess.uncropImages",
        "preprocess.findCroppingBoundaries",
        "preprocess.cropImages",
        "atlas.computeTransformedTemplate",
        "preprocess.updateTransformedTemplate(atlas.getTemplate())",
    ):
        assert call in excerpt_text
    assert b"alignShapes = true" in default_excerpt
    source_positions = [
        excerpt_text.index("if (alignShapes)"),
        excerpt_text.index("atlas.registerShapesPyramid();"),
        excerpt_text.index("preprocess.uncropImages();"),
        excerpt_text.index("preprocess.findCroppingBoundaries();"),
        excerpt_text.index("preprocess.cropImages();"),
        excerpt_text.rindex("atlas.setImageInfo("),
        excerpt_text.index("atlas.computeTransformedTemplate();"),
        excerpt_text.index("preprocess.updateTransformedTemplate(atlas.getTemplate());"),
    ]
    assert source_positions == sorted(source_positions)
    assert order.index("atlas.registerShapesPyramid") < order.index("preprocess.uncropImages")
    assert order.index("preprocess.cropImages") < order.index("atlas.setImageInfo after recrop")
    assert order.index("atlas.setImageInfo after recrop") < order.index("atlas.computeTransformedTemplate")
    assert order.index("atlas.computeTransformedTemplate") < order.index(
        "preprocess.updateTransformedTemplate(atlas.getTemplate())"
    )
