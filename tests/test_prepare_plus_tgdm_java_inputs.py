from scripts.prepare_plus_tgdm_java_inputs import INPUT_ORDER


def test_plus_tgdm_replay_preserves_workflow_boundary_and_tissue_input_order():
    assert INPUT_ORDER == (
        "segmentation",
        "toads_wm_mask",
        "multi_atlas_gm",
        "multi_atlas_wm",
        "multi_atlas_mask",
        "toads_gm",
        "toads_wm_fill",
        "toads_wm_mask",
    )
