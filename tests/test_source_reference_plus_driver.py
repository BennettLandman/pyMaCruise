from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "source-reference-run/matlab/run_plus_selected_mask_reference.m"


def test_plus_reference_driver_preserves_selected_stage_order_and_mask_contract():
    text = DRIVER.read_text(encoding="utf-8")
    selected_calls = (
        "FlairRegistration(",
        "prepare_raw_plus(",
        "cortical_seperation(",
        "run_N4Correct(T1_file",
        "run_N4Correct(Flair_file",
        "apply_mask(python_executable, mask_helper, t1_n4",
        "apply_mask(python_executable, mask_helper, flair_n4",
        "run_EM(",
        "run_RemoveDura(",
        "apply_mask(python_executable, mask_helper, flair_masked, dura_mask",
        "moveFile(t1_dura, opts, 'removeDura_final'",
        "moveFile(flair_dura_stripped, opts, 'removeDura_final'",
        "run_LesionToads(",
        "run_TransferToOriginSpace(",
        "run_CropBaseOnGM(",
        "run_GetMABoundary(",
        "update_WMandGM(",
        "run_CleanBoundaryByUpdateWM(",
        "run_CruisePlusSurf_pipeline(",
        "run_SegmentationRefine(",
    )
    positions = [text.index(call) for call in selected_calls]
    assert positions == sorted(positions)
    assert "opts.rigOnly = 1;" in text
    assert "opts.LesionToadsOutput = 'dura removal inputs';" in text
    assert "FlairRegistration" in text and "aladin_reg" not in text
    assert "run_Bet(" not in text
    assert "run_Spectre(" not in text
    assert "write_FreeViewBatch(" not in text
    assert "run_MakeLevelsetAndSurfaceMNI(" not in text
    assert "opts.freesurf_dir" not in text
    assert "transferVTKtoASC_cruise(" not in text

