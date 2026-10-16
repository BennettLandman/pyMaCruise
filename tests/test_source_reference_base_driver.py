from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "source-reference-run/matlab/run_base_selected_mask_reference.m"


def test_base_reference_driver_preserves_selected_stage_order_and_mask_contract():
    text = DRIVER.read_text(encoding="utf-8")
    selected_calls = (
        "prepare_raw(",
        "cortical_seperation(",
        "run_N4Correct(T1_file",
        "system(mask_cmd)",
        "run_Toads( ...\n    masked_n4_file, opts, 'Toads1stPass'",
        "run_RemoveDura(",
        "run_Toads( ...\n    dura_t1, opts, 'Toads2ndPass'",
        "transfer_orig = run_TransferToOriginSpace(",
        "run_CropBaseOnGM(",
        "run_GetMABoundary(",
        "run_CleanBoundaryByToads(",
        "update_WMandGM(",
        "run_CleanBoundaryByUpdateWM(",
        "run_maCruiseSurf_pipeline(",
        "run_SegmentationRefine(",
    )
    positions = [text.index(call) for call in selected_calls]
    assert positions == sorted(positions)
    assert text.rfind("opts.working_dir = working_dir;", 0, positions[7]) > positions[6]
    assert "opts.N3orN4 = 'N4';" in text
    assert "opts.ToadsOutput = 'dura removal inputs';" in text
    assert "opts.ToadsOutput = 'cruise inputs';" in text
    assert "'cruise-atlas-10obj-toads2008.txt'" in text
    assert "Missing selected BASE TOADS atlas" in text
    assert "run_Bet(" not in text
    assert "run_Spectre(" not in text
    assert "write_FreeViewBatch(" not in text
    assert "run_MakeLevelsetAndSurfaceMNI(" not in text
    assert "opts.freesurf_dir" not in text

