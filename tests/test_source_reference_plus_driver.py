from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "source-reference-run/matlab/run_plus_selected_mask_reference.m"
PACKET = ROOT / "source-reference-run/README.md"
RUN_INFO = ROOT / "source-reference-run/RUN_INFO.template.md"


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
        "run_plus_atlasem_reference(",
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
        "run_plus_thickness_reference(",
    )
    positions = [text.index(call) for call in selected_calls]
    assert positions == sorted(positions)
    assert "opts.rigOnly = 1;" in text
    assert "opts.LesionToadsOutput = 'dura removal inputs';" in text
    assert "opts.atlasem_descriptor = config.atlasem_descriptor;" in text
    assert "'ibsr8obj'" not in text
    assert "'cruise-atlas-12obj-lesiontoads2012.txt'" in text
    assert "Missing selected PLUS atlas" in text
    assert "'atlasem_gm', EM_gm" in text
    assert "'atlasem_wmfill', EM_wmfill" in text
    assert "'atlasem_wmmask', EM_wmmask" in text
    assert "FlairRegistration" in text and "aladin_reg" not in text
    assert "run_Bet(" not in text
    assert "run_Spectre(" not in text
    assert "write_FreeViewBatch(" not in text
    assert "run_MakeLevelsetAndSurfaceMNI(" not in text
    assert "opts.freesurf_dir" not in text
    assert "transferVTKtoASC_cruise(" not in text


def test_plus_source_run_packet_records_atlasem_as_a_distinct_handoff():
    driver = DRIVER.read_text(encoding="utf-8")
    packet = PACKET.read_text(encoding="utf-8")
    run_info = RUN_INFO.read_text(encoding="utf-8")
    for output in ("atlasem_gm", "atlasem_wmfill", "atlasem_wmmask"):
        assert f"'{output}'" in driver
    assert "initial Dura removal" in packet
    assert "do not combine AtlasEM and Dura" in packet
    assert "PLUS targets `ibsr11objV3`" in packet
    assert "source-config.json" in packet
    assert "PLUS AtlasEM descriptor" in run_info
    assert "PLUS AtlasEM GM, WM-fill, and WM-mask" in run_info

