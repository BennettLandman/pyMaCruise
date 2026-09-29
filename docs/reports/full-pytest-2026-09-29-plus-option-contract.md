# Full test suite: MaCRUISEplus option contract

- Date: 2026-09-29
- Tested implementation commit: `926d4672ad2e7aa79d9869db873e43b110bed7f5`
- Command: `.venv/bin/pytest -q`
- Result: **465 passed in 56.34 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The selected MaCRUISEplus entry-point overrides are represented by
`PlusPipelineOptions` and included in plus-profile manifest configuration.
The source-selected lesion output/correction, WM isolevel, GVF lambda,
enhancement, intermediate-output, CSF stage, boundary, and topology-refinement
values match the documented MATLAB settings. Invalid settings are rejected.
The base profile does not gain the plus-options field in its serialized
configuration. This is option-routing evidence only; it does not validate
LesionTOADS, ACE, topology correction, or surface outputs.

See the source-to-stage mapping in `docs/plus-pipeline-contract.md` and the
focused checks in `tests/test_plus_pipeline.py` and
`tests/test_workflow.py`.
