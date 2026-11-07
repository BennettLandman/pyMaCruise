"""Exercise kit preparation from bundled source data and its native layout."""
from __future__ import annotations
import importlib
import json
from pathlib import Path
import nibabel as nib
import numpy as np
import pytest


def test_kit_has_all_native_payloads_and_preserves_original_voxels(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / 'scripts'))
    prepare = importlib.import_module('prepare_plus_testing')
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    for name in ('T1', 'FLAIR', 'MaSeg', 'brain_mask'):
        image = nib.Nifti1Image(np.ones((3, 4, 5), dtype=np.uint8), np.eye(4))
        image.header.set_xyzt_units('mm')
        nib.save(image, inputs / f'{name}.nii.gz')
    output = tmp_path / 'kit'
    script = prepare.prepare(output, inputs)
    kit = json.loads((output / 'testing-kit.json').read_text())
    native = Path(kit['source_atlasem_descriptor'])
    parsed = prepare._source_assets(Path(kit['original_atlasem_descriptor']))
    rows = native.read_text().splitlines()
    # The selected Java parser only loops over consecutive Structure/File
    # pairs. Ensure every source class remains reachable in that grammar.
    first = next(i for i, line in enumerate(rows) if line.startswith('Structure:'))
    for i, (name, _, _) in enumerate(parsed['structures']):
        assert rows[first + 2*i] == f'Structure:\t{name}'
        assert rows[first + 2*i + 1].startswith('File:\t')
        path = native.parent / rows[first + 2*i + 1].split('\t')[1]
        assert prepare.hashlib.sha256(path.read_bytes()).hexdigest() == prepare.hashlib.sha256(parsed['shape_records'][name].read_bytes()).hexdigest()
    adapted = prepare._source_assets(Path(kit['original_atlasem_descriptor']))
    assert len(adapted['structures']) == 12
    config = json.loads((output / 'source-config.json').read_text())
    assert config['atlasem_descriptor'] == str(native)
    assert '--regrid-inputs' in script.read_text()
    assert '--no-regrid-inputs' not in script.read_text()
    with pytest.raises(FileExistsError):
        prepare.prepare(output, inputs)
