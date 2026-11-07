"""Tester tools preserve full defaults and collect diagnostics without images."""
from pathlib import Path
import importlib.util,json,zipfile
import sys
import gzip
import hashlib
import nibabel as nib
import numpy as np
import pytest

ROOT=Path(__file__).parents[1]/'scripts/tester'
def load(name):
    spec=importlib.util.spec_from_file_location('tester_'+name,ROOT/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_feedback_whitelist_excludes_scientific_data_and_external_links(tmp_path,monkeypatch):
    module=load('collect_feedback');run=tmp_path/'run';stage=run/'stages/core/attempt-01';stage.mkdir(parents=True)
    (run/'run.json').write_text('{"status":"failed"}')
    (stage/'stage.log').write_text('specific processing error\n')
    (stage/'provenance.json').write_text('{}')
    for name in ('T1.nii.gz','surface.vtk','secret.npy','screenshot.png'):(stage/name).write_bytes(b'not for feedback')
    outside=tmp_path/'outside.log';outside.write_text('outside content')
    (run/'tester-run.log').symlink_to(outside)
    monkeypatch.setenv('PRIVATE_TOKEN','must-not-be-in-feedback')
    path,summary=module.collect(run,verify=False)
    with zipfile.ZipFile(path) as stream:
        assert set(stream.namelist())=={'run.json','stages/core/attempt-01/stage.log','stages/core/attempt-01/provenance.json','feedback-environment.json'}
        assert b'must-not-be-in-feedback' not in stream.read('feedback-environment.json')
    assert summary['verification']['performed'] is False


def test_incomplete_run_still_produces_feedback(tmp_path):
    module=load('collect_feedback');run=tmp_path/'run';run.mkdir();(run/'run.json').write_text('{"status":"failed"}')
    checker=tmp_path/'verify-base-run.py';checker.write_text('import sys,json\nfrom pathlib import Path\nPath(sys.argv[sys.argv.index("--report")+1]).write_text(json.dumps({"artifact_and_grid_checks_passed":False,"error":"run incomplete"}))\nraise SystemExit(1)\n')
    archive,summary=module.collect(run,t1=tmp_path/'T1.nii.gz',bundle=tmp_path)
    assert summary['verification']['performed'] is True
    assert summary['verification']['passed'] is False
    with zipfile.ZipFile(archive) as stream:assert 'checks.json' in stream.namelist()


def test_launcher_uses_argument_list_full_defaults_and_explicit_n4_policy(tmp_path):
    module=load('run_base');inputs=tmp_path/'inputs with spaces';output=tmp_path/'new output';atlas=tmp_path/'atlas file.txt'
    args=module.command(inputs,output,atlas,resume=True)
    assert '--resume' in args
    assert args[args.index('--t1')+1]==str(inputs/'T1.nii.gz')
    assert args[args.index('--atlas-descriptor')+1]==str(atlas)
    assert args[args.index('--n4-estimation-mask')+1]=='all'
    for option in ('--toads-iterations','--tgdm-iterations','--gvf-iterations','--dura-iterations'):assert option not in args


def test_launcher_rejects_dependency_drift_before_processing(monkeypatch):
    module=load('run_base')
    monkeypatch.setattr(module.importlib.metadata,'version',lambda name:'unexpected-version')
    release={'tested_python':'.'.join(map(str,sys.version_info[:2])),
             'tested_scientific_dependencies':{'numpy':'2.2.6'}}
    with pytest.raises(RuntimeError,match='constraints-tested.txt'):
        module.check_environment(release)


def test_label_metrics_handle_sparse_large_label_ids_without_dense_joint_histogram():
    module=load('compare_base');a=np.array([0,44,44,999999999]);b=np.array([0,44,207,999999999])
    result=module.label_metrics(a,b)
    assert result['differing_voxels']==1
    assert result['per_label']['44']['dice']==pytest.approx(2/3)
    assert result['per_label']['207']['dice']==0
    assert result['per_label']['999999999']['dice']==1


@pytest.mark.parametrize('values',[np.array([.5,1]),np.array([-1,1]),np.array([np.nan,1])])
def test_comparison_rejects_invalid_original_labels(values):
    with pytest.raises(ValueError,match='finite, nonnegative integers'):load('compare_base').label_metrics(values,np.ones(2))


def test_geometry_distinguishes_header_rounding_from_physical_shift():
    module=load('compare_base');affine=np.diag([.89,.89,1.,1.]);rounded=affine.copy();rounded[0,0]+=6e-8
    a=nib.Nifti1Image(np.zeros((256,256,150),dtype=np.uint8),rounded)
    b=nib.Nifti1Image(np.zeros((256,256,150),dtype=np.uint8),affine)
    report=module.geometry_metrics(a,b)
    assert not report['affines_bit_identical']
    assert report['same_lattice_within_0_0001_mm']
    assert report['maximum_corner_displacement_mm']==pytest.approx(255*6e-8,rel=.001)
    shifted=affine.copy();shifted[0,3]=1.
    report=module.geometry_metrics(nib.Nifti1Image(np.zeros((2,2,2)),shifted),b)
    assert report['maximum_corner_displacement_mm']==pytest.approx(1)
    assert not report['same_lattice_within_0_0001_mm']


def test_feedback_records_pipeline_completion_separately_from_failed_verification(tmp_path):
    module=load('collect_feedback');run=tmp_path/'run';run.mkdir()
    (run/'run.json').write_text('{"status":"complete","converged":false}')
    (run/'tester-run-summary.json').write_text('{"returncode":0}')
    (tmp_path/'verify-base-run.py').write_text('import sys\nraise SystemExit(1)\n')
    archive,summary=module.collect(run,t1=tmp_path/'T1.nii.gz',bundle=tmp_path)
    with zipfile.ZipFile(archive) as stream:execution=json.loads(stream.read('tester-run-summary.json'))
    assert execution['pipeline_returncode']==0
    assert execution['returncode']==1
    assert execution['verification_passed'] is False
    assert summary['verification']['passed'] is False


def test_voxel_fingerprints_ignore_gzip_recompression_and_detect_changed_values(tmp_path):
    module=load('collect_feedback');data=np.arange(24,dtype=np.float32).reshape(2,3,4)
    a=tmp_path/'a.nii.gz';b=tmp_path/'b.nii.gz';c=tmp_path/'c.nii.gz'
    nib.save(nib.Nifti1Image(data,np.eye(4)),a)
    b.write_bytes(gzip.compress(gzip.decompress(a.read_bytes()),compresslevel=1,mtime=123))
    altered=data.copy();altered[0,0,0]=7.;nib.save(nib.Nifti1Image(altered,np.eye(4)),c)
    assert hashlib.sha256(a.read_bytes()).digest()!=hashlib.sha256(b.read_bytes()).digest()
    manifest={'stages':{'n4':{'outputs':[{'path':str(p)} for p in (a,b,c)]}}}
    report=module.image_diagnostics(tmp_path,manifest)['images']
    assert report['n4/a.nii.gz']['voxel_sha256']==report['n4/b.nii.gz']['voxel_sha256']
    assert report['n4/a.nii.gz']['voxel_sha256']!=report['n4/c.nii.gz']['voxel_sha256']
    assert report['n4/a.nii.gz']['mean']==pytest.approx(data.mean())
    assert not any(key in report['n4/a.nii.gz'] for key in ('data','values','samples','voxels'))


@pytest.fixture
def review_case(tmp_path):
    from pymacruise.surfaces import SurfaceMesh,write_legacy_vtk
    run=tmp_path/'run';source=tmp_path/'original';inputs=tmp_path/'inputs'
    run.mkdir();inputs.mkdir()
    for folder in ('Raw','SegRefine','Surfaces','Surfaces_FreeView'):(source/folder).mkdir(parents=True)
    affine=np.array([[.89,0,0,-2],[0,.89,0,-3],[0,0,1,-4],[0,0,0,1]],dtype=float)
    labels=np.arange(8,dtype=np.int16).reshape(2,2,2);t1=nib.Nifti1Image(labels,affine)
    nib.save(t1,inputs/'T1.nii.gz');nib.save(t1,source/'Raw/T1.nii.gz')
    nib.save(t1,inputs/'MaSeg.nii.gz');nib.save(t1,source/'Raw/Seg.nii.gz')
    rounded=affine.copy();rounded[0,0]+=6e-8
    nib.save(nib.Nifti1Image(labels,rounded),source/'SegRefine/SEG_refine.nii.gz')
    candidate=labels.copy();candidate[0,0,0]=12
    nib.save(nib.Nifti1Image(candidate,affine),run/'SEG_refine.nii.gz')
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
    faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
    world=np.einsum('ij,nj->ni',nib.load(source/'Raw/Seg.nii.gz').affine[:3,:3],vertices)+affine[:3,3]
    products=[{'path':str(run/'SEG_refine.nii.gz')}]
    for name in ('inner','central','outer'):
        write_legacy_vtk(SurfaceMesh(vertices,faces),source/f'Surfaces/target_image_GMimg_{name}Surf.vtk')
        with (source/f'Surfaces_FreeView/target_image_GMimg_{name}Surf.asc').open('w') as stream:
            stream.write('# physical\n4 4\n');np.savetxt(stream,world)
        path=run/(name+'_surface.vtk');write_legacy_vtk(SurfaceMesh(world+.5,faces),path)
        products.append({'path':str(path)})
    manifest={'status':'complete','converged':False,'inputs':{'t1':{'path':str(inputs/'T1.nii.gz')},
        'maseg':{'path':str(inputs/'MaSeg.nii.gz')}},'final_outputs':products}
    (run/'run.json').write_text(json.dumps(manifest))
    return run,source,labels,candidate,world,faces


def test_paired_review_preserves_labels_faces_and_physical_coordinates(review_case,tmp_path):
    from pymacruise.surfaces import read_legacy_vtk
    run,source,original,candidate,world,faces=review_case;output=tmp_path/'review'
    report=load('prepare_review').prepare(run,source,output)
    a=nib.load(output/'original_SEG_refine.nii.gz');b=nib.load(output/'python_SEG_refine.nii.gz')
    np.testing.assert_array_equal(a.affine,b.affine)
    np.testing.assert_array_equal(a.dataobj,original);np.testing.assert_array_equal(b.dataobj,candidate)
    old=read_legacy_vtk(output/'original_central_RAS_mm.vtk')
    new=read_legacy_vtk(output/'python_central_RAS_mm.vtk')
    np.testing.assert_array_equal(old.faces,faces);np.testing.assert_array_equal(new.faces,faces)
    np.testing.assert_allclose(old.vertices,world,rtol=0,atol=1e-9)
    np.testing.assert_allclose(new.vertices,world+.5,rtol=0,atol=1e-9)
    assert report['source_run_converged'] is False
    assert report['fitted_registration'] is False
    assert (run/'review-geometry.json').is_file()
    with pytest.raises(FileExistsError):load('prepare_review').prepare(run,source,output)


def test_paired_review_refuses_real_grid_shift(review_case,tmp_path):
    run,source,*_=review_case;p=source/'SegRefine/SEG_refine.nii.gz';im=nib.load(p)
    affine=im.affine.copy();affine[0,3]+=1
    nib.save(nib.Nifti1Image(np.asarray(im.dataobj),affine),p)
    with pytest.raises(ValueError,match='refusing header-only'):
        load('prepare_review').prepare(run,source,tmp_path/'review')
    assert not (tmp_path/'review').exists()


def test_paired_review_refuses_unknown_original_surface_convention(review_case,tmp_path):
    run,source,*_=review_case
    (source/'Surfaces_FreeView/target_image_GMimg_innerSurf.asc').write_text('# wrong frame\n4 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n')
    with pytest.raises(ValueError,match='does not map'):
        load('prepare_review').prepare(run,source,tmp_path/'review')
    assert not (tmp_path/'review').exists()
