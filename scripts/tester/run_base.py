#!/usr/bin/env python3
"""Run the packaged full BASE defaults and collect feedback automatically."""
from __future__ import annotations
import argparse
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

BUNDLE = Path(__file__).resolve().parent
THREADS = ('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')


def check_environment(release):
    expected_python=tuple(int(v) for v in release['tested_python'].split('.')[:2])
    if sys.version_info[:2]!=expected_python:
        raise RuntimeError('Use Python '+'.'.join(map(str,expected_python))+' for this constrained testing round.')
    differences=[]
    for name,expected in release['tested_scientific_dependencies'].items():
        try:actual=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:actual='missing'
        if actual!=expected:differences.append(f'{name}: expected {expected}, found {actual}')
    if differences:
        raise RuntimeError('Scientific dependencies differ from this tested release: '+'; '.join(differences)+'. Reinstall with constraints-tested.txt in the fresh environment.')


def command(inputs,output,atlas,*,resume=False,n4_estimation_mask='all'):
    args = [sys.executable,'-u','-c','from pymacruise.cli import main; raise SystemExit(main())',
        'base-run','--t1',str(inputs/'T1.nii.gz'),'--brain-mask',str(inputs/'brain_mask.nii.gz'),
        '--maseg',str(inputs/'MaSeg.nii.gz'),'--atlas-descriptor',str(atlas),
        '--output-dir',str(output),'--n4-estimation-mask',n4_estimation_mask]
    if resume:args.append('--resume')
    return args


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=Path('inputs'))
    parser.add_argument('--output-dir',type=Path,default=Path('outputs_pymacruise_0.1.6'))
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--n4-estimation-mask',choices=('all','legacy-otsu'),default='all')
    args = parser.parse_args(argv)
    inputs=args.input_dir.expanduser().resolve();output=args.output_dir.expanduser().resolve()
    for name in ('T1.nii.gz','brain_mask.nii.gz','MaSeg.nii.gz'):
        if not (inputs/name).is_file():parser.error(f'Missing input: {inputs/name}')
    if output.exists() and any(output.iterdir()) and not args.resume:
        parser.error('Use a fresh output folder; --resume is only for the same package/settings run.')
    if args.resume and not (output/'run.json').is_file():parser.error('No run.json exists to resume.')
    release=json.loads((BUNDLE/'release.json').read_text())
    for key in THREADS:os.environ[key]='1'
    os.environ.pop('PYTHONPATH',None);os.environ['PYTHONNOUSERSITE']='1';os.environ['PYTHONUNBUFFERED']='1'
    output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();code=1;error=None;child=None
    log_mode='a' if args.resume else 'w'
    try:
        installed=importlib.metadata.version('pymacruise')
        if installed!=release['version']:raise RuntimeError(f'Expected PyMaCRUISE {release["version"]}; found {installed}. Install the included wheel in the fresh environment.')
        check_environment(release)
        from pymacruise.base_workflow import _implementation_signature
        if _implementation_signature()!=release['implementation_signature']:
            raise RuntimeError('Installed Python source differs from this tested package. Reinstall the included wheel.')
        import ants  # Check the N4 extension before starting the full run.
        from pymacruise.atlas_resources import get_bundled_atlas_descriptor
        atlas=get_bundled_atlas_descriptor('cruise-atlas-10obj-toads2008.txt')
        args_list=command(inputs,output,atlas,resume=args.resume,n4_estimation_mask=args.n4_estimation_mask)
        with (output/'tester-run.log').open(log_mode,encoding='utf-8') as log:
            heading=f'PyMaCRUISE {installed}; full BASE defaults; N4 policy={args.n4_estimation_mask}; numerical threads=1\n'
            print(heading,end='',flush=True);log.write(heading);log.flush()
            child=subprocess.Popen(args_list,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                text=True,bufsize=1,env=os.environ.copy())
            for line in child.stdout:
                print(line,end='',flush=True);log.write(line);log.flush()
            code=child.wait()
    except KeyboardInterrupt:
        error='Interrupted by user; completed stage checkpoints are preserved.';code=130
        if child is not None and child.poll() is None:
            child.terminate();child.wait(timeout=30)
    except Exception as exception:
        error=str(exception);print(error,file=sys.stderr,flush=True)
    summary={'package_version':release['version'],'pipeline_returncode':code,'returncode':code,'error':error,
        'elapsed_seconds':time.monotonic()-started,'resume':args.resume,
        'n4_estimation_mask':args.n4_estimation_mask,'numerical_threads':1}
    (output/'tester-run-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    spec=importlib.util.spec_from_file_location('tester_feedback',BUNDLE/'collect-feedback.py')
    collector=importlib.util.module_from_spec(spec);spec.loader.exec_module(collector)
    archive,feedback=collector.collect(output,t1=inputs/'T1.nii.gz',bundle=BUNDLE)
    print('Feedback archive:',archive,flush=True)
    print('Run checks:','passed' if feedback['verification']['passed'] else 'failed or incomplete; diagnostics included',flush=True)
    return code if code else (0 if feedback['verification']['passed'] else 1)

if __name__=='__main__':raise SystemExit(main())
