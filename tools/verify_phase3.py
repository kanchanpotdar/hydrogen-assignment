"""Execute Phase 3 regression without running the Phase 4 contour workflow."""
import csv
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import hyd as h
from run_metadata import prepare_output, record_run
from tools.baseline_regression import compare_csv


def main(out):
    out=prepare_output(out,'phase3')
    with record_run(out,'phase3',dict(scope='default numeric regression and presentation corrections, not Phase 4')) as manifest:
        manifest['commands']=[]
        # Execute from a different working directory to test caller-relative --out.
        relative=os.path.relpath(out/'model',ROOT.parent)
        commands=[('hyd.py',['--validate','--out',relative],ROOT.parent),
                  ('compare_h_length.py',['--out',str(out/'h_length_sensitivity')],ROOT)]
        for script,args,cwd in commands:
            cmd=[sys.executable,str(ROOT/script),*args]
            with (out/(script+'.log')).open('w') as log:
                p=subprocess.run(cmd,cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
            manifest['commands'].append(dict(argv=cmd,cwd=str(cwd),returncode=p.returncode))
            if p.returncode:raise RuntimeError(script+' failed')
        comparisons=[]
        for name,p in [('smr_summary.csv',out/'model/smr_summary.csv'),
                       ('validation_table.csv',out/'model/fixed_reference/validation_table.csv'),
                       ('validation_equilibrium.csv',out/'model/fixed_reference/validation_equilibrium.csv'),
                       ('h_length_comparison.csv',out/'h_length_sensitivity/h_length_comparison.csv')]:
            old=subprocess.check_output(['git','show','v3-baseline:results/'+name],cwd=ROOT)
            comparisons.append(compare_csv(name,old,p.read_bytes()))
        checks=[]
        for cat in ['HMMC','57-4Q','25-4Q']:
            res=h.run(cat)
            with np.load(ROOT/f'phase1_runs/2026-11-01-v3-reproduction/profiles_{cat}.npz') as old:
                equal=np.array_equal(res['F0'],old['F0'])
                for sec in ['PBR','PR1','PR2']:
                    equal=equal and np.array_equal(res[sec][0],old[sec+'_z_m']) and np.array_equal(res[sec][1],old[sec+'_state'])
            checks.append(dict(catalyst=cat,full_profile_array_equal=bool(equal)))
        # Protect every baseline input and historical output.
        names=subprocess.check_output(['git','ls-tree','-r','--name-only','v3-baseline'],cwd=ROOT,text=True).splitlines()
        for name in names:
            if name.startswith(('data/','results/')):
                assert (ROOT/name).read_bytes()==subprocess.check_output(['git','show','v3-baseline:'+name],cwd=ROOT),name
        report=dict(csv_comparisons=comparisons,profiles=checks,historical_inputs_outputs_unchanged=True,
                    tests_not_run=['full contours','convergence','independent solver','energy integral','all-catalyst seed sensitivity'])
        (out/'numeric_regression.json').write_text(json.dumps(report,indent=2)+'\n')
        if not all(x['passed'] for x in comparisons) or not all(x['full_profile_array_equal'] for x in checks):
            raise AssertionError('Numeric regression failed')
        print('Phase 3 default numeric regression: PASS')


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True)
    main(ap.parse_args().out)
