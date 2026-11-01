"""Isolated, non-overwriting V3 reproduction; never alters model or baseline outputs."""
import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = 'v3-baseline'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def tolerance(filename, column):
    if column == 'digitization_uncertainty' or column.startswith('n_'):
        return 0.0
    if filename == 'smr_summary.csv':
        return {'T_K': .02, 'X_CH4': 2e-4}.get(column, 5e-6)
    if filename == 'validation_table.csv':
        return .2 if column == 'rel_err_pct' else 2e-5
    if filename == 'h_length_comparison.csv':
        return {'U_inlet_W_m2K': 2., 'T_dip_K': .2, 'T_dip_paper_text_K': .2,
                'T_curve_RMSE_K': .2}.get(column, 2e-5)
    return 2e-5

def compare(filename, original, current):
    a = csv.DictReader(io.StringIO(original.decode()))
    b = csv.DictReader(io.StringIO(current.decode()))
    ar, br = list(a), list(b)
    issues = []
    maxima = {}
    if a.fieldnames != b.fieldnames or len(ar) != len(br):
        issues.append('schema or row count mismatch')
    for i, (x, y) in enumerate(zip(ar, br), 1):
        for k in x:
            if k not in y:
                issues.append(f'row {i}: missing {k}')
                continue
            try:
                u, v = float(x[k]), float(y[k])
            except (ValueError, TypeError):
                if x[k] != y[k]:
                    issues.append(f'row {i} {k}: {x[k]} != {y[k]}')
                continue
            d = abs(u - v)
            maxima[k] = max(maxima.get(k, 0.), d)
            if not math.isfinite(u) or not math.isfinite(v) or d > tolerance(filename, k) + 1e-14:
                issues.append(f'row {i} {k}: delta {d}')
    return dict(file=filename, byte_equal=original == current, rows=len(br),
                max_absolute_difference=maxima, issues=issues, passed=not issues)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out == ROOT or ROOT / 'results' == out or out.exists():
        raise SystemExit('Output must be a new directory, separate from preserved results.')
    out.mkdir(parents=True)
    manifest = dict(baseline_commit=git('rev-parse', BASE+'^{}').decode().strip(),
                    execution_commit=git('rev-parse', 'HEAD').decode().strip(),
                    python=sys.version, executable=sys.executable, platform=platform.platform(),
                    started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    configuration=dict(T_K=853., flow_ml_min_RTP=120., pressure_bar=1.01325,
                                       h_length='tube', KP_MODEL='D', solver='LSODA',
                                       rtol=1e-8, atol_species=1e-16, atol_T=1e-6, h2_seed=1e-6),
                    packages={d.metadata['Name']: d.version for d in importlib.metadata.distributions()},
                    working_tree=git('status', '--porcelain').decode(), files={}, commands=[], status='RUNNING')
    dump(out/'manifest.json', manifest)
    try:
        paths = git('ls-tree', '-r', '--name-only', BASE).decode().splitlines()
        for p in paths:
            if p.endswith('.py') or p.startswith(('data/', 'results/')):
                original = git('show', f'{BASE}:{p}')
                current = (ROOT/p).read_bytes()
                manifest['files'][p] = dict(baseline_sha256=sha(original), current_sha256=sha(current))
                if original != current:
                    raise RuntimeError(f'Baseline changed: {p}')
        env = dict(os.environ, MPLCONFIGDIR=str(out/'mplconfig'), PYTHONNOUSERSITE='1',
                   PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
        for script in ['hyd.py', 'validation.py', 'compare_h_length.py']:
            cmd = [sys.executable, str(ROOT/script), '--out', str(out)]
            start = time.monotonic()
            with (out/(script+'.log')).open('w') as log:
                proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            manifest['commands'].append(dict(argv=cmd, cwd=str(ROOT), returncode=proc.returncode,
                                              seconds=time.monotonic()-start))
            dump(out/'manifest.json', manifest)
            print(script, 'exit', proc.returncode, flush=True)
            if proc.returncode:
                raise RuntimeError(f'{script} failed; see retained log')
        comparisons = []
        for f in ['smr_summary.csv', 'validation_table.csv', 'validation_equilibrium.csv', 'h_length_comparison.csv']:
            comparisons.append(compare(f, git('show', f'{BASE}:results/{f}'), (out/f).read_bytes()))
        dump(out/'csv_comparison.json', comparisons)
        # Additional explicit model executions solely for saved-state sanity and detailed metrics.
        sys.path.insert(0, str(ROOT))
        import numpy as np
        import hyd as h
        import validation as v
        exp, curves = v.load_experiment(), v.load_paper_curves()
        sanity, metrics = {}, {}
        atoms = np.array([[1,0,1,1,0,0], [4,2,0,0,2,0], [0,1,1,2,0,0], [0,0,0,0,0,1]])
        for cat in v.CATS:
            res = h.run(catalyst=cat)
            states = np.hstack([res[s][1] for s in ['PBR', 'PR1', 'PR2']])
            initial = atoms @ res['F0']
            residual = np.max(np.abs((atoms @ states[:6] - initial[:,None])/initial[:,None]), axis=1)
            norm = float(np.max(np.abs(h.mole_fractions(states).sum(axis=0)-1)))
            minimum = float(states[:6].min())
            passed = bool(np.isfinite(states).all() and minimum >= -1e-14 and norm <= 1e-12 and residual.max() <= 1e-8)
            sanity[cat] = dict(passed=passed, min_raw_flow_mol_s=minimum, normalization_error=norm,
                               relative_C_H_O_He_errors=residual.tolist())
            np.savez_compressed(out/f'profiles_{cat}.npz', F0=res['F0'], **{s+'_'+k: res[s][i] for s in ['PBR','PR1','PR2'] for i,k in enumerate(['z_m','state'])})
            outlet = v.exit_composition(res)
            e = np.array([outlet[s]-exp[cat][s][0] for s in v.SPECIES])
            pct = np.array([100*abs(outlet[s]-exp[cat][s][0])/exp[cat][s][0] for s in v.SPECIES])
            z,y,T,b = h.profiles(res); mask=z<=b['PBR_end']
            rmse = {}
            for s in v.SPECIES+['T']:
                d=curves[(s,cat)]; d=d[d[:,0]>(.15 if s=='T' else .05)]
                pred=np.interp(d[:,0],z[mask],(T if s=='T' else y[v.IDX[s]])[mask])
                rmse[s]=dict(n=len(d), z_min_mm=float(d[:,0].min()), z_max_mm=float(d[:,0].max()), RMSE=float(np.sqrt(np.mean((pred-d[:,1])**2))))
            metrics[cat]=dict(experiment_MAE=float(abs(e).mean()), experiment_RMSE=float(np.sqrt(np.mean(e**2))), experiment_MAPE_pct=float(pct.mean()), PBR_curve_metrics=rmse)
        dump(out/'sanity_checks.json', sanity); dump(out/'validation_metrics.json', metrics)
        manifest['diagnostic_execution']='Three additional default h.run calls: all catalysts; conservation and detailed metrics.'
        passed = all(c['passed'] for c in comparisons) and all(c['passed'] for c in sanity.values())
        manifest['status']='PASS' if passed else 'FAIL'
        # Recheck protected files after executions.
        for p, hashes in manifest['files'].items():
            if sha((ROOT/p).read_bytes()) != hashes['baseline_sha256']:
                raise RuntimeError(f'Protected file altered: {p}')
        manifest['output_hashes']={p.name:sha(p.read_bytes()) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}
    except Exception as exc:
        manifest['status']='FAIL'; manifest['failure']=repr(exc)
        raise
    finally:
        manifest['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        dump(out/'manifest.json',manifest)
    print('Phase 1 numerical reproduction:', manifest['status'], flush=True)
    return 0 if manifest['status']=='PASS' else 1

if __name__ == '__main__':
    raise SystemExit(main())
