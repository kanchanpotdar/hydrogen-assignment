"""Final pipeline inventory, V3 regression and immutable manifest audit."""
import argparse,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from matplotlib.image import imread
from run_metadata import prepare_output,record_run,sha
from tools.baseline_regression import compare_csv


def main(workflow,out):
    base=Path(workflow).resolve();out=prepare_output(out,'artifact_audit')
    with record_run(out,'artifact_audit',dict(workflow=str(base))) as meta:
        manifests=list(base.rglob('manifest_*.json'))
        if len(manifests)!=5:raise AssertionError('Expected workflow plus four stage manifests')
        for p in manifests:
            d=json.loads(p.read_text());assert d['status']=='COMPLETED',p
            for name,digest in d['output_sha256'].items():assert sha(p.parent/name)==digest,(p,name)
        figures=list(base.rglob('*.png'));assert len(figures)==10
        assert len(list(base.rglob('*.csv')))==7
        assert len(list(base.rglob('*.npz')))==7
        assert len(list((base/'contours/physical_checks').glob('*.json')))==660
        meta['figures']=[str(p.relative_to(base)) for p in figures]
        inventory=[]
        for p in sorted(base.rglob('*')):
            if p.is_file():
                item=dict(path=str(p.relative_to(base)),bytes=p.stat().st_size,sha256=sha(p))
                if p.suffix=='.png':
                    pixels=imread(p);assert pixels.ndim==3 and min(pixels.shape[:2])>500
                    item['pixel_dimensions']=list(pixels.shape)
                inventory.append(item)
        (out/'file_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
        meta['inventory_files']=len(inventory)
        comparisons=[]
        for name,stage in [('smr_summary.csv','model'),('validation_table.csv','fixed_reference'),('validation_equilibrium.csv','fixed_reference'),('h_length_comparison.csv','h_length_sensitivity')]:
            old=subprocess.check_output(['git','show','v3-baseline:results/'+name],cwd=ROOT)
            comparisons.append(compare_csv(name,old,(base/stage/name).read_bytes()))
        meta['CSV_regression']=comparisons
        meta['profile_array_equal']={}
        for cat in ['HMMC','57-4Q','25-4Q']:
            with np.load(base/'model'/f'profiles_{cat}.npz') as f,np.load(ROOT/'phase1_runs/2026-11-01-v3-reproduction'/f'profiles_{cat}.npz') as g:
                equal=all(np.array_equal(f[k],g[k]) for k in g.files)
                meta['profile_array_equal'][cat]=equal
        with np.load(base/'contours/operating_grid.npz') as current,np.load(ROOT/'phase4_runs/2026-11-01-workflow/contours/operating_grid.npz') as old:
            meta['operating_grid_unchanged']=bool(np.array_equal(current['fields'],old['fields'],equal_nan=True))
        meta['axial_map_before_after']={}
        for cat in ['HMMC','57-4Q','25-4Q']:
            with np.load(base/'contours'/f'axial_grid_{cat}.npz') as current,np.load(ROOT/'phase4_runs/2026-11-01-workflow/contours'/f'axial_grid_{cat}.npz') as old:
                delta=abs(current['fields']-old['fields'])
                meta['axial_map_before_after'][cat]=dict(max_species_change=float(delta[:,:5,:].max()),max_T_change_K=float(delta[:,5,:].max()))
        names=subprocess.check_output(['git','ls-tree','-r','--name-only','v3-baseline'],cwd=ROOT,text=True).splitlines()
        for name in names:
            if name.startswith(('data/','results/')):
                assert (ROOT/name).read_bytes()==subprocess.check_output(['git','show','v3-baseline:'+name],cwd=ROOT)
        assert all(r['passed'] for r in comparisons) and all(meta['profile_array_equal'].values()) and meta['operating_grid_unchanged']
        meta['verification_status']='PASS';print('Artifact inventory/regression: PASS')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workflow',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();main(a.workflow,a.out)
