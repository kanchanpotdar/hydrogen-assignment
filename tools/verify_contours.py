"""Independent contour-node, interpolation and per-point evidence checks."""
import argparse
import csv
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import hyd as h
import hyd_contours as c
from numerical_checks import state_checks
from run_metadata import prepare_output,record_run


def main(workflow,out):
    directory=Path(workflow).resolve()/'contours';out=prepare_output(out,'contour_verification')
    with record_run(out,'contour_verification',dict(source_contours=str(directory),reference_tolerance_factor=.01)) as meta:
        rows=[];physical=[json.loads(p.read_text()) for p in (directory/'physical_checks').glob('*.json')]
        if len(physical)!=660 or not all(x['passed'] for x in physical):raise AssertionError('Missing or failed contour state checks')
        with np.load(directory/'operating_grid.npz') as f:
            T=f['T_K'];q=f['flow_ml_min_RTP'];fields=f['fields'];mask=f['CO_ratio_mask']
        assert fields.shape==(3,17,11,3)
        assert np.isfinite(fields[:,:,:,:2]).all()
        assert np.array_equal(mask,np.isnan(fields[:,:,:,2]))
        assert np.array_equal(mask,fields[:,:,:,0]<.05)
        meta['raw_conversion_min_max']=[float(fields[:,:,:,0].min()),float(fields[:,:,:,0].max())]
        meta['CO_mask_count']=int(mask.sum());meta['physical_cases']=len(physical)
        meta['max_C_H_O_He_errors']=np.max([x['relative_C_H_O_He_errors'] for x in physical],axis=0).tolist()
        meta['min_flow_mol_s']=min(x['min_flow_mol_s'] for x in physical)
        meta['temperature_range_K']=[min(x['min_T_K'] for x in physical),max(x['max_T_K'] for x in physical)]
        for i,cat in enumerate(c.CATS):
            for temp,flow in [(T[0],q[0]),(T[0],q[-1]),(T[-1],q[0]),(T[-1],q[-1]),(T[8],q[5]),(853.,120.)]:
                with patch.object(h,'TOL_SCALE',.01):res=h.run(cat,Tf=temp,flow=flow)
                s=h.summarize(res);F=res['PR2'][1][:6,-1];conv=s['X_CH4'];direct=np.array([conv,F[4]/F.sum(),F[2]/(res['F0'][0]-F[0])])
                interpolated=np.array([RegularGridInterpolator((T,q),fields[i,:,:,k])([[temp,flow]])[0] for k in range(3)])
                offgrid=not(temp in T and flow in q);limit=.005 if offgrid else 2e-5
                delta=float(np.max(abs(direct-interpolated)))
                rows.append(dict(catalyst=cat,kind='operating',T_K=temp,flow_ml_min_RTP=flow,max_species_or_ratio_difference=delta,max_T_difference_K=0.,criterion=limit,passed=bool(delta<=limit),description='off-grid bilinear interpolation' if offgrid else 'stored node vs tighter independent solve'))
            with np.load(directory/f'axial_grid_{cat}.npz') as f:ts=f['T_K'];z=f['z_mm'];raw=f['fields']
            assert raw.shape==(33,6,241) and np.isfinite(raw).all()
            for idx in [0,12,32]:
                temp=float(ts[idx]);base=h.run(cat,Tf=temp,flow=120,n=c.AXIAL_SAMPLE_POINTS);zz,Y=base['PBR'];y=h.mole_fractions(Y)
                independent=np.array([np.interp(z,zz*1000,values) for values in [y[0],y[4],y[2],y[3],y[1],Y[6]-temp]])
                d=float(abs(independent-raw[idx]).max())
                rows.append(dict(catalyst=cat,kind='axial',T_K=temp,flow_ml_min_RTP=120,max_species_or_ratio_difference=d,max_T_difference_K=float(abs(independent[-1]-raw[idx,-1]).max()),criterion=2e-5,passed=bool(d<=2e-5),description='stored map vs independently reconstructed V3 interpolation'))
        with (out/'node_comparisons.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        meta['verification_status']='PASS' if all(r['passed'] for r in rows) else 'FAIL'
        print('Contour verification:',meta['verification_status'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workflow',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();main(a.workflow,a.out)
