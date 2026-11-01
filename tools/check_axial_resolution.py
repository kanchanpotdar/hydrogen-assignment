"""Compare stored axial-map interpolation to independently evaluated Radau states."""
import argparse,csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.integrate import solve_ivp
import hyd as h
from numerical_checks import verification_exit_code
from run_metadata import prepare_output,record_run


def main(workflow,out):
    out=prepare_output(out,'axial_resolution')
    with record_run(out,'axial_resolution',dict(reference='Radau rtol1e-9, species atol1e-17, T atol1e-7')) as meta:
        rows=[]
        contour_manifest=json.loads((Path(workflow)/'contours/manifest_contours.json').read_text())
        stored_n=contour_manifest['settings'].get('PBR_axial_sample_points',400)
        for cat in ['HMMC','57-4Q','25-4Q']:
            with np.load(Path(workflow)/'contours'/f'axial_grid_{cat}.npz') as data:
                temperatures=data['T_K'];z=data['z_mm'];raw=data['fields']
            for idx in [0,12,32]:
                T=float(temperatures[idx]);r=h.SMR_REACTOR;y=np.array([.05,.2,0,0,1e-6,.75]);y/=y.sum()
                F=h.feed_moles(120,y)
                sol=solve_ivp(h.rhs_pbr,(0,h.area(r)*r.L_pbr),np.r_[F,T],args=(r,T,h.ZETA[cat]),method='Radau',rtol=1e-9,atol=np.r_[np.full(6,1e-17),1e-7],dense_output=True)
                if not sol.success:raise RuntimeError(sol.message)
                Y=sol.sol(z*.001*h.area(r));yy=h.mole_fractions(Y)
                ref=np.array([yy[0],yy[4],yy[2],yy[3],yy[1],Y[6]-T])
                for n in sorted(set([400,1600,3200,stored_n])):
                    if n==stored_n:estimate=raw[idx]
                    else:
                        res=h.run(cat,Tf=T,n=n);zz,B=res['PBR'];by=h.mole_fractions(B)
                        estimate=np.array([np.interp(z,zz*1000,v) for v in [by[0],by[4],by[2],by[3],by[1],B[6]-T]])
                    ds=float(abs(estimate[:5]-ref[:5]).max());dt=float(abs(estimate[5]-ref[5]).max())
                    rows.append(dict(catalyst=cat,T_K=T,PBR_saved_points=n,max_species_difference=ds,max_T_difference_K=dt,passed=bool(ds<=2e-5 and dt<=.1)))
                print(cat,T,'axial-resolution checked',flush=True)
        with (out/'axial_resolution.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        meta['verification_status']='PASS' if all(r['passed'] for r in rows if r['PBR_saved_points']==stored_n) else 'FAIL'
    return verification_exit_code(meta)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workflow',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();raise SystemExit(main(a.workflow,a.out))
