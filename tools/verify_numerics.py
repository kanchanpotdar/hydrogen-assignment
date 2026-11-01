"""Phase 4 independent numerical experiments, preserving all production physics."""
import argparse
import csv
import json
from pathlib import Path
import sys
import time
import warnings

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.optimize import minimize_scalar
import hyd as h
from numerical_checks import state_checks, verification_exit_code
from run_metadata import prepare_output, record_run


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def grid(L,n):
    return np.unique(np.r_[0,np.geomspace(L*1e-12,L,n//2),np.linspace(0,L,n//2)])


def solve_chain(cat,method='LSODA',factor=1.,seed=1e-6):
    r=h.SMR_REACTOR
    y=np.array([.05,.2,0,0,seed,.75]);y/=y.sum()
    F0=h.feed_moles(120,y);initial=np.r_[F0,853.]
    specs=[('PBR',h.rhs_pbr,r.L_pbr,(r,853.,h.ZETA[cat])),
           ('PR1',h.rhs_pr1,r.L_pr1,(r,853.)),('PR2',h.rhs_pr2,r.L_pr2,(r,523.))]
    out={'F0':F0};solutions={}
    for name,fun,L,args in specs:
        sol=solve_ivp(fun,(0,h.area(r)*L),initial,args=args,method=method,
                      rtol=1e-8*factor,atol=np.r_[np.full(6,1e-16),1e-6]*factor,dense_output=True)
        if not sol.success:raise RuntimeError(f'{cat}/{method}/{seed}/{name}: {sol.message}')
        zz=grid(L,4000);Y=sol.sol(zz*h.area(r));out[name]=(zz,Y)
        solutions[name]=(sol,L);initial=sol.y[:,-1]
    return out,solutions


def difference(a,b):
    dy=max(float(np.max(abs(h.mole_fractions(a[s][1])-h.mole_fractions(b[s][1])))) for s in ['PBR','PR1','PR2'])
    dt=max(float(np.max(abs(a[s][1][6]-b[s][1][6]))) for s in ['PBR','PR1','PR2'])
    dx=abs(h.summarize(a)['X_CH4']-h.summarize(b)['X_CH4'])
    return dict(max_species_difference=dy,max_T_difference_K=dt,conversion_difference=dx,
                passed=bool(dy<=2e-5 and dt<=.1 and dx<=1e-5))


def dense_minimum(solutions):
    sol,L=solutions['PBR'];z=grid(L,8000);T=sol.sol(z*h.area(h.SMR_REACTOR))[6];k=int(T.argmin())
    lo,hi=z[max(0,k-1)],z[min(len(z)-1,k+1)]
    fit=minimize_scalar(lambda zz:sol.sol(zz*h.area(h.SMR_REACTOR))[6],bounds=(lo,hi),method='bounded',options={'xatol':1e-15})
    return float(min(T.min(),fit.fun)),float(fit.x*1000)


def energy_terms(cat,solutions,n):
    """Independent thermochemical endpoint/line/quadrature checks, all terms W."""
    r=h.SMR_REACTOR;A=h.area(r);records=[]
    for sec,(sol,L) in solutions.items():
        zz=np.unique(np.r_[0,np.geomspace(L*1e-18,L,n//2),np.linspace(0,L,n//2)])
        V=zz*A;Y=sol.sol(V);F=Y[:6];T=Y[6]
        qwall=[];qreaction=[];Fcp=[]
        for k,temp in enumerate(T):
            flow=F[:,k];y=flow/flow.sum();p=y*h.P_BAR
            if sec=='PBR':
                U=h.U_pbr(flow,temp,r)[0];wall=853.
                rates=h.rates_catalyst(p,temp,h.ZETA[cat])*r.m_cat/(A*r.L_pbr)*1000/3600
                # Stoichiometry independent of the production NU array.
                sf=np.array([-rates[0]-rates[2],-rates[0]-rates[1]-2*rates[2],rates[0]-rates[1],rates[1]+rates[2],3*rates[0]+rates[1]+4*rates[2],0.])
            else:
                U=h.U_pipe(flow,temp,r,r.kw_glass if sec=='PR1' else r.kw_steel);wall=853. if sec=='PR1' else 523.
                rw=h.rate_wgs_gas(y,temp) if sec=='PR1' else h.rate_wgs_wall(p,temp)*1000/3600*(4/r.d)
                sf=np.array([0.,-rw,-rw,rw,rw,0.])
            qwall.append(U*(4/r.d)*(wall-temp))
            qreaction.append(-float(h.hf_i(temp)@sf))
            Fcp.append(float(flow@h.cp_i(temp)))
        qw=float(simpson(qwall,x=V));qr=float(simpson(qreaction,x=V))
        dh=float(F[:,-1]@h.hf_i(T[-1])-F[:,0]@h.hf_i(T[0]))
        # Independent temperature line integral: trapezoid in T, not the ODE dT/dV.
        thermal=float(np.sum(.5*(np.asarray(Fcp)[:-1]+np.asarray(Fcp)[1:])*np.diff(T)))
        scale=max(abs(qw),abs(qr),abs(dh),1e-6)
        records.append(dict(catalyst=cat,section=sec,n=len(V),wall_heat_W=qw,reaction_heat_W=qr,
                            delta_total_enthalpy_W=dh,sensible_line_integral_W=thermal,
                            total_energy_residual_W=dh-qw,thermal_residual_W=thermal-qw-qr,
                            normalization_W=scale,total_relative_residual=abs(dh-qw)/scale,
                            thermal_relative_residual=abs(thermal-qw-qr)/scale))
    return records


def verification_passed(comparisons, energy, sanity, failures):
    return bool(comparisons and energy and sanity and not failures and
                all(x['passed'] for x in comparisons+energy+sanity))


def main(out):
    out=prepare_output(out,'numerical_verification')
    with record_run(out,'numerical_verification',dict(configuration='v3_reproduction',protocol='docs/phase4/acceptance_criteria.md')) as manifest:
        comparisons=[];minima=[];sanity=[];energy=[];failures=[]
        for cat in ['HMMC','57-4Q','25-4Q']:
            cache={}
            cases=[('tol_1','LSODA',1.,1e-6),('tol_0.1','LSODA',.1,1e-6),('tol_0.01','LSODA',.01,1e-6),
                   ('energy_tol_0.001','LSODA',.001,1e-6),('Radau_0.1','Radau',.1,1e-6),('seed_1e-8','LSODA',.1,1e-8),('seed_1e-10','LSODA',.1,1e-10),('seed_1e-12','LSODA',.1,1e-12)]
            for label,method,factor,seed in cases:
                start=time.monotonic()
                try:
                    res,sol=solve_chain(cat,method,factor,seed);cache[label]=(res,sol)
                    check=state_checks(res);sanity.append(dict(catalyst=cat,case=label,**check))
                    value,loc=dense_minimum(sol)
                    minima.append(dict(catalyst=cat,case=label,dense_Tmin_K=value,location_mm=loc))
                    np.savez_compressed(out/f'profiles_{cat}_{label}.npz',F0=res['F0'],**{s+'_'+name:res[s][i] for s in ['PBR','PR1','PR2'] for i,name in enumerate(['z_m','state'])})
                    print(cat,label,'completed',round(time.monotonic()-start,2),'s',flush=True)
                except Exception as exc:
                    failures.append(dict(catalyst=cat,case=label,error=repr(exc)));print('FAILED',failures[-1],flush=True)
            if 'tol_0.01' in cache:
                for label in ['tol_1','tol_0.1','Radau_0.1']:
                    if label in cache:comparisons.append(dict(catalyst=cat,category='solver',case=label,reference='tol_0.01',**difference(cache[label][0],cache['tol_0.01'][0])))
            if 'seed_1e-12' in cache:
                for label in ['tol_0.1','seed_1e-8','seed_1e-10']:
                    if label in cache:comparisons.append(dict(catalyst=cat,category='seed',case=label,reference='seed_1e-12',**difference(cache[label][0],cache['seed_1e-12'][0])))
            if 'tol_1' in cache:
                res,sol=cache['tol_1'];true_min=dense_minimum(sol)[0]
                for n,n2 in [(400,2000),(800,4000),(1600,8000)]:
                    # Re-run unchanged production integration with refined output counts,
                    # including P-R2 (whose output count is otherwise fixed at 2000).
                    from unittest.mock import patch
                    original_integrate=h.integrate
                    def sampled_integrate(rhs,state0,L,args,n=400):
                        count=n2 if rhs is h.rhs_pr2 else n
                        return original_integrate(rhs,state0,L,args,n=count)
                    with patch.object(h,'integrate',sampled_integrate):
                        sampled=h.run(cat,n=n)
                    species_error=0.;temperature_error=0.
                    for section in ['PBR','PR1','PR2']:
                        z_values,Y_values=sampled[section]
                        reference=sol[section][0].sol(z_values*h.area(h.SMR_REACTOR))
                        species_error=max(species_error,float(abs(h.mole_fractions(Y_values)-h.mole_fractions(reference)).max()))
                        temperature_error=max(temperature_error,float(abs(Y_values[6]-reference[6]).max()))
                    min_error=float(sampled['PBR'][1][6].min()-true_min)
                    record=dict(catalyst=cat,category='output_grid',case=str(n)+'/'+str(n2),reference='dense log-inlet grid',max_species_difference=species_error,max_T_difference_K=temperature_error,conversion_difference=abs(h.summarize(sampled)['X_CH4']-h.summarize(res)['X_CH4']),passed=bool(species_error<=2e-5 and temperature_error<=.1 and min_error<=.1))
                    record['sampled_minimum_error_K']=min_error
                    comparisons.append(record)
                # Persist original h.run baseline comparison; grid sampling itself does not change the adaptive solve.
                old=np.load(ROOT/f'phase1_runs/2026-11-01-v3-reproduction/profiles_{cat}.npz')
                actual=h.run(cat)
                manifest.setdefault('baseline_exact_regression',{})[cat]=all(np.array_equal(actual[s][1],old[s+'_state']) for s in ['PBR','PR1','PR2'])
            if 'energy_tol_0.001' in cache:
                coarse=energy_terms(cat,cache['energy_tol_0.001'][1],8000)
                fine=energy_terms(cat,cache['energy_tol_0.001'][1],16000)
                for x,y in zip(coarse,fine):
                    y['quadrature_relative_change']=max(abs(y['wall_heat_W']-x['wall_heat_W']),abs(y['reaction_heat_W']-x['reaction_heat_W']),abs(y['sensible_line_integral_W']-x['sensible_line_integral_W']))/y['normalization_W']
                    y['passed']=bool(max(y['total_relative_residual'],y['thermal_relative_residual'])<=1e-4 and y['quadrature_relative_change']<=1e-5)
                    energy.append(y)
            # Save partial evidence after each catalyst, including any failures.
            
            for row in comparisons:row.setdefault('sampled_minimum_error_K',None)
            write_csv(out/'comparison_results.csv',comparisons);write_csv(out/'minimum_temperatures.csv',minima)
            write_csv(out/'energy_balance.csv',energy)
            (out/'state_checks.json').write_text(json.dumps(sanity,indent=2)+'\n');(out/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
        manifest['verification_status']='PASS' if verification_passed(comparisons,energy,sanity,failures) else 'FAIL'
        print('Numerical verification status:',manifest['verification_status'],flush=True)
    return verification_exit_code(manifest)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);raise SystemExit(main(ap.parse_args().out))
