"""Phase 2 controlled transport comparison; new directory required for every run."""
import argparse
import csv
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import hyd as h
import validation as v
from tools.reproduce_phase1 import compare

def save_json(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

# All columns are SMR columns on paper p.11; cp converted kJ/kg/K -> J/kg/K.
TABLE5=[dict(T_K=773.,rho=.0474,us=.200,mu=3.56e-5,cp=4520.,kf=.242,ks=1.99,
             Re_paper=1.600,Pr_paper=.665,Nu_paper=2.91,U_paper=584.),
        dict(T_K=1173.,rho=.0206,us=.459,mu=6.49e-5,cp=4710.,kf=.464,ks=3.39,
             Re_paper=.878,Pr_paper=.659,Nu_paper=2.83,U_paper=1019.)]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);a=ap.parse_args()
    out=Path(a.out).resolve()
    if out.exists() or out.is_relative_to(ROOT/'results') or out.is_relative_to(ROOT/'phase1_runs'):
        raise SystemExit('Use a new directory outside preserved baseline results.')
    out.mkdir(parents=True)
    os.environ['MPLCONFIGDIR']=str(out/'mplconfig');os.environ['MPLBACKEND']='Agg'
    meta=dict(status='RUNNING',started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
              command=sys.argv,python=sys.version,platform=platform.platform(),
              commit=git('rev-parse','HEAD').decode().strip(),cwd=str(Path.cwd()),
              source_hashes={str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'hyd.py',ROOT/'validation.py',ROOT/'compare_h_length.py',Path(__file__).resolve(),*sorted((ROOT/'data').glob('*.csv'))]},
              dependencies={n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib']},
              settings=dict(T_K=853,flow_ml_min_RTP=120,pressure_bar=h.P_BAR,seed=1e-6,solver=h.SOLVER,rtol=1e-8,atol_species=1e-16,atol_T=1e-6),
              configurations={name:asdict(h.Reactor(transport_config=name)) for name in h.TRANSPORT_CONFIGS},
              comparison_protocol_sha256=digest(ROOT/'docs/phase2/comparison_protocol.md'),
              commands=[])
    save_json(out/'manifest.json',meta)
    try:
        # Re-run existing interfaces without changing their numerical options.
        legacy=out/'legacy_regression';legacy.mkdir()
        for script in ['hyd.py','validation.py','compare_h_length.py']:
            cmd=[sys.executable,str(ROOT/script),'--out',str(legacy)]
            start=time.monotonic()
            with (legacy/(script+'.log')).open('w') as log:
                p=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            meta['commands'].append(dict(argv=cmd,returncode=p.returncode,seconds=time.monotonic()-start))
            if p.returncode:raise RuntimeError(script+' failed')
        comparisons=[]
        for name in ['smr_summary.csv','validation_table.csv','validation_equilibrium.csv','h_length_comparison.csv']:
            comparisons.append(compare(name,git('show','v3-baseline:results/'+name),(legacy/name).read_bytes()))
        save_json(out/'legacy_csv_regression.json',comparisons)
        exp,curves=v.load_experiment(),v.load_paper_curves()
        rows=[];transport=[];sanity=[];curve_rows=[];all_runs={}
        atoms=np.array([[1,0,1,1,0,0],[4,2,0,0,2,0],[0,1,1,2,0,0],[0,0,0,0,0,1]])
        for config in h.TRANSPORT_CONFIGS:
            reactor=h.Reactor(transport_config=config)
            for cat in v.CATS:
                res=h.run(cat,reactor=reactor);all_runs[config,cat]=res
                z,y,T,b=h.profiles(res);mask=z<=b['PBR_end']
                rms={}
                for s in v.SPECIES+['T']:
                    d=curves[s,cat];d=d[d[:,0]>(.15 if s=='T' else .05)]
                    pred=np.interp(d[:,0],z[mask],(T if s=='T' else y[v.IDX[s]])[mask])
                    rms[s]=float(np.sqrt(np.mean((pred-d[:,1])**2)))
                    curve_rows.append(dict(configuration=config,catalyst=cat,species=s,n=len(d),z_min_mm=d[:,0].min(),z_max_mm=d[:,0].max(),RMSE=rms[s],units='K' if s=='T' else 'mole_fraction'))
                d=curves['T',cat];d=d[(d[:,0]>.15)&(d[:,0]<5)]
                tr=float(np.sqrt(np.mean((np.interp(d[:,0],z[mask],T[mask])-d[:,1])**2)))
                outlet=v.exit_composition(res);e=np.array([outlet[s]-exp[cat][s][0] for s in v.SPECIES])
                mape=float(np.mean([100*abs(outlet[s]-exp[cat][s][0])/exp[cat][s][0] for s in v.SPECIES]))
                summary=h.summarize(res)
                row=dict(configuration=config,catalyst=cat,mean_species_RMSE=float(np.mean([rms[s] for s in v.SPECIES])),temperature_RMSE_K=rms['T'],temperature_RMSE_inlet_K=tr,experimental_MAE=float(abs(e).mean()),experimental_RMSE=float(np.sqrt(np.mean(e**2))),experimental_MAPE_pct=mape,X_CH4=summary['X_CH4'],PBR_Tmin_K=float(res['PBR'][1][6].min()),PBR_exit_H2=float(h.mole_fractions(res['PBR'][1])[4,-1]),outlet_CH4=outlet['CH4'],outlet_H2=outlet['H2'])
                rows.append(row)
                U,props=h.U_pbr(res['F0'],853,reactor)
                transport.append(dict(catalyst=cat,T_K=853.,**props))
                states=np.hstack([res[s][1] for s in ['PBR','PR1','PR2']]);initial=atoms@res['F0']
                errs=np.max(np.abs((atoms@states[:6]-initial[:,None])/initial[:,None]),axis=1)
                norm=float(abs(h.mole_fractions(states).sum(axis=0)-1).max());minimum=float(states[:6].min())
                good=bool(np.isfinite(states).all() and minimum>=-1e-14 and max(errs)<=1e-8 and norm<=1e-12)
                sanity.append(dict(configuration=config,catalyst=cat,passed=good,min_raw_flow_mol_s=minimum,C_error=errs[0],H_error=errs[1],O_error=errs[2],He_error=errs[3],normalization_error=norm))
                np.savez_compressed(out/f'profiles_{config}_{cat}.npz',F0=res['F0'],**{s+'_'+k:res[s][i] for s in ['PBR','PR1','PR2'] for i,k in enumerate(['z_m','state'])})
                print(config,cat,'completed',flush=True)
        # Differences relative to the same catalyst's preserved reproduction.
        base={r['catalyst']:r for r in rows if r['configuration']=='v3_reproduction'}
        for r in rows:
            for key in ['X_CH4','PBR_Tmin_K','outlet_CH4','outlet_H2']:
                r['delta_'+key]=r[key]-base[r['catalyst']][key]
        write_csv(out/'figure4_metrics.csv',rows);write_csv(out/'per_curve_metrics.csv',curve_rows)
        write_csv(out/'figure4_inlet_transport_SI.csv',transport);write_csv(out/'sanity_checks.csv',sanity)
        table=[]
        flow=24000*h.area(h.SMR_REACTOR)*h.SMR_REACTOR.L_pbr*1e6/60
        for pub in TABLE5:
            T=pub['T_K']
            for mode in ['tabulated_properties','correlated_properties']:
                if mode=='tabulated_properties':
                    values={k:pub[k] for k in ['rho','us','mu','cp','kf','ks']}
                else:
                    F=h.feed_moles(flow,[.05,.2,0,0,0,.75]);mu,kf,cp,rho=h.mix_props(F/F.sum(),T)
                    us=F.sum()*h.R*T/(h.P_PA*h.area(h.SMR_REACTOR))
                    values=dict(rho=rho,us=us,mu=mu,cp=cp,kf=kf,ks=h.k_solid(T))
                for config in h.TRANSPORT_CONFIGS:
                    U,props=h.transport_from_properties(h.Reactor(transport_config=config),**values)
                    table.append(dict(property_source=mode,T_K=T,correlated_flow_ml_min_RTP=flow if mode=='correlated_properties' else '',**props,**{k:pub[k] for k in ['Re_paper','Pr_paper','Nu_paper','U_paper']},U_difference_W_m2K=U-pub['U_paper'],U_difference_pct=100*(U/pub['U_paper']-1)))
        write_csv(out/'table5_comparison_SI.csv',table);save_json(out/'table5_published_SI.json',TABLE5)
        save_json(out/'units_and_domains.json',dict(Re='dimensionless',Pr='dimensionless',Nu='dimensionless',h='W/m2/K',U_W_m2K='W/m2/K',kr='W/m/K',kf='W/m/K',ks='W/m/K',cp='J/kg/K',rho='kg/m3',us='m/s',mu='Pa s',Re_length_m='m',Nu_length_m='m',species_domain='identical supplied PBR digitized points z>0.05 mm',temperature_domain='identical supplied PBR points z>0.15 mm; inlet metric also z<5 mm',experiment='five species, full 40 mm P-R2 outlet',table5='773/1173 K; correlated case GHSV 2.4e4/h inferred RTP flow; tabulated case uses listed properties directly'))
        # All four variants on common axes, no retuning or automatic winner.
        import matplotlib.pyplot as plt
        styles=['-','--','-.',':']
        for cat in v.CATS:
            fig,axs=plt.subplots(2,3,figsize=(14,8));axs=axs.ravel()
            for ax,s in zip(axs,v.SPECIES+['T']):
                d=curves[s,cat];ax.scatter(d[:,0],d[:,1],s=3,c='gray',label='paper PBR points')
                for config,ls in zip(h.TRANSPORT_CONFIGS,styles):
                    z,y,T,b=h.profiles(all_runs[config,cat]);ax.plot(z,T if s=='T' else y[v.IDX[s]],ls,label=config)
                ax.axvline(12,color='gray',lw=.6);ax.set_title(s);ax.set_xlabel('display coordinate [mm]; P-R1 omitted')
                ax.set_ylabel('K' if s=='T' else 'mole fraction')
                if s!='T':ax.plot([15],[exp[cat][s][0]],'ko',mfc='none',label='40 mm outlet experiment')
            handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=8)
            fig.suptitle(cat+': transport interpretations, 853 K / 120 mL/min')
            fig.tight_layout(rect=(0,.12,1,.96));fig.savefig(out/f'comparison_{cat}.png',dpi=130);plt.close(fig)
        if not all(r['passed'] and r['byte_equal'] for r in comparisons):raise RuntimeError('Legacy CSV regression failed')
        if not all(r['passed'] for r in sanity):raise RuntimeError('Physical sanity failure')
        meta['status']='PASS'
    except Exception as exc:
        meta['status']='FAIL';meta['failure']=repr(exc);raise
    finally:
        meta['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
        meta['output_hashes']={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and 'mplconfig' not in p.parts and p.name!='manifest.json'}
        save_json(out/'manifest.json',meta)
    return 0

if __name__=='__main__':raise SystemExit(main())
