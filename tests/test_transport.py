"""Equation, dimensional-scaling, legacy interface and full-profile regression checks."""
import dataclasses
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import hyd as h

# Independent SI fixture; expected choices intentionally not read from production map.
CHOICES = {
    'v3_reproduction': (.006, .006, 'gas'),
    'v3_h_particle_only': (.006, .0002, 'gas'),
    'particle_lengths_gas_pr': (.0002, .0002, 'gas'),
    'printed_equations': (.0002, .0002, 'bed'),
}

class TransportTests(unittest.TestCase):
    def test_equations_and_effective_parameters(self):
        for name, (lr, lh, pk) in CHOICES.items():
            with self.subTest(name=name), patch.object(h, 'k_radial_ZS', return_value=.4):
                r=h.Reactor(transport_config=name)
                U,p=h.transport_from_properties(r,rho=.1,us=.2,mu=3e-5,cp=4000,kf=.15,ks=2.)
                re=.1*.2*lr/3e-5
                pr=4000*3e-5/(.15 if pk=='gas' else .4)
                nu=2.67+.53*re**.77*pr**.53
                hh=nu*.4/lh
                self.assertAlmostEqual(U,1/(1/hh+.001/1.6),places=12)
                for k,v in [('Re',re),('Pr',pr),('Nu',nu),('h',hh),('Re_length_m',lr),('Nu_length_m',lh)]:
                    self.assertAlmostEqual(p[k],v,places=12)
                self.assertEqual(p['Pr_conductivity'],pk)
                self.assertEqual(p['configuration'],name)

    def test_dimensional_scaling_and_wall_resistance(self):
        # Scaling conductivity and viscosity by q preserves Re/Pr when rho also
        # scales by q; with zero wall resistance, h and U must scale by q.
        for name in CHOICES:
            r=h.Reactor(transport_config=name,s=0)
            U,p=h.transport_from_properties(r,rho=.1,us=.2,mu=3e-5,cp=4000,kf=.15,ks=2.)
            V,q=h.transport_from_properties(r,rho=.3,us=.2,mu=9e-5,cp=4000,kf=.45,ks=6.)
            for k in ['Re','Pr','Nu']:self.assertAlmostEqual(p[k],q[k],places=12)
            self.assertAlmostEqual(V/U,3.,places=12)
            Uw,_=h.transport_from_properties(dataclasses.replace(r,s=.001),rho=.1,us=.2,mu=3e-5,cp=4000,kf=.15,ks=2.)
            self.assertLess(Uw,U)

    def test_legacy_aliases_and_explicit_precedence(self):
        F=h.feed_moles(120,[.05,.2,0,0,0,.75])
        for old,new in [('tube','v3_reproduction'),('particle','v3_h_particle_only')]:
            a=h.U_pbr(F,853,h.Reactor(h_length=old))
            b=h.U_pbr(F,853,h.Reactor(transport_config=new))
            self.assertEqual(a,b)
        self.assertEqual(h.effective_transport(h.Reactor(h_length='particle',transport_config='v3_reproduction'))['configuration'],'v3_reproduction')
        for r in [h.Reactor(h_length='bad'),h.Reactor(transport_config='bad')]:
            with self.assertRaises(ValueError):h.effective_transport(r)

    def test_cli_choices_and_conflicts(self):
        for args in [['--transport-config','printed_equations','--h-length','tube'],
                     ['--transport-config','printed_equations','--validate']]:
            p=subprocess.run([sys.executable,str(ROOT/'hyd.py'),*args],capture_output=True,text=True)
            self.assertEqual(p.returncode,2)
        p=subprocess.run([sys.executable,str(ROOT/'hyd.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0)
        for name in CHOICES:self.assertIn(name,p.stdout)

    def test_default_profiles_and_csv_match_baseline(self):
        prior=ROOT/'phase1_runs/2026-11-01-v3-reproduction'
        results={}
        for cat in ['HMMC','57-4Q','25-4Q']:
            res=h.run(cat);results[cat]=res
            with np.load(prior/f'profiles_{cat}.npz') as old:
                np.testing.assert_allclose(res['F0'],old['F0'],rtol=1e-10,atol=1e-16)
                for sec in ['PBR','PR1','PR2']:
                    np.testing.assert_array_equal(res[sec][0],old[sec+'_z_m'])
                    np.testing.assert_allclose(res[sec][1][:6],old[sec+'_state'][:6],rtol=1e-10,atol=1e-16)
                    np.testing.assert_allclose(res[sec][1][6],old[sec+'_state'][6],rtol=0,atol=1e-7)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'summary.csv';h.write_summary(results,p)
            expected=subprocess.check_output(['git','show','v3-baseline:results/smr_summary.csv'],cwd=ROOT)
            self.assertEqual(p.read_bytes(),expected)

if __name__=='__main__':unittest.main()
