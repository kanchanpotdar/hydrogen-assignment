import ast
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import hyd as h
import validation as v
import run_all
from run_metadata import prepare_output, record_run, sha


class ReproducibilityTests(unittest.TestCase):
    def test_custom_validation_rejected_before_output(self):
        for args in [['--T','900'],['--flow','80'],['--catalyst','HMMC'],
                     ['--h-length','particle'],['--transport-config','printed_equations']]:
            with tempfile.TemporaryDirectory() as td:
                out=Path(td)/'not_created'
                p=subprocess.run([sys.executable,str(ROOT/'hyd.py'),'--validate','--out',str(out),*args],capture_output=True,text=True)
                self.assertEqual(p.returncode,2,p.stderr)
                self.assertFalse(out.exists())
                self.assertIn('fixed Figure 4 benchmark',p.stderr)

    def test_output_safety_and_caller_relative_paths(self):
        for out in [ROOT,ROOT/'results',ROOT/'results/new',ROOT/'phase1_runs/new',ROOT/'phase2_runs/new']:
            with self.assertRaises(ValueError):prepare_output(out,'test')
        with tempfile.TemporaryDirectory() as td:
            before=Path.cwd()
            try:
                os.chdir(td)
                out=prepare_output('relative','test')
                self.assertEqual(out,Path(td).resolve()/'relative')
                (out/'existing.txt').write_text('preserve')
                with self.assertRaises(ValueError):prepare_output(out,'test')
            finally:os.chdir(before)
        # Generate default paths without leaving test directories behind.
        a=prepare_output(None,'test');b=prepare_output(None,'test')
        self.assertNotEqual(a,b);a.rmdir();b.rmdir()

    def test_success_and_failure_manifests(self):
        with tempfile.TemporaryDirectory() as td:
            out=prepare_output(Path(td)/'success','test')
            with record_run(out,'test',{'configuration':'v3_reproduction'}):
                (out/'artifact.txt').write_text('result')
            data=json.loads((out/'manifest_test.json').read_text())
            self.assertEqual(data['status'],'COMPLETED')
            self.assertEqual(data['output_sha256']['artifact.txt'],sha(out/'artifact.txt'))
            self.assertIn('hyd.py',data['source_data_sha256'])
            self.assertIn('numpy',data['dependencies'])
            bad=prepare_output(Path(td)/'failure','test')
            with self.assertRaisesRegex(RuntimeError,'injected'):
                with record_run(bad,'test',{}):raise RuntimeError('injected')
            self.assertEqual(json.loads((bad/'manifest_test.json').read_text())['status'],'FAILED')

    def test_workflow_uses_absolute_stage_directories(self):
        with tempfile.TemporaryDirectory() as td, patch.object(run_all,'subprocess') as child:
            call=child.run
            call.return_value=subprocess.CompletedProcess([],0)
            run_all.main(Path(td)/'workflow',workers=2)
            self.assertEqual(call.call_count,4)
            paths=[]
            for args in call.call_args_list:
                argv=args.args[0];p=Path(argv[argv.index('--out')+1]);paths.append(p)
                self.assertTrue(p.is_absolute())
                self.assertEqual(args.kwargs['cwd'],ROOT)
            self.assertEqual(len(set(paths)),4)
            self.assertEqual(json.loads((Path(td)/'workflow/manifest_workflow.json').read_text())['status'],'COMPLETED')

    def test_display_alias_and_reference_data(self):
        self.assertEqual(h.PR2_DISPLAY_MM,3.)
        self.assertEqual(h.PR2_VALIDATION_MM,h.PR2_DISPLAY_MM)
        curves=v.load_paper_curves();d=curves['T','57-4Q']
        self.assertEqual(len(d),35);self.assertEqual(d[:,0].max(),1.7946)
        self.assertEqual(len(curves),18)
        exp=v.load_experiment()
        self.assertEqual(exp['HMMC']['H2'][1],.002)

    def test_plot_section_labels_and_heatmap_visibility(self):
        import matplotlib.pyplot as plt
        states=np.zeros((7,4));states[:6]=np.array([.01,.1,.01,.02,.1,.7])[:,None];states[6]=[853,820,850,853]
        res={'PBR':(np.linspace(0,.012,4),states),
             'PR1':(np.linspace(0,.050,4),states),
             'PR2':(np.linspace(0,.040,100),np.repeat(states[:,-1:],100,axis=1))}
        with tempfile.TemporaryDirectory() as td, patch.object(plt,'close'):
            h.plot_results({'57-4Q':res},str(Path(td)/'profiles.png'),with_pr1=True)
            fig=plt.gcf()
            labels=[t.get_text() for ax in fig.axes for t in ax.texts]
            self.assertIn('P-R1',labels);self.assertIn('P-R2',labels)
            h.plot_axial_heatmaps(res,str(Path(td)/'fields.png'),'57-4Q')
            fig=plt.gcf();fig.canvas.draw();renderer=fig.canvas.get_renderer()
            for ax in fig.axes:
                if ax.get_ylabel() in ['CH4','H2','CO','CO2','H2O','He','Temperature']:
                    self.assertGreaterEqual(ax.yaxis.label.get_window_extent(renderer).x0,0)
        plt.close('all')

    def test_contour_annotations_with_synthetic_fixture(self):
        import hyd_contours as c
        import matplotlib.pyplot as plt
        class FixturePool:
            def map(self,func,jobs,**kwargs):
                return [(.3+.0001*T,.12+.00001*q,.2+.0001*q) for cat,T,q in jobs]
        with tempfile.TemporaryDirectory() as td, patch.object(plt,'close'):
            c.fig_operating(FixturePool(),td)
            title=plt.gcf()._suptitle.get_text()
            self.assertIn('120 mL/min only',title)
            self.assertNotIn('Validated only at the red stars',title)
        plt.close('all')

if __name__=='__main__':unittest.main()
