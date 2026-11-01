import unittest
import numpy as np
import hyd as h
from numerical_checks import state_checks,require_valid

class RawStateTests(unittest.TestCase):
    def fixture(self):
        F=h.feed_moles(120,[.05,.2,0,0,0,.75]);Y=np.repeat(np.r_[F,853.][:,None],3,axis=1)
        return {'F0':F,**{s:(np.linspace(0,L,3),Y.copy()) for s,L in [('PBR',.012),('PR1',.05),('PR2',.04)]}}
    def test_constant_state_conservation(self):
        self.assertTrue(state_checks(self.fixture())['passed'])
    def test_reject_negative_or_nonfinite_states(self):
        for row,value in [(0,-1e-8),(6,float('nan')),(6,-1.)]:
            r=self.fixture();r['PBR'][1][row,1]=value
            with self.assertRaises(RuntimeError):require_valid(r)
    def test_reject_helium_and_carbon_leak(self):
        for row in [0,5]:
            r=self.fixture();r['PR2'][1][row,-1]*=1.01
            self.assertFalse(state_checks(r)['passed'])

class VerificationVerdictTests(unittest.TestCase):
    def test_failed_sampling_check_stays_failed(self):
        from tools.verify_numerics import verification_passed
        yes=[{'passed':True}];no=[{'passed':False}]
        self.assertFalse(verification_passed(no,yes,yes,[]))
        self.assertFalse(verification_passed(yes,yes,yes,[{'error':'solver failed'}]))
        self.assertFalse(verification_passed([],yes,yes,[]))
        self.assertTrue(verification_passed(yes,yes,yes,[]))

if __name__=='__main__':unittest.main()
