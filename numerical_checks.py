"""Saved-state physical checks; never projects, clips or changes reactor states."""
import numpy as np


def verification_exit_code(manifest):
    """Execution completion alone must never turn a failed verification into success."""
    return 0 if manifest.get('verification_status') == 'PASS' else 1

ATOM_MATRIX = np.array([[1,0,1,1,0,0], [4,2,0,0,2,0],
                        [0,1,1,2,0,0], [0,0,0,0,0,1]])


def state_checks(res):
    initial=ATOM_MATRIX @ res['F0']
    states=np.hstack([res[s][1] for s in ['PBR','PR1','PR2']])
    F=states[:6];total=F.sum(axis=0)
    errors=np.max(np.abs((ATOM_MATRIX@F-initial[:,None])/initial[:,None]),axis=1)
    norm=float(np.max(np.abs((F/total).sum(axis=0)-1)))
    conv=1-F[0]/res['F0'][0]
    result=dict(finite=bool(np.isfinite(states).all()),min_flow_mol_s=float(F.min()),
                min_T_K=float(states[6].min()),max_T_K=float(states[6].max()),
                max_normalization_error=norm,relative_C_H_O_He_errors=errors.tolist(),
                min_CH4_conversion=float(conv.min()),max_CH4_conversion=float(conv.max()))
    result['passed']=bool(result['finite'] and result['min_flow_mol_s']>=-1e-14 and
                           result['min_T_K']>0 and np.all(total>0) and norm<=1e-12 and
                           max(errors)<=1e-8 and conv.min()>=-1e-8 and conv.max()<=1+1e-8)
    return result


def require_valid(res):
    result=state_checks(res)
    if not result['passed']:
        raise RuntimeError(f'Invalid raw reactor states: {result}')
    return result
