"""Run the default workflow into isolated stage directories with execution manifests."""
import argparse
from pathlib import Path
import subprocess
import sys
from run_metadata import prepare_output, record_run


def main(out=None, workers=4):
    out = prepare_output(out, 'workflow')
    here = Path(__file__).resolve().parent
    stages = [('hyd.py', ['--seed-study'], 'model'), ('validation.py', [], 'fixed_reference'),
              ('compare_h_length.py', [], 'h_length_sensitivity'),
              ('hyd_contours.py', ['--workers', str(workers)], 'contours')]
    with record_run(out, 'workflow', dict(configuration='v3_reproduction', workers=workers)) as record:
        record['commands'] = []
        for script, extra, folder in stages:
            command = [sys.executable, str(here/script), '--out', str(out/folder), *extra]
            print(f'Running {script}', flush=True)
            result = subprocess.run(command, cwd=here)
            record['commands'].append(dict(argv=command, returncode=result.returncode))
            result.check_returncode()
    print('All done. Results in', out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    main(args.out, args.workers)
