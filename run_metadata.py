"""Output isolation and execution provenance; contains no reactor calculations."""
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def prepare_output(out, stage):
    """Resolve relative paths against the caller, refusing existing nonempty runs."""
    if out is None:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        path = ROOT / 'runs' / f'{stamp}-{stage}'
    else:
        path = Path(out).expanduser().resolve()
    protected = [ROOT/'results', ROOT/'phase1_runs', ROOT/'phase2_runs']
    if path == ROOT or any(path == p or path.is_relative_to(p) for p in protected):
        raise ValueError(f'Protected output location: {path}; use a new run directory')
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f'Output is not empty: {path}; use a new run directory')
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@contextmanager
def record_run(out, stage, settings):
    """Record execution completion/failure, not scientific verification status."""
    out = Path(out)
    manifest = out / f'manifest_{stage}.json'
    def git(*args):
        try:
            return subprocess.check_output(['git', *args], cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except (OSError, subprocess.CalledProcessError):
            return None
    sources = sorted([*ROOT.glob('*.py'), *ROOT.glob('tools/*.py'), *ROOT.glob('data/*.csv')])
    data = dict(stage=stage, status='RUNNING', started_utc=datetime.now(timezone.utc).isoformat(),
                argv=sys.argv, cwd=str(Path.cwd()), output_directory=str(out),
                python=sys.version, executable=sys.executable, platform=platform.platform(),
                dependencies={n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib']},
                git_commit=git('rev-parse','HEAD'), git_status=git('status','--porcelain'),
                settings=settings, source_data_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources})
    def write():
        tmp=manifest.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        tmp.replace(manifest)
    write()
    try:
        yield data
    except BaseException as exc:
        data['status']='FAILED';data['error']=f'{type(exc).__name__}: {exc}'
        raise
    else:
        data['status']='COMPLETED'
    finally:
        data['finished_utc']=datetime.now(timezone.utc).isoformat()
        data['output_sha256']={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*'))
                               if p.is_file() and p!=manifest and not p.name.endswith('.tmp')}
        write()
