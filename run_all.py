"""run_all.py - runs the complete workflow in the right order (several minutes, mostly the contour stage).
   1) hyd.py              model run, results plots, summary csv
   2) validation.py       Fig. 4 overlay, error table, equilibrium check
   3) compare_h_length.py tube- vs particle-diameter comparison
   4) hyd_contours.py     operating-window and axial contour figures
   All outputs go to ./results (change with --out)."""
import argparse, os, subprocess, sys

ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results")
ap.add_argument("--workers", type=int, default=4, help="parallel workers of the contour stage")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
for script, extra in (("hyd.py", ["--seed-study"]), ("validation.py", []),
                      ("compare_h_length.py", []), ("hyd_contours.py", ["--workers", str(a.workers)])):
    print(f"\n=== {script} ===", flush=True)
    subprocess.run([sys.executable, os.path.join(here, script), "--out", a.out] + extra, check=True, cwd=here)
print("\nAll done. Results in", os.path.abspath(a.out))
