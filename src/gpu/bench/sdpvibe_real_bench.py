"""SDPviBe (Antoine Tilloy, MIT; used unchanged, as a library) on REAL Schur blocks dumped by
sdpvibe_dump_schur.py: its batched Cholesky F.Chol(S) and its L^-1 B (solve_panels), timed on the same GPU as
the backend's real_block_check, so that the two are matched on the same data.
  python sdpvibe_real_bench.py --sdpvibe <clone> --dir <dump> [--chunk 0] [--panels 4] [--reps 2]
"""
import argparse, json, os, sys, time
ap = argparse.ArgumentParser()
ap.add_argument("--sdpvibe", required=True)
ap.add_argument("--dir", required=True)
ap.add_argument("--chunk", type=int, default=0)
ap.add_argument("--panels", type=int, default=4)
ap.add_argument("--reps", type=int, default=2)
a = ap.parse_args()
sys.path.insert(0, a.sdpvibe)
import numpy as np
import torch
from sdpvibe import mpf as F
from sdpvibe.solver import solve_panels

meta = json.load(open(os.path.join(a.dir, "meta.json")))
L, ch = meta["L"], meta["chunks"][a.chunk]

def load(name, shape):
    m = np.fromfile(os.path.join(a.dir, f"{name}{a.chunk}.m.bin"), dtype="<i4").reshape(L, *shape)
    e = np.fromfile(os.path.join(a.dir, f"{name}{a.chunk}.e.bin"), dtype="<i8").reshape(*shape)
    return F.MPF(torch.from_numpy(m).cuda(), torch.from_numpy(e).cuda())

S = load("S", ch["S"])
B = load("B", ch["B"])
print(f"device {torch.cuda.get_device_name(0)}; SDPviBe on real blocks {a.dir}: S {ch['S']}, B {ch['B']}, L={L}",
      flush=True)

def timed(name, f, reps):
    f(); torch.cuda.synchronize()                      # warm-up (kernel compilation)
    best = 1e30
    for r in range(reps):
        torch.cuda.synchronize(); t0 = time.time(); out = f(); torch.cuda.synchronize(); t = time.time() - t0
        best = min(best, t)
        print(f"  {name} rep {r}: {t:.4f} s, peak {torch.cuda.max_memory_allocated() / 1e9:.1f} GB", flush=True)
        del out
    return best

tc = timed("Chol(S)", lambda: F.Chol(S), a.reps)
cs = F.Chol(S)
ts = timed("solve_panels(L^-1 B)", lambda: solve_panels(cs, B, a.panels), a.reps)
print(f"SDPviBe real blocks: Chol {tc:.4f} s + L^-1 B {ts:.4f} s = {tc + ts:.4f} s", flush=True)
