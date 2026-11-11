"""SDPviBe (Antoine Tilloy, MIT) on the nmax14 dim-5 Schur-block shape, for matched timings with
sdpb_gpu::Rns_Gemm on the same GPU:
  * its exact batched product F.gemm(X, B): (nb, P, P) x (nb, P, N), the operation Rns_Gemm performs;
  * its triangular solve L^-1 B as its solver runs it (Chol(S), then solve_panels(cs, B, 4));
  * its batched Cholesky of the S blocks.
  python sdpvibe_block_bench.py --sdpvibe <clone> --nb 36 --P 945 --N 3406 [--L 28] [--reps 2]
"""
import argparse, sys, time
ap = argparse.ArgumentParser()
ap.add_argument("--sdpvibe", required=True)
ap.add_argument("--nb", type=int, default=36)
ap.add_argument("--P", type=int, default=945)
ap.add_argument("--N", type=int, default=3406)
ap.add_argument("--L", type=int, default=28)
ap.add_argument("--reps", type=int, default=2)
ap.add_argument("--panels", type=int, default=4)
a = ap.parse_args()
sys.path.insert(0, a.sdpvibe)
import torch
from sdpvibe import mpf as F
from sdpvibe.solver import _rand_mpf, solve_panels

def timed(name, f, reps):
    f(); torch.cuda.synchronize()                      # warm-up (kernel compilation)
    best = 1e30
    for r in range(reps):
        torch.cuda.synchronize(); t0 = time.time(); out = f(); torch.cuda.synchronize(); t = time.time() - t0
        best = min(best, t); print(f"  {name} rep {r}: {t:.4f} s, peak {torch.cuda.max_memory_allocated() / 1e9:.1f} GB", flush=True)
        del out
    return best

gen = torch.Generator(device="cuda"); gen.manual_seed(1)
print(f"device {torch.cuda.get_device_name(0)}; SDPviBe, nb={a.nb} blocks P={a.P}, N={a.N}, L={a.L}", flush=True)
X = _rand_mpf((a.nb, a.P, a.P), a.L, gen)
B = _rand_mpf((a.nb, a.P, a.N), a.L, gen)
t = timed("gemm X B", lambda: F.gemm(X, B), a.reps)
mpmac = a.nb * a.P * a.P * a.N
print(f"SDPviBe exact batched product: best {t:.4f} s, {mpmac:.3g} MP-MAC -> {mpmac / t:.3g} MP-MAC/s", flush=True)
S = F.gemm(X, X, transB=True)                           # X X^T: SPD
del X; torch.cuda.empty_cache()
tc = timed("Chol(S)", lambda: F.Chol(S), a.reps)
print(f"SDPviBe batched Cholesky of S: best {tc:.4f} s", flush=True)
cs = F.Chol(S); del S
ts = timed("solve_panels(L^-1 B)", lambda: solve_panels(cs, B, a.panels), a.reps)
print(f"SDPviBe L^-1 B (solve_panels, {a.panels} panels): best {ts:.4f} s, {mpmac / 2:.3g} MP-MAC (triangular) -> "
      f"{mpmac / 2 / ts:.3g} MP-MAC/s", flush=True)
