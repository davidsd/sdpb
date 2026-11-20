"""SDPviBe's Cholesky of Q (one N x N SPD matrix, as its solver does at nmax14: F.Chol(copy_lower_to_upper(Q)),
or BlockChol when Q is large), timed on a random SPD matrix of the same size, for a matched comparison with
sdpb_gpu::Block_Factor at batch 1. SDPviBe by Antoine Tilloy (MIT), used unchanged.
  python sdpvibe_cholq_bench.py --sdpvibe <clone> --N 3406 [--L 28] [--reps 2] [--block 0|2048]"""
import argparse, sys, time
ap = argparse.ArgumentParser()
ap.add_argument("--sdpvibe", required=True)
ap.add_argument("--N", type=int, default=3406)
ap.add_argument("--L", type=int, default=28)
ap.add_argument("--reps", type=int, default=2)
ap.add_argument("--block", type=int, default=0)
a = ap.parse_args()
sys.path.insert(0, a.sdpvibe)
import torch
from sdpvibe import mpf as F
from sdpvibe.solver import _rand_mpf
gen = torch.Generator(device="cuda"); gen.manual_seed(3)
M = _rand_mpf((1, a.N, a.N), a.L, gen)
Q = F.gemm(M, M, transB=True)
Q = F.MPF(Q.m[:, 0].contiguous(), Q.e[0].contiguous())
print(f"device {torch.cuda.get_device_name(0)}; SDPviBe Cholesky of an N={a.N} SPD matrix, L={a.L}", flush=True)
def run():
    if a.block:
        from sdpvibe.blockchol import BlockChol
        return BlockChol([Q.clone()], a.block)
    return F.Chol(Q)
run(); torch.cuda.synchronize()
best = 1e30
for r in range(a.reps):
    torch.cuda.synchronize(); t0 = time.time(); c = run(); torch.cuda.synchronize(); t = time.time() - t0
    best = min(best, t); print(f"  rep {r}: {t:.3f} s, peak {torch.cuda.max_memory_allocated() / 1e9:.1f} GB", flush=True)
    del c
print(f"SDPviBe Chol(Q) N={a.N}{' BlockChol ' + str(a.block) if a.block else ''}: best {best:.3f} s")
