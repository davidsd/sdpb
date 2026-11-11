"""SDPviBe's GramAcc (sdpvibe/mpf.py, Antoine Tilloy, MIT), timed on the shapes of
gram_check, for a matched comparison with sdpb_gpu::Gram_Accumulator on the same GPU.

  python sdpvibe_gram_bench.py --sdpvibe <clone at 0d02884> --N 3406 --rows 17520 --chunks 4 [--L 28] [--reps 3]

As in SDPviBe's solver (lean mode): Q accumulated chunk by chunk with GramAcc.add, one result() (CRT) at
the end; the reference exponents come from a first pass (the solver carries them between iterations, so
no flush happens in steady state).  One random chunk (solver._rand_mpf) is reused for every chunk.
"""
import argparse, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--sdpvibe", required=True)
ap.add_argument("--N", type=int, required=True)
ap.add_argument("--rows", type=int, required=True)
ap.add_argument("--chunks", type=int, default=1)
ap.add_argument("--L", type=int, default=28)          # 832 bits: ceil((832 + 2) / 32) + 1
ap.add_argument("--reps", type=int, default=3)
a = ap.parse_args()
sys.path.insert(0, a.sdpvibe)
import torch
from sdpvibe import mpf as F
from sdpvibe.solver import _rand_mpf

gen = torch.Generator(device="cuda"); gen.manual_seed(0)
P = _rand_mpf((a.rows, a.N), a.L, gen)
torch.cuda.synchronize()
print(f"device {torch.cuda.get_device_name(0)}; SDPviBe GramAcc, L={a.L} limbs, N={a.N}, rows/chunk={a.rows}, "
      f"chunks={a.chunks}; RNS primes {F.ctx(a.L + 1).n} (4 INT8 products each)", flush=True)
t0 = time.time()
acc = F.GramAcc(a.N, a.L)
for c in range(a.chunks): acc.add(P)
q = acc.result(); torch.cuda.synchronize()
print(f"first pass (kernel compilation, reference exponents): {time.time() - t0:.2f} s, flushes {acc.flushes}", flush=True)
ref = acc.next_reference()
del q, acc
best = 1e30
for rep in range(a.reps):
    torch.cuda.synchronize(); t0 = time.time()
    acc = F.GramAcc(a.N, a.L, ref)
    for c in range(a.chunks): acc.add(P)
    t_add = time.time()
    torch.cuda.synchronize(); t_add = time.time() - t0
    q = acc.result(); torch.cuda.synchronize()
    t = time.time() - t0
    best = min(best, t)
    print(f"rep {rep}: {t:.4f} s total (add {t_add:.4f} s, CRT {t - t_add:.4f} s), flushes {acc.flushes}, "
          f"peak {torch.cuda.max_memory_allocated() / 1e9:.1f} GB", flush=True)
    del q, acc
mpmac = a.rows * a.chunks * a.N * (a.N + 1) / 2
print(f"best {best:.4f} s; {mpmac:.3g} multiprecision multiply-adds (lower triangle) -> {mpmac / best:.3g} MP-MAC/s")
