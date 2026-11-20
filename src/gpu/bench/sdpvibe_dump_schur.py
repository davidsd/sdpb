"""Dump real Schur blocks with SDPviBe (Antoine Tilloy, MIT; used unchanged, as a library): load an SDP and a
checkpoint, run SDPviBe's iteration up to the Schur complement, and write, for the chunk(s) of blocks of
dimension --dim, the Schur blocks S (nb, P, P), the constraint rows B (nb, P, N) and SDPviBe's own
P = L^-1 B (Chol + solve_panels), as raw little-endian arrays in SDPviBe's MPF layout:
  <out>/<name>.m.bin  int32 (L, nb, R, C)     <out>/<name>.e.bin  int64 (nb, R, C)     <out>/meta.json
  python sdpvibe_dump_schur.py --sdpvibe <clone> --sdp <sdp dir> --ck <checkpoint dir> --out <dir> [--dim 5]
SDPVIBE_NO_CACHE=1 is set so that SDPviBe never writes its parse cache into the SDP directory.
"""
import argparse, json, os, sys, time
os.environ["SDPVIBE_NO_CACHE"] = "1"
ap = argparse.ArgumentParser()
ap.add_argument("--sdpvibe", required=True)
ap.add_argument("--sdp", required=True)
ap.add_argument("--ck", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dim", type=int, default=5)
ap.add_argument("--precision", type=int, default=832)
a = ap.parse_args()
sys.path.insert(0, a.sdpvibe)
import numpy as np
import torch
from sdpvibe import mpf as F
from sdpvibe import solver as SOL
from sdpvibe import checkpoint as CK, state as STATE

L = (a.precision + 2 + 31) // 32 + 1
t0 = time.time()
sdp = SOL.SDP(a.sdp, L)
print(f"SDP loaded in {time.time() - t0:.0f} s: N={sdp.N}, {sdp.nblocks} blocks; {sdp.read_info}", flush=True)
solver = SOL.Solver(sdp, L, precision=a.precision, maxComplementarity="1e130")
st = CK.load(a.ck, sdp, L)
assert st is not None, "no checkpoint"
STATE.set_state(solver, st)
print(f"checkpoint {st['source']}", flush=True)

cap = []
orig = SOL.Solver.schur
def schur(self, g, AXI, AY):
    S = orig(self, g, AXI, AY)
    if g.dim == a.dim:
        cap.append((g, S))
    return S
SOL.Solver.schur = schur
class Stop(Exception):
    pass
def q_block(self, Q):
    raise Stop()
SOL.Solver._q_block = q_block
try:
    solver.run(max_iter=1, verbose=False)
except Stop:
    pass
print(f"captured {len(cap)} chunk(s) of dim-{a.dim} Schur blocks", flush=True)
os.makedirs(a.out, exist_ok=True)

def save(name, A):
    A.m.contiguous().cpu().numpy().astype("<i4").tofile(os.path.join(a.out, name + ".m.bin"))
    A.e.contiguous().cpu().numpy().astype("<i8").tofile(os.path.join(a.out, name + ".e.bin"))

meta = {"L": L, "N": sdp.N, "chunks": []}
for c, (g, S) in enumerate(cap):
    gi = sdp.groups.index(g)
    B = sdp.store.get(("B", gi))
    t1 = time.time()
    cs = F.Chol(S)
    P = SOL.solve_panels(cs, B, sdp.layout.panels)
    torch.cuda.synchronize()
    print(f"chunk {c}: nb={g.nb} P={g.P} blocks {g.idx[:3]}..; SDPviBe Chol + solve {time.time() - t1:.2f} s", flush=True)
    save(f"S{c}", S); save(f"B{c}", B); save(f"P{c}", P)
    meta["chunks"].append({"nb": g.nb, "P": g.P, "blocks": g.idx, "S": list(S.e.shape), "B": list(B.e.shape)})
json.dump(meta, open(os.path.join(a.out, "meta.json"), "w"), indent=1)
print("done", flush=True)
