"""unchanged SDPviBe on the frozen Q; full factor export.
SDPviBe by Antoine Tilloy (MIT) is the backend's prototype.
"""
import argparse, os, sys, time, json
p=argparse.ArgumentParser();p.add_argument('--sdpvibe',required=True);p.add_argument('--dir',required=True);p.add_argument('--out',required=True);p.add_argument('--reps',type=int,default=7);p.add_argument('--extended-dir');p.add_argument('--columns',type=int,default=128)
a=p.parse_args();sys.path.insert(0,a.sdpvibe)
import numpy as np
import torch
from sdpvibe import mpf as F, solver as S
from sdpvibe.limbops import bitlen_words
meta=json.load(open(a.dir+'/meta.json')); n=meta['Q_N']; L=meta['L']
m=torch.from_numpy(np.fromfile(a.dir+'/S0.m.bin',dtype='<i4').reshape(L,n,n)).cuda()
e=torch.from_numpy(np.fromfile(a.dir+'/S0.e.bin',dtype='<i8').reshape(n,n)).cuda()
Q=F.MPF(m,e); c=F.Chol(Q);torch.cuda.synchronize()
ts=[]
for r in range(a.reps):
 t=time.perf_counter();c=F.Chol(Q);torch.cuda.synchronize();ts.append(time.perf_counter()-t)
print('TIME_SDPVIBE',n,ts,'median',np.median(ts),flush=True)
idx=torch.arange(n,device='cuda'); d=c.Lt.limbs[:,0,idx,idx].contiguous()
# Normalize every diagonal limb before reading a float. _top inspects all limbs.
diag=F.normalize(d,torch.full((n,),c.Lt.exp,dtype=torch.int64,device='cuda'),L)
actual_top=F._top(diag)+c.k[0]
print('DIAG_SDPVIBE original_estimate',float(S.chol_log2_cond(c)[0]),'full_limb_min_top',int(actual_top.min()),'max_top',int(actual_top.max()),'range_estimate',2*int(actual_top.max()-actual_top.min()),flush=True)
# Normalize all fixed-point L entries, then undo D row scaling, keeping padding.
npad=c.npad; le=torch.full((npad*npad,),c.Lt.exp,dtype=torch.int64,device='cuda')
LL=F.normalize(c.Lt.limbs.reshape(L,-1),le,L).reshape(npad,npad)
k=torch.zeros(npad,dtype=torch.int64,device='cuda');k[:n]=c.k[0]
LL=LL.shift(k[:,None]);os.makedirs(os.path.dirname(a.out),exist_ok=True)
# backend layout is column-major, unlike the input row-major matrix
LL.m.transpose(-1,-2).contiguous().cpu().numpy().tofile(a.out+'.m.bin')
LL.e.T.contiguous().cpu().numpy().tofile(a.out+'.e.bin')
print('EXPORTED',a.out,npad,flush=True)

if a.extended_dir:
 out=a.extended_dir;os.makedirs(out,exist_ok=True)
 targets=[j if j<32 else 32+(n-33)*(j-32)//95 for j in range(a.columns)]
 one=F.from_mp([1],L,(1,)); bm=torch.zeros((L,n,a.columns),dtype=torch.int32,device='cuda');be=torch.full((n,a.columns),F.EZ,dtype=torch.int64,device='cuda')
 for j,i in enumerate(targets): bm[:,i,j]=one.m[:,0];be[i,j]=one.e[0]
 B=F.MPF(bm,be);P=c.solve_lower(B);torch.cuda.synchronize()
 for name,A in [('B0',B),('P0',P)]:
  A.m.contiguous().cpu().numpy().tofile(out+'/'+name+'.m.bin');A.e.contiguous().cpu().numpy().tofile(out+'/'+name+'.e.bin')
 for name in ['S0.m.bin','S0.e.bin']:
  path=out+'/'+name
  if not os.path.lexists(path):os.symlink(os.path.relpath(os.path.abspath(a.dir+'/'+name), out), path)
 json.dump({'N':a.columns,'Q_N':n,'L':L,'targets':targets},open(out+'/meta.json','w'),indent=2)
 print('EXTENDED_COLS',out,targets,flush=True)
