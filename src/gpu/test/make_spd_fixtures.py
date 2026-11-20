"""Exact dyadic SPD fixtures, no floating Gram formation.
M has independent integer entries and column scales spanning 2^200. A modular
rank check proves M nonsingular, hence Q=MM^T is SPD over the rationals.
Optional row scales reproduce the real Q's large diagonal range.
"""
import argparse, array, json, pathlib, random, time
p=argparse.ArgumentParser();p.add_argument('out');p.add_argument('--sizes',default='64,128,256');a=p.parse_args()
L=28;ez=-(1<<40);out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
def full_rank_mod(M,p=65521):
 A=[[x%p for x in row] for row in M];n=len(A)
 for j in range(n):
  k=next((i for i in range(j,n) if A[i][j]),None)
  if k is None:return False
  A[j],A[k]=A[k],A[j];iv=pow(A[j][j],-1,p)
  for i in range(j+1,n):
   f=A[i][j]*iv%p
   A[i][j:]=[(x-f*y)%p for x,y in zip(A[i][j:],A[j][j:])]
 return True
for n in map(int,a.sizes.split(',')):
 t=time.monotonic();rng=random.Random(325+n)
 M=[[rng.randrange(-(1<<30),1<<30) for j in range(n)] for i in range(n)]
 assert full_rank_mod(M),'regenerate with another seed'
 M=[[x<<(200-200*j//(n-1)) for j,x in enumerate(row)] for row in M]
 Q=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1):Q[i][j]=Q[j][i]=sum(x*y for x,y in zip(M[i],M[j]))
 for span in [0,100]:
  d=out/f'n{n}-row{span}';d.mkdir(exist_ok=True)
  m=array.array('I',[0])*(L*n*n);e=array.array('q',[ez])*(n*n)
  for i in range(n):
   for j in range(n):
    z=Q[i][j];q=i*n+j
    if not z:continue
    sh=894-abs(z).bit_length();assert sh>=0
    z=(z<<sh)%(1<<896);e[q]=-462-span*i//(n-1)-span*j//(n-1)-sh
    for k in range(L):m[k*n*n+q]=(z>>(32*k))&0xffffffff
  with open(d/'S0.m.bin','wb') as f:m.tofile(f)
  with open(d/'S0.e.bin','wb') as f:e.tofile(f)
  json.dump({'L':L,'N':32,'Q_N':n,'column_span':200,'row_span':span,'seed':325+n,'rank_mod_prime':65521,'exact_gram':True},open(d/'meta.json','w'),indent=2)
 print(n,'two fixtures, exact and modular nonsingularity checked, seconds',time.monotonic()-t,flush=True)
