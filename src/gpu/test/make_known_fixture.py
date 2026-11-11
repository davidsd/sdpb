"""Exact triangular known answer, including a zero correction.
L has diagonal 2 and every strict-lower entry -1. Both L and L^-1 are dyadic,
so Q=LL^T, B=I and P=L^-1 fit the stored precision exactly at n=64.
P0 is an analytic answer in this fixture, not an SDPviBe computation.
"""
import argparse,array,json,pathlib
from fractions import Fraction
p=argparse.ArgumentParser();p.add_argument('out');p.add_argument('--n',type=int,default=64);a=p.parse_args()
n=a.n;out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
L=28;ez=-(1<<40)
T=[[Fraction(2 if i==j else -1 if i>j else 0) for j in range(n)] for i in range(n)]
X=[[Fraction(0) for j in range(n)] for i in range(n)]
for j in range(n):
 for i in range(j,n):X[i][j]=(int(i==j)-sum(T[i][k]*X[k][j] for k in range(j,i)))/T[i][i]
Q=[[sum(T[i][k]*T[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
I=[[Fraction(i==j) for j in range(n)] for i in range(n)]
def write(name,A,column_major=False):
 m=array.array('I',[0])*(L*n*n);e=array.array('q',[ez])*(n*n)
 for i in range(n):
  for j in range(n):
   x=A[i][j];q=j*n+i if column_major else i*n+j
   if not x:continue
   den=x.denominator;assert den&(den-1)==0
   sh=894-abs(x.numerator).bit_length();assert sh>=0
   z=(x.numerator<<sh)%(1<<896);e[q]=-(den.bit_length()-1)-sh
   for k in range(L):m[k*n*n+q]=(z>>(32*k))&0xffffffff
 with open(out/(name+'.m.bin'),'wb') as f:m.tofile(f)
 with open(out/(name+'.e.bin'),'wb') as f:e.tofile(f)
for name,A in [('S0',Q),('B0',I),('P0',X)]:write(name,A)
write('known.L',T,True);write('known.X',X,True)
json.dump(dict(L=L,N=n,Q_N=n,reference='exact dyadic analytic triangular inverse'),open(out/'meta.json','w'),indent=2)
print('Exact known-answer fixture',out,'dimension',n)
