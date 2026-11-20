"""Known integer factor and permuted-column cache control (one CPU, seconds)."""
import argparse, os, pathlib, struct, subprocess
p=argparse.ArgumentParser();p.add_argument('binary');p.add_argument('out');p.add_argument('--compare-binary');a=p.parse_args()
o=pathlib.Path(a.out);o.mkdir(parents=True,exist_ok=True)
n=8;L=28;K=4
T=[[((i+j)%3-1 if j<i else (2 if i==j else 0)) for j in range(n)] for i in range(n)]
Y=[[i+2*c-3 for c in range(K)] for i in range(n)]
S=[[sum(T[i][k]*T[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
B=[[sum(T[i][k]*Y[k][c] for k in range(n)) for c in range(K)] for i in range(n)]
def words(x):return [(x>>(32*k))&0xffffffff for k in range(L)]
(o/'S0.m.bin').write_bytes(b''.join(struct.pack('<I',words(S[i][j])[t]) for t in range(L) for i in range(n) for j in range(n)))
(o/'S0.e.bin').write_bytes(struct.pack('<'+'q'*(n*n),*([0]*(n*n))))
def cols(order,perturb=False):
 b=b'QCOL'+struct.pack('<5i',1,n,K,L,K)+struct.pack('<4i',*order)+bytes(4*K)
 for c in order:
  for which,mat in enumerate((B,Y,Y)):
   for i in range(n):
    x=mat[i][c]+int(perturb and which==1 and i==3 and c==2)
    b+=struct.pack('<28I',*words(x))+struct.pack('<q',0)
 return b
(o/'orig.bin').write_bytes(cols([0,1,2,3]));(o/'permuted.bin').write_bytes(cols([3,1,0,2]));(o/'wrong.bin').write_bytes(cols([3,1,0,2],True))
def run(name,*opts,expected=0):
 r=subprocess.run([a.binary,str(o),'0',str(o/name),'1280',*map(str,opts)],env=dict(os.environ,OMP_NUM_THREADS='1'),text=True,capture_output=True)
 assert r.returncode==expected,(name,r)
 print(r.stdout);return r.stdout
base=run('orig.bin','--save-ref',o/'ref.txt');other=run('permuted.bin','--load-ref',o/'ref.txt');wrong=run('wrong.bin','--load-ref',o/'ref.txt','--require-no-worse-sv',1,expected=1)
assert 'worst backend 2^-1000000000.00' in base and 'worst backend 2^-1000000000.00' in other
assert 'worst backend 2^-1000000000.00' not in wrong
print('CACHE, PERMUTATION, AND PERTURBATION CONTROLS PASSED')
if a.compare_binary:
 ref=(o/'ref.txt').read_text()
 cases={'empty.txt':('',2),'duplicate.txt':(ref+ref,2),
        'malformed.txt':(ref+'oops\n',2),'nonfinite.txt':('0\n@NaN@\n',2)}
 changed=ref.splitlines();changed[1]='0';cases['changed-reference.txt']=('\n'.join(changed)+'\n',1)
 cases['exact-reference.txt']=(ref,0)
 for name,(contents,expected) in cases.items():
  (o/name).write_text(contents)
  r=subprocess.run([a.compare_binary,str(n),str(o/'ref.txt'),str(o/name),'-1000'],text=True,capture_output=True)
  print(name,'exit',r.returncode,'expected',expected);print(r.stdout+r.stderr)
  assert r.returncode==expected,(name,r)
 print('PRECISION-CONVERGENCE POSITIVE AND NEGATIVE CONTROLS PASSED')
