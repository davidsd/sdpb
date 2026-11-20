"""Independent references and gates for exact SPD fixtures.
Small fixtures only; capped at one CPU. Run larger tests on a compute node.
"""
import argparse, json, os, pathlib, re, subprocess, time
p=argparse.ArgumentParser();p.add_argument('data');p.add_argument('results');a=p.parse_args()
root=pathlib.Path(a.data);out=pathlib.Path(a.results);out.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS='1');summary=[];start=time.monotonic()
for d in sorted((root/'synthetic').glob('n*-row*')):
 n=json.load(open(d/'meta.json'))['Q_N']
 for prec in [1024,1280]:
  cache=d/f'ref{prec}.txt'
  for mode in ['full','adaptive','original']:
   args=[str(root/'ref_solve'),str(d),'0',str(d/f'{mode}.cols.bin'),str(prec),
         '--save-ref' if mode=='full' else '--load-ref',str(cache)]
   if mode!='original':args+=['--require-no-worse-sv','1']
   s=subprocess.check_output(args,env=env,text=True)
   (out/f'{d.name}-{mode}-mpfr{prec}.log').write_text(s)
   cols=[tuple(map(float,m)) for m in re.findall(r'^COL \d+ \d+ \d+ ([-0-9.]+) ([-0-9.]+)',s,re.M)]
   assert len(cols)==32
   if mode!='original':assert all(b<=v for b,v in cols),(d,mode,cols)
   summary.append(dict(fixture=d.name,precision=prec,mode=mode,backend=max(b for b,v in cols),sdpvibe=max(v for b,v in cols)))
 s=subprocess.check_output([str(root/'compare_references'),str(n),str(d/'ref1024.txt'),str(d/'ref1280.txt'),'-580'],text=True)
 (out/f'{d.name}-convergence.log').write_text(s)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));print('All corrected columns no worse than SDPviBe; precision-convergence gate passed. Seconds:',time.monotonic()-start)
