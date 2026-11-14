#!/usr/bin/env python3
"""Tool-independent mask of the 178-bp centromere satellite (CEN180) in Col-CEN:
1-kb windows in which >=60% of bases equal the base 176..180 bp downstream (best lag).
Unrelated sequence scores about 0.3."""
import numpy as np, sys
name=None; chunks=[]; seqs={}
for l in open('col.fa'):
    if l[0]=='>':
        if name: seqs[name]=''.join(chunks)
        name=l[1:].split()[0]; chunks=[]
    else: chunks.append(l.strip())
seqs[name]=''.join(chunks)
W=1000; out=open('cen180.mask.bed','w'); tot=0
for n,s in seqs.items():
    a=np.frombuffer(s.upper().encode(),dtype=np.uint8); best=None
    for lag in range(176,181):
        m=(a[:-lag]==a[lag:]).astype(np.int32); nw=len(m)//W
        f=m[:nw*W].reshape(nw,W).mean(axis=1)
        best=f if best is None else np.maximum(best[:len(f)],f[:len(best)])
    on=best>=0.6; i=0; c=0
    while i<len(on):
        if on[i]:
            j=i
            while j<len(on) and on[j]: j+=1
            out.write(f"{n}\t{i*W}\t{j*W}\n"); c+=(j-i)*W; i=j
        else: i+=1
    tot+=c; print(n,len(s),'CEN180-like bp:',c, 'background median match frac %.2f'%np.median(best))
print('total',tot)
