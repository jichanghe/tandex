#!/usr/bin/env python3
"""cen_eval.py name=calls.tsv ... [region=Chr:start-end]: coverage of the CEN180 mask by each tool."""
import sys, bisect
from collections import defaultdict
region=None; files=[]
for a in sys.argv[1:]:
    k,v=a.split('=')
    if k=='region': c,r=v.split(':'); region=(c,int(r.split('-')[0]),int(r.split('-')[1]))
    else: files.append((k,v))
mask=defaultdict(list)
for l in open('cen180.mask.bed'):
    c,s,e=l.split(); s=int(s); e=int(e)
    if region:
        if c!=region[0]: continue
        s=max(s,region[1]); e=min(e,region[2])
        if e<=s: continue
    mask[c].append((s,e))
mbp=sum(e-s for v in mask.values() for s,e in v)
def union(iv):
    o=[]
    for s,e in sorted(iv):
        if o and s<=o[-1][1]: o[-1][1]=max(o[-1][1],e)
        else:o.append([s,e])
    return o
def inter(a,b):
    i=j=0;t=0
    while i<len(a) and j<len(b):
        s=max(a[i][0],b[j][0]);e=min(a[i][1],b[j][1])
        if e>s:t+=e-s
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return t
print(f"CEN180 mask: {mbp} bp" + (f" in {region}" if region else ""))
print("tool\tcalls in mask\tmask covered by any call\tmask covered by calls with period ~178 (176-180)\tby period = n x 178 (+-2%)\tlargest single call")
for name,f in files:
    allc=defaultdict(list); mono=defaultdict(list); mult=defaultdict(list); n=0; big=0
    for l in open(f):
        if l[0]=='#' or l.startswith('SeqID'): continue
        x=l.split('\t')
        try: c,s,e,p=x[0],int(x[1]),int(x[2]),int(float(x[3]))
        except ValueError: continue
        if c not in mask: continue
        allc[c].append((s,e))
        k=round(p/178)
        if 176<=p<=180: mono[c].append((s,e))
        if k>=1 and abs(p-k*178)<=0.02*p+1: mult[c].append((s,e))
    ca=cm=cx=0
    for c in mask:
        m=union(mask[c]); ca+=inter(m,union(allc[c])); cm+=inter(m,union(mono[c])); cx+=inter(m,union(mult[c]))
        for s,e in allc[c]:
            if inter(m,[[s,e]])>0.5*(e-s): n+=1; big=max(big,e-s)
    print(f"{name}\t{n}\t{100*ca/mbp:.1f}%\t{100*cm/mbp:.1f}%\t{100*cx/mbp:.1f}%\t{big}")
