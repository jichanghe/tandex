import sys, time, pytrf
name=None; seq=[]
for l in open(sys.argv[1]):
    if l[0]=='>': name=l[1:].split()[0]
    else: seq.append(l.strip())
seq=''.join(seq).upper()
t=time.time(); out=open(sys.argv[2],'w')
for r in pytrf.ATRFinder(name, seq, min_motif=1, max_motif=int(sys.argv[3]), min_identity=70):
    out.write(f"{name}\t{r.start-1}\t{r.end}\t{len(r.motif)}\n")
for r in pytrf.GTRFinder(name, seq, min_motif=1, max_motif=int(sys.argv[3]), min_repeat=2, min_length=25):
    out.write(f"{name}\t{r.start-1}\t{r.end}\t{len(r.motif)}\n")
print('pytrf seconds: %.1f'%(time.time()-t), file=sys.stderr)
