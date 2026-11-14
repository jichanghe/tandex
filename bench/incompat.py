import bisect, sys, collections
sys.argv=[sys.argv[0]]+sys.argv[1:]
src=open('catalog_eval.py').read().split("bins = [")[0]; exec(src)
C=load(sys.argv[1],10); st=[c[0] for c in C]; mx=max(c[1]-c[0] for c in C)
R=load('chr21.trf.bed',5); rs=[c[0] for c in R]; rmx=max(c[1]-c[0] for c in R)
lo,hi=int(sys.argv[2]),int(sys.argv[3]); n=0; kinds=collections.Counter()
for s,e,m in cat:
    if not lo<=e-s<=hi: continue
    def hits(C,st,mx):
        i=bisect.bisect_left(st,s-mx); o=[]
        while i<len(C) and C[i][0]<e:
            if min(e,C[i][1])-max(s,C[i][0])>=0.5*(e-s): o.append(C[i])
            i+=1
        return o
    h=hits(C,st,mx)
    if not h or any(c[2]==len(m) or compatible(m,c[3]) for c in h): continue
    t=hits(R,rs,rmx); tok=any(c[2]==len(m) or compatible(m,c[3]) for c in t)
    ps=sorted(c[2] for c in h); k='shorter' if ps[-1]<len(m) else 'longer' if ps[0]>len(m) else 'mixed'
    kinds[(k,'TRF ok' if tok else 'TRF also off')]+=1
    n+=1
    if n<=int(sys.argv[4]): print(f"cat {s}-{e} len={e-s} motif({len(m)})={m[:40]} | tandex {[(c[0],c[1],c[2],c[3][:24]) for c in h][:2]} | TRF periods {sorted(set(c[2] for c in t))[:6]}")
print(n,kinds)
