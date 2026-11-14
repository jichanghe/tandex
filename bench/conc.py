import bisect,sys
def load(f,sc):
    return sorted((int(x[1]),int(x[2]),int(x[3]),float(x[sc])) for x in (l.split('\t') for l in open(f) if l[0]!='#'))
def union(iv):
    o=[]
    for s,e in sorted(iv):
        if o and s<=o[-1][1]: o[-1][1]=max(o[-1][1],e)
        else:o.append([s,e])
    return o
def cov(u,us,s,e):
    c=0;i=max(0,bisect.bisect_right(us,s)-1)
    while i<len(u) and u[i][0]<e:
        c+=max(0,min(e,u[i][1])-max(s,u[i][0]));i+=1
    return c
T=load('chr21.tandex.tsv',5); R=load('chr21.trf.bed',4)
ut=union([(s,e) for s,e,p,_ in T]); us=[u[0] for u in ut]; ur=union([(s,e) for s,e,p,_ in R]); ur_s=[u[0] for u in ur]
bt=sum(e-s for s,e in ut); br=sum(e-s for s,e in ur); sh=sum(cov(ut,us,s,e) for s,e in ur)
print('tandex: %d calls, %d bp | TRF: %d records (%d merged loci), %d bp | shared %d bp = %.1f%% of TRF bp'%(len(T),bt,len(R),len(ur),br,sh,100*sh/br))
Rm=[r for r in R if cov(ut,us,r[0],r[1])<0.5*(r[1]-r[0])]
print('TRF records not found by tandex: %d (%.1f%%); TRF score <60: %d, 60-99: %d, >=100: %d'%(len(Rm),100*len(Rm)/len(R),sum(r[3]<60 for r in Rm),sum(60<=r[3]<100 for r in Rm),sum(r[3]>=100 for r in Rm)))
Tm=[t for t in T if cov(ur,ur_s,t[0],t[1])<0.5*(t[1]-t[0])]
print('tandex calls not found by TRF: %d (%d bp); tandex score <50: %d, 50-99: %d, >=100: %d; period >2000: %d'%(len(Tm),sum(t[1]-t[0] for t in Tm),sum(t[3]<50 for t in Tm),sum(50<=t[3]<100 for t in Tm),sum(t[3]>=100 for t in Tm),sum(t[2]>2000 for t in Tm)))
if len(sys.argv)>1:
    for r in sorted(Rm,key=lambda r:-r[3])[:8]: print('   TRF-only',r)
    for t in sorted(Tm,key=lambda r:-r[3])[:8]: print('   tandex-only',t)
