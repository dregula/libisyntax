import random
def fwd1(x):                      # reversible 5/3, even=low, whole-sample symmetric extension, len(x) even
    N=len(x); X=lambda i: x[-i] if i<0 else (x[2*(N-1)-i] if i>=N else x[i])
    d=[X(2*n+1)-((X(2*n)+X(2*n+2))>>1) for n in range(N//2)]
    D=lambda n: d[-n-1] if n<0 else d[n]          # d[-1]=d[0]
    s=[X(2*n)+((D(n-1)+D(n)+2)>>2) for n in range(N//2)]
    return s,d
def levels(x,J):
    out=[]; cur=x
    for _ in range(J):
        s,d=fwd1(cur); out.append((s,d)); cur=s
    return out
random.seed(7)
def sig(n): 
    v=128; o=[]
    for i in range(n): v=max(0,min(255,v+random.randint(-9,9))); o.append(v)
    return o
for F,J in ((256,5),(512,5),(256,3)):
    G=4096; x=sig(G); a=G//2 - F//2; a-=a%64     # frame start multiple of 64
    gl=levels(x,J); fl=levels(x[a:a+F],J)
    print('\nframe',F,'px, %d levels  (frame start %d)'%(J,a))
    tot=0;diff=0
    for j in range(J):
        sg,dg=gl[j]; sl,dl=fl[j]; k0=a>>(j+1); n=F>>(j+1)
        ms=[i for i in range(n) if sg[k0+i]!=sl[i]]; md=[i for i in range(n) if dg[k0+i]!=dl[i]]
        # report distance of mismatches from nearest edge
        def span(m): 
            L=[i for i in m if i<n//2]; R=[n-1-i for i in m if i>=n//2]
            return (max(L)+1 if L else 0, max(R)+1 if R else 0)
        print(' level %d: %3d coeffs/band | LOW mismatches %3d, widest from (left,right) edge %s | HIGH mismatches %3d, widest %s'%(j+1,n,len(ms),span(ms),len(md),span(md)))
