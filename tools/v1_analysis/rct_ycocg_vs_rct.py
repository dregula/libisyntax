import random, itertools
fix=lambda x: x//2 if x>=0 else -((-x)//2)
fl=lambda x: x>>1
def fwd(R,G,B,h):
    Co=R-B; t=B+h(Co); Cg=G-t; return t+h(Cg),Co,Cg
def inv(Y,Co,Cg,h):
    t=Y-h(Cg); G=t+Cg; B=t-h(Co); return B+Co,G,B
vals=list(range(0,256,3))+[255]
cube=list(itertools.product(vals,vals,vals)); n=len(cube)
print('RGB sample points:',n)
for name,h in (('forward = fix (Philips compression samples; scalar C inverse)',fix),('forward = floor >>1 (Philips viewing samples; SIMD C inverse)',fl)):
    ok=all(inv(*fwd(R,G,B,h),h)==(R,G,B) for R,G,B in cube)
    print('\n##',name,'| same-rounding inverse exact:',ok)
    for qn,q in (('q=fix',fix),('q=floor',fl)):
        eY=eCb=eCr=0; mY=mCb=mCr=0; dv=set()
        for R,G,B in cube:
            Y,Co,Cg=fwd(R,G,B,h)
            Yr=(R+2*G+B)>>2; Cb=B-G; Cr=R-G
            cb=-q(Co)-Cg; cr=cb+Co
            if Y!=Yr: eY+=1; dv.add(Yr-Y)
            if cb!=Cb: eCb+=1; mCb=max(mCb,abs(cb-Cb))
            if cr!=Cr: eCr+=1; mCr=max(mCr,abs(cr-Cr))
        print("  direct matrix, %s: Y'=Y wrong %.1f%% (Y'-Y in %s) | Cb wrong %.1f%% (max %d) | Cr wrong %.1f%% (max %d)"%(qn,100*eY/n,sorted(dv),100*eCb/n,mCb,100*eCr/n,mCr))
random.seed(1); N=300000; d=0; mx=0
for _ in range(N):
    Y=random.randint(50,255);Co=random.randint(-60,60);Cg=random.randint(-60,60)
    a=inv(Y,Co,Cg,fix); b=inv(Y,Co,Cg,fl)
    if a!=b: d+=1; mx=max(mx,max(abs(x-y) for x,y in zip(a,b)))
print('\n## arbitrary decoded-like triples: trunc vs floor inverse differ in %.1f%% of pixels, max|dRGB| %d'%(100*d/N,mx))
