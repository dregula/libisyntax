import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import *
path,csvp=sys.argv[1],sys.argv[2]
g,rows=open_file(path,csvp); top,pad=g['top'],g['pad']; Hs,LLs={}, {}
for c in range(3): Hs[c],LLs[c]=decode_coeffs(path,rows,g,c)
decs=[Windowed(Hs[c],LLs[c],top,g['Hp'],g['Wp']) for c in range(3)]
i0y,i0x,F=64*157+pad,64*237+pad,256
pix=[d.out_region(0,i0y,i0y+F,i0x,i0x+F) for d in decs]
R,G,B=inv_rgb(pix[0],pix[1],pix[2],fl); Yr,Cbr,Crr=rct_fwd(R,G,B)
print('exact pixel-domain Cb range',Cbr.min(),Cbr.max())
for D in (1,2,3,4,5):
    py=[full_pyramid(p,D) for p in pix]
    for qn,q in (('floor',fl),):
        (LLy,lvy),cbp,crp=convert_pyramid(py,q); Cbp=inv_pyramid(*cbp); Crp=inv_pyramid(*crp)
        e=Cbp-Cbr; e2=Crp-Crr
        print('D=%d q=%s: Cb error mean %+.3f std %.3f min %d max %d | Cr error mean %+.3f std %.3f | err==0 in %.1f%%' % (D,qn,e.mean(),e.std(),e.min(),e.max(),e2.mean(),e2.std(),100*np.mean(e==0)))
# sanity: random tiny image through the same machinery
rng=np.random.default_rng(2); Rr=rng.integers(0,256,(64,64)); Gr=rng.integers(0,256,(64,64)); Br=rng.integers(0,256,(64,64))
def fwdc(R,G,B,h):
    Co=R-B; t=B+h(Co); Cg=G-t; return t+h(Cg),Co,Cg
Yy,Cc,Gg=fwdc(Rr,Gr,Br,fl); Yq,Cbq,Crq=rct_fwd(Rr,Gr,Br)
py=[full_pyramid(p,3) for p in (Yy,Cc,Gg)]; (a,b),cbp,crp=convert_pyramid(py,fl)
Cbp=inv_pyramid(*cbp); print('random image D=3: Cb error mean %+.3f std %.3f max %d' % ((Cbp-Cbq).mean(),(Cbp-Cbq).std(),np.abs(Cbp-Cbq).max()))
