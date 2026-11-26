import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import *
path,csvp=sys.argv[1],sys.argv[2]
g,rows=open_file(path,csvp); top,pad=g['top'],g['pad']; Hs,LLs={}, {}
for c in range(3): Hs[c],LLs[c]=decode_coeffs(path,rows,g,c)
decs=[Windowed(Hs[c],LLs[c],top,g['Hp'],g['Wp']) for c in range(3)]
D=5; F=512
res={}
for name,i0y,i0x in (('tissue',64*157+pad,64*237+pad),('tissue edge',64*137+pad,64*137+pad),('glass',64*37+pad,64*37+pad),('scan edge',64*197+pad,64*385+pad)):
    pix=[d.out_region(0,i0y,i0y+F,i0x,i0x+F) for d in decs]
    py=[full_pyramid(p,D) for p in pix]
    R,G,B=inv_rgb(pix[0],pix[1],pix[2],fl); Yr,Cbr,Crr=rct_fwd(R,G,B)
    (LLy,lvy),cbw,crw=convert_pyramid(py,fl)
    Cbp=inv_pyramid(*cbw); Crp=inv_pyramid(*crw); Yp=pix[0]
    eb=Cbp-Cbr; er=Crp-Crr
    res[name]=(round(eb.mean()),round(er.mean()))
    print('%-12s D=5 F=%d: raw bias Cb %+.2f Cr %+.2f' % (name,F,eb.mean(),er.mean()))
# one global constant per component (median over frames), applied by adding to LL_D only (exact DC shift of the whole component)
Kb=int(np.median([v[0] for v in res.values()])); Kr=int(np.median([v[1] for v in res.values()]))
print('global DC correction chosen from the frames: Cb %+d, Cr %+d (added to the LL_%d band => exact constant shift of the decoded plane)'%(-Kb,-Kr,D))
for name,i0y,i0x in (('tissue',64*157+pad,64*237+pad),('tissue edge',64*137+pad,64*137+pad),('glass',64*37+pad,64*37+pad),('scan edge',64*197+pad,64*385+pad)):
    pix=[d.out_region(0,i0y,i0y+F,i0x,i0x+F) for d in decs]; py=[full_pyramid(p,D) for p in pix]
    R,G,B=inv_rgb(pix[0],pix[1],pix[2],fl); Yr,Cbr,Crr=rct_fwd(R,G,B)
    (LLy,lvy),(LLb,lvb),(LLr,lvr)=convert_pyramid(py,fl)
    Cbp=inv_pyramid(LLb-Kb,lvb); Crp=inv_pyramid(LLr-Kr,lvr); RGB=rct_inv(pix[0],Cbp,Crp)
    eb=Cbp-Cbr; er=Crp-Crr; dR=np.stack([a-b for a,b in zip(RGB,(R,G,B))])
    print('%-12s after DC fix: Cb err mean %+.2f std %.2f range [%d,%d] exact %.1f%% | Cr err mean %+.2f std %.2f range [%d,%d] | RGB |d|<=1: %.1f%%, <=2: %.1f%%, max %d' % (
        name,eb.mean(),eb.std(),eb.min(),eb.max(),100*np.mean(eb==0),er.mean(),er.std(),er.min(),er.max(),100*np.mean(np.abs(dR).max(0)<=1),100*np.mean(np.abs(dR).max(0)<=2),np.abs(dR).max()))
