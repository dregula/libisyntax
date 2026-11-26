import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colour_ebcot_probe import *
path,csvp=sys.argv[1],sys.argv[2]
g,rows=open_file(path,csvp); top,pad=g['top'],g['pad']; Hs,LLs={}, {}
for c in range(3): Hs[c],LLs[c]=decode_coeffs(path,rows,g,c)
decs=[Windowed(Hs[c],LLs[c],top,g['Hp'],g['Wp']) for c in range(3)]
def h0(a):
    v,cnt=np.unique(a.ravel(),return_counts=True); p=cnt/a.size; return float(-(p*np.log2(p)).sum())
def lowbit(a):
    o=int(np.bitwise_or.reduce(np.abs(a).ravel())); return (o&-o).bit_length()-1 if o else '-'
for name,i0y,i0x,F in (('tissue 1024',64*157+pad,64*237+pad,1024),('tissue edge 1024',64*137+pad,64*137+pad,1024)):
    pix=[d.out_region(0,i0y,i0y+F,i0x,i0x+F) for d in decs]; D=5
    py=[full_pyramid(p,D) for p in pix]
    R,G,B=inv_rgb(pix[0],pix[1],pix[2],fl); Yr,Cbr,Crr=rct_fwd(R,G,B)         # exact pixel-domain chroma (inverse=floor)
    ex={'Cb':full_pyramid(Cbr,D),'Cr':full_pyramid(Crr,D)}
    (LLy,lvy),cbw,crw=convert_pyramid(py,fl); wd={'Cb':cbw,'Cr':crw}
    print('\n%s: per level (1 = finest): non-zero %% / H0 bits per coefficient / empty LSB planes, detail bands only' % name)
    print(' level | Co stored        Cg stored        | Cb wavelet-dom   Cb pixel-dom(exact) | Cr wavelet-dom   Cr pixel-dom(exact)')
    for j in range(D):
        def st(lv): a=np.concatenate([x.ravel() for x in lv[j]]); return '%6.2f%% %5.2fb L%s' % (100*np.mean(a!=0),h0(a),lowbit(a))
        print(' %5d | %-16s %-16s | %-16s %-19s | %-16s %-16s' % (j+1,st(py[1][1]),st(py[2][1]),st(wd['Cb'][1]),st(ex['Cb'][1]),st(wd['Cr'][1]),st(ex['Cr'][1])))
