# Savings from omitting tiles, from the per-tile CSVs written by omission_size_estimate.py.   python3 -B -I omission_summary.py [dir]
# Strategy A : omit only the exact filler tiles (YCoCg = 255, 0, 0 everywhere; synthetic fill, not scanned glass).
# Strategy B(T): additionally omit scanned-glass tiles whose pixels are nearly constant: robust RGB range (99.5 % minus 0.5 % percentile, worst channel) <= T
#                and no pixel outside the decoded array.  B(T) includes A.  "fill spread" = how far the median colour of the omitted GLASS tiles lies from one
#                common fill colour (their median): a single recommended absent-pixel colour would be off by this much.
# Sizes are encoded codestream bytes only (wavelet-domain colour route); DICOM per-frame metadata is not included (assumed about 250 bytes/frame: assumption).
import sys, os, csv, glob, numpy as np
d = sys.argv[1] if len(sys.argv) > 1 else '../transcoding_review'
TS = [2, 3, 4, 5, 6, 8, 16]; SS = [8, 16, 32, 64]; META = 250       # SS: strict (max - min) RGB range limits, no outlier pixel allowed beyond it
def load(F, k):
    r = list(csv.DictReader(open(os.path.join(d, 'tile_stats_F%d_L%d.csv' % (F, k)))))
    return dict(strict=np.array([float(x['rgb_range']) for x in r]), exact=np.array([int(x['exact_filler']) for x in r], bool), rob=np.array([float(x['rgb_range_robust']) for x in r]), out=np.array([float(x['outside_frac']) for x in r]),
                med=np.array([[float(x['med_r']), float(x['med_g']), float(x['med_b'])] for x in r]), b=np.array([int(x['bytes']) for x in r]), n=len(r))
for F in (256, 512, 1024):
    tot_all = 0; sav = {'A': 0}; sav.update({T: 0 for T in TS}); cnt = {'A': 0}; cnt.update({T: 0 for T in TS}); nframes = 0; ssav = {S: 0 for S in SS}; scnt = {S: 0 for S in SS}
    print('\n=== frame size %d x %d ===' % (F, F))
    print(' layer | tiles | all tiles MB | A: omit filler only: tiles / MB saved (%%) | B(T): + near-constant glass, tiles / MB saved (%%) [fill spread: 95th pct / max, RGB levels]')
    for k in range(6):
        t = load(F, k); tot = t['b'].sum(); tot_all += tot; nframes += t['n']; A = t['exact']; line = ' %5d | %5d | %7.3f | A %4d / %6.3f (%4.1f%%)' % (k, t['n'], tot / 1e6, A.sum(), t['b'][A].sum() / 1e6, 100 * t['b'][A].sum() / tot)
        sav['A'] += t['b'][A].sum(); cnt['A'] += A.sum(); parts = []
        for S in SS: ms = A | ((t['strict'] <= S) & (t['out'] == 0)); ssav[S] += t['b'][ms].sum(); scnt[S] += ms.sum()
        for T in TS:
            m = A | ((t['rob'] <= T) & (t['out'] == 0)); sav[T] += t['b'][m].sum(); cnt[T] += m.sum(); gl = m & ~A
            if gl.sum() > 1:
                dev = np.abs(t['med'][gl] - np.median(t['med'][gl], axis=0)).max(1); sp = '%.1f/%.1f' % (np.percentile(dev, 95), dev.max())
            else: sp = '-'
            parts.append('T=%-2d %4d / %6.3f (%4.1f%%) [%s]' % (T, m.sum(), t['b'][m].sum() / 1e6, 100 * t['b'][m].sum() / tot, sp))
        if k <= 1: print(line + ' | ' + ' | '.join(parts))
    print(' ALL 6 LAYERS: %d frames, %.2f MB encoded' % (nframes, tot_all / 1e6))
    print('   A (filler only): %d frames omitted, saves %.2f MB = %.2f%% (+ about %.2f MB metadata at %d B/frame)' % (cnt['A'], sav['A'] / 1e6, 100 * sav['A'] / tot_all, cnt['A'] * META / 1e6, META))
    for T in TS:
        t0_ = load(F, 0); gl = (t0_['rob'] <= T) & (t0_['out'] == 0) & ~t0_['exact']
        st = t0_['strict'][gl]; er = ('layer-0 glass tiles: strict RGB range (bound on the worst pixel error vs the tile median) median %.0f, 95th pct %.0f, max %.0f' % (np.median(st), np.percentile(st, 95), st.max())) if gl.sum() else '-'
        print('   B(T=%2d): %d frames omitted, saves %.2f MB = %.2f%% (+ about %.2f MB metadata) | %s' % (T, cnt[T], sav[T] / 1e6, 100 * sav[T] / tot_all, cnt[T] * META / 1e6, er))
    for S in SS: print('   strict range <= %2d (no speck beyond S levels): %d frames omitted, saves %.2f MB = %.2f%%' % (S, scnt[S], ssav[S] / 1e6, 100 * ssav[S] / tot_all))
