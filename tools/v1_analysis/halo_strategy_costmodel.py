# Core-second model of the halo strategies, from the logs of halo_strategy_footprint.py (decodes, bytes, transform samples) and measured constants.
#   python3 -B -I halo_strategy_costmodel.py <footprint.log> <iSyntax file bytes> [coded_pyramid_ratio=2.3]
# Constants (1.isyntax, build-ssse3 = Linux Meson -mssse3, one core; see Halo_Strategy_Comparison_Claude.md):
#   decode   = 30.9 us per non-empty codeblock + 7.7 ms per MB compressed (least-squares fit of the six per-level C timings, within 5 %)
#   lifting/forward transforms = 3 ns per sample (cold libisyntax_tile_read minus decode, includes RGB conversion, so slightly pessimistic)
#   Tier-1+2 = 6.5 MB/s of coded output per core (tier1_profile.py, noisy); coded pyramid = ratio x file bytes (2.3 measured on 1.isyntax only)
#   strip sharing of IDWT cells: the 'union' saving measured on 1.isyntax (65 % of the lifted samples) is assumed for every file (UNION=0 on the big ones)
import sys, re
log = open(sys.argv[1]).read(); fbytes = float(sys.argv[2]); ratio = float(sys.argv[3]) if len(sys.argv) > 3 else 2.3
A, B, NS, T1 = 30.9e-6, 7.7e-3, 3e-9, 6.5e6
def grab(pat):
    m = re.search(pat, log); return [float(v) for v in m.groups()] if m else None
def dec(n, mb): return A * n + B * mb
strip_layer = grab(r'STRIP per layer.*?: (\d+) decodes, ([\d.]+) MB'); strip_shared = grab(r'STRIP shared.*?: (\d+) decodes, ([\d.]+) MB')
thin0 = grab(r'THIN, per-frame task, no cache\s+: (\d+) decodes.*?, ([\d.]+) MB'); rec0 = grab(r'RECRUIT whole tiles, no cache\s+: (\d+) decodes.*?, ([\d.]+) MB')
lru = {c: grab(r'THIN with LRU cache of\s+%d blocks\s+: (\d+) decodes.*?, ([\d.]+) MB' % c) for c in (256, 1024, 4096)}
w = grab(r'THIN: lifted ([\d.]+) M output samples \+ strip forward ([\d.]+) M samples ; RECRUIT: (\d+) tile inverse transforms = ([\d.]+) M samples')
frames = int(re.search(r'over (\d+) interior frames', log).group(1)); lift, fwd, rtiles, rlift = w[0] * 1e6, w[1] * 1e6, w[2], w[3] * 1e6
frame_fwd = frames * 256 * 256 * 1.33                                  # RECRUIT still runs a full forward transform per frame
t1 = fbytes * ratio / T1
def tr(samples_per_colour): return 3 * samples_per_colour * NS
rows = [('STRIP (per-layer windows, shared IDWT cells)', dec(*strip_layer), tr(lift * 0.35 + fwd)),
        ('STRIP, all layers in lockstep', dec(*strip_shared), tr(lift * 0.35 + fwd)),
        ('THIN per-frame tasks + LRU cache (smallest sufficient)', dec(*(lru[4096] or lru[1024])), tr(lift + fwd)),
        ('THIN per-frame tasks, no cache', dec(*thin0), tr(lift + fwd)),
        ('RECRUIT whole tiles, no cache', dec(*rec0), tr(rlift + frame_fwd))]
base = rows[0][1] + rows[0][2] + t1
print('Tier-1+2 (hypothesis: coded pyramid = %.1f x file) : %7.0f core-s' % (ratio, t1))
print('%-58s %9s %11s %9s %7s' % ('strategy', 'decode s', 'transform s', 'total s', 'vs strip'))
for n, d, t in rows: print('%-58s %9.1f %11.1f %9.0f %6.2fx' % (n, d, t, d + t + t1, (d + t + t1) / base))
print('strip total per MB of file: %.2f core-s/MB' % (base / (fbytes / 1e6)))
