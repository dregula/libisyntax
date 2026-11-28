// iSyntax (v1) -> DICOM-WSI JPEG 2000 frame coefficients: STRIP WINDOW implementation (added in this fork; new file, decoder untouched).
//
// For every layer k of the DICOM pyramid and every F x F frame (raster order) it produces the reversible 5/3 coefficients that a JPEG 2000 encoder would hold
// for that frame, WITHOUT leaving the wavelet domain: interior detail coefficients are copied from the stored iSyntax H blocks, a ring of depth h along the
// frame edges is recomputed from thin pixel strips obtained by a windowed top-down IDWT (see iSyntax_to_DICOM_WSI_Transcoding_Claude.md section 4 and
// tools/v1_analysis/windowed_idwt_halo.py, which this file ports). Frames are scheduled in raster rows over ONE sliding window of decoded codeblocks per layer
// (Halo_Strategy_Comparison_Claude.md): the window holds the blocks the current and the next frame row need and drops the rest, so each block is decoded
// once per layer; frame tasks of a row run in parallel over the read-only window.
//   Wavelet-domain colour conversion (Pipeline_Plan_v2 section 2.2): Y' = Y, Cb = -q(Co) - Cg, Cr = Cb + Co (q = floor or truncating halving), then -128 on LL of Y'.
//   Omission (layers 0 and 1 only): frames whose RGB strict range (truncating inverse, clamped) is <= S are omitted; a coefficient prefilter avoids the pixel
//   decode of obvious tissue.  Fast path: a frame whose whole footprint consists of empty blocks and a constant top LL is a constant frame (exact).
//   Tier-1: patched OpenJPEG 2.5.4 (opj_set_coefficient_input), see patches/ and tools/transcode_tests/build_patched_openjpeg.sh.
// Build: tools/transcode_tests/build_strip_transcoder.sh      Run: isyntax_strip_transcode FILE [options]  (see usage())
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <chrono>
#include <algorithm>
#include "libisyntax.h"
#include "openjpeg.h"
#include "dicom_layer_writer.h"
#ifndef DEFAULT_ICC
#define DEFAULT_ICC "src/transcode/srgb_lcms.icc"
#endif

typedef int32_t i32; typedef int16_t i16; typedef int64_t i64; typedef uint64_t u64;
[[noreturn]] static void die(const char* m, long a = 0, long b = 0, long c = 0) { fprintf(stderr, "FATAL: %s (%ld %ld %ld)\n", m, a, b, c); fflush(stderr); abort(); }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ------------------------------------------------------------------------------------------------ matrices and lifting
struct Mat {
    int r = 0, c = 0; std::vector<i32> v;
    Mat() {}
    Mat(int r_, int c_) : r(r_), c(c_), v((size_t)r_ * c_, 0) {}
    i32& at(int y, int x) { return v[(size_t)y * c + x]; }
    const i32& at(int y, int x) const { return v[(size_t)y * c + x]; }
};
static Mat crop(const Mat& m, int r0, int r1, int c0, int c1) {
    Mat o(r1 - r0, c1 - c0); for (int y = r0; y < r1; y++) memcpy(&o.v[(size_t)(y - r0) * o.c], &m.v[(size_t)y * m.c + c0], (size_t)(c1 - c0) * 4); return o;
}
static void paste(Mat& d, int r0, int c0, const Mat& s) { for (int y = 0; y < s.r; y++) memcpy(&d.v[(size_t)(r0 + y) * d.c + c0], &s.v[(size_t)y * s.c], (size_t)s.c * 4); }
// iSyntax inverse lifting (reconstruct_np.py lift_axis): L, H same shape; 2N outputs, even <- H path, odd <- L path, mirrored (duplicated) at the window ends.
static Mat lift_cols(const Mat& L, const Mat& H) {
    int N = L.c; Mat out(L.r, 2 * N); std::vector<i32> xl(N);
    for (int y = 0; y < L.r; y++) {
        const i32* l = &L.v[(size_t)y * N]; const i32* h = &H.v[(size_t)y * N]; i32* o = &out.v[(size_t)y * 2 * N];
        for (int n = 0; n < N; n++) { i32 hn = n + 1 < N ? h[n + 1] : h[N - 1]; xl[n] = l[n] - ((h[n] + hn + 2) >> 2); }
        for (int n = 0; n < N; n++) { i32 p = n ? xl[n - 1] : xl[0]; o[2 * n] = h[n] + ((p + xl[n]) >> 1); o[2 * n + 1] = xl[n]; }
    }
    return out;
}
static Mat lift_rows(const Mat& L, const Mat& H) {
    int N = L.r, C = L.c; Mat out(2 * N, C), xl(N, C);
    for (int n = 0; n < N; n++) { const i32* l = &L.v[(size_t)n * C]; const i32* h = &H.v[(size_t)n * C]; const i32* hn = &H.v[(size_t)(n + 1 < N ? n + 1 : N - 1) * C]; i32* x = &xl.v[(size_t)n * C];
        for (int i = 0; i < C; i++) x[i] = l[i] - ((h[i] + hn[i] + 2) >> 2); }
    for (int n = 0; n < N; n++) { const i32* h = &H.v[(size_t)n * C]; const i32* p = &xl.v[(size_t)(n ? n - 1 : 0) * C]; const i32* x = &xl.v[(size_t)n * C]; i32* o0 = &out.v[(size_t)(2 * n) * C]; i32* o1 = &out.v[(size_t)(2 * n + 1) * C];
        for (int i = 0; i < C; i++) { o0[i] = h[i] + ((p[i] + x[i]) >> 1); o1[i] = x[i]; } }
    return out;
}
// standard JPEG 2000 reversible 5/3 (T.800), mirrored at the ends; n even
static void fwd_rows(const Mat& x, Mat& s, Mat& d) {
    int n2 = x.r / 2, C = x.c; s = Mat(n2, C); d = Mat(n2, C);
    for (int n = 0; n < n2; n++) { const i32* ev = &x.v[(size_t)(2 * n) * C]; const i32* od = &x.v[(size_t)(2 * n + 1) * C]; const i32* evn = &x.v[(size_t)(2 * (n + 1 < n2 ? n + 1 : n2 - 1)) * C]; i32* dd = &d.v[(size_t)n * C];
        for (int i = 0; i < C; i++) dd[i] = od[i] - ((ev[i] + evn[i]) >> 1); }
    for (int n = 0; n < n2; n++) { const i32* ev = &x.v[(size_t)(2 * n) * C]; const i32* d0 = &d.v[(size_t)n * C]; const i32* dp = &d.v[(size_t)(n ? n - 1 : 0) * C]; i32* ss = &s.v[(size_t)n * C];
        for (int i = 0; i < C; i++) ss[i] = ev[i] + ((dp[i] + d0[i] + 2) >> 2); }
}
static void fwd_cols(const Mat& x, Mat& s, Mat& d) {
    int n2 = x.c / 2, R = x.r; s = Mat(R, n2); d = Mat(R, n2);
    for (int y = 0; y < R; y++) { const i32* row = &x.v[(size_t)y * x.c]; i32* sr = &s.v[(size_t)y * n2]; i32* dr = &d.v[(size_t)y * n2];
        for (int n = 0; n < n2; n++) { i32 evn = row[2 * (n + 1 < n2 ? n + 1 : n2 - 1)]; dr[n] = row[2 * n + 1] - ((row[2 * n] + evn) >> 1); }
        for (int n = 0; n < n2; n++) { i32 dp = dr[n ? n - 1 : 0]; sr[n] = row[2 * n] + ((dp + dr[n] + 2) >> 2); } }
}
struct Q4 { Mat LL, HL, LH, HH; };
static Q4 fwd2(const Mat& x) { Mat lo, hi; fwd_rows(x, lo, hi); Q4 q; fwd_cols(lo, q.LL, q.HL); fwd_cols(hi, q.LH, q.HH); return q; }   // vertical, then horizontal
static Mat inv_cols(const Mat& s, const Mat& d) {
    int N = s.c; Mat out(s.r, 2 * N);
    for (int y = 0; y < s.r; y++) { const i32* sr = &s.v[(size_t)y * N]; const i32* dr = &d.v[(size_t)y * N]; i32* o = &out.v[(size_t)y * 2 * N]; std::vector<i32> ev(N);
        for (int n = 0; n < N; n++) { i32 dp = dr[n ? n - 1 : 0]; ev[n] = sr[n] - ((dp + dr[n] + 2) >> 2); }
        for (int n = 0; n < N; n++) { i32 evn = n + 1 < N ? ev[n + 1] : ev[N - 1]; o[2 * n] = ev[n]; o[2 * n + 1] = dr[n] + ((ev[n] + evn) >> 1); } }
    return out;
}
static Mat inv_rows(const Mat& s, const Mat& d) {
    int N = s.r, C = s.c; Mat out(2 * N, C), ev(N, C);
    for (int n = 0; n < N; n++) { const i32* sr = &s.v[(size_t)n * C]; const i32* dr = &d.v[(size_t)n * C]; const i32* dp = &d.v[(size_t)(n ? n - 1 : 0) * C]; i32* e = &ev.v[(size_t)n * C];
        for (int i = 0; i < C; i++) e[i] = sr[i] - ((dp[i] + dr[i] + 2) >> 2); }
    for (int n = 0; n < N; n++) { const i32* dr = &d.v[(size_t)n * C]; const i32* e = &ev.v[(size_t)n * C]; const i32* en = &ev.v[(size_t)(n + 1 < N ? n + 1 : N - 1) * C]; i32* o0 = &out.v[(size_t)(2 * n) * C]; i32* o1 = &out.v[(size_t)(2 * n + 1) * C];
        for (int i = 0; i < C; i++) { o0[i] = e[i]; o1[i] = dr[i] + ((e[i] + en[i]) >> 1); } }
    return out;
}
static Mat inv2(const Mat& LL, const Mat& HL, const Mat& LH, const Mat& HH) { return inv_rows(inv_cols(LL, HL), inv_cols(LH, HH)); }   // horizontal, then vertical

struct Pyr { Mat LL; std::vector<std::array<Mat, 3>> lv; };      // lv[j] = (HL, LH, HH) of J2K level j (0 = finest); LL = LL_D
static Pyr full_pyramid(const Mat& fr, int D) {
    Pyr p; Mat cur = fr;
    for (int j = 0; j < D; j++) { Q4 q = fwd2(cur); p.lv.push_back({q.HL, q.LH, q.HH}); cur = q.LL; }
    p.LL = cur; return p;
}
static Mat inv_pyramid(const Pyr& p) { Mat cur = p.LL; for (int j = (int)p.lv.size() - 1; j >= 0; j--) cur = inv2(cur, p.lv[j][0], p.lv[j][1], p.lv[j][2]); return cur; }
static bool same(const Mat& a, const Mat& b) { return a.r == b.r && a.c == b.c && a.v == b.v; }
static bool same(const Pyr& a, const Pyr& b) { if (!same(a.LL, b.LL) || a.lv.size() != b.lv.size()) return false; for (size_t j = 0; j < a.lv.size(); j++) for (int t = 0; t < 3; t++) if (!same(a.lv[j][t], b.lv[j][t])) return false; return true; }

// ------------------------------------------------------------------------------------------------ geometry and the window of decoded blocks
static int need(i64 n, int s) { return (int)((n + ((i64)1 << (s + 1)) - 1) >> (s + 1)) + 2; }
struct Geo {
    int L = 0, top = 0; i64 Wp = 0, Hp = 0, pad = 0; std::vector<int> limR, limC, nbr, nbc;
    int R(int s) const { return limR[s]; } int C(int s) const { return limC[s]; }
};
struct Blk { i16 d[9 * 16384]; bool zero = false; };            // kind 1 (H): [colour][band HL,LH,HH][128][128];  kind 0 (LL): [colour][128][128]
static inline u64 bkey(int s, int kind, int bx, int by) { return ((u64)s << 56) | ((u64)kind << 54) | ((u64)bx << 27) | (u64)by; }
struct Window {
    isyntax_t* is = nullptr; const Geo* g = nullptr; std::unordered_map<u64, std::unique_ptr<Blk>> m; std::atomic<long> decodes{0}; size_t peak = 0;
    const Blk* find(int s, int kind, int bx, int by) const {
        if (bx < 0 || by < 0 || bx >= g->nbc[s] || by >= g->nbr[s]) return nullptr;                 // outside the tile grid: zeros
        auto it = m.find(bkey(s, kind, bx, by)); if (it == m.end()) die("window miss (planner bug)", s, bx, by); return it->second.get();
    }
};
// H coefficients (band 0..2) of colour c at scale s, rows [r0,r1), cols [c0,c1) -> Mat; kind 0 = LL (top only)
static Mat fetch(const Window& w, int c, int kind, int band, int s, int r0, int r1, int c0, int c1) {
    if (r0 < 0 || c0 < 0 || r1 > w.g->R(s) || c1 > w.g->C(s) || r1 <= r0 || c1 <= c0) die("fetch outside the coefficient array", s, r0, r1);
    Mat o(r1 - r0, c1 - c0);
    for (int by = r0 / 128; by <= (r1 - 1) / 128; by++) for (int bx = c0 / 128; bx <= (c1 - 1) / 128; bx++) {
        const Blk* b = w.find(s, kind, bx, by); if (!b) continue;
        const i16* d = b->d + (size_t)(kind ? c * 3 + band : c) * 16384;
        int ya = std::max(r0, by * 128), yb = std::min(r1, by * 128 + 128), xa = std::max(c0, bx * 128), xb = std::min(c1, bx * 128 + 128);
        for (int y = ya; y < yb; y++) for (int x = xa; x < xb; x++) o.at(y - r0, x - c0) = d[(y - by * 128) * 128 + (x - bx * 128)];
    }
    return o;
}

// ------------------------------------------------------------------------------------------------ planner: which blocks does a frame touch? (indices only, no data)
struct Plan {
    std::set<u64> blocks; std::vector<std::array<int, 6>> llrects;      // top-scale LL rectangles read (r0,r1,c0,c1) ; clipped: any window touched an array limit
    bool clipped = false;                                              // a window reached index 0: the zero shift L'[0] = 0 enters, so no constant shortcut
};
struct Planner {
    const Geo& g; Plan& p; Planner(const Geo& g_, Plan& p_) : g(g_), p(p_) {}
    void touch(int s, int kind, int r0, int r1, int c0, int c1) { for (int by = r0 / 128; by <= (r1 - 1) / 128; by++) for (int bx = c0 / 128; bx <= (c1 - 1) / 128; bx++) p.blocks.insert(bkey(s, kind, bx, by)); }
    void cur_region(int s, int r0, int r1, int c0, int c1) {
        if (!(0 <= r0 && r0 < r1 && r1 <= g.R(s) && 0 <= c0 && c0 < c1 && c1 <= g.C(s))) die("planner: cur_region outside", s, r0, c0);
        if (s == g.top) { touch(s, 0, r0, r1, c0, c1); p.llrects.push_back({s, r0, r1, c0, c1, 0}); return; }
        int pr = r0 == 0, pc = c0 == 0; if (pr || pc) p.clipped = true;
        if (r1 - pr <= 0 || c1 - pc <= 0) return;
        out_region(s + 1, r0 - 1 + pr, r1 - 1, c0 - 1 + pc, c1 - 1);
    }
    void out_region(int s, int r0, int r1, int c0, int c1) {
        int R = g.R(s), C = g.C(s); if (!(0 <= r0 && r0 < r1 && r1 <= 2 * R && 0 <= c0 && c0 < c1 && c1 <= 2 * C)) die("planner: out_region outside", s, r0, r1);
        int nr0 = std::max(0, (r0 >> 1) - 1), nr1 = std::min(R, ((r1 - 1) >> 1) + 3), nc0 = std::max(0, (c0 >> 1) - 1), nc1 = std::min(C, ((c1 - 1) >> 1) + 3);
        cur_region(s, nr0, nr1, nc0, nc1); touch(s, 1, nr0, nr1, nc0, nc1);
    }
};
template <class G> static void rect_fill_calls(int n, G getter, int W) {      // the calls rect_fill makes (rect_fill itself is below, with data)
    if (n <= 2 * W) { getter(0, n, 0, n); return; }
    getter(0, n, 0, W); getter(0, n, n - W, n); getter(0, W, 0, n); getter(n - W, n, 0, n);
}
static void plan_interior_frame(const Geo& g, Plan& p, int k, int i0y, int i0x, int F, int D, int h) {
    Planner pl(g, p); int W = 2 * h + 4, qy = i0y + 1, qx = i0x + 1;
    rect_fill_calls(F, [&](int r0, int r1, int c0, int c1) { pl.out_region(k, i0y + r0, i0y + r1, i0x + c0, i0x + c1); }, W); int n = F;
    for (int j = 0; j < D; j++) {
        int s = k + j, nn = n / 2; if (s > g.top) { n = nn; continue; }
        int hx = ((qx + 2) >> (j + 1)) - 1, hy = ((qy + 2) >> (j + 1)) - 1;
        pl.touch(s, 1, hy - 1, hy + nn - 1, hx, hx + nn); pl.touch(s, 1, hy, hy + nn, hx - 1, hx + nn - 1); pl.touch(s, 1, hy, hy + nn, hx, hx + nn);
        int by = ((qy + 2) >> (j + 1)) - 2, bx = ((qx + 2) >> (j + 1)) - 2;
        if (j == D - 1 || s == g.top || nn <= 2 * W) pl.cur_region(s, by, by + nn, bx, bx + nn);
        else rect_fill_calls(nn, [&](int r0, int r1, int c0, int c1) { pl.cur_region(s, by + r0, by + r1, bx + c0, bx + c1); }, W);
        n = nn;
    }
}
static void plan_pixel_frame(const Geo& g, Plan& p, int k, int i0y, int i0x, int F) {
    int R = 2 * g.R(k), C = 2 * g.C(k), y0 = std::max(i0y, 0), x0 = std::max(i0x, 0), y1 = std::min(i0y + F, R), x1 = std::min(i0x + F, C);
    if (y1 > y0 && x1 > x0) { Planner pl(g, p); pl.out_region(k, y0, y1, x0, x1); }
}

// ------------------------------------------------------------------------------------------------ engine: windowed IDWT and the frame builder (data)
struct Engine {
    const Window& w; const Geo& g; int c;
    Engine(const Window& w_, int c_) : w(w_), g(*w_.g), c(c_) {}
    Mat cur_region(int s, int r0, int r1, int c0, int c1) {
        if (!(0 <= r0 && r0 < r1 && r1 <= g.R(s) && 0 <= c0 && c0 < c1 && c1 <= g.C(s))) die("cur_region outside", s, r0, c0);
        if (s == g.top) return fetch(w, c, 0, 0, s, r0, r1, c0, c1);
        int pr = r0 == 0, pc = c0 == 0;
        if (r1 - pr <= 0 || c1 - pc <= 0) return Mat(r1 - r0, c1 - c0);
        Mat sub = out_region(s + 1, r0 - 1 + pr, r1 - 1, c0 - 1 + pc, c1 - 1), o(r1 - r0, c1 - c0); paste(o, pr, pc, sub); return o;
    }
    Mat out_region(int s, int r0, int r1, int c0, int c1) {
        int R = g.R(s), C = g.C(s); if (!(0 <= r0 && r0 < r1 && r1 <= 2 * R && 0 <= c0 && c0 < c1 && c1 <= 2 * C)) die("out_region outside", s, r0, r1);
        int nr0 = std::max(0, (r0 >> 1) - 1), nr1 = std::min(R, ((r1 - 1) >> 1) + 3), nc0 = std::max(0, (c0 >> 1) - 1), nc1 = std::min(C, ((c1 - 1) >> 1) + 3);
        Mat L = cur_region(s, nr0, nr1, nc0, nc1);
        Mat HL = fetch(w, c, 1, 0, s, nr0, nr1, nc0, nc1), LH = fetch(w, c, 1, 1, s, nr0, nr1, nc0, nc1), HH = fetch(w, c, 1, 2, s, nr0, nr1, nc0, nc1);
        Mat out = lift_rows(lift_cols(L, HL), lift_cols(LH, HH));
        return crop(out, r0 - 2 * nr0, r1 - 2 * nr0, c0 - 2 * nc0, c1 - 2 * nc0);
    }
};
template <class G> static void rect_fill(Mat& a, G getter, int W) {
    int n = a.r;
    if (n <= 2 * W) { paste(a, 0, 0, getter(0, n, 0, n)); return; }
    paste(a, 0, 0, getter(0, n, 0, W)); paste(a, 0, n - W, getter(0, n, n - W, n)); paste(a, 0, 0, getter(0, W, 0, n)); paste(a, n - W, 0, getter(n - W, n, 0, n));
}
// port of windowed_idwt_halo.py build_windowed(): J2K coefficients (YCoCg domain) of the interior frame of layer k whose first pixel is out_k[i0y, i0x]
static Pyr build_frame(Engine& e, int k, int i0y, int i0x, int F, int D, int h) {
    int W = 2 * h + 4, qy = i0y + 1, qx = i0x + 1, top = e.g.top; if (((qy + 2) & ((1 << D) - 1)) || ((qx + 2) & ((1 << D) - 1))) die("frame origin not aligned", i0y, i0x);
    Mat curf(F, F);
    rect_fill(curf, [&](int r0, int r1, int c0, int c1) { return e.out_region(k, i0y + r0, i0y + r1, i0x + c0, i0x + c1); }, W);
    Pyr pr;
    for (int j = 0; j < D; j++) {
        int s = k + j, n = curf.r, nn = n / 2;
        if (s > top) { Q4 q = fwd2(curf); pr.lv.push_back({q.HL, q.LH, q.HH}); curf = q.LL; continue; }
        int hx = ((qx + 2) >> (j + 1)) - 1, hy = ((qy + 2) >> (j + 1)) - 1, by = ((qy + 2) >> (j + 1)) - 2, bx = ((qx + 2) >> (j + 1)) - 2;
        Mat HLc = fetch(e.w, e.c, 1, 0, s, hy - 1, hy + nn - 1, hx, hx + nn), LHc = fetch(e.w, e.c, 1, 1, s, hy, hy + nn, hx - 1, hx + nn - 1), HHc = fetch(e.w, e.c, 1, 2, s, hy, hy + nn, hx, hx + nn);
        Mat LLc(nn, nn);
        if (j == D - 1 || s == top || nn <= 2 * W) LLc = e.cur_region(s, by, by + nn, bx, bx + nn);
        else rect_fill(LLc, [&](int r0, int r1, int c0, int c1) { return e.cur_region(s, by + r0, by + r1, bx + c0, bx + c1); }, W);
        if (nn <= 2 * h) { Q4 q = fwd2(curf); LLc = q.LL; HLc = q.HL; LHc = q.LH; HHc = q.HH; }
        else {
            auto put = [&](const Q4& q, int rs, int cs, int tr, int tc, int nrw, int ncl) {
                Mat* dst[4] = {&LLc, &HLc, &LHc, &HHc}; const Mat* src[4] = {&q.LL, &q.HL, &q.LH, &q.HH};
                for (int t = 0; t < 4; t++) paste(*dst[t], rs, cs, crop(*src[t], tr, tr + nrw, tc, tc + ncl));
            };
            Q4 q = fwd2(crop(curf, 0, n, 0, W));      put(q, 0, 0, 0, 0, nn, h);
            q = fwd2(crop(curf, 0, n, n - W, n));     put(q, 0, nn - h, 0, 2, nn, h);
            q = fwd2(crop(curf, 0, W, 0, n));         put(q, 0, 0, 0, 0, h, nn);
            q = fwd2(crop(curf, n - W, n, 0, n));     put(q, nn - h, 0, 2, 0, h, nn);
        }
        pr.lv.push_back({HLc, LHc, HHc}); curf = LLc;
    }
    pr.LL = curf; return pr;
}
static Mat pixel_frame(Engine& e, int k, int i0y, int i0x, int F) {                    // frame_from_decode(): decoded part, zero filler outside the array
    int R = 2 * e.g.R(k), C = 2 * e.g.C(k), y0 = std::max(i0y, 0), x0 = std::max(i0x, 0), y1 = std::min(i0y + F, R), x1 = std::min(i0x + F, C);
    Mat fr(F, F); if (e.c == 0) std::fill(fr.v.begin(), fr.v.end(), 255);                    // outside the array: white (Y = 255, Co = Cg = 0), as the iSyntax margin filler
    if (y1 > y0 && x1 > x0) paste(fr, y0 - i0y, x0 - i0x, e.out_region(k, y0, y1, x0, x1)); return fr;
}
// Layer k+1 from layer k: out_{k+1}[m] = 5/3 low-pass (T.800) of out_k centred on sample 2m + 3 (cross-layer relation i = 2 i' + 3, Transcoding doc 8a), vertical then horizontal, mirrored at the ends.
static inline int refl(int i, int N) { if (N == 1) return 0; int per = 2 * (N - 1); i = i < 0 ? -i : i; i %= per; return i > N - 1 ? per - i : i; }
static Mat lp_down(const Mat& x) {
    auto lp = [](auto at, int N, int m) { int i = 2 * m + 3; auto X = [&](int j) { return at(refl(j, N)); };
        i32 dm = X(i - 1) - ((X(i - 2) + X(i)) >> 1), dp = X(i + 1) - ((X(i) + X(i + 2)) >> 1); return X(i) + ((dm + dp + 2) >> 2); };
    int Rn = std::max(1, (x.r - 2) / 2), Cn = std::max(1, (x.c - 2) / 2); Mat t(Rn, x.c), o(Rn, Cn);
    for (int m = 0; m < Rn; m++) for (int c = 0; c < x.c; c++) t.at(m, c) = lp([&](int j) { return x.at(j, c); }, x.r, m);
    for (int r = 0; r < Rn; r++) for (int m = 0; m < Cn; m++) o.at(r, m) = lp([&](int j) { return t.at(r, j); }, x.c, m);
    return o;
}
static Pyr const_pyramid(int F, int D, i32 v) { Pyr p; for (int j = 0; j < D; j++) { int n = F >> (j + 1); p.lv.push_back({Mat(n, n), Mat(n, n), Mat(n, n)}); } p.LL = Mat(F >> D, F >> D); std::fill(p.LL.v.begin(), p.LL.v.end(), v); return p; }

// ------------------------------------------------------------------------------------------------ colour, Mallat layout
static inline i32 qhalf(i32 x, bool trunc_) { return trunc_ ? (x >= 0 ? x >> 1 : -((-x) >> 1)) : (x >> 1); }
static void colour_convert(Pyr& Y, Pyr& Co, Pyr& Cg, bool trunc_) {
    auto cv = [&](Mat& a, Mat& b, Mat& c3) { for (size_t i = 0; i < a.v.size(); i++) { i32 co = b.v[i], cg = c3.v[i]; i32 cb = -qhalf(co, trunc_) - cg; b.v[i] = cb; c3.v[i] = cb + co; } };
    cv(Y.LL, Co.LL, Cg.LL); for (size_t j = 0; j < Y.lv.size(); j++) for (int t = 0; t < 3; t++) cv(Y.lv[j][t], Co.lv[j][t], Cg.lv[j][t]);
    for (auto& v : Y.LL.v) v -= 128;                                       // DC level shift on LL of Y' only (see Pipeline_Plan_v2 2.1)
}
static void mallat(const Pyr& p, int F, int D, i32* A) {                  // OpenJPEG layout, row stride F
    memset(A, 0, (size_t)F * F * 4); int hh = F >> D;
    for (int y = 0; y < hh; y++) for (int x = 0; x < hh; x++) A[y * F + x] = p.LL.at(y, x);
    for (int j = D - 1; j >= 0; j--) { int n = F >> (j + 1); const Mat& HL = p.lv[j][0]; const Mat& LH = p.lv[j][1]; const Mat& HH = p.lv[j][2];
        for (int y = 0; y < n; y++) for (int x = 0; x < n; x++) { A[y * F + n + x] = HL.at(y, x); A[(n + y) * F + x] = LH.at(y, x); A[(n + y) * F + n + x] = HH.at(y, x); } }
}
static Pyr unmallat(const i32* A, int F, int D) {
    Pyr p; int hh = F >> D; p.LL = Mat(hh, hh); for (int y = 0; y < hh; y++) for (int x = 0; x < hh; x++) p.LL.at(y, x) = A[y * F + x];
    p.lv.resize(D);
    for (int j = 0; j < D; j++) { int n = F >> (j + 1); Mat HL(n, n), LH(n, n), HH(n, n);
        for (int y = 0; y < n; y++) for (int x = 0; x < n; x++) { HL.at(y, x) = A[y * F + n + x]; LH.at(y, x) = A[(n + y) * F + x]; HH.at(y, x) = A[(n + y) * F + n + x]; }
        p.lv[j] = {HL, LH, HH}; }
    return p;
}

// ------------------------------------------------------------------------------------------------ Tier-1 through the patched OpenJPEG (coefficient input)
struct MS { std::vector<uint8_t>* b; size_t pos; };
static OPJ_SIZE_T ms_write(void* p, OPJ_SIZE_T n, void* u) { MS* m = (MS*)u; if (m->pos + n > m->b->size()) m->b->resize(m->pos + n); memcpy(m->b->data() + m->pos, p, n); m->pos += n; return n; }
static OPJ_SIZE_T ms_read(void* p, OPJ_SIZE_T n, void* u) { MS* m = (MS*)u; if (m->pos >= m->b->size()) return (OPJ_SIZE_T)-1; size_t k = std::min((size_t)n, m->b->size() - m->pos); memcpy(p, m->b->data() + m->pos, k); m->pos += k; return k; }
static OPJ_OFF_T ms_skip(OPJ_OFF_T n, void* u) { MS* m = (MS*)u; i64 np = (i64)m->pos + n; if (np < 0) return -1; if ((size_t)np > m->b->size()) m->b->resize(np); m->pos = np; return n; }
static OPJ_BOOL ms_seek(OPJ_OFF_T p, void* u) { MS* m = (MS*)u; if (p < 0) return OPJ_FALSE; if ((size_t)p > m->b->size()) m->b->resize(p); m->pos = p; return OPJ_TRUE; }
static void opj_quiet(const char*, void*) {}
static bool encode_frame(const i32* coef3, int F, std::vector<uint8_t>& out) {          // coef3 = 3 planes of F*F int32 (Mallat layout, colour converted, shifted)
    opj_cparameters_t p; opj_set_default_encoder_parameters(&p);
    p.tcp_numlayers = 1; p.cp_disto_alloc = 1; p.tcp_rates[0] = 0; p.irreversible = 0; p.tcp_mct = 1; p.numresolution = 6; p.cblockw_init = 64; p.cblockh_init = 64; p.cod_format = 0;
    opj_image_cmptparm_t cp[3]; memset(cp, 0, sizeof cp); for (int c = 0; c < 3; c++) { cp[c].dx = cp[c].dy = 1; cp[c].w = cp[c].h = F; cp[c].prec = 8; cp[c].sgnd = 0; }
    opj_image_t* im = opj_image_create(3, cp, OPJ_CLRSPC_SRGB); im->x0 = im->y0 = 0; im->x1 = im->y1 = F;
    opj_codec_t* cd = opj_create_compress(OPJ_CODEC_J2K); opj_set_error_handler(cd, opj_quiet, nullptr); opj_set_warning_handler(cd, opj_quiet, nullptr);
    bool ok = opj_setup_encoder(cd, &p, im) && opj_set_coefficient_input(cd, 1);
    out.clear(); MS ms{&out, 0}; opj_stream_t* st = opj_stream_create(1 << 16, OPJ_FALSE);
    opj_stream_set_write_function(st, ms_write); opj_stream_set_skip_function(st, ms_skip); opj_stream_set_seek_function(st, ms_seek); opj_stream_set_user_data(st, &ms, nullptr);
    ok = ok && opj_start_compress(cd, im, st) && opj_write_tile(cd, 0, (OPJ_BYTE*)coef3, (OPJ_UINT32)((size_t)F * F * 3 * 4), st) && opj_end_compress(cd, st);
    opj_stream_destroy(st); opj_destroy_codec(cd); opj_image_destroy(im); return ok;
}
static bool decode_codestream(const std::vector<uint8_t>& cs, int F, std::vector<i32>& planes) {   // stock decoder -> 3 planes of F*F
    std::vector<uint8_t> copy = cs; MS ms{&copy, 0}; opj_stream_t* st = opj_stream_create(1 << 16, OPJ_TRUE);
    opj_stream_set_read_function(st, ms_read); opj_stream_set_skip_function(st, ms_skip); opj_stream_set_seek_function(st, ms_seek); opj_stream_set_user_data(st, &ms, nullptr); opj_stream_set_user_data_length(st, copy.size());
    opj_dparameters_t dp; opj_set_default_decoder_parameters(&dp); opj_codec_t* cd = opj_create_decompress(OPJ_CODEC_J2K); opj_set_error_handler(cd, opj_quiet, nullptr); opj_set_warning_handler(cd, opj_quiet, nullptr);
    opj_setup_decoder(cd, &dp); opj_image_t* im = nullptr; bool ok = opj_read_header(st, cd, &im) && opj_decode(cd, st, im) && im && im->numcomps == 3;
    if (ok) { planes.resize((size_t)F * F * 3); for (int c = 0; c < 3; c++) for (size_t i = 0; i < (size_t)F * F; i++) planes[c * (size_t)F * F + i] = im->comps[c].data[i]; }
    opj_stream_destroy(st); opj_destroy_codec(cd); if (im) opj_image_destroy(im); return ok;
}

// ------------------------------------------------------------------------------------------------ thread pool
struct Pool {                                                  // persistent workers; every run() is its own Job (own counters), so a late-waking worker can never touch the next job's counters
    struct Job { std::function<void(int, int)> fn; int n = 0; std::atomic<int> next{0}, done{0}; };
    int T; std::vector<std::thread> th; std::mutex mu; std::condition_variable cv, cvd; std::shared_ptr<Job> job; long gen = 0; bool stop = false;
    explicit Pool(int t) : T(t) {
        for (int i = 0; i < T; i++) th.emplace_back([this, i] {
            long my = 0;
            for (;;) {
                std::shared_ptr<Job> j;
                { std::unique_lock<std::mutex> lk(mu); cv.wait(lk, [&] { return stop || gen != my; }); if (stop) break; my = gen; j = job; }
                for (;;) { int k = j->next++; if (k >= j->n) break; j->fn(k, i); if (++j->done == j->n) { std::lock_guard<std::mutex> lk(mu); cvd.notify_all(); } }
            }
            libisyntax_thread_shutdown();
        });
    }
    void run(int count, std::function<void(int, int)> f) {
        if (count <= 0) return;
        auto j = std::make_shared<Job>(); j->fn = std::move(f); j->n = count;
        { std::lock_guard<std::mutex> lk(mu); job = j; gen++; }
        cv.notify_all();
        { std::unique_lock<std::mutex> lk(mu); cvd.wait(lk, [&] { return j->done.load() >= count; }); }
    }
    ~Pool() { { std::lock_guard<std::mutex> lk(mu); stop = true; } cv.notify_all(); for (auto& t : th) t.join(); }
};

// ------------------------------------------------------------------------------------------------ the driver
struct Opt {
    std::string file, csv, dumpdir; int F = 256, D = 5, h = 4, threads = 8, S = 16, prefilter = -1, verify_every = 0, layer_lo = 0, layer_hi = 99, omit_layers = 2, max_rows = 0;
    std::string dicom_dir, icc = DEFAULT_ICC; double orient[6] = {0, -1, 0, -1, 0, 0}; double clamp_um = 0;
    bool trunc_ = true, encode = true, fastpath = true, audit = false, roundtrip = false, oob = false; std::vector<std::array<int, 3>> dumps;
};
struct Stat { std::atomic<long> const_omit{0}, white{0}, masked{0}, frames{0}, interior{0}, edge{0}, fast{0}, omitted{0}, prefilter_pass{0}, encoded{0}, verified{0}, verify_bad{0}, rt_checked{0}, rt_bad{0}, audit_false_reject{0}, audit_omit{0}, fast_verified{0}, fast_bad{0}; std::atomic<long> bytes{0};
    std::atomic<long> t_plan_ns{0}, t_decode_ns{0}, t_build_ns{0}, t_pixel_ns{0}, t_omit_ns{0}, t_encode_ns{0}, t_verify_ns{0}, t_fast_ns{0}; std::atomic<int> max_audit_abs[3] = {{0}, {0}, {0}}; };
static void usage() { fprintf(stderr, "usage: isyntax_strip_transcode FILE [--layers a-b] [--threads N] [--colour floor|trunc] [--omit S] [--no-fastpath] [--no-encode] [--verify-every K] [--roundtrip] [--audit]\n"
    "         [--csv out.csv] [--dump k ty tx dir] [--prefilter T] [--max-rows N]\n         [--dicom OUTDIR] [--icc file.icc] [--orientation a,b,c,d,e,f]   (DICOM: one file per layer, omitted glass frames as constant-colour frames)\n"); exit(2); }
struct FrameRec { int k, ty, tx, i0y, i0x; char kind, path; int omitted; long bytes; long oob; std::shared_ptr<std::vector<uint8_t>> cs; };

int main(int argc, char** argv) {
    if (argc < 2) usage();
    Opt o; o.file = argv[1];
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i]; auto nx = [&]() -> const char* { if (i + 1 >= argc) usage(); return argv[++i]; };
        if (a == "--layers") { int lo, hi; if (sscanf(nx(), "%d-%d", &lo, &hi) != 2) usage(); o.layer_lo = lo; o.layer_hi = hi; }
        else if (a == "--threads") o.threads = atoi(nx()); else if (a == "--colour") o.trunc_ = std::string(nx()) == "trunc"; else if (a == "--omit") o.S = atoi(nx());
        else if (a == "--no-fastpath") o.fastpath = false; else if (a == "--no-encode") o.encode = false; else if (a == "--verify-every") o.verify_every = atoi(nx());
        else if (a == "--roundtrip") o.roundtrip = true; else if (a == "--audit") o.audit = true; else if (a == "--oob") o.oob = true; else if (a == "--csv") o.csv = nx(); else if (a == "--prefilter") o.prefilter = atoi(nx());
        else if (a == "--max-rows") o.max_rows = atoi(nx());
        else if (a == "--origin-clamp-um") o.clamp_um = atof(nx()); else if (a == "--dicom") o.dicom_dir = nx(); else if (a == "--icc") o.icc = nx();
        else if (a == "--orientation") { if (sscanf(nx(), "%lf,%lf,%lf,%lf,%lf,%lf", o.orient, o.orient + 1, o.orient + 2, o.orient + 3, o.orient + 4, o.orient + 5) != 6) usage(); }
        else if (a == "--dump") { int k = atoi(nx()), ty = atoi(nx()), tx = atoi(nx()); o.dumpdir = nx(); o.dumps.push_back({k, ty, tx}); }
        else usage();
    }
    if (o.prefilter < 0) o.prefilter = 3 * o.S;                       // audit: max |H| of an omittable frame was 17 / 32 / 20 on files 1 / 5 / 3 (S = 16)
    double t0 = now(); libisyntax_init(); isyntax_t* is = nullptr;
    if (libisyntax_open(o.file.c_str(), (libisyntax_open_flags_t)0, &is) != LIBISYNTAX_OK) { fprintf(stderr, "cannot open %s\n", o.file.c_str()); return 1; }
    const isyntax_image_t* wsi = libisyntax_get_wsi_image(is); Geo g; g.L = libisyntax_image_get_level_count(wsi); g.top = g.L - 1; g.pad = (3 << g.L) - 3;
    g.Wp = libisyntax_level_get_width(libisyntax_image_get_level(wsi, 0)) + 2 * g.pad; g.Hp = libisyntax_level_get_height(libisyntax_image_get_level(wsi, 0)) + 2 * g.pad;
    for (int s = 0; s < g.L; s++) { const isyntax_level_t* lv = libisyntax_image_get_level(wsi, s); g.limR.push_back(need(g.Hp, s)); g.limC.push_back(need(g.Wp, s)); g.nbr.push_back(libisyntax_level_get_height_in_tiles(lv)); g.nbc.push_back(libisyntax_level_get_width_in_tiles(lv)); }
    int Wv = (int)(g.Wp - 2 * g.pad), Hv = (int)(g.Hp - 2 * g.pad);
    printf("file %s: %d levels, padded %lldx%lld, valid %dx%d (%.0f Mpx); F=%d D=%d h=%d colour=%s omit S=%d (layers 0..%d) fastpath=%d threads=%d\n", o.file.c_str(), g.L, (long long)g.Wp, (long long)g.Hp, Wv, Hv, (double)Wv * Hv / 1e6, o.F, o.D, o.h, o.trunc_ ? "trunc" : "floor", o.S, o.omit_layers - 1, o.fastpath, o.threads);
    for (int s = 0; s < g.L; s++) if (g.limR[s] > g.nbr[s] * 128 || g.limC[s] > g.nbc[s] * 128) die("tile grid smaller than the coefficient array", s);
    std::vector<std::vector<uint8_t>> EX(g.L);
    for (int s = 0; s < g.L; s++) { EX[s].assign((size_t)g.nbr[s] * g.nbc[s], 0); for (int ty = 0; ty < g.nbr[s]; ty++) for (int tx = 0; tx < g.nbc[s]; tx++) EX[s][(size_t)ty * g.nbc[s] + tx] = libisyntax_tile_exists(is, s, tx, ty) ? 1 : 0; }
    auto tile_exists = [&](int s, i64 ty, i64 tx) { return !(ty < 0 || tx < 0 || ty >= g.nbr[s] || tx >= g.nbc[s]) && EX[s][(size_t)ty * g.nbc[s] + tx] != 0; };
    auto fdiv = [](int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); };
    Window win; win.is = is; win.g = &g; Pool pool(o.threads); Stat st; std::vector<FrameRec> recs; std::mutex recmu; std::atomic<long> peak_blocks{0};
    const int F = o.F, D = o.D, h = o.h;
    {   long n0 = 0; for (auto v : EX[0]) n0 += v != 0; printf("tiles that exist in the file at scale 0: %ld of %zu (non-existent tiles are drawn white by SlideScape and written white here)\n", n0, EX[0].size()); }
    auto exist_state = [&](int k, int i0y, int i0x) -> int {              // 0: no iSyntax tile under the frame exists, 1: all, 2: some
        if (k > g.top) return 1; int n = 0, ne = 0;
        for (int ty = fdiv(i0y, 256); ty <= fdiv(i0y + F - 1, 256); ty++) for (int tx = fdiv(i0x, 256); tx <= fdiv(i0x + F - 1, 256); tx++) { n++; ne += tile_exists(k, ty, tx); }
        return ne == 0 ? 0 : (ne == n ? 1 : 2);
    };
    auto masked_pixel_frame = [&](Engine& e, int c, int k, int i0y, int i0x, int Fs) {      // decoded pixels, white (Y = 255, Co = Cg = 0) outside the array and on non-existent tiles
        Mat m = pixel_frame(e, k, i0y, i0x, Fs);
        if (k <= g.top) for (int ty = fdiv(i0y, 256); ty <= fdiv(i0y + Fs - 1, 256); ty++) for (int tx = fdiv(i0x, 256); tx <= fdiv(i0x + Fs - 1, 256); tx++) if (!tile_exists(k, ty, tx)) {
            int ya = std::max(ty * 256, i0y), yb = std::min(ty * 256 + 256, i0y + Fs), xa = std::max(tx * 256, i0x), xb = std::min(tx * 256 + 256, i0x + Fs);
            for (int y = ya; y < yb; y++) for (int x = xa; x < xb; x++) m.at(y - i0y, x - i0x) = c == 0 ? 255 : 0; }
        return m;
    };
    std::map<int, std::shared_ptr<std::vector<uint8_t>>> ccache; std::mutex ccmu;
    auto const_codestream = [&](const i32* rgb) {                                   // a frame of one colour: RCT of (R,G,B), LL constant, all details 0 (exact for the 5/3 transform)
        int key = (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
        { std::lock_guard<std::mutex> lk(ccmu); auto it = ccache.find(key); if (it != ccache.end()) return it->second; }
        i32 val[3] = {((rgb[0] - 128) + 2 * (rgb[1] - 128) + (rgb[2] - 128)) >> 2, rgb[2] - rgb[1], rgb[0] - rgb[1]};
        std::vector<i32> coef((size_t)F * F * 3); for (int c = 0; c < 3; c++) { Pyr pc = const_pyramid(F, D, val[c]); mallat(pc, F, D, &coef[(size_t)c * F * F]); }
        auto cs = std::make_shared<std::vector<uint8_t>>(); if (!encode_frame(coef.data(), F, *cs)) die("constant frame encode failed");
        std::lock_guard<std::mutex> lk(ccmu); ccache[key] = cs; return cs;
    };
    // Frame grids of all layers (dicom_tile_tests.py axis_grid, changed): every layer's matrix origin must lie on a common position, otherwise viewers (wsidicom: origin threshold 0.01 mm)
    // treat the layers as different pyramids. Sample i of layer k sits at 2^k (i + 3) - 2 (padded level-0 units) and the J2K alignment needs (i_base + 3) % 32 == 0, so the layer
    // origins can be put on the common grid 2^(top+5) M; the remaining spread is only the sample-centre offsets (<= 2^(top-1) units = 4 um at top = 5).
    struct Ax { int i_base = 0, tiles = 0; i64 cols = 0; };
    std::vector<Ax> GY, GX;
    auto compute_grids = [&](int Kmax, std::vector<Ax>& GYo, std::vector<Ax>& GXo) {
        GYo.assign(Kmax + 1, Ax()); GXo.assign(Kmax + 1, Ax());
        for (int axis = 0; axis < 2; axis++) {
            int nvalid = axis ? Wv : Hv; std::vector<i64> ilo(Kmax + 1), ihi(Kmax + 1); i64 M = -1;
            for (int k = 0; k <= Kmax; k++) {
                auto p = [&](i64 i) { return (double)((i64)1 << k) * (double)(i + 3) - 2.0; }; double half = (double)((i64)1 << k) / 2; i64 lo = 0; while (p(lo) + half <= (double)g.pad - 0.5) lo++;
                i64 hi = lo; while (p(hi + 1) - half < (double)(g.pad + nvalid) - 0.5) hi++; ilo[k] = lo; ihi[k] = hi;
                i64 q = (lo + 3) / ((i64)32 << (Kmax - k)); M = M < 0 ? q : std::min(M, q);
            }
            for (int k = 0; k <= Kmax; k++) { Ax a; a.i_base = (int)(((i64)32 << (Kmax - k)) * M - 3); a.cols = ihi[k] - a.i_base + 1; a.tiles = (int)((a.cols + F - 1) / F); (axis ? GXo : GYo)[k] = a; }
        }
    };
    // DICOM-WSI wants layers until the whole image fits one tile: K = first layer whose matrix is <= F x F (the number of layers therefore varies with the slide)
    int K = 0; for (; K < 30; K++) { compute_grids(K, GY, GX); if (GY[K].cols <= F && GX[K].cols <= F) break; } if (K >= 30) die("no single-tile layer found");
    printf("pyramid layers 0..%d (iSyntax stores scales 0..%d%s)\n", K, g.top, K > g.top ? "; the layers above are derived by the 5/3 low-pass" : "");
    const double mpp_um = (libisyntax_get_is_mpp_known(is) && libisyntax_get_mpp_x(is) > 0) ? libisyntax_get_mpp_x(is) : 0.25, unit_mm = mpp_um / 1000.0;
    if (std::fabs(mpp_um - 0.25) > 0.01) printf("NOTE: pixel spacing %.4f um is not the usual 40x 0.25 um\n", mpp_um);
    auto finish_layer = [&](int k, FILE* spool, const std::string& spool_path, std::vector<uint64_t>& offs, std::vector<uint64_t>& lens, long cs_total, int tix, int tiy, i64 mcx, i64 mcy, int ibx, int iby) {
        fclose(spool);
        LayerMeta m; m.layer = k; m.rows = m.cols = F; m.frames_x = tix; m.frames_y = tiy; m.matrix_cols = (long)mcx; m.matrix_rows = (long)mcy; m.spacing_mm = unit_mm * (double)((i64)1 << k); m.resampled = k > 0;
        double ox = libisyntax_image_get_offset_x(wsi), oy = libisyntax_image_get_offset_y(wsi);            // level-0 units (scale factor of the file); sample i of layer k is centred at ox + 2^k (i + 3) - 2
        auto corner = [&](int kk, int ib_x, int ib_y, double& cx, double& cy) { double pxl = (double)((i64)1 << kk); cx = (ox + pxl * (ib_x + 3) - 2.0 - pxl / 2) * unit_mm; cy = (oy + pxl * (ib_y + 3) - 2.0 - pxl / 2) * unit_mm; };
        double Bx, By, B0x, B0y; corner(k, ibx, iby, Bx, By); corner(0, GX[0].i_base, GY[0].i_base, B0x, B0y);
        if (o.clamp_um > 0) { double c = o.clamp_um / 1000.0; Bx = B0x + std::max(-c, std::min(c, Bx - B0x)); By = B0y + std::max(-c, std::min(c, By - B0y)); }      // optional: keep layer origins within a viewer's tolerance (wsidicom: 0.01 mm)
        // Placeholder mapping (owner confirmed the default orientation looks the same in SlideScape and OmniHisto): slide position of pixel (r, c) = (X0 + o0 c s + o3 r s, Y0 + o1 c s + o4 r s);
        // X0, Y0 put the layer-0 matrix in the lower-left of its slide bounding box.
        double W0 = (double)GX[0].cols * unit_mm, H0 = (double)GY[0].cols * unit_mm; const double* oo = o.orient;
        double minx = std::min(0.0, oo[0] * W0) + std::min(0.0, oo[3] * H0), miny = std::min(0.0, oo[1] * W0) + std::min(0.0, oo[4] * H0);
        m.origin_x_mm = Bx - minx; m.origin_y_mm = By - miny; for (int q = 0; q < 6; q++) m.orientation[q] = o.orient[q];
        std::string base = o.file; size_t sl = base.find_last_of("/\\"); if (sl != std::string::npos) base = base.substr(sl + 1); size_t dt = base.rfind('.'); if (dt != std::string::npos) base = base.substr(0, dt);
        m.source_name = base; m.study_uid = make_uid(base, "study"); m.series_uid = make_uid(base, "series"); m.for_uid = make_uid(base, "frame-of-reference"); m.specimen_uid = make_uid(base, "specimen"); m.dimorg_uid = make_uid(base, "dimorg");
        m.sop_uid = make_uid(base, "sop-layer-" + std::to_string(k));
        auto sv = [](const char* v, const char* d) { return std::string(v && *v ? v : d); };
        m.manufacturer = sv(libisyntax_get_manufacturer(is), "UNKNOWN"); m.model = sv(libisyntax_get_manufacturers_model_name(is), "UNKNOWN"); m.serial = sv(libisyntax_get_device_serial_number(is), "UNKNOWN");
        m.software = "libisyntax strip transcoder (fork) 0.1"; m.acquisition_datetime = sv(libisyntax_get_acquisition_datetime(is), ""); m.lossy_method = sv(libisyntax_get_lossy_image_compression_method(is), "");
        m.source_ratio = libisyntax_get_lossy_image_compression_ratio(is); m.our_ratio = (double)m.matrix_cols * m.matrix_rows * 3 / std::max<long>(1, cs_total); m.icc_path = o.icc;
        std::string outp = o.dicom_dir + "/" + base + "_L" + std::to_string(k) + ".dcm", err;
        if (!write_layer_file(m, outp, spool_path, offs, lens, err)) die(("DICOM write failed: " + err).c_str());
        remove(spool_path.c_str()); printf("  wrote %s (%zu frames, %.1f MB of codestreams)\n", outp.c_str(), offs.size(), cs_total / 1e6); fflush(stdout);
    };
    auto spool_row = [&](FILE* spool, std::vector<FrameRec>& local, std::vector<uint64_t>& offs, std::vector<uint64_t>& lens, uint64_t& spool_pos, long& cs_total) {     // frames of one row, in order, as encapsulated items
        for (size_t t2 = 0; t2 < local.size(); t2++) {
            FrameRec& r = local[t2]; if (!r.cs) die("frame without codestream", r.k, r.ty, (long)t2); uint32_t len = (uint32_t)((r.cs->size() + 1) & ~(size_t)1);
            unsigned char hdr[8] = {0xFE, 0xFF, 0x00, 0xE0, (unsigned char)(len & 255), (unsigned char)((len >> 8) & 255), (unsigned char)((len >> 16) & 255), (unsigned char)(len >> 24)};
            fwrite(hdr, 1, 8, spool); fwrite(r.cs->data(), 1, r.cs->size(), spool); if (r.cs->size() & 1) fputc(0, spool);
            offs.push_back(spool_pos); lens.push_back(len); spool_pos += 8 + len; cs_total += (long)r.cs->size(); }
    };
    auto ensure_blocks = [&](const std::set<u64>& need_set) {                  // keep exactly these blocks; decode the missing ones in parallel
        double t = now();
        for (auto it = win.m.begin(); it != win.m.end();) { if (!need_set.count(it->first)) it = win.m.erase(it); else ++it; }
        std::vector<u64> missing; for (u64 key : need_set) { int by = (int)(key & ((1u << 27) - 1)), bx = (int)((key >> 27) & ((1u << 27) - 1)), s2 = (int)(key >> 56);
            if (bx < g.nbc[s2] && by < g.nbr[s2] && !win.m.count(key)) { missing.push_back(key); win.m[key] = std::make_unique<Blk>(); } }
        pool.run((int)missing.size(), [&](int i, int) { u64 key = missing[i]; int by = (int)(key & ((1u << 27) - 1)), bx = (int)((key >> 27) & ((1u << 27) - 1)), kind = (int)((key >> 54) & 3), s2 = (int)(key >> 56);
            Blk* b = win.m[key].get(); int32_t z = 0;
            if (libisyntax_tile_read_coefficients(is, s2, bx, by, kind ? LIBISYNTAX_COEFFICIENTS_H : LIBISYNTAX_COEFFICIENTS_LL, b->d, &z) != LIBISYNTAX_OK) die("block decode failed", s2, bx, by);
            b->zero = z != 0; });
        win.decodes += (long)missing.size(); win.peak = std::max(win.peak, win.m.size()); st.t_decode_ns += (long)((now() - t) * 1e9);
    };
    Mat plane[3]; int plane_layer = -1;                                          // pixel plane of the last derived layer
    for (int k = std::max(0, o.layer_lo); k <= std::min(K, o.layer_hi); k++) {
        int iby = GY[k].i_base, tiy = GY[k].tiles, ibx = GX[k].i_base, tix = GX[k].tiles; i64 mcy = GY[k].cols, mcx = GX[k].cols;
        if (k > g.top) {                                                       // ---- layers above the top stored scale: 5/3 low-pass of the layer below, whole plane in memory (small)
            if (plane_layer != k - 1) {
                if (k - 1 != g.top) die("derived plane chain broken", k);
                int Rt = 2 * g.R(g.top), Ct = 2 * g.C(g.top);
                { Plan pl; plan_pixel_frame(g, pl, g.top, 0, 0, std::max(Rt, Ct)); ensure_blocks(std::set<u64>(pl.blocks.begin(), pl.blocks.end())); }
                Mat tp[3]; for (int c = 0; c < 3; c++) { Engine e(win, c); tp[c] = e.out_region(g.top, 0, Rt, 0, Ct); }
                if (o.verify_every && g.top >= 1) {                              // check: LP of layer top-1 reproduces the stored/decoded layer top in the interior
                    int Rp = 2 * g.R(g.top - 1), Cp = 2 * g.C(g.top - 1); Plan p2; plan_pixel_frame(g, p2, g.top - 1, 0, 0, std::max(Rp, Cp)); ensure_blocks(std::set<u64>(p2.blocks.begin(), p2.blocks.end()));
                    long bad = 0, tot = 0; for (int c = 0; c < 3; c++) { Engine e(win, c); Mat lp = lp_down(e.out_region(g.top - 1, 0, Rp, 0, Cp));
                        for (int y = 8; y < std::min(lp.r, tp[c].r) - 8; y++) for (int x = 8; x < std::min(lp.c, tp[c].c) - 8; x++) { tot++; bad += lp.at(y, x) != tp[c].at(y, x); } }
                    printf("  LP check layer %d -> %d (interior): %ld of %ld samples differ\n", g.top - 1, g.top, bad, tot); }
                for (int c = 0; c < 3; c++) for (int ty = 0; ty * 256 < Rt; ty++) for (int tx = 0; tx * 256 < Ct; tx++) if (!tile_exists(g.top, ty, tx))
                    for (int y = ty * 256; y < std::min(Rt, ty * 256 + 256); y++) for (int x = tx * 256; x < std::min(Ct, tx * 256 + 256); x++) tp[c].at(y, x) = c == 0 ? 255 : 0;
                for (int c = 0; c < 3; c++) plane[c] = tp[c];
                plane_layer = g.top; win.m.clear();
            }
            for (int c = 0; c < 3; c++) plane[c] = lp_down(plane[c]);
            plane_layer = k; int Rk = plane[0].r, Ck = plane[0].c;
            printf("layer %d (derived): plane %d x %d, %d x %d frames (origin %d,%d)\n", k, Ck, Rk, tix, tiy, ibx, iby); fflush(stdout);
            FILE* spool = nullptr; std::string spool_path; std::vector<uint64_t> offs, lens; uint64_t spool_pos = 0; long cs_total = 0;
            if (!o.dicom_dir.empty()) { spool_path = o.dicom_dir + "/.spool_L" + std::to_string(k) + ".tmp"; spool = fopen(spool_path.c_str(), "wb"); if (!spool) die("cannot create spool"); }
            for (int ty = 0; ty < tiy; ty++) {
                std::vector<FrameRec> local(tix);
                pool.run(tix, [&](int tx, int) {
                    int i0y = iby + ty * F, i0x = ibx + tx * F; FrameRec fr{k, ty, tx, i0y, i0x, 'D', 'D', 0, 0, -1}; st.frames++;
                    int y0 = std::max(i0y, 0), x0 = std::max(i0x, 0), y1 = std::min(i0y + F, Rk), x1 = std::min(i0x + F, Ck);
                    if (y1 <= y0 || x1 <= x0) { i32 w[3] = {255, 255, 255}; fr.cs = const_codestream(w); st.white++; }
                    else {
                        Pyr py[3]; for (int c = 0; c < 3; c++) { Mat pix(F, F); if (c == 0) std::fill(pix.v.begin(), pix.v.end(), 255);
                            for (int y = y0; y < y1; y++) memcpy(&pix.v[(size_t)(y - i0y) * F + (x0 - i0x)], &plane[c].v[(size_t)y * Ck + x0], (size_t)(x1 - x0) * 4); py[c] = full_pyramid(pix, D); }
                        colour_convert(py[0], py[1], py[2], o.trunc_); std::vector<i32> coef((size_t)F * F * 3); for (int c = 0; c < 3; c++) mallat(py[c], F, D, &coef[(size_t)c * F * F]);
                        auto cs = std::make_shared<std::vector<uint8_t>>(); if (!encode_frame(coef.data(), F, *cs)) die("encode failed", k, ty, tx); fr.cs = cs; st.encoded++; }
                    fr.bytes = (long)fr.cs->size(); st.bytes += fr.bytes; local[tx] = fr; });
                if (spool) spool_row(spool, local, offs, lens, spool_pos, cs_total);
                for (auto& r : local) { r.cs.reset(); recs.push_back(r); }
            }
            if (spool) finish_layer(k, spool, spool_path, offs, lens, cs_total, tix, tiy, mcx, mcy, ibx, iby);
            continue;
        }
        int R = 2 * g.R(k), C = 2 * g.C(k);
        printf("layer %d: %d x %d frames (origin %d,%d)\n", k, tix, tiy, ibx, iby); fflush(stdout);
        FILE* spool = nullptr; std::string spool_path; std::vector<uint64_t> offs, lens; uint64_t spool_pos = 0; long cs_total = 0;
        if (!o.dicom_dir.empty()) { if (!o.encode) die("--dicom needs encoding"); spool_path = o.dicom_dir + "/.spool_L" + std::to_string(k) + ".tmp"; spool = fopen(spool_path.c_str(), "wb"); if (!spool) die("cannot create spool in --dicom directory"); }
        struct FP { int ty, tx, i0y, i0x; bool edge; int state = 1; Plan plan; };
        auto plan_row = [&](int ty, std::vector<FP>& row) {
            row.assign(tix, FP()); double t = now();
            pool.run(tix, [&](int tx, int) { FP& f = row[tx]; f.ty = ty; f.tx = tx; f.i0y = iby + ty * F; f.i0x = ibx + tx * F;
                f.edge = !(f.i0y >= 8 && f.i0x >= 8 && f.i0y + F <= R - 8 && f.i0x + F <= C - 8);
                f.state = exist_state(k, f.i0y, f.i0x); if (f.state == 2) f.edge = true;                // frames partly over non-existent tiles use the pixel path with a white mask
                if (f.state == 0) { /* wholly over non-existent tiles: white, nothing to decode */ }
                else if (f.edge) plan_pixel_frame(g, f.plan, k, f.i0y, f.i0x, F); else plan_interior_frame(g, f.plan, k, f.i0y, f.i0x, F, D, h); });
            st.t_plan_ns += (long)((now() - t) * 1e9);
        };
        auto ensure = [&](const std::vector<FP>& row) { std::set<u64> ns; for (auto& f : row) ns.insert(f.plan.blocks.begin(), f.plan.blocks.end()); ensure_blocks(ns); };
        std::vector<FP> cur, nxt; int rows = o.max_rows ? std::min(tiy, o.max_rows) : tiy; plan_row(0, cur); ensure(cur);
        for (int ty = 0; ty < rows; ty++) {
            std::vector<FrameRec> local(tix);
            pool.run(tix, [&](int tx, int) {
                FP& f = cur[tx]; FrameRec fr{k, ty, tx, f.i0y, f.i0x, f.edge ? 'E' : 'I', 'G', 0, 0, -1}; st.frames++; (f.edge ? st.edge : st.interior)++;
                if (f.state == 0) { i32 w[3] = {255, 255, 255}; fr.path = 'W'; st.white++; if (o.encode) { fr.cs = const_codestream(w); fr.bytes = (long)fr.cs->size(); st.bytes += fr.bytes; } local[tx] = fr; return; }
                if (f.state == 2) st.masked++;
                Engine eng[3] = {Engine(win, 0), Engine(win, 1), Engine(win, 2)}; Pyr py[3]; bool have_py = false; bool omit = false; i32 ocol[3] = {255, 255, 255};
                double ta = now(); bool is_dump = false; for (auto& d : o.dumps) if (d[0] == k && d[1] == ty && d[2] == tx) is_dump = true;
                i32 cvv[3] = {0, 0, 0};
                auto try_constant = [&](i32* v) -> bool {                              // exact: all footprint H blocks empty and the top LL constant => the frame is constant
                    if (f.edge || !o.fastpath || f.plan.clipped) return false;
                    for (u64 key : f.plan.blocks) { if (((key >> 54) & 3) != 1) continue; int s2 = (int)(key >> 56), by = (int)(key & ((1u << 27) - 1)), bx = (int)((key >> 27) & ((1u << 27) - 1)); const Blk* b = win.find(s2, 1, bx, by); if (b && !b->zero) return false; }
                    for (int c = 0; c < 3; c++) { bool first = true; for (auto& rc : f.plan.llrects) { Mat m = fetch(win, c, 0, 0, rc[0], rc[1], rc[2], rc[3], rc[4]); for (i32 x : m.v) { if (first) { v[c] = x; first = false; } else if (x != v[c]) return false; } } }
                    return true;
                };
                // ---- omission test (layers 0 and 1)
                bool omit_layer = k < o.omit_layers && o.S >= 0; bool candidate = false; Mat pix[3]; bool have_pix = false;
                if (omit_layer && !f.edge && try_constant(cvv)) { omit = true; candidate = true; st.prefilter_pass++; st.const_omit++; if (o.audit) st.audit_omit++;
                    { i32 t2 = cvv[0] - qhalf(cvv[2], true), G2 = t2 + cvv[2], B2 = t2 - qhalf(cvv[1], true), R2 = B2 + cvv[1]; ocol[0] = std::min(255, std::max(0, R2)); ocol[1] = std::min(255, std::max(0, G2)); ocol[2] = std::min(255, std::max(0, B2)); } }   // constant frame: strict range 0, no pixel decode
                else if (omit_layer && !f.edge) {
                    int mx[3] = {0, 0, 0};                                            // interior H maxima over the frame's own scales
                    for (int c = 0; c < 3; c++) for (int j = 0; j < D && k + j <= g.top; j++) { int s = k + j, nn = F >> (j + 1), hx = ((f.i0x + 1 + 2) >> (j + 1)) - 1, hy = ((f.i0y + 1 + 2) >> (j + 1)) - 1;
                        Mat a = fetch(win, c, 1, 0, s, hy - 1, hy + nn - 1, hx, hx + nn), b = fetch(win, c, 1, 1, s, hy, hy + nn, hx - 1, hx + nn - 1), d = fetch(win, c, 1, 2, s, hy, hy + nn, hx, hx + nn);
                        for (auto* m : {&a, &b, &d}) for (i32 v : m->v) mx[c] = std::max(mx[c], std::abs(v)); }
                    candidate = mx[0] <= o.prefilter && mx[1] <= o.prefilter && mx[2] <= o.prefilter; if (candidate) st.prefilter_pass++;
                    if (o.audit && !candidate) { /* exact check below for the audit */ }
                    if (candidate || o.audit) {
                        double tb = now(); for (int c = 0; c < 3; c++) pix[c] = masked_pixel_frame(eng[c], c, k, f.i0y, f.i0x, F); have_pix = true;
                        int lo[3] = {1 << 30, 1 << 30, 1 << 30}, hi[3] = {-(1 << 30), -(1 << 30), -(1 << 30)}; long sm[3] = {0, 0, 0};
                        for (size_t i = 0; i < pix[0].v.size(); i++) { i32 Y = pix[0].v[i], Co = pix[1].v[i], Cg = pix[2].v[i]; i32 t = Y - qhalf(Cg, true), G = t + Cg, B = t - qhalf(Co, true), Rr = B + Co;
                            i32 rgbv[3] = {std::min(255, std::max(0, Rr)), std::min(255, std::max(0, G)), std::min(255, std::max(0, B))}; for (int q = 0; q < 3; q++) { lo[q] = std::min(lo[q], rgbv[q]); hi[q] = std::max(hi[q], rgbv[q]); sm[q] += rgbv[q]; } }
                        int range = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]}); bool om = range <= o.S;
                        if (o.audit) { if (om) { st.audit_omit++; if (!candidate) st.audit_false_reject++; for (int c = 0; c < 3; c++) { int cur_max = st.max_audit_abs[c]; while (mx[c] > cur_max && !st.max_audit_abs[c].compare_exchange_weak(cur_max, mx[c])) {} } } }
                        omit = om; { long n = (long)pix[0].v.size(); for (int q = 0; q < 3; q++) ocol[q] = (int)((sm[q] + n / 2) / n); } st.t_omit_ns += (long)((now() - tb) * 1e9);
                    }
                }
                if (omit) { fr.omitted = 1; st.omitted++; fr.path = 'O'; if (!is_dump) { if (o.encode) { fr.cs = const_codestream(ocol); fr.bytes = (long)fr.cs->size(); st.bytes += fr.bytes; } local[tx] = fr; return; } }
                // ---- coefficients: fast path / windowed builder / pixel path
                double tc = now();
                if (!is_dump && try_constant(cvv)) {
                    for (int c = 0; c < 3; c++) py[c] = const_pyramid(F, D, cvv[c]); have_py = true; fr.path = 'F'; st.fast++; st.t_fast_ns += (long)((now() - tc) * 1e9);
                    if (o.verify_every && (st.fast % o.verify_every == 0)) { st.fast_verified++; for (int c = 0; c < 3; c++) { Pyr ref = build_frame(eng[c], k, f.i0y, f.i0x, F, D, h); if (!same(ref, py[c])) { st.fast_bad++; break; } } }
                }
                if (!have_py && f.edge && f.plan.blocks.empty()) { for (int c = 0; c < 3; c++) py[c] = const_pyramid(F, D, c == 0 ? 255 : 0); have_py = true; fr.path = 'F'; st.fast++; }   // padding frame: white
                if (!have_py) {
                    if (f.edge || have_pix) { double tp = now(); if (!have_pix) for (int c = 0; c < 3; c++) pix[c] = masked_pixel_frame(eng[c], c, k, f.i0y, f.i0x, F); for (int c = 0; c < 3; c++) py[c] = full_pyramid(pix[c], D); fr.path = 'P'; st.t_pixel_ns += (long)((now() - tp) * 1e9); }
                    else { double tb2 = now(); for (int c = 0; c < 3; c++) py[c] = build_frame(eng[c], k, f.i0y, f.i0x, F, D, h); st.t_build_ns += (long)((now() - tb2) * 1e9); }
                    have_py = true;
                }
                if (o.verify_every && fr.path == 'G' && ((long)(st.interior) % o.verify_every == 0)) {       // exactness: builder == full forward transform of the decoded pixels
                    double tv = now(); bool bad = false; for (int c = 0; c < 3; c++) { Mat px = have_pix ? pix[c] : masked_pixel_frame(eng[c], c, k, f.i0y, f.i0x, F); if (!same(full_pyramid(px, D), py[c])) bad = true; }
                    st.verified++; if (bad) st.verify_bad++; st.t_verify_ns += (long)((now() - tv) * 1e9);
                }
                colour_convert(py[0], py[1], py[2], o.trunc_);
                std::vector<i32> coef((size_t)F * F * 3); for (int c = 0; c < 3; c++) mallat(py[c], F, D, &coef[(size_t)c * F * F]);
                if (o.oob) {                                                   // pixels the stock decoder would have to clamp (outside 0..255 after the inverse of our coefficients)
                    Mat Yp = inv_pyramid(unmallat(&coef[0], F, D)), Cb = inv_pyramid(unmallat(&coef[(size_t)F * F], F, D)), Cr = inv_pyramid(unmallat(&coef[(size_t)2 * F * F], F, D)); long bad = 0;
                    for (size_t i = 0; i < (size_t)F * F; i++) { i32 G = Yp.v[i] - ((Cb.v[i] + Cr.v[i]) >> 2), Rr = Cr.v[i] + G, B = Cb.v[i] + G; i32 e[3] = {Rr + 128, G + 128, B + 128}; bool b = false; for (int c = 0; c < 3; c++) b |= e[c] < 0 || e[c] > 255; bad += b; }
                    fr.oob = bad; }
                if (is_dump) { char fn[1024]; snprintf(fn, sizeof fn, "%s/cpp_L%d_%d_%d.coef", o.dumpdir.c_str(), k, ty, tx); FILE* ff = fopen(fn, "wb"); fwrite(coef.data(), 4, coef.size(), ff); fclose(ff); }
                if (o.encode) { double te = now(); std::vector<uint8_t> cs; if (!encode_frame(coef.data(), F, cs)) die("encode failed", k, ty, tx); fr.bytes = (long)cs.size(); fr.cs = std::make_shared<std::vector<uint8_t>>(cs); st.bytes += fr.bytes; st.encoded++; st.t_encode_ns += (long)((now() - te) * 1e9);
                    if (o.roundtrip && (st.encoded % std::max(1, o.verify_every ? o.verify_every : 50) == 0)) {          // stock decoder vs the inverse of our own coefficients (catches Tier-1 range overflow)
                        std::vector<i32> dec; st.rt_checked++; bool bad = !decode_codestream(cs, F, dec);
                        if (!bad) { Mat Yp = inv_pyramid(unmallat(&coef[0], F, D)), Cb = inv_pyramid(unmallat(&coef[(size_t)F * F], F, D)), Cr = inv_pyramid(unmallat(&coef[(size_t)2 * F * F], F, D));
                            for (size_t i = 0; i < (size_t)F * F && !bad; i++) { i32 G = Yp.v[i] - ((Cb.v[i] + Cr.v[i]) >> 2), Rr = Cr.v[i] + G, B = Cb.v[i] + G; i32 e[3] = {Rr, G, B};
                                for (int c = 0; c < 3; c++) { i32 ex = std::min(255, std::max(0, e[c] + 128)); if (dec[(size_t)c * F * F + i] != ex) { bad = true; break; } } } }
                        if (bad) st.rt_bad++; } }
                local[tx] = fr; (void)ta;
            });
            if (spool) spool_row(spool, local, offs, lens, spool_pos, cs_total);
            for (auto& r : local) r.cs.reset();
            { std::lock_guard<std::mutex> lk(recmu); for (auto& r : local) recs.push_back(r); }
            if (ty + 1 < rows) { plan_row(ty + 1, nxt); ensure(nxt); cur.swap(nxt); }
            if (ty % 50 == 49 || ty + 1 == rows) { printf("  layer %d row %d/%d: %.0fs, window %zu blocks (%.0f MB), decodes %ld, frames %ld, bytes %.1f MB\n", k, ty + 1, rows, now() - t0, win.m.size(), win.m.size() * sizeof(Blk) / 1e6, win.decodes.load(), st.frames.load(), st.bytes.load() / 1e6); fflush(stdout); }
        }
        if (spool) { if (rows != tiy) die("--dicom needs all frame rows (no --max-rows)"); finish_layer(k, spool, spool_path, offs, lens, cs_total, tix, tiy, mcx, mcy, ibx, iby); }
        win.m.clear();
    }
    double tt = now() - t0;
    long rss_kb = 0; { FILE* f = fopen("/proc/self/status", "r"); char line[256]; while (f && fgets(line, sizeof line, f)) if (!strncmp(line, "VmHWM:", 6)) rss_kb = atol(line + 6); if (f) fclose(f); }
    printf("\nSUMMARY %s (%s, %d threads)\n", o.file.c_str(), o.trunc_ ? "trunc" : "floor", o.threads);
    printf("  frames %ld (interior %ld, edge %ld): fast-path constant %ld (+%ld omitted without a pixel decode), omitted %ld (prefilter passed %ld), encoded %ld, coded bytes %.2f MB\n", st.frames.load(), st.interior.load(), st.edge.load(), st.fast.load(), st.const_omit.load(), st.omitted.load(), st.prefilter_pass.load(), st.encoded.load(), st.bytes.load() / 1e6);
    printf("  frames over non-existent iSyntax tiles written white: %ld; frames partly over them (pixel path with white mask): %ld\n", st.white.load(), st.masked.load());
    printf("  block decodes (blocks, 3 colours each): %ld, peak window %zu blocks (%.0f MB), process peak RSS %.2f GB\n", win.decodes.load(), win.peak, win.peak * sizeof(Blk) / 1e6, rss_kb / 1e6);
    printf("  exactness: builder vs full forward transform of decoded pixels: %ld frames checked, %ld mismatches; fast-path vs builder: %ld checked, %ld mismatches; Tier-1 round trip: %ld checked, %ld bad\n", st.verified.load(), st.verify_bad.load(), st.fast_verified.load(), st.fast_bad.load(), st.rt_checked.load(), st.rt_bad.load());
    if (o.audit) printf("  omission audit (all layer-0/1 interior frames exact-tested): %ld omittable, %ld of them rejected by the prefilter (T=%d); max |H| seen in omittable frames per colour: %d %d %d\n", st.audit_omit.load(), st.audit_false_reject.load(), o.prefilter, st.max_audit_abs[0].load(), st.max_audit_abs[1].load(), st.max_audit_abs[2].load());
    printf("  wall %.1f s; thread-seconds: plan %.1f, block decode (wall in phase) %.1f, builder %.1f, pixel path %.1f, omission pixel test %.1f, fast path %.2f, Tier-1 %.1f, verify %.1f\n", tt, st.t_plan_ns / 1e9, st.t_decode_ns / 1e9, st.t_build_ns / 1e9, st.t_pixel_ns / 1e9, st.t_omit_ns / 1e9, st.t_fast_ns / 1e9, st.t_encode_ns / 1e9, st.t_verify_ns / 1e9);
    if (!o.csv.empty()) { FILE* f = fopen(o.csv.c_str(), "w"); fprintf(f, "layer,ty,tx,i0y,i0x,kind,path,bytes,oob\n"); std::sort(recs.begin(), recs.end(), [](const FrameRec& a, const FrameRec& b) { return std::tie(a.k, a.ty, a.tx) < std::tie(b.k, b.ty, b.tx); });
        for (auto& r : recs) { fprintf(f, "%d,%d,%d,%d,%d,%c,%c,%ld,%ld\n", r.k, r.ty, r.tx, r.i0y, r.i0x, r.kind, r.path, r.bytes, r.oob); } fclose(f); }
    return (st.verify_bad || st.fast_bad || st.rt_bad || st.audit_false_reject) ? 3 : 0;
}
