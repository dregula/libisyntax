// Example (added in this fork): read the stored wavelet coefficients of tiles with libisyntax_tile_read_coefficients(),
// print per-band statistics, optionally dump a tile as raw int16, and check that reading many tiles from several threads
// gives exactly the same result as reading them one by one.
//
//   isyntax_coefficients FILE LEVEL TILE_X TILE_Y [h|ll] [dump.raw]      one tile: stats (+ raw int16 dump)
//   isyntax_coefficients FILE LEVEL all THREADS                           every H tile of LEVEL, serial vs THREADS threads
//
// Raw dump layout: h = int16 [color Y/Co/Cg][band HL/LH/HH][block_h][block_w]; ll = int16 [color][block_h][block_w].
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "libisyntax.h"

static const char* color_name[3] = {"Y", "Co", "Cg"};
static const char* band_name[3] = {"HL", "LH", "HH"};

static void stats(const int16_t* p, size_t n, long* nz, int* minabs, int* maxabs) {
    *nz = 0; *minabs = 1 << 20; *maxabs = 0;
    for (size_t i = 0; i < n; ++i) {
        if (p[i] == 0) continue;
        int a = p[i] < 0 ? -p[i] : p[i];
        (*nz)++; if (a < *minabs) *minabs = a; if (a > *maxabs) *maxabs = a;
    }
    if (*nz == 0) *minabs = 0;
}

typedef struct { isyntax_t* isyntax; int level; int64_t w, h; int64_t first, step; uint64_t* sums; int failed; } job_t;

static uint64_t fnv(const int16_t* p, size_t n) {
    uint64_t h = 1469598103934665603ull;
    for (size_t i = 0; i < n; ++i) { h ^= (uint16_t)p[i]; h *= 1099511628211ull; }
    return h;
}

static void* worker(void* arg) {
    job_t* j = (job_t*)arg;
    size_t n = 9 * (size_t)libisyntax_get_block_width(j->isyntax) * libisyntax_get_block_height(j->isyntax);
    int16_t* buf = malloc(n * sizeof(int16_t));
    for (int64_t t = j->first; t < j->w * j->h; t += j->step) {
        int32_t zero;
        if (libisyntax_tile_read_coefficients(j->isyntax, j->level, t % j->w, t / j->w, LIBISYNTAX_COEFFICIENTS_H, buf, &zero) != LIBISYNTAX_OK) j->failed++;
        j->sums[t] = fnv(buf, n);
    }
    free(buf);
    return NULL;
}

static void* thread_main(void* arg) {
    worker(arg);
    libisyntax_thread_shutdown();   // release this thread's 16 MB scratch memory (not on the main thread)
    return NULL;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s FILE LEVEL TILE_X TILE_Y [h|ll] [dump.raw]\n       %s FILE LEVEL all THREADS\n", argv[0], argv[0]); return 1; }
    libisyntax_init();
    isyntax_t* isyntax;
    if (libisyntax_open(argv[1], 0, &isyntax) != LIBISYNTAX_OK) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    int level = atoi(argv[2]);
    const isyntax_image_t* wsi = libisyntax_get_wsi_image(isyntax);
    const isyntax_level_t* lvl = libisyntax_image_get_level(wsi, level);
    int64_t w = libisyntax_level_get_width_in_tiles(lvl), h = libisyntax_level_get_height_in_tiles(lvl);
    int bw = libisyntax_get_block_width(isyntax), bh = libisyntax_get_block_height(isyntax);
    printf("level %d: %lld x %lld tiles, block %d x %d, levels %d\n", level, (long long)w, (long long)h, bw, bh, libisyntax_image_get_level_count(wsi));
    size_t plane = (size_t)bw * bh;
    if (strcmp(argv[3], "all") == 0) {
        int threads = argc > 4 ? atoi(argv[4]) : 4; if (threads < 1) threads = 1;
        uint64_t* serial = calloc(w * h, sizeof(uint64_t)); uint64_t* par = calloc(w * h, sizeof(uint64_t));
        job_t js = {isyntax, level, w, h, 0, 1, serial, 0};
        worker(&js);                                           // serial reference (on the main thread)
        pthread_t* tid = malloc(threads * sizeof(pthread_t)); job_t* jobs = malloc(threads * sizeof(job_t));
        for (int i = 0; i < threads; ++i) { jobs[i] = (job_t){isyntax, level, w, h, i, threads, par, 0}; pthread_create(&tid[i], NULL, thread_main, &jobs[i]); }
        for (int i = 0; i < threads; ++i) pthread_join(tid[i], NULL);
        long bad = 0, failed = js.failed; for (int i = 0; i < threads; ++i) failed += jobs[i].failed;
        for (int64_t t = 0; t < w * h; ++t) bad += serial[t] != par[t];
        printf("%lld tiles, %d threads: %ld tiles differ between serial and parallel, %ld read errors\n", (long long)(w * h), threads, bad, failed);
        return (bad || failed) ? 2 : 0;
    }
    int64_t tx = atoll(argv[3]), ty = argc > 4 ? atoll(argv[4]) : 0;
    int want_ll = argc > 5 && strcmp(argv[5], "ll") == 0;
    size_t n = (want_ll ? 3 : 9) * plane;
    int16_t* buf = malloc(n * sizeof(int16_t)); int32_t zero = 0;
    isyntax_error_t err = libisyntax_tile_read_coefficients(isyntax, level, tx, ty, want_ll ? LIBISYNTAX_COEFFICIENTS_LL : LIBISYNTAX_COEFFICIENTS_H, buf, &zero);
    if (err != LIBISYNTAX_OK) { fprintf(stderr, "libisyntax_tile_read_coefficients failed: %d\n", err); return 1; }
    printf("tile (%lld,%lld) %s: all-zero/empty=%d\n", (long long)tx, (long long)ty, want_ll ? "LL" : "H", zero);
    for (int c = 0; c < 3; ++c) for (int b = 0; b < (want_ll ? 1 : 3); ++b) {
        long nz; int mn, mx; stats(buf + ((size_t)c * (want_ll ? 1 : 3) + b) * plane, plane, &nz, &mn, &mx);
        printf("  %-2s %-2s: nonzero %6ld of %zu (%.3f%%), min|c| %d, max|c| %d\n", color_name[c], want_ll ? "LL" : band_name[b], nz, plane, 100.0 * nz / plane, mn, mx);
    }
    if (argc > 6) { FILE* f = fopen(argv[6], "wb"); if (f) { fwrite(buf, sizeof(int16_t), n, f); fclose(f); printf("wrote %zu bytes to %s\n", n * 2, argv[6]); } }
    free(buf); libisyntax_close(isyntax);
    return 0;
}
