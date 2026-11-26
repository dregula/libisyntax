/* Test of the OpenJPEG coefficient-input patch (ext/openjpeg-2.5.4, patches/openjpeg-2.5.4-coefficient-input.patch).
 *
 *   coeff_input_test <dir> <name> [<name>...]
 *
 * For each <name> it reads <dir>/<name>_rgb.raw (256x256x3 uint8) and <dir>/<name>_coef.raw (3 x 256x256 int32: DC-shifted,
 * RCT-transformed, 5-level 5/3 coefficients in OpenJPEG's Mallat layout; see tools/v1_analysis/make_coeff_test_frames.py), then
 *   A: encodes the RGB with stock OpenJPEG settings (raw J2K codestream, reversible 5/3, RCT, 6 resolutions, 64x64 code-blocks),
 *   B: encodes the coefficients with opj_set_coefficient_input(), same parameters,
 * and requires the two codestreams to be byte-identical, then decodes B with the stock decoder and requires the RGB back.
 * Exit code 0 only if every frame passes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openjpeg.h"
#define F 256
static unsigned char *slurp(const char *p, size_t *n) {
    FILE *f = fopen(p, "rb"); if (!f) { perror(p); exit(2); }
    fseek(f, 0, SEEK_END); *n = (size_t)ftell(f); fseek(f, 0, SEEK_SET);
    unsigned char *b = malloc(*n); if (fread(b, 1, *n, f) != *n) { perror("read"); exit(2); } fclose(f); return b;
}
static void quiet(const char *m, void *d) { (void)d; fputs(m, stderr); }
/* encode with the tile-write API; coef != NULL selects the coefficient-input mode. Returns codestream (malloc) and its size. */
static unsigned char *encode(const unsigned char *rgb, const int *coef, size_t *out_n, const char *tmp) {
    opj_cparameters_t p; opj_set_default_encoder_parameters(&p);
    p.tcp_numlayers = 1; p.cp_disto_alloc = 1; p.tcp_rates[0] = 0; p.irreversible = 0; p.tcp_mct = 1;
    p.numresolution = 6; p.cblockw_init = 64; p.cblockh_init = 64; p.cod_format = 0;
    opj_image_cmptparm_t cp[3]; memset(cp, 0, sizeof cp);
    for (int c = 0; c < 3; c++) { cp[c].dx = cp[c].dy = 1; cp[c].w = cp[c].h = F; cp[c].prec = 8; cp[c].sgnd = 0; }
    opj_image_t *im = opj_image_create(3, cp, OPJ_CLRSPC_SRGB); im->x0 = im->y0 = 0; im->x1 = im->y1 = F;
    opj_codec_t *cd = opj_create_compress(OPJ_CODEC_J2K); opj_set_error_handler(cd, quiet, NULL);
    if (!opj_setup_encoder(cd, &p, im)) { fprintf(stderr, "setup failed\n"); exit(2); }
    if (coef && !opj_set_coefficient_input(cd, 1)) { fprintf(stderr, "coefficient input failed\n"); exit(2); }
    opj_stream_t *st = opj_stream_create_default_file_stream(tmp, 0);
    if (!opj_start_compress(cd, im, st)) { fprintf(stderr, "start failed\n"); exit(2); }
    size_t npix = (size_t)F * F; unsigned char *data; size_t dn;
    if (coef) { dn = npix * 3 * 4; data = malloc(dn); memcpy(data, coef, dn); }
    else {  /* planar bytes: component after component (the tile-write API wants planar data) */
        dn = npix * 3; data = malloc(dn);
        for (size_t i = 0; i < npix; i++) for (int c = 0; c < 3; c++) data[c * npix + i] = rgb[i * 3 + c];
    }
    if (!opj_write_tile(cd, 0, data, (OPJ_UINT32)dn, st)) { fprintf(stderr, "write_tile failed\n"); exit(2); }
    if (!opj_end_compress(cd, st)) { fprintf(stderr, "end failed\n"); exit(2); }
    opj_stream_destroy(st); opj_destroy_codec(cd); opj_image_destroy(im); free(data);
    return slurp(tmp, out_n);
}
static int decode_check(const unsigned char *cs, size_t n, const unsigned char *rgb, const char *tmp) {
    FILE *f = fopen(tmp, "wb"); fwrite(cs, 1, n, f); fclose(f);
    opj_dparameters_t dp; opj_set_default_decoder_parameters(&dp);
    opj_stream_t *st = opj_stream_create_default_file_stream(tmp, 1); opj_codec_t *cd = opj_create_decompress(OPJ_CODEC_J2K);
    opj_setup_decoder(cd, &dp); opj_image_t *im = NULL;
    if (!opj_read_header(st, cd, &im) || !opj_decode(cd, st, im) || im->numcomps != 3) { fprintf(stderr, "decode failed\n"); return 0; }
    long bad = 0, shown = 0;
    for (size_t i = 0; i < (size_t)F * F; i++) for (int c = 0; c < 3; c++) {
        if (im->comps[c].data[i] != rgb[i * 3 + c]) { bad++; if (shown++ < 5) fprintf(stderr, "  mismatch comp %d y %zu x %zu: decoded %d, expected %d\n", c, i / F, i % F, im->comps[c].data[i], rgb[i * 3 + c]); }
    }
    if (bad) fprintf(stderr, "  %ld of %d samples differ\n", bad, F * F * 3);
    opj_stream_destroy(st); opj_destroy_codec(cd); opj_image_destroy(im); return bad == 0;
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <dir> <name>...\n", argv[0]); return 2; }
    const char *tmpd = getenv("TMPDIR"); char tmp[1024]; snprintf(tmp, sizeof tmp, "%s/coeff_test.j2k", tmpd ? tmpd : "/tmp");
    int fail = 0;
    for (int a = 2; a < argc; a++) {
        char p1[1024], p2[1024]; size_t n1, n2;
        snprintf(p1, sizeof p1, "%s/%s_rgb.raw", argv[1], argv[a]); snprintf(p2, sizeof p2, "%s/%s_coef.raw", argv[1], argv[a]);
        unsigned char *rgb = slurp(p1, &n1); int *coef = (int *)slurp(p2, &n2);
        if (n1 != (size_t)F * F * 3 || n2 != (size_t)F * F * 3 * 4) { fprintf(stderr, "%s: bad file size\n", argv[a]); return 2; }
        size_t na, nb; unsigned char *A = encode(rgb, NULL, &na, tmp), *B = encode(rgb, coef, &nb, tmp);
        int same = na == nb && memcmp(A, B, na) == 0, rt = decode_check(B, nb, rgb, tmp);
        printf("%-8s stock %zu bytes, coefficient-input %zu bytes: %s; decoded back to the RGB: %s\n", argv[a], na, nb, same ? "IDENTICAL" : "DIFFERENT", rt ? "yes" : "NO");
        fail |= !(same && rt); free(A); free(B); free(rgb); free(coef);
    }
    printf(fail ? "FAIL\n" : "PASS\n"); return fail;
}
