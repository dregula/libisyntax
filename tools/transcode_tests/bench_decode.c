#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "libisyntax.h"
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
int main(int argc,char**argv){
  libisyntax_init(); isyntax_t* s; if(libisyntax_open(argv[1],0,&s)!=LIBISYNTAX_OK) return 1;
  const isyntax_image_t* wsi=libisyntax_get_wsi_image(s); int L=libisyntax_image_get_level_count(wsi);
  int bw=libisyntax_get_block_width(s), bh=libisyntax_get_block_height(s); int16_t* buf=malloc(9*bw*bh*2); long sink=0; double tot=0;
  for(int lv=0;lv<L;lv++){
    const isyntax_level_t* l=libisyntax_image_get_level(wsi,lv); int w=libisyntax_level_get_width_in_tiles(l), h=libisyntax_level_get_height_in_tiles(l);
    double t0=now(); long nz_tiles=0;
    for(int pass=0;pass<3;pass++) for(int y=0;y<h;y++) for(int x=0;x<w;x++){ int32_t z; libisyntax_tile_read_coefficients(s,lv,x,y,LIBISYNTAX_COEFFICIENTS_H,buf,&z); sink+=buf[5]; if(pass==0&&!z) nz_tiles++; }
    double dt=(now()-t0)/3; tot+=dt; printf("level %d: %d tiles, %ld with data, decode(+read) %.3f s per pass\n",lv,w*h,nz_tiles,dt);
  }
  printf("all levels, H coefficients, one core: %.3f s (sink %ld)\n",tot,sink); return 0; }
