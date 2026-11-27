#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "libisyntax.h"
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
int main(int argc,char**argv){
  libisyntax_init(); isyntax_t* s; if(libisyntax_open(argv[1],0,&s)!=LIBISYNTAX_OK) return 1; isyntax_cache_t* c; libisyntax_cache_create("b",20000,&c); libisyntax_cache_inject(c,s);
  const isyntax_image_t* wsi=libisyntax_get_wsi_image(s); const isyntax_level_t* l=libisyntax_image_get_level(wsi,0);
  uint32_t* px=malloc(256*256*4); int x0=40,x1=80,y0=30,y1=50; int n=0; double t0=now(); long sink=0;
  for(int y=y0;y<y1;y++) for(int x=x0;x<x1;x++){ libisyntax_tile_read(s,c,0,x,y,px,LIBISYNTAX_PIXEL_FORMAT_RGBA); sink+=px[1000]; n++; }
  double dt=now()-t0; printf("level 0, %d tiles (cache 20000): %.3f s total, %.3f ms per tile (decode of the touched blocks + IDWT chain + RGB), sink %ld\n",n,dt,dt/n*1e3,sink); return 0; }
