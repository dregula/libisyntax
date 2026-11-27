#pragma once

#include "isyntax.h"
#include "libisyntax.h"
#include "platform_mutex.h"

typedef struct isyntax_tile_list_t {
    isyntax_tile_t* head;
    isyntax_tile_t* tail;
    int count;
    const char* dbg_name;
} isyntax_tile_list_t;

typedef struct isyntax_cache_t {
    isyntax_tile_list_t cache_list;
    platform_mutex_t mutex;
    // TODO(avirodov): int refcount;
    int target_cache_size;
    block_allocator_t* ll_coeff_block_allocator;
    block_allocator_t* h_coeff_block_allocator;
	bool is_block_allocator_owned;
    int allocator_block_width;
    int allocator_block_height;
} isyntax_cache_t;

// TODO(avirodov): can this ever fail?
void isyntax_tile_read(isyntax_t* isyntax, isyntax_cache_t* cache, int scale, int tile_x, int tile_y,
                       uint32_t* pixels_buffer, enum isyntax_pixel_format_t pixel_format);

// Raw (stored) wavelet coefficients of one tile, WITHOUT any inverse transform and WITHOUT colour conversion.
// Added in this fork. Reads and Hulsken-decodes the codeblocks of the three colour channels (Y, Co, Cg) straight from the
// file into the caller's buffer. Uses no cache, no shared mutable state and no lock, so any number of threads may call it
// concurrently on the same isyntax_t (the calling thread needs thread memory: see libisyntax_thread_init()).
//   want_ll == false: H coefficients, layout out[color][band][block_height][block_width], band 0 = HL, 1 = LH, 2 = HH
//                     (3 * 3 * block_w * block_h i16). Exists for every scale.
//   want_ll == true:  LL coefficients, layout out[color][block_height][block_width] (3 * block_w * block_h i16). Only the
//                     top scale (wsi->max_scale) stores LL; every other scale's LL is the output of its parent's inverse transform.
// *out_is_all_zero is set when every codeblock of the tile is an 'empty' block (<= 8 bytes) or the tile does not exist.
// Returns false on a decode or I/O error, or if want_ll is requested for a scale that does not store LL.
bool isyntax_read_tile_coefficients(isyntax_t* isyntax, int scale, int tile_x, int tile_y, bool want_ll,
                                    i16* out_coefficients, bool* out_is_all_zero);

void tile_list_init(isyntax_tile_list_t* list, const char* dbg_name);
void tile_list_remove(isyntax_tile_list_t* list, isyntax_tile_t* tile);
