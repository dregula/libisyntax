/*
  BSD 2-Clause License

  Copyright (c) 2019-2026, Pieter Valkema

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
     list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// =====================================================================================================
// ORIENTATION NOTES (added to this fork during a bit-level study of the format; COMMENTS ONLY, no code
// was changed). Companion documents in the repository root: CLAUDE.md, Example_iSyntax_File_Claude.md,
// Example_FIC_Files_Claude.md, Philips_Compression_Code_Samples_Claude.md, Image_Format_Code_Samples_Claude.md,
// Viewing_Code_Samples_Claude.md, Third_Party_Documentation.md. Statements marked [real file] were measured on
// example_isyntax_files/1.isyntax (data model v1, 25600 x 21504 padded, 6 levels, 33867 codeblocks).
//
// FILE LAYOUT
//   XML header | EOT 0D 0A 04 | seektable (v1 only) | FFFE,E0DD delimiter + empty FFFE,E000 item | codeblocks
//   - Every codeblock is [FFFE,E000 + u32 size = 8-byte DICOM item header][payload], back to back, so the
//     codeblock section is one contiguous run that ends exactly at end-of-file [real file]. Offsets and sizes
//     in the tables always refer to the payload, not to the 8-byte item header.
//   - Data model v1 (UFS interface version < 100, e.g. 5.0): the XML holds UFS_IMAGE_BLOCK_HEADER_TABLE
//     (0x301D,0x2014; base64; coordinates, template id and, in the "full" form, offset and size) and the file
//     holds a seektable (tag 0x301D,0x2015; 40-byte entries) with the offset and size of every codeblock,
//     indexed by a "block id" computed from the block coordinates (see isyntax_open()).
//   - Data model v2 (>= 100, e.g. 100.4): no seektable. The XML holds UFS_IMAGE_CLUSTER_HEADER_TEMPLATES
//     (0x2016) and UFS_IMAGE_CLUSTER_HEADER_TABLE (0x301D,0x201F; base64) which lists, per cluster ("chunk"),
//     its coordinates, template id, data offset and size, and the relative offset and size of each codeblock.
//   - The v1 seektable and the v2 cluster table describe identical data: for the same scan, all 33867
//     (colour, scale, coefficient, offset, size) tuples are equal [real file vs 1.fic].
//   - Philips' own Python/MATLAB samples understand only v1; v2 is documented only by this code.
//
// GEOMETRY (all derived from the 5/3 filter support; verified against Philips' MATLAB and on the real file)
//   - scale s = 0 is full resolution. A codeblock holds 128 x 128 coefficients, so it covers a tile of
//     256 * 2^s pixels per side after the inverse transform. Block coordinates are in level-0 pixels.
//   - Padded image coordinates: the image is zero-padded by 3 * (2^L - 1) = (3 << L) - 3 pixels on every side
//     (L = number of levels), and further on the right/bottom so the grid is a multiple of the block size.
//   - First valid coefficient position at scale s: H = (3 << s) - 2 (1, 4, 10, 22, 46, 94, ...), LL = H + (1 << s).
//     At scale s, H coefficient n sits at position H + 2^(s+1) * n and LL coefficient n at LL + 2^(s+1) * n.
//     Hence at every level the stream starts with an H sample and the low band lies on odd positions of the
//     interleaved output ("cas = 1" in isyntax_dwt.c).
//   - Valid size = padded size - 2 * padding. The valid-data envelope vertices in the XML start at the padding
//     (189 for 6 levels, 765 for 8, 1533 for 9).
//
// WAVELET: integer reversible 5/3 (LeGall, "rbio2.2" integer-to-integer in Philips' MATLAB), identical in its
//   integer rounding to the JPEG 2000 reversible 5/3 (see isyntax_dwt.c). Philips' encoder does not use
//   symmetric extension at image edges: it zero-pads and keeps only "valid" coefficients, which is why
//   neighbouring tiles must overlap by a few coefficients.
//
// QUANTIZATION AND DEADZONE: ENCODER-ONLY (the decoder never needs them, see the stubs in isyntax.c)
//   The XML carries DP_WAVELET_QUANTIZER (0x101B) and DP_WAVELET_DEADZONE (0x101C) per colour group and per
//   level. They describe what the scanner's encoder did to the wavelet coefficients BEFORE Hulsken coding:
//       for each HL, LH, HH coefficient c of a quantized level (LL is never touched):
//           if abs(c) <= deadzone:  c = 0                                     (deadzone first)
//           c = fix(c / 2^quantizer) * 2^quantizer                            (then truncate toward zero)
//   Nothing is divided out: values keep their original scale, they merely have their lowest `quantizer`
//   magnitude bits cleared. That is why decoding needs no dequantization. In the bitplane domain (see
//   isyntax_hulsken_decompress) the effect is: the lowest `quantizer` magnitude planes are always empty and are
//   simply not stored (their bit in the per-coefficient bitmask is 0), and no nonzero value is smaller than the
//   deadzone (after truncation the smallest survivor is exactly the deadzone value, e.g. 17 -> 16).
//   Parameter mapping: PerColor entry 1 = Y, entry 2 = Co and Cg; PerLevel entries are level 0 first, where
//   level 0 = the first (finest) DWT level = scale 0 in this code; levels beyond the table are unquantized.
//   Philips "Quality=2" ("Q2") table, as found in every sample header:
//                        scale 0 (Quantizer, Deadzone)   scale 1    scale 2 and up
//         Y                       (1, 2)                  (0, 0)     unquantized
//         Co, Cg                  (4, 16)                 (2, 8)     unquantized
//   [real file] Checked on all 17626 non-empty codeblocks: in compressor-v1 mask numbering (bit 0 = sign,
//   bit k = magnitude bit k-1) no block ever sets mask bit 1 for Y scale 0, bits 1-4 for chroma scale 0, or
//   bits 1-2 for chroma scale 1; the smallest nonzero |c| in sampled blocks is 2, 16 and 8 respectively. The spec's
//   wording for the deadzone ("greater or equal") is wrong: the encoder zeroes values less than or equal to it.
//
// COLOUR: lossy iSyntax stores YCoCg-R (Y, Co, Cg). See ycocg_to_rgb() in isyntax.c for the two rounding
//   conventions that exist (they change 64.5% of the pixels of the real file by up to 2 levels).
// =====================================================================================================
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "libisyntax.h"
#include "block_allocator.h"
#include "work_queue.h"
#include "mathutils.h"
#include "platform.h"

#include "yxml.h"

// Coefficients are held in 16 bits. The stored Hulsken data has 16 bit planes per coefficient (sign + 15
// magnitude bits), but in practice values are small: the largest magnitude seen in any intermediate array of an
// independent full reconstruction of the real file was 282 (9 bits), so neither the 16-bit storage nor the
// 16-bit SIMD lifting (which wraps instead of widening) overflows for ordinary 8-bit slides.
#define DWT_COEFF_BITS 16
#if (DWT_COEFF_BITS==16)
typedef i16 icoeff_t;
#else
typedef i32 icoeff_t;
#endif

// Overlap needed for one tile's inverse transform: each of the four subband quadrants (LL, HL, LH, HH) of a
// tile is extended by 4 coefficients on each side with data from the neighbouring tiles of the same scale
// (PAD_L on the left/top, PAD_R on the right/bottom). The 5/3 synthesis needs at most 2 neighbouring
// coefficients per level, so 4 leaves a safety margin. ISYNTAX_IDWT_FIRST_VALID_PIXEL is the index, in the
// (2 * (128 + 8))-wide inverse-transform buffer, of the first output pixel that is not contaminated by the
// buffer edges. Numerically it equals 2 * PAD_L - 1 = 7; the margin itself comes from the reach of the 5/3 synthesis filters.
#define ISYNTAX_IDWT_PAD_L 4
#define ISYNTAX_IDWT_PAD_R 4
#define ISYNTAX_IDWT_FIRST_VALID_PIXEL 7

#define ISYNTAX_ADJ_TILE_TOP_LEFT 0x100
#define ISYNTAX_ADJ_TILE_TOP_CENTER 0x80
#define ISYNTAX_ADJ_TILE_TOP_RIGHT 0x40
#define ISYNTAX_ADJ_TILE_CENTER_LEFT 0x20
#define ISYNTAX_ADJ_TILE_CENTER 0x10
#define ISYNTAX_ADJ_TILE_CENTER_RIGHT 8
#define ISYNTAX_ADJ_TILE_BOTTOM_LEFT 4
#define ISYNTAX_ADJ_TILE_BOTTOM_CENTER 2
#define ISYNTAX_ADJ_TILE_BOTTOM_RIGHT 1


enum isyntax_image_type_enum {
	ISYNTAX_IMAGE_TYPE_NONE = 0,
	ISYNTAX_IMAGE_TYPE_MACROIMAGE = 1,
	ISYNTAX_IMAGE_TYPE_LABELIMAGE = 2,
	ISYNTAX_IMAGE_TYPE_WSI = 3,
};

enum isyntax_node_type_enum {
	ISYNTAX_NODE_NONE = 0,
	ISYNTAX_NODE_LEAF = 1, // ex. <Attribute Name="DICOM_MANUFACTURER" Group="0x0008" Element="0x0070"PMSVR="IString">PHILIPS</ Attribute>
	ISYNTAX_NODE_BRANCH = 2, // ex. <DataObject ObjectType="DPScannedImage"> (leaf nodes) </DataObject>
	ISYNTAX_NODE_ARRAY = 3, // <Array> (contains one or more similar type of leaf/branch nodes)
};

// NOTE: Most of these have DICOM group 0x301D. Currently there seem to be no element ID collisions.
// The wavelet quantizer elements (0x1019 .. 0x101C) are parsed as empty stubs in isyntax.c because the
// decoder does not need them: quantization and deadzone are applied only by the encoder (see the notes at the
// top of this file). Nesting in the XML (inside the WSI DPScannedImage, after DP_COLOR_MANAGEMENT):
//   DP_WAVELET_QUANTIZER_SETTINGS_PER_COLOR (0x1019)  -> array of 2 DPWaveletQuantizerSettingsPerColor
//      [0] = luma (Y), [1] = the other colour components (Co and Cg)
//     DP_WAVELET_QUANTIZER_SETTINGS_PER_LEVEL (0x101A) -> array, level 0 (= scale 0, finest) first
//       DP_WAVELET_QUANTIZER (0x101B): coefficients are truncated to a multiple of 2^value
//       DP_WAVELET_DEADZONE  (0x101C): coefficients with abs(c) <= value are set to 0
enum isyntax_group_data_object_dicom_element_enum {
	// Group 0x301D
	PIM_DP_SCANNED_IMAGES                   = 0x1003, // DPScannedImage
	DP_IMAGE_POST_PROCESSING                = 0x1014, // DPImagePostProcessing
	DP_WAVELET_QUANTIZER_SETTINGS_PER_COLOR = 0x1019, // DPWaveletQuantizerSeetingsPerColor
	DP_WAVELET_QUANTIZER_SETTINGS_PER_LEVEL = 0x101a, // DPWaveletQuantizerSeetingsPerLevel
	UFS_IMAGE_GENERAL_HEADERS               = 0x2000, // UFSImageGeneralHeader
	UFS_IMAGE_DIMENSIONS                    = 0x2003, // UFSImageDimension
	UFS_IMAGE_BLOCK_HEADER_TEMPLATES        = 0x2009, // UFSImageBlockHeaderTemplate
	UFS_IMAGE_DIMENSION_RANGES              = 0x200a, // UFSImageDimensionRange
	DP_COLOR_MANAGEMENT                     = 0x200b, // DPColorManagement
	UFS_IMAGE_BLOCK_HEADERS                 = 0x200d, // UFSImageBlockHeader              // new in iSyntax v2
	UFS_IMAGE_CLUSTER_HEADER_TEMPLATES      = 0x2016, // UFSImageClusterHeaderTemplate    // new in iSyntax v2
	UFS_IMAGE_VALID_DATA_ENVELOPES          = 0x2023, // UFSImageValidDataEnvelope        // new in iSyntax v2
	UFS_IMAGE_OPP_EXTREME_VERTICES          = 0x2024, // UFSImageOppExtremeVertex         // new in iSyntax v2
	// Group 8B01
	PIIM_PIXEL_DATA_REPRESENTATION_SEQUENCE = 0x1001, // PixelDataRepresentation
};

enum isyntax_data_object_flag_enum {
	ISYNTAX_OBJECT_DPUfsImport = 1,
	ISYNTAX_OBJECT_DPScannedImage = 2,
	ISYNTAX_OBJECT_UFSImageGeneralHeader = 4,
	ISYNTAX_OBJECT_UFSImageBlockHeaderTemplate = 8,
	ISYNTAX_OBJECT_UFSImageDimension = 0x10,
	ISYNTAX_OBJECT_UFSImageDimensionRange = 0x20,
	ISYNTAX_OBJECT_DPColorManagement = 0x40,
	ISYNTAX_OBJECT_DPImagePostProcessing = 0x80,
	ISYNTAX_OBJECT_DPWaveletQuantizerSeetingsPerColor = 0x100,
	ISYNTAX_OBJECT_DPWaveletQuantizerSeetingsPerLevel = 0x200,
	ISYNTAX_OBJECT_PixelDataRepresentation = 0x400,
	ISYNTAX_OBJECT_UFSImageBlockHeader = 0x800,             // new in iSyntax v2
	ISYNTAX_OBJECT_UFSImageClusterHeaderTemplate = 0x1000,  // new in iSyntax v2
	ISYNTAX_OBJECT_UFSImageValidDataEnvelope = 0x2000,      // new in iSyntax v2
	ISYNTAX_OBJECT_UFSImageOppExtremeVertex = 0x4000,       // new in iSyntax v2
};


#pragma pack(push, 1)
typedef struct isyntax_dicom_tag_header_t {
	u16 group;
	u16 element;
	u32 size;
} isyntax_dicom_tag_header_t;

typedef struct isyntax_partial_block_header_t {
	isyntax_dicom_tag_header_t sequence_element_header;
	isyntax_dicom_tag_header_t block_coordinates_header;
	u32 x_coordinate;
	u32 y_coordinate;
	u32 color_component;
	u32 scale;
	u32 coefficient;
	/* [MISSING] dicom_tag_header_t block_data_offset_header; */
	/* [MISSING] u64 block_data_offset; */
	/* [MISSING] dicom_tag_header_t block_size_header; */
	/* [MISSING] u64 block_size; */
	isyntax_dicom_tag_header_t block_header_template_id_header;
	u32 block_header_template_id;
} isyntax_partial_block_header_t;

typedef struct isyntax_full_block_header_t {
	isyntax_dicom_tag_header_t sequence_element_header;
	isyntax_dicom_tag_header_t block_coordinates_header;
	u32 x_coordinate;
	u32 y_coordinate;
	u32 color_component;
	u32 scale;
	u32 coefficient;
	isyntax_dicom_tag_header_t block_data_offset_header;
	u64 block_data_offset;
	isyntax_dicom_tag_header_t block_size_header;
	u64 block_size;
	isyntax_dicom_tag_header_t block_header_template_id_header;
	u32 block_header_template_id;
} isyntax_full_block_header_t;

typedef struct isyntax_seektable_codeblock_header_t {
	isyntax_dicom_tag_header_t start_header;
	isyntax_dicom_tag_header_t block_data_offset_header;
	u64 block_data_offset;
	isyntax_dicom_tag_header_t block_size_header;
	u64 block_size;
} isyntax_seektable_codeblock_header_t;
#pragma pack(pop)

typedef struct isyntax_image_dimension_range_t {
	i32 start;
	i32 step;
	i32 end;
	i32 numsteps;
} isyntax_image_dimension_range_t;

// One template per (colour, scale, band group). Observed in v2 headers: 3 * (L + 1) templates (21 for 6
// levels, 27 for 8, 30 for 9): H group for every scale plus one LL template per colour at the top scale.
typedef struct isyntax_block_header_template_t {
	u32 block_width;     // e.g. 128
	u32 block_height;    // e.g. 128
	u8 color_component;  // 0=Y 1=Co 2=Cg
	u8 scale;            // range 0-8
	u8 waveletcoeff;     // 1 = the HL+LH+HH group (a codeblock with coefficient == 1), 3 = LL (coefficient == 0); set from the dimension range start (0 -> 1, otherwise 3). NOTE: this comment used to say the opposite
} isyntax_block_header_template_t;

typedef struct isyntax_cluster_block_header_t {
	u32 x_coordinate;
	u32 y_coordinate;
	u32 color_component;
	u32 scale;
	u32 coefficient;
} isyntax_cluster_block_header_t;

typedef struct isyntax_cluster_relative_coords_t {
	u32 raw_coords[5];
	u32 block_header_template_id;
	u32 x;
	u32 y;
	u32 color_component;
	u32 scale;
	u32 waveletcoeff;
} isyntax_cluster_relative_coords_t;

// Observed cluster sizes (blocks per cluster, all colours): 63 = 3 colours * (16 + 4 + 1) for three scales
// without LL, 66 for a cluster that also holds the top-level LL (3 * 22), and 18 for a top cluster of only two
// scales plus LL (3 * (4 + 1 + 1)). 70 is therefore a comfortable upper bound for the files seen so far.
#define MAX_CODEBLOCKS_PER_CLUSTER 70 // NOTE: what is the actual maximum possible?

typedef struct isyntax_cluster_header_template_t {
	u32 base_x;
	u32 base_y;
	u8 base_scale;
	u8 base_waveletcoeff;
	u8 base_color_component;
	isyntax_cluster_relative_coords_t relative_coords_for_codeblock_in_cluster[MAX_CODEBLOCKS_PER_CLUSTER];
	i32 codeblock_in_cluster_count;
	i32 dimension_order[5];
	u8 dimension_count;
} isyntax_cluster_header_template_t;

typedef struct isyntax_valid_data_envelope_t {
	v2i vertices[64];
	i32 vertex_count;
} isyntax_valid_data_envelope_t;

// One compressed block of coefficients: one colour x one scale x one tile position.
//   coefficient == 1: the three high-pass subbands HL, LH, HH of that tile merged into one block (3 * 128 * 128
//                     values, decoded in that order); coefficient == 0: the LL subband (1 * 128 * 128 values),
//                     which is stored only at the top scale. (The v2 header says "waveletcoeff 3 = LL" and
//                     the code above maps it to coefficient 0.)
//   block_data_offset / block_size: payload position in the file (the 8-byte DICOM item header precedes it).
//                     A size of 8 bytes or less is an "empty" block (all coefficients zero). In real files these
//                     are tiles outside the tissue AND tiles zeroed by the encoder's deadzone: on the real file 81.7%
//                     of Co and 92.6% of Cg blocks at scale 0 are empty.
//   x_coordinate / y_coordinate: absolute position in level-0 pixels of the block's first valid coefficient,
//                     i.e. image origin + (3 << scale) - 2 for H (plus (1 << scale) for LL) + k * (256 << scale).
typedef struct isyntax_codeblock_t {
	u32 x_coordinate;
	u32 y_coordinate;
	u32 color_component;
	u32 scale;
	u32 coefficient;
	u64 block_data_offset;
	u64 block_size;
	u32 block_header_template_id;
	i32 x_adjusted;
	i32 y_adjusted;
	i32 block_x;
	i32 block_y;
	u64 block_id;
} isyntax_codeblock_t;

// A chunk (called a cluster in the v2 header) is a group of codeblocks stored contiguously in the file:
// three consecutive scales for all three colours. Order inside a chunk [real file]: colour-major (Y, Co, Cg);
// within a colour the coarsest scale's block first, then its 2 x 2 children (x fastest), then the 4 x 4
// grandchildren (x fastest); the top-level LL block, when present, comes last for that colour. That is
// 1 + 4 + 16 = 21 blocks per colour (63 per chunk), or 22 per colour (66) when LL is included.
// The real file has 525 chunks of 63 blocks (scales 0..2) followed by 12 chunks of 66 (scales 3..5 + LL).
typedef struct isyntax_data_chunk_t {
	i64 offset;
	u32 size;
	i32 top_codeblock_index;
	i32 codeblock_count_per_color;
	i32 scale;
	i32 level_count;
	u8* data;
} isyntax_data_chunk_t;

typedef struct isyntax_tile_channel_t {
	icoeff_t* coeff_h;
	icoeff_t* coeff_ll;
	u32 neighbors_loaded;
} isyntax_tile_channel_t;

typedef struct isyntax_tile_t {
	u32 codeblock_index;
	u32 codeblock_chunk_index;
	u32 data_chunk_index;
	isyntax_tile_channel_t color_channels[3];
	u32 ll_invalid_edges;
	bool exists;
	bool has_ll;
	bool has_h;
	bool is_submitted_for_h_coeff_decompression;
	bool is_submitted_for_loading;
	bool is_loaded;

    // Cache management.
    // TODO(avirodov): need to rethink this, maybe an external struct that points to isyntax_tile_t. The benefit
    //   is that the cache is usually smaller than the number of tiles. The con is that I'll need to manage list memory
    //   (probably another allocator for small objects - list nodes).
    bool cache_marked;
    struct isyntax_tile_t* cache_next;
    struct isyntax_tile_t* cache_prev;

    // Note(avirodov): this is needed for isyntax_reader. It is very convenient to be able to compute neighbors
    // from the tile itself, although at the cost of additional memory (3 ints) per tile.
    // TODO(avirodov): reconsider this as part of moving out cache_* fields, if applicable.
    // TODO(avirodov): tile_x and tile_y can be computed by O(1) pointer arithmetic given scale. Scale can be
    //  computed as well, but in O(L) where L is number of levels, and that will be computed often.
    int tile_scale;
    int tile_x;
    int tile_y;
} isyntax_tile_t;

typedef struct isyntax_level_t {
	i32 scale;
	i32 width_in_tiles;
	i32 height_in_tiles;
	i32 width;
	i32 height;
	float downsample_factor;
	float um_per_pixel_x;
	float um_per_pixel_y;
	float x_tile_side_in_um;
	float y_tile_side_in_um;
	u64 tile_count;
	i32 origin_offset_in_pixels;
	v2f origin_offset;
	isyntax_tile_t* tiles;
	bool is_fully_loaded;
} isyntax_level_t;

// width_including_padding / height_including_padding are the padded image size from the XML dimension ranges;
// width / height remove the 2 * level0_padding border, where level0_padding = (3 << level_count) - 3.
// codeblocks / data_chunks come from the v1 block header table + seektable, or from the v2 cluster table.
// compressor_version: 1 for UFS_IMAGE_BLOCK_COMPRESSION_METHOD 16, 2 for method 19. All six sample headers and
// the real file use method 16. The container data model (v1/v2) and the compressor version are independent.
typedef struct isyntax_image_t {
	u32 image_type;
    i64 base64_encoded_jpg_file_offset;
    size_t base64_encoded_jpg_len;
	i32 width_including_padding;
	i32 height_including_padding;
	i32 width;
	i32 height;
	i32 level0_padding;
	i32 offset_x;
	i32 offset_y;
	i32 level_count;
	i32 max_scale;
	isyntax_level_t levels[16];
	i32 compressor_version;
	bool compression_is_lossy;
	i32 lossy_image_compression_ratio;
	i32 number_of_blocks;
	i32 codeblock_count;
	isyntax_codeblock_t* codeblocks;
	i32 data_chunk_count;
	isyntax_data_chunk_t* data_chunks;
	bool header_codeblocks_are_partial;
	bool first_load_complete;
	bool first_load_in_progress;
	i64 base64_encoded_icc_profile_file_offset;
	size_t base64_encoded_icc_profile_len;
} isyntax_image_t;

typedef struct isyntax_parser_node_t {
	u32 node_type; // leaf, branch, or array
	bool has_children;
	bool has_base64_content;
	u16 group;
	u16 element;
} isyntax_parser_node_t;

#define ISYNTAX_MAX_NODE_DEPTH 16

typedef struct isyntax_xml_parser_t {
	yxml_t* x;
	isyntax_image_t* current_image;
	i32 running_image_index;
	u32 current_image_type;
	char* attrbuf;
	char* attrbuf_end;
	char* attrcur;
	size_t attrlen;
	size_t attrbuf_capacity;
	char* contentbuf;
	char* contentcur;
	size_t contentlen;
	size_t contentbuf_capacity;
    i64 content_file_offset;
	char current_dicom_attribute_name[256];
	u32 current_dicom_group_tag;
	u32 current_dicom_element_tag;
	i32 attribute_index;
	u32 current_node_type;
	bool current_node_has_children;
	isyntax_parser_node_t node_stack[ISYNTAX_MAX_NODE_DEPTH];
	i32 node_stack_index;
	isyntax_parser_node_t data_object_stack[ISYNTAX_MAX_NODE_DEPTH];
	i32 data_object_stack_index;
	u32 data_object_flags;
	i32 block_header_template_index;
	i32 cluster_header_template_index;
	i32 block_header_index_for_cluster;
	i32 dimension_index;
	i32 valid_data_envelope_index;
	bool initialized;
} isyntax_xml_parser_t;

// data_model_major_version: < 100 (usually 5) = v1 layout, >= 100 (e.g. 100.4) = v2 layout (see the notes at the top).
// dicom_derivation_description carries e.g. "PHILIPS UFS V1.8.6614 | Quality=2 | DWT=1 | Compressor=16":
// Quality=2 selects Philips' "Q2" quantizer/deadzone table, DWT=1 the 5/3 wavelet, Compressor=16 the Hulsken v1 coder.
typedef struct isyntax_t {
	enum libisyntax_open_flags_t open_flags;
	i64 filesize;
	file_handle_t file_handle;
	isyntax_image_t images[16];
	i32 image_count;
	isyntax_block_header_template_t block_header_templates[64];
	i32 block_header_template_count;
	isyntax_cluster_header_template_t cluster_header_templates[8];
	i32 cluster_header_template_count;
	isyntax_valid_data_envelope_t valid_data_envelopes[16];
	i32 valid_data_envelope_count;
	i32 macro_image_index;
	i32 label_image_index;
	i32 wsi_image_index;
	isyntax_xml_parser_t parser;
	float mpp_x;
	float mpp_y;
	bool is_mpp_known;
	i32 block_width;
	i32 block_height;
	i32 tile_width;
	i32 tile_height;
	icoeff_t* black_dummy_coeff;
	icoeff_t* white_dummy_coeff;
	block_allocator_t* ll_coeff_block_allocator;
	block_allocator_t* h_coeff_block_allocator;
    bool is_block_allocator_owned;
	float loading_time;
	float total_rgb_transform_time;
	i32 data_model_major_version; // <100 (usually 5) for iSyntax format v1, >= 100 for iSyntax format v2
	i32 data_model_minor_version;
	char barcode[64];
	bool is_barcode_read;
	isyntax_cache_t* cache;
	thread_pool_t* work_submission_pool;
	volatile i32 refcount;
	char dicom_acquisition_datetime[33]; // e.g. "20210609111602.000000"
	char dicom_manufacturer[65]; // e.g. "PHILIPS"
	char dicom_manufacturers_model_name[65]; // e.g. "UFS Scanner"
	char dicom_device_serial_number[65]; // e.g. "FMT0586"
	char (*dicom_software_versions)[65]; // e.g. ""1.8.6824" "20180906_R51""
	i32 dicom_software_versions_count;
	char dicom_derivation_description[1025]; // e.g. PHILIPS UFS V1.8.6824 | Quality=2 | DWT=1 | Compressor=16
	char (*dicom_date_of_last_calibration)[9];
	i32 dicom_date_of_last_calibration_count;
	char (*dicom_time_of_last_calibration)[7];
	i32 dicom_time_of_last_calibration_count;
	bool dicom_lossy_image_compression;
	double dicom_lossy_image_compression_ratio;
	char dicom_lossy_image_compression_method[65]; // e.g. "PHILIPS_DP_1_0"
	char image_dimension_unit[65];
} isyntax_t;

// function prototypes
// Decodes one codeblock (compressed payload) into coeff_count * block_width * block_height i16 values in
// two's complement, where coeff_count is 3 for coefficient == 1 (order HL, LH, HH) and 1 for LL.
// See the long description at the definition in isyntax.c.
bool isyntax_hulsken_decompress(u8 *compressed, size_t compressed_size, i32 block_width, i32 block_height, i32 coefficient, i32 compressor_version, i16* out_buffer);
void isyntax_set_thread_pool(isyntax_t* isyntax, thread_pool_t* thread_pool);
bool isyntax_open(isyntax_t* isyntax, const char* filename, enum libisyntax_open_flags_t flags);
void isyntax_destroy(isyntax_t* isyntax);
void isyntax_idwt(icoeff_t* idwt, i32 quadrant_width, i32 quadrant_height, bool output_steps_as_png, const char* png_name);
void isyntax_load_tile(isyntax_t* isyntax, isyntax_image_t* wsi, i32 scale, i32 tile_x, i32 tile_y, block_allocator_t* ll_coeff_block_allocator,
                       u32* out_buffer_or_null, enum isyntax_pixel_format_t pixel_format);
u32 isyntax_get_adjacent_tiles_mask(isyntax_level_t* level, i32 tile_x, i32 tile_y);
u32 isyntax_get_adjacent_tiles_mask_only_existing(isyntax_level_t* level, i32 tile_x, i32 tile_y);
u32 isyntax_idwt_tile_for_color_channel(isyntax_t* isyntax, isyntax_image_t* wsi, i32 scale, i32 tile_x, i32 tile_y, i32 color, icoeff_t* dest_buffer);
void isyntax_decompress_codeblock_in_chunk(isyntax_codeblock_t* codeblock, i32 block_width, i32 block_height, u8* chunk, u64 chunk_base_offset, i32 compressor_version, i16* out_buffer);
i32 isyntax_get_chunk_codeblocks_per_color_for_level(i32 level, bool has_ll);
u8* isyntax_get_associated_image_pixels(isyntax_t* isyntax, isyntax_image_t* image, enum isyntax_pixel_format_t pixel_format);
u8* isyntax_get_associated_image_jpeg(isyntax_t* isyntax, isyntax_image_t* image, u32* jpeg_size);
u8* isyntax_get_icc_profile(isyntax_t* isyntax, isyntax_image_t* image, u32* icc_profile_size);


#ifdef __cplusplus
}
#endif
