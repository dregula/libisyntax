// Writes one resolution layer of the iSyntax -> DICOM-WSI transcoding as a DICOM file (VL Whole Slide Microscopy Image Storage, JPEG 2000 Lossless Only,
// YBR_RCT, TILED_FULL). Metadata through DCMTK dcmdata (3.6.9); Pixel Data streamed by hand (encapsulated, empty Basic Offset Table, one fragment per frame,
// Extended Offset Table patched in after the frames are known). Most Type 1/2 tags carry dummy values; the geometry tags carry real values.
// Written for Halo_Strategy / Pipeline_Plan_v2 (new file, decoder untouched).
#pragma once
#include <string>
#include <vector>
#include <cstdint>
struct LayerMeta {
    int layer = 0, rows = 256, cols = 256;                    // frame size
    int frames_x = 0, frames_y = 0;                           // frame grid; NumberOfFrames = frames_x * frames_y (raster order)
    long matrix_cols = 0, matrix_rows = 0;                    // Total Pixel Matrix
    double spacing_mm = 0.00025;                              // pixel spacing of this layer (0.25 um x 2^layer)
    double origin_x_mm = 0, origin_y_mm = 0;                  // slide-coordinate position of the top-left corner of the matrix (see --orientation: mapping unverified)
    double orientation[6] = {1, 0, 0, 0, 1, 0};              // Image Orientation (Slide): row direction cosines, then column direction cosines
    std::string study_uid, series_uid, for_uid, sop_uid, specimen_uid, dimorg_uid;
    std::string source_name, manufacturer, model, serial, software, acquisition_datetime, lossy_method;
    double source_ratio = 0, our_ratio = 0;                   // Lossy Image Compression Ratio: source (iSyntax, 0 = unknown) and this file's own
    std::string icc_path;                                     // sRGB profile for the optical path
    bool resampled = false;                                   // layers 1.. are RESAMPLED
};
// spool: file holding, for every frame in order, the item (FFFE,E000) header + even-padded codestream; offsets/lengths per frame as defined for the Extended Offset Table
bool write_layer_file(const LayerMeta& m, const std::string& out_path, const std::string& spool_path, const std::vector<uint64_t>& offsets,
                      const std::vector<uint64_t>& lengths, std::string& err);
std::string make_uid(const std::string& seed, const std::string& purpose);     // 2.25.<128-bit decimal>, deterministic from the seed text
