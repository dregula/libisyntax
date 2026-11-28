#include "dicom_layer_writer.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcvrov.h"

std::string make_uid(const std::string& seed, const std::string& purpose) {              // two splitmix64 words from an FNV-1a hash of seed + purpose
    uint64_t h = 1469598103934665603ull; for (unsigned char c : seed + "|" + purpose) { h ^= c; h *= 1099511628211ull; }
    auto sm = [](uint64_t& x) { uint64_t z = (x += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); };
    unsigned __int128 v = ((unsigned __int128)sm(h) << 64) | sm(h); v &= (((unsigned __int128)1) << 124) - 1;   // < 2^124: keeps the decimal at 38 digits at most
    std::string d; if (v == 0) d = "0"; while (v) { d.insert(d.begin(), (char)('0' + (int)(v % 10))); v /= 10; }
    return "2.25." + d;
}
static std::string ds(double v) { char b[64]; snprintf(b, sizeof b, "%.10g", v); std::string s = b; if (s.size() > 16) { snprintf(b, sizeof b, "%.8g", v); s = b; } return s; }
static void code_item(DcmItem* seq_item, const char* value, const char* scheme, const char* meaning) {
    seq_item->putAndInsertString(DCM_CodeValue, value); seq_item->putAndInsertString(DCM_CodingSchemeDesignator, scheme); seq_item->putAndInsertString(DCM_CodeMeaning, meaning);
}
#define CHK(x) do { OFCondition _c = (x); if (_c.bad()) { err = std::string(#x) + ": " + _c.text(); return false; } } while (0)

bool write_layer_file(const LayerMeta& m, const std::string& out_path, const std::string& spool_path, const std::vector<uint64_t>& offsets, const std::vector<uint64_t>& lengths, std::string& err) {
    DcmFileFormat ff; DcmDataset* d = ff.getDataset(); size_t nframes = offsets.size();
    if (nframes != (size_t)m.frames_x * m.frames_y || lengths.size() != nframes) { err = "frame count mismatch"; return false; }
    std::string date = "20000101", time = "000000";
    if (m.acquisition_datetime.size() >= 14) { date = m.acquisition_datetime.substr(0, 8); time = m.acquisition_datetime.substr(8, 6); }
    // ---- SOP common, patient, study, series (dummy values)
    CHK(d->putAndInsertString(DCM_SOPClassUID, UID_VLWholeSlideMicroscopyImageStorage)); CHK(d->putAndInsertString(DCM_SOPInstanceUID, m.sop_uid.c_str()));
    CHK(d->putAndInsertString(DCM_PatientName, "Anonymous^iSyntaxTranscode")); CHK(d->putAndInsertString(DCM_PatientID, m.source_name.c_str())); CHK(d->putAndInsertString(DCM_PatientBirthDate, "")); CHK(d->putAndInsertString(DCM_PatientSex, ""));
    CHK(d->putAndInsertString(DCM_StudyInstanceUID, m.study_uid.c_str())); CHK(d->putAndInsertString(DCM_StudyDate, date.c_str())); CHK(d->putAndInsertString(DCM_StudyTime, time.c_str()));
    CHK(d->putAndInsertString(DCM_ReferringPhysicianName, "")); CHK(d->putAndInsertString(DCM_StudyID, "1")); CHK(d->putAndInsertString(DCM_AccessionNumber, ""));
    CHK(d->putAndInsertString(DCM_Modality, "SM")); CHK(d->putAndInsertString(DCM_SeriesInstanceUID, m.series_uid.c_str())); CHK(d->putAndInsertString(DCM_SeriesNumber, "1"));
    CHK(d->putAndInsertString(DCM_FrameOfReferenceUID, m.for_uid.c_str())); CHK(d->putAndInsertString(DCM_PositionReferenceIndicator, "SLIDE_CORNER"));
    CHK(d->putAndInsertString(DCM_Manufacturer, m.manufacturer.c_str())); CHK(d->putAndInsertString(DCM_ManufacturerModelName, m.model.c_str()));
    CHK(d->putAndInsertString(DCM_DeviceSerialNumber, m.serial.c_str())); CHK(d->putAndInsertString(DCM_SoftwareVersions, m.software.c_str()));
    CHK(d->putAndInsertString(DCM_InstanceNumber, std::to_string(m.layer + 1).c_str())); CHK(d->putAndInsertString(DCM_ContentDate, date.c_str())); CHK(d->putAndInsertString(DCM_ContentTime, time.c_str()));
    CHK(d->putAndInsertString(DCM_AcquisitionDateTime, (date + time).c_str()));
    // ---- image module
    std::string itype = std::string("DERIVED\\PRIMARY\\VOLUME\\") + (m.resampled ? "RESAMPLED" : "NONE");
    CHK(d->putAndInsertString(DCM_ImageType, itype.c_str())); CHK(d->putAndInsertString(DCM_VolumetricProperties, "VOLUME"));
    CHK(d->putAndInsertString(DCM_SamplesPerPixel, "3")); CHK(d->putAndInsertString(DCM_PhotometricInterpretation, "YBR_RCT")); CHK(d->putAndInsertUint16(DCM_PlanarConfiguration, 0));
    CHK(d->putAndInsertUint16(DCM_Rows, (Uint16)m.rows)); CHK(d->putAndInsertUint16(DCM_Columns, (Uint16)m.cols)); CHK(d->putAndInsertUint16(DCM_BitsAllocated, 8)); CHK(d->putAndInsertUint16(DCM_BitsStored, 8));
    CHK(d->putAndInsertUint16(DCM_HighBit, 7)); CHK(d->putAndInsertUint16(DCM_PixelRepresentation, 0)); CHK(d->putAndInsertString(DCM_NumberOfFrames, std::to_string(nframes).c_str()));
    CHK(d->putAndInsertString(DCM_BurnedInAnnotation, "NO")); CHK(d->putAndInsertString(DCM_SpecimenLabelInImage, "NO")); CHK(d->putAndInsertString(DCM_FocusMethod, "AUTO")); CHK(d->putAndInsertString(DCM_ExtendedDepthOfField, "NO"));
    CHK(d->putAndInsertUint32(DCM_NumberOfOpticalPaths, 1)); CHK(d->putAndInsertUint32(DCM_TotalPixelMatrixFocalPlanes, 1));
    CHK(d->putAndInsertString(DCM_LossyImageCompression, "01"));
    { std::string meth = m.lossy_method.empty() ? "ISO_15444_1" : m.lossy_method + "\\ISO_15444_1", ratio = m.source_ratio > 0 ? ds(m.source_ratio) + "\\" + ds(m.our_ratio) : ds(m.our_ratio);
      CHK(d->putAndInsertString(DCM_LossyImageCompressionMethod, meth.c_str())); CHK(d->putAndInsertString(DCM_LossyImageCompressionRatio, ratio.c_str())); }
    // ---- whole-slide geometry
    CHK(d->putAndInsertUint32(DCM_TotalPixelMatrixColumns, (Uint32)m.matrix_cols)); CHK(d->putAndInsertUint32(DCM_TotalPixelMatrixRows, (Uint32)m.matrix_rows));
    CHK(d->putAndInsertFloat32(DCM_ImagedVolumeWidth, (Float32)(m.matrix_cols * m.spacing_mm))); CHK(d->putAndInsertFloat32(DCM_ImagedVolumeHeight, (Float32)(m.matrix_rows * m.spacing_mm))); CHK(d->putAndInsertFloat32(DCM_ImagedVolumeDepth, 0.003f));
    { std::string o; for (int i = 0; i < 6; i++) o += (i ? "\\" : "") + ds(m.orientation[i]); CHK(d->putAndInsertString(DCM_ImageOrientationSlide, o.c_str())); }
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_TotalPixelMatrixOriginSequence, it, -2)); CHK(it->putAndInsertString(DCM_XOffsetInSlideCoordinateSystem, ds(m.origin_x_mm).c_str())); CHK(it->putAndInsertString(DCM_YOffsetInSlideCoordinateSystem, ds(m.origin_y_mm).c_str())); }
    // ---- specimen, optical path
    CHK(d->putAndInsertString(DCM_ContainerIdentifier, "SLIDE1")); CHK(d->insert(new DcmSequenceOfItems(DCM_IssuerOfTheContainerIdentifierSequence)));
    CHK(d->insert(new DcmSequenceOfItems(DCM_AcquisitionContextSequence)));
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_ContainerTypeCodeSequence, it, -2)); code_item(it, "433466003", "SCT", "Microscope slide"); }
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_SpecimenDescriptionSequence, it, -2)); CHK(it->putAndInsertString(DCM_SpecimenIdentifier, "SPECIMEN1")); CHK(it->insert(new DcmSequenceOfItems(DCM_IssuerOfTheSpecimenIdentifierSequence))); CHK(it->putAndInsertString(DCM_SpecimenUID, m.specimen_uid.c_str()));
      CHK(it->insert(new DcmSequenceOfItems(DCM_SpecimenPreparationSequence))); }
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_OpticalPathSequence, it, -2)); CHK(it->putAndInsertString(DCM_OpticalPathIdentifier, "1"));
      std::ifstream f(m.icc_path, std::ios::binary); std::vector<char> icc((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>()); if (icc.empty()) { err = "ICC profile not readable: " + m.icc_path; return false; }
      CHK(it->putAndInsertUint8Array(DCM_ICCProfile, (const Uint8*)icc.data(), icc.size()));
      DcmItem* c = nullptr; CHK(it->findOrCreateSequenceItem(DCM_IlluminationColorCodeSequence, c, -2)); code_item(c, "414298005", "SCT", "Full Spectrum");
      CHK(it->findOrCreateSequenceItem(DCM_IlluminationTypeCodeSequence, c, -2)); code_item(c, "111744", "DCM", "Brightfield illumination"); }
    // ---- multi-frame dimension (tiled full) and functional groups
    CHK(d->putAndInsertString(DCM_DimensionOrganizationType, "TILED_FULL"));
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_DimensionOrganizationSequence, it, -2)); CHK(it->putAndInsertString(DCM_DimensionOrganizationUID, m.dimorg_uid.c_str())); }
    { DcmItem* it = nullptr; CHK(d->findOrCreateSequenceItem(DCM_DimensionIndexSequence, it, -2)); CHK(it->putAndInsertTagKey(DCM_DimensionIndexPointer, DCM_ColumnPositionInTotalImagePixelMatrix));
      CHK(it->putAndInsertTagKey(DCM_FunctionalGroupPointer, DCM_PlanePositionSlideSequence)); CHK(it->putAndInsertString(DCM_DimensionOrganizationUID, m.dimorg_uid.c_str()));
      CHK(d->findOrCreateSequenceItem(DCM_DimensionIndexSequence, it, -2)); CHK(it->putAndInsertTagKey(DCM_DimensionIndexPointer, DCM_RowPositionInTotalImagePixelMatrix));
      CHK(it->putAndInsertTagKey(DCM_FunctionalGroupPointer, DCM_PlanePositionSlideSequence)); CHK(it->putAndInsertString(DCM_DimensionOrganizationUID, m.dimorg_uid.c_str())); }
    { DcmItem* sh = nullptr; CHK(d->findOrCreateSequenceItem(DCM_SharedFunctionalGroupsSequence, sh, -2)); DcmItem* x = nullptr;
      CHK(sh->findOrCreateSequenceItem(DCM_PixelMeasuresSequence, x, -2)); std::string ps = ds(m.spacing_mm) + "\\" + ds(m.spacing_mm); CHK(x->putAndInsertString(DCM_PixelSpacing, ps.c_str())); CHK(x->putAndInsertString(DCM_SliceThickness, "0.003"));
      CHK(sh->findOrCreateSequenceItem(DCM_WholeSlideMicroscopyImageFrameTypeSequence, x, -2)); CHK(x->putAndInsertString(DCM_FrameType, itype.c_str()));
      CHK(sh->findOrCreateSequenceItem(DCM_OpticalPathIdentificationSequence, x, -2)); CHK(x->putAndInsertString(DCM_OpticalPathIdentifier, "1")); }
    // ---- Extended Offset Table placeholders (patched after the frames are appended)
    { std::vector<Uint64> zero(nframes, 0); DcmOther64bitVeryLong* a = new DcmOther64bitVeryLong(DCM_ExtendedOffsetTable); CHK(a->putUint64Array(zero.data(), nframes)); CHK(d->insert(a));
      DcmOther64bitVeryLong* b = new DcmOther64bitVeryLong(DCM_ExtendedOffsetTableLengths); CHK(b->putUint64Array(zero.data(), nframes)); CHK(d->insert(b)); }
    // ---- write metadata part
    CHK(ff.saveFile(out_path.c_str(), EXS_JPEG2000LosslessOnly, EET_ExplicitLength, EGL_withoutGL));
    std::vector<char> head; { std::ifstream f(out_path, std::ios::binary); head.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>()); }
    auto find_hdr = [&](uint16_t grp, uint16_t el) -> long { unsigned char pat[8] = {(unsigned char)(grp & 255), (unsigned char)(grp >> 8), (unsigned char)(el & 255), (unsigned char)(el >> 8), 'O', 'V', 0, 0};
        auto it = std::search(head.begin(), head.end(), (char*)pat, (char*)pat + 8); return it == head.end() ? -1 : (long)(it - head.begin()) + 12; };       // offset of the value
    long po = find_hdr(0x7FE0, 0x0001), pl = find_hdr(0x7FE0, 0x0002); if (po < 0 || pl < 0) { err = "extended offset table placeholder not found in the written metadata"; return false; }
    for (size_t i = 0; i < nframes; i++) { memcpy(&head[po + 8 * i], &offsets[i], 8); memcpy(&head[pl + 8 * i], &lengths[i], 8); }                                   // little endian host
    FILE* o = fopen(out_path.c_str(), "wb"); if (!o) { err = "cannot write " + out_path; return false; }
    fwrite(head.data(), 1, head.size(), o);
    const unsigned char pix_hdr[12] = {0xE0, 0x7F, 0x10, 0x00, 'O', 'B', 0, 0, 0xFF, 0xFF, 0xFF, 0xFF}, bot[8] = {0xFE, 0xFF, 0x00, 0xE0, 0, 0, 0, 0}, seq_end[8] = {0xFE, 0xFF, 0xDD, 0xE0, 0, 0, 0, 0};
    fwrite(pix_hdr, 1, 12, o); fwrite(bot, 1, 8, o);
    FILE* sp = fopen(spool_path.c_str(), "rb"); if (!sp) { fclose(o); err = "cannot read spool"; return false; }
    std::vector<char> buf(1 << 20); size_t n; while ((n = fread(buf.data(), 1, buf.size(), sp)) > 0) fwrite(buf.data(), 1, n, o);
    fclose(sp); fwrite(seq_end, 1, 8, o); fclose(o); return true;
}
