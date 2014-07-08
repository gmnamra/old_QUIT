#include "Header.h"
#include "nifti1.h"

using namespace std;

namespace Nifti {

/*
 * Returns the DataType enum for a particular code
 *
 *\sa NIFTI1_TYPE group in nifti1.h
 */
DataType DataTypeForCode(const int code) {
	static const map<int, DataType> c2dt{
		{NIFTI_TYPE_UINT8, DataType::UINT8},
		{NIFTI_TYPE_UINT16, DataType::UINT16},
		{NIFTI_TYPE_UINT32, DataType::UINT32},
		{NIFTI_TYPE_UINT64, DataType::UINT64},
		{NIFTI_TYPE_INT8, DataType::INT8},
		{NIFTI_TYPE_INT16, DataType::INT16},
		{NIFTI_TYPE_INT32, DataType::INT32},
		{NIFTI_TYPE_INT64, DataType::INT64},
		{NIFTI_TYPE_FLOAT32, DataType::FLOAT32},
		{NIFTI_TYPE_FLOAT64, DataType::FLOAT64},
		{NIFTI_TYPE_FLOAT128, DataType::FLOAT128},
		{NIFTI_TYPE_COMPLEX64, DataType::COMPLEX64},
		{NIFTI_TYPE_COMPLEX128, DataType::COMPLEX128},
		{NIFTI_TYPE_COMPLEX256, DataType::COMPLEX256},
		{NIFTI_TYPE_RGB24, DataType::RGB24},
		{NIFTI_TYPE_RGBA32, DataType::RGBA32},
	};
	auto dt = c2dt.find(code);
	if (dt != c2dt.end())
		return dt->second;
	else
		throw(std::invalid_argument("Unsupported data format code: " + to_string(code)));
}


/*	Internal map of datatype properties
 *
 *  The map is declared here because making it a static member of Nifti1 was
 *  causing problems with looking up the datatype in close() when called by
 *  ~Nifti1. It's possible for C++ to destruct static members even when
 *  objects still exist in another translation unit.
 */
const DataTypeInfo &TypeInfo(const DataType dt) {
	static const map<DataType, DataTypeInfo> DTInfo{
		{DataType::UINT8,      {DataType::UINT8,      NIFTI_TYPE_UINT8,       1,  0, "UINT8"} },
		{DataType::UINT16,     {DataType::UINT16,     NIFTI_TYPE_UINT16,      2,  2, "UINT16"} },
		{DataType::UINT32,     {DataType::UINT32,     NIFTI_TYPE_UINT32,      4,  4, "UINT32"} },
		{DataType::UINT64,     {DataType::UINT64,     NIFTI_TYPE_UINT64,      8,  8, "UINT64"} },
		{DataType::INT8,       {DataType::INT8,       NIFTI_TYPE_INT8,        1,  0, "INT8"} },
		{DataType::INT16,      {DataType::INT16,      NIFTI_TYPE_INT16,       2,  2, "INT16"} },
		{DataType::INT32,      {DataType::INT32,      NIFTI_TYPE_INT32,       4,  4, "INT32"} },
		{DataType::INT64,      {DataType::INT64,      NIFTI_TYPE_INT64,       8,  8, "INT64"} },
		{DataType::FLOAT32,    {DataType::FLOAT32,    NIFTI_TYPE_FLOAT32,     4,  4, "FLOAT32"} },
		{DataType::FLOAT64,    {DataType::FLOAT64,    NIFTI_TYPE_FLOAT64,     8,  8, "FLOAT64"} },
		{DataType::FLOAT128,   {DataType::FLOAT128,   NIFTI_TYPE_FLOAT128,   16, 16, "FLOAT128"} },
		{DataType::COMPLEX64,  {DataType::COMPLEX64,  NIFTI_TYPE_COMPLEX64,   8,  4, "COMPLEX64"} },
		{DataType::COMPLEX128, {DataType::COMPLEX128, NIFTI_TYPE_COMPLEX128, 16,  8, "COMPLEX128"} },
		{DataType::COMPLEX256, {DataType::COMPLEX256, NIFTI_TYPE_COMPLEX256, 32, 16, "COMPLEX256"} },
		{DataType::RGB24,      {DataType::RGB24,      NIFTI_TYPE_RGB24,       3,  0, "RGB24"} },
		{DataType::RGBA32,     {DataType::RGBA32,     NIFTI_TYPE_RGBA32,      4,  0, "RGBA32"} }
	};
	auto info = DTInfo.find(dt);
	if (info != DTInfo.end())
		return info->second;
	else
		throw(std::invalid_argument("Missing type information, contact libNifti1 author."));
}

/*
 * Returns the string representation of a NIfTI transform code.
 *
 *\sa NIFTI1_XFORM_CODES group in nifti1.h
 */
const string XFormName(const XForm c) {
	switch (c) {
		case XForm::Unknown: return "Unknown"; break;
		case XForm::ScannerAnatomy: return "Scanner Anatomy"; break;
		case XForm::AlignedAnatomy: return "Aligned Anatomy"; break;
		case XForm::Talairach: return "Talairach"; break;
		case XForm::MNI_152: return "MNI 152"; break;
	}
}

/*
 * Converts a NIfTI transform enum into the corresponding integer code
 *
 *\sa NIFTI1_XFORM_CODES group in nifti1.h
 */
int XFormCode(const XForm c) {
	switch (c) {
		case XForm::Unknown: return NIFTI_XFORM_UNKNOWN; break;
		case XForm::ScannerAnatomy: return NIFTI_XFORM_SCANNER_ANAT; break;
		case XForm::AlignedAnatomy: return NIFTI_XFORM_ALIGNED_ANAT; break;
		case XForm::Talairach: return NIFTI_XFORM_TALAIRACH; break;
		case XForm::MNI_152: return NIFTI_XFORM_MNI_152; break;
	}
}

/*
 * Returns an integer code into the corresponding NIfTI transform enum
 *
 *\sa NIFTI1_XFORM_CODES group in nifti1.h
 */
XForm XFormForCode(const int c) {
	switch (c) {
		case NIFTI_XFORM_UNKNOWN: return XForm::Unknown; break;
		case NIFTI_XFORM_SCANNER_ANAT: return XForm::ScannerAnatomy; break;
		case NIFTI_XFORM_ALIGNED_ANAT: return XForm::AlignedAnatomy; break;
		case NIFTI_XFORM_TALAIRACH: return XForm::Talairach; break;
		case NIFTI_XFORM_MNI_152: return XForm::MNI_152; break;
		default:
			throw(std::invalid_argument("Invalid transform code: " + to_string(c)));
	}
}
} // End namespace Nifti
