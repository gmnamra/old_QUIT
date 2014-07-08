/** \file Header.h
 \brief Declaration for Header class
 - Written by Tobias Wood, IoP KCL
 */
#ifndef NIFTI_HEADER_H
#define NIFTI_HEADER_H

#include <string>
#include <map>

namespace Nifti {
enum class DataType {
	UINT8, UINT16, UINT32, UINT64, INT8, INT16, INT32, INT64,
	FLOAT32, FLOAT64, FLOAT128, COMPLEX64, COMPLEX128, COMPLEX256,
	RGB24, RGBA32
};
DataType DataTypeForCode(const int code);

struct DataTypeInfo {
	DataType type;
	size_t code, size, swapsize;
	std::string name;
}; //!< Contains all the information needed to read/write a Nifti datatype
const DataTypeInfo &TypeInfo(const DataType dt);

enum class XForm {
	Unknown, ScannerAnatomy, AlignedAnatomy, Talairach, MNI_152
};
const std::string XFormName(const XForm t);
int XFormCode(const XForm c);
XForm XFormForCode(const int c);

class Header {

};

} // End namespace Nifti

#endif // NIFTI_HEADER_H
