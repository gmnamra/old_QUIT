/** \file File.h
 \brief Declaration for File class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef NIFTI_IMAGE
#define NIFTI_IMAGE

#include <string>
#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>
#include <list>
#include <map>
#include <limits>
#include <exception>

#include "Eigen/Geometry"

#include "nifti1.h" // NIFTI-1 header specification
#include "nifti_analyze.h" // NIFTI version of the ANALYZE 7.5 header
#include "ZipFile.h"
#include "Extension.h"

using namespace std;
using namespace Eigen;

#pragma mark Start namespace Nifti
namespace Nifti {

#pragma mark Enums, structs and typedefs
struct DataType {
	int code, size, swapsize;
	string name;
};
typedef map<int, DataType> DTMap;
static void printDTypeList();
typedef map<int, string> StringMap;

/*
 *  Used when opening a File to specify read or write. Files
 *  can only be open for either reading or writing at any one time, and
 *  must be closed before re-opening. READ_HEADER is a special mode that
 *  will just read the header and then close the file.
 *
 */
enum class Modes : char {
	Closed = 0,
	Read = 'r',
	ReadHeader = 'h',
	ReadSkipExt = 's',
	Write = 'w',
	WriteSkipExt = 'x'
};

#pragma mark NIfTI File Class
class File {
	private:		
		static const DTMap &DataTypes();
		
		Array<size_t, 7, 1> m_dim;   //!< Number of voxels in each dimension. Note that here we do NOT store the rank in dim[0], so only 7 elements required.
		Array<float, 7, 1> m_voxdim; //!< Size of each voxel. As above, only 7 elements because the rank is not stored.
		Affine3f m_qform, m_sform;   //!< Tranformation matrices from voxel indices to physical co-ords.
		
		string m_basepath;            //!< Path to file without extension.
		bool m_nii, m_gz;
		Modes m_mode;                 //!< Whether the file is closed or open for reading/writing.
		ZipFile m_file;
		DataType m_datatype;          //!< Datatype on disk.
		int m_voxoffset;              //!< Offset to start of voxel data.
		int m_swap;                   //!< True if byte order on disk is different to CPU.
		
		list<Extension> m_extensions;
		
		static int needs_swap(short dim0, int hdrsize); //!< Check if file endianism matches host endianism.
		static float fixFloat(const float f); //!< Converts invalid floats to 0 to ensure a marginally sane header
		
		static void SwapBytes(size_t n, int siz, void *ar);
		static void SwapNiftiHeader(struct nifti_1_header *h);
		static void SwapAnalyzeHeader(nifti_analyze75 *h);
		
		void readHeader();      //!< Attempts to read a header structure from the currently open file.
		void readExtensions();  //!< Attempts to read any extensions
		void writeHeader();     //!< Attempts to write a header structure to the currently open file.
		void writeExtensions(); //!< Attempts to write extensions
		int totalExtensionSize(); //!< Counts the total number of bytes for all extensions.
		char *readBytes(size_t start, size_t length, char *buffer);
		void writeBytes(size_t start, size_t length, char *buffer);
		
		template<typename T> void convertFromBytes(const vector<char> &bytes, const size_t nEl, vector<T> &data);
		template<typename T> void convertFromBytes(const vector<complex<T>> &bytes, const size_t nEl, vector<T> &data);
		template<typename T> vector<char> convertToBytes(const vector<T> &data);
		template<typename T> vector<char> convertToBytes(const vector<complex<T>> &data);
		
	#pragma mark Public Methods
	public:
		~File();
		File();                              //!< Default constructor. Initialises an empty header, size 1 in all dimensions.
		File(const File &other);             //!< Copy constructor. Copies all elements, and if the original is open then also opens new file handles.
		File &operator=(const File &other);  //!< Assignment. Copies all elements except file handles, and marks destination as Closed.
		File(File &&other) noexcept;         //!< Move constructor. Copies all elements, including the file handles, and marks the original as Closed.
		
		File(const int nx, const int ny, const int nz, const int nt,
			 const float dx, const float dy, const float dz, const float dt,
			 const int datatype = NIFTI_TYPE_FLOAT32, const Affine3f &transform = Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		File(const Array<size_t, Dynamic, 1> &dim, const ArrayXf &voxdim,
			 const int datatype = NIFTI_TYPE_FLOAT32, const Affine3f &transform = Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		File(const File &other, const size_t nt, const int datatype = NIFTI_TYPE_FLOAT32);               //!< Copies only basic geometry information from other, then sets the datatype and number of volumes. Does not copy scaling information etc.
		File(const string &filename, const Modes &mode);
		
		void open(const string &filename, const Modes &mode); //!< Attempts to open a NIfTI file. Throws runtime_error or invalid_argument on failure.
		void close();                                         //!< Closes the file
		bool isOpen();                                        //!< Returns true if file is currently open for reading or writing.
		
		const string basePath() const;
		const string imagePath() const;
		const string headerPath() const;
		char *readRawVolume(const int vol);
		char *readRawAllVolumes();
		
		size_t dimensions() const;                              //!< Get the number of dimensions (rank) of this image.
		size_t dim(const size_t d) const;                       //!< Get the size (voxel count) of a dimension. Valid dimensions are 1-7.
		void setDim(const size_t d, const size_t n);            //!< Set the size (voxel count) of a dimension. Valid dimensions are 1-7.
		const Array<size_t, Dynamic, 1> dims() const;           //!< Get all dimension sizes.
		void setDims(const Array<size_t, Dynamic, 1> &newDims); //!< Set all dimension sizes.
		size_t voxelsPerSlice() const;                          //!< Voxel count for a whole slice (dim1 x dim2).
		size_t voxelsPerVolume() const;                         //!< Voxel count for a volume (dim1 x dim2 x dim3).
		size_t voxelsTotal() const;                             //!< Voxel count for whole image (all dimensions).
		
		float voxDim(const size_t d) const;                     //!< Get the voxel size along dimension d. Valid dimensions are 1-7.
		void setVoxDim(const size_t d, const float f);          //!< Set the voxel size along dimension d. Valid dimensions are 1-7.
		const ArrayXf voxDims() const;                          //!< Get all voxel sizes.
		void setVoxDims(const ArrayXf &newVoxDims);             //!< Set all voxel sizes.
		
		const int &datatype() const;
		const string &dtypeName() const;
		const int &bytesPerVoxel() const;
		void setDatatype(const int dt);
		
		float scaling_slope;
		float scaling_inter;
		float calibration_min;
		float calibration_max;
		
		void setTransform(const Affine3f &t, const int transform_code = NIFTI_XFORM_SCANNER_ANAT); //!< Set the qform and sform from a 4x4 general matrix. The qform will be set to closest matching linear transform, the sform will be an exact copy.
		const Affine3f &transform() const;           //!< Return the transform with the highest priority.
		const Affine3f &qform() const;               //!< Return just the qform.
		const Affine3f &sform() const;               //!< Return just the sform.
		int qform_code;
		int sform_code;
		bool matchesSpace(const File &other) const;  //!< Check if voxel dimensions, data size and transform match
		bool matchesVoxels(const File &other) const; //!< Looser check if voxel dimensions and data size match
		
		int freq_dim ;                //!< Index of the frequency encode direction (1-3)
		int phase_dim;                //!< Index of the phase encode direction (1-3)
		int slice_dim;                //!< Index of the slice direction (1-3)
		
		int   slice_code;             //!< code for slice timing pattern
		int   slice_start;            //!< index for start of slices
		int   slice_end;              //!< index for end of slices
		float slice_duration;         //!< time between individual slices
		float toffset;                //!< time coordinate offset
	
		int xyz_units;                //!< dx,dy,dz units: NIFTI_UNITS_* code
		int time_units;               //!< dt       units: NIFTI_UNITS_* code
		int   intent_code ;           //!< statistic type (or something)
		float intent_p1 ;             //!< intent parameters
		float intent_p2 ;             //!< intent parameters
		float intent_p3 ;             //!< intent parameters
		string intent_name;           //!< optional description of intent data
		string description;           //!< optional text to describe dataset
		string aux_file;              //!< auxiliary filename
		
		static const string &TransformName(const int code);
		const string &qformName() const;
		const string &sformName() const;
		const string &spaceUnits() const;
		const string &timeUnits() const;
		const string &intentName() const;
		const string &sliceName() const;
		
		void addExtension(const int code, const vector<char> &data);
		void addExtension(const Extension &e);
		const list<Extension> &extensions() const;
		
		template<typename T> void readVolume(const size_t &vol, vector<T> &buffer);
		template<typename T> vector<T> readVolume(const size_t &vol);
		template<typename T> void readAllVolumes(vector<T> &buffer);
		template<typename T> vector<T> readAllVolumes();
		template<typename T> void readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
		                                        const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
												vector<T> &buffer);
		template<typename T> void readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
		                                        const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et);
		template<typename T> void writeVolume(const size_t vol, const vector<T> &data);
		template<typename T> void writeAllVolumes(const vector<T> &data);
		template<typename T> void writeSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
												 const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
												 const vector<T> &data);
};

#include "Nifti-inl.h"

}; // End namespace Nifti
#endif // NIFTI_IMAGE
