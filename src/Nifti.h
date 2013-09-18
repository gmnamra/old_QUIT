/** \file File.h
 \brief Declaration for File class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef NIFTI_IMAGE
#define NIFTI_IMAGE

#include <zlib.h>

#include <string>
#include <iostream>
#include <cstdio>
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
enum class Modes : char
{
	Closed = 0,
	Read = 'r',
	ReadHeader = 'h',
	ReadSkipExt = 's',
	Write = 'w',
	WriteSkipExt = 'x'
};



#pragma mark ZipFile
/*! Utility class that wraps unzipped and zipped files into one object */
// zlib 1.2.5 and above support a "Transparent" mode that would remove the need for this,
// but Mac OS is stuck on 1.2.1
class ZipFile {
	private:
		FILE *m_plainFile;
		gzFile m_gzipFile;
		
	public:
		ZipFile();
		bool open(const string &path, const string &mode, const bool zip);
		void close();
		size_t read(void *buff, unsigned size);   //!< Attempts to reads size bytes from the image file to buff. Returns actual number read.
		size_t write(const void *buff, int size); //!< Attempts to write size bytes from buff to the image file. Returns actual number written.
		bool seek(long offset, int whence);       //!< Seeks to the specified position in the file. Returns true if successful.
		long tell() const;                        //!< Returns the current position in the file
		void flush();                             //!< Flushes unwritten buffer contents
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

		/**
		  *   Internal function to convert the internal NIfTI data to the desired dataype.
		  *
		  *   Converts a sequence of bytes from a NIfTI image to the templated datatype.
		  *   If the file contains complex data but a real datatype is requested then the magnitude is returned.
		  *   If the file contains real data but a complex datatype is requested then the imaginary part will be 0.
		  *   
		  *   @param T Desired datatype. Valid (scalar) conversions must exist.
		  *   @param bytes std::vector of byte data
		  *   @param nEl Number of elements (not bytes) expected in the data
		  *   @param data std::vector& to converted data in. Will be resized to ensure enough space.
		  */
		template<typename T> void convertFromBytes(const vector<char> &bytes, const size_t nEl, vector<T> &data) {
			assert(nEl == (bytes.size() / bytesPerVoxel()));
			data.resize(nEl);
			for (size_t i = 0; i < nEl; i++) {
				switch (m_datatype.code) {
					case NIFTI_TYPE_INT8:      data[i] = static_cast<T>(reinterpret_cast<const char *>(bytes.data())[i]); break;
					case NIFTI_TYPE_UINT8:     data[i] = static_cast<T>(reinterpret_cast<const unsigned char *>(bytes.data())[i]); break;
					case NIFTI_TYPE_INT16:     data[i] = static_cast<T>(reinterpret_cast<const short *>(bytes.data())[i]); break;
					case NIFTI_TYPE_UINT16:    data[i] = static_cast<T>(reinterpret_cast<const unsigned short *>(bytes.data())[i]); break;
					case NIFTI_TYPE_INT32:     data[i] = static_cast<T>(reinterpret_cast<const int *>(bytes.data())[i]); break;
					case NIFTI_TYPE_UINT32:    data[i] = static_cast<T>(reinterpret_cast<const unsigned int *>(bytes.data())[i]); break;
					case NIFTI_TYPE_FLOAT32:   data[i] = static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]); break;
					case NIFTI_TYPE_FLOAT64:   data[i] = static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]); break;
					case NIFTI_TYPE_INT64:     data[i] = static_cast<T>(reinterpret_cast<const long *>(bytes.data())[i]); break;
					case NIFTI_TYPE_UINT64:    data[i] = static_cast<T>(reinterpret_cast<const unsigned long *>(bytes.data())[i]); break;
					case NIFTI_TYPE_FLOAT128:  data[i] = static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]); break;
					// NOTE: C++11 specifies that C++ 'complex<type>' and C 'type complex'
					// should be interchangeable even at pointer level
					case NIFTI_TYPE_COMPLEX64:  data[i] = static_cast<T>(abs(reinterpret_cast<const complex<float> *>(bytes.data())[i])); break;
					case NIFTI_TYPE_COMPLEX128: data[i] = static_cast<T>(abs(reinterpret_cast<const complex<double> *>(bytes.data())[i])); break;
					case NIFTI_TYPE_COMPLEX256: data[i] = static_cast<T>(abs(reinterpret_cast<const complex<long double> *>(bytes.data())[i])); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						throw(runtime_error("RGB/RGBA datatypes not supported.")); break;
				}
			}
		}
		
		template<typename T> void convertFromBytes(const vector<complex<T>> &bytes, const size_t nEl, vector<T> &data) {
			assert(nEl == (bytes.size() / bytesPerVoxel()));
			data.resize(nEl);
			for (size_t i = 0; i < nEl; i++) {
				switch (m_datatype.code) {
					case NIFTI_TYPE_INT8:      data[i] = complex<T>(static_cast<T>(reinterpret_cast<const char *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_UINT8:     data[i] = complex<T>(static_cast<T>(reinterpret_cast<const unsigned char *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_INT16:     data[i] = complex<T>(static_cast<T>(reinterpret_cast<const short *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_UINT16:    data[i] = complex<T>(static_cast<T>(reinterpret_cast<const unsigned short *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_INT32:     data[i] = complex<T>(static_cast<T>(reinterpret_cast<const int *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_UINT32:    data[i] = complex<T>(static_cast<T>(reinterpret_cast<const unsigned int *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_FLOAT32:   data[i] = complex<T>(static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_FLOAT64:   data[i] = complex<T>(static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_INT64:     data[i] = complex<T>(static_cast<T>(reinterpret_cast<const long *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_UINT64:    data[i] = complex<T>(static_cast<T>(reinterpret_cast<const unsigned long *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_FLOAT128:  data[i] = complex<T>(static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]), 0.); break;
					case NIFTI_TYPE_COMPLEX64:  data[i] = static_cast<complex<T> >(reinterpret_cast<const complex<float> *>(bytes.data())[i]); break;
					case NIFTI_TYPE_COMPLEX128: data[i] = static_cast<complex<T> >(reinterpret_cast<const complex<double> *>(bytes.data())[i]); break;
					case NIFTI_TYPE_COMPLEX256: data[i] = static_cast<complex<T> >(reinterpret_cast<const complex<long double> *>(bytes.data())[i]); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						throw(runtime_error("RGB/RGBA datatypes not supported.")); break;
				}
			}
		}
		
		/**
		  *   Internal function to convert data to the internal NIfTI datatype.
		  *
		  *   Converts from the templated type to a sequence of bytes suitable for writing to a NIfTI image.
		  *   @param data std::vector of the data.
		  *   @param T Desired datatype. Valid (scalar) conversions must exist.
		  *   @return std::vector of the data stored in a byte array.
		  */
		template<typename T> vector<char> convertToBytes(const vector<T> &data) {
			vector<char> bytes(data.size() * bytesPerVoxel());
			for (size_t i = 0; i < data.size(); i++) {
				switch (m_datatype.code) {
					case NIFTI_TYPE_INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<char>(data[i]); break;
					case NIFTI_TYPE_UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<unsigned char>(data[i]); break;
					case NIFTI_TYPE_INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<short>(data[i]); break;
					case NIFTI_TYPE_UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<unsigned short>(data[i]); break;
					case NIFTI_TYPE_INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int>(data[i]); break;
					case NIFTI_TYPE_UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<unsigned int>(data[i]); break;
					case NIFTI_TYPE_FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(data[i]); break;
					case NIFTI_TYPE_FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(data[i]); break;
					case NIFTI_TYPE_INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<long>(data[i]); break;
					case NIFTI_TYPE_UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<unsigned long>(data[i]); break;
					case NIFTI_TYPE_FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(data[i]); break;
					case NIFTI_TYPE_COMPLEX64:  reinterpret_cast<complex<float> *>(bytes.data())[i] = complex<float>(static_cast<float>(data[i])); break;
					case NIFTI_TYPE_COMPLEX128: reinterpret_cast<complex<double> *>(bytes.data())[i] = complex<double>(static_cast<double>(data[i])); break;
					case NIFTI_TYPE_COMPLEX256: reinterpret_cast<complex<long double> *>(bytes.data())[i] = complex<long double>(static_cast<long double>(data[i])); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						throw(runtime_error("RGB/RGBA datatypes not supported.")); break;				}
			}
			return bytes;
		}
		
		template<typename T> vector<char> convertToBytes(const vector<complex<T>> &data) {
			vector<char> bytes(data.size() * bytesPerVoxel());
			for (int i = 0; i < data.size(); i++) {
				switch (m_datatype.code) {
					case NIFTI_TYPE_INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<char>(abs(data[i])); break;
					case NIFTI_TYPE_UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<unsigned char>(abs(data[i])); break;
					case NIFTI_TYPE_INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<short>(abs(data[i])); break;
					case NIFTI_TYPE_UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<unsigned short>(abs(data[i])); break;
					case NIFTI_TYPE_INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int>(abs(data[i])); break;
					case NIFTI_TYPE_UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<unsigned int>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(abs(data[i])); break;
					case NIFTI_TYPE_INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<long>(abs(data[i])); break;
					case NIFTI_TYPE_UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<unsigned long>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(abs(data[i])); break;
					case NIFTI_TYPE_COMPLEX64:  reinterpret_cast<complex<float> *>(bytes.data())[i] = static_cast<complex<float> >(data[i]); break;
					case NIFTI_TYPE_COMPLEX128: reinterpret_cast<complex<double> *>(bytes.data())[i] = static_cast<complex<double> >(data[i]); break;
					case NIFTI_TYPE_COMPLEX256: reinterpret_cast<complex<long double> *>(bytes.data())[i] = static_cast<complex<long double> >(data[i]); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						throw(runtime_error("RGB/RGBA datatypes not supported.")); break;
				}
			}
			return bytes;
		}
	
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
		
		template<typename T> void readVolume(const size_t &vol, vector<T> &buffer) {
			size_t bytesPerVolume = voxelsPerVolume() * bytesPerVoxel();
			vector<char> raw(bytesPerVolume);
			readBytes(vol * bytesPerVolume, bytesPerVolume, raw.data());
			convertFromBytes(raw, voxelsPerVolume(), buffer);
		}
		
		template<typename T> vector<T> readVolume(const size_t &vol) {
			vector<T> buffer;
			readVolume(vol, buffer);
			return buffer;
		}
		
		template<typename T> void readAllVolumes(vector<T> &buffer) {
			vector<char> raw(voxelsTotal() * bytesPerVoxel());
			readBytes(0, voxelsTotal() * bytesPerVoxel(), raw.data());
			return convertFromBytes<T>(raw, voxelsTotal(), buffer);
		}
		
		template<typename T> vector<T> readAllVolumes() {
			vector<T> buffer;
			readAllVolumes(buffer);
			return buffer;
		}
		
		/**
		  *   Reads a subvolume
		  *
		  *   Converts to the templated type from the NIfTI datatype.
		  *   @param sx, sy, sz, st - Start voxels for each dimension.
		  *   @param ex, ey, ez, et - End voxels for each dimension. Values are clamped to dimension sizes.
		  *   @param buffer std::vector to store the data.
		  *
		  */
		template<typename T> void readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
		                                        const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
												vector<T> &buffer) {
			size_t lx, ly, lz, lt, total, toRead;
			lx = ((ex > dim(1)) ? dim(1) : ex) - sx;
			ly = ((ey > dim(2)) ? dim(2) : ey) - sy;
			lz = ((ez > dim(3)) ? dim(3) : ez) - sz;
			lt = ((et > dim(4)) ? dim(4) : et) - st;
			total = lx * ly * lz * lt;
			
			if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
				throw(out_of_range("Invalid subvolume read dimensions: " + imagePath()));
			}
			
			// Collapse successive full dimensions into a single compressed read
			toRead = lx * bytesPerVoxel();
			if (lx == dim(1)) {
				toRead *= ly;
				if (ly == dim(2)) {
					toRead *= lz;
					if (lz == dim(3)) {
						// If we've got to here we're actual reading the whole image
						toRead *= lt;
						lt = 1;
					}
					lz = 1;
				}
				ly = 1;
			}
						
			vector<char> raw(total * bytesPerVoxel());
			char *nextRead = raw.data();
			for (size_t t = st; t < st+lt; t++) {
				size_t tOff = t * voxelsPerVolume();
				for (size_t z = sz; z < sz+lz; z++) {
					size_t zOff = z * voxelsPerSlice();
					for (size_t y = sy; y < sy+ly; y++) {
						size_t yOff = y * dim(1);
						if (readBytes((tOff + zOff + yOff) * bytesPerVoxel(), toRead, nextRead))
							nextRead += toRead;
					}
				}
			}
			convertFromBytes<T>(raw, total, buffer);
		}
		
		template<typename T> void readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
		                                        const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et) {
			vector<T> buffer;
			readSubvolume(sx, sy, sz, st, ex, ey, ez, et, buffer);
			return buffer;
		}
		
		template<typename T> void writeVolume(const size_t vol, const vector<T> &data) {
			if (data.size() != voxelsPerVolume()) {
				throw(invalid_argument("Insufficient data to write volume for file: " + imagePath()));
			}
			vector<char> converted = convertToBytes(data);
			writeBytes(vol * voxelsPerVolume() * bytesPerVoxel(), converted.size(), converted.data());
		}
		
		template<typename T> void writeAllVolumes(const vector<T> &data) {
			if (data.size() != voxelsTotal()) {
				throw(invalid_argument("Insufficient data to write all volumes for file: " + imagePath()));
			}
			vector<char> converted = convertToBytes(data);
			writeBytes(0, converted.size(), converted.data());
		}
		
		/**
		  *   Writes a subvolume
		  *
		  *   Converts from the templated type to a sequence of bytes suitable for writing to a NIfTI image.
		  *   @param sx, sy, sz, st - Start voxels for each dimension.
		  *   @param ex, ey, ez, et - End voxels for each dimension. Values are clamped to dimension sizes.
		  *   @param data std::vector of the data.
		  *
		  */
		template<typename T> void writeSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
												 const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
												 const vector<T> &data) {
			size_t lx, ly, lz, lt, total, toWrite;
			lx = ((ex > dim(1)) ? dim(1) : ex) - sx;
			ly = ((ey > dim(2)) ? dim(2) : ey) - sy;
			lz = ((ez > dim(3)) ? dim(3) : ez) - sz;
			lt = ((et > dim(4)) ? dim(4) : et) - st;
			total = lx * ly * lz * lt;
			
			if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
				throw(out_of_range("Invalid subvolume read dimensions: " + imagePath()));
			}
			
			// Collapse successive full dimensions into a single write
			toWrite = lx;
			if (lx == dim(1)) {
				toWrite *= ly;
				if (ly == dim(2)) {
					toWrite *= lz;
					if (lz == dim(3)) {
						// If we've got to here we're actually writing the whole image
						toWrite *= lt;
						lt = 1;
					}
					lz = 1;
				}
				ly = 1;
			}
			if (toWrite != data.size()) {
				throw(invalid_argument("Insufficient data to write subvolume for file: " + imagePath()));
			}
			toWrite *= bytesPerVoxel();
			
			vector<char> raw = convertToBytes(data);
			char *nextWrite = raw.data();
			for (size_t t = st; t < st+lt; t++) {
				size_t tOff = t * voxelsPerVolume();
				for (size_t z = sz; z < sz+lz; z++) {
					size_t zOff = z * voxelsPerSlice();
					for (size_t y = sy; y < sy+ly; y++) {
						size_t yOff = y * dim(1);
						writeBytes((tOff + zOff + yOff) * bytesPerVoxel(), toWrite, nextWrite);
						nextWrite += toWrite;
					}
				}
			}
		}
};

}; // End namespace Nifti
#endif // NIFTI_IMAGE
