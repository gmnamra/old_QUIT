/** \file NiftiImage.h
 \brief Declaration for NiftiImage class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef NIFTI_IMAGE
#define NIFTI_IMAGE

#include <stdio.h>
#include <zlib.h>

#include <string>
#include <iostream>
#include <algorithm>
#include <complex>
#include <map>

#include "Eigen/Geometry"

#include "nifti1.h" // NIFTI-1 header specification
#include "nifti_analyze.h" // NIFTI version of the ANALYZE 7.5 header

using namespace std;
using namespace Eigen;

// Convenience macros for printing errors. Note that err is NOT encased in ()
// so that NIFTI_ERROR( "string" << number ); works
#define NIFTI_ERROR( err ) do { cerr << __PRETTY_FUNCTION__ << ": " << err << endl; } while(0)
#define NIFTI_FAIL( err ) do { NIFTI_ERROR( err ); exit(EXIT_FAILURE); } while(0)

/*! Utility class that wraps unzipped and zipped files into one object */
// zlib 1.2.5 and above support a "Transparent" mode that would remove the need for this,
// but Mac OS is stuck on 1.2.1
class ZipFile {
	private:
		FILE *_unzipped;
		gzFile _zipped;
		bool _zip;
		static const int MaxZippedBytes = 1 << 30;
		
	public:
		bool open(const string &path, const string &mode, const bool zip);
		void close();
		int read(void *buff, int size);        //!< Attempts to reads size bytes from the image file to buff. Returns actual number read.
		int write(const void *buff, int size); //!< Attempts to write size bytes from buff to the image file. Returns actual number written.
		long seek(long offset, int whence);    //!< Seeks to the specified position in the file
		long tell();                           //!< Returns the current position in the file
		void flush();                          //!< Flushes unwritten buffer contents
};
	
/*! NIfTI header class */
class NiftiImage {
	private:

		struct DataType {
			int code, size, swapsize;
			string name;
		};
		typedef map<int, DataType> DTMap;
		typedef map<int, string> StringMap;
		
		Array<int, 7, 1> _dim;     //!< Number of voxels = nx*ny*nz*...*nw
		Array<float, 7, 1> _voxdim;//!< Dimensions of each voxel
		Affine3f _qform, _sform;   //!< Tranformation matrices from voxel indices to physical co-ords
		
		string _basepath;          //!< Path to file without extension
		bool _nii, _gz;
		char _mode;                //!< Whether the file is closed or open for reading/writing
		ZipFile _file;
		DataType _datatype;        //!< Datatype on disk
		int _voxoffset;            //!< Offset to start of voxel data
		int _swap;                 //!< True if byte order on disk is different to CPU
		
		static int needs_swap(short dim0, int hdrsize); //!< Check if file endianism matches host endianism.
		static float fixFloat(const float f); //!< Converts invalid floats to 0 to ensure a marginally sane header
		
		static void SwapBytes(size_t n, int siz, void *ar);
		static void SwapNiftiHeader(struct nifti_1_header *h);
		static void SwapAnalyzeHeader(nifti_analyze75 *h);
		
		bool readHeader();  //!< Attempts to read a header structure from the currently open file. Returns true on success, false on failure.
		bool writeHeader(); //!< Attempts to write a header structure to the currently open file. Returns true on success, false on failure.
		char *readBytes(size_t start, size_t length, char *buffer = NULL);
		void writeBytes(char *buffer, size_t start, size_t length);

		/**
		  *   Internal function to convert the internal NIfTI data to the desired dataype.
		  *
		  *   Converts a sequence of bytes from a NIfTI image to the templated datatype.
		  *   
		  *   @param T Desired datatype. Valid (scalar) conversions must exist.
		  *   @param bytes Pointer to the NIfTI data
		  *   @param nEl Number of elements (not bytes) in the data
		  *   @param data Pointer to the storage for the converted data. If NULL then
		  *          storage is allocated internally and a pointer returned.
		  *   @return Pointer to the converted data.
		  */
		template<typename T> T* convertFromBytes(char *const bytes, const size_t nEl, T *data = NULL)
		{
			if (!data)
				data = new T[nEl];
			for (int i = 0; i < nEl; i++) {
				switch (_datatype.code) {
					case NIFTI_TYPE_INT8:      data[i] = static_cast<T>(reinterpret_cast<char *>(bytes)[i]); break;
					case NIFTI_TYPE_UINT8:     data[i] = static_cast<T>(reinterpret_cast<unsigned char *>(bytes)[i]); break;
					case NIFTI_TYPE_INT16:     data[i] = static_cast<T>(reinterpret_cast<short *>(bytes)[i]); break;
					case NIFTI_TYPE_UINT16:    data[i] = static_cast<T>(reinterpret_cast<unsigned short *>(bytes)[i]); break;
					case NIFTI_TYPE_INT32:     data[i] = static_cast<T>(reinterpret_cast<int *>(bytes)[i]); break;
					case NIFTI_TYPE_UINT32:    data[i] = static_cast<T>(reinterpret_cast<unsigned int *>(bytes)[i]); break;
					case NIFTI_TYPE_FLOAT32:   data[i] = static_cast<T>(reinterpret_cast<float *>(bytes)[i]); break;
					case NIFTI_TYPE_FLOAT64:   data[i] = static_cast<T>(reinterpret_cast<double *>(bytes)[i]); break;
					case NIFTI_TYPE_INT64:     data[i] = static_cast<T>(reinterpret_cast<long *>(bytes)[i]); break;
					case NIFTI_TYPE_UINT64:    data[i] = static_cast<T>(reinterpret_cast<unsigned long *>(bytes)[i]); break;
					case NIFTI_TYPE_FLOAT128:  data[i] = static_cast<T>(reinterpret_cast<long double *>(bytes)[i]); break;
					// NOTE: C++11 specifies that C++ 'complex<type>' and C 'type complex'
					// should be interchangeable even at pointer level
					case NIFTI_TYPE_COMPLEX64:  data[i] = static_cast<T>(abs(reinterpret_cast<complex<float> *>(bytes)[i])); break;
					case NIFTI_TYPE_COMPLEX128: data[i] = static_cast<T>(abs(reinterpret_cast<complex<double> *>(bytes)[i])); break;
					case NIFTI_TYPE_COMPLEX256: data[i] = static_cast<T>(abs(reinterpret_cast<complex<long double> *>(bytes)[i])); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						NIFTI_FAIL("RGB/RGBA datatypes not supported.");
						break;
				}
			}
			return data;
		}
		
		template<typename T> complex<T> *convertFromBytes(char *const bytes, const size_t nEl, complex<T> *data)
		{
			if (!data)
				data = new complex<T>[nEl];
			for (int i = 0; i < nEl; i++) {
				switch (_datatype.code) {
					case NIFTI_TYPE_INT8:      data[i].real(static_cast<T>(reinterpret_cast<char *>(bytes)[i])); break;
					case NIFTI_TYPE_UINT8:     data[i].real(static_cast<T>(reinterpret_cast<unsigned char *>(bytes)[i])); break;
					case NIFTI_TYPE_INT16:     data[i].real(static_cast<T>(reinterpret_cast<short *>(bytes)[i])); break;
					case NIFTI_TYPE_UINT16:    data[i].real(static_cast<T>(reinterpret_cast<unsigned short *>(bytes)[i])); break;
					case NIFTI_TYPE_INT32:     data[i].real(static_cast<T>(reinterpret_cast<int *>(bytes)[i])); break;
					case NIFTI_TYPE_UINT32:    data[i].real(static_cast<T>(reinterpret_cast<unsigned int *>(bytes)[i])); break;
					case NIFTI_TYPE_FLOAT32:   data[i].real(static_cast<T>(reinterpret_cast<float *>(bytes)[i])); break;
					case NIFTI_TYPE_FLOAT64:   data[i].real(static_cast<T>(reinterpret_cast<double *>(bytes)[i])); break;
					case NIFTI_TYPE_INT64:     data[i].real(static_cast<T>(reinterpret_cast<long *>(bytes)[i])); break;
					case NIFTI_TYPE_UINT64:    data[i].real(static_cast<T>(reinterpret_cast<unsigned long *>(bytes)[i])); break;
					case NIFTI_TYPE_FLOAT128:  data[i].real(static_cast<T>(reinterpret_cast<long double *>(bytes)[i])); break;
					case NIFTI_TYPE_COMPLEX64:  data[i] = static_cast<complex<T> >(reinterpret_cast<complex<float> *>(bytes)[i]); break;
					case NIFTI_TYPE_COMPLEX128: data[i] = static_cast<complex<T> >(reinterpret_cast<complex<double> *>(bytes)[i]); break;
					case NIFTI_TYPE_COMPLEX256: data[i] = static_cast<complex<T> >(reinterpret_cast<complex<long double> *>(bytes)[i]); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						NIFTI_FAIL("RGB/RGBA datatypes not supported.");
						break;
				}
			}
			return data;
		}
		
		/**
		  *   Internal function to convert data to the internal NIfTI datatype.
		  *
		  *   Converts from the templated type to a sequence of bytes suitable for writing to a NIfTI image.
		  *   @param data Pointer to the data.
		  *   @param T Desired datatype. Valid (scalar) conversions must exist.
		  *   @param nEl Number of elements (not bytes) in the data
		  *   @param bytes Pointer to the NIfTI data. If NULL then storage is
		  *          allocated internally and a pointer returned.
		  *   @return Pointer to the data stored in a byte array.
		  */
		template<typename T> char *convertToBytes(const T *const data, const size_t nEl, char *bytes = NULL)
		{
			if (!bytes)
				bytes = new char[nEl * _datatype.size];
			for (int i = 0; i < nEl; i++) {
				switch (_datatype.code) {
					case NIFTI_TYPE_INT8:              reinterpret_cast<char *>(bytes)[i] = static_cast<char>(data[i]); break;
					case NIFTI_TYPE_UINT8:    reinterpret_cast<unsigned char *>(bytes)[i] = static_cast<unsigned char>(data[i]); break;
					case NIFTI_TYPE_INT16:            reinterpret_cast<short *>(bytes)[i] = static_cast<short>(data[i]); break;
					case NIFTI_TYPE_UINT16:  reinterpret_cast<unsigned short *>(bytes)[i] = static_cast<unsigned short>(data[i]); break;
					case NIFTI_TYPE_INT32:              reinterpret_cast<int *>(bytes)[i] = static_cast<int>(data[i]); break;
					case NIFTI_TYPE_UINT32:    reinterpret_cast<unsigned int *>(bytes)[i] = static_cast<unsigned int>(data[i]); break;
					case NIFTI_TYPE_FLOAT32:          reinterpret_cast<float *>(bytes)[i] = static_cast<float>(data[i]); break;
					case NIFTI_TYPE_FLOAT64:         reinterpret_cast<double *>(bytes)[i] = static_cast<double>(data[i]); break;
					case NIFTI_TYPE_INT64:             reinterpret_cast<long *>(bytes)[i] = static_cast<long>(data[i]); break;
					case NIFTI_TYPE_UINT64:   reinterpret_cast<unsigned long *>(bytes)[i] = static_cast<unsigned long>(data[i]); break;
					case NIFTI_TYPE_FLOAT128:   reinterpret_cast<long double *>(bytes)[i] = static_cast<long double>(data[i]); break;
					case NIFTI_TYPE_COMPLEX64:  reinterpret_cast<complex<float> *>(bytes)[i].real(static_cast<float>(data[i])); break;
					case NIFTI_TYPE_COMPLEX128: reinterpret_cast<complex<double> *>(bytes)[i].real(static_cast<double>(data[i])); break;
					case NIFTI_TYPE_COMPLEX256: reinterpret_cast<complex<long double> *>(bytes)[i].real(static_cast<long double>(data[i])); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						NIFTI_FAIL("RGB/RGBA datatypes not supported."); break;
				}
			}
			return bytes;
		}
		
		template<typename T> char *convertToBytes(const complex<T> *const data, const size_t nEl, char *bytes = NULL)
		{
			if (!bytes)
				bytes = new char[nEl * _datatype.size];
			for (int i = 0; i < nEl; i++) {
				switch (_datatype.code) {
					case NIFTI_TYPE_INT8:              reinterpret_cast<char *>(bytes)[i] = static_cast<char>(abs(data[i])); break;
					case NIFTI_TYPE_UINT8:    reinterpret_cast<unsigned char *>(bytes)[i] = static_cast<unsigned char>(abs(data[i])); break;
					case NIFTI_TYPE_INT16:            reinterpret_cast<short *>(bytes)[i] = static_cast<short>(abs(data[i])); break;
					case NIFTI_TYPE_UINT16:  reinterpret_cast<unsigned short *>(bytes)[i] = static_cast<unsigned short>(abs(data[i])); break;
					case NIFTI_TYPE_INT32:              reinterpret_cast<int *>(bytes)[i] = static_cast<int>(abs(data[i])); break;
					case NIFTI_TYPE_UINT32:    reinterpret_cast<unsigned int *>(bytes)[i] = static_cast<unsigned int>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT32:          reinterpret_cast<float *>(bytes)[i] = static_cast<float>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT64:         reinterpret_cast<double *>(bytes)[i] = static_cast<double>(abs(data[i])); break;
					case NIFTI_TYPE_INT64:             reinterpret_cast<long *>(bytes)[i] = static_cast<long>(abs(data[i])); break;
					case NIFTI_TYPE_UINT64:   reinterpret_cast<unsigned long *>(bytes)[i] = static_cast<unsigned long>(abs(data[i])); break;
					case NIFTI_TYPE_FLOAT128:   reinterpret_cast<long double *>(bytes)[i] = static_cast<long double>(abs(data[i])); break;
					case NIFTI_TYPE_COMPLEX64:  reinterpret_cast<complex<float> *>(bytes)[i] = static_cast<complex<float> >(data[i]); break;
					case NIFTI_TYPE_COMPLEX128: reinterpret_cast<complex<double> *>(bytes)[i] = static_cast<complex<double> >(data[i]); break;
					case NIFTI_TYPE_COMPLEX256: reinterpret_cast<complex<long double> *>(bytes)[i] = static_cast<complex<long double> >(data[i]); break;
					case NIFTI_TYPE_RGB24: case NIFTI_TYPE_RGBA32:
						NIFTI_FAIL("RGB/RGBA datatypes not supported."); break;
				}
			}
			return bytes;
		}
		
	public:
		/*
		 *  Used when opening a NiftiImage to specify read or write. NiftiImages
		 *  can only be open for either reading or writing at any one time, and
		 *  must be closed before re-opening. READ_HEADER is a special mode that
		 *  will just read the header and then close the file.
		 *
		 */
		enum IMAGE_MODES
		{
			CLOSED = 0,
			READ = 'r',
			WRITE = 'w',
			READ_HEADER = 'h'
		};
		
		~NiftiImage();
		NiftiImage();
		NiftiImage(const int nx, const int ny, const int nz, const int nt,
		           const float dx, const float dy, const float dz, const float dt,
				   const int datatype);
		NiftiImage(const ArrayXi &dim, const ArrayXf &voxdim, const int &datatype,
                   const Matrix4f &qform = Matrix4f::Identity(), const Matrix4f &sform = Matrix4f::Identity());
		NiftiImage(const string &filename, const char &mode);
		NiftiImage &operator=(const NiftiImage &other);
		static void printDTypeList();
		
		bool open(const string &filename, const char &mode); //!< Attempts to open a NIfTI file. Returns true on success, false on failure.
		void close();
		const string basePath() const;
		const string imagePath() const;
		const string headerPath() const;
		char *readRawVolume(const int vol);
		char *readRawAllVolumes();
		
		int dimensions() const;                //!< Get the number of dimensions (rank) of this image
		void setDimensions(const ArrayXi &dims, const ArrayXf &voxDims); //!< Set all dimension information in one go
		
		int dim(const int d) const;             //!< Get the size (voxel count) of a dimension
		void setDim(const int d, const int n);  //!< Set the size (voxel count) of a dimension d
		const ArrayXi dims() const;            //!< Get all dimension sizes
		void setDims(const ArrayXi &newDims);   //!< Set all dimension sizes
		int voxelsPerSlice() const;             //!< Voxel count for a whole slice (dim1 x dim2)
		int voxelsPerVolume() const;            //!< Voxel count for a volume (dim1 x dim2 x dim3)
		int voxelsTotal() const;                //!< Voxel count for whole image (all dimensions)
		
		float voxDim(const int d) const;            //!< Get the voxel size along dimension d
		void setVoxDim(const int d, const float f); //!< Set the voxel size along dimension d
		const ArrayXf voxDims() const;             //!< Get all voxel sizes
		void setVoxDims(const ArrayXf &newVoxDims); //!< Set all voxel sizes
		
		const int &datatype() const;
		const string &dtypeName() const;
		const int &bytesPerVoxel() const;
		void setDatatype(const int dt);
		bool matchesSpace(const NiftiImage &other) const; //!< Check if voxel dimensions, data size and transform match
		bool matchesVoxels(const NiftiImage &other) const; //!< Looser check if voxel dimensions and data size match
		
		float scaling_slope;
		float scaling_inter;
		float calibration_min;
		float calibration_max;
		
		int qform_code;
		int sform_code;
		const Matrix4f &qform() const;
		const Matrix4f &sform() const;
		const Matrix4f &ijk_to_xyz() const;
		const Matrix4f &xyz_to_ijk() const;
		
		int freq_dim ;                //!< indexes (1,2,3, or 0) for MRI
		int phase_dim;                //!< directions in dim[]/pixdim[]
		int slice_dim;                //!< directions in dim[]/pixdim[]
		
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
		string intent_name;      //!< optional description of intent data
		string description;      //!< optional text to describe dataset
		string aux_file;         //!< auxiliary filename
		
		static const string &TransformName(const int code);
		const string &qformName() const;
		const string &sformName() const;
		const string &spaceUnits() const;
		const string &timeUnits() const;
		const string &intentName() const;
		const string &sliceName() const;
		
		template<typename T> T *readVolume(const int &vol, T *converted = NULL)
		{
			size_t bytesPerVolume = voxelsPerVolume() * _datatype.size;
			char *raw = readBytes(vol * bytesPerVolume, bytesPerVolume);
			converted = convertFromBytes<T>(raw, voxelsPerVolume(), converted);
			delete[] raw;
			return converted;
		}
		
		template<typename T> T *readAllVolumes()
		{
			char *raw =	readBytes(0, voxelsTotal() * _datatype.size);
			T *converted = convertFromBytes<T>(raw, voxelsTotal());
			delete[] raw;
			return converted;
		}
		
		template<typename T> T *readSubvolume(const int &sx, const int &sy, const int &sz, const int &st,
		                                      const int &ex, const int &ey, const int &ez, const int &et,
											  T *converted = NULL)
		{
			size_t lx, ly, lz, lt, total, toRead;
			lx = ((ex == -1) ? dim(1) : ex) - sx;
			ly = ((ey == -1) ? dim(2) : ey) - sy;
			lz = ((ez == -1) ? dim(3) : ez) - sz;
			lt = ((et == -1) ? dim(4) : et) - st;
			total = lx * ly * lz * lt;
			
			if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
				NIFTI_ERROR("Invalid subvolume read of dimensions " <<
							lx << "," << ly << "," << lz << "," << lt << " requested. Nothing read.");
				return NULL;
			}
			
			// Collapse successive full dimensions into a single compressed read
			toRead = lx * _datatype.size;
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
						
			char *raw = new char[total * _datatype.size];
			char *nextRead = raw;
			for (int t = st; t < st+lt; t++)
			{
				size_t tOff = t * voxelsPerVolume();
				for (int z = sz; z < sz+lz; z++)
				{
					size_t zOff = z * voxelsPerSlice();
					for (int y = sy; y < sy+ly; y++)
					{
						size_t yOff = y * dim(1);
						if (readBytes((tOff + zOff + yOff) * _datatype.size, toRead, nextRead))
							nextRead += toRead;
						else
							NIFTI_FAIL("failed to read subvolume from file.");
					}
				}
			}
			converted = convertFromBytes<T>(raw, total, converted);
			delete[] raw;
			return converted;
		}
		
		template<typename T> void writeVolume(const int vol, const T *data)
		{
			size_t bytesPerVolume = voxelsPerVolume() * _datatype.size;
			char *converted = convertToBytes(data, voxelsPerVolume());
			writeBytes(converted, vol * bytesPerVolume, bytesPerVolume);
			delete[] converted;
		}
		
		template<typename T> void writeAllVolumes(const T *data)
		{
			char *converted = convertToBytes(data, voxelsTotal());
			writeBytes(converted, 0, voxelsTotal() * _datatype.size);
			delete[] converted;
		}
		
		template<typename T> void writeSubvolume(const int &sx, const int &sy, const int &sz, const int &st,
												 const int &ex, const int &ey, const int &ez, const int &et,
											     const T *data)
		{
			size_t lx, ly, lz, lt, total, toWrite;
			lx = ((ex == -1) ? dim(1) : ex) - sx;
			ly = ((ey == -1) ? dim(2) : ey) - sy;
			lz = ((ez == -1) ? dim(3) : ez) - sz;
			lt = ((et == -1) ? dim(4) : et) - st;
			total = lx * ly * lz * lt;
			
			if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
				NIFTI_ERROR("Invalid subvolume write of dimensions " <<
							lx << "," << ly << "," << lz << "," << lt << " requested. Nothing written.");
				return;
			}
			
			// Collapse successive full dimensions into a single write
			toWrite = lx * _datatype.size;
			if (lx == dim(1)) {
				toWrite *= ly;
				if (ly == dim(2)) {
					toWrite *= lz;
					if (lz == dim(3)) {
						// If we've got to here we're actual writing the whole image
						toWrite *= lt;
						lt = 1;
					}
					lz = 1;
				}
				ly = 1;
			}
						
			char *raw = convertToBytes(data, total);
			char *nextWrite = raw;
			for (int t = st; t < st+lt; t++)
			{
				size_t tOff = t * voxelsPerVolume();
				for (int z = sz; z < sz+lz; z++)
				{
					size_t zOff = z * voxelsPerSlice();
					for (int y = sy; y < sy+ly; y++)
					{
						size_t yOff = y * dim(1);
						writeBytes(nextWrite, (tOff + zOff + yOff) * _datatype.size, toWrite);
						nextWrite += toWrite;
					}
				}
			}
			delete[] raw;
		}
};

#endif // NIFTI_IMAGE
