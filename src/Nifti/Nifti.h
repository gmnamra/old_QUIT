/** \file Nifti.h
 \brief Declaration for File class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef NIFTI_NIFTI
#define NIFTI_NIFTI

#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <vector>
#include <list>
#include <map>
#include <limits>
#include <exception>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "ZipFile.h"

#pragma mark NIfTI File Class
class Nifti {
	public:
		typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXs;
		enum class Mode : char {
			Closed = 0, Read = 'r', ReadHeader = 'h', ReadSkipExt = 's', Write = 'w', WriteSkipExt = 'x'
		};

		enum class DataType {
			UINT8, UINT16, UINT32, UINT64, INT8, INT16, INT32, INT64,
			FLOAT32, FLOAT64, FLOAT128, COMPLEX64, COMPLEX128, COMPLEX256,
			RGB24, RGBA32
		};
		
		struct DataTypeInfo {
			DataType type;
			size_t code, size, swapsize;
			std::string name;
		}; //!< Contains all the information needed to read/write a Nifti datatype
		static const DataTypeInfo &TypeInfo(const DataType dt);
		
		enum class XForm {
			Unknown, ScannerAnatomy, AlignedAnatomy, Talairach, MNI_152
		};
		static const std::string XFormName(const XForm t);
		
		/*
		 *  Nifti Extension Class.
		 *
		 *  Provides a minimal way to read and write Nifti extensions as
		 *  vectors of bytes.
		 */
		class Extension {
			private:
				int m_code;          //!< Extension code, one of the NIFTI_ECODE_ values
				std::vector<char> m_data; //!< Raw data, with no byte swapping (length is esize-8)
			
			public:
				static const std::string &CodeName(const int code);
				
				Extension(int code, std::vector<char> data);
				Extension(int size, int code, char *data);
				const int rawSize() const;
				const int size() const;
				const int padding() const;
				const int code() const;
				const std::string &codeName() const;
				void setCode(int code);
				
				const std::vector<char> &data() const;
				void setData(const std::vector<char> &data);
		};

	private:
		static const DataType DataTypeForCode(const int code);
		static const XForm XFormForCode(const int code);
		static const int XFormCode(const XForm t);
		
		Eigen::Array<size_t, 7, 1> m_dim;      //!< Number of voxels in each dimension. Note that here we do NOT store the rank in dim[0], so only 7 elements required.
		Eigen::Array<size_t, 7, 1> m_strides;  //!< Strides into the data on disk.
		Eigen::Array<float, 7, 1> m_voxdim;    //!< Size of each voxel. As above, only 7 elements because the rank is not stored.
		Eigen::Affine3f m_qform, m_sform;      //!< Tranformation matrices from voxel indices to physical co-ords.
		XForm m_qcode, m_scode;                //!< Codes to define what the transformations represent.
		std::string m_basepath;                //!< Path to file without extension.
		bool m_nii, m_gz;
		Mode m_mode;                    //!< Whether the file is closed or open for reading/writing.
		ZipFile m_file;
		DataTypeInfo m_typeinfo;        //!< Informatio for datatype on disk.
		int m_voxoffset;                //!< Offset to start of voxel data.
		int m_swap;                     //!< True if byte order on disk is different to CPU.
		
		std::list<Extension> m_extensions;
				
		void readHeader();      //!< Attempts to read a header structure from the currently open file.
		void readExtensions();  //!< Attempts to read any extensions
		void writeHeader();     //!< Attempts to write a header structure to the currently open file.
		void writeExtensions(); //!< Attempts to write extensions
		int totalExtensionSize(); //!< Counts the total number of bytes for all extensions.
		char *readBytes(size_t start, size_t length, char *buffer);
		void readBytes(std::vector<char> &data);
		void writeBytes(size_t start, size_t length, char *buffer);
		void writeBytes(const std::vector<char> & data);
		void calcStrides();
		void seekToVoxel(const ArrayXs &target);
		
		template<typename T> void convertFromBytes(const std::vector<char> &bytes, const typename std::vector<T>::iterator begin, const typename std::vector<T>::iterator end);
		template<typename T> void convertFromBytes(const std::vector<char> &bytes, const typename std::vector<std::complex<T>>::iterator begin, const typename std::vector<std::complex<T>>::iterator end);
		template<typename T> void convertToBytes(const typename std::vector<T>::iterator begin, const typename std::vector<T>::iterator end, std::vector<char> &bytes);
		template<typename T> void convertToBytes(const typename std::vector<std::complex<T>>::iterator begin, const typename std::vector<std::complex<T>>::iterator end, std::vector<char> &bytes);
				
	#pragma mark Public Class Methods
	public:
		~Nifti();
		Nifti();                               //!< Default constructor. Initialises an empty header, size 1 in all dimensions.
		Nifti(const Nifti &other);             //!< Copy constructor. Copies all elements, and if the original is open then also opens new file handles.
		Nifti &operator=(const Nifti &other);  //!< Assignment. Copies all elements except file handles, and marks destination as Closed.
		Nifti(Nifti &&other) noexcept;         //!< Move constructor. Copies all elements, including the file handles, and marks the original as Closed.
		
		Nifti(const int nx, const int ny, const int nz, const int nt,
			  const float dx, const float dy, const float dz, const float dt,
			  const DataType dtype = DataType::FLOAT32, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		Nifti(const ArrayXs &dim, const Eigen::ArrayXf &voxdim,
			  const DataType dtype = DataType::FLOAT32, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		Nifti(const Nifti &other, const size_t nt, const DataType dtype = DataType::FLOAT32);               //!< Copies only basic geometry information from other, then sets the datatype and number of volumes. Does not copy scaling information etc.
		Nifti(const std::string &filename, const Mode &mode);
		
		void open(const std::string &filename, const Mode &mode); //!< Attempts to open a NIfTI file. Throws runtime_error or invalid_argument on failure.
		void close();                                             //!< Closes the file
		bool isOpen();                                            //!< Returns true if file is currently open for reading or writing.
		
		const std::string basePath() const;
		const std::string imagePath() const;
		const std::string headerPath() const;
		char *readRawVolume(const int vol);
		char *readRawAllVolumes();
		
		size_t dimensions() const;                              //!< Get the number of dimensions (rank) of this image.
		size_t dim(const size_t d) const;                       //!< Get the size (voxel count) of a dimension. Valid dimensions are 1-7.
		void setDim(const size_t d, const size_t n);            //!< Set the size (voxel count) of a dimension. Valid dimensions are 1-7.
		const ArrayXs dims() const;                             //!< Get all dimension sizes.
		void setDims(const ArrayXs &newDims);                   //!< Set all dimension sizes.
		size_t voxelsPerSlice() const;                          //!< Voxel count for a whole slice (dim1 x dim2).
		size_t voxelsPerVolume() const;                         //!< Voxel count for a volume (dim1 x dim2 x dim3).
		size_t voxelsTotal() const;                             //!< Voxel count for whole image (all dimensions).
		
		float voxDim(const size_t d) const;                     //!< Get the voxel size along dimension d. Valid dimensions are 1-7.
		void setVoxDim(const size_t d, const float f);          //!< Set the voxel size along dimension d. Valid dimensions are 1-7.
		const Eigen::ArrayXf voxDims() const;                   //!< Get all voxel sizes.
		void setVoxDims(const Eigen::ArrayXf &newVoxDims);      //!< Set all voxel sizes.
		
		const DataType &datatype() const;
		void setDatatype(const DataType dt);
		
		float scaling_slope;
		float scaling_inter;
		float calibration_min;
		float calibration_max;
		
		void setTransform(const Eigen::Affine3f &t, const XForm tc = XForm::ScannerAnatomy); //!< Set the qform and sform from a 4x4 general matrix. The qform will be set to closest matching linear XForm, the sform will be an exact copy.
		const Eigen::Affine3f &transform() const;            //!< Return the XForm with the highest priority.
		const Eigen::Affine3f &qform() const;                //!< Return just the qform.
		const Eigen::Affine3f &sform() const;                //!< Return just the sform.
		const XForm qcode() const;               //!< Find out what transformation the qform represents.
		const XForm scode() const;               //!< Find out what transformation the sform represents.
		
		bool matchesSpace(const Nifti &other) const;  //!< Check if voxel dimensions, data size and XForm match
		bool matchesVoxels(const Nifti &other) const; //!< Looser check if voxel dimensions and data size match
		
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
		std::string intent_name;      //!< optional description of intent data
		std::string description;      //!< optional text to describe dataset
		std::string aux_file;         //!< auxiliary filename
		
		const std::string &spaceUnits() const;
		const std::string &timeUnits() const;
		const std::string &intentName() const;
		const std::string &sliceName() const;
		
		void addExtension(const int code, const std::vector<char> &data);
		void addExtension(const Extension &e);
		const std::list<Extension> &extensions() const;
		
		template<typename T> void readWriteVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size, std::vector<T> &data);
		template<typename T> void readVolumes(const size_t first, const size_t nvol, std::vector<T> &data);
		template<typename T> void writeVolumes(const size_t vol, const size_t nvol, const std::vector<T> &data);
};

#include "Nifti-inl.h"

#endif // NIFTI_IMAGE
