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
#include <stdexcept>
#include <cstdint>
#include <type_traits>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "ZipFile.h"

#pragma mark NIfTI File Class
class Nifti {
	public:
		typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXs;
		enum class Mode : char {
			Closed = 0, Read = 'r', ReadHeader = 'h', Write = 'w'
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
				size_t rawSize() const;
				int size() const;
				int padding() const;
				int code() const;
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
		
		void setPaths(const std::string &path); //!< Works out the basepath and file extensions.
		void readHeader();                 //!< Attempts to read a header structure from the currently open file.
		void readExtensions();             //!< Attempts to read any extensions
		void writeHeader();                //!< Attempts to write a header structure to the currently open file.
		void writeExtensions();            //!< Attempts to write extensions
		int totalExtensionSize();          //!< Counts the total number of bytes for all extensions.

		void calcStrides();
		
		void seekToVoxel(const ArrayXs &target);
		void readBytes(std::vector<char> &data);
		void writeBytes(const std::vector<char> & data);
		
		template<typename FromTp, typename ToTp>
		class Scale {
		public:
			static ToTp Forward(const FromTp val, const float slope, const float inter) {
				return static_cast<ToTp>((val * slope) + inter);
			}
			static ToTp Reverse(const FromTp val, const float slope, const float inter) {
				return static_cast<ToTp>((val - inter) / slope);
			}
		};
		
		template<typename FromTp, typename ToTp>
		class Scale<const std::complex<FromTp>, ToTp> {
		public:
			static ToTp Forward(const std::complex<FromTp> val, const float slope, const float inter) {
				return static_cast<ToTp>((std::abs(val) * slope) + inter);
			}
			static ToTp Reverse(const std::complex<FromTp> val, const float slope, const float inter) {
				return static_cast<ToTp>((std::abs(val) - inter) / slope);
			}
		};
		
		template<typename FromTp, typename ToTp>
		class Scale<const FromTp, std::complex<ToTp>> {
		public:
			static std::complex<ToTp> Forward(const FromTp val, const float slope, const float inter) {
				return std::complex<ToTp>((val * slope) + inter, 0.);
			}
			static std::complex<ToTp> Reverse(const FromTp val, const float slope, const float inter) {
				return std::complex<ToTp>((val - inter) / slope, 0.);
			}
		};
		
		template<typename FromTp, typename ToTp>
		class Scale<const std::complex<FromTp>, std::complex<ToTp>> {
		public:
			static std::complex<ToTp> Forward(const std::complex<FromTp> val, const float slope, const float inter) {
				return static_cast<std::complex<ToTp>>(val * static_cast<FromTp>(slope) + static_cast<FromTp>(inter));
			}
			static std::complex<ToTp> Reverse(const std::complex<FromTp> val, const float slope, const float inter) {
				return static_cast<std::complex<ToTp>>((val - static_cast<FromTp>(inter)) / static_cast<FromTp>(slope));
			}
		};
		
		template<typename ByteIter, typename Iter>
		void fromLoop(const ByteIter fromBegin, Iter toBegin, Iter toEnd) {
			typedef const typename std::remove_reference<decltype(*fromBegin)>::type f_t;
			typedef typename std::remove_reference<decltype(*toBegin)>::type t_t;
			ByteIter from(fromBegin);
			Iter to(toBegin);
			while (to != toEnd) {
				f_t temp = *from;
				*to = Scale<f_t, t_t>::Forward(temp, scaling_slope, scaling_inter);
				to++;
				from++;
			}
		}
		
		template<typename ByteIter, typename Iter>
		void toLoop(const Iter fromBegin, const Iter fromEnd, ByteIter toBegin) {
			typedef const typename std::remove_reference<decltype(*fromBegin)>::type f_t;
			typedef typename std::remove_reference<decltype(*toBegin)>::type t_t;
			Iter from(fromBegin);
			ByteIter to(toBegin);
			while (from != fromEnd) {
				*to = Scale<f_t, t_t>::Reverse(*from, scaling_slope, scaling_inter);
				to++;
				from++;
			}
		}
		
		template<typename Iter>
		void fromBytes(const std::vector<char> &bytes, Iter begin, Iter end);
		template<typename Iter>
		void toBytes(const Iter begin, const Iter end, std::vector<char> &bytes);
		
		template<typename Iter>
		void readWriteVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size,
						     Iter &begin, Iter &end);
		
	#pragma mark Public Class Methods
	public:
		~Nifti();
		Nifti();                               //!< Default constructor. Initialises an empty header, size 1 in all dimensions.
		Nifti(const Nifti &other);             //!< Copy constructor. Copies all elements, and if the original is open then also opens new file handles.
		Nifti &operator=(const Nifti &other);  //!< Copy Assignment. Copies all elements except file handles, and marks destination as Closed.
		Nifti(Nifti &&other) noexcept;         //!< Move constructor. Copies all elements, including the file handles, and marks the original as Closed.
		Nifti &operator=(Nifti &&other);       //!< Move assignment. Copies all elements, including the file handles, and marks the original as Closed.
		Nifti(const int nx, const int ny, const int nz, const int nt,
			  const float dx, const float dy, const float dz, const float dt,
			  const DataType dtype = DataType::FLOAT32, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		Nifti(const ArrayXs &dim, const Eigen::ArrayXf &voxdim,
			  const DataType dtype = DataType::FLOAT32, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity()); //!< Constructs a header with the specified dimension and voxel sizes.
		Nifti(const Nifti &other, const size_t nt, const DataType dtype = DataType::FLOAT32);                        //!< Copies only basic geometry information from other, then sets the datatype and number of volumes. Does not copy scaling information etc.
		Nifti(const std::string &filename, const Mode &mode);
		
		void open(const std::string &filename, const Mode &mode); //!< Attempts to open a NIfTI file. Throws runtime_error or invalid_argument on failure.
		void close();                                             //!< Closes the file
		bool isOpen();                                            //!< Returns true if file is currently open for reading or writing.
		
		const std::string &basePath() const;
		std::string imagePath() const;
		std::string headerPath() const;
		
		const DataType &datatype() const;
		void setDatatype(const DataType dt);
		
		size_t rank() const;                                    //!< Get the rank (number of dimensions) of the image.
		size_t dim(const size_t d) const;                       //!< Get the size (voxel count) of a dimension. Valid dimensions are 1-7.
		ArrayXs dims() const;                                   //!< Get all dimension sizes.
		void setDim(const size_t d, const size_t n);            //!< Set the size (voxel count) of a dimension. Valid dimensions are 1-7.
		void setDims(const ArrayXs &newDims);                   //!< Set all dimension sizes.
		
		float voxDim(const size_t d) const;                     //!< Get the voxel size along dimension d. Valid dimensions are 1-7.
		Eigen::ArrayXf voxDims() const;                  //!< Get all voxel sizes.
		void setVoxDim(const size_t d, const float f);          //!< Set the voxel size along dimension d. Valid dimensions are 1-7.
		void setVoxDims(const Eigen::ArrayXf &newVoxDims);      //!< Set all voxel sizes.
		
		void setTransform(const Eigen::Affine3f &t, const XForm tc = XForm::ScannerAnatomy); //!< Set the qform and sform from a 4x4 general matrix. The qform will be set to closest matching linear XForm, the sform will be an exact copy.
		const Eigen::Affine3f &transform() const;            //!< Return the XForm with the highest priority.
		const Eigen::Affine3f &qform() const;                //!< Return just the qform.
		const Eigen::Affine3f &sform() const;                //!< Return just the sform.
		const XForm &qcode() const;               //!< Find out what transformation the qform represents.
		const XForm &scode() const;               //!< Find out what transformation the sform represents.
		bool matchesSpace(const Nifti &other) const;  //!< Check if voxel dimensions, data size and XForm match
		bool matchesVoxels(const Nifti &other) const; //!< Looser check if voxel dimensions and data size match
		
		template<typename T> void readVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size,
						                     typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end);
		template<typename T> void readVolumes(const size_t first, const size_t nvol,
						                      typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end);
		template<typename T> void writeVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size,
						                      typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end);
		template<typename T> void writeVolumes(const size_t vol, const size_t nvol,
						                       typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end);
		void addExtension(const int code, const std::vector<char> &data);
		void addExtension(const Extension &e);
		const std::list<Extension> &extensions() const;
		
		#pragma mark Information bits of the NIfTI header
		float scaling_slope;          //!< Slope of scaling between data on disk and in memory.
		float scaling_inter;          //!< Intercept of scaling between data on disk and in memory.
		float calibration_min;        //!< Suggested minimum for display.
		float calibration_max;        //!< Suggested maximum for display.
		
		int freq_dim ;                //!< Index of the frequency encode direction (1-3).
		int phase_dim;                //!< Index of the phase encode direction (1-3).
		int slice_dim;                //!< Index of the slice direction (1-3).
		
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
};

#include "Nifti-inl.h"

#endif // NIFTI_NIFTI
