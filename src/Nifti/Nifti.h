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

#include "ZipFile.h"
#include "Extension.h"
#include "Header.h"

namespace Nifti {

enum class Mode : char { Closed = 0, Read = 'r', ReadHeader = 'h', Write = 'w' };

class File {
	private:
		std::string m_basepath;            //!< Path to file without extension.
		ZipFile m_file;                    //!< Actual disk file
		Header m_header;                   //!< The header object for this file
		std::list<Extension> m_extensions; //!< List of Nifti extensions
		Mode m_mode;                       //!< Whether the file is closed or open for reading/writing.
		bool m_nii;                        //!< Is this a .nii file?
		bool m_gz;                         //!< Is this a gzipped file?
		Version m_nifti_version;           //!< Is this a nifti1 or a nifti2 file?
		int m_swap;                        //!< True if byte order on disk is different to CPU.
		
		void setPaths(const std::string &path); //!< Works out the basepath and file extensions.
		void readHeader();                 //!< Attempts to read a header structure from the currently open file.
		void readExtensions();             //!< Attempts to read any extensions
		void writeHeader();                //!< Attempts to write a header structure to the currently open file.
		void writeExtensions();            //!< Attempts to write extensions
		int totalExtensionSize();          //!< Counts the total number of bytes for all extensions.

		void seekToVoxel(const IndexArray &target);
		void readBytes(std::vector<char> &data);
		void writeBytes(const std::vector<char> &data);

		template<typename FromTp, typename ToTp> class Scale; //!< Templated class to scale and cast between types. Definitions in Nifti-inl.h
		template<typename Iter> void readWriteVoxels(const IndexArray &start, const IndexArray &size, Iter &begin, Iter &end); //!< Core IO routine. Can read/write to any storage that supports iterators/pointers

	#pragma mark Public Class Methods
	public:
		~File();
		File();                               //!< Default constructor. Initialises an empty header, size 1 in all dimensions.
		File(const File &other);             //!< Copy constructor. Copies all elements, and if the original is open then also opens new file handles.
		File(File &&other) noexcept;         //!< Move constructor. Copies all elements, including the file handles, and marks the original as Closed.
		File(const Header &hdr, const std::string &filename); //!< Creates a new File and opens it for writing with the specified header and filename.
		File(const File &other, const size_t nt, const DataType dtype = DataType::FLOAT32);                        //!< Copies only basic geometry information from other, then sets the datatype and number of volumes. Does not copy scaling information etc.
		File(const std::string &filename, const Mode &mode);
		
		void open(const std::string &filename, const Mode &mode); //!< Attempts to open a NIfTI file. Throws runtime_error or invalid_argument on failure.
		void close();                                             //!< Closes the file
		bool isOpen();                                            //!< Returns true if file is currently open for reading or writing.
		
		const std::string &basePath() const;
		std::string imagePath() const;
		std::string headerPath() const;
		
		const Header &header() const;                           //!< Get the header information for this image.
		size_t rank() const;                                    //!< Get the rank (number of dimensions) of the image.
		size_t dim(const size_t d) const;                       //!< Get the size (voxel count) of a dimension. Valid dimensions are 1-7.
		IndexArray dims() const;                                //!< Get all dimension sizes.

		template<typename IterTp> void readVoxels(IterTp begin, IterTp end, const Eigen::Ref<IndexArray> &start, const Eigen::Ref<IndexArray> &size);
		template<typename IterTp> void readVolumes(IterTp begin, IterTp end, const size_t first = 0, const size_t nvol = 0);
		template<typename IterTp> void readAll(IterTp begin, IterTp end);
		template<typename IterTp> void writeVoxels(IterTp begin, IterTp end, const Eigen::Ref<IndexArray> &start, const Eigen::Ref<IndexArray> &size);
		template<typename IterTp> void writeVolumes(IterTp begin, IterTp end, const size_t first = 0, const size_t nvol = 0);
		template<typename IterTp> void writeAll(IterTp begin, IterTp end);
		void addExtension(const int code, const std::vector<char> &data);
		void addExtension(const Extension &e);
		const std::list<Extension> &extensions() const;
		void setExtensions(const std::list<Extension> &es);
};

#include "Nifti-inl.h"

} // End namespace Nifti

#endif // NIFTI_NIFTI
