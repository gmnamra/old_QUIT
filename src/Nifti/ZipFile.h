/** \file ZipFile.h
 \brief Declaration for a wrapper around gzipped or unzipped files
 - Written by Tobias Wood, IoP KCL
 - Based on znzlib (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */


#ifndef NIFTI_ZIPFILE
#define NIFTI_ZIPFILE

#include <cstdio>
#include <string>
#include <limits>

using std::string;
using std::numeric_limits;

#include <zlib.h>

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

#endif // NIFTI_ZIPFILE
