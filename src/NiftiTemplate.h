/** \file Nifti-inl.h
 \brief Definitions of templated IO code for Nifti::File
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 - This file should NOT be included directly in projects. It is just to 
 - separate the interface from the definition a bit better.
 */

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
template<typename T> void File::convertFromBytes(const vector<char> &bytes, const size_t nEl, vector<T> &data) {
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

/**
  *   Internal function to convert the internal NIfTI data to the desired complex dataype.
  *
  *   Converts a sequence of bytes from a NIfTI image to the templated complex datatype.
  *   If the file contains complex data but a real datatype is requested then the magnitude is returned.
  *   If the file contains real data but a complex datatype is requested then the imaginary part will be 0.
  *   
  *   @param T Desired datatype. Valid (scalar) conversions must exist.
  *   @param bytes std::vector of byte data
  *   @param nEl Number of elements (not bytes) expected in the data
  *   @param data std::vector& to converted data in. Will be resized to ensure enough space.
  */
template<typename T> void File::convertFromBytes(const vector<complex<T>> &bytes, const size_t nEl, vector<T> &data) {
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
template<typename T> vector<char> File::convertToBytes(const vector<T> &data) {
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

template<typename T> vector<char> File::convertToBytes(const vector<complex<T>> &data) {
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

template<typename T> void File::readVolume(const size_t &vol, vector<T> &buffer) {
	size_t bytesPerVolume = voxelsPerVolume() * bytesPerVoxel();
	vector<char> raw(bytesPerVolume);
	readBytes(vol * bytesPerVolume, bytesPerVolume, raw.data());
	convertFromBytes(raw, voxelsPerVolume(), buffer);
}

template<typename T> vector<T> File::readVolume(const size_t &vol) {
	vector<T> buffer;
	readVolume(vol, buffer);
	return buffer;
}

template<typename T> void File::readAllVolumes(vector<T> &buffer) {
	vector<char> raw(voxelsTotal() * bytesPerVoxel());
	readBytes(0, voxelsTotal() * bytesPerVoxel(), raw.data());
	return convertFromBytes<T>(raw, voxelsTotal(), buffer);
}

template<typename T> vector<T> File::readAllVolumes() {
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
template<typename T> void File::readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
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

template<typename T> void File::readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
										const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et) {
	vector<T> buffer;
	readSubvolume(sx, sy, sz, st, ex, ey, ez, et, buffer);
	return buffer;
}

template<typename T> void File::writeVolume(const size_t vol, const vector<T> &data) {
	if (data.size() != voxelsPerVolume()) {
		throw(invalid_argument("Insufficient data to write volume for file: " + imagePath()));
	}
	vector<char> converted = convertToBytes(data);
	writeBytes(vol * voxelsPerVolume() * bytesPerVoxel(), converted.size(), converted.data());
}

template<typename T> void File::writeAllVolumes(const vector<T> &data) {
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
template<typename T> void File::writeSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
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
