/** \file Nifti-inl.h
 \brief Definitions of templated IO code for Nifti
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 - This file should NOT be included directly in projects. It is just to 
 - separate the interface from the definition a bit better.
 */

#ifndef NIFTI_NIFTI_INL
#define NIFTI_NIFTI_INL

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
template<typename T> void Nifti::convertFromBytes(const std::vector<char> &bytes, const size_t nEl, std::vector<T> &data) {
	assert(nEl == (bytes.size() / m_typeinfo.size));
	data.resize(nEl);
	for (size_t i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:      data[i] = static_cast<T>(reinterpret_cast<const char *>(bytes.data())[i]); break;
			case DataType::UINT8:     data[i] = static_cast<T>(reinterpret_cast<const unsigned char *>(bytes.data())[i]); break;
			case DataType::INT16:     data[i] = static_cast<T>(reinterpret_cast<const short *>(bytes.data())[i]); break;
			case DataType::UINT16:    data[i] = static_cast<T>(reinterpret_cast<const unsigned short *>(bytes.data())[i]); break;
			case DataType::INT32:     data[i] = static_cast<T>(reinterpret_cast<const int *>(bytes.data())[i]); break;
			case DataType::UINT32:    data[i] = static_cast<T>(reinterpret_cast<const unsigned int *>(bytes.data())[i]); break;
			case DataType::FLOAT32:   data[i] = static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]); break;
			case DataType::FLOAT64:   data[i] = static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]); break;
			case DataType::INT64:     data[i] = static_cast<T>(reinterpret_cast<const long *>(bytes.data())[i]); break;
			case DataType::UINT64:    data[i] = static_cast<T>(reinterpret_cast<const unsigned long *>(bytes.data())[i]); break;
			case DataType::FLOAT128:  data[i] = static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]); break;
			// NOTE: C++11 specifies that C++ 'complex<type>' and C 'type complex'
			// should be interchangeable even at pointer level
			case DataType::COMPLEX64:  data[i] = static_cast<T>(abs(reinterpret_cast<const std::complex<float> *>(bytes.data())[i])); break;
			case DataType::COMPLEX128: data[i] = static_cast<T>(abs(reinterpret_cast<const std::complex<double> *>(bytes.data())[i])); break;
			case DataType::COMPLEX256: data[i] = static_cast<T>(abs(reinterpret_cast<const std::complex<long double> *>(bytes.data())[i])); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
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
template<typename T> void Nifti::convertFromBytes(const std::vector<std::complex<T>> &bytes, const size_t nEl, std::vector<T> &data) {
	assert(nEl == (bytes.size() / m_typeinfo.size));
	data.resize(nEl);
	for (size_t i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:      data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const char *>(bytes.data())[i]), 0.); break;
			case DataType::UINT8:     data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const unsigned char *>(bytes.data())[i]), 0.); break;
			case DataType::INT16:     data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const short *>(bytes.data())[i]), 0.); break;
			case DataType::UINT16:    data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const unsigned short *>(bytes.data())[i]), 0.); break;
			case DataType::INT32:     data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const int *>(bytes.data())[i]), 0.); break;
			case DataType::UINT32:    data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const unsigned int *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT32:   data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT64:   data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]), 0.); break;
			case DataType::INT64:     data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const long *>(bytes.data())[i]), 0.); break;
			case DataType::UINT64:    data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const unsigned long *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT128:  data[i] = std::complex<T>(static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]), 0.); break;
			case DataType::COMPLEX64:  data[i] = static_cast<std::complex<T> >(reinterpret_cast<const std::complex<float> *>(bytes.data())[i]); break;
			case DataType::COMPLEX128: data[i] = static_cast<std::complex<T> >(reinterpret_cast<const std::complex<double> *>(bytes.data())[i]); break;
			case DataType::COMPLEX256: data[i] = static_cast<std::complex<T> >(reinterpret_cast<const std::complex<long double> *>(bytes.data())[i]); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
	}
}

/**
  *   Internal function to convert data to the internal NIfTI datatype.
  *
  *   Converts from the templated type to a sequence of bytes suitable for writing to a NIfTI image.
  *   @param data std::std::vector of the data.
  *   @param T Desired datatype. Valid (scalar) conversions must exist.
  *   @return std::std::vector of the data stored in a byte array.
  */
template<typename T> std::vector<char> Nifti::convertToBytes(const std::vector<T> &data) {
	std::vector<char> bytes(data.size() * m_typeinfo.size);
	for (size_t i = 0; i < data.size(); i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<char>(data[i]); break;
			case DataType::UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<unsigned char>(data[i]); break;
			case DataType::INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<short>(data[i]); break;
			case DataType::UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<unsigned short>(data[i]); break;
			case DataType::INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int>(data[i]); break;
			case DataType::UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<unsigned int>(data[i]); break;
			case DataType::FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(data[i]); break;
			case DataType::FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(data[i]); break;
			case DataType::INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<long>(data[i]); break;
			case DataType::UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<unsigned long>(data[i]); break;
			case DataType::FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(data[i]); break;
			case DataType::COMPLEX64:  reinterpret_cast<std::complex<float> *>(bytes.data())[i] = std::complex<float>(static_cast<float>(data[i])); break;
			case DataType::COMPLEX128: reinterpret_cast<std::complex<double> *>(bytes.data())[i] = std::complex<double>(static_cast<double>(data[i])); break;
			case DataType::COMPLEX256: reinterpret_cast<std::complex<long double> *>(bytes.data())[i] = std::complex<long double>(static_cast<long double>(data[i])); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
	}
	return bytes;
}

template<typename T> std::vector<char> Nifti::convertToBytes(const std::vector<std::complex<T>> &data) {
	std::vector<char> bytes(data.size() * m_typeinfo.size);
	for (int i = 0; i < data.size(); i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<char>(abs(data[i])); break;
			case DataType::UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<unsigned char>(abs(data[i])); break;
			case DataType::INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<short>(abs(data[i])); break;
			case DataType::UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<unsigned short>(abs(data[i])); break;
			case DataType::INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int>(abs(data[i])); break;
			case DataType::UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<unsigned int>(abs(data[i])); break;
			case DataType::FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(abs(data[i])); break;
			case DataType::FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(abs(data[i])); break;
			case DataType::INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<long>(abs(data[i])); break;
			case DataType::UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<unsigned long>(abs(data[i])); break;
			case DataType::FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(abs(data[i])); break;
			case DataType::COMPLEX64:  reinterpret_cast<std::complex<float> *>(bytes.data())[i] = static_cast<std::complex<float> >(data[i]); break;
			case DataType::COMPLEX128: reinterpret_cast<std::complex<double> *>(bytes.data())[i] = static_cast<std::complex<double> >(data[i]); break;
			case DataType::COMPLEX256: reinterpret_cast<std::complex<long double> *>(bytes.data())[i] = static_cast<std::complex<long double> >(data[i]); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
	}
	return bytes;
}

template<typename T> void Nifti::readVolume(const size_t &vol, std::vector<T> &buffer) {
	size_t bytesPerVolume = voxelsPerVolume() * m_typeinfo.size;
	std::vector<char> raw(bytesPerVolume);
	readBytes(vol * bytesPerVolume, bytesPerVolume, raw.data());
	convertFromBytes(raw, voxelsPerVolume(), buffer);
}

template<typename T> std::vector<T> Nifti::readVolume(const size_t &vol) {
	std::vector<T> buffer;
	readVolume(vol, buffer);
	return buffer;
}

template<typename T> void Nifti::readAllVolumes(std::vector<T> &buffer) {
	std::vector<char> raw(voxelsTotal() * m_typeinfo.size);
	readBytes(0, voxelsTotal() * m_typeinfo.size, raw.data());
	return convertFromBytes<T>(raw, voxelsTotal(), buffer);
}

template<typename T> std::vector<T> Nifti::readAllVolumes() {
	std::vector<T> buffer;
	readAllVolumes(buffer);
	return buffer;
}

/**
  *   Reads a subvolume
  *
  *   Converts to the templated type from the NIfTI datatype.
  *   @param sx, sy, sz, st - Start voxels for each dimension.
  *   @param ex, ey, ez, et - End voxels for each dimension. Values are clamped to dimension sizes.
  *   @param buffer std::std::vector to store the data.
  *
  */
template<typename T> void Nifti::readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
										      const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
										      std::vector<T> &buffer) {
	size_t lx, ly, lz, lt, total, toRead;
	lx = ((ex > dim(1)) ? dim(1) : ex) - sx;
	ly = ((ey > dim(2)) ? dim(2) : ey) - sy;
	lz = ((ez > dim(3)) ? dim(3) : ez) - sz;
	lt = ((et > dim(4)) ? dim(4) : et) - st;
	total = lx * ly * lz * lt;
	
	if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
		throw(std::out_of_range("Invalid subvolume read dimensions: " + imagePath()));
	}
	
	// Collapse successive full dimensions into a single compressed read
	toRead = lx * m_typeinfo.size;
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
				
	std::vector<char> raw(total * m_typeinfo.size);
	char *nextRead = raw.data();
	for (size_t t = st; t < st+lt; t++) {
		size_t tOff = t * voxelsPerVolume();
		for (size_t z = sz; z < sz+lz; z++) {
			size_t zOff = z * voxelsPerSlice();
			for (size_t y = sy; y < sy+ly; y++) {
				size_t yOff = y * dim(1);
				if (readBytes((tOff + zOff + yOff) * m_typeinfo.size, toRead, nextRead))
					nextRead += toRead;
			}
		}
	}
	convertFromBytes<T>(raw, total, buffer);
}

template<typename T> void Nifti::readSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
										const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et) {
	std::vector<T> buffer;
	readSubvolume(sx, sy, sz, st, ex, ey, ez, et, buffer);
	return buffer;
}

template<typename T> void Nifti::writeVolume(const size_t vol, const std::vector<T> &data) {
	if (data.size() != voxelsPerVolume()) {
		throw(std::invalid_argument("Insufficient data to write volume for file: " + imagePath()));
	}
	std::vector<char> converted = convertToBytes(data);
	writeBytes(vol * voxelsPerVolume() * m_typeinfo.size, converted.size(), converted.data());
}

template<typename T> void Nifti::writeAllVolumes(const std::vector<T> &data) {
	if (data.size() != voxelsTotal()) {
		throw(std::invalid_argument("Insufficient data to write all volumes for file: " + imagePath()));
	}
	std::vector<char> converted = convertToBytes(data);
	writeBytes(0, converted.size(), converted.data());
}

/**
  *   Writes a subvolume
  *
  *   Converts from the templated type to a sequence of bytes suitable for writing to a NIfTI image.
  *   @param sx, sy, sz, st - Start voxels for each dimension.
  *   @param ex, ey, ez, et - End voxels for each dimension. Values are clamped to dimension sizes.
  *   @param data std::std::vector of the data.
  *
  */
template<typename T> void Nifti::writeSubvolume(const size_t &sx, const size_t &sy, const size_t &sz, const size_t &st,
										       const size_t &ex, const size_t &ey, const size_t &ez, const size_t &et,
										       const std::vector<T> &data) {
	size_t lx, ly, lz, lt, total, toWrite;
	lx = ((ex > dim(1)) ? dim(1) : ex) - sx;
	ly = ((ey > dim(2)) ? dim(2) : ey) - sy;
	lz = ((ez > dim(3)) ? dim(3) : ez) - sz;
	lt = ((et > dim(4)) ? dim(4) : et) - st;
	total = lx * ly * lz * lt;
	
	if (lx < 1 || ly < 1 || lz < 1 || lt < 1) { // There is nothing to write
		throw(std::out_of_range("Invalid subvolume read dimensions: " + imagePath()));
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
		throw(std::invalid_argument("Insufficient data to write subvolume for file: " + imagePath()));
	}
	toWrite *= m_typeinfo.size;
	
	std::vector<char> raw = convertToBytes(data);
	char *nextWrite = raw.data();
	for (size_t t = st; t < st+lt; t++) {
		size_t tOff = t * voxelsPerVolume();
		for (size_t z = sz; z < sz+lz; z++) {
			size_t zOff = z * voxelsPerSlice();
			for (size_t y = sy; y < sy+ly; y++) {
				size_t yOff = y * dim(1);
				writeBytes((tOff + zOff + yOff) * m_typeinfo.size, toWrite, nextWrite);
				nextWrite += toWrite;
			}
		}
	}
}

#endif // NIFTI_NIFTI_INL
