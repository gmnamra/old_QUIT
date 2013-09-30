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
  *   The number of elements is implicit in the size of the input byte array as the NIfTI stores the element size.
  *   The output array can be larger than required, in which case offset can be non-zero and data will be written
  *   from that point onwards.
  *   
  *   @param T Desired datatype. Valid (scalar) conversions must exist.
  *   @param bytes  - Raw bytes read from NIfTI file.
  *   @param data   - Storage for converted data. Must contain enough space.
  *   @param offset - Location in data to start writing.
  */
template<typename T> void Nifti::convertFromBytes(const std::vector<char> &bytes,
                                                  const typename std::vector<T>::iterator begin,
												  const typename std::vector<T>::iterator end) {
	size_t nEl = bytes.size() / m_typeinfo.size;
	assert(nEl == std::distance(begin, end));
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	for (size_t i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:       *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const int8_t *>(bytes.data())[i]); break;
			case DataType::INT16:      *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const int16_t *>(bytes.data())[i]); break;
			case DataType::INT32:      *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const int32_t *>(bytes.data())[i]); break;
			case DataType::INT64:      *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const int64_t *>(bytes.data())[i]); break;
			case DataType::UINT8:      *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const uint8_t *>(bytes.data())[i]); break;
			case DataType::UINT16:     *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const uint16_t *>(bytes.data())[i]); break;
			case DataType::UINT32:     *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const uint32_t *>(bytes.data())[i]); break;
			case DataType::UINT64:     *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const uint64_t *>(bytes.data())[i]); break;
			case DataType::FLOAT32:    *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]); break;
			case DataType::FLOAT64:    *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]); break;
			case DataType::FLOAT128:   *el = sc_in + sc_sl * static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]); break;
			case DataType::COMPLEX64:  *el = sc_in + sc_sl * static_cast<T>(abs(reinterpret_cast<const std::complex<float> *>(bytes.data())[i])); break;
			case DataType::COMPLEX128: *el = sc_in + sc_sl * static_cast<T>(abs(reinterpret_cast<const std::complex<double> *>(bytes.data())[i])); break;
			case DataType::COMPLEX256: *el = sc_in + sc_sl * static_cast<T>(abs(reinterpret_cast<const std::complex<long double> *>(bytes.data())[i])); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
		el++;
	}
}

template<typename T> void Nifti::convertFromBytes(const std::vector<char> &bytes,
                                                  const typename std::vector<std::complex<T>>::iterator begin,
												  const typename std::vector<std::complex<T>>::iterator end) {
	size_t nEl = bytes.size() / m_typeinfo.size;
	assert(nEl == std::distance(begin, end));
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	for (size_t i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:       *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const int8_t *>(bytes.data())[i]), 0.); break;
			case DataType::INT16:      *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const int16_t *>(bytes.data())[i]), 0.); break;
			case DataType::INT32:      *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const int32_t *>(bytes.data())[i]), 0.); break;
			case DataType::INT64:      *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const int64_t *>(bytes.data())[i]), 0.); break;
			case DataType::UINT8:      *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const uint8_t *>(bytes.data())[i]), 0.); break;
			case DataType::UINT16:     *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const uint16_t *>(bytes.data())[i]), 0.); break;
			case DataType::UINT32:     *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const uint32_t *>(bytes.data())[i]), 0.); break;
			case DataType::UINT64:     *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const uint64_t *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT32:    *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const float *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT64:    *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const double *>(bytes.data())[i]), 0.); break;
			case DataType::FLOAT128:   *el = sc_in + sc_sl * std::complex<T>(static_cast<T>(reinterpret_cast<const long double *>(bytes.data())[i]), 0.); break;
			case DataType::COMPLEX64:  *el = sc_in + sc_sl * static_cast<std::complex<T> >(reinterpret_cast<const std::complex<float> *>(bytes.data())[i]); break;
			case DataType::COMPLEX128: *el = sc_in + sc_sl * static_cast<std::complex<T> >(reinterpret_cast<const std::complex<double> *>(bytes.data())[i]); break;
			case DataType::COMPLEX256: *el = sc_in + sc_sl * static_cast<std::complex<T> >(reinterpret_cast<const std::complex<long double> *>(bytes.data())[i]); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
		el++;
	}
}

/**
  *   Internal function to convert from arbitrary data to the internal NIfTI dataype.
  *
  *   Converts a from the templated datatype to a sequence of bytes.
  *   If complex data is passed in but the file datatype is not complex, then the absolute magnitude is taken.
  *   If non-complex data is passed in but the file datatype is complex, then the imaginary part will be 0.
  *   
  *   The number of elements is implicit in the size of the input byte array as the NIfTI stores the element size.
  *   The output array can be larger than required, in which case offset can be non-zero and data will be written
  *   from that point onwards.
  *   
  *   @param T Desired datatype. Valid (scalar) conversions must exist.
  *   @param bytes  - Raw bytes read from NIfTI file.
  *   @param data   - Storage for converted data. Must contain enough space.
  *   @param offset - Location in data to start writing.
  */
template<typename T> void Nifti::convertToBytes(const typename std::vector<T>::iterator begin,
                                                const typename std::vector<T>::iterator end,
												std::vector<char> &bytes) {
	size_t nEl = std::distance(begin, end);
	assert(nEl == bytes.size() / m_typeinfo.size);
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	for (size_t i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<int8_t>(*el / sc_sl - sc_in); break;
			case DataType::INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<int16_t>(*el / sc_sl - sc_in); break;
			case DataType::INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int32_t>(*el / sc_sl - sc_in); break;
			case DataType::INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<int64_t>(*el / sc_sl - sc_in); break;
			case DataType::UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<uint8_t>(*el / sc_sl - sc_in); break;
			case DataType::UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<uint16_t>(*el / sc_sl - sc_in); break;
			case DataType::UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<uint32_t>(*el / sc_sl - sc_in); break;
			case DataType::UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<uint64_t>(*el / sc_sl - sc_in); break;
			case DataType::FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(*el / sc_sl - sc_in); break;
			case DataType::FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(*el / sc_sl - sc_in); break;
			case DataType::FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(*el / sc_sl - sc_in); break;
			case DataType::COMPLEX64:  reinterpret_cast<std::complex<float> *>(bytes.data())[i] = std::complex<float>(static_cast<float>(*el / sc_sl - sc_in)); break;
			case DataType::COMPLEX128: reinterpret_cast<std::complex<double> *>(bytes.data())[i] = std::complex<double>(static_cast<double>(*el / sc_sl - sc_in)); break;
			case DataType::COMPLEX256: reinterpret_cast<std::complex<long double> *>(bytes.data())[i] = std::complex<long double>(static_cast<long double>(*el / sc_sl - sc_in)); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
		el++;
	}
}

template<typename T> void Nifti::convertToBytes(const typename std::vector<std::complex<T>>::iterator begin,
												const typename std::vector<std::complex<T>>::iterator end,
												std::vector<char> &bytes) {
	size_t nEl = std::distance(begin, end);
	assert(nEl == bytes.size() / m_typeinfo.size);
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	for (int i = 0; i < nEl; i++) {
		switch (m_typeinfo.type) {
			case DataType::INT8:              reinterpret_cast<char *>(bytes.data())[i] = static_cast<int8_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::INT16:            reinterpret_cast<short *>(bytes.data())[i] = static_cast<int16_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::INT32:              reinterpret_cast<int *>(bytes.data())[i] = static_cast<int32_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::INT64:             reinterpret_cast<long *>(bytes.data())[i] = static_cast<int64_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::UINT8:    reinterpret_cast<unsigned char *>(bytes.data())[i] = static_cast<uint8_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::UINT16:  reinterpret_cast<unsigned short *>(bytes.data())[i] = static_cast<uint16_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::UINT32:    reinterpret_cast<unsigned int *>(bytes.data())[i] = static_cast<uint32_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::UINT64:   reinterpret_cast<unsigned long *>(bytes.data())[i] = static_cast<uint64_t>(abs(*el / sc_sl - sc_in)); break;
			case DataType::FLOAT32:          reinterpret_cast<float *>(bytes.data())[i] = static_cast<float>(abs(*el / sc_sl - sc_in)); break;
			case DataType::FLOAT64:         reinterpret_cast<double *>(bytes.data())[i] = static_cast<double>(abs(*el / sc_sl - sc_in)); break;
			case DataType::FLOAT128:   reinterpret_cast<long double *>(bytes.data())[i] = static_cast<long double>(abs(*el / sc_sl - sc_in)); break;
			case DataType::COMPLEX64:  reinterpret_cast<std::complex<float> *>(bytes.data())[i] = static_cast<std::complex<float> >(*el / sc_sl - sc_in); break;
			case DataType::COMPLEX128: reinterpret_cast<std::complex<double> *>(bytes.data())[i] = static_cast<std::complex<double> >(*el / sc_sl - sc_in); break;
			case DataType::COMPLEX256: reinterpret_cast<std::complex<long double> *>(bytes.data())[i] = static_cast<std::complex<long double> >(*el / sc_sl - sc_in); break;
			case DataType::RGB24: case DataType::RGBA32:
				throw(std::runtime_error("Unsupported datatype.")); break;
		}
		el++;
	}
}

/*
 *   Core IO routine. Depending on the mode of the file, reads/writes a
 *   contiguous subregion of the region and converts to/from the desired datatype.
 *
 *   start and size can be short of the full number of dimensions, e.g. if you
 *   only want a part of the first volume in a multi-volume image, then they can
 *   just be 3D. The only limitation is that they must not have more dimensions
 *   than the image.
 *
 *   @param start The voxel indices of the first desired voxel.
 *   @param size  The size of the desired subregion.
 *   @param data  Storage for the data to read/write. Must be sufficiently large.
 */
template<typename T> void Nifti::readWriteVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size, std::vector<T> &data) {
	if (start.rows() != size.rows()) throw(std::out_of_range("Start and size must have same dimension in image: " + imagePath()));
	if (start.rows() > m_dim.rows()) throw(std::out_of_range("Too many read dimensions specified in image: " + imagePath()));
	if ((size == 0).any()) throw(std::out_of_range("Requested a zero-length read in 1 or more dimensions of image: " + imagePath()));
	if (((start + size) > m_dim.head(start.rows())).any()) throw(std::out_of_range("Requested read was larger than image dimensions: " + imagePath()));
	if (size.prod() < data.size()) throw(std::out_of_range("Allocated memory is insufficient."));
	
	size_t firstDim = 0; // We can always read first dimension in one go
	size_t blockSize = size(firstDim);
	auto dataIt = data.begin();
	while ((size(firstDim) == m_dim(firstDim)) && (firstDim < size.rows() - 1)) {
		firstDim++;
		blockSize *= size(firstDim);
	}
	//cout << "firstDim " << firstDim << " blockSize " << blockSize << endl;
	std::vector<char> block(blockSize * m_typeinfo.size);
	ArrayXs blockStart = start;
	std::function<void (const size_t dim)> dimLoop;
	dimLoop = [&] (const size_t dim) {
		if (dim == firstDim) {
			//cout << "target " << target.transpose() << endl;
			seekToVoxel(blockStart);
			if (m_mode == Nifti::Mode::Read) {
				readBytes(block);
				convertFromBytes<T>(block, dataIt, dataIt + blockSize);
			} else {
				convertToBytes<T>(dataIt, dataIt + blockSize, block);
				writeBytes(block);
			}
			dataIt += blockSize;
		} else {
			//cout << "dim " << dim << " start " << start(dim) << " end " << start(dim) + size(dim) << endl;
			for (size_t v = start(dim); v < start(dim) + size(dim); v++) {
				blockStart(dim) = v;
				dimLoop(dim - 1);
			}
		}
	};
	dimLoop(start.rows() - 1);
}

template<typename T> void Nifti::readVolumes(const size_t first, const size_t nvol, std::vector<T> &data) {
	if (!((m_mode == Mode::Read) || (m_mode == Mode::ReadSkipExt)))
		throw(std::runtime_error("File must be opened for reading: " + basePath()));
	if (data.size() != (voxelsPerVolume() * nvol))
		throw(std::runtime_error("Insufficient storage allocated for read: " + basePath()));
	
	ArrayXs start, size;
	start << 0, 0, first;
	size << dim(1), dim(2), first + nvol;
	readWriteVoxels(start, size, data);
}

template<typename T> void Nifti::writeVolumes(const size_t first, const size_t nvol, const std::vector<T> &data) {
	if (!((m_mode == Mode::Write) || (m_mode == Mode::WriteSkipExt)))
		throw(std::runtime_error("File must be opened for writing: " + basePath()));
	if (data.size() != (voxelsPerVolume() * nvol))
		throw(std::runtime_error("Insufficient data for write: " + basePath()));
	
	ArrayXs start, size;
	start << 0, 0, first;
	size << dim(1), dim(2), first + nvol;
	readWriteVoxels(start, size, data);
}

#endif // NIFTI_NIFTI_INL
