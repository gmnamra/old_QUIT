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
template<typename T, typename Iter>
void Nifti::Converter<T, Iter>::FromBytes(const std::vector<char> &bytes, const Nifti::DataTypeInfo &tInfo,
                                    const float sc_sl, const float sc_in,
                                    const Iter &begin, const Iter &end) {
	typename Iter::difference_type nEl = bytes.size() / tInfo.size;
	assert(nEl == std::distance(begin, end));
	switch (tInfo.type) {
		case DataType::INT8:       { ConverterLoop<int8_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT16:      { ConverterLoop<int16_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT32:      { ConverterLoop<int32_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT64:      { ConverterLoop<int64_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT8:      { ConverterLoop<uint8_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT16:     { ConverterLoop<uint16_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT32:     { ConverterLoop<uint32_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT64:     { ConverterLoop<uint64_t, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT32:    { ConverterLoop<float, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT64:    { ConverterLoop<double, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT128:   { ConverterLoop<long double, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX64:  { ConverterLoop<std::complex<float>, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX128: { ConverterLoop<std::complex<double>, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX256: { ConverterLoop<std::complex<long double>, T, Iter>::From(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
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
template<typename T, typename Iter>
void Nifti::Converter<T, Iter>::ToBytes(std::vector<char> &bytes, const Nifti::DataTypeInfo &tInfo,
							      const float sc_sl, const float sc_in,
							      const Iter &begin, const Iter &end) {
	size_t nEl = std::distance(begin, end);
	assert(nEl == bytes.size() / tInfo.size);
	switch (tInfo.type) {
		case DataType::INT8:       { ConverterLoop<T, int8_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT16:      { ConverterLoop<T, int16_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT32:      { ConverterLoop<T, int32_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::INT64:      { ConverterLoop<T, int64_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT8:      { ConverterLoop<T, uint8_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT16:     { ConverterLoop<T, uint16_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT32:     { ConverterLoop<T, uint32_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::UINT64:     { ConverterLoop<T, uint64_t, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT32:    { ConverterLoop<T, float, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT64:    { ConverterLoop<T, double, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::FLOAT128:   { ConverterLoop<T, long double, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX64:  { ConverterLoop<T, std::complex<float>, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX128: { ConverterLoop<T, std::complex<double>, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::COMPLEX256: { ConverterLoop<T, std::complex<long double>, Iter>::To(bytes, begin, end, sc_sl, sc_in); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
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
 *   An entry of 0 in size indicates that we want to read the whole dimension.
 *
 *   @param start The voxel indices of the first desired voxel.
 *   @param size  The size of the desired subregion.
 *   @param data  Storage for the data to read/write. Must be sufficiently large.
 */
template<typename T, typename Iter>
void Nifti::readWriteVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &inSize,
                            const Iter &begin, const Iter &end) {
	ArrayXs size = inSize;
	for (ArrayXs::Index i = 0; i < size.rows(); i++)
		if (size(i) == 0) size(i) = m_dim(i);
	
	if (start.rows() != size.rows()) throw(std::out_of_range("Start and size must have same dimension in image: " + imagePath()));
	if (start.rows() > m_dim.rows()) throw(std::out_of_range("Too many read/write dimensions specified in image: " + imagePath()));
	if (((start + size) > m_dim.head(start.rows())).any()) throw(std::out_of_range("Read/write past image dimensions requested: " + imagePath()));
	if (size.prod() < static_cast<unsigned long>(end - begin)) throw(std::out_of_range("Storage size does not match requested read/write size in image: " + imagePath()));
	
	ArrayXs::Index firstDim = 0; // We can always read first dimension in one go
	ArrayXs::Index blockSize = size(firstDim);
	auto dataIt = begin;
	while ((size(firstDim) == m_dim(firstDim)) && (firstDim < size.rows() - 1)) {
		firstDim++;
		blockSize *= size(firstDim);
	}
	std::vector<char> block(blockSize * m_typeinfo.size);
	ArrayXs blockStart = start;
	std::function<void (const size_t dim)> dimLoop;
	dimLoop = [&] (const ArrayXs::Index dim) {
		if (dim == firstDim) {
			seekToVoxel(blockStart);
			if (m_mode == Nifti::Mode::Read) {
				readBytes(block);
				Converter<T, Iter>::FromBytes(block, m_typeinfo, scaling_slope, scaling_inter, dataIt, dataIt + blockSize);
			} else {
				Converter<T, Iter>::ToBytes(block, m_typeinfo, scaling_slope, scaling_inter, dataIt, dataIt + blockSize);
				writeBytes(block);
			}
			dataIt += blockSize;
		} else {
			for (size_t v = start(dim); v < start(dim) + size(dim); v++) {
				blockStart(dim) = v;
				dimLoop(dim - 1);
			}
		}
	};
	dimLoop(start.rows() - 1);
}

template<typename T>
void Nifti::readVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size,
					   typename std::vector<T>::iterator begin,
					   typename std::vector<T>::iterator end) {
	if (!(m_mode == Mode::Read))
		throw(std::runtime_error("File must be opened for reading: " + basePath()));
	readWriteVoxels<T>(start, size, begin, end);
}

template<typename T>
void Nifti::readVolumes(const size_t first, const size_t nvol,
                        typename std::vector<T>::iterator begin,
						typename std::vector<T>::iterator end) {
	if (!(m_mode == Mode::Read))
		throw(std::runtime_error("File must be opened for reading: " + basePath()));
	Eigen::Array<size_t, 4, 1> start{0, 0, 0, first};
	Eigen::Array<size_t, 4, 1> size{dim(1), dim(2), dim(3), nvol};
	readWriteVoxels<T>(start, size, begin, end);
}

template<typename T>
void Nifti::writeVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size,
                        typename std::vector<T>::iterator begin,
						typename std::vector<T>::iterator end) {
	if (!(m_mode == Mode::Write))
		throw(std::runtime_error("File must be opened for writing: " + basePath()));
	readWriteVoxels<T>(start, size, begin, end);
}

template<typename T>
void Nifti::writeVolumes(const size_t first, const size_t nvol,
                         typename std::vector<T>::iterator begin,
						 typename std::vector<T>::iterator end) {
	if (!(m_mode == Mode::Write))
		throw(std::runtime_error("File must be opened for writing: " + basePath()));
	Eigen::Array<size_t, 4, 1> start{0, 0, 0, first};
	Eigen::Array<size_t, 4, 1> size{dim(1), dim(2), dim(3), nvol};
	readWriteVoxels<T>(start, size, begin, end);
}

#endif // NIFTI_NIFTI_INL
