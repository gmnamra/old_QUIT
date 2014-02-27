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
	typename std::vector<T>::iterator::difference_type nEl = bytes.size() / m_typeinfo.size;
	assert(nEl == std::distance(begin, end));
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	#define DECL_PTR( BYTE_TYPE ) auto p = reinterpret_cast<const BYTE_TYPE *>(bytes.data())
	#define REAL_LOOP( BYTE_TYPE ) DECL_PTR(BYTE_TYPE); while (el != end) { *el = sc_in + sc_sl * static_cast<T>(*p); el++; p++; }
	#define COMP_LOOP( BYTE_TYPE ) DECL_PTR(std::complex<BYTE_TYPE>); while (el != end) { *el = sc_in + sc_sl * static_cast<T>(abs(*p)); el++; p++; }
	switch (m_typeinfo.type) {
		case DataType::INT8:       { REAL_LOOP(int8_t); }; break;
		case DataType::INT16:      { REAL_LOOP(int16_t); }; break;
		case DataType::INT32:      { REAL_LOOP(int32_t); }; break;
		case DataType::INT64:      { REAL_LOOP(int64_t); }; break;
		case DataType::UINT8:      { REAL_LOOP(uint8_t); }; break;
		case DataType::UINT16:     { REAL_LOOP(uint16_t); }; break;
		case DataType::UINT32:     { REAL_LOOP(uint32_t); }; break;
		case DataType::UINT64:     { REAL_LOOP(uint64_t); }; break;
		case DataType::FLOAT32:    { REAL_LOOP(float); }; break;
		case DataType::FLOAT64:    { REAL_LOOP(double); }; break;
		case DataType::FLOAT128:   { REAL_LOOP(long double); }; break;
		case DataType::COMPLEX64:  { COMP_LOOP(float); }; break;
		case DataType::COMPLEX128: { COMP_LOOP(double); }; break;
		case DataType::COMPLEX256: { COMP_LOOP(long double); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
	}
	#undef COMP_LOOP
	#undef REAL_LOOP
	#undef DECL_PTR
}

template<typename T> void Nifti::convertFromBytes(const std::vector<char> &bytes,
                                                  const typename std::vector<std::complex<T>>::iterator begin,
												  const typename std::vector<std::complex<T>>::iterator end) {
	size_t nEl = bytes.size() / m_typeinfo.size;
	assert(nEl == std::distance(begin, end));
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	#define DECL_PTR( BYTE_TYPE ) auto p = reinterpret_cast<const BYTE_TYPE *>(bytes.data())
	#define REAL_LOOP( BYTE_TYPE ) DECL_PTR(BYTE_TYPE); while (el != end) { *el = std::complex<T>(sc_in + sc_sl * static_cast<T>(*p), 0.); el++; p++; }
	#define COMP_LOOP( BYTE_TYPE ) DECL_PTR(std::complex<BYTE_TYPE>); while (el != end) { *el = sc_in + sc_sl * static_cast<std::complex<T>>(*p); el++; p++; }
	switch (m_typeinfo.type) {
		case DataType::INT8:       { REAL_LOOP(int8_t); }; break;
		case DataType::INT16:      { REAL_LOOP(int16_t); }; break;
		case DataType::INT32:      { REAL_LOOP(int32_t); }; break;
		case DataType::INT64:      { REAL_LOOP(int64_t); }; break;
		case DataType::UINT8:      { REAL_LOOP(uint8_t); }; break;
		case DataType::UINT16:     { REAL_LOOP(uint16_t); }; break;
		case DataType::UINT32:     { REAL_LOOP(uint32_t); }; break;
		case DataType::UINT64:     { REAL_LOOP(uint64_t); }; break;
		case DataType::FLOAT32:    { REAL_LOOP(float); }; break;
		case DataType::FLOAT64:    { REAL_LOOP(double); }; break;
		case DataType::FLOAT128:   { REAL_LOOP(long double); }; break;
		case DataType::COMPLEX64:  { COMP_LOOP(float); }; break;
		case DataType::COMPLEX128: { COMP_LOOP(double); }; break;
		case DataType::COMPLEX256: { COMP_LOOP(long double); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
	}
	#undef COMP_LOOP
	#undef REAL_LOOP
	#undef DECL_PTR
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
	#define DECL_PTR( BYTE_TYPE ) auto p = reinterpret_cast<BYTE_TYPE *>(bytes.data())
	#define REAL_LOOP( BYTE_TYPE ) DECL_PTR( BYTE_TYPE ); while (el != end) { *p = static_cast<BYTE_TYPE>(*el / sc_sl - sc_in); el++; p++; }
	#define COMP_LOOP( BYTE_TYPE ) DECL_PTR( std::complex<BYTE_TYPE> ); while (el != end) { *p = std::complex<BYTE_TYPE>(static_cast<BYTE_TYPE>(*el / sc_sl - sc_in), 0.); el++; p++; }
	switch (m_typeinfo.type) {
		case DataType::INT8:       { REAL_LOOP(int8_t); }; break;
		case DataType::INT16:      { REAL_LOOP(int16_t); }; break;
		case DataType::INT32:      { REAL_LOOP(int32_t); }; break;
		case DataType::INT64:      { REAL_LOOP(int64_t); }; break;
		case DataType::UINT8:      { REAL_LOOP(uint8_t); }; break;
		case DataType::UINT16:     { REAL_LOOP(uint16_t); }; break;
		case DataType::UINT32:     { REAL_LOOP(uint32_t); }; break;
		case DataType::UINT64:     { REAL_LOOP(uint64_t); }; break;
		case DataType::FLOAT32:    { REAL_LOOP(float); }; break;
		case DataType::FLOAT64:    { REAL_LOOP(double); }; break;
		case DataType::FLOAT128:   { REAL_LOOP(long double); }; break;
		case DataType::COMPLEX64:  { COMP_LOOP(float); }; break;
		case DataType::COMPLEX128: { COMP_LOOP(double); }; break;
		case DataType::COMPLEX256: { COMP_LOOP(long double); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
	}
	#undef COMP_LOOP
	#undef REAL_LOOP
	#undef DECL_PTR
}

template<typename T> void Nifti::convertToBytes(const typename std::vector<std::complex<T>>::iterator begin,
												const typename std::vector<std::complex<T>>::iterator end,
												std::vector<char> &bytes) {
	size_t nEl = std::distance(begin, end);
	assert(nEl == bytes.size() / m_typeinfo.size);
	T sc_sl = static_cast<T>(scaling_slope);
	T sc_in = static_cast<T>(scaling_inter);
	auto el = begin;
	#define DECL_PTR( BYTE_TYPE ) auto p = reinterpret_cast<BYTE_TYPE *>(bytes.data())
	#define REAL_LOOP( BYTE_TYPE ) DECL_PTR( BYTE_TYPE ); while (el != end) { *p = static_cast<BYTE_TYPE>(abs(*el) / sc_sl - sc_in); el++; p++; }
	#define COMP_LOOP( BYTE_TYPE ) DECL_PTR( std::complex<BYTE_TYPE> ); while (el != end) { *p = static_cast<std::complex<BYTE_TYPE>>(*el / sc_sl - sc_in); el++; p++; }
	switch (m_typeinfo.type) {
		case DataType::INT8:       { REAL_LOOP(int8_t); }; break;
		case DataType::INT16:      { REAL_LOOP(int16_t); }; break;
		case DataType::INT32:      { REAL_LOOP(int32_t); }; break;
		case DataType::INT64:      { REAL_LOOP(int64_t); }; break;
		case DataType::UINT8:      { REAL_LOOP(uint8_t); }; break;
		case DataType::UINT16:     { REAL_LOOP(uint16_t); }; break;
		case DataType::UINT32:     { REAL_LOOP(uint32_t); }; break;
		case DataType::UINT64:     { REAL_LOOP(uint64_t); }; break;
		case DataType::FLOAT32:    { REAL_LOOP(float); }; break;
		case DataType::FLOAT64:    { REAL_LOOP(double); }; break;
		case DataType::FLOAT128:   { REAL_LOOP(long double); }; break;
		case DataType::COMPLEX64:  { COMP_LOOP(float); }; break;
		case DataType::COMPLEX128: { COMP_LOOP(double); }; break;
		case DataType::COMPLEX256: { COMP_LOOP(long double); }; break;
		case DataType::RGB24: case DataType::RGBA32:
			throw(std::runtime_error("Unsupported datatype.")); break;
	}
	#undef COMP_LOOP
	#undef REAL_LOOP
	#undef DECL_PTR
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
template<typename T> void Nifti::readWriteVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &inSize, std::vector<T> &data) {
	ArrayXs size = inSize;
	for (ArrayXs::Index i = 0; i < size.rows(); i++)
		if (size(i) == 0) size(i) = m_dim(i);
	
	if (start.rows() != size.rows()) throw(std::out_of_range("Start and size must have same dimension in image: " + imagePath()));
	if (start.rows() > m_dim.rows()) throw(std::out_of_range("Too many read/write dimensions specified in image: " + imagePath()));
	if (((start + size) > m_dim.head(start.rows())).any()) throw(std::out_of_range("Read/write past image dimensions requested: " + imagePath()));
	if (size.prod() < data.size()) throw(std::out_of_range("Allocated memory is insufficient for read/write in image: " + imagePath()));
	
	ArrayXs::Index firstDim = 0; // We can always read first dimension in one go
	ArrayXs::Index blockSize = size(firstDim);
	auto dataIt = data.begin();
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
				convertFromBytes<T>(block, dataIt, dataIt + blockSize);
			} else {
				convertToBytes<T>(dataIt, dataIt + blockSize, block);
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

template<typename T> void Nifti::readVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size, std::vector<T> &data) {
	if (!(m_mode == Mode::Read))
		throw(std::runtime_error("File must be opened for reading: " + basePath()));
	readWriteVoxels(start, size, data);
}

template<typename T> void Nifti::readVolumes(const size_t first, const size_t nvol, std::vector<T> &data) {
	if (!(m_mode == Mode::Read))
		throw(std::runtime_error("File must be opened for reading: " + basePath()));
	if (data.size() != (m_dim.head(3).prod() * nvol))
		throw(std::runtime_error("Insufficient storage allocated for read: " + basePath()));
	
	Eigen::Array<size_t, 4, 1> start, size;
	start << 0, 0, 0, first;
	size << dim(1), dim(2), dim(3), nvol;
	readWriteVoxels(start, size, data);
}

template<typename T> void Nifti::writeVoxels(const Eigen::Ref<ArrayXs> &start, const Eigen::Ref<ArrayXs> &size, std::vector<T> &data) {
	if (!(m_mode == Mode::Write))
		throw(std::runtime_error("File must be opened for writing: " + basePath()));
	readWriteVoxels(start, size, data);
}

template<typename T> void Nifti::writeVolumes(const size_t first, const size_t nvol, std::vector<T> &data) {
	if (!(m_mode == Mode::Write))
		throw(std::runtime_error("File must be opened for writing: " + basePath()));
	if (data.size() != (m_dim.head(3).prod() * nvol))
		throw(std::runtime_error("Insufficient data for write: " + basePath()));
	
	Eigen::Array<size_t, 4, 1> start, size;
	start << 0, 0, 0, first;
	size << dim(1), dim(2), dim(3), nvol;
	readWriteVoxels(start, size, data);
}

#endif // NIFTI_NIFTI_INL
