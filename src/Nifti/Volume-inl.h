//
//  Volume-inl.h
//  Volume
//
//  Created by Tobias Wood on 10/03/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef VOLUME_VOLUME_INL
#define VOLUME_VOLUME_INL

template<typename Tp, size_t rank>
void VolumeBase<Tp, rank>::calcStrides() {
	m_strides[0] = 1;
	for (typename Indx::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase() :
	m_offset{0},
	m_dims{Indx::Zero()},
	m_strides{Indx::Zero()}
{
	
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase(const Indx &inDims) :
	m_offset{0},
	m_dims{inDims}
{
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase(const SliceIndx &inDims, const size_t finalDim) :
	m_offset{0}
{
	m_dims.head(rank - 1) = inDims;
	m_dims[rank - 1] = finalDim;
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase(Nifti &img) {
	readFrom(img);
}

template<typename Tp, size_t rank>
void VolumeBase<Tp, rank>::readFrom(Nifti &img) {
	assert(rank == img.rank());
	m_offset = 0;
	m_dims = img.dims().head(rank);
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
	img.readVolumes<Tp>(0, img.dim(4), m_ptr->begin() + m_offset, m_ptr->end());
}

template<typename Tp, size_t rank>
void VolumeBase<Tp, rank>::writeTo(Nifti &img) {
	// We might be a view, so work these out
	auto begin = m_ptr->begin() + m_offset;
	auto end   = begin + size();
	assert(end <= m_ptr->end());
	img.writeVolumes<Tp>(0, img.dim(4), begin, end);
}

template<typename Tp, size_t rank>
const typename VolumeBase<Tp, rank>::Indx &VolumeBase<Tp, rank>::dims() const {
	return m_dims;
}

template<typename Tp, size_t rank>
size_t VolumeBase<Tp, rank>::size() const {
	return m_dims.prod();
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::ConstTpRef VolumeBase<Tp, rank>::operator[](const Indx &vox) const {
	assert((vox < m_dims).all());
	return (*m_ptr)[m_offset + (vox * m_strides).sum()];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::TpRef VolumeBase<Tp, rank>::operator[](const Indx &vox) {
	assert((vox < m_dims).all());
	return (*m_ptr)[m_offset + (vox * m_strides).sum()];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::ConstTpRef VolumeBase<Tp, rank>::operator[](const size_t i) const {
	assert(i < size());
	return (*m_ptr)[m_offset + i];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::TpRef VolumeBase<Tp, rank>::operator[](const size_t i) {
	assert(i < size());
	return (*m_ptr)[m_offset + i];
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase(const Indx &dims, const Indx &strides, const size_t offset, const PtrTp &ptr) :
	m_dims{dims}, m_strides{strides}, m_offset(offset), m_ptr(ptr) {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	//std::cout << *this << std::endl;
}

template<typename Tp, size_t rank>
auto VolumeBase<Tp, rank>::viewSlice(const size_t i, const size_t d) -> SliceTp {
	SliceIndx newDims, newStrides;
	size_t to_dim = 0, from_dim = 0;
	while(from_dim < rank) {
		if (from_dim != (d-1)) {
			newDims[to_dim] = m_dims[from_dim];
			newStrides[to_dim] = m_strides[from_dim];
			to_dim++;
		} else {
		
		}
		from_dim++;
	}
	size_t sliceOffset = m_offset + m_strides[d-1] * i;
	return SliceTp{newDims, newStrides, sliceOffset, m_ptr};
}

template<typename Tp, size_t rank>
auto VolumeBase<Tp, rank>::line(const size_t i) const -> LineTp {
	assert(i < m_dims.head(rank-1).prod());
	auto first = m_ptr->data() + m_offset + i;
	const LineTp s(first, m_dims[rank-1], Eigen::InnerStride<>(m_strides[rank-1]));
	return s;
}

template<typename Tp, size_t rank>
auto VolumeBase<Tp, rank>::line(const SliceIndx &vox, const size_t lineD) const -> LineTp {
	size_t idx = 0, d1 = 0, d2 = 0;
	while (d1 < rank) {
		if (d1 != (lineD-1)) {
			if (vox[d2] > m_dims[d2]) {
				throw(std::out_of_range("Requested line outside of volume dimensions."));
			}
			idx += m_strides[d1] * vox[d2];
			d2++;
		}
		d1++;
	}
	auto first = m_ptr->data() + m_offset + idx;
	const LineTp s(first, m_dims[lineD-1], Eigen::InnerStride<>(m_strides[lineD-1]));
	return s;
}

#endif
