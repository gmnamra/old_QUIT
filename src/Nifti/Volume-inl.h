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
	for (typename IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase() :
	m_offset{0},
	m_dims{IndexArray::Zero()},
	m_strides{IndexArray::Zero()}
{
	
}

template<typename Tp, size_t rank>
VolumeBase<Tp, rank>::VolumeBase(const IndexArray &inDims) :
	m_offset{0},
	m_dims{inDims}
{
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
	img.writeVolumes<Tp>(0, img.dim(4), m_ptr->begin() + m_offset, m_ptr->end());
}

template<typename Tp, size_t rank>
const typename VolumeBase<Tp, rank>::IndexArray &VolumeBase<Tp, rank>::dims() const {
	return m_dims;
}

template<typename Tp, size_t rank>
size_t VolumeBase<Tp, rank>::size() const {
	return m_dims.prod();
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::ConstTpRef VolumeBase<Tp, rank>::operator[](const IndexArray &vox) const {
	assert((vox < m_dims).all());
	return (*m_ptr)[(vox * m_strides).sum()];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::TpRef VolumeBase<Tp, rank>::operator[](const IndexArray &vox) {
	assert((vox < m_dims).all());
	return (*m_ptr)[(vox * m_strides).sum()];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::ConstTpRef VolumeBase<Tp, rank>::operator[](const size_t i) const {
	assert(i < size());
	return (*m_ptr)[i];
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::TpRef VolumeBase<Tp, rank>::operator[](const size_t i) {
	assert(i < size());
	return (*m_ptr)[i];
}

/*
template<typename Tp, size_t rank>
VolumeBase<Tp, rank> VolumeBase<Tp, rank>::view(const IndexArray &start, const IndexArray &size, const IndexArray &viewStride) {
	Volume v;
	for (IndexArray::Index i = 0; i < size.size(); i++) {
		if (size[i] == MaxIndex) size[i] = m_dims[i];
	}
	assert((start < m_dims).all());
	assert(((start + size) < m_dims).all());
	
	v.m_ptr     = m_ptr;
	v.m_offset  = (start * m_strides).sum();
	v.m_dims    = size;
	v.m_strides = (m_strides * viewStride);
	return v;
}*/

/*
template<typename Tp>
void VolumeSeries<Tp>::readVolumesFrom(Nifti &img, size_t first, size_t n) {
	assert (img.dim(4) == n);
	auto b = m_data.begin() + m_strides[3]*first;
	auto e = b + m_strides[3]*n;
	img.readVolumes<Tp>(0, n, b, e);
}*/

/*
template<typename Tp>
void VolumeSeries<Tp>::writeVolumesTo(Nifti &img, size_t first, size_t n) {
	assert (img.dim(4) == n);
	auto b = m_data.begin() + m_strides[3]*first;
	auto e = b + m_strides[3]*n;
	img.writeVolumes<Tp>(0, n, b, e);
}*/

template<typename Tp, size_t rank>
const typename VolumeBase<Tp, rank>::SeriesTp VolumeBase<Tp, rank>::series(const size_t i) const {
	assert(rank > 3);
	assert(i < m_dims.head(3).prod());
	const SeriesTp s(m_ptr->data() + i, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::SeriesTp VolumeBase<Tp, rank>::series(const size_t i) {
	assert(rank > 3);
	assert(i < m_dims.head(3).prod());
	SeriesTp s(m_ptr->data() + i, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}


template<typename Tp, size_t rank>
const typename VolumeBase<Tp, rank>::SeriesTp VolumeBase<Tp, rank>::series(const typename VolumeBase<Tp, rank>::IndexArray &vox) const {
	assert(rank > 3);
	assert((vox < m_dims.head(3)).all());
	size_t idx = (vox * m_strides.head(3)).sum();
	const SeriesTp s(m_ptr->data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp, size_t rank>
typename VolumeBase<Tp, rank>::SeriesTp VolumeBase<Tp, rank>::series(const typename VolumeBase<Tp, rank>::IndexArray &vox) {
	assert(rank > 3);
	assert((vox.head(3) < m_dims.head(3)).all());
	size_t idx = (vox.head(3) * m_strides.head(3)).sum();
	SeriesTp s(m_ptr->data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

#endif
