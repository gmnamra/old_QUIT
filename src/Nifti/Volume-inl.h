//
//  Volume-inl.h
//  Volume
//
//  Created by Tobias Wood on 10/03/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef VOLUME_VOLUME_INL
#define VOLUME_VOLUME_INL

template<typename Tp>
void Volume<Tp>::calcStrides() {
	m_strides[0] = 1;
	for (IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp>
Volume<Tp>::Volume() :
	m_offset{0},
	m_dims{0, 0, 0},
	m_strides{0, 0, 0}
{
	
}

template<typename Tp>
Volume<Tp>::Volume(const IndexArray &inDims) :
	m_offset{0},
	m_dims{inDims}
{
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp>
Volume<Tp>::Volume(Nifti &img, const size_t vol) {
	readFrom(img, vol);
}

template<typename Tp>
Volume<Tp> Volume<Tp>::view(const IndexArray &start, const IndexArray &size, const IndexArray &viewStride) {
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
}

template<typename Tp>
void Volume<Tp>::readFrom(Nifti &img, const size_t vol) {
	m_dims = img.dims().head(3);
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
	img.readVolumes<Tp>(vol, 1, m_ptr->begin() + m_offset, m_ptr->end());
}

template<typename Tp>
void Volume<Tp>::writeTo(Nifti &img, const size_t vol) {
	img.writeVolumes<Tp>(vol, 1, m_ptr->begin() + m_offset, m_ptr->end());
}

template<typename Tp>
const typename Volume<Tp>::IndexArray &Volume<Tp>::dims() const {
	return m_dims;
}

template<typename Tp>
size_t Volume<Tp>::size() const {
	return m_dims.prod();
}

template<typename Tp>
typename Volume<Tp>::ConstTpRef Volume<Tp>::operator[](const IndexArray &vox) const {
	assert((vox < m_dims).all());
	return (*m_ptr)[(vox * m_strides).sum()];
}

template<typename Tp>
typename Volume<Tp>::TpRef Volume<Tp>::operator[](const IndexArray &vox) {
	assert((vox < m_dims).all());
	return (*m_ptr)[(vox * m_strides).sum()];
}

template<typename Tp>
typename Volume<Tp>::ConstTpRef Volume<Tp>::operator[](const size_t i) const {
	assert(i < size());
	return (*m_ptr)[i];
}

template<typename Tp>
typename Volume<Tp>::TpRef Volume<Tp>::operator[](const size_t i) {
	assert(i < size());
	return (*m_ptr)[i];
}


#pragma mark VolumeSeries

template<typename Tp>
VolumeSeries<Tp>::VolumeSeries() {
	
}

template<typename Tp>
VolumeSeries<Tp>::VolumeSeries(const IndexArray &inDims) {
	m_dims = inDims;
	calcStrides();
	m_data.resize(m_dims.prod());
}
	
template<typename Tp>
VolumeSeries<Tp>::VolumeSeries(const typename Volume<Tp>::IndexArray &dims, const size_t nt) {
	assert(nt > 0);
	m_dims.head(3) = dims;
	m_dims(3) = nt;
	calcStrides();
	m_data.resize(m_dims.prod());
}

template<typename Tp>
VolumeSeries<Tp>::VolumeSeries(Nifti &img) {
	readFrom(img);
}

template<typename Tp>
void VolumeSeries<Tp>::readFrom(Nifti &img) {
	if (img.rank() < 4) {
		m_dims.setOnes();
		m_dims.head(img.rank()) = img.dims();
	} else {
		m_dims = img.dims().head(4);
	}
	calcStrides();
	m_data.resize(m_dims.prod());
	img.readVolumes<Tp>(0, img.dim(4), m_data.begin(), m_data.end());
}

template<typename Tp>
void VolumeSeries<Tp>::readVolumesFrom(Nifti &img, size_t first, size_t n) {
	assert (img.dim(4) == n);
	auto b = m_data.begin() + m_strides[3]*first;
	auto e = b + m_strides[3]*n;
	img.readVolumes<Tp>(0, n, b, e);
}

template<typename Tp>
void VolumeSeries<Tp>::writeTo(Nifti &img) {
	assert(img.dim(4) == m_dims[3]);
	img.writeVolumes<Tp>(0, m_dims[3], m_data.begin(), m_data.end());
}

template<typename Tp>
void VolumeSeries<Tp>::writeVolumesTo(Nifti &img, size_t first, size_t n) {
	assert (img.dim(4) == n);
	auto b = m_data.begin() + m_strides[3]*first;
	auto e = b + m_strides[3]*n;
	img.writeVolumes<Tp>(0, n, b, e);
}

template<typename Tp>
void VolumeSeries<Tp>::calcStrides() {
	m_strides[0] = 1;
	for (IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp>
size_t VolumeSeries<Tp>::size() const {
	return m_dims.prod();
}

template<typename Tp>
const typename VolumeSeries<Tp>::SeriesTp VolumeSeries<Tp>::series(const size_t i) const {
	assert(i < size());
	const SeriesTp s(m_data.data() + i, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp>
typename VolumeSeries<Tp>::SeriesTp VolumeSeries<Tp>::series(const size_t i) {
	assert(i < size());
	SeriesTp s(m_data.data() + i, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}


template<typename Tp>
const typename VolumeSeries<Tp>::SeriesTp VolumeSeries<Tp>::series(const typename Volume<Tp>::IndexArray &vox) const {
	assert((vox < m_dims.head(3)).all());
	size_t idx = (vox * m_strides.head(3)).sum();
	const SeriesTp s(m_data.data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp>
typename VolumeSeries<Tp>::SeriesTp VolumeSeries<Tp>::series(const typename Volume<Tp>::IndexArray &vox) {
	assert((vox < m_dims.head(3)).all());
	size_t idx = (vox * m_strides.head(3)).sum();
	SeriesTp s(m_data.data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp>
const typename VolumeSeries<Tp>::IndexArray &VolumeSeries<Tp>::dims() const {
	return m_dims;
}

#endif
