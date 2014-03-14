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
Volume<Tp>::Volume() {
	
}

template<typename Tp>
Volume<Tp>::Volume(const IndexArray &inDims) {
	m_dims = inDims;
	calcStrides();
	m_data.resize(m_dims.prod());
}

template<typename Tp>
Volume<Tp>::Volume(Nifti &img, const size_t vol) {
	readFrom(img, vol);
}

template<typename Tp>
void Volume<Tp>::readFrom(Nifti &img, const size_t vol) {
	m_dims = img.dims().head(3);
	calcStrides();
	m_data.resize(m_dims.prod());
	img.readVolumes<Tp>(vol, 1, m_data.begin(), m_data.end());
}

template<typename Tp>
void Volume<Tp>::writeTo(Nifti &img, const size_t vol) {
	img.writeVolumes<Tp>(vol, 1, m_data.begin(), m_data.end());
}

template<typename Tp>
void Volume<Tp>::calcStrides() {
	m_strides[0] = 1;
	for (IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp>
const typename Volume<Tp>::IndexArray &Volume<Tp>::dims() const {
	return m_dims;
}

template<typename Tp>
typename Volume<Tp>::ConstTpRef Volume<Tp>::operator[](const IndexArray &vox) const {
	assert((vox < m_dims).all());
	return m_data[(vox * m_strides).sum()];
}

template<typename Tp>
typename Volume<Tp>::TpRef Volume<Tp>::operator[](const IndexArray &vox) {
	assert((vox < m_dims).all());
	return m_data[(vox * m_strides).sum()];
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
	m_dims = img.dims().head(4);
	calcStrides();
	m_data.resize(m_dims.prod());
	img.readVolumes<Tp>(0, img.dim(4), m_data.begin(), m_data.end());
}

template<typename Tp>
void VolumeSeries<Tp>::writeTo(Nifti &img) {
	assert(img.dim(4) == m_dims[3]);
	img.writeVolumes<Tp>(0, m_dims[3], m_data.begin(), m_data.end());
}

template<typename Tp>
void VolumeSeries<Tp>::calcStrides() {
	m_strides[0] = 1;
	for (IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
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
