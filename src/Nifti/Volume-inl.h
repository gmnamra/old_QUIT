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
Volume<Tp>::Volume(std::vector<size_t> inDims) {
	assert(inDims.size() <= 4);
	m_dims.resize(inDims.size());
	for (size_t i = 0; i < inDims.size(); i++)
		m_dims[i] = inDims[i];
	calcStrides();
	m_data.resize(m_dims.prod());
}
	
template<typename Tp>
Volume<Tp>::Volume(const Volume<Tp>::IndexArray inDims) {
	assert(inDims.size() <= 4);
	m_dims = inDims;
	calcStrides();
	m_data.resize(m_dims.prod());
}

template<typename Tp>
Volume<Tp>::Volume(Nifti &img) {
	readFrom(img);
}

template<typename Tp>
void Volume<Tp>::readFrom(Nifti &img) {
	assert(img.rank() <= 4);
	m_dims = img.dims();
	calcStrides();
	m_data.resize(m_dims.prod());
	img.readVolumes(0, img.dim(4), m_data);
}

template<typename Tp>
void Volume<Tp>::writeTo(Nifti &img) {
	assert(img.rank() <= 4);
	if (img.rank() == 4) {
		assert(img.dim(4) == m_dims[3]);
		img.writeVolumes<Tp>(0, m_dims[3], m_data);
	} else {
		img.writeVolumes<Tp>(0, 1, m_data);
	}
}

template<typename Tp>
void Volume<Tp>::calcStrides() {
	m_strides.resize(m_dims.size());
	m_strides[0] = 1;
	for (IndexArray::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp>
size_t Volume<Tp>::calcIndex(const std::vector<size_t> indices) const {
	IndexArray aind(indices.size());
	for (size_t i = 0; i < indices.size(); i++)
		aind[i] = indices[i];
	if (!(aind < m_dims.head(aind.size())).all())
		throw (std::out_of_range("Requested index outside of volume."));
	size_t index = (aind * m_strides.head(aind.size())).sum();
	return index;
}

template<typename Tp>
const Tp &Volume<Tp>::operator[](const std::vector<size_t> indices) const {
	return m_data[calcIndex(indices)];
}

template<typename Tp>
Tp &Volume<Tp>::operator[](const std::vector<size_t> indices) {
	return m_data[calcIndex(indices)];
}

template<typename Tp>
typename Volume<Tp>::VectorTp Volume<Tp>::series(const size_t index) const {
	VectorTp s(m_dims[3]);
	size_t idx = index;
	for (size_t i = 0; i < m_dims[3]; i++) {
		s[i] = m_data.at(idx);
		idx += m_strides[3];
	}
	return s;
}

template<typename Tp>
typename Volume<Tp>::VectorTp Volume<Tp>::series(const std::vector<size_t> indices) const {
	VectorTp s(m_dims[3]);
	size_t idx = calcIndex(indices);
	for (size_t i = 0; i < m_dims[3]; i++) {
		s[i] = m_data.at(idx);
		idx += m_strides[3];
	}
	return s;
}

template<typename Tp>
const typename Volume<Tp>::IndexArray &Volume<Tp>::dims() const {
	return m_dims;
}

#endif
