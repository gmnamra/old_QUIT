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
Volume<Tp>::Volume(const Eigen::Array<size_t, 3, 1> dims, const size_t nt) {
	assert(nt > 0);
	m_dims.head(3) = dims;
	m_dims(3) = nt;
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
	m_dims.setOnes();
	m_dims.head(img.dims().size()) = img.dims();
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
size_t Volume<Tp>::calcIndex(const IndexArray &aind) const {
	if (!(aind < m_dims.head(aind.size())).all())
		throw (std::out_of_range("Requested index outside of volume."));
	size_t index = (aind * m_strides.head(aind.size())).sum();
	return index;
}

template<typename Tp>
typename Volume<Tp>::ConstTpRef Volume<Tp>::operator[](const IndexArray &indices) const {
	return m_data[calcIndex(indices)];
}

template<typename Tp>
typename Volume<Tp>::TpRef Volume<Tp>::operator[](const IndexArray &indices) {
	return m_data[calcIndex(indices)];
}

template<typename Tp>
const typename Volume<Tp>::SeriesTp Volume<Tp>::series(const IndexArray &indices) const {
	size_t idx = calcIndex(indices);
	const SeriesTp s(m_data.data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp>
typename Volume<Tp>::SeriesTp Volume<Tp>::series(const IndexArray &indices) {
	size_t idx = calcIndex(indices);
	SeriesTp s(m_data.data() + idx, m_dims[3], Eigen::InnerStride<>(m_strides[3]));
	return s;
}

template<typename Tp>
const typename Volume<Tp>::IndexArray &Volume<Tp>::dims() const {
	return m_dims;
}

#endif
