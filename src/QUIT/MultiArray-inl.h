//
//  MultiArray-inl.h
//  MultiArray
//
//  Created by Tobias Wood on 10/03/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef VOLUME_VOLUME_INL
#define VOLUME_VOLUME_INL

template<typename Tp, size_t rank>
void MultiArray<Tp, rank>::calcStrides() {
	m_strides[0] = 1;
	for (typename Index::Index i = 1; i < m_dims.size(); i++)
		m_strides[i] = m_strides[i - 1] * m_dims[i - 1];
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray() :
	m_offset{0},
	m_dims{Index::Zero()},
	m_strides{Index::Zero()}
{
	
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const Index &inDims) :
	m_offset{0},
	m_dims{inDims}
{
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const SliceIndex &inDims, const size_t finalDim) :
	m_offset{0}
{
	m_dims.head(rank - 1) = inDims;
	m_dims[rank - 1] = finalDim;
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(Nifti &img) {
	readFrom(img);
}

template<typename Tp, size_t rank>
void MultiArray<Tp, rank>::readFrom(Nifti &img) {
	assert(rank == img.rank());
	m_offset = 0;
	m_dims = img.dims().head(rank);
	calcStrides();
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
	img.readVolumes(0, img.dim(4), m_ptr->begin() + m_offset, m_ptr->end());
}

template<typename Tp, size_t rank>
void MultiArray<Tp, rank>::writeTo(Nifti &img) {
	// We might be a view, so work these out
	auto begin = m_ptr->begin() + m_offset;
	auto end   = begin + size();
	assert(end <= m_ptr->end());
	img.writeVolumes(0, img.dim(4), begin, end);
}

template<typename Tp, size_t rank>
const typename MultiArray<Tp, rank>::Index &MultiArray<Tp, rank>::dims() const {
	return m_dims;
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::strides() const -> const Index & {
	return m_strides;
}

template<typename Tp, size_t rank>
size_t MultiArray<Tp, rank>::size() const {
	return m_dims.prod();
}

template<typename Tp, size_t rank>
typename MultiArray<Tp, rank>::const_reference MultiArray<Tp, rank>::operator[](const Index &vox) const {
	if ((vox >= m_dims).any()) {
		std::stringstream ss;
		ss << "Voxel " << vox.transpose() << " outside volume.\n" << print();
		throw(std::out_of_range(ss.str()));
	}
	return (*m_ptr)[m_offset + (vox * m_strides).sum()];
}

template<typename Tp, size_t rank>
typename MultiArray<Tp, rank>::reference MultiArray<Tp, rank>::operator[](const Index &vox) {
	return const_cast<reference>(static_cast<const MultiArray<Tp, rank> &>(*this).operator[](vox));
}

template<typename Tp, size_t rank>
typename MultiArray<Tp, rank>::const_reference MultiArray<Tp, rank>::operator[](const size_t i) const {
	if (i >= size()) {
		throw(std::out_of_range("Index " + std::to_string(i) + " out of range.\n" + print()));
	}
	return (*m_ptr)[m_offset + i];
}

template<typename Tp, size_t rank>
typename MultiArray<Tp, rank>::reference MultiArray<Tp, rank>::operator[](const size_t i) {
	return const_cast<reference>(static_cast<const MultiArray<Tp, rank> &>(*this).operator[](i));
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const Index &dims, const Index &strides, const size_t offset, const PtrTp &ptr) :
	m_dims{dims}, m_strides{strides}, m_offset(offset), m_ptr(ptr) {
	
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::viewSlice(const size_t i, const size_t d) -> SliceTp {
	SliceIndex newDims, newStrides;
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
auto MultiArray<Tp, rank>::line(const size_t i) const -> LineTp {
	assert(i < m_dims.head(rank-1).prod());
	auto first = m_ptr->data() + m_offset + i;
	const LineTp s(first, m_dims[rank-1], Eigen::InnerStride<>(m_strides[rank-1]));
	return s;
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::line(const SliceIndex &vox, const size_t lineD) const -> LineTp {
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

template<typename Tp, size_t rank>
std::string MultiArray<Tp, rank>::print() const {
	std::stringstream ss;
	ss << "Dims:    " << m_dims.transpose() << std::endl;
	ss << "Strides: " << m_strides.transpose() << std::endl;
	ss << "Offset:  " << m_offset << " Ptr Count: " << m_ptr.use_count() << std::endl;
	return ss.str();
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArrayIterator::MultiArrayIterator(MultiArray &array, Index start) :
	m_array(array), m_index(start) {

}

template<typename Tp, size_t rank>
Tp &MultiArray<Tp, rank>::iterator::operator*() {
	return m_array[m_index];
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::iterator::operator++() -> iterator & {
	size_t dim = 0;
	for (size_t dim = 0; dim < rank; dim++) {
		m_index[dim]++;
		if (m_index[dim] == m_array.dims()[dim]) {
			m_index[dim] = 0; // Reset this dim to zero, go and increment next dimension
		} else {
			break; // This dimension still has increments left
		}
	}
	return *this;
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::iterator::operator++(int) -> iterator {
	iterator tmp(*this);
	operator++(); // Call pre-increment operator
	return tmp;
}

template<typename Tp, size_t rank>
bool MultiArray<Tp, rank>::iterator::operator==(const iterator &other) const {
	if ((m_index == other.m_index).all()) {
		return true;
	} else {
		return false;
	}
}

template<typename Tp, size_t rank>
bool MultiArray<Tp, rank>::iterator::operator!=(const iterator &other) const {
	return !operator==(other);
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::begin() -> iterator {
	iterator b(*this, Index::Zero());
	return b;
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::end() -> iterator {
	iterator e(*this, m_dims);
	return e;
}

#endif
