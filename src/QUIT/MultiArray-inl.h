//
//  MultiArray-inl.h
//  MultiArray
//
//  Created by Tobias Wood on 10/03/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef MULTIARRAY_INL_H
#define MULTIARRAY_INL_H

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::CalcStrides(const Index &dims) -> Index {
	Index strides;
	strides[0] = 1;
	for (typename Index::Index i = 1; i < dims.size(); i++)
		strides[i] = strides[i - 1] * dims[i - 1];
	return strides;
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray() :
	m_offset{0},
	m_dims{Index::Zero()},
	m_strides{Index::Zero()},
	m_packed{true}
{
	
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const Index &inDims) :
	m_offset{0},
	m_dims{inDims},
	m_strides{CalcStrides(inDims)},
	m_ptr{std::make_shared<std::vector<Tp>>(inDims.prod())},
	m_packed{true}
{

}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const Eigen::Array<size_t, rank - 1, 1> &inDims, const size_t finalDim) :
	m_offset{0},
	m_packed{true}
{
	m_dims.head(rank - 1) = inDims;
	m_dims[rank - 1] = finalDim;
	m_strides = CalcStrides(m_dims);
	m_ptr = std::make_shared<std::vector<Tp>>(m_dims.prod());
}

template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArray(const Index &dims, const PtrTp &ptr, const Index &strides, const size_t offset) :
	m_dims{dims},
	m_strides{strides},
	m_offset{offset},
	m_ptr{ptr}
{
	if ((m_strides == Index::Zero()).all()) {
		m_strides = CalcStrides(m_dims);
	}
	m_packed = (CalcStrides(m_dims) == m_strides).all();
	//std::cout << "Created new MultiArray" << std::endl;
	//std::cout << *this << std::endl;
}

template<typename Tp, size_t rank> auto MultiArray<Tp, rank>::dims()     const -> const Index & { return m_dims; }
template<typename Tp, size_t rank> auto MultiArray<Tp, rank>::strides()  const -> const Index & { return m_strides; }
template<typename Tp, size_t rank> size_t MultiArray<Tp, rank>::size()   const { return m_dims.prod(); }
template<typename Tp, size_t rank> bool MultiArray<Tp, rank>::isPacked() const { return m_packed; }

template<typename Tp, size_t rank> void MultiArray<Tp, rank>::resize(const Index &newDims) {
	*this = MultiArray<Tp, rank>(newDims);
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
template<size_t newRank>
MultiArray<Tp, newRank> MultiArray<Tp, rank>::reshape(const typename MultiArray<Tp, newRank>::Index &newDims) {
	if (newDims.prod() != m_dims.prod()) {
		throw(std::logic_error("Reshape cannot change the number of voxels."));
	}
	return MultiArray<Tp, newRank>(newDims, m_ptr);
}


template<typename Tp, size_t rank>
template<size_t newRank>
MultiArray<Tp, newRank> MultiArray<Tp, rank>::slice(const Index &start, const Index &inSize) const {
	typename MultiArray<Tp, newRank>::Index newDims, newStrides;

	// Replace any "ALL" dimensions with actual size
	Index size;
	for (size_t i = 0; i < rank; i++) {
		if (inSize[i] == size_t(-1)) {
			size[i] = m_dims[i];
		} else {
			size[i] = inSize[i];
		}
	}
	if (((start + size) > m_dims).any()) {
		std::stringstream mesg; mesg << "Requested rank " << std::to_string(newRank) << " slice (start: " << start.transpose()
			<< ", size: " << size.transpose() << ") exceeds array dimensions: " << m_dims.transpose();
		throw(std::out_of_range(mesg.str()));
	}
	// Check that we have the correct number of 0 dimensions
	size_t zero_dim = 0;
	for (size_t d = 0; d < rank; d++) {
		if (size[d] == 0) zero_dim++;
	}
	if (newRank != (rank - zero_dim)) {
		throw(std::out_of_range("Incorrect number of zero dimensions for slice of rank " + std::to_string(newRank)));
	}
	// Now go through and copy over dimensions/strides for slice
	size_t to_dim = 0, from_dim = 0;
	while(from_dim < rank) {
		if (size[from_dim] > 0) {
			newDims[to_dim] = m_dims[from_dim];
			newStrides[to_dim] = m_strides[from_dim];
			to_dim++;
		}
		from_dim++;
	}
	size_t newOffset = m_offset + (m_strides*start).sum();
	MultiArray<Tp, newRank> slice(newDims, m_ptr, newStrides, newOffset);
	return slice;
}

// Can't partially specialize this, just not allowed :-(
template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::asArray() const -> MapTp {
	auto ptr = m_ptr->data() + m_offset;
	if (rank == 1) {
		// Outer and inner strides are reversed in Eigen constructor
		const Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> strides(0, m_strides[0]);
		MapTp array(ptr, m_dims[0], 1, strides);
		return array;
	} else if (rank == 2) {
		Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> strides(m_strides[1], m_strides[0]);
		MapTp array(ptr, m_dims[0], m_dims[1], strides);
		return array;
	} else {
		throw(std::logic_error("Cannot convert to an Eigen::Array if more than 2 dimensions."));
	}
}

template<typename Tp, size_t rank>
std::string MultiArray<Tp, rank>::print() const {
	std::stringstream ss;
	if (m_ptr) {
		ss << "MultiArray @" << m_ptr->data() << "+" << m_offset << ", share count: " << m_ptr.use_count() << std::endl;
		ss << "Dims:    " << m_dims.transpose() << std::endl;
		ss << "Strides: " << m_strides.transpose();
	} else {
		ss << "Unitialised MultiArray.";
	}
	return ss.str();
}

/******************************************************************************
 *
 * Iterator Methods
 *
 *****************************************************************************/
template<typename Tp, size_t rank>
MultiArray<Tp, rank>::MultiArrayIterator::MultiArrayIterator(MultiArray &array, Index start) :
	m_array(array),
	m_index(start),
	m_packedIndex(start.prod())
{

}

template<typename Tp, size_t rank>
Tp &MultiArray<Tp, rank>::iterator::operator*() {
	if (m_array.isPacked()) {
		return m_array[m_packedIndex];
	} else {
		return m_array[m_index];
	}
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::iterator::operator++() -> iterator & {
	if (m_array.isPacked()) {
		m_packedIndex++;
	} else {
		for (size_t dim = 0; dim < rank; dim++) {
			m_index[dim]++;
			if (m_index[dim] == m_array.dims()[dim]) {
				m_index[dim] = 0; // Reset this dim to zero, go and increment next dimension
			} else {
				break; // This dimension still has increments left
			}
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
	if (m_array.isPacked() && (m_packedIndex == other.m_packedIndex)) {
		return true;
	} else if ((m_index == other.m_index).all()) {
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
	iterator b(*this);
	return b;
}

template<typename Tp, size_t rank>
auto MultiArray<Tp, rank>::end() -> iterator {
	iterator e(*this, m_dims);
	return e;
}

#endif // MULTIARRAY_INL_H
