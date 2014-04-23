//
//  MultiArray.h
//  MultiArray
//
//  Created by Tobias Wood on 10/03/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef VOLUME_VOLUME
#define VOLUME_VOLUME

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <exception>
#include <stdexcept>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Nifti/Nifti.h"

template<typename Tp, size_t rank>
class MultiArray {
	public:
		typedef Eigen::Array<size_t, rank, 1> Index;
		typedef Eigen::Array<size_t, rank - 1, 1> SliceIndex;
		typedef std::vector<Tp> StorageTp;
		typedef std::shared_ptr<StorageTp> PtrTp;
		typedef MultiArray<Tp, rank - 1> SliceTp;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> LineTp;
		// These typedefs are for STL Iterator compatibility
		typedef Tp        value_type;
		typedef size_t    size_type;
		typedef ptrdiff_t difference_type;
		typedef typename StorageTp::const_reference const_reference;
		typedef typename StorageTp::reference       reference;
		
		class MultiArrayIterator {
			public:
				typedef std::forward_iterator_tag iterator_category;
				typedef Tp        value_type;
				typedef size_t    size_type;
				typedef ptrdiff_t difference_type;
				typedef Tp        *pointer;
				typedef const Tp  *const_pointer;
				typedef Tp        &reference;
				typedef const Tp  &const_reference;

			private:
				MultiArray &m_array;
				Index m_index;

			public:
				MultiArrayIterator(MultiArray &array, Index start);
				Tp &operator*();

				MultiArrayIterator &operator++();
				MultiArrayIterator operator++(int);

				bool operator==(const MultiArrayIterator &other) const;
				bool operator!=(const MultiArrayIterator &other) const;
		};
		typedef MultiArrayIterator iterator;

		static const size_t MaxIndex{std::numeric_limits<size_t>::max()};
	protected:
		PtrTp   m_ptr;
		size_t  m_offset;
		Index    m_dims, m_strides;
		
		void calcStrides();
	public:
		MultiArray();
		MultiArray(const Index &dims);
		MultiArray(const Index &dims, const Index &strides, const size_t offset, const PtrTp &ptr);
		MultiArray(const SliceIndex &dims, const size_t finalDim);
		MultiArray(Nifti &img);
		
		void readFrom(Nifti &img);
		void writeTo(Nifti &img);
		
		const Index &dims() const;
		const Index &strides() const;
		size_t size() const;
		
		SliceTp viewSlice(const size_t i, const size_t d=rank);

		LineTp line(const SliceIndex &vox, const size_t d=rank) const;
		LineTp line(const size_t i) const;

		// STL-like interface
		const_reference operator[](const size_t i) const;
		const_reference operator[](const Index &vox) const;
		reference operator[](const size_t i);
		reference operator[](const Index &vox);
		
		iterator begin();
		iterator end();

		std::string print() const;
		friend std::ostream &operator<<(std::ostream &os, const MultiArray &v) {
			os << v.print();
			return os;
		}
};

template<typename Tp> using Volume = MultiArray<Tp, 3>;
template<typename Tp> using Series = MultiArray<Tp, 4>;

#include "MultiArray-inl.h"

#endif //VOLUME_VOLUME
