//
//  Volume.h
//  Volume
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
class VolumeBase {
	public:
		typedef Eigen::Array<size_t, rank, 1> IndexArray;
		typedef Eigen::Array<size_t, rank - 1, 1> ViewIndexArray;
		typedef std::vector<Tp> StorageTp;
		typedef std::shared_ptr<StorageTp> PtrTp;
		typedef typename StorageTp::const_reference ConstTpRef;
		typedef typename StorageTp::reference TpRef;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> SeriesTp;
		
		static const size_t MaxIndex{std::numeric_limits<size_t>::max()};
	protected:
		PtrTp       m_ptr;
		size_t      m_offset;
		IndexArray  m_dims, m_strides;
		
		void calcStrides();
	public:
		VolumeBase();
		VolumeBase(const IndexArray &dims);
		VolumeBase(const IndexArray &dims, const size_t offset, const PtrTp &ptr);
		VolumeBase(const ViewIndexArray &dims, const size_t finalDim);
		VolumeBase(Nifti &img);
		
		VolumeBase view(const IndexArray &start, const IndexArray &end = IndexArray{MaxIndex, MaxIndex, MaxIndex}, const IndexArray &stride = IndexArray{1, 1, 1});
		VolumeBase copy(const IndexArray &start, const IndexArray &end, const IndexArray &stride);
		VolumeBase<Tp, rank - 1> view(const size_t i);
		
		void readFrom(Nifti &img);
		void writeTo(Nifti &img);
		
		const IndexArray &dims() const;
		size_t size() const;
		
		ConstTpRef operator[](const size_t i) const;
		TpRef operator[](const size_t i);
		
		ConstTpRef operator[](const IndexArray &vox) const;
		TpRef operator[](const IndexArray &vox);
		
		const SeriesTp series(const size_t i) const;
		SeriesTp series(const size_t i);
		
		const SeriesTp series(const ViewIndexArray &vox) const;
		SeriesTp series(const ViewIndexArray &vox);
};

template<typename Tp> using Volume = VolumeBase<Tp, 3>;
template<typename Tp> using Series = VolumeBase<Tp, 4>;

#include "Volume-inl.h"

#endif //VOLUME_VOLUME
