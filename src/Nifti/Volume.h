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
		typedef Eigen::Array<size_t, rank, 1> Indx;
		typedef Eigen::Array<size_t, rank - 1, 1> SliceIndx;
		typedef std::vector<Tp> StorageTp;
		typedef std::shared_ptr<StorageTp> PtrTp;
		typedef typename StorageTp::const_reference ConstTpRef;
		typedef typename StorageTp::reference TpRef;
		typedef VolumeBase<Tp, rank - 1> SliceTp;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> LineTp;
		
		static const size_t MaxIndex{std::numeric_limits<size_t>::max()};
	protected:
		PtrTp   m_ptr;
		size_t  m_offset;
		Indx    m_dims, m_strides;
		
		void calcStrides();
	public:
		VolumeBase();
		VolumeBase(const Indx &dims);
		VolumeBase(const Indx &dims, const Indx &strides, const size_t offset, const PtrTp &ptr);
		VolumeBase(const SliceIndx &dims, const size_t finalDim);
		VolumeBase(Nifti &img);
		
		void readFrom(Nifti &img);
		void writeTo(Nifti &img);
		
		const Indx &dims() const;
		size_t size() const;
		
		ConstTpRef operator[](const size_t i) const;
		TpRef operator[](const size_t i);
		
		ConstTpRef operator[](const Indx &vox) const;
		TpRef operator[](const Indx &vox);
		
		SliceTp viewSlice(const size_t i, const size_t d=rank);
		
		LineTp line(const SliceIndx &vox, const size_t d=rank) const;
		LineTp line(const size_t i) const;
		
		friend std::ostream &operator<<(std::ostream &os, const VolumeBase &v) {
			os << "Dims:    " << v.m_dims.transpose() << std::endl;
			os << "Strides: " << v.m_strides.transpose() << std::endl;
			os << "Offset:  " << v.m_offset << " Ptr Count: " << v.m_ptr.use_count() << std::endl;
			return os;
		}
};

template<typename Tp> using Volume = VolumeBase<Tp, 3>;
template<typename Tp> using Series = VolumeBase<Tp, 4>;

#include "Volume-inl.h"

#endif //VOLUME_VOLUME
