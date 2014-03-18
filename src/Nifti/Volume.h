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
		typedef typename std::vector<Tp>::const_reference ConstTpRef;
		typedef typename std::vector<Tp>::reference TpRef;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> SeriesTp;
		
		static const size_t MaxIndex{std::numeric_limits<size_t>::max()};
	protected:
		std::shared_ptr<std::vector<Tp>> m_ptr;
		size_t      m_offset;
		IndexArray  m_dims, m_strides;
		
		void calcStrides();
	public:
		VolumeBase();
		VolumeBase(const VolumeBase &ov);
		VolumeBase(const VolumeBase &&ov);
		VolumeBase(const IndexArray &dims);
		VolumeBase(Nifti &img);
		
		VolumeBase view(const IndexArray &start, const IndexArray &end = IndexArray{MaxIndex, MaxIndex, MaxIndex}, const IndexArray &stride = IndexArray{1, 1, 1});
		VolumeBase copy(const IndexArray &start, const IndexArray &end, const IndexArray &stride);
		
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
		
		const SeriesTp series(const IndexArray &vox) const;
		SeriesTp series(const IndexArray &vox);
};

typedef VolumeBase<float, 3> Volumef;
typedef VolumeBase<bool, 3>  Volumeb;
typedef VolumeBase<float, 4> Seriesf;
typedef VolumeBase<std::complex<float>, 4> Seriescf;

template<typename Tp> class Volume : public VolumeBase<Tp, 3> {};
template<typename Tp> class Series : public VolumeBase<Tp, 4> {};
#include "Volume-inl.h"

#endif //VOLUME_VOLUME
