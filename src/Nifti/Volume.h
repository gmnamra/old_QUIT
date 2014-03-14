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
#include <exception>
#include <stdexcept>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Nifti/Nifti.h"

template<typename Tp>
class Volume {
	public:
		typedef Eigen::Array<size_t, 3, 1> IndexArray;
		typedef typename std::vector<Tp>::const_reference ConstTpRef;
		typedef typename std::vector<Tp>::reference TpRef;
	protected:
		std::vector<Tp> m_data;
		IndexArray      m_dims, m_strides;
		
		void calcStrides();
	public:
		Volume();
		Volume(const IndexArray &dims);
		Volume(Nifti &img, const size_t vol=0);
		
		void readFrom(Nifti &img, const size_t vol=0);
		void writeTo(Nifti &img, const size_t vol=0);
		
		const IndexArray &dims() const;
		
		ConstTpRef operator[](const IndexArray &vox) const;
		TpRef operator[](const IndexArray &vox);
};

template<typename Tp>
class VolumeSeries {
	public:
		typedef Eigen::Array<size_t, 4, 1> IndexArray;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> SeriesTp;
	protected:
		std::vector<Tp> m_data;
		IndexArray      m_dims, m_strides;
		
		void calcStrides();
	public:
		VolumeSeries();
		VolumeSeries(const IndexArray &dims);
		VolumeSeries(const typename Volume<Tp>::IndexArray &dims, const size_t nt);
		VolumeSeries(Nifti &img);
		
		void readFrom(Nifti &img);
		void writeTo(Nifti &img);
		
		const IndexArray &dims() const;
		
		const SeriesTp series(const typename Volume<Tp>::IndexArray &vox) const;
		SeriesTp series(const typename Volume<Tp>::IndexArray &vox);
};

#include "Volume-inl.h"

#endif //VOLUME_VOLUME
