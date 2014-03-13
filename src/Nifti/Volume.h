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
		typedef Eigen::Array<size_t, 4, 1> IndexArray;
		typedef Eigen::Map<Eigen::Array<Tp, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<>> SeriesTp;
		typedef typename std::vector<Tp>::const_reference ConstTpRef;
		typedef typename std::vector<Tp>::reference TpRef;
	private:
		std::vector<Tp> m_data;
		IndexArray      m_dims, m_strides;
		
		void calcStrides();
		size_t calcIndex(const IndexArray &indices) const;
	public:
		Volume();
		Volume(const Eigen::Array<size_t, 3, 1> dims, const size_t nt);
		Volume(const IndexArray &dims);
		//Volume(const IndexArray dims);
		Volume(Nifti &img);
		
		void readFrom(Nifti &img);
		void writeTo(Nifti &img);
		
		ConstTpRef operator[](const IndexArray &indices) const;
		TpRef operator[](const IndexArray &indices);
		
		const SeriesTp series(const IndexArray &indices) const;
		SeriesTp series(const IndexArray &indices);
		
		size_t size() const;
		const IndexArray &dims() const;
};

#include "Volume-inl.h"

#endif //VOLUME_VOLUME
