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
		typedef Eigen::Array<size_t, Eigen::Dynamic, 1> IndexArray;
		typedef Eigen::Array<Tp, Eigen::Dynamic, 1> VectorTp;
	private:
		std::vector<Tp> m_data;
		IndexArray      m_dims, m_strides;
		
		void calcStrides();

	public:
		Volume();
		Volume(const std::vector<size_t> dims);
		Volume(Nifti &img);
		
		void readFrom(Nifti &img);
		
		Tp at(const size_t index) const;
		Tp at(const std::vector<size_t> indices) const;
		
		VectorTp series(const size_t index) const;
		VectorTp series(const std::vector<size_t> indices) const;
		
		size_t size() const;
		const IndexArray &dims() const;
};

#include "Volume-inl.h"

#endif //VOLUME_VOLUME
