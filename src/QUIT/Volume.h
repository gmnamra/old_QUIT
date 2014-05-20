//
//  MultiArray.h
//  MultiArray
//
//  Created by Tobias Wood on 23/04/2014.
//  Copyright (c) 2014 Tobias Wood. All rights reserved.
//

#ifndef VOLUME_H
#define VOLUME_H

#include "Nifti/Nifti.h"
#include "QUIT/MultiArray.h"

enum class Extrap { Zero, Clamp, Error };
enum class Interp { NN };

template<typename Tp, size_t rank>
class VolumeBase {
	public:
		typedef MultiArray<Tp, rank> StorageTp;
		typedef typename StorageTp::const_reference const_reference;
		typedef typename StorageTp::reference reference;
		typedef typename StorageTp::Index Index;
		typedef typename StorageTp::SmallIndex SmallIndex;

	private:
		StorageTp m_storage;
		Eigen::Affine3f m_xform;
		Extrap m_extrap;

	public:
		VolumeBase();
		VolumeBase(const Index &dims, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity(), const Extrap ext = Extrap::Zero);
		VolumeBase(const SmallIndex &dims, const size_t finalDim, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity(), const Extrap ext = Extrap::Zero);
		VolumeBase(const MultiArray<Tp, rank> &a, const Eigen::Affine3f &xform = Eigen::Affine3f::Identity(), const Extrap ext = Extrap::Zero);
		VolumeBase(Nifti &img);

		void readFrom(Nifti &img);
		void writeTo(Nifti &img);

		Index dims() const;
		size_t size() const;
		Eigen::Affine3f transform() const;

		const Extrap &extrapolate() const;
		void setExtrapolate(const Extrap ext);

		const_reference operator[](const size_t i) const;
		const_reference operator[](const Index &vox) const;
		//const_reference operator[](const Eigen::Vector3f &point) const;
		reference operator[](const size_t i);
		reference operator[](const Index &vox);
		//reference operator[](const Eigen::Vector3f &point);
		StorageTp &data() { return m_storage; }
};

template<typename Tp> using Volume = VolumeBase<Tp, 3>;
template<typename Tp> using Series = VolumeBase<Tp, 4>;

#include "Volume-inl.h"

#endif
