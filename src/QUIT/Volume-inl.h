/*
 *  Volume-inl.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_VOLUMEINL
#define QUIT_VOLUMEINL

template<typename Tp, size_t rank> VolumeBase<Tp, rank>::VolumeBase() { }

template<typename Tp, size_t rank> VolumeBase<Tp, rank>::VolumeBase(const Index &dims, const Eigen::Affine3f &xform, const Extrap ext) :
	m_storage(dims), m_xform(xform), m_extrap(ext)
{ }

template<typename Tp, size_t rank> VolumeBase<Tp, rank>::VolumeBase(const SmallIndex &dims, const size_t finalDim, const Eigen::Affine3f &xform, const Extrap ext) :
	m_storage(dims, finalDim), m_xform(xform), m_extrap(ext)
{ }

template<typename Tp, size_t rank> VolumeBase<Tp, rank>::VolumeBase(const MultiArray<Tp, rank> &a, const Eigen::Affine3f &xform, const Extrap ext) :
	m_storage(a), m_xform(xform), m_extrap(ext)
{ }

template<typename Tp, size_t rank> VolumeBase<Tp, rank>::VolumeBase(Nifti::Nifti1 &img) :
	m_storage(img.dims().head(rank)),
	m_xform(img.transform()),
	m_extrap(Extrap::Zero)
{
	readFrom(img);
}

template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::dims() const -> Index { return m_storage.dims(); }
template<typename Tp, size_t rank> size_t VolumeBase<Tp, rank>::size() const        { return m_storage.size(); }
template<typename Tp, size_t rank> Eigen::Affine3f VolumeBase<Tp, rank>::transform() const { return m_xform; }
template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::extrapolate() const -> const Extrap & { return m_extrap; }
template<typename Tp, size_t rank> void VolumeBase<Tp, rank>::setExtrapolate(const Extrap ext)      { m_extrap = ext; }

template<typename Tp, size_t rank> void VolumeBase<Tp, rank>::readFrom(Nifti::Nifti1 &img) { img.readVolumes(0, img.dim(4), m_storage.begin(), m_storage.end()); }
template<typename Tp, size_t rank> void VolumeBase<Tp, rank>::writeTo(Nifti::Nifti1 &img) { img.writeVolumes(0, img.dim(4), m_storage.begin(), m_storage.end()); }

template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::operator[](const Index &inVox) const -> const_reference {
	static Tp zero;
	Index vox(inVox);
	for (size_t i = 0; i < rank; i++) {
		if (vox[i] < 0) {
			switch (m_extrap) {
				case Extrap::Zero : return zero; break;
				case Extrap::Clamp : vox[i] = 0; break;
				case Extrap::Error : throw(std::out_of_range("Tried to access outside a volume.")); break;
			}
		} else if (vox[i] >= m_storage.dims()[i]) {
			switch (m_extrap) {
				case Extrap::Zero : return zero; break;
				case Extrap::Clamp : vox[i] = m_storage.dims()[i] - 1; break;
				case Extrap::Error : throw(std::out_of_range("Tried to access outside a volume.")); break;
			}
		}
	}
	return m_storage[vox];
}

template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::operator[](const Index &vox) -> reference {
	return const_cast<reference>(static_cast<const VolumeBase<Tp, rank> &>(*this).operator[](vox));
}

template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::operator[](const size_t i) const -> const_reference { return m_storage[i]; }
template<typename Tp, size_t rank> auto VolumeBase<Tp, rank>::operator[](const size_t i) -> reference {
	return const_cast<reference>(static_cast<const VolumeBase<Tp, rank> &>(*this).operator[](i));
}


#endif // VOLUMEINL_H
