/*
 *  VolumeFunctions.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_VOLUMEFUNCTIONS_H
#define QUIT_VOLUMEFUNCTIONS_H

#include "Volume.h"

namespace QUIT {

enum class DiffType { Forward, Backward, Central };

template<typename Tp>
void VolumeDerivative(const Volume<Tp> &vol, Volume<Tp> g, Volume<Eigen::Vector3f> d, const DiffType &type = DiffType::Forward) {
	typedef typename Volume<Tp>::Index Index;
	Index dx{1, 0, 0}, dy{0, 1, 0}, dz{0, 0, 1};
	for (size_t k = 0; k < vol.dims()[2]; k++) {
		for (size_t j = 0; j < vol.dims()[1]; j++) {
			for (size_t i = 0; i < vol.dims()[0]; i++) {
				Index ijk{i, j, k};
				switch (type) {
					case DiffType::Forward : {
						d[ijk] = Eigen::Vector3f{(vol[ijk + dx] - vol[ijk]),
						                         (vol[ijk + dy] - vol[ijk]),
						                         (vol[ijk + dz] - vol[ijk])};
					} break;
					case DiffType::Backward : {
						d[ijk][0] = (vol[ijk] - vol[ijk - dx]);
						d[ijk][1] = (vol[ijk] - vol[ijk - dy]);
						d[ijk][2] = (vol[ijk] - vol[ijk - dz]);
					} break;
					case DiffType::Central : {
						d[ijk][0] = (vol[ijk + dx] - vol[ijk - dx]) / 2.;
						d[ijk][1] = (vol[ijk + dy] - vol[ijk - dy]) / 2.;
						d[ijk][2] = (vol[ijk + dz] - vol[ijk - dz]) / 2.;

					} break;
				}
				g[ijk] = d[ijk].norm();
			}
		}
	}
}

} // End namespace QUIT

#endif // QUIT_VOLUMEFUNCTIONS_H
