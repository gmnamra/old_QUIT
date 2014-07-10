//
//  Internal-inl.h
//  NiftiImage
//
//  Created by Tobias Wood on 18/09/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#ifndef NIFTI_INTERNAL_INL
#define NIFTI_INTERNAL_INL

inline float FixFloat(const float f) {
	if (isfinite(f))
		return f;
	else
		return 0.;
}

#endif
