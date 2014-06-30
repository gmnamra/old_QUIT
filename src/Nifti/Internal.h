//
//  Internal.h
//  NiftiImage
//
//  Created by Tobias Wood on 18/09/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#ifndef NIFTI_INTERNAL
#define NIFTI_INTERNAL

#include <cstddef>
#include <cmath>

namespace Nifti {

// Include these here to prevent namespace pollution

#include "nifti1.h" // NIFTI-1 header specification
#include "nifti_analyze.h" // NIFTI version of the ANALYZE 7.5 header

void swapBytes(size_t n, size_t siz, void *ar);
void swapNiftiHeader(struct nifti_1_header *h);
void swapAnalyzeHeader(nifti_analyze75 *h);

inline float fixFloat(const float f); //!< Converts invalid floats to 0 to ensure a marginally sane header

#include "Nifti/Internal-inl.h"

} // End namespace Nifti

#endif // NIFTI_INTERNAL
