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
#include <string>

#include "Nifti/Enum.h"

namespace Nifti {

// Include these here to prevent namespace pollution

#include "nifti1.h" // NIfTI-1 header specification
#include "nifti2.h" // NIfIT-2 header specification
#include "nifti_analyze.h" // NIfTI version of the ANALYZE 7.5 header

void SwapBytes(size_t n, size_t siz, void *ar);
void SwapNiftiHeader(struct nifti_1_header *h);
void SwapAnalyzeHeader(nifti_analyze75 *h);

inline float FixFloat(const float f); //!< Converts invalid floats to 0 to ensure a marginally sane header

#include "Nifti/Internal-inl.h"

int CodeForIntent(const Intent i);
Intent IntentForCode(const int c);
std::string IntentName(Intent i);

} // End namespace Nifti

#endif // NIFTI_INTERNAL
