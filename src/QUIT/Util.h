#ifndef UTIL_H
#define UTIL_H

#include "Agilent/procpar.h"
#include "Nifti/Nifti.h"
#include "Nifti/ExtensionCodes.h"

bool ReadPP(const Nifti &nii, Agilent::ProcPar &pp);

#endif // UTIL_H
