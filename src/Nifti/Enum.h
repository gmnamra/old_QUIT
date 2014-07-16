/** \file Header.h
 \brief Enums to replace NIFTI_... definitions
 - Written by Tobias Wood, IoP KCL
 - This code is released to the public domain. Do with it what you will.
 */
#ifndef LIBNIFTI_ENUM_H
#define LIBNIFTI_ENUM_H

namespace Nifti {

enum class Intent {
	None, Correlation, TTest, FTest, ZScore, ChiSquared, Beta, Binomial,
	Gamma, Poisson, Normal, FTestNonCentral, ChiSquaredNonCentral,
	Logistic, Laplace, Uniform, TTestNonCentral, Weibull, Chi,
	InverseGuassian, ExtremeValue, PValue, LogPValue, Log10PValue,
	Estimate, Label, Neuroname, MatrixGeneral, MatrixSymmetric,
	VectorDisplacement, Vector, Pointset, Triangle, Quaternion,
	Dimensionless,
	Timeseries, NodeIndex, RGBVector, RGBAVector, Shape
};

}

#endif // ENUM_H
