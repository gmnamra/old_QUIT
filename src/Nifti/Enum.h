#ifndef ENUM_H
#define ENUM_H

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


#endif // ENUM_H
