/*
 *  mathsOptimisers.c
 *
 *  Created by Tobias Wood on 01/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include "mathsOptimisers.h"
//******************************************************************************
#pragma mark Residuals
//******************************************************************************
// Single evaulation function residuals
inline double calcResiduals(double *parameters, double *constants,
                            double *dataX, double *dataY, size_t nD,
				            eval_type *function, double *residuals, bool norm)
{
	double sum = 0.;
	double funcY[nD];
	for (int d = 0; d < nD; d++)
		funcY[d] = function(dataX[d], parameters, constants);
	if (norm)
		arrayScale(funcY, funcY, 1. /arrayMean(funcY, nD), nD);
	for (int d = 0; d < nD; d++)
	{
		double res = dataY[d] - funcY[d];
		if (residuals)
			residuals[d] = res;
		sum += res*res;
	}
	return sqrt(sum);
}

// One array function residuals
double calcAResiduals(double *params, double *consts,
                      double *dataX, double *dataY, size_t nD,
	 	              eval_array_type *function, double *residuals, bool norm)
{
	double sum = 0.;
	double funcY[nD];
	function(dataX, params, consts, funcY, nD);
	if (norm)
		arrayScale(funcY, funcY, 1. /arrayMean(funcY, nD), nD);
	for (int d = 0; d < nD; d++)
	{
		double res = dataY[d] - funcY[d];
		if (residuals)
			residuals[d] = res;
		sum += res*res;
	}
	return sqrt(sum);
}

// Evaluate multiple function residuals
inline double calcMResiduals(double *params, double **consts, size_t nM,
							 double **dataX, double **dataY, size_t *nD,
							 eval_type **funcs, double **residuals, bool norm)
{
	double sum = 0.;
	for (int m = 0; m < nM; m++)
	{
		double *funcY = malloc(nD[m] * sizeof(double));
		for (int d = 0; d < nD[m]; d++)
			funcY[d] = funcs[m](dataX[m][d], params, consts[m]);
		if (norm)
			arrayScale(funcY, funcY, 1. / arrayMean(funcY, nD[m]), nD[m]);
		double thisRes = 0.;
		for (int d = 0; d < nD[m]; d++)
		{
			double res = dataY[m][d] - funcY[d];
			if (residuals)
				residuals[m][d] = res;
			thisRes += res*res;
		}
		sum += thisRes;
		free(funcY);
	}
	sum = sqrt(sum);
	return sum;
}

// Evaluate multiple array function residuals
inline double calcMAResiduals(double *params, double **consts, size_t nM,
							  double **dataX, double **dataY, size_t *nD,
							  eval_array_type **funcs, double **residuals,
							  bool norm)
{
	double sum = 0.;
	for (int m = 0; m < nM; m++)
	{
		double *funcY = malloc(nD[m] * sizeof(double));
		funcs[m](dataX[m], params, consts[m], funcY, nD[m]);
		
		if (norm)
			arrayScale(funcY, funcY, 1. / arrayMean(funcY, nD[m]), nD[m]);
		double thisRes = 0.;
		for (int d = 0; d < nD[m]; d++)
		{
			double res = dataY[m][d] - funcY[d];
			if (residuals)
				residuals[m][d] = res;
			thisRes += res*res;
		}
		sum += thisRes;
		free(funcY);
	}
	sum = sqrt(sum);
	return sum;
}
//******************************************************************************
#pragma mark Distributions
//******************************************************************************
double uniform(double lo, double hi);
inline double uniform(double lo, double hi)
{
	double range = hi - lo;
	double r = ((1. * rand()) / RAND_MAX);
	return (r * range) + lo;
}

//******************************************************************************
// Basic least squares fitting
//******************************************************************************
void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res)
{
	double sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < nD; i++)
	{
		double x = X[i];
		double y = Y[i];
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	*slope = (nD * sumXY - (sumX * sumY)) / (nD * sumXX - (sumX * sumX));
	*inter = (sumY - (*slope) * sumX) / nD;
	
	if (res)
	{
		*res = 0.;
		double m = *slope; double c = *inter;
		for (int i = 0; i < nD; i++)
			*res += pow(Y[i] - (m*X[i] + c), 2.);
	}
}

//******************************************************************************
#pragma mark Actual Optimisers
//******************************************************************************

//******************************************************************************
// Golden Ratio section search
//******************************************************************************
int goldenSection(double *parameters, int nP, int P,
                  double loP, double hiP, double *constants,
                  double *dataX, double *dataY, int nD,
				  eval_type *function, double *finalResidue)
{
	// Golden Section Search to find B1	
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	double R = 0.61803399; // Golden ratio - 1
	double C = 1 - R;
	double precision = 0.001;	

	double res1, res2, residuals[nD];
	double P0 = loP, P1, P2, P3 = hiP;
	
	int iterations = 0;
	parameters[P] = P0;
	res1 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);
	parameters[P] = P3;
	res2 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);

	if (res1 < res2)
	{
		P1 = P0 + 0.2;
		P2 = P1 + C * (P3 - P1);
	}
	else
	{
		P2 = P3 - 0.2;
		P1 = P2 - C * (P2 - P0);
	}

	parameters[P] = P1;
	res1 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);
	parameters[P] = P2;
	res2 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);
	
	while ( fabs(P3 - P0) > precision * (fabs(P1) + fabs(P2)))
	{
		iterations++;
		if (res2 < res1)
		{
			P0 = P1; P1 = P2;
			P2 = R * P1 + C * P3;
			res1 = res2;
			parameters[P] = P2;
			res2 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);
		}
		else
		{
			P3 = P2; P2 = P1;
			P1 = R * P2 + C * P0;
			res2 = res1;
			parameters[P] = P1;
			res1 = calcResiduals(parameters, constants, dataX, dataY, nD, function, residuals, false);
		}
	}
	
	// Best value for B1
	if (res1 < res2)
	{
		*finalResidue = res1;
		parameters[P] = P1;
	}
	else
	{
		*finalResidue = res2;
		parameters[P] = P2;
	}
	return iterations;
}

//******************************************************************************
// Non-linear least squares fitting using Levenberg-Marquardt
//******************************************************************************
int levenbergMarquardt(double *parameters, size_t nP, double **constants,
                       double **dataX, double **dataY, eval_array_type **funcs,
		               jacob_type **jacFuncs, size_t *nD, size_t nF,
		               double *loBounds, double *hiBounds,
		               bool normalise, double *finalResidue)
{
	// Set up variables
	// nP is number of parameters, nD is number of data points
	size_t evaluations = 0, MAX_EVALUATIONS = 100;
	size_t totalD = 0;
	
	if (nF == 0)
	{
		fprintf(stderr, "Invalid parameters for Levenberg-Marquardt, must have at least 1 evaluation function.\n");
		return 0;
	}
	
	for (size_t f = 0; f < nF; f++)
		totalD += nD[f];
	
	double *Jacob = matrixAlloc(totalD, nP);
	double *Jacobt = matrixAlloc(nP, totalD);
	double *JtJ = matrixAlloc(nP, nP);
	double *JtJAug = matrixAlloc(nP, nP);
	double *residuals[nF], allRes[totalD];
	for (size_t f = 0; f < nF; f++)
		residuals[f] = arrayAlloc(nD[f]);
	double *JtRes = arrayAlloc(nP);
	double *delta = arrayAlloc(nP);
	double *newPars = arrayAlloc(nP);
	
	double oldSum = 0., lambda = 1., tol = sqrt(EPS); // Starting values
	double sumResidues = calcMAResiduals(parameters, constants, nF, dataX, dataY, nD, funcs, residuals, normalise);
	
	if (MATHS_DEBUG > 0)
		fprintf(stdout, "\nStarting Levenberg-Marquardt. Initial residue = %g\n", sumResidues);
	bool outer = true;
	while (outer && (evaluations < MAX_EVALUATIONS) && (fabs(sumResidues - oldSum) > tol))
	{	// Outer loop (Levenberg)
		if (MATHS_DEBUG > 0)
			fprintf(stdout, "\nOuter Loop. Evaluations = %zu, lambda = %f, residual = %g\n", evaluations, lambda, sumResidues);
		oldSum = sumResidues;
		
		// Calculate Jacobians
		size_t jacOffset = 0;
		for (size_t f = 0; f < nF; f++)
		{
			arrayCopy(allRes + jacOffset, residuals[f], nD[f]);
			
			if (jacFuncs) // We have analytical gradients
				jacFuncs[f](dataX[f], nD[f], parameters, constants[f], Jacob + (jacOffset * nP));
			else // Numerical gradients
			{
				double *tempP = arrayAlloc(nP), *valAtP = arrayAlloc(nD[f]),
				 	  *tempVal = arrayAlloc(nD[f]), *grad = arrayAlloc(nD[f]);
				funcs[f](dataX[f], parameters, constants[f], valAtP, nD[f]);
				arrayCopy(tempP, parameters, nP);
				for (size_t p = 0; p < nP; p++)
				{
					double dP = sqrt(EPS) * tempP[p];
					if (dP < sqrt(EPS))
						dP = sqrt(EPS);
					tempP[p] += dP;
					funcs[f](dataX[f], tempP, constants[f], tempVal, nD[f]);
					arraySub(grad, tempVal, valAtP, nD[f]);
					arrayScale(grad, grad, 1. / dP, nD[f]);
					for (size_t d = 0; d < nD[f]; d++)
						Jacob[(jacOffset * nP) + (d * nP) + p] = grad[d];
					tempP[p] = parameters[p]; // Restore original value
				}
				free(tempP); free(valAtP); free(tempVal); free(grad);
			}
			
			jacOffset += nD[f];
		}
		// Calculate J'J and J'res
		matrixTranspose(Jacobt, Jacob, totalD, nP);
		matrixMult(JtRes, Jacobt, allRes, nP, totalD, 1);
		matrixMult(JtJ, Jacobt, Jacob, nP, totalD, nP);
		
		if (MATHS_DEBUG > 1)
		{
			MAT_D(Jacob, totalD, nP);
			MAT_D(Jacobt, nP, totalD);
			ARR_D(allRes, totalD);
			MAT_D(JtRes, nP, 1);			
			MAT_D(JtJ, nP, nP);
		}
				
		double deltaTolerance = tol * matrixInfNorm(parameters, nP, 1);
		bool inner = true;
		while(inner && (evaluations < MAX_EVALUATIONS))
		{	// Inner loop (Marquardt)
			evaluations++;
			// Augment with lambda * diag(JtJ)
			arrayCopy(JtJAug, JtJ, nP * nP);
			for (int i = 0; i < nP; i++)
				JtJAug[i * nP + i] += lambda * JtJ[i * nP + i];
			matrixSolve(delta, JtJAug, JtRes, nP, 1);
			// New parameter estimate
			arrayAdd(newPars, parameters, delta, nP);
			
			// Get in bounds
			for (size_t p = 0; p < nP; p++)
			{
				if (loBounds && (newPars[p] < loBounds[p]))
				{
					double ratio = fabs(parameters[p] - loBounds[p]) / fabs(delta[p]);
					arrayAddScale(newPars, parameters, 1., delta, ratio, nP);
				}
				if (hiBounds && (newPars[p] > hiBounds[p]))
				{
					double ratio = fabs(parameters[p] - hiBounds[p]) / fabs(delta[p]);
					arrayAddScale(newPars, parameters, 1., delta, ratio, nP);
				}
			}
			
			sumResidues = calcMAResiduals(newPars, constants, nF, dataX, dataY, nD, funcs, residuals, normalise);
			if (MATHS_DEBUG > 0)
			{
				fprintf(stdout, "\nInner Loop. Evaluations = %zu, lambda = %f, residual = %g\n", evaluations, lambda, sumResidues);
				ARR_D(delta, nP);
				ARR_D(newPars, nP);
			}
			if (!isfinite(sumResidues))
				inner = false;			
			if ((matrixInfNorm(delta, nP, 1) < deltaTolerance) || // Change in parameters is small, back to outer loop
			    (sumResidues < oldSum))                                // Improvement in estimate, accept
			{
				arrayCopy(parameters, newPars, nP);
				lambda /= 10.;
				inner = false; 
			}
			else
			{
				lambda *= 10.;
				if (!isfinite(lambda))
				{	// No further improvement is possible
					inner = false;
					outer = false;
				}
			}
		}
	}
	if (MATHS_DEBUG > 0)
	{
		fprintf(stdout, "Finished Levenberg-Marquardt. Final evaluations = %zu, residual = %g\n", evaluations, sumResidues);
		ARR_D(parameters, nP);
	}
	for (size_t f = 0; f < nF; f++)
		free(residuals[f]);
	free(Jacob);
	free(Jacobt);
	free(JtJ);
	free(JtJAug);
	free(JtRes);
	free(delta);
	free(newPars);
	*finalResidue = sqrt(sumResidues);
	return evaluations;
}

//******************************************************************************
// Basic Simplex/Amoeba Search Optimiser (Nelder-Mead)
//******************************************************************************
int MATHS_SIMPLEX_DEBUG = 0;
int simplex(double *params, size_t nP, double **consts, size_t nM,
		    double **dataX, double **dataY, size_t *nD,
		    eval_type **funcs, double **initial, double *finalResidue)
{
	size_t nV = nP + 1, best[nV], MAX_EVALUATIONS = 250, evaluations = 0;
	double simplex[nV][nP], centroid[nP], reflection[nP], expansion[nP], contraction[nP], residuals[nV];
	
	double alpha = 1., gamma = 2.5, rho = 0.5, sigma = 0.25;
	double tol = sqrt(EPS), resDiff = 0.;

	if (initial)
	{
		for (int v = 0; v < nV; v++)
			arrayCopy(simplex[v], initial[v], nP);
	}
	else
	{	// Randomly peturb starting params to produce initial simplex
		//srand(clock());
		for (int v = 0; v < nV; v++)
		{
			arrayCopy(simplex[v], params, nP);
			for (int p = 0; p < nP; p++)
			{
				double r = uniform(0.8, 1.2);
				simplex[v][p] *= r;
			}
		}
	}
	if (MATHS_SIMPLEX_DEBUG > 0)
	{
		fprintf(stdout, "Start params: "); ARR_D( params, nP); fprintf(stdout, "\n");
	}
	
	for (size_t v = 0; v < nV; v++)
		residuals[v] = calcMResiduals(simplex[v], consts, nM, dataX, dataY, nD, funcs, NULL, false);
	arrayIndexSort(residuals, best, SORT_ASCEND, nV);
	resDiff = fabs(residuals[best[0]] - residuals[best[nV - 1]]);

	while ((evaluations < MAX_EVALUATIONS) && (resDiff > tol))
	{
		if (MATHS_SIMPLEX_DEBUG > 0)
		{
			fprintf(stdout, "Iteration %ld. Simplex: \n", evaluations);
			for (size_t v = 0; v < nV; v++)
			{
				ARR_D( simplex[v], nP );
				fprintf(stdout, "Res: %f\n", residuals[v]);
			}
			fprintf(stdout, "Res Diff: %f\n", resDiff);
		}

		// Ignore worst vertex, calculate centroid
		arrayZero(centroid, nP);
		for (size_t v = 0; v < (nV - 1); v++)
			arrayAdd(centroid, centroid, simplex[v], nP);
		arrayScale(centroid, centroid, 1. / nP, nP);
		
		// Calculate reflection and compare
		arrayAddScale(reflection, centroid, 1. + alpha, simplex[best[nV -1]], -alpha, nP);
		double resRef = calcMResiduals(reflection, consts, nM, dataX, dataY, nD, funcs, NULL, false);
		if ((resRef < residuals[best[nV - 2]]) && (resRef > residuals[best[0]]))
			arrayCopy(simplex[best[nV - 1]], reflection, nP);
		else if (resRef < residuals[best[0]])
		{	// Reflection is better than any point in simplex
			arrayAddScale(expansion, centroid, 1. + gamma, simplex[best[nV - 1]], -gamma, nP);
			double resExp = calcMResiduals(expansion, consts, nM, dataX, dataY, nD, funcs, NULL, false);
			if (resExp < resRef)
				arrayCopy(simplex[best[nV - 1]], expansion, nP);
			else 
				arrayCopy(simplex[best[nV - 1]], reflection, nP);
		}
		else
		{	// Reflection not better than 2nd worst point in simplex
			arrayAddScale(contraction, centroid, 1. + rho, simplex[best[nV - 1]], -rho, nP);
			double resCon = calcMResiduals(contraction, consts, nM, dataX, dataY, nD, funcs, NULL, false);
			if (resCon < residuals[best[nV - 1]])
				arrayCopy(simplex[best[nV - 1]], contraction, nP);
			else
			{	// Still not better, just shrink the simplex
				double bestCopy[nP]; arrayCopy(bestCopy, simplex[best[0]], nP);
				for (int v = 0; v < nV; v++)
				{
					arrayAddScale(simplex[v], bestCopy, 1 - sigma, simplex[v], sigma, nP);
				}
			}
		}
		
		for (size_t v = 0; v < nV; v++)
			residuals[v] = calcMResiduals(simplex[v], consts, nM, dataX, dataY, nD, funcs, NULL, false);
		arrayIndexSort(residuals, best, SORT_ASCEND, nV);
		resDiff = fabs(residuals[best[0]] - residuals[best[nV - 1]]);
		evaluations++;		
	}
	arrayCopy(params, simplex[best[0]], nP);
	*finalResidue = residuals[best[0]];
	return evaluations;
}

//******************************************************************************
// Stochastic Region Contraction
// nS = number of points to sample in parameter space at each step
// nR = number of best points to retain and use for contracting the bounds
//******************************************************************************
void regionContraction(double *params, size_t nP, double **consts, size_t nM,
					   double **dataX, double **dataY, size_t *nD, bool norm,
					   eval_array_type **funcs, double **initBounds,
					   bool **constrained, size_t nS, size_t nR,
					   size_t maxContractions, double thresh, double expand,
					   double *finalResidue)
{

	double **samples = malloc(nS * sizeof(double *));
	double *sampleRes = malloc(nS * sizeof(double));
	size_t *sortedRes = malloc(nS * sizeof(size_t));
	for (int s = 0; s < nS; s++)
		samples[s] = malloc(nP * sizeof(double));
	double bounds[2][nP];
	arrayCopy(bounds[0], initBounds[0], nP);
	arrayCopy(bounds[1], initBounds[1], nP);
	
	size_t c;
	srand(time(NULL));
	for (c = 0; c < maxContractions; c++)
	{
		for (int s = 0; s < nS; s++)
		{
			for (int p = 0; p < nP; p++)
				samples[s][p] = uniform(bounds[0][p], bounds[1][p]);
			sampleRes[s] = calcMAResiduals(samples[s], consts, nM, dataX, dataY, nD, funcs, NULL, norm);
		}
		arrayIndexSort(sampleRes, sortedRes, SORT_ASCEND, nS);
		
		// Find the min and max for each parameter in the top nR samples
		arraySet(bounds[0],  NUM_MAX, nP);
		arraySet(bounds[1], -NUM_MAX, nP);
		
		for (int r = 0; r < nR; r++)
		{
			for (int p = 0; p < nP; p++)
			{
				if (samples[sortedRes[r]][p] < bounds[0][p])
					bounds[0][p] = samples[sortedRes[r]][p];
				if (samples[sortedRes[r]][p] > bounds[1][p])
					bounds[1][p] = samples[sortedRes[r]][p];
			}
		}
		
		double regionSize[nP];
		arraySub(regionSize, bounds[1], bounds[0], nP);
		// Terminate if ALL the distances between bounds are under the threshold
		bool all = true;
		for (int p = 0; p < nP; p++)
		{
			if ((regionSize[p] / bounds[1][p]) > thresh)
			{
				all = false;
				break;
			}
		}
		if (all)
			break;
		
		// Expand the boundaries back out in case we just missed a minima,
		// but don't go past initial boundaries if constrained
		for (int p = 0; p < nP; p++)
		{
			bounds[0][p] -= regionSize[p] * expand;
			if (constrained[0][p] && (bounds[0][p] < initBounds[0][p]))
				bounds[0][p] = initBounds[0][p];
			
			bounds[1][p] += regionSize[p] * expand;
			if (constrained[1][p] && (bounds[1][p] > initBounds[1][p]))
				bounds[1][p] = initBounds[1][p];
		}
		
	}
	// Return the best evaluated solution so far
	for (int p = 0; p < nP; p++)
		params[p] = samples[sortedRes[0]][p];
	*finalResidue = sampleRes[sortedRes[0]];
	
	free(sortedRes); free(sampleRes);
	for (int s = 0; s < nS; s++)
		free(samples[s]);
	free(samples);
	if (MATHS_DEBUG > 0)
	{
		ARR_D(bounds[0], nP);
		ARR_D(bounds[1], nP);
		fprintf(stdout, "   "); ARR_D(params, nP);
	}
}