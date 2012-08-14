/*
 *  mathsOptimisers.c
 *
 *  Created by Tobias Wood on 01/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include "mathsOptimisers.h"

extern int MATHS_DEBUG;
int LEV_DEBUG = 0;
//******************************************************************************
#pragma mark Residuals
//******************************************************************************
// Single evaulation function residuals
double calcResiduals(gsl_vector *params, void *consts,
					 gsl_vector *dataX, gsl_vector *dataY,
					 eval_type  *function, gsl_vector *residuals)
{
	bool alloc = false;
	gsl_vector *temp;
	if (!residuals)
	{
		alloc = true;
		temp = gsl_vector_alloc(dataX->size);
	}
	else
		temp = residuals;
	function(dataX, params, consts, temp);
	gsl_vector_sub(temp, dataY);
	gsl_vector_mul(temp, temp);
	double sumsq = vector_sum(temp);
	if (alloc)
		gsl_vector_free(temp);
	return sqrt(sumsq);
}

// Evaluate multiple array function residuals
double calcMResiduals(gsl_vector *params, size_t nF, void **consts,
                      gsl_vector **dataX, gsl_vector **dataY,
					  eval_type  **funcs, gsl_vector **residuals)
{
	bool alloc = false;
	gsl_vector **temp;
	if (!residuals)
	{
		alloc = true;
		temp = malloc(nF * sizeof(gsl_vector *));
		for (int i = 0; i < nF; i++)
			temp[i] = gsl_vector_alloc(dataX[i]->size);
	}
	else
		temp = residuals;
	
	double sumsq = 0;
	for (int m = 0; m < nF; m++)
	{
		funcs[m](dataX[m], params, consts[m], temp[m]);
		gsl_vector_sub(temp[m], dataY[m]);
		gsl_vector_mul(temp[m], temp[m]);
		sumsq += vector_sum(temp[m]);
		if (alloc)
			gsl_vector_free(temp[m]);
	}
	if (alloc)
		free(temp);
	return sqrt(sumsq);
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
int goldenSection(gsl_vector *parameters, int P,
                  double loP, double hiP, void *constants,
                  gsl_vector *dataX, gsl_vector *dataY,
				  eval_type *function, double *finalResidue)
{
	// Golden Section Search to find B1	
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	double R = 0.61803399; // Golden ratio - 1
	double C = 1 - R;
	double precision = 0.001;	

	double res1, res2;
	gsl_vector *residuals = gsl_vector_alloc(dataX->size);
	double P0 = loP, P1, P2, P3 = hiP;
	
	int iterations = 0;
	gsl_vector_set(parameters, P, P0);
	res1 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);
	gsl_vector_set(parameters, P, P3);
	res2 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);

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

	gsl_vector_set(parameters, P, P1);
	res1 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);
	gsl_vector_set(parameters, P, P2);
	res2 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);
	
	while ( fabs(P3 - P0) > precision * (fabs(P1) + fabs(P2)))
	{
		iterations++;
		if (res2 < res1)
		{
			P0 = P1; P1 = P2;
			P2 = R * P1 + C * P3;
			res1 = res2;
			gsl_vector_set(parameters, P, P2);
			res2 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);
		}
		else
		{
			P3 = P2; P2 = P1;
			P1 = R * P2 + C * P0;
			res2 = res1;
			gsl_vector_set(parameters, P, P1);
			res1 = calcResiduals(parameters, constants, dataX, dataY, function, residuals);
		}
	}
	
	// Best value for B1
	if (res1 < res2)
	{
		*finalResidue = res1;
		gsl_vector_set(parameters, P, P1);
	}
	else
	{
		*finalResidue = res2;
		gsl_vector_set(parameters, P, P2);
	}
	return iterations;
}

//******************************************************************************
// Non-linear least squares fitting using Levenberg-Marquardt
//******************************************************************************
int levenbergMarquardt(gsl_vector *parameters, size_t nF, void **constants,
                       gsl_vector **dataX, gsl_vector **dataY,
					   eval_type  **funcs, jacob_type **jacFuncs, double *finalResidue)
{
	// Set up variables
	// nP is number of parameters, nD is number of data points
	size_t evaluations = 0, MAX_EVALUATIONS = 100;
	size_t totalD = 0;
	size_t nP = parameters->size;
	if (nF == 0)
	{
		fprintf(stderr, "Invalid parameters for Levenberg-Marquardt, must have at least 1 evaluation function.\n");
		return 0;
	}
	
	for (size_t f = 0; f < nF; f++)
		totalD += dataX[f]->size;
	
	gsl_matrix *Jacob = gsl_matrix_alloc(totalD, nP);
	gsl_matrix *Jacobt = gsl_matrix_alloc(nP, totalD);
	gsl_matrix *JtJ = gsl_matrix_alloc(nP, nP);
	gsl_matrix *JtJAug = gsl_matrix_alloc(nP, nP);
	gsl_vector *allRes = gsl_vector_alloc(totalD);;
	gsl_vector **residuals = malloc(nF * sizeof(gsl_vector *));
	for (size_t f = 0; f < nF; f++)
		residuals[f] = gsl_vector_alloc(dataX[f]->size);
	gsl_vector *JtRes = gsl_vector_alloc(nP);
	gsl_vector *delta = gsl_vector_alloc(nP);
	gsl_vector *tolerance = gsl_vector_alloc(nP);
	gsl_vector *newPars = gsl_vector_alloc(nP);
	
	double oldSum = 0., lambda = 1., tol = sqrt(DBL_EPSILON); // Starting values
	double sumResidues = calcMResiduals(parameters, nF, constants, dataX, dataY, funcs, residuals);
	evaluations++;
	if (LEV_DEBUG > 0)
		fprintf(stdout, "\nStarting Levenberg-Marquardt. Initial residue = %g\n", sumResidues);
	bool outer = true;
	while (outer && (evaluations < MAX_EVALUATIONS))
	{	// Outer loop (Levenberg)
		if (fabs(sumResidues - oldSum) < (tol * oldSum))
		{
			if (LEV_DEBUG > 0)
				fprintf(stdout, "Improvement in residue (%f) is below tolerance (%f). Terminating.\n",
				        fabs(sumResidues - oldSum), (tol * oldSum));
			break;
		}
		if (LEV_DEBUG > 0)
			fprintf(stdout, "\nOuter Loop. Evaluations = %zu, lambda = %f, residual = %g\n", evaluations, lambda, sumResidues);
		
		oldSum = sumResidues;
		// Calculate Jacobians
		size_t jacOffset = 0;
		for (size_t f = 0; f < nF; f++)
		{
			size_t nX = dataX[f]->size;
			gsl_vector_view temp = gsl_vector_subvector(allRes, jacOffset, nX);
			gsl_vector_memcpy(&temp.vector, residuals[f]);
			if (jacFuncs) // We have analytical gradients
				jacFuncs[f](dataX[f], parameters, constants[f], Jacob + (jacOffset * nP));
			else // Numerical gradients
			{
				gsl_vector *tempP = gsl_vector_alloc(nP);
				gsl_vector_memcpy(tempP, parameters);
				gsl_vector *Y = gsl_vector_alloc(nX);
				gsl_vector *dYdX = gsl_vector_alloc(nX);
				funcs[f](dataX[f], parameters, constants[f], Y);
				for (size_t p = 0; p < nP; p++)
				{
					double dP = sqrt(DBL_EPSILON) * gsl_vector_get(tempP, p);
					if (dP < sqrt(DBL_EPSILON))
						dP = sqrt(DBL_EPSILON);
					tempP->data[p * tempP->stride] += dP;
					funcs[f](dataX[f], tempP, constants[f], dYdX); // Calc Y + dY
					gsl_vector_sub(dYdX, Y);                       // Sub to get dY
					gsl_vector_scale(dYdX, 1. / dP);               // Div to get dYdX
					gsl_matrix_view Jview = gsl_matrix_submatrix(Jacob, jacOffset, 0, nX, nP);
					gsl_matrix_set_col(&Jview.matrix, p, dYdX);
					// Restore original value
					tempP->data[p*tempP->stride] = gsl_vector_get(parameters, p); 
				}
				gsl_vector_free(tempP);
				gsl_vector_free(Y);
				gsl_vector_free(dYdX);
			}
			jacOffset += nX;
		}
		// Calculate J'J and J'res
		gsl_matrix_transpose_memcpy(Jacobt, Jacob);
		matrix_mulv(JtRes, Jacobt, allRes);
		matrix_mult(JtJ, Jacobt, Jacob);
		
		gsl_vector_memcpy(tolerance, parameters);
		gsl_vector_scale(tolerance, tol);
		
		if (LEV_DEBUG > 0)
		{
			fprintf(stdout, "Outer Loop Calculations Finished.\n");
			MAT_PRINT(Jacob);
			MAT_PRINT(Jacobt);
			MAT_PRINT(JtJ);
			VEC_PRINT(tolerance);
			fprintf(stdout, "Starting inner loop.\n\n");
		}
		bool inner = true;
		while(inner && (evaluations < MAX_EVALUATIONS))
		{	// Inner loop (Marquardt)
			// Augment with lambda * diag(JtJ)
			gsl_matrix_memcpy(JtJAug, JtJ);
			for (int i = 0; i < nP; i++)
				JtJAug->data[i * JtJAug->tda + i] += lambda * JtJ->data[i * JtJ->tda + i];
			gsl_vector_memcpy(delta, JtRes);
			matrix_solvev(JtJAug, delta);
			// New parameter estimate
			if (!isfinite(vector_norm(delta)))
			{
				if (LEV_DEBUG > 0)
					fprintf(stdout, "Change in parameters is too large. Back to outer loop.\n");
				inner = false;
				break; // Avoid a residuals calculation with an infinite parameter
			}			
			gsl_vector_sub(delta, tolerance);
			if (gsl_vector_isneg(delta))
			{
				if (LEV_DEBUG)
					fprintf(stdout, "Change in parameters is small, accept and back to outer loop.\n");
				gsl_vector_add(parameters, delta);
				lambda /= 10.;
				inner = false;
				break;
			}
			sumResidues = calcMResiduals(parameters, nF, constants, dataX, dataY, funcs, residuals);
			evaluations++;
			if(sumResidues < oldSum)
			{
				if (LEV_DEBUG)
					fprintf(stdout, "Estimate has improved. Back to outer loop.\n");
				lambda /= 10.;
				inner = false; 
			}
			else
			{
				lambda *= 100.;
				if (!isfinite(lambda))
				{
					if (LEV_DEBUG)
						fprintf(stdout, "Lambda is infinite. No further improvement possible. Terminating.\n");
					inner = false;
					outer = false;
				}
			}
		}
	}
	if (LEV_DEBUG > 0)
	{
		fprintf(stdout, "Finished Levenberg-Marquardt. Final evaluations = %zu, residual = %g\n", evaluations, sumResidues);
		VEC_PRINT(parameters);
	}
	gsl_matrix_free(Jacob);
	gsl_matrix_free(Jacobt);
	gsl_matrix_free(JtJ);
	gsl_matrix_free(JtJAug);
	gsl_vector_free(JtRes);
	gsl_vector_free(delta);
	gsl_vector_free(newPars);
	gsl_vector_free(tolerance);
	gsl_vector_free(allRes);
	for (size_t f = 0; f < nF; f++)
		gsl_vector_free(residuals[f]);
	free(residuals);
	*finalResidue = sqrt(sumResidues);
	return evaluations;
}

//******************************************************************************
// Basic Simplex/Amoeba Search Optimiser (Nelder-Mead)
//******************************************************************************
int MATHS_SIMPLEX_DEBUG = 0;
int simplex(gsl_vector *params, size_t nF, void **consts,
		    gsl_vector **dataX, gsl_vector **dataY,
		    eval_type  **funcs, gsl_vector **initial, double *finalResidue)
{
	size_t nP = params->size, nV = params->size + 1, MAX_EVALUATIONS = 250, evaluations = 0;
	gsl_vector *centroid = gsl_vector_alloc(nP),
			   *reflection = gsl_vector_alloc(nP),
			   *expansion = gsl_vector_alloc(nP),
			   *contraction = gsl_vector_alloc(nP),
			   *residuals = gsl_vector_alloc(nV),
	           *simplex[nV];
	gsl_permutation *sort_perm = gsl_permutation_alloc(nV);
	for (size_t i = 0; i < nV; i++)
		simplex[i] = gsl_vector_alloc(nP);
	
	double alpha = 1., gamma = 2.5, rho = 0.5, sigma = 0.25;
	double tol = sqrt(DBL_EPSILON), resDiff = 0.;

	if (initial)
	{
		for (int v = 0; v < nV; v++)
			gsl_vector_memcpy(simplex[v], initial[v]);
	}
	else
	{	// Randomly peturb starting params to produce initial simplex
		//srand(clock());
		for (int v = 0; v < nV; v++)
		{
			gsl_vector_memcpy(simplex[v], params);
			for (int p = 0; p < nP; p++)
			{
				double r = uniform(0.8, 1.2);
				simplex[v]->data[p * simplex[v]->stride] *= r;
			}
		}
	}
	if (MATHS_SIMPLEX_DEBUG > 0)
	{
		//fprintf(stdout, "Start params: "); ARR_D( params, nP); fprintf(stdout, "\n");
	}
	
	for (size_t v = 0; v < nV; v++)
		gsl_vector_set(residuals, v, calcMResiduals(simplex[v], nF, consts, dataX, dataY, funcs, NULL));
	
	
	gsl_sort_vector_index(sort_perm, residuals);
	resDiff = fabs(gsl_vector_get(residuals, gsl_permutation_get(sort_perm, 0)) -
	               gsl_vector_get(residuals, gsl_permutation_get(sort_perm, nV - 1)));

	while ((evaluations < MAX_EVALUATIONS) && (resDiff > tol))
	{
		if (MATHS_SIMPLEX_DEBUG > 0)
		{
			fprintf(stdout, "Iteration %ld. Simplex: \n", evaluations);
			for (size_t v = 0; v < nV; v++)
			{
				//ARR_D( simplex[v], nP );
				fprintf(stdout, "Res: %f\n", gsl_vector_get(residuals, v));
			}
			fprintf(stdout, "Res Diff: %f\n", resDiff);
		}

		// Ignore worst vertex, calculate centroid
		gsl_vector_set_zero(centroid);
		for (size_t v = 0; v < (nV - 1); v++)
			gsl_vector_add(centroid, simplex[v]);
		gsl_vector_scale(centroid, 1. / nP);
		
		// Calculate reflection and compare
		gsl_vector *best = simplex[gsl_permutation_get(sort_perm, 0)];
		gsl_vector *worst = simplex[gsl_permutation_get(sort_perm, nV - 1)];
		
		vector_add_scale(reflection, centroid, 1. + alpha, worst, -alpha);
		double resRef = calcMResiduals(reflection, nF, consts, dataX, dataY, funcs, NULL);
		if ((resRef < gsl_vector_get(residuals, gsl_permutation_get(sort_perm, nV - 2))) &&
			(resRef > gsl_vector_get(residuals, gsl_permutation_get(sort_perm, 0))))
			gsl_vector_memcpy(worst, reflection);
		else if (resRef < gsl_vector_get(residuals, gsl_permutation_get(sort_perm, 0)))
		{	// Reflection is better than any point in simplex
			vector_add_scale(expansion, centroid, 1. + gamma, worst, -gamma);
			double resExp = calcMResiduals(expansion, nF, consts, dataX, dataY, funcs, NULL);
			if (resExp < resRef)
				gsl_vector_memcpy(worst, expansion);
			else 
				gsl_vector_memcpy(worst, reflection);
		}
		else
		{	// Reflection not better than 2nd worst point in simplex
			vector_add_scale(contraction, centroid, 1. + rho, worst, -rho);
			double resCon = calcMResiduals(contraction, nF, consts, dataX, dataY, funcs, NULL);
			if (resCon < gsl_vector_get(residuals, gsl_permutation_get(sort_perm, nV - 1)))
				gsl_vector_memcpy(worst, contraction);
			else
			{	// Still not better, just shrink the simplex
				gsl_vector *bestCopy = gsl_vector_alloc(nP);
				gsl_vector_memcpy(bestCopy, best);
				for (int v = 0; v < nV; v++)
					vector_add_scale(simplex[v], bestCopy, 1 - sigma, simplex[v], sigma);
				gsl_vector_free(bestCopy);
			}
		}
		
		for (size_t v = 0; v < nV; v++)
			gsl_vector_set(residuals, v, calcMResiduals(simplex[v], nF, consts, dataX, dataY, funcs, NULL));
		gsl_sort_vector_index(sort_perm, residuals);
		resDiff = fabs(gsl_vector_get(residuals, gsl_permutation_get(sort_perm, 0)) -
					   gsl_vector_get(residuals, gsl_permutation_get(sort_perm, nV - 1)));
		evaluations++;		
	}
	gsl_vector_memcpy(params, simplex[gsl_permutation_get(sort_perm, 0)]);
	*finalResidue = gsl_vector_get(residuals, gsl_permutation_get(sort_perm, 0));
	return evaluations;
}

//******************************************************************************
// Stochastic Region Contraction
// nS = number of points to sample in parameter space at each step
// nR = number of best points to retain and use for contracting the bounds
//******************************************************************************
void regionContraction(gsl_vector *params, size_t nF, void **consts,
					   gsl_vector **dataX, gsl_vector **dataY,
					   eval_type  **funcs, gsl_vector *initBounds[2],
					   bool **constrained, size_t nS, size_t nR,
					   size_t maxContractions, double thresh, double expand,
					   double *finalResidue)
{
	size_t nP = params->size;
	gsl_vector **samples = malloc(nS * sizeof(gsl_vector *));
	for (int s = 0; s < nS; s++)
		samples[s] = gsl_vector_alloc(nP);
	gsl_vector *sampleRes = gsl_vector_alloc(nS);
	gsl_vector *bounds[2];
	bounds[0] = gsl_vector_alloc(nP);
	bounds[1] = gsl_vector_alloc(nP);
	gsl_vector *regionSize = gsl_vector_alloc(nP);
	gsl_vector_memcpy(bounds[0], initBounds[0]);
	gsl_vector_memcpy(bounds[1], initBounds[1]);
	gsl_permutation *retained = gsl_permutation_alloc(nR);
	size_t c;
	srand(time(NULL));
	for (c = 0; c < maxContractions; c++)
	{
		for (int s = 0; s < nS; s++)
		{
			for (int p = 0; p < nP; p++)
			{
				double rval = uniform(gsl_vector_get(bounds[0], p),
							           gsl_vector_get(bounds[1], p));
				gsl_vector_set(samples[s], p, rval);
				               
			}
			gsl_vector_set(sampleRes, s,
			               calcMResiduals(samples[s], nF, consts, dataX, dataY, funcs, NULL));
		}
		gsl_sort_vector_smallest_index(retained->data, nR, sampleRes);
		// Find the min and max for each parameter in the top nR samples
		gsl_vector_set_all(bounds[0], DBL_MAX);
		gsl_vector_set_all(bounds[1],-DBL_MAX);
		
		for (int r = 0; r < nR; r++)
		{
			for (int p = 0; p < nP; p++)
			{
				double pval = gsl_vector_get(samples[gsl_permutation_get(retained, r)], p);
				if (pval < gsl_vector_get(bounds[0], p))
					gsl_vector_set(bounds[0], p, pval);

				if (pval > gsl_vector_get(bounds[1], p))
					gsl_vector_set(bounds[1], p, pval);
			}
		}
		
		// Terminate if ALL the distances between bounds are under the threshold
		vector_add_scale(regionSize, bounds[1], 1., bounds[0], -1.);
		gsl_vector_div(regionSize, bounds[1]);
		double max_size = gsl_vector_max(regionSize);
		if (max_size < thresh)
			break;
		
		// Expand the boundaries back out in case we just missed a minima,
		// but don't go past initial boundaries if constrained
		vector_add_scale(regionSize, bounds[1], expand, bounds[0], expand);
		gsl_vector_sub(bounds[0], regionSize);
		gsl_vector_add(bounds[1], regionSize);
		for (int p = 0; p < nP; p++)
		{
			if (constrained[0][p] && (gsl_vector_get(bounds[0], p) < gsl_vector_get(initBounds[0], p)))
				gsl_vector_set(bounds[0], p, gsl_vector_get(initBounds[0], p));
			
			if (constrained[1][p] && (gsl_vector_get(bounds[0], p) > gsl_vector_get(initBounds[1], p)))
				gsl_vector_set(bounds[1], p, gsl_vector_get(initBounds[1], p));
		}
		
	}
	// Return the best evaluated solution so far
	gsl_vector_memcpy(params, samples[gsl_permutation_get(retained, 0)]);
	*finalResidue = gsl_vector_get(sampleRes, gsl_permutation_get(retained, 0));
	
	for (int s = 0; s < nS; s++)
		gsl_vector_free(samples[s]);
	free(samples);

	gsl_vector_free(sampleRes);
	gsl_vector_free(bounds[0]);
	gsl_vector_free(bounds[1]);
	gsl_vector_free(regionSize);
	gsl_permutation_free(retained);

}