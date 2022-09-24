#include <stdlib.h>
#include <math.h>

#include "Decs.h"
#include "Constants.h"

//compares two recip lattice point structs based on their Q values
int cmpr(const void* a, const void* b){
	rlpt* a_ = (rlpt*)a;
	rlpt* b_ = (rlpt*)b;
	
	if(a_->Q < b_->Q) return -1;
	if(a_->Q > b_->Q) return +1;
	return 0;
}


//Fills an array of recip lat point structures given the max
//magnitude of h, k, and l, updates it's size
//Also sorts the array by Q value and sets the structure factor to 0
void MakeRecipPointArr(uint*restrict len, rlpt*restrict*restrict arr, 
					   const double B[DIM][DIM],
					   const int hmax, const int kmax, const int lmax,
					   const double wavelength){
	///Make array
	*len = (uint)((2*hmax + 1)*(2*kmax + 1)*(2*lmax + 1));
	*arr = malloc(*len*sizeof(rlpt));

	///Fill array
	uint counter = 0;
	for(int h = -hmax; h <= hmax; ++h){
	for(int k = -kmax; k <= kmax; ++k){
	for(int l = -lmax; l <= lmax; ++l){
		rlpt pt;
		pt.hkl[0] = h; pt.hkl[1] = k; pt.hkl[2] = l;
		SetG(WORK_ARR_D, B, h, k, l);
		Dot_vv(&pt.Q, WORK_ARR_D, WORK_ARR_D);
		pt.Q = sqrt(pt.Q);

#if APPLY_LP
		double arg = pt.Q*wavelength/4.0/PI;
		if(0.0 < arg && arg < 1.0){ ///(should be) no need to check negatives
			double theta = asin(arg);
			double cos2theta = cos(2.0*theta);
			pt.lp = (1.0 + cos2theta*cos2theta) / (arg*arg * cos(theta));
		}
		else pt.lp = -1.0;
#endif

		pt.F.real = 0.0;
		pt.F.imag = 0.0;

		(*arr)[counter] = pt;
		counter++;
	}
	}
	}

	///Organize array
	qsort(*arr, *len, sizeof(rlpt), cmpr);

	return;
}


///Calculates the phi in the scattering factor equation f_j e^(i phi)
///i.e. phi = 2pi (hkl) dot (abc)
///WARNING: IF THIS GETS MULTITHREADED, USE THE BOTTOM FUNCTION
double CalcPhi(const int hkl[DIM], const double abc[DIM]){
#if 1 ///if runtime really becomes a problem, you may try setting this false
	double res;
	for(uint i = 0; i < 3; ++i){
		WORK_ARR_D[i] = (double)hkl[i];
	}
	Dot_vv(&res, WORK_ARR_D, abc);

	return 2.0*PI*res;
#else
	return 2.0*PI*((((double)hkl[0])*abc[0]) + 
				   (((double)hkl[1])*abc[1]) + 
				   (((double)hkl[2])*abc[2]));
#endif
}


//Given xMin, xMax, # of x sample pts, and a (theoretically) continuous 
//function f(x), gives x & y:
//x is the x samples and y are the gaussian smoothed f(x) interpolated to x
//f(x) is given as an array of rlpts: x = rlpt.Q, y = rlpt.F.real
void MakeDiffGraph(const uint nSamPts, double** x, double** y, 
				   const uint nRecipPts, const rlpt* f, 
				   const double xMin, const double xMax,
				   const double sigma, const uint cutoff){

	///Make x, y_
	double dx = (xMax - xMin) / (double)(nSamPts - 1u);
	*x = malloc(nSamPts*sizeof(double));
	double* y_ = calloc(nSamPts, sizeof(double));
	for(uint i = 0; i < nSamPts; ++i){
		(*x)[i] = xMin + (double)i*dx;
	}
	for(uint i = 0; i < nRecipPts; ++i){
		int ind = (int)((f[i].Q - xMin) / dx);
		if(0 <= ind && ind < nSamPts){
			y_[ind] += f[i].F.real;
		}
	}


	///Make gaussian profile
	uint gProfLen; double* gProf;
	GaussianProfile(&gProfLen, &gProf, sigma, cutoff);


	///Smooth things over
	*y = calloc(nSamPts, sizeof(double));
	int lowInd;
	for(uint i = 0; i < nSamPts; ++i){
		////make sure we don't overstep any boundaries
		lowInd = i - ((int)gProfLen - 1)/2;
		if(lowInd < 0) lowInd = 0;
		else if(lowInd >= nSamPts - gProfLen) lowInd = (int)nSamPts - (int)gProfLen; 

		(*y)[i] = CConvS(gProfLen, y_ + lowInd, gProf);
	}

	free(y_);
	free(gProf);

}