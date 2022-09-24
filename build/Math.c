#include <stdlib.h>
#include <math.h>

#include "Constants.h"
#include "Decs.h"

//3D vector math, limited complex value support, etc.

//Vector Math------------------------------------------------------------------
//Accum the result of a + b into res
void Add_vv(double res[DIM], const double a[DIM], const double b[DIM]){
	for(uint i = 0; i < DIM; ++i){
		res[i] = a[i] + b[i];
	}

	return;
}

//Accum the result of scalar multiplication of scalar s into vector b
void Mul_iv(double res[DIM], const int a, const double b[DIM]){
	for(uint i = 0; i < DIM; ++i){
		res[i] = (double)a * b[i];
	}

	return;
}

//Accum the result of scalar multiplication of scalar s into vector b
void Mul_dv(double res[DIM], const double a, const double b[DIM]){
	for(uint i = 0; i < DIM; ++i){
		res[i] = a * b[i];
	}

	return;
}

//Accum the result of a dot b into res
void Dot_vv(double* res, const double a[DIM], const double b[DIM]){
	*res = 0.0;
	for(uint i = 0; i < DIM; ++i){
		*res += a[i] * b[i];
	}

	return;
}

//Accum the result of a cross b into res
void Crs_vv(double res[DIM], const double a[DIM], const double b[DIM]){
	res[0] = +(a[1]*b[2] - a[2]*b[1]);
	res[1] = -(a[0]*b[2] - a[2]*b[0]);
	res[2] = +(a[0]*b[1] - a[1]*b[0]);

	return;
}


//Complex Math-----------------------------------------------------------------
//Accum the result of a + b into res
void Add_c(cpxd res, const cpxd a, const cpxd b){
	res.real = a.real + b.real;
	res.imag = a.imag + b.imag;

	return;
}

//Accum the result of ae^(ix) into res, given a real (double)x
//and an already initialized complex number res
void ExpAcc_c(cpxd* res, const double a, const double x){
	res->real += a*cos(x);
	res->imag += a*sin(x);

	return;
}

//Accum the result of a* into res
void Conj_c(cpxd res, const cpxd a){
	res.real = a.real;
	res.imag = -a.imag;

	return;
}

//Accum the square magnitude of a, |a|^2 = a*a, into res
void Mag2_c(double* res, const cpxd a){
	*res = a.real*a.real + a.imag*a.imag;

	return;
}

//Gaussian Smearing------------------------------------------------------------
//Fills arr with a gaussian profile, updates its length
//Based on dx, the distance between two adjacent indices in ANY (not
//(just x) direction.  Cutoff = number of standard deviations before
//the profile is cut (3 to 4 is good)
void GaussianProfile(uint* len, double **arr, 
					 const double sigma, const uint cutoff){
	*len = 2u*(uint)(ceil((double)cutoff*sigma) - 1.0);
	if(*len%2u == 0u){
		(*len)++;
	}

	int mp = ((int)*len - 1)/2; ///the midpoint of the array
	double twoSigSqd = 2.0*sigma*sigma;

	*arr = malloc(*len*sizeof(double));
	for(int i = 0; i < *len; ++i){
		(*arr)[i] = 1.0/sqrt(PI*twoSigSqd) * 
			        exp(-(double)((i-mp)*(i-mp))/twoSigSqd);
	}
}

//Computes the cross-convolution between x and p(rofile), both of len
//len (should be an odd number), at x's midpoint
//Assumes that x and gProf are completly filled with appropriate 
//extrapolations, etc.  This is done in Grid.c
double CConvS(const uint len, const double* x, const double* p){
	double ret = 0.0;
	for(uint i = 0; i < len; ++i){ 
			ret += p[i]*(x[i]);
	}

	return ret;
}


