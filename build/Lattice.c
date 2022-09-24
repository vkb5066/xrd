#include "Decs.h"
#include "Constants.h"




//Changes direct coordinates of a site s.t. 0 <= site coord i < 1
void MoveToCell(site** s){
	for(uint i = 0; i < DIM; ++i){
		while((*s)->crdsD[i] < 0.0) {
			(*s)->crdsD[i] += 1.0;
		}
		while((*s)->crdsD[i] >= 1.0) {
			(*s)->crdsD[i] -= 1.0;
		}
	}

	return;
}

//Computes the reciprocal lattice of A, where A is a double array of
// {{a0}, {a1}, {a2}}, and puts them in B = {{b0}, {b1}, {b2}}
void RecipLattice(double B[DIM][DIM], const double A[DIM][DIM]){

	///real-space volume is a1 dot (a2 cross a3)
	double V;
	Crs_vv(WORK_ARR_D, A[1], A[2]); 
	Dot_vv(&V, A[0], WORK_ARR_D);

	double scale = 2.0*PI/V;
	Crs_vv(WORK_ARR_D + 0u*DIM, A[1], A[2]);
	Mul_dv(WORK_ARR_D + 0u*DIM, scale, WORK_ARR_D + 0u*DIM);
	Crs_vv(WORK_ARR_D + 1u*DIM, A[2], A[0]);
	Mul_dv(WORK_ARR_D + 1u*DIM, scale, WORK_ARR_D + 1u*DIM);
	Crs_vv(WORK_ARR_D + 2u*DIM, A[0], A[1]);
	Mul_dv(WORK_ARR_D + 2u*DIM, scale, WORK_ARR_D + 2u*DIM);

	for(uint i = 0; i < DIM; ++i){
		for(uint j = 0; j < DIM; ++j){
			B[j][i] = WORK_ARR_D[j*DIM + i];
		}
	}

	return;
}

//Computes G = h*b1 + k*b2 + l*b3 and stores it in res 
void SetG(double res[DIM], const double B[DIM][DIM],
		  const int h, const int k, const int l){

	Mul_iv(WORK_ARR_D + 0u*DIM, h, B[0]);
	Mul_iv(WORK_ARR_D + 1u*DIM, k, B[1]);
	Mul_iv(WORK_ARR_D + 2u*DIM, l, B[2]);

	for(uint i = 0; i < DIM; ++i){
		res[i] = WORK_ARR_D[0u*DIM + i] + ///b1 
				 WORK_ARR_D[1u*DIM + i] + ///b2
			     WORK_ARR_D[2u*DIM + i];  ///b3
	}
	
	return;
}
