#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Constants.h"
#include "Decs.h"

extern int WORK_ARR_I[WARR_SIZE] = {0};
extern double WORK_ARR_D[WARR_SIZE] = {0.0};

int main(int argc, char* argv[]){

	//Read any params from the command line to overwrite defaults
	char* infile; char* outfileInt; char* outfileGraph;
	int hmax; int kmax; int lmax;
	double qMin; double qMax; uint nAFFPts; uint nQPts; double sigma; uint cutoff;
	double wavelength;
	CLAHandler(argc, argv, &infile, &outfileInt, &outfileGraph,
			   &hmax, &kmax, &lmax, 
			   &qMin, &qMax, &nAFFPts, &nQPts, &sigma, &cutoff, &wavelength); 
	double dqAff = (qMax - qMin) / (double)(nAFFPts - 1u);


	//Get the lattice(s)
	double A[DIM][DIM]; double B[DIM][DIM];
	uint nSpecies; uint* species; uint nSites; site** sites;
	ReadLatticeFile(infile, A, &nSpecies, &species, &nSites, &sites);
	RecipLattice(B, A);
	

	//Set up tables: mapping from z values to small array values ...
	uint map[ZMAX]; BuildZMap(nSpecies, species, &map);
	// ... and the matrix of aff values to interp between later
	double** affMatrix;
	BuildAffMatrix(nSpecies, nAFFPts, &affMatrix, species, qMin, dqAff);

	
	//Build an array of reciprocal lattice points
	uint nRlpts; rlpt* rlpts;
	MakeRecipPointArr(&nRlpts, &rlpts, B, hmax, kmax, lmax, wavelength);


	//Initialize interpolation params ...
	double qLo; double qHi;
	double interpParams[ZMAX][2]; /// = {{affLo 0, slope Aff 0}, { ... }}
	for(uint i = 0; i < nSpecies; ++i) UpdateQParams(qMin, qMin, dqAff, 
										affMatrix[i], &qLo, &qHi, 
										&interpParams[species[i]][0], 
										&interpParams[species[i]][1]);
	// ... and then begin main part of the program
	double fj; double phi;
	for(uint i = 0; i < nRlpts; ++i){
		///No need to do anything if we're at Q < qMin
		if(rlpts[i].Q < qMin) continue;

		///Compute F(G) through sum of all unit cell atoms
		for(uint j = 0; j < nSites; ++j){
			fj = SampleAffArr(rlpts[i].Q, qLo, 
							  interpParams[sites[j]->species][0], 
							  interpParams[sites[j]->species][1]);

			phi = CalcPhi(rlpts[i].hkl, sites[j]->crdsD);
			ExpAcc_c(&(rlpts[i].F), fj*sites[j]->occ, phi);
		}

		///Store the |F|^2 into the real part, apply LP correction
		Mag2_c(&(rlpts[i].F.real), rlpts[i].F);
#if APPLY_LP
		rlpts[i].F.real *= rlpts[i].lp;
#endif

		///We only have interpolated AFF values up to qMax - break loop if
		///we're past that point since our array is sorted
		if(rlpts[i].Q > qMax) break;

		///Update interpolation params if necessary
		if(rlpts[i].Q > qHi){
			for(uint j = 0; j < nSpecies; ++j){
				UpdateQParams(rlpts[i].Q, qMin, dqAff, affMatrix[j], &qLo, &qHi, 
							  &(interpParams[species[j]][0]), 
							  &(interpParams[species[j]][1]));
			}
		}
	}


	//Plotting Loops: apply profile function and interpolate Q values to
	//array indices
	WriteInts(outfileInt, nRlpts, rlpts);
	double* x; double* y;
	MakeDiffGraph(nQPts, &x, &y, nRlpts, rlpts, qMin, qMax, sigma, cutoff);
	WriteGraph(outfileGraph, nQPts, x, y);


	//Clean up, end
	free(infile); free(outfileInt); free(outfileGraph);
	free(species);
	free(sites);
	for(uint i = 0; i < nSpecies; ++i) free(affMatrix[i]);
	free(affMatrix);
	free(rlpts);
	free(x); free(y);

	return 0;
}