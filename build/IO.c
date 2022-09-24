#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "Decs.h"
#include "Constants.h"

//Command line arguments-------------------------------------------------------


void CLAHandler(int argc, char* argv[],
				char** infileName, 
				char** outLisName, char** outGraphName,
				int* absMaxH, int* absMaxK, int* absMaxL,
				double* qMin, double* qMax, 
				uint* nAffSamps, uint* nQSamps, double* sigma, uint* cutoff,
				double* wavelength){

	///Set defaults
	*infileName = malloc(LINESIZE*sizeof(char)); 
	strcpy(*infileName, INFILE_LATTICE_DEFAULT);
	*outLisName = malloc(LINESIZE*sizeof(char));
	strcpy(*outLisName, OUTFILE_INT_DEFAULT);
	*outGraphName = malloc(LINESIZE*sizeof(char));
	strcpy(*outGraphName, OUTFILE_GRAPH_DEFAULT);
	*absMaxH = MAX_MAG_H; *absMaxK = MAX_MAG_K; *absMaxL = MAX_MAG_L;
	*qMin = QMIN_DEFAULT;
	*qMax = QMAX_DEFAULT;
	*nAffSamps = N_AFF_SAMPLES_DEFAULT;
	*nQSamps = N_Q_SAMPLES_DEFAULT;
	*sigma = SIGMA_DEFAULT;
	*cutoff = CUTOFF_DEFAULT;
	*wavelength = WVLNGTH_DEFAULT;

	if(argc < 2) goto fend;

	for(int i = 1; i < argc; ++i){
		///lattice input, intensities & graph output files
		if(strstr(argv[i], "-li")){
			for(uint j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					(*infileName)[j] = '\0';
					break;
				}
				(*infileName)[j] = argv[i + 1][j];
			}
			continue;
		}
		if(strstr(argv[i], "-io")){
			for(uint j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					(*outLisName)[j] = '\0';
					break;
				}
				(*outLisName)[j] = argv[i + 1][j];
			}
			continue;
		}
		if(strstr(argv[i], "-go")){
			for(uint j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					(*outGraphName)[j] = '\0';
					break;
				}
				(*outGraphName)[j] = argv[i + 1][j];
			}
			continue;
		}

		///max miller indices
		if(strstr(argv[i], "-hkl")){
			*absMaxH = strtod(argv[i + 1], NULL);
			*absMaxK = strtod(argv[i + 2], NULL);
			*absMaxL = strtod(argv[i + 3], NULL);
			continue;
		}

		///qmin, qmax
		if(strstr(argv[i], "-qm")){
			*qMin = strtod(argv[i + 1], NULL);
			continue;
		}
		if(strstr(argv[i], "-qM")){
			*qMax = strtod(argv[i + 1], NULL);
			continue;
		}
		///number of samples
		if(strstr(argv[i], "-nas")){
			*nAffSamps = strtoul(argv[i + 1], NULL, BASE);
			continue;
		}
		if(strstr(argv[i], "-nqs")){
			*nQSamps = strtoul(argv[i + 1], NULL, BASE);
			continue;
		}
		///g-smearing sigma
		if(strstr(argv[i], "-s")){
			*sigma = strtod(argv[i + 1], NULL);
			continue;
		}

		///wavelength
		if(strstr(argv[i], "-w")){
			*wavelength = strtod(argv[i + 1], NULL);
			continue;
		}
	}


	fend: return;
}


//Input------------------------------------------------------------------------
//Reads A, lattice from file
//HAS A FIXED FORMAT:
//a0x a0y a0z
//a1x a1y a1z
//a2x a2y a2z
//number of species, number of sites 
//Z a b c f (Z = atomic number, a, b, c = dir crds, f = 0.0 <= frac occ <= 1.0
//etc.
void ReadLatticeFile(const char* fileName,
					 double A[DIM][DIM],
					 uint*restrict nSpecies, 
					 uint*restrict*restrict species,
					 uint*restrict nSites, 
					 site*restrict*restrict*restrict sites){

	FILE* infile;
	infile = fopen(fileName, "r");
	if(!infile){
		printf("\nERR: unable to open %s\n", fileName);
		exit(1);
	}


	char line[LINESIZE];
	char* next;

	//Unit / supercell definition
	for(uint i = 0; i < DIM; ++i){
		fgets(line, LINESIZE, infile);
		A[i][0] = strtod(line, &next);
		for(uint j = 1; j < DIM; ++j){
			A[i][j] = strtod(next, &next);
		}
	}

	//Number of species, sites
	fgets(line, LINESIZE, infile);
	*nSpecies = strtoul(line, &next, BASE);
	*species = malloc(*nSpecies*sizeof(uint));
	*nSites = strtoul(next, NULL, BASE);
	*sites = malloc(*nSites*sizeof(site*));

	//Read sites, species
	uint specCounter = 0;
	uint lastSpecies = -1;
	for(uint i = 0; i < *nSites; ++i){
		fgets(line, LINESIZE, infile);

		site* thisSite = malloc(sizeof(site));
		thisSite->species = strtoul(line, &next, BASE);
		for(uint j = 0; j < DIM; ++j){
			thisSite->crdsD[j] = strtod(next, &next);
		}
		thisSite->occ = strtod(next, NULL);
		MoveToCell(&thisSite);

		(*sites)[i] = thisSite;

		if(thisSite->species != lastSpecies){
			(*species)[specCounter] = thisSite->species;
			lastSpecies = thisSite->species;
			specCounter++;
		}
	}


	fclose(infile);
	return;
}


//Writes diffraction profile in the following format:
//nEntries
//I0 Q0 h0 k0 l0
//I1 Q1 h1 k1 l1
// ...
//I(N-1) Q(N-1) h(N-1) k(N-1) l(N-1)
//sorted from low to high by Q values
void WriteInts(const char* outfileName,
			   const uint len, const rlpt* pts){

	FILE* outfile;
	outfile = fopen(outfileName, "w");
	if(!outfile){ ///this should really never happen, but just in case ...
		printf("\nERR: unable to open %s\n", outfileName);
		exit(1);
	}

	fprintf(outfile, "%u\n", len);
	for(uint i = 0; i < len; ++i){
		fprintf(outfile, "%f %f %i %i %i\n", pts[i].F.real, pts[i].Q, 
				pts[i].hkl[0], pts[i].hkl[1], pts[i].hkl[2]);
	}
		

	fclose(outfile);
	return;
}

//Writes graph to file in the following format:
//N Q0 Q1 Q2 ... Q(N-1) I0 I1 I2 ... I(N-1)
//with N the number of sample points, Q in inverse angstroms and
//I the UNSCALED intensity
void WriteGraph(const char* outfileName,
				const uint len, const double* x, const double* y){

	FILE* outfile;
	outfile = fopen(outfileName, "w");
	if(!outfile){ ///this should really never happen, but just in case ...
		printf("\nERR: unable to open %s\n", outfileName);
		exit(1);
	}

	fprintf(outfile, "%u", len);
	for(uint i = 0; i < len; ++i) fprintf(outfile, " %f", x[i]);
	for(uint i = 0; i < len; ++i) fprintf(outfile, " %f", y[i]);

	fclose(outfile);
	return;
}