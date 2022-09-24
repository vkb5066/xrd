#include "Constants.h"



//Functions--------------------------------------------------------------------
///AFF
void BuildZMap(const uint nSpecies, const uint* specArr, uint map[ZMAX]);
void BuildAffMatrix(const uint nSpec, const uint nAffSamps, 
					double*restrict*restrict*restrict matrix, 
					const uint*restrict species,
					const double qMin, const double dq);
double SampleAffArr(const double sampleQ, 
				    const double q0, const double aff0, const double dAff);
void UpdateQParams(const double Q, const double qMin, const double dq,
				   const double*restrict affArr, 
				   double*restrict qLo, double*restrict qHi,
				   double*restrict affLo, double*restrict dAff);

///IO
void CLAHandler(int argc, char* argv[],
				char** infileName, 
				char** outLisName, char** outGraphName,
				int* absMaxH, int* absMaxK, int* absMaxL,
				double* qMin, double* qMax, 
				uint* nAffSamps, uint* nQSamps, double* sigma, uint* cutoff,
				double* wavelength);
void ReadLatticeFile(const char* fileName,
					 double A[DIM][DIM],
					 uint*restrict nSpecies, 
					 uint*restrict*restrict species,
					 uint*restrict nSites, 
					 site*restrict*restrict*restrict sites);
void WriteInts(const char* outfileName,
			   const uint len, const rlpt* pts);
void WriteGraph(const char* outfileName,
				const uint len, const double* x, const double* y);


///Lattice
void MoveToCell(site** s);
void RecipLattice(double B[DIM][DIM], const double A[DIM][DIM]);
void SetG(double res[DIM], const double B[DIM][DIM],
		  const int h, const int k, const int l);

///Math
void Add_vv(double res[DIM], const double a[DIM], const double b[DIM]);
void Mul_iv(double res[DIM], const int a, const double b[DIM]);
void Mul_dv(double res[DIM], const double a, const double b[DIM]);
void Dot_vv(double* res, const double a[DIM], const double b[DIM]);
void Crs_vv(double res[DIM], const double a[DIM], const double b[DIM]);
void ExpAcc_c(cpxd* res, const double a, const double x);
void Mag2_c(double* res, const cpxd a);
void GaussianProfile(uint* len, double **arr, 
					 const double sigma, const uint cutoff);
double CConvS(const uint len, const double* x, const double* p);


///XRD
int cmpr(const void* a, const void* b);
void MakeRecipPointArr(uint*restrict len, rlpt*restrict*restrict arr, 
					   const double B[DIM][DIM],
					   const int hmax, const int kmax, const int lmax,
					   const double wavelength);
double CalcPhi(const int hkl[DIM], const double abc[DIM]);
void MakeDiffGraph(const uint nSamPts, double** x, double** y, 
				   const uint nRecipPts, const rlpt* f, 
				   const double xMin, const double xMax,
				   const double sigma, const uint cutoff);

//Structures-------------------------------------------------------------------
///Complex number
struct cpxd{
	double real;
	double imag;
};

///Site structure
struct site{
	uint species; ///id # for this element (Z = 0 -> 98 for periodic table)
	double crdsD[DIM]; ///crystal coordinates a, b, c

	double occ;
};

///Variable containing h, k, l, and their respective |G| = Q
struct rlpt{
	int hkl[DIM]; ///{h, k, l}

	double Q; /// = | h*b1 + k*b2 + l*b3 |
	double lp; /// = (1 + cos^2(2theta)) / (sin^2(theta) cos(theta))

	cpxd F;
};

