//General
typedef unsigned int uint;
#ifdef _WIN32
#define restrict __restrict ///no, typdefs don't work here...
#endif

typedef struct site site;
typedef struct cpxd cpxd;
typedef struct rlpt rlpt;

#define NOP ;
#define min(X,Y) (((X) < (Y)) ? (X) : (Y))
#define max(X,Y) (((X) > (Y)) ? (X) : (Y))

//Math
#define PI 3.14159265358979323846
#define DIM 3u

//Work arrays
#define WARR_SIZE 16u
extern int WORK_ARR_I[WARR_SIZE];
extern double WORK_ARR_D[WARR_SIZE];

//Default runparams, IO name(s)
#define INFILE_LATTICE_DEFAULT "lattice.xrd\0"
#define OUTFILE_INT_DEFAULT "intensities.xrd\0"
#define OUTFILE_GRAPH_DEFAULT "graph.xrd\0"
#define BASE 10u
#define LINESIZE 128u
#define INFILE_LINE_MAX 1024u ///lines to read in infile before aborting
///Reciprocal space sphere max params
#define MAX_MAG_H 10
#define MAX_MAG_K 10
#define MAX_MAG_L 10
///Interpolation, graphing constants
#define QMIN_DEFAULT 0.0 ///for both AFF interpolation and graph writing
#define QMAX_DEFAULT 8.0 ///units: inverse angstroms

#define N_AFF_SAMPLES_DEFAULT 10u 
#define N_Q_SAMPLES_DEFAULT 1000u
#define SIGMA_DEFAULT 1.5
#define CUTOFF_DEFAULT 4u;

//Max number in periodic table
//If you add any special aff params, you need to increment ZMAX
//If adding more, I suggest putting them at the bottom and then setting a fake
//z number in the lattice file so that you don't screw up the rest of the 
//lookup table
#define ZMAX 99u ///up to Cf (index 0 = vacancy)
#define NCOEFFS 9u ///ai, bi (1 <= i <= 4), c

//Apply an LP correction to the intensities:
//If true, mark intensities w/ Q values corresponding to 0 > 2theta > 180
//as negative in 
#define APPLY_LP 0
#define WVLNGTH_DEFAULT 1.5406 ///angstroms, eg Cu K-alpha = 1.5406


