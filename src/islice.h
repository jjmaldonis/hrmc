#define BW (2.0F/3.0F)  /* bandwidth limit */
#define ABERR 1.0e-4    /* max error for a,b */

#define NSMAX 1000      /* max number of slices */
#define NPARAM  64      /* number of parameters */

#define NCMAX 256       /* max characters in file names */

#define NZMAX   103     /* max atomic number Z */

#define long32 int      /*  for 32/64 bit machines, copied from tiffsubs.h */


/* define functions at end of this file (i.e. so main can be 1st) */

void trlayer( const float x[], const float y[], const float occ[],
        const int Znum[], const int natom, 
        const float ax, const float by, const float kev,
        float **transr, float **transi, const long nx, const long ny,
        double *phirms, long *nbeams, const float k2max );

double vzatomLUT( int Z, double r );

//float *islice(int npc, int natom, int *Znum, float *x, float *y, float *z, 
float *islice(int natom, int *Znum, float *x, float *y, float *z, 
	      float *occ, float *wobble,  float ax, float by, float cz, 
	   int nx, int ny, float v0, double deltaz, float apert, double px, double py);

double chi( double **aber, double wl, double kx, double ky, double dx, double dy);

int probe(float **pixr, float **pixi, int nx, int ny, float ax, 
	  float by, double keV, float apert, double dx, double dy);

void writeModel(char *n, int natom, int *Znum, float *x, float *y, float *z, 
		float *occ, float *wobble,  float ax, float by, float cz);
