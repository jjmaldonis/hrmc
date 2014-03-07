/* islice.c

Functions to implement multislice calculation of the diffracted
intensity inside the femsim Fortran framework.

Based on autoslice.c by Earl Kirkland, version of 07/03/2008

first version 07-15-10, pmv
added free() for all memory allocated in islice() or probe() 08-10-10 pmv
*/

#include <stdio.h>      /* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /* FFT routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "islice.h"

/* global data for trlayer() */
float *kx2, *ky2;

/* parameters:
natom - number of atoms
Znun[] - atom atomic numbers
x[] - x atom positions
y[] - y atom positions
z[] - z atom positions
occ[] - atom occupancies
wobble[] - atom thermal vibration amplitudes
ax - supercell size in x
by - supercell size in y
cz - supercell size in z
nx - # of wave function pixels in x
ny - # of wave function pixels in y
v0 - incident beam energy
deltaz - slice thickness
res - Rayleigh criterion probe diameter in Angstroms
px - probe x position
py - probe y position
waver - real part of the propagated wave function.  must be an unallocated pointer
wavei - imaginary part of the progagated wave function.  must be an unallocated pointer
*/
//int WriteDMFixedFormat(FILE *fp, float **in_data, long32 size_x, long32 size_y);

//12/26/2010
//in probe function, when cycle through all pixels
//add  chi0 = chi( aber, wavlen, kx, ky, dx, dy ) to the 
// if ( ( ismoth != 0) && 
//                        ( fabs(k2-k2max) <= pixel) )
//end of this if condition by FY on 12/26/2010

//use i_aav=(wr*wr + wi*wi) instead of sqrt(wr*wr +wi*wi) by FY on 12/29/2010

//float *islice(int npc, int natom, int *Znum, float *x, float *y, float *z, float *occ, float *wobble,
//	  float ax, float by, float cz, int nx, int ny, float v0, double deltaz, float res, double px, double py)  
float *islice(int natom, int *Znum, float *x, float *y, float *z, float *occ, float *wobble,
	  float ax, float by, float cz, int nx, int ny, float v0, double deltaz, float res, double px, double py)  
{
    int lstart=0, lwobble=0;
    int ix, iy, iz, i, islice, kr, ixmid, iymid;
    int  istart, na;
    long nbeams;
    long  ltime;
    unsigned long iseed;
    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax, apert;
    float *kx, *ky, *xpos, *ypos, *i_aav, *npix;
    float dk, k, k2, k2max, scale, wavlen, rx, ry, pi, rmin, rmax, aimin, aimax, temperature;
    float tr, ti, wr, wi;
    float **transr, **transi, **waver, **wavei, **wavem, *propxr, *propxi, *propyr, *propyi, **tempr, **tempi, **pix;
    double sum, phi, t, zslice, phirms;
    long32 nxl, nyl;
    char fname[255];
    //This is for test image output
    char *image_gfx_name;
    FILE *image_file;
    int n_image_gfx_name;
    //float **image_r, **image_i; //for output image
    float **DP_image; //for output image
    int DP_size; //DP image size
    int image_x, image_y;
    int DP_x, DP_y;
    int verbose = 0; // set verbose to true or false
    /*  echo version date and get input file name */

    //for output diffraction image
    DP_size = (int) ax; //up to 1 A^-1 in DP
    DP_image = (float**) malloc2D( 2*DP_size, 2*DP_size, sizeof(float), "DP_image" );

    //printf("perform multislice with automatic slicing for pixel %d\n\n", npc);
    //sprintf(fname, "model_%d.xyz", npc);

    pi = (float) (4.0 * atan( 1.0 ));
    ixmid = nx/2;
    iymid = ny/2;
    dk = 1.0/ax;

    /*  get simulation options */
    if( deltaz < 1.0 ) {
        printf("WARNING: this slice thickness is probably too thin for autoslice to work properly.\n");
    }


    lwobble = 0;
    if( lwobble == 1 ) {
        if(verbose) printf( "Type the temperature in degrees K:\n");
        scanf( "%g", &temperature );
        /* get random number seed from time if available otherwise ask for a seed */
        if( lwobble == 1 ) {
            ltime = (long) time( NULL );
            iseed = (unsigned) ltime;
            if( ltime == -1 ) {
                printf("Type initial seed for random number generator:\n");
                scanf("%ld", &iseed);
            }
        else
            if(verbose) printf( "Random number seed initialized to %ld\n", iseed );
        }
    } else temperature = 0.0F;


    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( v0 );
    if(verbose) printf("electron wavelength = %g Angstroms\n", wavlen);

    /*  read in specimen coordinates and scattering factors */
    if(verbose) printf("%d atomic coordinates read in\n", natom );
    if(verbose) printf("Size in pixels Nx, Ny= %d x %d = %d beams\n", nx,ny, nx*ny);
    if(verbose) printf("Lattice constant a,b = %12.4f, %12.4f\n", ax,by);

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    if(verbose) printf("Total specimen range is\n %g to %g in x\n %g to %g in y\n %g to %g in z\n", xmin, xmax, ymin, ymax, zmin, zmax );
    if(verbose){ if( lwobble == 1 ) printf("Range of thermal rms displacements (300K) = %g to %g\n", wmin, wmax ); }

    /*  calculate spatial frequencies and positions for future use */
    rx = 1.0F/ax;
    ry = 1.0F/by;
    nxl = nx;
    nyl = ny;

    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
    freqn( ky, ky2, ypos, ny, by );

    /*  allocate some more arrays and initialize wavefunction */
    transr = (float**) malloc2D( nx, ny, sizeof(float), "transr" );
    transi = (float**) malloc2D( nx, ny, sizeof(float), "transi" );

    /* calculate incident probe wave function */
    waver = (float**) malloc2D( 2*nx, ny, sizeof(float), "waver" );
    wavei = waver + nx;
    apert = 1000.0*0.61*wavlen / res;
    if(verbose) printf("wavelength = %g kV.  aperture size = %g mrad.\n", wavlen, apert);
    probe(waver, wavei, nx, ny, ax, by, v0, apert, px, py);

    sum = 0.0;
    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++)
    sum += waver[ix][iy]*waver[ix][iy] + wavei[ix][iy]*wavei[ix][iy];

    if(verbose) printf("in islice, after probe, waver / wavei sum is: %g\n", sum);

    //debug
    //write probe image in gfx format
    /*image_gfx_name = (char*) malloc1D(30, sizeof(char), "probe_image");
    n_image_gfx_name = sprintf(image_gfx_name, "%s_%d_%d%s", "Probe_image", index_rot, index_pix,".gfx");

    printf("Probe image name is: %s\n", image_gfx_name);
    if (!(image_file = fopen(image_gfx_name, "w+b")))
    {
    printf("Cannot open the output probe file. Exiting.\n");
    exit(1);
    }

    WriteDMFixedFormat(image_file, waver, 2*nx, ny);
    fclose(image_file);
    free(image_gfx_name);
    */
    //finish writing probe image

    /*  calculate propagator function  */
    k2max = nx/(2.0F*ax);
    k2max = BW * k2max;
    if(verbose) printf("Bandwidth limited to a real space resolution of %f Angstroms\n", 1.0F/k2max);
    if(verbose) printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n", wavlen*k2max*1000.0F);
    k2max = k2max*k2max;

    propxr = (float*) malloc1D( nx, sizeof(float), "propxr" );
    propxi = (float*) malloc1D( nx, sizeof(float), "propxi" );
    propyr = (float*) malloc1D( ny, sizeof(float), "propyr" );
    propyi = (float*) malloc1D( ny, sizeof(float), "propyi" );

    scale = pi * ((float)deltaz);

    for( ix=0; ix<nx; ix++) {
        t = scale * ( kx2[ix]*wavlen );
        propxr[ix] = (float)  cos(t);
        propxi[ix] = (float) -sin(t);
    }
    for( iy=0; iy<ny; iy++) {
        t = scale * ( ky2[iy]*wavlen );
        propyr[iy] = (float)  cos(t);
        propyi[iy] = (float) -sin(t);
    }


    /*  iterate the multislice algorithm proper

    NOTE: zero freg is in the bottom left corner and
    expands into all other corners - not in the center
    this is required for the FFT - don't waste time rearranging */

    /* ---- start coherent method below ----------------
    (remember that waver,i[][] was initialize above) */

    /*  add random thermal displacements scaled by temperature if requested 
    remember that initial wobble is at 300K for each direction */
    if( lwobble == 1 ){
        scale = (float) sqrt(temperature/300.0) ;
        for( i=0; i<natom; i++) {
            x[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
            y[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
            z[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
        }
    }

    if(verbose) printf( "Sorting atoms by depth...\n");
    sortByZ( x, y, z, occ, Znum, natom );

    if( lwobble == 1 ){
        zmin = z[0];            /* reset zmin/max after wobble */
        zmax = z[natom-1];
        if(verbose) printf("Thickness range with thermal displacements is %g to %g (in z)\n", zmin, zmax );
    }

    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /* start a little before top of unit cell */
    istart = 0;
    islice = 1;

    while( (istart < natom) && ( zslice < (zmax+deltaz) ) ) {
        /* find range of atoms for current slice */
        na = 0;
        for(i=istart; i<natom; i++) 
            if( z[i] < zslice ) na++; else break;

        /* calculate transmission function, skip if layer empty */
        if( na > 0 ) {
            trlayer( &x[istart], &y[istart], &occ[istart],
               &Znum[istart], na, ax, by, v0, transr, transi,
               nxl, nyl, &phirms, &nbeams, k2max );
            //??? printf("average atompot comparison = %g\n", phirms/(wavlen*mm0) );
            transmit( waver, wavei, transr, transi, nx, ny );
        }

        /*  bandwidth limit */
        fft2d( waver, wavei, nxl, nyl, +1);

        /* remember: prop needed here to get anti-aliasing right */
        propagate( waver,  wavei, propxr, propxi, propyr, propyi, kx2,  ky2,  k2max, nx, ny );
        fft2d( waver, wavei, nxl, nyl, -1);

        sum = 0.0;
        for( ix=0; ix<nx; ix++)
            for( iy=0; iy<ny; iy++)
                sum += waver[ix][iy]*waver[ix][iy] +
                       wavei[ix][iy]*wavei[ix][iy];
        sum = sum * scale;

        if(verbose) printf("z= %f A, %ld beams, %d coord., \naver. phase= %f, total intensity = %f\n",
        zslice, nbeams, na, phirms, sum );

        zslice += deltaz;
        istart += na;
        islice++;
    } /* end while(istart<natom..) */

    //debug
    //write auto image in gfx format
    /*image_gfx_name = (char*) malloc1D(30, sizeof(char), "auto_image");
    n_image_gfx_name = sprintf(image_gfx_name, "%s_%d_%d%s", "Auto_image", index_rot, index_pix,".gfx");

    printf("Auto image name is: %s\n", image_gfx_name);
    if (!(image_file = fopen(image_gfx_name, "w+b")))
    {
        printf("Cannot open the output probe file. Exiting.\n");
        exit(1);
    }

    WriteDMFixedFormat(image_file, waver, 2*nx, ny);
    fclose(image_file);
    free(image_gfx_name);
    */
    //finish writing auto image



    /* end of multislice algorithm. waver and wavei now contain real-space exit-face wave function */

    rmin  = waver[0][0];
    rmax  = rmin;
    aimin = wavei[0][0];
    aimax = aimin;

    i_aav = (float *)malloc1D(ixmid/2, sizeof(float), "i_aav");
    npix = (float *)malloc1D(ixmid/2, sizeof(float), "npix");
    for(ix=0; ix<ixmid/2; ix++) {
        i_aav[ix] = 0.0;
        npix[ix] = 0.0;
    }

    /* Back to diffraction.  slice algorithm ends in real space. */
    fft2d ( waver, wavei, nx, ny, +1);
    //invert2D( pixr, nx, ny); //for output image, then need to recenter it
    //invert2D( pixi, nx, ny);

    //invert2D( waver, nx, ny); //for output image, then need to recenter it
    //invert2D( wavei, nx, ny);

    /*for output image after autoslice and image function in multislice code */
    //image_r = (float**) malloc2D(2*nx, ny, sizeof(float), "image_r");
    //image_i = image_r + nx;	

    /* calculate diffracted intensity, and diffracted intensity annular average */
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            wr = waver[ix][iy];
            wi = wavei[ix][iy];
            if( wr < rmin ) rmin = wr;
            if( wr > rmax ) rmax = wr;
            if( wi < aimin ) aimin = wi;
            if( wi > aimax ) aimax = wi;
            k = sqrt( kx2[ix] + ky2[iy] );
            kr = floor( k / (2.0*dk) );
            if (kr < ixmid/2) {
                //i_aav[kr] += sqrt( wr*wr + wi*wi );
                i_aav[kr] += ( wr*wr + wi*wi );
                npix[kr] += 1;
            }

            /* maybe call invert2D instead of this code */
            //output image, orignal waver, the zero is at the corner instead at the center
            //need to adjust it
            /*if(ix<nx/2)
                {image_x = nx/2-1-ix;}
            else
                {image_x = nx -1- ix + nx/2; }

            if(iy<ny/2)
                {image_y = ny/2-1-iy;}
            else
                {image_y = ny -1- iy + ny/2; }

            image_r[image_x][image_y]= wr;
            image_i[image_x][image_y]= wi;
            */

            //DP image, part of diffraction image
            if (ix<DP_size)
                { DP_x = DP_size - ix; }
            else if ((nx-ix) < DP_size )
                {DP_x = nx - ix + DP_size; }
            if (iy<DP_size)
                { DP_y = DP_size - iy; }
            else if ((ny-iy) < DP_size )
                {DP_y = ny - iy + DP_size; }
            if (kr < ixmid/2)
                { DP_image [DP_x][DP_y] = (wr*wr +wi*wi); }
        } //ix = 0


    //debug
    //write  image image in gfx format
    /*image_gfx_name = (char*) malloc1D(30, sizeof(char), "image_image");
    n_image_gfx_name = sprintf(image_gfx_name, "%s_%d_%d%s", "Image_image", index_rot, index_pix,".gfx");

    printf("Image image name is: %s\n", image_gfx_name);
    if (!(image_file = fopen(image_gfx_name, "w+b")))
    {
        printf("Cannot open the output probe file. Exiting.\n");
        exit(1);
    }

    WriteDMFixedFormat(image_file, waver, 2*nx, ny);
    fclose(image_file);
    free(image_gfx_name);
    */

    /*
    image_gfx_name = (char*) malloc1D(30, sizeof(char), "DP_image");
    n_image_gfx_name = sprintf(image_gfx_name, "%s_%d_%d%s", "DP_Ni_image", index_rot, index_pix,".gfx"); //DP for 99235 DRP
    //DP_Ni for AlLaNi50000

    //printf("Image image name is: %s\n", image_gfx_name);
    if (!(image_file = fopen(image_gfx_name, "w+b")))
    {
        printf("Cannot open the output probe file. Exiting.\n");
        exit(1);
    }

    WriteDMFixedFormat(image_file,DP_image, 2*DP_size, 2*DP_size);		
    fclose(image_file);
    free(image_gfx_name);
    free2D((void **) DP_image, 2*DP_size);
    */
    //free2D((void **)image_r, 2*nx);
    //finish writing probe image

    for(ix=0; ix < ixmid/2; ix++) {
        i_aav[ix] /= npix[ix];
        i_aav[ix] *= scale;
        if(verbose) printf("islice i_aav[%d] = %g\n", ix, i_aav[ix]);
    }

    if(verbose) printf( "pix range %g to %g real,\n          %g to %g imag\n",  rmin,rmax,aimin,aimax );

    /* free allocated memory */
    free(kx);
    free(kx2);
    free(xpos);
    free(ky);
    free(ky2);
    free(ypos);
    free2D((void **)transr, nx);
    free2D((void **)transi, nx);
    free2D((void **)waver, 2*nx);
    free(propxr);
    free(propxi);
    free(propyr);
    free(propyi);
    free(npix);


    return i_aav;

} /* end intensity  */


/*--------------------- trlayer() -----------------------*/
/*
  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  transr  = 2D array to get real part of specimen
                transmission function
  transi  = 2D array to get imag. part of specimen
                transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

*/
void trlayer( const float x[], const float y[], const float occ[],
     const int Znum[], const int natom, 
     const float ax, const float by, const float kev,
     float **transr, float **transi, const long nx, const long ny,
     double *phirms, long *nbeams, const float k2max )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rminsq, sum, scale, scalex, scaley;
    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */
    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rminsq = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<ny; iy++)
            transr[ix][iy] = 0.0F;

    for( i=0; i<natom; i++) {
        ixo = (int) ( x[i]/scalex );
        iyo = (int) ( y[i]/scaley );
        nx1 = ixo - idx;
        nx2 = ixo + idx;
        ny1 = iyo - idy;
        ny2 = iyo + idy;

        /* add proj. atomic potential at a local region near its center
        taking advantage of small range of atomic potential */
        for( ix=nx1; ix<=nx2; ix++) {
            rx2 = x[i] - ((double)ix)*scalex;
            rx2 = rx2 * rx2;
            ixw = ix;
            while( ixw < 0 ) ixw = ixw + nx;
            ixw = ixw % nx;
            for( iy=ny1; iy<=ny2; iy++) {
                rsq = y[i] - ((double)iy)*scaley;
                rsq = rx2 + rsq*rsq;
                if( rsq <= rmax2 ) {
                    iyw = iy;
                    while( iyw < 0 ) iyw = iyw + ny;
                    iyw = iyw % ny;
                    if( rsq < rminsq ) rsq = rminsq;
                    /* r = sqrt( r );
                    vz = occ[i] * scale * vzatom( Znum[i], r ); slow */
                    vz = occ[i] * vzatomLUT( Znum[i], rsq );
                    transr[ixw][iyw] += (float) vz;
                }
            } /* end for(iy... */
        } /* end for(ix... */
    } /* end for(i=0... */

    /* convert phase to a complex transmission function */
    sum = 0;
    for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<ny; iy++) {
            vz = scale * transr[ix][iy];
            sum += vz;
            transr[ix][iy] = (float) cos( vz );
            transi[ix][iy] = (float) sin( vz );
        }

    *phirms = sum / ( ((double)nx)*((double)ny) );

    /* bandwidth limit the transmission function */
    *nbeams = 0;
    fft2d( transr, transi, nx, ny, +1);
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
        k2 = ky2[iy] + kx2[ix];
        if (k2 < k2max) *nbeams += 1;
        else transr[ix][iy] = transi[ix][iy] = 0.0F;
    }
    fft2d( transr, transi, nx, ny, -1);
    
    return;
 }  /* end trlayer() */
 


/*------------------------ chi() ---------------------*/
/*
        calculate the aberration function chi
        with 3rd and 5th order Cs plus astigmatism plus offset

        put in a subroutine so I only have to get it right once

  input values:

  Aberration parameters are:
	C1, defocus, in nm
	A1, 2-fold astigmatism, in nm
	A2, 3-fold astigmatism, in nm
	B2, axial coma, in nm
	C3, primary spherical aberration, in um
	A3, 4-fold astigmasitm, in um
	S3, star aberration, in um
	A4, 5-fold astigmatism, in um
	D4, 3-lobe aberration, in um
	B4, axial coma, in um
	C5, 5th order spherical aberration, in mm
	A5, 5th order spherical aberration, in mm
  
In the aberrations array, the first column contains the aberration coefficient, the
second column contains the rotation angle for that aberration, the third column
contains the sine of the angle, and the fourth column contains the cosine.

        wl = wavelength (in Ang.)
        kx, ky = spatial freq (in 1/Ang)
        dx,dy = offset (in Ang)

*/
double chi( double **aber, double wl, double kx, double ky, double dx, double dy)
{
    double c;
    double kxr, kyr;

    kxr = kx;
    kyr = ky;

    c = 0.0;
      
    /* Evaluate the phase shifts from the various aberrations */
    c += -1.0*(dx*kx + dy*ky);                                                              /* probe center shift */
    c += 0.5*wl*aber[0][0]*(kx*kx + ky*ky);                                                     /* C1, defocus */

    kxr = kx*aber[1][3] + ky*aber[1][2];  kyr = ky*aber[1][3] - kx*aber[1][2];
    c += 0.5*wl*aber[1][0]*(kxr*kxr - kyr*kyr);                                                 /* A1, 2-fold astigmatism */
    kxr = kx*aber[2][3] + ky*aber[2][2];  kyr = ky*aber[2][3] - kx*aber[2][2];
    c += (1.0/3.0)*wl*wl*aber[2][0]*(pow(kxr,3) - 3*kxr*kyr*kyr);                               /* A2, 3-fold astigmatism */
    kxr = kx*aber[3][3] + ky*aber[3][2];  kyr = ky*aber[3][3] - kx*aber[3][2];
    c += wl*wl*aber[3][0]*(pow(kxr, 3) + kxr*kyr*kyr);                                          /* B2, axial coma */

    c += 0.25*pow(wl,3)*aber[4][0]*(pow(kx, 4) + 2*pow(kx,2)*pow(ky,2) + pow(ky, 4));           /* C3, primary spherical aberration */

    kxr = kx*aber[5][3] + ky*aber[5][2];  kyr = ky*aber[5][3] - kx*aber[5][2];
    c += 0.25*pow(wl,3)*aber[5][0]*(pow(kxr, 4) - 6*pow(kxr,2)*pow(kyr, 2) + pow(kyr,4));	/* A3, 4-fold astigmasitm */
    kxr = kx*aber[6][3] + ky*aber[6][2];  kyr = ky*aber[6][3] - kx*aber[6][2];
    c += pow(wl,3)*aber[6][0]*(pow(kxr,4) - pow(kyr,4));                                        /* S3, star aberration */ 
    kxr = kx*aber[7][3] + ky*aber[7][2];  kyr = ky*aber[7][3] - kx*aber[7][2];
    c += 0.2*pow(wl,4)*aber[7][0]*(pow(kxr, 5) - 10*pow(kxr,3)*pow(kyr,2) + 5*kxr*pow(kyr,4));	/* A4, 5-fold astigmatism */
    kxr = kx*aber[8][3] + ky*aber[8][2];  kyr = ky*aber[8][3] - kx*aber[8][2];
    c += pow(wl,4)*aber[8][0]*(pow(kxr, 5) - 2*pow(kxr,3)*pow(kyr,2) - 3*kxr*pow(kyr, 4));	/* D4, 3-lobe aberration */
    kxr = kx*aber[9][3] + ky*aber[9][2];  kyr = ky*aber[9][3] - kx*aber[9][2];
    c += pow(wl,4)*aber[9][0]*(pow(kxr,5) + 2*pow(kxr,3)*pow(kyr,2) + kxr*pow(kyr,4));	        /* B4, axial coma */

    c += (1.0/6.0)*pow(wl,5)*aber[10][0]*(pow(kx,6) + 3*pow(kx,4)*pow(ky,2) 
					  + 3*pow(kx,2)*pow(ky,4) + pow(ky,6));                 /* C5, 5th order spherical aberration */

    kxr = kx*aber[11][3] + ky*aber[11][2];  kyr = ky*aber[11][3] - kx*aber[11][2];
    c += (1.0/6.0)*pow(wl,5)*aber[11][0]*(pow(kxr,6) - 15*pow(kxr,4)*pow(kyr,2) 
					  + 15*pow(kxr,2)*pow(kyr,4) - pow(kyr,6));             /* A5, 5th order spherical aberration */

    c *= 2.0*3.1415927;

    /*phi = atan2( ky, kx );
    
    //        c = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
    //               + dfa2*sin( 2.0*(phi-dfa2phi) ) 
    //                + 2.0*dfa3*wavlen*sqrt(k2)* sin( 3.0*(phi-dfa3phi) )/3.0 )
    //                        - ( (dx2p*kx) + (dy2p*ky) ); */
    return( c );
}

int probe(float **pixr, float **pixi, int nx, int ny, float ax, 
	float by, double keV, float apert, double dx, double dy)
{
    int ix, iy, ixmid, iymid, i, ismoth, npixels;
    double  kx, ky, ky2, k2, k2max, wavlen, rx, ry, rx2, ry2, scale, pixel,sum, chi0;
	double **aber;
    int verbose = 0; // set verbose to true or false

	aber = (double**)malloc2D(12, 4, sizeof(double), "Aberrations array");
    for(i=0; i<12; i++) {
        aber[i][0] = 0.0;
        aber[i][1] = 0.0;
        aber[i][2] = 0.0;
        aber[i][3] = 1.0;
	}

    /*  Echo version date etc.  */
    if(verbose) printf( "probe version dated 3-jul-2008 ejk\n");
    if(verbose) printf( "calculate focused probe wave function\n\n");

    if( (nx != powerof2(nx)) || (ny != powerof2(ny)) ) {
        if(verbose) printf("Nx=%d, Ny=%d must be a power of 2,\ntry again.\n", nx, ny);
        exit( 0 );
    }
	ismoth = 1;

    /*  Calculate misc constants  */
    rx  = 1.0/ax;
    rx2 = rx * rx;
    ry  = 1.0/by;
    ry2 = ry * ry;
    ixmid = nx/2;
    iymid = ny/2;
    wavlen = wavelength( keV );
    if(verbose) printf("electron wavelength = %g Angstroms\n", wavlen);
    k2max = apert*0.001/wavlen;
    k2max = k2max * k2max;

    /*Calculate MTF 
    NOTE zero freq. is in the bottom left corner and
    expands into all other corners - not in the center
    this is required for FFT

    PIXEL = diagonal width of pixel squared
    if a pixel is on the aperture boundary give it a weight
    of 1/2 otherwise 1 or 0 */
    pixel = ( rx2 + ry2 );
    npixels = 0;

    for( iy=0; iy<ny; iy++) {
        ky = (double) iy;
        if( iy > iymid ) ky = (double) (iy-ny);
        ky = ky*ry;
        ky2 = ky*ky;
        for( ix=0; ix<nx; ix++) {
            kx = (double) ix;
            if( ix > ixmid ) kx = (double) (ix-nx);
            kx = kx*rx;
            k2 = kx*kx + ky2;
            if ( ( ismoth != 0) && 
                ( fabs(k2-k2max) <= pixel) ) {
                /*  old chi= chi1*k2* (chi2*k2-df) -
                2.0F*pi*( (dx*kx/ax) + (dy*ky/by) ); */
                chi0 = chi( aber, wavlen, kx, ky, dx, dy ); //added by Feng Yi on 12/26/2010
                pixr[ix][iy]= (float) ( 0.5 * cos(chi0));
                pixi[ix][iy]= (float) (-0.5 * sin(chi0));
                if(verbose) printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy );
            } else if ( k2 <= k2max ) {
                chi0 = chi( aber, wavlen, kx, ky, dx, dy );
                pixr[ix][iy]= (float)  cos(chi0);
                pixi[ix][iy]= (float) -sin(chi0);
                npixels++;
            } else {
                pixr[ix][iy] = pixi[ix][iy] = 0.0F;
            }
       }
    }
    if(verbose) printf("There were %d pixels inside the aperture\n", npixels );
    fft2d( pixr, pixi, nx, ny, -1);

    /*  Normalize probe intensity to unity  */

    sum = 0.0;
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) 
        sum +=  pixr[ix][iy]*pixr[ix][iy] +
                pixi[ix][iy]*pixi[ix][iy];

    scale = 1.0 / sum;
    scale = scale * ((double)nx) * ((double)ny);
    scale = (double) sqrt( scale );

	if(verbose) printf("probe sum = %g.  probe scale = %g.\n", sum, scale);

    for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] *= (float) scale;
            pixi[ix][iy] *= (float) scale;
        }

	/* free allocated memory */
	free2D((void **)aber, 12);

	if(verbose) printf("\n\n");
    return 0;

}  /* end probe() */


int WriteDMFixedFormat(FILE *fp, float **in_data, long32 size_x, long32 size_y) {

/* fixed format data types enum from Gatan web site
   http://www.gatan.com/~software/ScriptingLanguage/ScriptingWebPages/importexport.html */
enum { NULL_DATA, SIGNED_INT16_DATA, REAL4_DATA, COMPLEX8_DATA, OBSELETE_DATA,
	 PACKED_DATA, UNSIGNED_INT8_DATA, SIGNED_INT32_DATA, RGB_DATA, SIGNED_INT8_DATA, 
	 UNSIGNED_INT16_DATA, UNSIGNED_INT32_DATA , REAL8_DATA, COMPLEX16_DATA, BINARY_DATA, 
	 RGBA_FLOAT32_DATA, RGB_UINT16_DATA , RGB_FLOAT64_DATA, RGBA_FLOAT64_DATA, 
	 RGBA_UINT16_DATA, RGB_UINT8_DATA , RGBA_UINT8_DATA, LAST_DATA, 
	 OS_RGBA_UINT8_DATA = RGB_DATA }; 


    unsigned long32 endian = 0x0000FFFF;	/* indicates non-byte swapped ordering */
    long32 size[3];				/* x, y, z dimensions of the image */
    long32 depth = 4;			/* byte depth of the data (redundant in this case) */
    long32 type = REAL4_DATA;		/* data type */
	float *data;

	int i,j;

	/*  write the fixed format heaers */
	fwrite(&endian, sizeof(unsigned long32), 1, fp);
	//printf("1\n");
	size[0] = size_x;
	size[1] = size_y;
	size[2] = 1;
	fwrite(size, sizeof(long32), 3, fp);
	//printf("2\n");
	fwrite(&depth, sizeof(long32), 1, fp);
	//printf("3\n");
	fwrite(&type, sizeof(long32), 1, fp);
	//printf("4\n");

	/* get the data in one big block for an fwrite */
	//debug
	//printf("size_x is: %d,\n, size_y is: %d\n", size_x, size_y);
	data = (float*) malloc1D((size_x*size_y), sizeof(float), "output data block");
	for(j=0; j<size_y; j++) {
		for(i=0; i<size_x; i++) {
		    data[size_x*j + i] = in_data[i][(size_y-1) - j];  /*y goes from size_y-1 to 0.*/
		}
	}

	/* write the image data */
	fwrite(data, sizeof(float), (size_t) size_x*size_y, fp);

	free(data);
	return 1;
}
