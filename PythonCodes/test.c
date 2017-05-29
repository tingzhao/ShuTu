//test program. 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include "libNeuTu.h"

//a test program for understanding how MYFFT's Real_FFT_nf work, especially how fft values are packed. 
//relative to our indexing in getindex3DZfirst, we need to switch x, y, when calling Real_FFT_nf.
//Also, the elements at i % (nxx/2) == 0 are combinations of two elements with no imaginary parts. 
void testMYFFT() {
	int i,j,k;
	float *im1;
	int nxx, nyy, nzz;
	nxx = 4;
	nyy = 8; 
	nzz = 2;
	im1 = (float *) malloc(nxx*nyy*nzz*sizeof(float));
	srand(1000);
	for (k=0; k<nzz; k++) {
		for (i=0; i<nxx; i++) {
			for (j=0; j<nyy; j++) {
				im1[getIndex3DZfirst(i,j,k,nxx,nyy,nzz)] = rand()*1.0/RAND_MAX;
			}
		}
	}
	printf("Fourier:\n");
	float PI=3.1415926;
	int p,q,s;
	for (s=0; s<nzz; s++) {
		printf("s=%d\n",s);
		for (p=0; p<nxx; p++) {
			for (q=0; q<nyy; q++) {
				float rr, im;
				rr = 0; im = 0;
				for (k=0; k<nzz; k++) {
					for (i=0; i<nxx; i++) {
						for (j=0; j<nyy; j++) {
							float ang;
							ang = 2*PI*(i*p*1.0/nxx+j*q*1.0/nyy+k*s*1.0/nzz);
							rr += im1[getIndex3DZfirst(i,j,k,nxx,nyy,nzz)]*cos(ang);
							im += im1[getIndex3DZfirst(i,j,k,nxx,nyy,nzz)]*sin(ang);
						}
					}
				}
				printf("(%5.3f,%5.3f) ",rr,im);
			}
			printf("\n");
		}
	}
	
	int dims[3];
	dims[0] = nyy; 
	dims[1] = nxx;	//note the switch of nxx, nyy!
	dims[2] = nzz;	
	Real_FFT_nf(3,dims,im1);
	printf("\nMyer FFT:\n");
	for (k=0; k<nzz; k++) {
		printf("k=%d:\n",k);
		for (i=0; i<nxx; i++) {
			for (j=0; j<nyy; j++) {
				printf("%5.3f ",im1[getIndex3DZfirst(i,j,k,nxx,nyy,nzz)]);
			}
			printf("\n");
		}
	}

	free(im1);
}

int main()
{
	int nx = 10, ny = 20, nz = 3;
	float *imR, *imG, *imB;
	int ntot, i, j, k;
	long ii;
	
	ntot = nx * ny * nz;
	imR = (float *) malloc(ntot * sizeof(float));
	imG = (float *) malloc(ntot * sizeof(float));
	imB = (float *) malloc(ntot * sizeof(float));
	
	for (k=0; k<nz; k++) {
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				imR[ii] = rand()*1.0/RAND_MAX * (k+1)*10;
				imG[ii] = rand()*1.0/RAND_MAX ;
				imB[ii] = rand()*1.0/RAND_MAX * (k+1)*10;
			}
		}
	}
	createRGBTiffStackFromFloatArray(nx,ny,nz,imR,imG,imB,"temp.tif");
	
	free(imR);
	free(imG);
	free(imB);
	
	return 0;
}

