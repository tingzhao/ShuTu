/*
These are C-functions used in ShuTu.py and ShuTuAutoTrace.c, ShuTuAutoTraceOneStack.c
 
Copyright (C) 2014-, Dezhe Z. Jin (dzj2@psu.edu)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <float.h>
#include <stdint.h> 
#include "image.h"
#include "libShuTu.h"

//flags for dianosis
#define growFromEndPoints_DIAG 0				//set to 1 to analyze the grow process. 
#define checkValidityOfLinkedPoint_DIAG 0		//set to 1 to analyze the performance of the validty criteria. 
#define adjustZOfLinkedPoint_DIAG 0				//set to 1 if need to analyze this process. The curves are saved tp temp3.py.

//1d index of 2d index (i,j)
long int getIndex(int i, int j, int nx, int ny)
{
	return i * ny + j;
}

//1d index of 3d index (i,j,k)
long int getIndex3D(int i, int j, int k, int nx, int ny, int nz)
{
	return i * ny * nz + j * nz + k;
}

//save image arrays. 
void saveImageData(char *filename, int nx, int ny, int nz, float *im3d, float *imFlat, float *imFlatRGB) {
	FILE *fpt;
	fpt = fopen(filename,"wb");
	fwrite(im3d,sizeof(float),nx*ny*nz,fpt);
	fwrite(imFlat,sizeof(float),nx*ny,fpt);
	fwrite(imFlatRGB,sizeof(float),nx*ny*3,fpt);
	fclose(fpt);
}

//read image arrays. 
void readImageData(char *filename, int nx, int ny, int nz, float *im3d, float *imFlat, float *imFlatRGB) {
	FILE *fpt;
	fpt = fopen(filename,"rb");
	fread(im3d,sizeof(float),nx*ny*nz,fpt);
	fread(imFlat,sizeof(float),nx*ny,fpt);
	fread(imFlatRGB,sizeof(float),nx*ny*3,fpt);
	fclose(fpt);
}

//probe tiff stack properties using Ting's program. 
void getTiffAttribute(char *filename, int *tiffinfo)
{
	Tiff *t;
	int maxDepth = 1000;
	int  j, type=-1, width=0, height=0, depth=0, nchan=0;
	t = Open_Tiff(filename,"r");
	if (t) {
		depth = 0;
		while ( ! Tiff_EOF(t)){ 
			depth += 1;
			if (depth > maxDepth) {
				printf("ERROR in getTiffAttribute, incorrect tif file %s\n",filename);
				tiffinfo[0] = height;
				tiffinfo[1] = width;
				tiffinfo[2] = depth;
				tiffinfo[3] = nchan;
				tiffinfo[4] = type;
				return;
			}		
			Advance_Tiff(t);
		}
		Rewind_Tiff(t);
		Get_IFD_Shape(t,&width,&height,&nchan);
		type = Get_IFD_Channel_Type(t,0);
		Close_Tiff(t);		
	}	
	tiffinfo[0] = height;
	tiffinfo[1] = width;
	tiffinfo[2] = depth;
	tiffinfo[3] = nchan;
	tiffinfo[4] = type;
}

//create gray scale tiff stack from float arrays
//im3d is the gray scale ndarray in python generated in NeuTu.py
void createGrayTiffStackFromFloatArray(int nx, int ny, int nz, float *im3d, char *outFilename)
{
	int i, ntot; 
	Array *image;
	uint8 *dat;
	float mmax = -1e5,mmin=1e5;
	
	ntot = nx*ny*nz;

	//get the maximum and mininum for scaling.
	for (i=0;i < ntot; i++) {
		if (im3d[i] > mmax) mmax = im3d[i];
		if (im3d[i] < mmin) mmin = im3d[i];
	}		

	image = Make_Array_With_Shape(PLAIN_KIND,UINT8_TYPE,Coord3(nz,nx,ny));
	dat = AUINT8(image);
	for (i=0; i<ntot; i++) {
		dat[i] = (uint8) ((im3d[i] - mmin)/(mmax-mmin) * 254.0);	// scale the values.
	}
	Write_Image(outFilename, image, DONT_PRESS);
	Kill_Array(image);
}

//create RGB tiff stack from float arrays
void createRGBTiffStackFromFloatArray(int nx, int ny, int nz, float *im3dR, float *im3dG, float *im3dB, char *outFilename)
{
	int i, ntot; 
	Array *image;
	uint8 *dat;
	float mmax = FLT_MIN,mmin=FLT_MAX;
	ntot = nx*ny*nz;
	//get the maximum and mininum for scaling.
	for (i=0;i < ntot; i++) {
		if (im3dR[i] > mmax) mmax = im3dR[i];
		if (im3dR[i] < mmin) mmin = im3dR[i];
		if (im3dG[i] > mmax) mmax = im3dG[i];
		if (im3dG[i] < mmin) mmin = im3dG[i];
		if (im3dB[i] > mmax) mmax = im3dB[i];
		if (im3dB[i] < mmin) mmin = im3dB[i];
	}	
	image = Make_Array_With_Shape(RGB_KIND,UINT8_TYPE,Coord3(nz,nx,ny));
	dat = AUINT8(image);
	for (i=0; i<ntot; i++) {
		dat[i] = (uint8) ((im3dR[i] - mmin)/(mmax-mmin) * 255.0);;
		dat[i+ntot] = (uint8) ((im3dG[i] - mmin)/(mmax-mmin) * 255.0);
		dat[i+ntot*2] = (uint8) ((im3dB[i] - mmin)/(mmax-mmin) * 255.0);
	}	
	Write_Image(outFilename, image, DONT_PRESS);
	Kill_Array(image);
}

//create tiff stack from slices using mylib.
void createTiffStackFromSlices(int nfiles, char **sliceFilenames, char *outFilename)
{
	int i, j, k, width, height, depth, nchan, type, kind;
	Array *plane;
	Tiff *stack, *slice;	
	
	stack = Open_Tiff(outFilename,"w");
	if (stack == NULL) {
		printf("ERROR IN createTiffStackFromSlices, failed to create file %s\n",outFilename);
		return;
	}
	//go through the slices and add to the output tiff.
	for (i=0; i<nfiles; i++) {
		slice = Open_Tiff(sliceFilenames[i],"r");
		if (slice) {
			if (!Get_IFD_Shape(slice,&width,&height,&nchan)) {
				type = Get_IFD_Channel_Type(slice,0);
				for (j=0;j<nchan;j++) {
					plane = Make_Array_With_Shape(PLAIN_KIND,type,Coord2(height,width));
					Get_IFD_Channel(slice,j,plane);
					kind = Get_IFD_Channel_Kind(slice,j);
					Add_IFD_Channel(stack,plane,kind);
					Free_Array(plane);
				}
				Update_Tiff(stack,DONT_PRESS);
			} else {
				printf("WARNING IN createTiffStackFromSlices, failed to read slice file %s\n",sliceFilenames[i]);
			}
			Close_Tiff(slice);
		} else {
			printf("WARNING IN createTiffStackFromSlices, failed to read slice file %s\n",sliceFilenames[i]);
		}
	}
	Close_Tiff(stack);
}

//read tiff stack using Ting's program. 
#define _READ_TIFF_STACK(dst) \
for (i=0; i<depth;i++) {\
	for (j = 0; j < nchan; j++) { \
		Get_IFD_Channel(t,j,plane);\
		switch (Get_IFD_Channel_Kind(t, j)) {\
		case PLAIN_CHAN:\
			im3d = im3dR;\
			break;\
		case MAPPED_CHAN:\
			im3d = im3dR;\
			break;\
		case RED_CHAN:\
			im3d = im3dR;\
			break;\
		case GREEN_CHAN:\
			im3d = im3dG;\
			break;\
		case BLUE_CHAN:\
			im3d = im3dB;\
			break;\
		default:\
			im3d = im3dR;\
			break;\
		}\
		for (k=0;k<nxy;k++) {\
			im3d[i*nxy+k] = (float) dst[k];\
		}\
	}\
	Advance_Tiff(t);\
}

void readTiffStack(char *filename, int *tiffinfo, float *im3dR, float *im3dG, float *im3dB)
{
	int i, j, k, width, height, depth, nchan, type, nxy;
	Array *plane;
	Tiff *t;
	uint8 *duint8;
	uint16 *duint16;
	uint32 *duint32;
	uint64 *duint64;
	int8 *dint8;
	int16 *dint16;
	int32 *dint32;
	int64 *dint64;
	float32 *dfloat32;
	float64 *dfloat64;
	float *im3d;
	
	height = tiffinfo[0]; 
	width  = tiffinfo[1];
	depth  = tiffinfo[2];
	nchan  = tiffinfo[3];
	type   = tiffinfo[4];
	nxy    = height * width;
	
	t = Open_Tiff(filename,"r");
	if (t) {
		plane = Make_Array_With_Shape(PLAIN_KIND,type,Coord2(height,width));
		switch (type) {
		case UINT8_TYPE:
			duint8 = ((uint8 *) (plane)->data);
			_READ_TIFF_STACK(duint8);
			break;
		case UINT16_TYPE:
			duint16 = ((uint16 *) (plane)->data);
			_READ_TIFF_STACK(duint16);
			break;
		case UINT32_TYPE:
			duint32 = ((uint32 *) (plane)->data);		 
			_READ_TIFF_STACK(duint32);
			break;
		case UINT64_TYPE:
			duint64 = ((uint64 *) (plane)->data);		
			_READ_TIFF_STACK(duint64);
			break;
		case INT8_TYPE: 
			dint8 = ((int8 *) (plane)->data);
			_READ_TIFF_STACK(dint8);
			break;	
		case INT16_TYPE: 
			dint16 = ((int16 *) (plane)->data);
			_READ_TIFF_STACK(dint16);
			break;
		case INT32_TYPE:
			dint32 = ((int32 *) (plane)->data);
			_READ_TIFF_STACK(dint32);
			break;	 
		case INT64_TYPE:
			dint64 = ((int64 *) (plane)->data);
			_READ_TIFF_STACK(dint64);
			break;
		case FLOAT32_TYPE:
			dfloat32 = ((float32 *) (plane)->data);
			_READ_TIFF_STACK(dfloat32);
			break;
		case FLOAT64_TYPE:
			dfloat64 = ((float64 *) (plane)->data);
			_READ_TIFF_STACK(dfloat64);
			break;
		}	
		Kill_Array(plane);	
		Close_Tiff(t);
	}	
}

//get tiff file dimensions. 
//tiffinfo is int tiffinfo[5]. 
void probeFileDimensions(char *filename, int *nxx, int *nyy, int *nzz, int *tiffinfo) 
{
	int nx,ny,nz,nchan,type;	
	getTiffAttribute(filename, tiffinfo);
	nx = tiffinfo[0];
	ny = tiffinfo[1];
	nz = tiffinfo[2];
	nchan = tiffinfo[3];
	type = tiffinfo[4];    
	printf("Image dimension=(%d %d %d) type=%d number of channels=%d\n",nx,ny,nz,type,nchan);
	(*nxx) = nx; (*nyy) = ny; (*nzz) = nz;
}	

//read tiff stack, return im3d.
//tiffinfo is int tiffinfo[5]. 
//must first call probeFileDimensions and allocate memory for im3d, size nx*ny*nz 
int readImageAndReturnIm3d(char *filename, int nx, int ny, int nz, float *im3d, int *tiffinfo) 
{
	int ntot, i, j, k, nxy, Wx, Wy;
	float maxI, minI, sigma, *im3dG, *im3dB, *imFlat;
	int nchan;
	float mmd;
	int itype;

	//read image. 

	ntot = nx * ny * nz;	
	nxy = nx * ny;
	nchan = tiffinfo[3];

	if (nchan == 1) {
		printf("One chanel stack. ");
		readTiffStack(filename, tiffinfo, im3d, NULL, NULL);
	} else if (nchan == 3) {
		printf("Color stack. ");
		im3dG = (float *) malloc(ntot*sizeof(float));	
		im3dB = (float *) malloc(ntot*sizeof(float));			
		readTiffStack(filename, tiffinfo, im3d, im3dG, im3dB);
		for (i=0; i<ntot; i++) {
			im3d[i] = 0.2126 * im3d[i] + 0.7152 * im3dG[i] + 0.0722 * im3dB[i];
		}
		free(im3dG);
		free(im3dB);
	} else {
		printf("ERROR: The number of channels in the image is currently assumed to be 1 or 3. \n");
		return 1;
	}	
	
	//smooth each plane in the tiff stack
	printf("Smoothing planes. ");
	sigma = 1.0;				
	Wx = (int) fmax(5,5*sigma); 
	Wy = Wx;	
	imFlat = (float *) malloc(nxy * sizeof(float));
	for (k=0; k<nz; k++) {
		i = k * nxy;
		gaussianFilter2D(nx, ny, Wx, Wy, im3d+i, imFlat, sigma);
		for (i=0;i<nxy;i++) im3d[k*nxy+i] = imFlat[i];
	}
	
	/*
	double x2, y2, z, sx, sy, Ix, Iy, dn;
	double alpha, beta, gm, xc, yc;
	long ii;
	printf("...removing tilts in planes.");
	//remove systematic linear changes of intensity in each plane.
	x2 = 0; y2 = 0; z=0; Ix = 0; Iy = 0; sx = 0; sy = 0;
	xc = nx/2.0; yc = ny/2.0;
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			x2 += (i - xc) * (i - xc);
			y2 += (j - yc) * (j - yc);
			z  += (i - xc) * (j - yc);
			sx += (i - xc);
			sy += (j - yc);
		}
	}
	dn = x2 * y2 - z * z;
	for (k=0; k<nz; k++) {
		Ix = 0; Iy = 0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				Ix += im3d[ii] * (i - xc);
				Iy += im3d[ii] * (j - yc);
			}
		}
		alpha = (Ix * y2 - Iy * z)/dn;
		beta  = (Iy * x2 - Ix * z)/dn;
		gm = - (alpha * sx + beta * sy)/nxy;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				im3d[ii] -= alpha * (i - xc) + beta * (j - yc) + gm;
			}
		}		
	}
	*/
	//reschale. 
	printf("Rescale. ");
	maxI = -FLT_MAX;
	minI = FLT_MAX;
	for (i=0; i<ntot; i++) {
		if (maxI < im3d[i]) maxI = im3d[i];
		if (minI > im3d[i]) minI = im3d[i];
	}
	if (PARAMS.imageType == 0) {
		itype = 0;	//bright field.
		printf("Bright field image. ");
	} else if (PARAMS.imageType == 1) {
		itype = 1; //dark field. 
		printf("Dark field image. ");
	} else {
		printf("Automatically detecting image type...");
		//detect background type. 
		mmd = getSparseThreshold(nx,ny,im3d,0.5);	//get the median value from the first plane. 
		if (maxI - mmd < mmd - minI) {	//bright field. 
			printf("Bright field image. ");
			itype = 0;
		} else {
			printf("Dark field image. ");
			itype = 1;
		}
	}
	for (i=0; i<ntot; i++) {
		if (itype == 1) {	//dark fild image, invert to bright field
			im3d[i] = (maxI - im3d[i])/(maxI - minI + 1e-10); //invert to bright field.
		} else { 
			im3d[i] = (im3d[i] - minI)/(maxI - minI + 1e-10); //bright field.
		} 
	}
	free(imFlat);
	return 0;
}	

//make minimum intensity projection of im3d, z index first.  
void minimumIntensityProjectionZfirst(int nx,int ny,int nz,float *im3d,float *imFlat)
{
	int i,j,k;
	long ii;
	float minI;
	//maximum intensity projection. 
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			minI = FLT_MAX;
			for (k=0; k<nz; k++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				if (minI > im3d[ii]) minI = im3d[ii];
			}
			ii = getIndex(i,j,nx,ny);
			imFlat[ii] = minI;
		}
	}
}

//make minimum intensity projection of im3d, regular x, y, z indexing. 
void minimumIntensityProjection(int nx,int ny,int nz,float *im3d,float *imFlat)
{
	int i,j,k;
	long ii;
	float minI;
	//maximum intensity projection. 
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			minI = FLT_MAX;
			for (k=0; k<nz; k++) {
				ii = getIndex3D(i,j,k,nx,ny,nz);
				if (minI > im3d[ii]) minI = im3d[ii];
			}
			ii = getIndex(i,j,nx,ny);
			imFlat[ii] = minI;
		}
	}
}

//process image. Creates .Proj.tif inputs. 
//converts dark field image to bright field and moves original one to .org.tif
int processImage(char *filenameBase, int imageType) 
{
	int ntot, i, k, nxy, nx, ny, nz;
	float maxI, minI, *im3d, *im3dG, *im3dB, *imFlatRGB, *imFlat;
	int nchan;
	int tiffinfo[5];
	char filenameOut[1000],filename[1000];
	int nRGB;
	//float sigma = 1.0;
	//int Wx = (int) fmax(5,5*sigma); 
	//int Wy = Wx;	
	
	//read image. 
	sprintf(filename,"%s.tif",filenameBase);
	probeFileDimensions(filename, &nx, &ny, &nz, tiffinfo);
	ntot = nx * ny * nz;	
	nxy = nx * ny;
	nchan = tiffinfo[3];
	if (nchan == 3 && imageType == 0) {
		nRGB = 3;
		imFlatRGB = (float *) malloc(nxy * 3 * sizeof(float));
	} else {
		nRGB = 1;
		imFlatRGB = (float *) malloc(nxy * sizeof(float));
	}	
	im3d = (float *) malloc(ntot * sizeof(float));

	if (nchan == 1) {
		printf("One chanel stack. ");
		readTiffStack(filename, tiffinfo, im3d, NULL, NULL);			
		//for (k=0; k<nz; k++) {//smooth each plane
		//	gaussianFilter2D(nx,ny,Wx,Wy,im3d+k*nxy,im3d+k*nxy,sigma);
		//}
	} else if (nchan == 3) {
		printf("Color stack. ");
		imFlat = (float *) malloc(nxy * sizeof(float));
		im3dG = (float *) malloc(ntot*sizeof(float));	
		im3dB = (float *) malloc(ntot*sizeof(float));			
		readTiffStack(filename, tiffinfo, im3d, im3dG, im3dB);
		//for (k=0; k<nz; k++) {//smooth each plane
		//	gaussianFilter2D(nx,ny,Wx,Wy,im3d+k*nxy,im3d+k*nxy,sigma);
		//	gaussianFilter2D(nx,ny,Wx,Wy,im3dG+k*nxy,im3dG+k*nxy,sigma);
		//	gaussianFilter2D(nx,ny,Wx,Wy,im3dB+k*nxy,im3dB+k*nxy,sigma);
		//}	
		if (imageType == 0) {//bright field, retain color.
			minimumIntensityProjectionZfirst(nx,ny,nz,im3d,imFlat);
			for (i=0; i<nxy; i++) imFlatRGB[i] = imFlat[i];
			minimumIntensityProjectionZfirst(nx,ny,nz,im3dG,imFlat);
			for (i=0; i<nxy; i++) imFlatRGB[i+nxy] = imFlat[i];
			minimumIntensityProjectionZfirst(nx,ny,nz,im3dB,imFlat);
			for (i=0; i<nxy; i++) imFlatRGB[i+nxy*2] = imFlat[i];
		}				
		for (i=0; i<ntot; i++) {
			im3d[i] = 0.2126 * im3d[i] + 0.7152 * im3dG[i] + 0.0722 * im3dB[i];	//make sure the coefficients add to 1.
		}
		free(im3dG);
		free(im3dB);
		free(imFlat);
	} else {
		printf("ERROR: The number of channels in the image is currently assumed to be 1 or 3. \n");
		return 1;
	}
	if (imageType == 1) {
		//convert to bright field, gray scale. 	
		maxI = FLT_MIN;
		for (i=0; i<ntot; i++) if (maxI < im3d[i]) maxI = im3d[i];
		for (i=0; i<ntot; i++) {
			im3d[i] = maxI - im3d[i];
		} 
		sprintf(filenameOut,"%s.org.tif",filenameBase);
		printf("Dark field image. Convert to bright field. Move original tif to %s.\n",filenameOut);
		rename(filename,filenameOut);
		createGrayTiffStackFromFloatArray(nx,ny,nz,im3d,filename);
	} 
	if (nchan != 3 || imageType != 0) {	
		minimumIntensityProjectionZfirst(nx,ny,nz,im3d,imFlatRGB);
	}
	sprintf(filenameOut,"%s.Proj.tif",filenameBase);
	printf("Save projection to %s\n",filenameOut);
	if (nRGB == 3) {
		createRGBTiffStackFromFloatArray(nx,ny,1,imFlatRGB,imFlatRGB+nxy,imFlatRGB+2*nxy,filenameOut);
	} else {
		createGrayTiffStackFromFloatArray(nx,ny,1,imFlatRGB,filenameOut);
	}
	free(im3d);
	free(imFlatRGB);
	return 0;
}	


//stitch tiff in filename1 to tiff in filename2. dx, dy, dz are the shifts of 2 relative to 1. 
//offsetDirection, 1 up, 2 down, 3 left, 4 right. 
//set imageType = 0, bright field, 1, dark field. 
//ret is a float array of length 4, returns dx, dy, dz, maxcorr.
#define DirectionalPairwiseStitch_DIAG 0	//set to 1 if need to analyze this function. 
void DirectionalPairwiseStitchC(int imageType, char *filename1,char *filename2,float overlapFract,int offsetDirection, float *ret)
{
	int nx, ny, nz, ntot, nxy, tiffinfo[5], i, j, k;
	char *filename[2];
	float *im3d[2], *imFlat[2], *imF1;
	int nxx, nyy, nzz, nx2, ny2;
	char outName[1000];
	float *im1, *im2;
	int i11,i12,j11,j12,k11,k12,i21,i22,j21,j22,k21,k22,is,js,ks;
	int dims[3];
	Complex_f *s;
	float nm, mmax, aa, bb, bg;
	long ii;
	int dx, dy, dz;
	float maxcorr;
	float sigma = 3.0;				
	int Wx = (int) fmax(5,5*sigma); 
	int Wy = Wx;	

	PARAMS.imageType = imageType; //this affects reading the tiff. After reading, loaded images are converted to bright field if needed. 

	filename[0] = filename1;
	filename[1] = filename2;
	printf("Stitching %s and %s. ",filename1,filename2);
	for (i=0; i<2; i++) {
		//read tiff stacks. 
		probeFileDimensions(filename[i],&nx,&ny,&nz,tiffinfo);
		ntot = nx * ny * nz;	
		nxy = nx * ny;	
		im3d[i] = (float *) malloc(ntot*sizeof(float));	
		readImageAndReturnIm3d(filename[i],nx,ny,nz,im3d[i],tiffinfo); //this returns normalized im3d
		#if DirectionalPairwiseStitch_DIAG
		imFlat[i] = (float *) malloc(nxy * sizeof(float));
		minimumIntensityProjectionZfirst(nx,ny,nz,im3d[i],imFlat[i]);
		sprintf(outName,"temp.%d.tif",i);
		printf("Saving imflat to %s\n",outName);
		createGrayTiffStackFromFloatArray(nx,ny,1,imFlat[i],outName);			
		#endif 
	}		
	//get the overlapping regions of the images. 		
	//set the dimension of the array. 
	if (overlapFract <0 || overlapFract > 1) {	//for bad values, set to default.
		printf("Bad value for overlapFrac. Set to 0.2. "); 
		overlapFract = 0.2;
	}
	if (offsetDirection == 1 || offsetDirection == 2) {//up, down.
		nx2 = (int) (nx * overlapFract);
		ny2 = ny;
	} else { //left, right. 
		ny2 = (int) (ny * overlapFract);
		nx2 = nx;		
	}
	//make the dimensions 2^n power, make it bigger to pad background.  
	nxx = (int)pow(2,(int)log2f(nx2)+1);
	nyy = (int)pow(2,(int)log2f(ny2)+1); 
	nzz = (int)pow(2,(int)log2f(nz)+1);
	ntot = nxx * nyy * nzz;
	printf("Dimensions of overlapping stack nxx=%d nyy=%d nzz=%d. ",nxx,nyy,nzz);		
	
	//get the sub regions in to arrays.
	im1 = (float *) malloc(ntot * sizeof(float));
	im2 = (float *) malloc(ntot * sizeof(float));
	//initialize with background. 
	if (PARAMS.imageType == 0) {
		bg = getSparseThreshold(nx,ny,im3d[0],0.8);
		//bg = 1.0; //bright field background
	} else {
		bg = getSparseThreshold(nx,ny,im3d[0],0.2);
		//bg = 0.0; //dark field background
	}
	for (ii=0; ii<ntot; ii++) {
		im1[ii] = bg;
		im2[ii] = bg;
	}
			
	//set ranges of indices for copying the image. 
	if (offsetDirection == 1) {//up
		i11 = 0; i12 = nx2;
		i21 = nx - nx2; i22 = nx;
		j11 = 0; j12 = j11 + ny;
		j21 = j11; j22 = j12;
	} else if (offsetDirection == 2) {//down
		i21 = 0; i22 = nx2;
		i11 = nx - nx2; i12 = nx;
		j11 = 0; j12 = j11 + ny;
		j21 = j11; j22 = j12;
	} else if (offsetDirection == 3) {//left
		j11 = 0; j12 = ny2;
		j21 = ny - ny2; j22 = ny;
		i11 = 0; i12 = nx;
		i21 = i11; i22 = i12;
	} else if (offsetDirection == 4) {//right
		j21 = 0; j22 = ny2;
		j11 = ny - ny2; j12 = ny;
		i11 = 0; i12 = nx;
		i21 = i11; i22 = i12;
	}
	//copy data to the array
	is = (nxx - nx2)/2;
	js = (nyy - ny2)/2;
	ks = (nzz - nz)/2;
	for (k=0; k<nz; k++) {
		for (i=i11; i<i12; i++) {
			for (j=j11; j<j12; j++) {
				ii = getIndex3DZfirst(is+i-i11,js+j-j11,ks+k,nxx,nyy,nzz);
				im1[ii] = im3d[0][getIndex3DZfirst(i,j,k,nx,ny,nz)];
			}
		}
	}
	for (k=0; k<nz; k++) {
		for (i=i21; i<i22; i++) {
			for (j=j21; j<j22; j++) {
				ii = getIndex3DZfirst(is+i-i21,js+j-j21,ks+k,nxx,nyy,nzz);
				im2[ii] = im3d[1][getIndex3DZfirst(i,j,k,nx,ny,nz)];
			}
		}
	}
	//smooth planes to remove boundary effects. 	
	for (k=0; k<nzz; k++) {
		gaussianFilter2D(nxx,nyy,Wx,Wy,im1+k*nxx*nyy,im1+k*nxx*nyy,sigma);
		gaussianFilter2D(nxx,nyy,Wx,Wy,im2+k*nxx*nyy,im2+k*nxx*nyy,sigma);
	}
	imF1 = (float*) malloc(nxx*nyy*sizeof(float));
	
	#if DirectionalPairwiseStitch_DIAG
	float  *imF2;
	minimumIntensityProjectionZfirst(nxx,nyy,nzz,im1,imF1);
	sprintf(outName,"temp.1sub.tif");
	printf("Saving imflat of sub slab to %s\n",outName);
	createGrayTiffStackFromFloatArray(nxx,nyy,1,imF1,outName);			
	imF2 = (float*) malloc(nxx*nyy*sizeof(float));
	minimumIntensityProjectionZfirst(nxx,nyy,nzz,im2,imF2);
	sprintf(outName,"temp.2sub.tif");
	printf("Saving imflat of sub slab to %s\n",outName);
	createGrayTiffStackFromFloatArray(nxx,nyy,1,imF2,outName);			
	#endif
	
	//phase correlation. 
	printf("Computing phase correlation. ");
	dims[0] = nyy;	//note the swith of x and y! 
	dims[1] = nxx;
	dims[2] = nzz;	
	
	Real_FFT_nf(3,dims,im1);
	Real_FFT_nf(3,dims,im2);
	//Real_Correlation_nf(3,dims,(Complex_f *)im1,(Complex_f *)im2);
	//compuet noramlized correlation. Note how the fft values are packed in im1.
	long ii1, ii2;
	for (k=0; k<nzz; k++) {
		for (i=0; i<nxx; i++) {
			for (j=0; j<nyy; j+=2) {
				ii1 = getIndex3DZfirst(i,j,k,nxx,nyy,nzz);
				ii2 = getIndex3DZfirst(i,j+1,k,nxx,nyy,nzz);
				if (j==0 && (i % (nxx/2) == 0)) {	//special pair that packs together two FF coefficients that are always real. 
					im1[ii1] = 1;
					im1[ii2] = 1;				
				} else {	//im1[ii1], im1[ii2] are the real and imaginary parts. 
					aa = im1[ii1] * im2[ii1] + im1[ii2] * im2[ii2];	//multiply C1 C2^*
					bb = im1[ii2] * im2[ii1] - im1[ii1] * im2[ii2];
					nm = sqrt(aa * aa + bb * bb)+1e-5;
					im1[ii1] = aa/nm;
					im1[ii2] = bb/nm;
				}	
			}
		}
	}
	Real_FFT_Inverse_nf(3,dims,(Complex_f *)im1);		
	
	//locate the maximum.
	//smooth the planes.
	printf("Smoothing planes. ");
	nxy = nxx*nyy;
	for (k=0; k<nzz; k++) {
		i = k * nxy;
		gaussianFilter2D(nxx, nyy, Wx, Wy, im1+i, imF1, sigma);
		for (i=0;i<nxy;i++) im1[k*nxy+i] = imF1[i];
	}
	
	mmax = FLT_MIN;
	for (k=0; k<nzz; k++) {
		for (i=0; i<nxx; i++) {
			for (j=0; j<nyy; j++) {
				aa = im1[getIndex3DZfirst(i,j,k,nxx,nyy,nzz)]; 
				if ( aa > mmax) {
					mmax = aa;
					dx = i;
					dy = j;
					dz = k;
				}
			}
		}
	}
	if (dx > nxx/2) dx = dx - nxx;
	if (dy > nyy/2) dy = dy - nyy;
	if (dz > nzz/2) dz = dz - nzz;
	maxcorr = mmax; 
	printf("Shifts dx=%d dy=%d dz=%d maxcorr=%f. ",dx,dy,dz,maxcorr);

	#if DirectionalPairwiseStitch_DIAG
	//save stack. 
	printf("Saving the phase correlation stack to temp.tif.\n");
	createGrayTiffStackFromFloatArray(nxx,nyy,nzz,im1,"temp.tif");				
	//construct composition figure. 
	float *imFC, *img1, *img2;
	int nxc, nyc, dxx, dyy;
	if (offsetDirection == 1) {
		img1 = imFlat[1];
		img2 = imFlat[0];
		dxx = nx - nx2 - dx;
		dyy = -dy;
	} else if (offsetDirection == 2) {
		img1 = imFlat[0];
		img2 = imFlat[1];
		dxx = nx - nx2 + dx;
		dyy = dy;
	} else if (offsetDirection == 3) {
		img1 = imFlat[1];
		img2 = imFlat[0];
		dyy = ny - ny2 - dy;
		dxx = -dx;		
	}  else if (offsetDirection == 4) {
		img1 = imFlat[0];
		img2 = imFlat[1];
		dyy = ny - ny2 + dy;
		dxx = dx;
	}
	i11 = (int)fmin(dxx,0);
	j11 = (int)fmin(dyy,0);
	i22 = (int)fmax(dxx+nx,nx);
	j22 = (int)fmax(dyy+ny,ny);
	nxc = i22 - i11;
	nyc = j22 - j11;
	imFC = (float *) malloc(nxc*nyc*sizeof(float));
	for (ii=0; ii<nxc*nyc; ii++) imFC[ii] = bg;
	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++) {
			imFC[getIndex(i-i11,j-j11,nxc,nyc)] = img1[getIndex(i,j,nx,ny)];
		}
	}
	for (i=dxx;i<dxx+nx;i++) {
		for (j=dyy;j<dyy+ny;j++) {
			imFC[getIndex(i-i11,j-j11,nxc,nyc)] = img2[getIndex(i-dxx,j-dyy,nx,ny)];
		}
	}
	sprintf(outName,"temp.4sub.tif");
	printf("Saving imFlatC to %s\n",outName);
	createGrayTiffStackFromFloatArray(nxc,nyc,1,imFC,outName);			
	free(imFC);
	for (ii=0; ii<ntot; ii++) im1[ii] = -im1[ii];
	minimumIntensityProjectionZfirst(nxx,nyy,nzz,im1,imF1);	
	sprintf(outName,"temp.3sub.tif");
	printf("Saving imflat of correlation to %s\n",outName);
	createGrayTiffStackFromFloatArray(nxx,nyy,1,imF1,outName);			

	free(imF2);
	free(imFlat[0]);
	free(imFlat[1]);
	#endif				
	
	free(im1);
	free(im2);
	free(imF1);
	free(im3d[0]);
	free(im3d[1]);
	
	printf("Convert shifts to offsets.\n");
	if (offsetDirection == 1) {
		dx = -nx + nx2 + dx;
	} else if (offsetDirection == 2) {
		dx =  nx - nx2 + dx;
	} else if (offsetDirection == 3) {
		dy = -ny + ny2 + dy;
	}  else if (offsetDirection == 4) {
		dy =  ny - ny2 + dy;
	}
		
	ret[0] = dx;
	ret[1] = dy;
	ret[2] = dz;
	ret[3] = maxcorr;
}


//Modified heap data structure code from http://www.thelearningpoint.net/computer-science/data-structures-heaps-with-c-program-source-code
/*Declaring heap globally so that we do not need to pass it as an argument every time*/
/* Heap used here is Min Heap */

/*Initialize Heap*/
void Init()
{
       heapSize = 0;
       heap[0].value = -FLT_MAX;
       heap[0].ind = -1;
}
/*Insert an element into the heap */
//if the element index exists, relpalce the old value and update the heap.
//assuming the updated value is smaller than before.  
void InsertOrUpdate(struct valInd element)
{
	int i, now, flag=0;
	for (i=1;i<=heapSize;i++) {
		if (heap[i].ind == element.ind) {
			flag = 1;
			break;
		}
	}
	if (flag == 0) {	//new index, increase the heap and adjust. 		
        heapSize++;
        heap[heapSize] = element; /*Insert in the last place*/
        now = heapSize;
    } else {	//existing index, replace value and adjust the heap. 
		now = i;   
	}
    /*Adjust its position*/
    while(heap[now/2].value > element.value) 
    {
         heap[now] = heap[now/2];
         now /= 2;
     }
     heap[now] = element;
}

struct valInd DeleteMin()
{
	/* heap[1] is the minimum element. So we remove heap[1]. Size of the heap is decreased. 
           Now heap[1] has to be filled. We put the last element in its place and see if it fits.
           If it does not fit, take minimum element among both its children and replaces parent with it.
           Again See if the last element fits in that place.*/
    struct valInd minElement,lastElement;
    int child,now;
    minElement = heap[1];
    lastElement = heap[heapSize--];
    /* now refers to the index at which we are now */
    for(now = 1; now*2 <= heapSize ;now = child) {
		/* child is the index of the element which is minimum among both the children */ 
        /* Indexes of children are i*2 and i*2 + 1*/
        child = now*2;
        /*child!=heapSize beacuse heap[heapSize+1] does not exist, which means it has only one 
        child */
        if(child != heapSize && heap[child+1].value < heap[child].value ) {
			child++;
        }
        /* To check if the last element fits ot not it suffices to check if the last element
           is less than the minimum element among both the children*/
        if(lastElement.value > heap[child].value) {
			heap[now] = heap[child];
        }
		else /* It fits there */
        {
			break;
        }
	}
    heap[now] = lastElement;
    return minElement;
}

float pixelDistance(int i, int j, float img, int ii, int jj, float imgFrom)
{
	float dd = exp(PARAMS.alphaDistance*img) * sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj));
	return dd;	
}

float pixelDistanceZ(int i, int j, float img, int ii, int jj, float imgFrom, float zfact)
{
	//double dd = (exp(alpha*abs(img-imgFrom)) + exp(alpha*img))/2.0 * sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj));
	float dd = exp(PARAMS.alphaDistance*img) * sqrt(zfact*zfact*(i-ii)*(i-ii)+(j-jj)*(j-jj));
	//double dd = img * sqrt(zfact*zfact*(i-ii)*(i-ii)+(j-jj)*(j-jj));
	return dd;	
}

//compute the shortest distance using Dijstra's algorithm. Computes shortest distance from one point to all other points. 
void dijstraComputeDists(int is, int js, int nx, int ny, float *img, float *dists) 
{
	int *visited, i, j, k, ii, jj, kk, ntot;
	float dd;
	struct valInd element;
	
	ntot = nx * ny;
	for (i=0;i<ntot;i++) dists[i] = FLT_MAX;
	visited = (int *) malloc(ntot * sizeof(int));
	for (i=0;i<ntot;i++) visited[i] = 0;

	k = getIndex(is,js,nx,ny);
	dists[k] = 0.0;
	visited[k] = 1;
	//allocate heat
	if (heap == NULL) heap = (struct valInd  *) malloc(nx *ny * sizeof(struct valInd));
	Init(); //heap structure for finding minimum.
	while (1) {
		//find the index of 
		ii = k / ny; 
		jj = k % ny;
		for (i=ii-1; i<ii+2; i++) {
			if (i < 0 || i >= nx) continue;
			for (j=jj-1; j<jj+2; j++) {
				if (j < 0 || j >= ny) continue;
				kk = getIndex(i,j,nx,ny);
				if (kk == k || visited[kk] == 1) continue;
				dd = pixelDistance(i,j,img[kk],ii,jj,img[k]) + dists[k];
				dists[kk] = fmin(dists[kk],dd);
				element.value = dists[kk];
				element.ind = kk;
				InsertOrUpdate(element);
			}
		}
		element = DeleteMin();
		k = element.ind;
		if (heapSize == 0) break;
		visited[k] = 1;
	}
	free(visited);
	if (heap != NULL) {
		free(heap); heap = NULL;
	} 
}

//compute the shortest path from (is,js) to (ie,je) given the distances from (is,js)
int dijstraComputePath(int ie, int je, int nx, int ny, float *dists, float *bw) 
{
	int i, j, ii, jj, k, kk, im, lp, ntot;
	float mmin;
	
	ntot = nx * ny;
	for (i=0;i<ntot;i++) bw[i] = 0;
	k = getIndex(ie,je,nx,ny);
	lp = 0;
	bw[k] = 1;
	while (lp < ntot) {
		ii = k / ny; 
		jj = k % ny;
		mmin = FLT_MAX;
		for (i=ii-1; i<ii+2; i++) {
			if (i < 0 || i >= nx) continue;
			for (j=jj-1; j<jj+2; j++) {
				if (j < 0 || j >= ny) continue;
				kk = getIndex(i,j,nx,ny);
				if (kk == k) continue;
				if (dists[kk] < mmin) {
					mmin = dists[kk];
					im = kk;
				}
			}
		}	
		if (mmin == FLT_MAX) break;
		k = im;
		bw[k] = 1;
		++lp;
		if (mmin == 0.0) break;
	}
	return lp;
}
     
float pixelDistance3D(int i, int j, int k, float img, int ii, int jj, int kk, float imgFrom, float zfact)
{
	float alpha = 20.0; //exponent in the distant measure
	//double dd = (exp(alpha*abs(img-imgFrom)) + exp(alpha*img))/2.0 * sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj)+(k-kk)*(k-kk));
	float dd = exp(alpha*img) * sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj)+zfact*zfact*(k-kk)*(k-kk));
	return dd;	
}

//dijstra shortest distance in 3D image, starting from point (is, js, ks), set heapAllocated = 1 if 
//the memory for heap is already allocated. useful if this function is repeated called. 
void dijstraComputeDists3D(int is, int js, int ks, int nx, int ny, int nz, float zfact, float *img3d, float *dists3d, int heapAllocated) 
{
	int *visited, i, j, k, ii, jj, kk, m, mm, irest, ntot;
	float dd;
	struct valInd element;
	ntot = nx * ny * nz;
	for (i=0;i<ntot;i++) dists3d[i] = FLT_MAX;
	visited = (int *) malloc(ntot * sizeof(int));
	for (i=0;i<ntot;i++) visited[i] = 0;

	m = getIndex3D(is,js,ks,nx,ny, nz);
	dists3d[m] = 0.0;
	visited[m] = 1;
	
	if (heapAllocated != 1 || heap == NULL) heap = (struct valInd *) malloc(ntot*27*sizeof(struct valInd));
	Init(); //heap structure for finding minimum.
	while (1) {
		//find the index of 
		ii = m / (ny * nz); 
		irest = m % (ny * nz);
		jj = irest / nz;
		kk = irest % nz;
		for (i=ii-1; i<ii+2; i++) {
			if (i < 0 || i >= nx) continue;
			for (j=jj-1; j<jj+2; j++) {
				if (j < 0 || j >= ny) continue;
				for (k=kk-1; k<kk+2; k++) {
					if (k< 0 || k >= nz) continue;
					mm = getIndex3D(i,j,k,nx,ny,nz);
					if (mm == m || visited[mm] == 1) continue;
					dd = pixelDistance3D(i,j,k,img3d[mm],ii,jj,kk,img3d[m],zfact) + dists3d[m];
					dists3d[mm] = fmin(dists3d[mm],dd);
					element.value = dists3d[mm];
					element.ind = mm;
					InsertOrUpdate(element);
				}
			}		
		}
		element = DeleteMin();
		m = element.ind;
		if (heapSize == 0) break;
		visited[m] = 1;
	}
	free(visited);
	if (heapAllocated != 1 && heap != NULL) {
		free(heap); heap = NULL;
	}
}

//compute the shortest path from (is,js,ks) to (ie,je,ke) given the distances from (is,js,ks)
//note that dists3d contains shortest distances to all points from (is,js,ks), and is computed using dijstraComputeDists3D.
//the path is stored in x, y, z. This array must be allocated when callig this function. 
int dijstraComputePath3D(int ie, int je, int ke, int nx, int ny, int nz, float *dists3d, float *x, float *y, float *z) 
{
	int i, j, k, ii, jj, kk, m, mm, im, lp, ntot, irest;
	float mmin;
	ntot = nx * ny * nz;
	m = getIndex3D(ie,je,ke,nx,ny,nz);
	lp = 0;
	x[lp] = ie; y[lp] = je; z[lp] = ke; lp++;
	while (lp < ntot) {
		ii = m / (ny * nz); 
		irest = m % (ny * nz);
		jj = irest / nz;
		kk = irest % nz;
		mmin = FLT_MAX;
		for (i=ii-1; i<ii+2; i++) {
			if (i < 0 || i >= nx) continue;
			for (j=jj-1; j<jj+2; j++) {
				if (j < 0 || j >= ny) continue;
				for (k=kk-1; k<kk+2; k++) {
					if (k < 0 || k >= nz) continue;
					mm = getIndex3D(i,j,k,nx,ny,nz);
					if (mm == m) continue;
					if (dists3d[mm] < mmin) {
						mmin = dists3d[mm];
						im = mm;
					}
				}	
			}
		}	
		if (mmin == FLT_MAX) break;
		m = im;
		ii = m / (ny * nz); 
		irest = m % (ny * nz);
		jj = irest / nz;
		kk = irest % nz;
		x[lp] = ii; y[lp] = jj; z[lp] = kk; 
		lp++;
		if (mmin == 0.0) break;
	}
	return lp;
}

//This function builds distance matrix starting from one point on the left.
void computeDists(int ix, int nx, int ny, float zfact, float *img,float *dists) {
	int i,j,iy,i1,i2,ii,jj,istart,iend;
	float ddmin,dd;
	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++) {
			dists[ny*i+j] = FLT_MAX;
		}
	}
	dists[ix*ny] = img[ix*ny];
	istart = ix;
	iend = ix+1;
	for (iy=1;iy<ny;iy++){
		istart -= 1;
		iend += 1;
		istart = fmax(istart,0);
		iend = fmin(iend,nx);
		for (ii=istart;ii<iend;ii++) {
			ddmin =FLT_MAX;
			i1 = fmax(0,ii-1); 
			i2 = fmin(nx,ii+2);
			for (jj=i1;jj<i2;jj++) {
				dd = pixelDistanceZ(ii,iy,img[ii*ny+iy],jj,iy-1,img[jj*ny+iy-1],zfact) + dists[jj*ny+iy-1];				
				ddmin = fmin(ddmin,dd);
			}
			dists[ii*ny+iy] = ddmin;
		}
	}
}

//This program computes the shortest path from all points from y=0 line to y = ny line. The path steps forward. from left to right
void shortestPathImageLeftToRight(int nx, int ny, float *img, float* ypath, float zfact)
{
    int i,ix,iy,ii;
    int iMin;
    float *minDists,ddmin,dd,*dists;

    dists = (float *)malloc(nx * ny * sizeof(float));
    minDists = (float *)malloc(nx * sizeof(float));

    //printf("Computing z using shortest distance...\n");
    for (ix=0;ix<nx;ix++) {
        computeDists(ix,nx,ny,zfact,img,dists);
        ddmin = FLT_MAX;
        for (ii=0; ii<nx; ii++) {
            if (dists[ii*ny+ny-1] < ddmin) {
                ddmin = dists[ii*ny+ny-1];
            }
        }
        minDists[ix] = ddmin;
    }

    //find the min distance path
    ddmin = FLT_MAX;
    for (ix=0;ix<nx;ix++) {
        if (ddmin > minDists[ix]) {
            ddmin = minDists[ix];
            iMin = ix;
        }
    }
    computeDists(iMin,nx,ny,zfact,img,dists);   
    for (iy=ny-1;iy >=0; iy--){
        ddmin = FLT_MAX;
        for (ix=0;ix<nx;ix++) {
            if (ddmin > dists[ix*ny+iy]) {
                ddmin = dists[ix*ny+iy];
                ii = ix;
            }
        }
        ypath[iy] = ii;
    }
    free(dists);
    free(minDists);
}

void MakeSet(int n, Node *nodes)
{
	int i;
	for (i=0;i<n;i++) {
		(nodes+i)->parent = nodes+i;	//pointing to itself.
		(nodes+i)->rank =0;				//initial rank 0. 
	}
}

Node* Find(Node *node) 					//find the root. 
{
	if (node->parent != node) {
		node->parent = Find(node->parent);
	}
	return node->parent;
}

void Union(Node *x, Node *y) 
{
	Node *rx, *ry;
	rx = Find(x);
	ry = Find(y);
	if (rx == ry) {
		return;
	}
	if (rx->rank < ry->rank) {		
		rx->parent = ry;
	} else if (rx->rank > ry->rank) {
		ry->parent = rx;
	} else {
		ry->parent = rx;
		rx->rank += 1;
	}	
}

//kruskal's minimum spanning tree algorithm
void kruskal(int n, float distThreshold, float *D, int *E)
{
	//n, number of verticies
	//distThreshold, if the distance is beyond this, the edge is not considered in the graph (open edge)
	//D, array of size n*n, distance matrix
	//E, array of size 2*n+1, edges of the spanning tree, the last number holds the number edges.
	
	int i,j,k,l,ip,ne;
	struct valInd edge;		//value holds the cost of the edge, ind holds the index of the edge.
							//index k = n * v + u, where the edge if between (v, u). 
	int heapAllocated = 0;						
							
	Node *nodes; 	//nodes represeting the verticies. 						
	nodes = (Node *)malloc(n * sizeof(struct Node));
	MakeSet(n,nodes);			//make initia disjoint trees. 
							
	//get the list of the ordered edges using heap data structure. 
	if (heap == NULL) {
		heap = (struct valInd *) malloc(n*n*sizeof(struct valInd));
		heapAllocated = 1;
	}
	Init(); //heap structure for finding minimum.
	for (i=0;i<n;i++) {
		for (j=i+1;j<n;j++) {
			k = i * n + j;
			if (D[k] > distThreshold) {
				continue;
			}
			edge.value = D[k];
			edge.ind = k;
			InsertOrUpdate(edge);
		}
	}
	
	ne = 0;
	while (heapSize > 0) {
		edge = DeleteMin();
		k = edge.ind;
		i = k / n;
		j = k % n;
		if (Find(nodes+i) != Find(nodes+j)) {	//add edge to the tree and merge the sets. this can be optimzed in future. 
			E[2*ne] = i;
			E[2*ne+1] = j;
			ne++;
			Union(nodes+i,nodes+j);	
		}
	}
	E[2*n] = ne;
	printf("Number of edges in the minimum spanning tree = %d\n",ne);	
	free(nodes); 
	if (heapAllocated == 1 && heap != NULL) {
		free(heap);
		heap = NULL;
	}
}


//Gradient vector flow functions. 

//1D case. 
//parameters, n, dimension of the vectors
//v, vector of n, gradient vector field
//Ix, gradient of the image. 
//mu, parameter for controlling the smoothness of v
//dt, time step, dx, lattice space, maxIter, maximum iteration number. 
void gvf1d(int n, float *v, float *Ix, float mu, float dt, float dx, int maxIter)
{
	int i,iter;
	float *v2, r, md;
	r = dt * mu/(dx *dx);
	v2 = (float *) malloc(n * sizeof(float));
	//printf("Last Ix = %f\n",Ix[n-1]);
	//printf("mu=%f dt=%f dx=%f maxIter=%d\n",mu,dt,dx,maxIter);
	for (i=0; i<n; i++) v[i] = 0.0;
	for (iter=0; iter<maxIter; iter++) {
		md = 0.0;
		for (i=1;i<n-1;i++) {
			v2[i] = r * (v[i-1] + v[i+1] - 2 * v[i]) - dt * Ix[i] *Ix[i] * (v[i] - Ix[i]);
			if (fabs(v2[i]) > md) md = fabs(v2[i]);
		}	
		for (i=1;i<n-1;i++) v[i] += v2[i];
		//printf("%f ",md);	
	}
	printf("Final error = %f\n",md);
	free(v2);
}

//2D case. 
//parameters, m, n, dimensions of the images
//v, m x n array, gradient vector field, x component
//u, m x n array, gradient vector field, y component
//Ix, m x n array, gradient of the image, x compoment 
//Iy, m x n array, gradient of the image, y compoment 
//mu, parameter for controlling the smoothness of v
//dt, time step, dx, dy lattice space, maxIter, maximum iteration number. 
//Important, I should be normalize such that its maximum value is close to 1. Otherwise the convergence can suffer. 
void gvf2d(int m, int n, float *v, float *u, float *Ix, float *Iy, float mu, float dt, float dx, float dy, int maxIter)
{
	int i, j, iter, nt, i1;
	float *v2, *u2, rx, ry, md, cm, gm;
	float tol = 1e-3;	//if relative error goes below tol the iteration stops. 
	nt = m * n;
	rx = dt * mu/(dx * dx);
	ry = dt * mu/(dy * dy);
	v2 = (float *) malloc(nt * sizeof(float));
	u2 = (float *) malloc(nt * sizeof(float));
	//printf("Last Ix = %f\n",Ix[nt-1]);
	//printf("Last Iy = %f\n",Iy[nt-1]);
	//printf("mu=%f dt=%f dx=%f dy=%f  maxIter=%d\n",mu,dt,dx,dy,maxIter);
	for (i=0;i<nt;i++) {
		v[i] = 0.0;
		u[i] = 0.0;
	}
	for (iter=0; iter<maxIter; iter++) {
		md = 0.0;
		for (i=1;i<m-1;i++) {
		for (j=1;j<n-1;j++) {	
			i1 = getIndex(i,j,m,n);
			gm = dt * (Ix[i1] * Ix[i1] + Iy[i1] * Iy[i1]);	
			cm =  rx * (v[getIndex(i-1,j,m,n)] + v[getIndex(i+1,j,m,n)] - 2 * v[i1]) 
				+ ry * (v[getIndex(i,j-1,m,n)] + v[getIndex(i,j+1,m,n)] - 2 * v[i1]); 
			v2[i1] = cm - gm * (v[i1] - Ix[i1]);
			cm =  rx * (u[getIndex(i-1,j,m,n)] + u[getIndex(i+1,j,m,n)] - 2 * u[i1]) 
				+ ry * (u[getIndex(i,j-1,m,n)] + u[getIndex(i,j+1,m,n)] - 2 * u[i1]); 
			u2[i1] = cm - gm * (u[i1] - Iy[i1]);			   
			if (fabs(v2[i1]) > md) md = fabs(v2[i1]);
			if (fabs(u2[i1]) > md) md = fabs(u2[i1]);		
		}}	
		cm = 0.0;
		for (i=1;i<m-1;i++) {
		for (j=1;j<n-1;j++) {
			i1 = getIndex(i,j,m,n);
			v[i1] += v2[i1];
			u[i1] += u2[i1];
			if (fabs(v[i1]) > cm) cm = fabs(v[i1]);
			if (fabs(u[i1]) > cm) cm = fabs(u[i1]);					
		}}
		if (md/cm < tol) break;
		//printf("%f ",md);	
	}
	printf("rel error = %f\n",md/cm);

	free(v2);
	free(u2);
}

float Dirac(float x, float sigma){
	float f, b;
	f=0.5/sigma*(1+cos(3.1415926*x/sigma));
	b = (x<=sigma) & (x>=-sigma);
	f = f*b;
	return f;
}

//single linked list of integers
LinkedList* GetLastInList(LinkedList *list)
{
	LinkedList *p;
	p = list;
	while (p->next != NULL) {
		p = p->next;
	}
	return p;
}

void AppendToList(LinkedList **list, int val)
{
	LinkedList *new,*last;
	new = (LinkedList *) malloc(sizeof(LinkedList));
	new->val = val;
	new->next = NULL;
	if ((*list) == NULL) {
		(*list) = new;
	} else {
		last = GetLastInList((*list));
		last->next = new;
	}
}

void DeleteList(LinkedList *list)
{
	LinkedList *p, *next;
	p = list;
	while (p != NULL) {
		next = p->next;
		free(p);
		p = next;
	}
	list = NULL;
}

void DeleteFirstElem(LinkedList *list)
{
	LinkedList *p;
	p = list;
	if (list == NULL) return;
	list = list->next;
	free(p);
}

//double linked list of integers
DLinkedList* GetLastInDList(DLinkedList *dlist)
{
	DLinkedList *p;
	p = dlist;
	while (p->next != NULL) {
		p = p->next;
	}
	return p;
}

void AppendToDList(DLinkedList *dlist, int val)
{
	DLinkedList *new,*last;
	new = (DLinkedList *) malloc(sizeof(DLinkedList));
	new->val = val;
	new->next = NULL;
	new->prev = NULL;
	last = GetLastInDList(dlist);
	last->next = new;
	new->prev = last;
}

void DeleteDList(DLinkedList *dlist)
{
	DLinkedList *p, *next;
	p = dlist->next;
	while (p != NULL) {
		next = p->next;
		free(p);
		p = next;
	}
	dlist = NULL;
}

void DeleteFromDList(DLinkedList *dlist, DLinkedList *pdel)
{
	if (pdel->prev == NULL && pdel->next == NULL) {
		dlist = NULL;
	} else {
		if (pdel->prev == NULL) {
			pdel->next->prev = NULL;
		} else if (pdel->next == NULL) {
			pdel->prev->next = NULL;
		} else {
			pdel->prev->next = pdel->next;
			pdel->next->prev = pdel->prev;
		}
	}		
	free(pdel);
}

//exact squared Euclidean distance transformation. 
//bw contains the binary image, and edt is the distance map. 
void sedt(int nx, int ny, float *bw, float *edt)
{
	int i, j, k, u, v, np, mp, iid, iid2, flag;
	float *h;
	int *ip, *jp; 
	LinkedList **bp, *list;
	//printf("Computing the square Euclidean distance transformation...\n"); 
	h = (float *) malloc(nx * ny * sizeof(float));
	ip = (int *) malloc(nx * ny * sizeof(int));
	bp = (LinkedList **) malloc(nx * sizeof(LinkedList *));
	for (i=0; i<nx; i++) {	// initialize the list for the boundary points. 
		bp[i] = NULL;
	}	
	//get the pixel points and the boundary points. 
	np = 0;		//pixel points
	for (i=1; i<nx-1; i++) {		//note here we ignore the edge of the image for simplicity of the code. 
		for (j=1; j<ny-1; j++) {
			iid = i * ny + j;
			h[iid] = 0;			//initialize h
			if (bw[iid] > 0) {
				ip[np] = iid;
				++np;
			} else {	//see if this is a boundary point in j direction
				flag = 0;
				for (v=j-1; v <= j+1; v++) {
					if (v != j && bw[i*ny + v] > 0) {
						flag = 1;
						break;
					}
					if (flag == 1) break; 
				}
				if (flag == 1) {	// boundary point, add to the list. 
					AppendToList(&bp[i],iid);
				}	
			}	
		}
	}
	//sweep i - index
	for (k=0; k<np; k++) {
		iid = ip[k];
		i = iid/ny;
		j = iid - i*ny;
		h[iid] = 1e10;
		list = bp[i];
		while (list != NULL) {
			iid2 = list->val - (list->val/ny)*ny;
			h[iid] = fmin(h[iid],(j - iid2 )*(j - iid2));
			list = list->next;
		}
	}
	//sweep j - index
	for (k=0; k<np; k++) {
		iid = ip[k];
		i = iid/ny;
		j = iid - i*ny;
		edt[iid] = 1e10;
		for (u=1; u<nx-1; u++) {
			iid2 = u*ny + j;
			edt[iid] = fmin(edt[iid],h[iid2] + (u - i)*(u - i));
		}
	}	
	//convert to distance. 
	for (k=0; k<nx*ny; k++) edt[k] = sqrt(edt[k]);
	
	//delete memory
	free(h);
	free(ip);
	for (i=0; i<nx; i++) {	
		DeleteList(bp[i]);
	}	
	free(bp);
}

//gaussian filter in 2d, modified from Vaa3D plugin gaussianfilter. 
void gaussianFilter2D(int nx, int ny, unsigned int Wx, unsigned int Wy, float *imgIn, float *imgOut, float sigma)
{
	// for filter kernel
	float sigma_s2 = 0.5/(sigma*sigma); // 1/(2*sigma*sigma)
	float pi_sigma = 1.0/(sqrt(2*3.1415926)*sigma); // 1.0/(sqrt(2*pi)*sigma)
	float min_val = FLT_MAX, max_val = 0;
	int i,Weight,ix,iy;
	float  *WeightsX = 0,*WeightsY=0;
	float Half,x,y,k,sum;
	float  *extension_bufferX = 0,*extension_bufferY = 0;
	float  *extStop,*extIter,*stop_line,*arrLeft,*extLeft,*extRight,*arrRight,*resIter;
	float  *weightIter, *End,*arrIter;  
	unsigned int offset;
	float *img;

	//make Wx, Wy odd numbers
	Wx = Wx/2 * 2 + 1;
	Wy = Wy/2 * 2 + 1;

	img = (float *) malloc(nx*ny*sizeof(float));
	for (i=0; i<nx*ny; i++) img[i] = imgIn[i];
    //create Gaussian kernel
    WeightsX = (float *) malloc(Wx*sizeof(float));

    // Gaussian filter equation:
    // http://en.wikipedia.org/wiki/Gaussian_blur
    //   for (unsigned int Weight = 0; Weight < Half; ++Weight)
    //   {
    //        const float  x = Half* float (Weight) / float (Half);
    //         WeightsX[(int)Half - Weight] = WeightsX[(int)Half + Weight] = pi_sigma * exp(-x * x *sigma_s2); // Corresponding symmetric WeightsX
    //    }
	Half = (float)(Wx-1)/2.0;
    for (Weight = 0; Weight <= Half; ++Weight){
		x = Weight -Half;
        WeightsX[Weight]= pi_sigma * exp(-(x * x *sigma_s2)); // Corresponding symmetric WeightsX
		WeightsX[Wx-Weight-1] = WeightsX[Weight];
    }

    k = 0.;
    for (Weight = 0; Weight < Wx; ++Weight) k += WeightsX[Weight];
    for (Weight = 0; Weight < Wx; ++Weight) WeightsX[Weight] /= k;		 

    //   Allocate 1-D extension array
    extension_bufferX = (float *) malloc((nx + (Wx<<1))*sizeof(float)); //size, nx + 2*Wx.
    offset = Wx>>1;			//Half. 

    //	along x
    extStop = extension_bufferX + nx + offset;	
         
	for(iy = 0; iy < ny; iy++) {				//column
		extIter = extension_bufferX + Wx;		//copy the values of the column, starting point.
        for(ix = 0; ix < nx; ix++) {
			*(extIter++) = img[ix*ny + iy]; //copy img to extension_bufferX, starting from extension_bufferX + Wx
        }

        //   Extend image
        stop_line = extension_bufferX - 1;
        extLeft = extension_bufferX + Wx - 1;   //=extension_bufferX + Wx - 1
        arrLeft = extLeft + 2;					//=extension_bufferX + Wx + 1
        extRight = extLeft + nx + 1;			//=extension_bufferX + Wx + nx
        arrRight = extRight - 2;				//=extension_bufferX + Wx + nx - 2
        while (extLeft > stop_line){
			*(extLeft--) = *(arrLeft++);		//reflect around extension_bufferX + Wx
            *(extRight++) = *(arrRight--);		//reflect around extension_bufferX + Wx + nx - 1
		}

        //	Filtering
        extIter = extension_bufferX + Wx - offset;	//=extension_bufferX + Wx/2, original code had a bug!!
		resIter = &(img[iy]);
        while (extIter < extStop) {
			sum = 0.;
            weightIter = WeightsX;
            End = WeightsX + Wx;
            arrIter = extIter;
            while (weightIter < End)
				sum += *(weightIter++) * 1.0 * (*(arrIter++));
            extIter++;
            *(resIter) = sum;
            resIter += ny;

            //for rescale
            if(max_val<*arrIter) max_val = *arrIter;
            if(min_val>*arrIter) min_val = *arrIter;
        }
       
	 }
     //de-alloc
     free(WeightsX); WeightsX=0;
     free(extension_bufferX); extension_bufferX=0;

     //create Gaussian kernel
     WeightsY = (float *)malloc(Wy * sizeof(float));
	 Half = (float)(Wy-1)/2.0;
     for (Weight = 0; Weight <= Half; ++Weight) {
		y = Weight-Half;
        WeightsY[Weight] =  pi_sigma * exp(-(y * y *sigma_s2)); // Corresponding symmetric WeightsY
        WeightsY[Wy-Weight-1] = WeightsY[Weight];
     }

	 k = 0.;
     for (Weight = 0; Weight < Wy; ++Weight) k += WeightsY[Weight];
     for (Weight = 0; Weight < Wy; ++Weight) WeightsY[Weight] /= k;

     //	along y
     extension_bufferY = (float *) malloc((ny + (Wy<<1))*sizeof(float));
     offset = Wy>>1;
     extStop = extension_bufferY + ny + offset;

     for(ix = 0; ix < nx; ix++) {
		extIter = extension_bufferY + Wy;
        for(iy = 0; iy < ny; iy++) {
			*(extIter++) = img[ix*ny + iy];
        }

       //   Extend image
       stop_line = extension_bufferY - 1;
       extLeft = extension_bufferY + Wy - 1;
       arrLeft = extLeft + 2;
       extRight = extLeft + ny + 1;
       arrRight = extRight - 2;

       while (extLeft > stop_line) {
		 *(extLeft--) = *(arrLeft++);
         *(extRight++) = *(arrRight--);
       }

      //	Filtering
      extIter = extension_bufferY + Wy - offset;		//original code had a bug!!

      resIter = &(img[ix*ny]);

      while (extIter < extStop){
		sum = 0.;
        weightIter = WeightsY;
        End = WeightsY + Wy;
        arrIter = extIter;
        while (weightIter < End)
			sum += *(weightIter++) * 1.0 * (*(arrIter++));
        extIter++;
        *(resIter++) = sum;

        //for rescale
        if(max_val<*arrIter) max_val = *arrIter;
        if(min_val>*arrIter) min_val = *arrIter;
      }
               
	}
	for (i=0;i<nx*ny;i++) imgOut[i] = img[i];
    //de-alloc
    free(WeightsY); WeightsY=0;
    free(extension_bufferY); extension_bufferY=0;
    free(img);
}

//double linked lists
PT *pt_create(long x, long y, long idx){
  PT* newpt = (PT*)malloc(sizeof(PT));
  if(newpt == NULL) return NULL;

  newpt->x = x;
  newpt->y = y;
  newpt->idx = idx;
  newpt->prev = NULL;
  newpt->next = NULL;
  return newpt;
}

LL *ll_create(){
  LL *newll = (LL*)malloc(sizeof(LL));
  if(newll == NULL) return NULL;
  newll->head = NULL;
  newll->curr = NULL;
  newll->length = 0;
  return newll;
}

void ll_push(LL *list, PT *add){
  if(add == NULL || list == NULL) return;
  add->next = list->head;
  add->prev = NULL;
  if(add->next != NULL){
    add->next->prev = add;
  }
  list->head = add;
  list->length++;
}

void ll_pushnew(LL *list, long x, long y, long idx){
  if(list == NULL) return;
  PT* add = pt_create(x,y,idx);
  if(add == NULL) return;
  ll_push(list,add);
}

void ll_destroy(LL *list){
  if(list==NULL) return;
  while(list->head != NULL){
    ll_pop_free(list);
  }
  free(list);
}

void ll_remcurr_free(LL *list){
  PT* p = ll_remcurr(list);
  if(p != NULL) free(p);
}

PT *ll_remcurr(LL *list){
  if(list == NULL) return NULL;
  PT* out = list->curr;
  if(out == list->head){
    return ll_pop(list);
  }
  else
  {
    if(out != NULL)
    {
      if(out->next != NULL){
        out->next->prev = out->prev;
      }
      if(out->prev != NULL){
        out->prev->next = out->next;
      }
      list->curr = out->next;
      list->length--;
    }
    return out;
  }
}

void ll_pop_free(LL *list){
  PT* p = ll_pop(list);
  if(p != NULL) free(p);
}

PT *ll_pop(LL *list){
  if(list == NULL) return NULL;
  PT *out = list->head;
  if(out != NULL){
    list->head = out->next;
    if(list->curr == out) list->curr = list->head;
    if(list->head != NULL) list->head->prev = NULL;
    list->length--;
  }
  return out;
}

void ll_init(LL *list){
  if(list == NULL) return;
  list->curr = list->head;
}

void ll_step(LL *list){
  if(list == NULL) return;
  if(list->curr != NULL){
    list->curr = list->curr->next;
  }
}

//Chan-Vese levelset method. 
/*
	This file implements the sparse field chan vese argorithm for level set segmentation. 
	The implementation follows 
	"Sparse Field Methods - Technical Report", Shawn Lankton, July 6, 2009

	assume:
	img, imgae arrary
	init, initial mask, 1 foreground, 0 background
*/
//get the neiboring point
void getNp(int i, int j, long nx, long ny, int *np) {
	if (i==0) {
		np[0] = -1;
	} else {
		np[0] = getIndex(i-1,j,nx,ny);
	}
	if (i==nx-1) {
		np[1] = -1;
	} else {
		np[1] = getIndex(i+1,j,nx,ny);
	}
	if (j==0) {
		np[2] = -1;
	} else {
		np[2] = getIndex(i,j-1,nx,ny);
	}
	if (j==ny-1) {
		np[3] = -1;
	} else {
		np[3] = getIndex(i,j+1,nx,ny);
	}
}

//Chan-Vese levelset segmentation. 
//nx,ny - dimensions of img, initial mask bw. Returns with binary bw updated.
//lambda parameter for minimizing length
//dt time step
void sparseFieldChanVese(int nx, int ny, float *img, float *bw, int iter, float lambda)
{
  float *label;
  LL *Lz, *Ln1, *Ln2, *Lp1, *Lp2;
  LL *Sz, *Sn1, *Sn2, *Sp1, *Sp2;
  int i, j, k, m, n, idx, it;
  int np[4];
  short flag;
  short offx[4] ={-1,1,0,0},offy[4]={0,0,-1,1}; //corresponding offsets in getNp
  const float eps = 1.0e-10;
  float mu1,mu2,a1,a2,n1,n2;
  long ntot = nx * ny;
  float *phi, *F, kappa, Fm, dx,dxx,dx2,dy,dyy,dy2,dxy;
  int nout = (int) iter * 0.1;
			  
  //create linked lists
  Lz  = ll_create();
  Ln1 = ll_create();
  Ln2 = ll_create();
  Lp1 = ll_create();
  Lp2 = ll_create();
 
  Sz  = ll_create();
  Sn1 = ll_create();
  Sn2 = ll_create();
  Sp1 = ll_create();
  Sp2 = ll_create();

  label = (float *) malloc(ntot * sizeof(float));
  phi = (float *) malloc(ntot * sizeof(float));

  //Procedure 1 Initialization
  //Pre-condition labelmap and phi
  printf("In sparseFieldChanVese...");
  for (i=0;i<nx;i++) { 
  for (j=0;j<ny;j++) {
	idx  = getIndex(i,j,nx,ny);
    if (bw[idx] == 0) {
    	label[idx] = 3.0; 
        phi[idx] = 3.0;
	} else if (bw[idx] == 1) {
		label[idx] = -3.0;
		phi[idx] = -3.0;
	}
  }}	

  // Find the zero-level set
  for (i=0;i<nx;i++) { 
  for (j=0;j<ny;j++) {
	idx  = getIndex(i,j,nx,ny);
	if (bw[idx] == 1) {
		//check neighbors, if any bw =0 in neighbors.
		getNp(i,j,nx,ny,np);
		for (k=0; k<4;k++) {
			if (np[k] == -1) { continue; }
			if (bw[np[k]] == 0) {
				ll_pushnew(Lz,i,j,idx);
				label[idx] = 0.0;
				phi[idx] = 0.0;
				break;
			}
		}
	}
  }}	

  // Find the +1 and -1 level set
  //scan Lz to create Ln1 and Lp1
  ll_init(Lz);
  while(Lz->curr != NULL){
  	i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== -3) {
			ll_pushnew(Ln1,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=-1; 
			phi[np[k]]=-1;
		} else if (label[np[k]] == 3) {
			ll_pushnew(Lp1,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=1; 
			phi[np[k]]=1;
		}
	}	
	ll_step(Lz);
  }

  //Find the +2 and -2 level set
  ll_init(Ln1);
  while(Ln1->curr != NULL){
  	i = Ln1->curr->x; j = Ln1->curr->y; idx = Ln1->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== -3) {
			ll_pushnew(Ln2,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=-2; 
			phi[np[k]]=-2;
		}
	}
	ll_step(Ln1);
  }

  ll_init(Lp1);
  while(Lp1->curr != NULL){
  	i = Lp1->curr->x; j = Lp1->curr->y; idx = Lp1->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== 3) {
			ll_pushnew(Lp2,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=2; 
			phi[np[k]]=2;
		}
	}
	ll_step(Lp1);
  }

  //iterate, move the level set using the force. 
  //compute mu1, mu2, as the average of the intensity inside or outside of the zero level
	
  //fprintf(stderr,"Moving the level set...\n");
  for (it=0;it<iter;it++) { 	 	  
	//compute mu1, mu2
	a1 = 0; n1 = 0; a2 = 0; n2 = 0;
	for (idx=0;idx<ntot;idx++) { 
		if (phi[idx] < 0) {
			a1 += img[idx];
			n1 += 1;
		} else {
			a2 += img[idx];
			n2 += 1;
		}
	}
	mu1 = a1/(n1+eps);
	mu2 = a2/(n2+eps);
		
	// compute force along the zero level
	// allocate space for F
  	F = (float*)malloc(Lz->length*sizeof(float)); 
  	 
	Fm = eps;
  	ll_init(Lz);
	k = 0;
  	while(Lz->curr != NULL){
  		i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;

		// compute curvature
		dx = 0; dxx = 0;
		if ((i > 0) && (i < nx -1)) {
			m = getIndex(i+1,j,nx,ny);
			n = getIndex(i-1,j,nx,ny);
			dx = (phi[m] - phi[n])/2.0; 
			dxx = phi[n] + phi[m] - 2.0 * phi[idx];
		}
		dx2 = dx*dx;
		dy = 0; dyy = 0;
		if ((j > 0) && (j < ny -1)) {
			m = getIndex(i,j+1,nx,ny);
			n = getIndex(i,j-1,nx,ny);
			dy = (phi[m] - phi[n])/2.0;
			dyy = phi[n] + phi[m] - 2.0 * phi[idx];
		}
		dy2 = dy*dy;
		dxy = 0;
		if ((i > 0) && (i < nx -1) && (j > 0) && (j < ny -1)) {
			dxy = (phi[getIndex(i-1,j-1,nx,ny)] + phi[getIndex(i+1,j+1,nx,ny)] - phi[getIndex(i+1,j-1,nx,ny)] - phi[getIndex(i-1,j+1,nx,ny)])/4.0;
		} 
		kappa = ((dxx + dyy)*(dx2+dy2) - (dx2 *dxx + dy2*dyy + 2*dx*dy*dxy))/pow((dx2+dy2+eps),1.5); //curvature.
		// compute force 
		F[k] = (img[idx] - mu1)*(img[idx] - mu1) - (img[idx] - mu2)*(img[idx] - mu2) + lambda * kappa;
		if (fabs(F[k]) > Fm) {
			Fm = fabs(F[k]);
		}	
		++k;
		ll_step(Lz);
	}
	//rescale if Fm > 0.5.
	for (k=0;k<Lz->length;k++) {
		if (Fm > 0.5) {
			F[k] = 0.5 * F[k]/Fm;
		} 	
	}

	// Update the zero level set
  	ll_init(Lz);
	k = 0;
  	while(Lz->curr != NULL){
  		i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;
		phi[idx] += F[k];	
		if (phi[idx] > 0.5) {
			ll_push(Sp1, ll_remcurr(Lz));
		} else if (phi[idx] < -0.5) {
			ll_push(Sn1, ll_remcurr(Lz));
		} else {
			ll_step(Lz);
		}
		++k;
	}
	free(F);
	
	// Update -1 and +1 level sets
  	ll_init(Ln1);
  	while(Ln1->curr != NULL){
  		i = Ln1->curr->x; j = Ln1->curr->y; idx = Ln1->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = -1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 0) flag = 1;
			if (label[np[k]] >=0) {
				if (Fm < phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			ll_push(Sn2, ll_remcurr(Ln1));
		} else {
			phi[idx] = Fm - 1.0;
			if (phi[idx] >= -0.5) {
				ll_push(Sz, ll_remcurr(Ln1));
			} else if (phi[idx] < -1.5) {
				ll_push(Sn2, ll_remcurr(Ln1));
			} else {
				ll_step(Ln1);
			}
		}
	}

  	ll_init(Lp1);
  	while(Lp1->curr != NULL){
  		i = Lp1->curr->x; j = Lp1->curr->y; idx = Lp1->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = 1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 0) flag = 1;
			if (label[np[k]] <= 0) {
				if (Fm > phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			ll_push(Sp2, ll_remcurr(Lp1));
		} else {
			phi[idx] = Fm + 1.0;
			if (phi[idx] <= 0.5) {
				ll_push(Sz, ll_remcurr(Lp1));
			} else if (phi[idx] > 1.5) {
				ll_push(Sp2, ll_remcurr(Lp1));
			} else {
				ll_step(Lp1);
			}
		}
	}
	
	// Update -2 and +2 level sets
  	ll_init(Ln2);
  	while(Ln2->curr != NULL){
  		i = Ln2->curr->x; j = Ln2->curr->y; idx = Ln2->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = -1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == -1) flag = 1;
			if (label[np[k]] >= -1) {
				if (Fm < phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			phi[idx] = -3.0;
			label[idx] = -3.0;
			ll_remcurr_free(Ln2);
		} else {
			phi[idx] = Fm - 1.0;
			if (phi[idx] >= -1.5) {
				ll_push(Sn1, ll_remcurr(Ln2));
			} else if (phi[idx] < -2.5) {
				phi[idx] = -3.0;
				label[idx] = -3.0;
				ll_remcurr_free(Ln2);
			} else {
				ll_step(Ln2);
			}
		}
	}

  	ll_init(Lp2);
  	while(Lp2->curr != NULL){
  		i = Lp2->curr->x; j = Lp2->curr->y; idx = Lp2->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = 1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 1) flag = 1;
			if (label[np[k]] <= 1) {
				if (Fm > phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			phi[idx] = 3.0;
			label[idx] = 3.0;
			ll_remcurr_free(Lp2);
		} else {
			phi[idx] = Fm + 1.0;
			if (phi[idx] <= 1.5) {
				ll_push(Sp1, ll_remcurr(Lp2));
			} else if (phi[idx] > 2.5) {
				phi[idx] = 3.0;
				label[idx] = 3.0;
				ll_remcurr_free(Lp2);
			} else {
				ll_step(Lp2);
			}
		}
	}
	
    // Move points into zero level set.
    ll_init(Sz);
	while (Sz->curr != NULL){
    	idx= Sz->curr->idx;
    	ll_push(Lz,ll_remcurr(Sz));
    	label[idx] = 0;
  	}
  	
  	// Move points into -1 and +1 level sets, and ensure -2, +2 neighbors
    ll_init(Sn1);
	while (Sn1->curr != NULL){
  		i = Sn1->curr->x; j = Sn1->curr->y; idx = Sn1->curr->idx;
  		label[idx] = -1;
  		ll_push(Ln1,ll_remcurr(Sn1));
		getNp(i,j,nx,ny,np);
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (phi[np[k]] == -3.0) {
				phi[np[k]] = phi[idx] - 1.0;
				ll_pushnew(Sn2,i+offx[k],j+offy[k],np[k]);
			}
		}
  	}

    ll_init(Sp1);
	while (Sp1->curr != NULL){
  		i = Sp1->curr->x; j = Sp1->curr->y; idx = Sp1->curr->idx;
  		label[idx] = 1;
  		ll_push(Lp1,ll_remcurr(Sp1));
		getNp(i,j,nx,ny,np);
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (phi[np[k]] == 3.0) {
				phi[np[k]] = phi[idx] + 1.0;
				ll_pushnew(Sp2,i+offx[k],j+offy[k],np[k]);
			}
		}
  	}
 
 	//Move points into -2 and +2 level sets
    ll_init(Sn2);
	while (Sn2->curr != NULL){
  		idx = Sn2->curr->idx;
  		label[idx] = -2;
  		ll_push(Ln2,ll_remcurr(Sn2));
	}

    ll_init(Sp2);
	while (Sp2->curr != NULL){
  		idx = Sp2->curr->idx;
  		label[idx] = 2;
  		ll_push(Lp2,ll_remcurr(Sp2));
	}	

	//if ((it/nout)*nout == it) {
	//	fprintf(stderr, "At it=%d\n",it);
	//}
			
  }	

  for (i=0; i<ntot; i++) {
	if (phi[i] < 0) {
		bw[i] = 1.0;
	} else {
		bw[i] = 0.0;
	}
  }

  free(phi);
  free(label);

  //destroy linked lists
  ll_destroy(Lz);
  ll_destroy(Ln1);
  ll_destroy(Ln2);
  ll_destroy(Lp1);
  ll_destroy(Lp2);

  ll_destroy(Sz);
  ll_destroy(Sn1);
  ll_destroy(Sn2);
  ll_destroy(Sp1);
  ll_destroy(Sp2);
  return;
}

void NeumannBoundCond(int nrow, int ncol, float *g)
{
    // Make a function satisfy Neumann boundary condition
	int i;
	
    g[getIndex(0,0,nrow,ncol)] = g[getIndex(2,2,nrow,ncol)];
    g[getIndex(0,ncol-1,nrow,ncol)] = g[getIndex(2,ncol-3,nrow,ncol)];
    g[getIndex(nrow-1,0,nrow,ncol)] = g[getIndex(nrow-3,2,nrow,ncol)];
    g[getIndex(nrow-1,ncol-1,nrow,ncol)] = g[getIndex(nrow-3,ncol-3,nrow,ncol)];

	for (i=1;i<ncol-1;i++) {
		g[getIndex(0,i,nrow,ncol)] = g[getIndex(2,i,nrow,ncol)];
		g[getIndex(nrow-1,i,nrow,ncol)] = g[getIndex(nrow-3,i,nrow,ncol)];
	}
	for (i=1;i<nrow-1;i++) {
		g[getIndex(i,0,nrow,ncol)] = g[getIndex(i, 2,nrow,ncol)];
		g[getIndex(i,ncol-1,nrow,ncol)] = g[getIndex(i,ncol-3,nrow,ncol)];
    }
}

//nx,ny - dimensions of matrices.  
//g - gradient indicator. 
//bw - initial mask. Returns with the final mask.
//iter - number of iterations.
//mu - coefficient for curvature.  
void sparseFieldLevelSet(int nx, int ny, float *g, float *bw, int iter, float mu)
{
  float *label;
  LL *Lz, *Ln1, *Ln2, *Lp1, *Lp2;
  LL *Sz, *Sn1, *Sn2, *Sp1, *Sp2;
  int i, j, k, m, n, idx, it;
  int np[4];
  short flag;
  short offx[4] ={-1,1,0,0},offy[4]={0,0,-1,1}; //corresponding offsets in getNp
  const float eps = 1.0e-5;
  long ntot = nx * ny;
  float *F, kappa, *reg, gx, gy, Fm, dx,dxx,dx2,dy,dyy,dy2,dxy;
  int nout = (int) iter * 0.1;
  float *phi,norm,mmin,mmax;
  int flagConverged = 0;
  float epsConverge = 0.01;	//If the maximum mmovements of levelset is 0.01, judged as converged. 
  
  //rescale g to the interval 0, 1.
  mmin = FLT_MAX;
  mmax = FLT_MIN;
  for (i=0; i<ntot; i++) {
	  if (g[i] > mmax) mmax = g[i];
	  if (g[i] < mmin) mmin = g[i];
  }
  for (i=0; i<ntot; i++) {
	  g[i] = (g[i] - mmin)/(mmax - mmin + eps);
  }
  
  //create linked lists
  Lz  = ll_create();
  Ln1 = ll_create();
  Ln2 = ll_create();
  Lp1 = ll_create();
  Lp2 = ll_create();
 
  Sz  = ll_create();
  Sn1 = ll_create();
  Sn2 = ll_create();
  Sp1 = ll_create();
  Sp2 = ll_create();

  label = (float *) malloc(ntot*sizeof(float));
  phi = (float *) malloc(ntot*sizeof(float));
  
  //Procedure 1 Initialization
  //Pre-condition labelmap and phi
  printf("In sparseFieldLevelset...");
  for (i=0;i<nx;i++) { 
  for (j=0;j<ny;j++) {
	idx  = getIndex(i,j,nx,ny);
    if (bw[idx] == 1.0) {
    	label[idx] = 3.0; 
        phi[idx] = 3.0;
	} else if (bw[idx] == 0.0) {
		label[idx] = -3.0;
		phi[idx] = -3.0;
	}
  }}	

  // Find the zero-level set
  for (i=0;i<nx;i++) { 
  for (j=0;j<ny;j++) {
	idx  = getIndex(i,j,nx,ny);
	if (bw[idx] == 1) {
		//check neighbors, if any init =0 in neighbors.
		getNp(i,j,nx,ny,np);
		for (k=0; k<4;k++) {
			if (np[k] == -1) { continue; }
			if (bw[np[k]] == 0) {
				ll_pushnew(Lz,i,j,idx);
				label[idx] = 0.0;
				phi[idx] = 0.0;
				break;
			}
		}
	}
  }}	

  // Find the +1 and -1 level set
  //scan Lz to create Ln1 and Lp1
  ll_init(Lz);
  while(Lz->curr != NULL){
  	i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== -3) {
			ll_pushnew(Ln1,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=-1; 
			phi[np[k]]=-1;
		} else if (label[np[k]] == 3) {
			ll_pushnew(Lp1,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=1; 
			phi[np[k]]=1;
		}
	}	
	ll_step(Lz);
  }

  //Find the +2 and -2 level set
  ll_init(Ln1);
  while(Ln1->curr != NULL){
  	i = Ln1->curr->x; j = Ln1->curr->y; idx = Ln1->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== -3) {
			ll_pushnew(Ln2,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=-2; 
			phi[np[k]]=-2;
		}
	}
	ll_step(Ln1);
  }

  ll_init(Lp1);
  while(Lp1->curr != NULL){
  	i = Lp1->curr->x; j = Lp1->curr->y; idx = Lp1->curr->idx;
 	getNp(i,j,nx,ny,np);
	for (k=0; k<4;k++) {
		if (np[k] == -1) continue;
 		if (label[np[k]]== 3) {
			ll_pushnew(Lp2,i+offx[k],j+offy[k],np[k]);
			label[np[k]]=2; 
			phi[np[k]]=2;
		}
	}
	ll_step(Lp1);
  }

  //iterate, move the level set using the force. 
	
  //fprintf(stderr,"Moving the level set...\n");
  for (it=0;it<iter;it++) {
		
	NeumannBoundCond(nx,ny,phi);	//make sure the correct boundary condition is met. 	
		
	// compute force along the zero level
	// allocate space for F
  	F = (float*)malloc(Lz->length*sizeof(float)); 
  	 
	Fm = FLT_MIN;
  	ll_init(Lz);
	k = 0;
  	while(Lz->curr != NULL){
  		i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;

		// compute curvature and gradiant.
		dx = 0; dxx = 0; gx = 0; 
		if ((i > 0) && (i < nx -1)) {
			m = getIndex(i+1,j,nx,ny);
			n = getIndex(i-1,j,nx,ny);
			dx = (phi[m] - phi[n])/2.0; 
			dxx = phi[n] + phi[m] - 2.0 * phi[idx];
			gx = (g[m] - g[n])/2.0;
		}
		dx2 = dx*dx;
		dy = 0; dyy = 0; gy = 0;
		if ((j > 0) && (j < ny -1)) {
			m = getIndex(i,j+1,nx,ny);
			n = getIndex(i,j-1,nx,ny);
			dy = (phi[m] - phi[n])/2.0;
			dyy = phi[n] + phi[m] - 2.0 * phi[idx];
			gy = (g[m] - g[n])/2.0;
		}
		
		dy2 = dy*dy;
		dxy = 0;
		if ((i > 0) && (i < nx -1) && (j > 0) && (j < ny -1)) {
			dxy = (phi[getIndex(i-1,j-1,nx,ny)] + phi[getIndex(i+1,j+1,nx,ny)] 
			     - phi[getIndex(i+1,j-1,nx,ny)] - phi[getIndex(i-1,j+1,nx,ny)])/4.0;
		}  
		kappa = ((dxx + dyy)*(dx2+dy2) - (dx2 *dxx + dy2*dyy + 2*dx*dy*dxy))/pow((dx2+dy2+eps),1.5); //curvature.
		norm = sqrt(dx2 + dy2)+eps;
		// compute force 
		F[k] = (gx * dx + gy * dy)/norm + (mu + g[idx]) * kappa;
		if (fabs(F[k]) > Fm) {
			Fm = fabs(F[k]);
		}	
		++k;
		ll_step(Lz);
	}
	
	//rescale if Fm > 0.5.
	for (k=0;k<Lz->length;k++) {
		if (Fm > 0.5) {
			F[k] = 0.5 * F[k]/Fm;
		} 	
	}

	// Update the zero level set
  	ll_init(Lz);
	k = 0;
  	while(Lz->curr != NULL){
  		i = Lz->curr->x; j = Lz->curr->y; idx = Lz->curr->idx;
		phi[idx] += F[k];	
		if (phi[idx] > 0.5) {
			ll_push(Sp1, ll_remcurr(Lz));
		} else if (phi[idx] < -0.5) {
			ll_push(Sn1, ll_remcurr(Lz));
		} else {
			ll_step(Lz);
		}
		++k;
	}
	free(F);
	
	if (Fm < epsConverge) {
		printf("Converged. \n");
		break;
	}
	
	// Update -1 and +1 level sets
  	ll_init(Ln1);
  	while(Ln1->curr != NULL){
  		i = Ln1->curr->x; j = Ln1->curr->y; idx = Ln1->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = -1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 0) flag = 1;
			if (label[np[k]] >=0) {
				if (Fm < phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			ll_push(Sn2, ll_remcurr(Ln1));
		} else {
			phi[idx] = Fm - 1.0;
			if (phi[idx] >= -0.5) {
				ll_push(Sz, ll_remcurr(Ln1));
			} else if (phi[idx] < -1.5) {
				ll_push(Sn2, ll_remcurr(Ln1));
			} else {
				ll_step(Ln1);
			}
		}
	}

  	ll_init(Lp1);
  	while(Lp1->curr != NULL){
  		i = Lp1->curr->x; j = Lp1->curr->y; idx = Lp1->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = 1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 0) flag = 1;
			if (label[np[k]] <= 0) {
				if (Fm > phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			ll_push(Sp2, ll_remcurr(Lp1));
		} else {
			phi[idx] = Fm + 1.0;
			if (phi[idx] <= 0.5) {
				ll_push(Sz, ll_remcurr(Lp1));
			} else if (phi[idx] > 1.5) {
				ll_push(Sp2, ll_remcurr(Lp1));
			} else {
				ll_step(Lp1);
			}
		}
	}
	
	// Update -2 and +2 level sets
  	ll_init(Ln2);
  	while(Ln2->curr != NULL){
  		i = Ln2->curr->x; j = Ln2->curr->y; idx = Ln2->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = -1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == -1) flag = 1;
			if (label[np[k]] >= -1) {
				if (Fm < phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			phi[idx] = -3.0;
			label[idx] = -3.0;
			ll_remcurr_free(Ln2);
		} else {
			phi[idx] = Fm - 1.0;
			if (phi[idx] >= -1.5) {
				ll_push(Sn1, ll_remcurr(Ln2));
			} else if (phi[idx] < -2.5) {
				phi[idx] = -3.0;
				label[idx] = -3.0;
				ll_remcurr_free(Ln2);
			} else {
				ll_step(Ln2);
			}
		}
	}

  	ll_init(Lp2);
  	while(Lp2->curr != NULL){
  		i = Lp2->curr->x; j = Lp2->curr->y; idx = Lp2->curr->idx;
		getNp(i,j,nx,ny,np);
		flag = 0;
		Fm = 1e10;
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (label[np[k]] == 1) flag = 1;
			if (label[np[k]] <= 1) {
				if (Fm > phi[np[k]]) Fm = phi[np[k]];
			}
		}
		if (flag == 0) {
			phi[idx] = 3.0;
			label[idx] = 3.0;
			ll_remcurr_free(Lp2);
		} else {
			phi[idx] = Fm + 1.0;
			if (phi[idx] <= 1.5) {
				ll_push(Sp1, ll_remcurr(Lp2));
			} else if (phi[idx] > 2.5) {
				phi[idx] = 3.0;
				label[idx] = 3.0;
				ll_remcurr_free(Lp2);
			} else {
				ll_step(Lp2);
			}
		}
	}
	
    // Move points into zero level set.
    ll_init(Sz);
	while (Sz->curr != NULL){
    	idx= Sz->curr->idx;
    	ll_push(Lz,ll_remcurr(Sz));
    	label[idx] = 0;
  	}
  	
  	// Move points into -1 and +1 level sets, and ensure -2, +2 neighbors
    ll_init(Sn1);
	while (Sn1->curr != NULL){
  		i = Sn1->curr->x; j = Sn1->curr->y; idx = Sn1->curr->idx;
  		label[idx] = -1;
  		ll_push(Ln1,ll_remcurr(Sn1));
		getNp(i,j,nx,ny,np);
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (phi[np[k]] == -3.0) {
				phi[np[k]] = phi[idx] - 1.0;
				ll_pushnew(Sn2,i+offx[k],j+offy[k],np[k]);
			}
		}
  	}

    ll_init(Sp1);
	while (Sp1->curr != NULL){
  		i = Sp1->curr->x; j = Sp1->curr->y; idx = Sp1->curr->idx;
  		label[idx] = 1;
  		ll_push(Lp1,ll_remcurr(Sp1));
		getNp(i,j,nx,ny,np);
		for (k=0;k<4;k++) {
			if (np[k] == -1) continue;
			if (phi[np[k]] == 3.0) {
				phi[np[k]] = phi[idx] + 1.0;
				ll_pushnew(Sp2,i+offx[k],j+offy[k],np[k]);
			}
		}
  	}
 
 	//Move points into -2 and +2 level sets
    ll_init(Sn2);
	while (Sn2->curr != NULL){
  		idx = Sn2->curr->idx;
  		label[idx] = -2;
  		ll_push(Ln2,ll_remcurr(Sn2));
	}

    ll_init(Sp2);
	while (Sp2->curr != NULL){
  		idx = Sp2->curr->idx;
  		label[idx] = 2;
  		ll_push(Lp2,ll_remcurr(Sp2));
	}	
  }	
  
  for (i=0; i<ntot; i++) {
	if (phi[i] > 0) {
		bw[i] = 1.0;
	} else {
		bw[i] = 0.0;
	}
  }

  free(label);
  free(phi);

  //destroy linked lists
  ll_destroy(Lz);
  ll_destroy(Ln1);
  ll_destroy(Ln2);
  ll_destroy(Lp1);
  ll_destroy(Lp2);

  ll_destroy(Sz);
  ll_destroy(Sn1);
  ll_destroy(Sn2);
  ll_destroy(Sp1);
  ll_destroy(Sp2);
  return;
}


//compar function for sorting, for inverse sorting, used in getSparseThreshold.
int comp (const void * elem1, const void * elem2) 
{
    float f = *((float*)elem1);
    float s = *((float*)elem2);
    if (f > s) return  -1;
    if (f < s) return 1;
    return 0;
}

float getSparseThreshold(int nx, int ny,float *bw,float sparse)
{
	int i,j,nxy, nn=100,mj;
	float *Is, thrs[100],p;
	
	nxy = nx * ny;
	Is = (float *) malloc(nxy * sizeof(float));
	for (i=0; i<nxy; i++) Is[i] = bw[i];
	qsort (Is, nxy, sizeof(float), comp);
	for (i=0; i<nn; i++) {
		thrs[i] = Is[0] - (Is[0] - Is[nxy-1])/nn*(i+1);
	}
	j = 0;
	mj = 0;
	for (i=0; i<nxy; i++) {
		if (Is[i] < thrs[j]) {
			p = i*1.0/nxy;
			if (p > sparse) {
				mj = j;
				break;
			}
			j++;
		}
	}
	free(Is);
	return thrs[mj];
}	

float getSparseThresholdAdaptive(int nx, int ny,float *bw,float sparseUp)
{
	int i,j,nxy, nn=100,mj,mu, flag;
	float thr, *Is, thrs[100], p[100],ps[100],mmax,dd[100];
	
	nxy = nx * ny;
	Is = (float *) malloc(nxy * sizeof(float));
	for (i=0; i<nxy; i++) Is[i] = bw[i];
	qsort (Is, nxy, sizeof(float), comp);
	for (i=0; i<nn; i++) {
		thrs[i] = Is[0] - (Is[0] - Is[nxy-1])/nn*(i+1);
	}
	for (i=0; i<nn; i++) {
		p[i] = 0.0;
		ps[i] = 0.0;
	}
	j = 0;
	mu = -1;
	for (i=0; i<nxy; i++) {
		if (Is[i] < thrs[j]) {
			p[j] = i*1.0/nxy;
			if (mu==-1 & p[j] > sparseUp) mu = i;
			j++;
		}
	}
	//smooth
	gaussianFilter2D(1,nn,1,30,p,ps,10.0);
	//get second derivative.
	for (i=1;i<nn;i++) {
		dd[i] = (ps[i+1] + ps[i-1] - 2 * ps[i])/2.0;
	}
	//smooth
	gaussianFilter2D(1,nn,1,10,dd,ps,2.0);
	mmax = -FLT_MAX;
	flag = 1;
	for (i=1;i<nn;i++) {
		if (flag == 1 & ps[i] > mmax) {
			mmax = ps[i];
		} else {
			mj = i;
			break;
		}
	}		
	if (p[mj] > sparseUp) {
		thr = thrs[mu];
	} else {
		thr = thrs[mj];
	}
	free(Is);
	return thr;
}

//binary image operations. 
//labeling connected pixels in binary image using two-pass algorithm based on Wikipedia description.
//http://en.wikipedia.org/wiki/Connected-component_labeling
int labelObjectsBinaryImage(int nx, int ny, float *bw, float *labels)
{
	int i,j,ii,jj,nxy,lb,k,L[8],nl,ml;
	Node *nodes; 				//nodes represeting each pixel's label. 
	int *LB;
	
	nxy = nx * ny; 						
	nodes = (Node *)malloc(nxy * sizeof(struct Node));
	MakeSet(nxy,nodes);			//make initia disjoint trees. 
	
	for (i=0; i<nxy; i++) labels[i] = 0.0;	//initialize. 
	//First pass
	lb = 0;
	for (i=0; i<nx; i++) {
	   for (j=0; j<ny; j++) {
		   ii = getIndex(i,j,nx,ny);
		   if (bw[ii] > 0.0) {
			   //check neighbors. 
			   nl = 0;
			   ml = nxy + 1;
			   if (j>0) {
				   jj = getIndex(i,j-1,nx,ny);
				   if (bw[jj] > 0.0) {
					   L[nl] = jj;
					   nl++;
					   if (ml > labels[jj]) ml = labels[jj];
				   }
				}
				if (i>0) {				   
				   jj = getIndex(i-1,j,nx,ny);
				   if (bw[jj] > 0.0) {
					   L[nl] = jj;
					   nl++;
					   if (ml > labels[jj]) ml = labels[jj];
				   }
				}
				if (i > 0 && j > 0) {
				   jj = getIndex(i-1,j-1,nx,ny);
				   if (bw[jj] > 0.0) {
					   L[nl] = jj;
					   if (ml > labels[jj]) ml = labels[jj];
					   nl++;
				   }
				}	
				if (i > 0 && j < ny-1) {
				   jj = getIndex(i-1,j+1,nx,ny);
				   if (bw[jj] > 0.0) {
					   L[nl] = jj;
					   if (ml > labels[jj]) ml = labels[jj];
					   nl++;
				   }
				}	
				
				if (nl ==0) {               //neighbors is empty, new label.
					lb += 1;
					labels[ii] = lb;
					(nodes+ii)->data = lb;
				} else {
					labels[ii] = ml;
					(nodes+ii)->data = ml;
					for (k=0;k<nl;k++) {
						Union(nodes+ii,nodes+L[k]);
					}
				}	
			}
		}	
	}		
  
	LB = (int *) malloc(lb * sizeof(int));
    //Second pass
    nl = 0;
    for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			ii = getIndex(i,j,nx,ny);
			if (bw[ii] > 0.0) {
				lb = Find(nodes+ii)->data;
				//check if this is a new label.
				if (nl == 0) {
					nl = 1;
					LB[0] = lb;
					k = 1;
				} else {
					for (k=0; k<nl; k++) {
						if (LB[k] == lb) {
							break;
						}
					}
					if (k==nl) {
						nl += 1;
						LB[nl-1] = lb;
						k = nl;
					} else {
						k = k+1;
					} 	
				}
				labels[ii] = k;
			}
		}
	}	
	
    free(nodes);
    free(LB);
    return nl;
}

//skeletonization. 
//Modified from the following source. 
/*
* Code for thinning a binary image using Zhang-Suen algorithm.
* Author: Nash (nash [at] opencv-code [dot] com)
* Website: http://opencv-code.com
*/
/*
* Perform one thinning iteration.
* Normally you wouldn't call this function directly from your code.
*
* Parameters:
* im Binary image with range = [0,1]
* iter 0=even, 1=odd
*/
void thinningIteration(int nx, int ny, float *img, int iter)
{
	int x, y, i, nxy, A, B, m1, m2;
	float *pAbove;
	float *pCurr;
	float *pBelow;
	float *nw, *no, *ne; // north (pAbove)
	float *we, *me, *ea;
	float *sw, *so, *se; // south (pBelow)
	float *pDst;
	float *marker;

	nxy = nx * ny;
	marker = (float *) malloc(nxy * sizeof(float));
	for (i=0; i<nxy; i++) marker[i] = 0.0;
	
	// initialize row pointers
	pAbove = NULL;
	pCurr = img;
	pBelow = img+ny;

	for (x = 1; x < nx-1; x++) {
		// shift the rows up by one
		pAbove = pCurr;
		pCurr = pBelow;
		pBelow = img + (x+1)*ny;
		pDst = marker + x*ny;
		// initialize col pointers
		no = &(pAbove[0]);
		ne = &(pAbove[1]);
		me = &(pCurr[0]);
		ea = &(pCurr[1]);
		so = &(pBelow[0]);
		se = &(pBelow[1]);

		for (y = 1; y < ny-1; y++) {
			// shift col pointers left by one (scan left to right)
			nw = no;
			no = ne;
			ne = &(pAbove[y+1]);
			we = me;
			me = ea;
			ea = &(pCurr[y+1]);
			sw = so;
			so = se;
			se = &(pBelow[y+1]);
			A = (*no == 0 && *ne == 1) + (*ne == 0 && *ea == 1) +
				(*ea == 0 && *se == 1) + (*se == 0 && *so == 1) +
				(*so == 0 && *sw == 1) + (*sw == 0 && *we == 1) +
				(*we == 0 && *nw == 1) + (*nw == 0 && *no == 1);
			B = *no + *ne + *ea + *se + *so + *sw + *we + *nw;
			m1 = iter == 0 ? (*no * *ea * *so) : (*no * *ea * *we);
			m2 = iter == 0 ? (*ea * *so * *we) : (*no * *so * *we);
			if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0) pDst[y] = 1;
		}
	}
	for (i=0; i<nxy; i++) {
		if (marker[i] > 0.0) img[i] = 0.0;
	}		
	free(marker);
}

/*
* Function for thinning the given binary image
*
* Parameters:
* src The source image, binary with range = [0,1.0]
* dst The destination image
*/
void skeletonizeZS(int nx, int ny, float *src, float *dst)
{
	int nxy, i, flag;
	float *prev, diff;
	
	nxy = nx * ny;
	prev = (float *) malloc(nxy * sizeof(float));	
	for (i=0; i<nxy; i++) {
		dst[i] = src[i];
		prev[i] = 0.0;
	}	
	while (1) {
		thinningIteration(nx, ny, dst, 0);
		thinningIteration(nx, ny, dst, 1);
		flag = 0;
		for (i=0; i<nxy; i++) {
			diff = abs(dst[i] - prev[i]);
			if (diff > 0.0) {
				flag = 1;
				break;
			}
		}
		for (i=0; i<nxy; i++) prev[i] = dst[i];
		if (flag == 0) break;
	}
	free(prev);
}

//prune small branches in a skeleton. 
void pruneSmallBranches(int nx, int ny, float *bw, float smallLen)
{
	int i,j,i1,j1,i2,j2,j3,ii,jj,kk,nxy,np,ne,iRemoved,nnp,jjp[8],fl;
	float *flags;
	int *x, *y, *endPoints;
	
	nxy = nx * ny;
	flags = (float *) malloc(nxy * sizeof(float));
	x = (int *) malloc(nxy * sizeof(int));
	y = (int *) malloc(nxy * sizeof(int));
	endPoints = (int *) malloc(nxy * sizeof(int));
	while(1) {
		iRemoved = 0;
		ne = 0;
		//get end points.
		for (ii=0; ii<nxy; ii++) {
			flags[ii] = bw[ii];
			if (bw[ii] > 0.0) {
				i = ii/ny;
				j = ii - i*ny;
				nnp = 0;		
				for (i1 = i-1; i1 < i+2; i1++) {
				for (j1 = j-1; j1 < j+2; j1++) {
					jj = getIndex(i1,j1,nx,ny);
					if (jj == ii || jj < 0 || jj >nxy-1) continue;
					if (bw[jj] > 0.0) {
						jjp[nnp++] = jj;
					}
				}}
				fl = 0;
				if (nnp <= 1) {
					fl = 1;
				} else if (nnp == 2) {
					i1 = jjp[0]/ny; j1 = jjp[0] - i1 * ny;
					i2 = jjp[1]/ny; j2 = jjp[1] - i2 * ny;
					if (sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2)) <= sqrt(2.0)) fl = 1;
				}
				if (fl == 1) {
					endPoints[ne++] = ii;
				}
			}
		}
		// construct path from the end points to the next
		for (j3=0;j3<ne; j3++) {
			ii = endPoints[j3];
			if (bw[ii] == 0.0) {//this point is already erased. 
				continue;
			}
			np = 0;
			i = ii/ny; j = ii - i*ny;
			x[np] = i; y[np] = j; np++;
			while (1) { //follow the skeleton.
				flags[ii] = 0;
				//select next point
				nnp = 0;		
				for (i1 = i-1; i1 < i+2; i1++) {
				for (j1 = j-1; j1 < j+2; j1++) {
					jj = getIndex(i1,j1,nx,ny);
					if (jj == ii || jj < 0 || jj >nxy-1) continue;
					if (flags[jj] > 0.0) {
						jjp[nnp++] = jj;
					}
				}}
				fl = 0;
				if (nnp == 1) {
					kk = jjp[0];	//this is the next point. 
					fl = 1;
				} else if (nnp == 2) {
					i1 = jjp[0]/ny; j1 = jjp[0] - i1 * ny;
					i2 = jjp[1]/ny; j2 = jjp[1] - i2 * ny;
					if (sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2)) == 1.0) {
						fl = 1;
						if (i1 == i || j1 == j) {
							kk = jjp[0];
						} else if (i2 == i || j2 == j) {
							kk = jjp[1];
						}
					}
				}
				if (fl == 1) {	//get to the next point. 
					ii = kk;
					i = ii/ny; j = ii - i*ny;
					x[np] = i; y[np] = j; np++;
					if (np > smallLen) break;	// this is a long branch. 
				} else {	//end point or cross point, stop. 
					break;
				}
			}
			if (np > 0 && np < smallLen) { // erase this branch. 
				iRemoved = 1;
				for (i1=0; i1<np; i1++) {
					jj = getIndex(x[i1],y[i1],nx,ny);
					bw[jj] = 0.0;
				}
			}
		}
		if (iRemoved == 0) break;
	}
	free(flags);
	free(x);
	free(y);
	free(endPoints);
}

//functions for autotracing neuron, NeuTuAutoTrace

//1d index of 3d index (i,j,k), used for dealing with tif image read in. 
long int getIndex3DZfirst(int i, int j, int k, int nx, int ny, int nz)
{
	return k * nx * ny + i * ny + j;
}

//set default parameters. Useful when testing functions called from python. 
void setDefaultParameters()
{
	printf("WARNING: Using default parameters. CHECK IF setDefaultParameters SHOULD BE CALLED!!!\n");
	//Parameter for image type. 
	PARAMS.imageType = 0;	//0 for bright field, 1 for dark field. 
	//Parameters for speciying pixel dimensions in the tiff stack, magnifications. These are determined when taking the images. 
	PARAMS.xyDist = 0.065; // distance between pixels in xy plane
	PARAMS.zDist = 0.5;  //distance in z between planes. 

	//Parameters for auto trace. Setting up
	PARAMS.nSplit = 8; //IMPORTANT, number of subdivision of a tiff stack in z (sub-slabs). 
	PARAMS.maxFracTotPoints = 0.2; //the maximum number of SWC points as a fraction of the pixel numbers in a stack. 
	PARAMS.zext = 10; //the extension of the z dimension in the sub-slabs. Useful for avoiding overlapping SWC points due to subdivisions. 

	//Parameters for creating 2D mask from 2D projections from the slabs. 
	PARAMS.sigmaFilter = 2.0;	//sigmaa for the valley detector filter. Typically set to the half width of the smallest widths of neurites. Increase if the image is noisy. 
	PARAMS.sparse = 0.2; // range from 0 to 1. IMPOARTANT, upper limit of sparseness in determining the threshold with Hessian matrix. 
	PARAMS.sigmaBack = 100.0; //length scale sigma for Gaussian smoothing to get the smooth background. 
	PARAMS.smallArea = 100.0; //area of small patches to be removed in the mask. Useful for removing noise. 
	PARAMS.lambdaRatioThr = 2.0; //ranges from 1 to 10, threshold for the ratio lambda1/lambda2 in Hessian to supress blob. 
	PARAMS.levelSetIter = 1500; //number of iterations for levelset smoothing. 
	PARAMS.levelSetMu = 0.0; //larger than 0, levelset parameter for controlling the smoothness of the mask. 

	//Parameters for skeletonization from the mask. 
	PARAMS.smallLen = 20;  //length smaller than this in skeleton are pruned. 

	//Parameters for getting linked points, which is the basis for generating the SWC file.
	PARAMS.sigmaSmooth = 5.0; //for Gaussian smoothing the skeleton segments when creating the linked point segments. 	PARAMS.zJumpFact = 2.0; //IMPORTANT, factor for determing whether two points should be disconnected because of too large z-jump. 
	PARAMS.distFactor = 2.0; //IMPORTNAT, factor for disconnecting two points. If the 2D distance between two points is larger than 
	PARAMS.angle = 60.0; //in degrees. If the turn of the new points is more than this angle, the new point is disconnected. 
	PARAMS.drFactor = 0.1; //greater than 0, a factor for reducing the distance between linked points for large radius. 

	//Parameters for growing from end points to fill in and bridges gaps. 
	PARAMS.factSearchMax = 3.0; //this factor multiplied by the radius is the maximum range for searching for the next point.
	PARAMS.factSearchMin = 2.5; //this factor multiplied by the radius is the minimum range for searching for the next point.
	PARAMS.zBoxSize = 5; //the z dimension of a volume around an end point for searching the next linked point for growing, 

	//Parameters for determining the validity of the linked points. These are crucial parameters for 
	//judging whether a linked point is at the neurite, and whether it has the correct radius. Also used to correct the positions of 
	//the linked points. 
	PARAMS.minRange = 50; //minimum range for computing the extent of the inverse peak. Typically set at 4*r if this is bigger than the minRange. 
	PARAMS.sigmaSmoothCurve = 5.0; //for smoothing the inverse peak. Imporant for judging the position and width of the peak. 
	PARAMS.sparseThreshold = 0.2; //range from 0 to 1, the fraction of top values in the curves. This is used  to determine the starting point of the background judging the validity of the linked point. Increasing this value makes the judgement more strict. 
	PARAMS.factSigmaThreshold = 8.0; //factor for determing the threshold for the background from the inverse peak curves. 
	PARAMS.factSmallDeriv = 0.5; //factor for determining the smallness in changes in the derivative. The threshold is set at this factor multiplied by sigma from the curves.  
	PARAMS.minRadius = 2; //minimum radius of the linked point. Points smaller than this are considered noise. 
	PARAMS.maxRadius = 100; //maximum radius of the linked point. Points smaller than this are considered noise. 
	PARAMS.factShift = 1.0; //this factor multiplied by r is the maximum shift allowed in the xy position of the linked point. 
	PARAMS.factAdjustRadius = 0.8; //factor for adjusting the radius measured. Often the width of the inverse peak can be consistently 

	//Parameters for marking occupied pixels. IMPORTANT!
	PARAMS.factMarkOccXY = 1.0; //pixels around an accepted linked point are marked occupied. This factor times the radius is the range in xy.
	PARAMS.zOcc = 10; //pixels around an accepted linked point are marked occupied. This is the z extent. 

	//Parameters for connecting end points. 
	PARAMS.distFactConn = 5.0; //sets the minimum distance for connecting points, distFactConn *(r1 + r2), where r1, r2 are the radii of the points to be connected.  

	//Parameters for creating SWC structure. 
	PARAMS.minNumPointsBr = 3; //minimum number of points in a branch for accepting. 
	PARAMS.minNumPointsBrIso = 10; //minimum number of points in an isolated branch for accepting. 

	PARAMS.zfact = PARAMS.zDist/PARAMS.xyDist; //ratio between z distance and the xy pixel distance, zDist/xyDist. 
	
	PARAMS.somaLengthScale = 50;			//length scale of the soma and thick dendrites.
	PARAMS.somaSparseThr = 0.05;			//fraction of top darkest pixels for getting the mask for the soma.  

	PARAMS.alphaDistance = 20.0;			//factor for weighted shortest distance. Used in computing left right shortest path and in extending from ends. 
}

void readParametersFromFile(char *filenameParam)
{
	int ind;
	size_t len;
	char *line;
	FILE *fpt;
	ssize_t read;	
	float numIn;

	printf("Reading parameters from file %s\n",filenameParam);
	fpt = fopen(filenameParam,"r");
	ind = 0;
	line = 0;
	len = 0;
	while ((read = getline(&line, &len, fpt)) != -1) {
		if (len > 0 && line[0] == '#') {	//this is a comment.
			continue;
		}
		if (len > 0 && line[0] == '\n') {	//this is an empty line
			continue;
		}
		if (len == 0) continue;
		sscanf(line,"%f",&numIn); 
		ind++;	
		if (ind == 1) {
			PARAMS.imageType = (int) numIn; // image type, 0 for bright field, 1 for dark field. 
		} else if (ind == 2) {
			PARAMS.xyDist = numIn; // distance between pixels in xy plane
		} else if (ind == 3) {
			PARAMS.zDist = numIn;  //distance in z between planes. 
		} else if (ind == 4) {
			PARAMS.nSplit = (int) numIn; //IMPORTANT, number of subdivision of a tiff stack in z (sub-slabs). 
		} else if (ind == 5) {	
			PARAMS.maxFracTotPoints = numIn; //the maximum number of SWC points as a fraction of the pixel numbers in a stack. 
		} else if (ind == 6) {
			PARAMS.zext = (int)fmax(1,numIn/PARAMS.zDist); //the extension of the z dimension in the sub-slabs. Useful for avoiding overlapping SWC points due to subdivisions. 
		} else if (ind == 7) {
			PARAMS.sigmaFilter = numIn/PARAMS.xyDist; //sigmaa for the valley detector filter. Typically set to the half width of the smallest widths of neurites. Increase if the image is noisy. 
		} else if (ind == 8) {	
			PARAMS.sparse = numIn; // range from 0 to 1. IMPOARTANT, upper limit of sparseness in determining the threshold with Hessian matrix. 
		} else if (ind == 9) {
			PARAMS.sigmaBack = numIn/PARAMS.xyDist; //length scale sigma for Gaussian smoothing to get the smooth background. 
		} else if (ind == 10) {
			PARAMS.smallArea = (int)numIn/(PARAMS.xyDist * PARAMS.xyDist); //area of small patches to be removed in the mask. Useful for removing noise. 
		} else if (ind == 11) {
			PARAMS.lambdaRatioThr = numIn; //ranges from 1 to 10, threshold for the ratio lambda1/lambda2 in Hessian to supress blob. 
		} else if (ind == 12) {
			PARAMS.levelSetIter = (int) numIn; //number of iterations for levelset smoothing. 
		} else if (ind == 13) {
			PARAMS.levelSetMu = numIn; //larger than 0, levelset parameter for controlling the smoothness of the mask. 
		} else if (ind == 14) {
			PARAMS.smallLen = (int) numIn/PARAMS.xyDist;  //length smaller than this in skeleton are pruned. 
		} else if (ind == 15) {
			PARAMS.sigmaSmooth = numIn/PARAMS.xyDist;  //for Gaussian smoothing the skeleton segments when creating the linked point segments. 	
		} else if (ind == 16) {	
			PARAMS.zJumpFact = numIn; //IMPORTANT, factor for determing whether two points should be disconnected because of too large z-jump. 
		} else if (ind == 17) {
			PARAMS.distFactor = numIn; //IMPORTNAT, factor for disconnecting two points. If the 2D distance between two points is larger than 
		} else if (ind == 18) {
			PARAMS.angle = numIn; //in degrees. If the turn of the new points is more than this angle, the new point is disconnected. 
		} else if (ind == 19) {
			PARAMS.drFactor = numIn; //greater than 0, a factor for reducing the distance between linked points for large radius. 
		} else if (ind == 20) {
			PARAMS.factSearchMax = numIn; //this factor multiplied by the radius is the maximum range for searching for the next point.
		} else if (ind == 21) {
			PARAMS.factSearchMin = numIn; //this factor multiplied by the radius is the minimum range for searching for the next point.
		} else if (ind == 22) {
			PARAMS.zBoxSize = (int) fmax(1,numIn/PARAMS.zDist); //the z dimension of a volume around an end point for searching the next linked point for growing, 
		} else if (ind == 23) {
			PARAMS.minRange = (int) fmax(1,numIn/PARAMS.xyDist); //minimum range for computing the extent of the inverse peak. Typically set at 4*r if this is bigger than the minRange. 
		} else if (ind == 24) {
			PARAMS.sigmaSmoothCurve = numIn/PARAMS.xyDist; //for smoothing the inverse peak. Imporant for judging the position and width of the peak. 
		} else if (ind == 25) {
			PARAMS.sparseThreshold = numIn; //range from 0 to 1, the fraction of top values in the curves. This is used  to determine the starting point of the background judging the validity of the linked point. Increasing this value makes the judgement more strict. 
		} else if (ind == 26) {
			PARAMS.factSigmaThreshold = numIn; //factor for determing the threshold for the background from the inverse peak curves. 
		} else if (ind == 27) {
			PARAMS.factSmallDeriv = numIn; //factor for determining the smallness in changes in the derivative. The threshold is set at this factor multiplied by sigma from the curves.  
		} else if (ind == 28) {
			PARAMS.minRadius = (int) fmax(1,numIn/PARAMS.xyDist); //minimum radius of the linked point. Points smaller than this are considered noise. 
		} else if (ind == 29) {
			PARAMS.maxRadius = (int) fmax(1,numIn/PARAMS.xyDist); //maximum radius of the linked point. Points smaller than this are considered noise. 
		} else if (ind == 30) {
			PARAMS.factShift = numIn; //this factor multiplied by r is the maximum shift allowed in the xy position of the linked point. 
		} else if (ind == 31) {
			PARAMS.factAdjustRadius = numIn; //factor for adjusting the radius measured. Often the width of the inverse peak can be consistently 
		} else if (ind == 32) {
			PARAMS.factMarkOccXY = numIn; //pixels around an accepted linked point are marked occupied. This factor times the radius is the range in xy.
		} else if (ind == 33) {
			PARAMS.zOcc = (int) fmax(1,numIn/PARAMS.zDist); //pixels around an accepted linked point are marked occupied. This is the z extent. 
		} else if (ind == 34) {
			PARAMS.distFactConn = numIn; //sets the minimum distance for connecting points, distFactConn *(r1 + r2), where r1, r2 are the radii of the points to be connected.  
		} else if (ind == 35) {
			PARAMS.minNumPointsBr = (int) numIn; //minimum number of points in a branch for accepting. 
		} else if (ind == 36) {
			PARAMS.minNumPointsBrIso = (int) numIn; //minimum number of points in an isolated branch for accepting. 
		} else if (ind == 37) {
			PARAMS.somaLengthScale = (int) fmax(1,numIn/PARAMS.xyDist); //length scale of the soma. 
		} else if (ind == 38) {
			PARAMS.somaSparseThr = numIn; //fraction of top darkest pixels for getting the mask for the soma. 
		} else if (ind == 39) {
			PARAMS.alphaDistance = numIn; //fraction of top darkest pixels for getting the mask for the soma. 
		} 
	}
	if (line) free(line);
	fclose(fpt);	
	PARAMS.zfact = PARAMS.zDist/PARAMS.xyDist; //ratio between z distance and the xy pixel distance, zDist/xyDist. 	
}

void printParameters()
{
	//print parameters. 
	printf("The parameters are:\n");
	printf ("1. xyDist		=%6.3f micron\n",PARAMS.xyDist);
	printf ("2. zDist		=%6.3f micron\n",PARAMS.zDist); 
	printf ("3. nSplit		=%d\n",PARAMS.nSplit);
	printf ("4. maxFracTotPoints	=%6.3f\n",PARAMS.maxFracTotPoints);
	printf ("5. zext			=%d plane\n",PARAMS.zext);
	printf ("6. sigmaFilter		=%6.3f pixel\n",PARAMS.sigmaFilter);
	printf ("7. sparse		=%6.3f\n",PARAMS.sparse);	
	printf ("8. sigmaBack		=%6.3f pixel\n",PARAMS.sigmaBack);
	printf ("9. smallArea		=%d pixel^2\n",PARAMS.smallArea);
	printf ("10. lambdaRatioThr	=%6.3f\n",PARAMS.lambdaRatioThr);
	printf ("11. levelSetIter	=%d\n",PARAMS.levelSetIter);
	printf ("12. levelSetMu		=%6.3f\n",PARAMS.levelSetMu);
	printf ("13. smallLen		=%d pixel\n",PARAMS.smallLen);
	printf ("14. sigmaSmooth		=%6.3f pixel\n",PARAMS.sigmaSmooth);
	printf ("15. zJumpFact		=%6.3f\n",PARAMS.zJumpFact);
	printf ("16. distFactor		=%6.3f\n",PARAMS.distFactor);
	printf ("17. angle		=%6.3f degree\n",PARAMS.angle);
	printf ("18. drFactor		=%6.3f\n",PARAMS.drFactor);
	printf ("19. factSearchMax	=%6.3f\n",PARAMS.factSearchMax);
	printf ("20. factSearchMin	=%6.3f\n",PARAMS.factSearchMin);
	printf ("21. zBoxSize		=%d plane\n",PARAMS.zBoxSize);
	printf ("22. minRange		=%d pixel\n",PARAMS.minRange);
	printf ("23. sigmaSmoothCurve	=%6.3f pixel\n",PARAMS.sigmaSmoothCurve);
	printf ("24. sparseThreshold	=%6.3f\n",PARAMS.sparseThreshold);
	printf ("25. factSigmaThreshold	=%6.3f\n",PARAMS.factSigmaThreshold);
	printf ("26. factSmallDeriv	=%6.3f\n",PARAMS.factSmallDeriv);
	printf ("27. minRadius		=%d pixel\n",PARAMS.minRadius);
	printf ("28. maxRadius		=%d pixel\n",PARAMS.maxRadius);
	printf ("29. factShift		=%6.3f\n",PARAMS.factShift);
	printf ("30. factAdjustRadius	=%6.3f\n",PARAMS.factAdjustRadius); 
	printf ("31. factMarkOccXY	=%6.3f\n",PARAMS.factMarkOccXY);
	printf ("32. zOcc		=%d plane\n",PARAMS.zOcc);
	printf ("33. distFactConn	=%6.3f\n",PARAMS.distFactConn);
	printf ("34. minNumPointsBr	=%d\n",PARAMS.minNumPointsBr);	
	printf ("35. minNumPointsBrIso	=%d\n",PARAMS.minNumPointsBrIso);	
	printf ("36. somaLengthScale	=%d pixel\n",PARAMS.somaLengthScale);	
	printf ("37. somaSparseThr	=%f\n",PARAMS.somaSparseThr);	
	printf ("38. alphaDistance	=%f\n",PARAMS.alphaDistance);	
	if (PARAMS.imageType == 0) {
		printf("Bright field image. \n");
	} else if (PARAMS.imageType == 1) {
		printf("Dark field image. \n");
	} else {
		printf("Image type unspecified. Will attempt to detect automatically. \n");
	}
}

#define createMask2D_DIAG 0						//set to zero unless for outputing images for diagnosis. 
#if createMask2D_DIAG
int counter = 0; char filenametemp[1000];	//for outputing images during creating maskes. 
#endif 
//create mask bw from image I.    
void createMask2D(int nx, int ny, float *I, float *bw)
{
	int i,j,nxy,ii,jj,kk,smallArea,i2;
	float mu, sigmap, sigmaf,dt,sigmaBack;
	float *Is, *Ix, *Iy, *Ixy, *v, *u, *vx, *uy, dx=1.0, dy=1.0;
	int Wxp, Wyp, Wxf, Wyf, WxBack, WyBack;
	float sparse,thrSparse,lambdaRatioThr;
	int nthrs = 100,iter;
	float thrs[100];
	int *labelCounts, *delLabs;
	float mmin,mmax;
	int maxiter;
	
	#if createMask2D_DIAG	
	readParametersFromFile("NeuTu.Parameters.dat");	//set up parameters. 
	#endif

	printf("Creating mask....\n");

	//SET UP PARAMETERS. ONLY USED WHEN THE FUNCTION IS CALLED FROM PYTHON FOR TESTING. COMMENT THIS OUT IN PRODUCTON 
	//setDefaultParameters();
	
	sigmap = 1.0;	
	sigmaf = PARAMS.sigmaFilter;
	Wxp = (int) fmax(5,5*sigmap); 
	Wyp = Wxp;			
	Wxf = (int) fmax(5,5*sigmaf); 
	Wyf = Wxf;			
	sparse = PARAMS.sparse;			
	lambdaRatioThr = PARAMS.lambdaRatioThr;	
	smallArea = PARAMS.smallArea;	
	iter = PARAMS.levelSetIter;		
	mu = PARAMS.levelSetMu;
	sigmaBack  = PARAMS.sigmaBack;
	WxBack = (int) fmax(6, 5*sigmaBack);
	WyBack = WxBack;
		
	nxy = nx * ny;
	Is = (float *)malloc(nxy * sizeof(float));
	Ix = (float *)malloc(nxy * sizeof(float));
	Iy = (float *)malloc(nxy * sizeof(float));
	Ixy = (float *)malloc(nxy * sizeof(float));
	u = (float *)malloc(nxy * sizeof(float));
	v = (float *)malloc(nxy * sizeof(float));
	uy = (float *)malloc(nxy * sizeof(float));
	vx = (float *)malloc(nxy * sizeof(float));
	
	//smooth
	gaussianFilter2D(nx, ny, Wxp, Wyp, I, Is, sigmap);	

	//get the background. 
	gaussianFilter2D(nx, ny, WxBack, WyBack, I, Ixy, sigmaBack);
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.proj.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,Is,filenametemp);
	sprintf(filenametemp,"temp.background.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,Ixy,filenametemp);
	#endif 	
	
	//get rid of background and normalize I. 
	mmin = FLT_MAX;
	mmax = FLT_MIN;
	for (i=0; i<nxy; i++) {
		Is[i] -= Ixy[i];
		if (Is[i] > mmax) mmax = Is[i];
		if (Is[i] < mmin) mmin = Is[i];
	}
	for (i=0; i<nxy; i++) {
		Is[i] = (Is[i] - mmin)/(mmax - mmin + 1e-5);
	}
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.projAdjusted.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,Is,filenametemp);
	#endif 	
	
	// compute Ix and Iy. 
	for (i=1; i<nx-1; i++) {
		for (j=0; j<ny; j++) {
			ii = getIndex(i,  j,nx,ny);
			jj = getIndex(i+1,j,nx,ny);
			kk = getIndex(i-1,j,nx,ny);
			Ix[ii] = (Is[jj] - Is[kk])/(2.0*dx);
		}
	}
	for (j=0; j<ny; j++){	//boundary
		ii = getIndex(0,j,nx,ny);
		jj = getIndex(1,j,nx,ny);
		kk = getIndex(0,j,nx,ny);
		Ix[ii] = (Is[jj] - Is[kk])/(1.0*dx);
		ii = getIndex(nx-1,j,nx,ny);
		jj = getIndex(nx-1,j,nx,ny);
		kk = getIndex(nx-2,j,nx,ny);
		Ix[ii] = (Is[jj] - Is[kk])/(1.0*dx);
	}	
	for (j=1; j<ny-1; j++) {
		for (i=0; i<nx; i++) {
			ii = getIndex(i,j,  nx,ny);
			jj = getIndex(i,j+1,nx,ny);
			kk = getIndex(i,j-1,nx,ny);			
			Iy[ii] = (Is[jj] - Is[kk])/(2.0*dy);
		}
	}	
	for (i=0; i<nx; i++) {	//boundary.
		ii = getIndex(i,0,nx,ny);
		jj = getIndex(i,1,nx,ny);
		kk = getIndex(i,0,nx,ny);
		Iy[ii] = (Is[jj] - Is[kk])/(1.0*dy);
		ii = getIndex(i,ny-1,nx,ny);
		jj = getIndex(i,ny-1,nx,ny);
		kk = getIndex(i,ny-2,nx,ny);
		Iy[ii] = (Is[jj] - Is[kk])/(1.0*dy);
	}	

	//smooth the derivatives. 
	gaussianFilter2D(nx, ny, Wxp, Wyp, Ix, vx, sigmap);	
	gaussianFilter2D(nx, ny, Wxp, Wyp, Iy, uy, sigmap);	
			
	// calculate the Hessian matrix
	for (i=1; i<nx-1; i++) {
		for (j=0; j<ny; j++){
			ii = getIndex(i,  j,nx,ny);
			jj = getIndex(i+1,j,nx,ny);
			kk = getIndex(i-1,j,nx,ny);			
			Ix[ii] = (vx[jj] - vx[kk])/(2.0*dx);
			Ixy[ii] =(uy[jj] - uy[kk])/(2.0*dx)/2.0;
		}
	}
	for (j=0; j<ny; j++){	//boundary. 
		ii = getIndex(0,j,nx,ny);
		jj = getIndex(1,j,nx,ny);
		kk = getIndex(0,j,nx,ny);
		Ix[ii] = (vx[jj] - vx[kk])/(1.0*dx);
		Ixy[ii] = (uy[jj] - uy[kk])/(1.0*dx)/2.0;
		ii = getIndex(nx-1,j,nx,ny);
		jj = getIndex(nx-1,j,nx,ny);
		kk = getIndex(nx-2,j,nx,ny);
		Ix[ii] = (vx[jj] - vx[kk])/(1.0*dx);
		Ixy[ii] = (uy[jj] - uy[kk])/(1.0*dx)/2.0;
	}	
	for (j=1; j<ny-1; j++) {
		for (i=0; i<nx; i++) {
			ii = getIndex(i,j,  nx,ny);
			jj = getIndex(i,j+1,nx,ny);
			kk = getIndex(i,j-1,nx,ny);						
			Iy[ii] = (uy[jj] - uy[kk])/(2.0*dy);		
			Ixy[ii] += (vx[jj] - vx[kk])/(2.0*dy)/2.0; 		
		}
	}
	for (i=0; i<nx; i++) {	//boundary
		ii = getIndex(i,0,nx,ny);
		jj = getIndex(i,1,nx,ny);
		kk = getIndex(i,0,nx,ny);
		Iy[ii] = (uy[jj] - uy[kk])/(1.0*dy);
		Ixy[ii] += (vx[jj] - vx[kk])/(1.0*dy)/2.0; 		
		ii = getIndex(i,ny-1,nx,ny);
		jj = getIndex(i,ny-1,nx,ny);
		kk = getIndex(i,ny-2,nx,ny);
		Iy[ii] = (uy[jj] - uy[kk])/(1.0*dy);
		Ixy[ii] += (vx[jj] - vx[kk])/(1.0*dy)/2.0; 		
	}	
	
	//smooth the second derivatives with the filter. 	
	gaussianFilter2D(nx, ny, Wxf, Wyf, Ix, Ix, sigmaf);	
	gaussianFilter2D(nx, ny, Wxf, Wyf, Iy, Iy, sigmaf);	
	gaussianFilter2D(nx, ny, Wxf, Wyf, Ixy,Ixy, sigmaf);	
	
	//get rid of the boundary effects
	for (i=0; i<3*sigmaf+1; i++) {
	for (j=0; j<ny; j++){
		ii = getIndex(i,j,nx,ny);
		Ix[ii] = 0.0;
		Ixy[ii] = 0.0;
		vx[ii] = 0;
		ii = getIndex(nx-1-i,j,nx,ny);
		Ix[ii] = 0.0;
		Ixy[ii] = 0.0;
		vx[ii] = 0;
	}}	
	for (j=0; j<3*sigmaf; j++) {
	for (i=0; i<nx; i++) {
		ii = getIndex(i,j,nx,ny);
		Iy[ii] = 0.0;
		Ixy[ii] = 0.0;		
		uy[ii] = 0;
		ii = getIndex(i,ny-1-j,nx,ny);
		Iy[ii] = 0.0;
		Ixy[ii] = 0.0;		
		uy[ii] = 0;
	}}	
		
	// compute the eigen values of the Hessian matrix. 	
	for (i=0; i<nxy; i++) {
		v[i] = 0.5 * (Ix[i] + Iy[i] + sqrt((Ix[i] - Iy[i])*(Ix[i] - Iy[i]) + 4 * Ixy[i] * Ixy[i])); //lambda1
		u[i] = 0.5 * (Ix[i] + Iy[i] - sqrt((Ix[i] - Iy[i])*(Ix[i] - Iy[i]) + 4 * Ixy[i] * Ixy[i]));			//lambda2
	}
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.lambda.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,v,filenametemp);
	#endif 	
	gaussianFilter2D(nx, ny, Wxp, Wyp, v, v, sigmap);	
	gaussianFilter2D(nx, ny, Wxp, Wyp, u, u, sigmap);
	
	mmax = FLT_MIN;
	mmin = FLT_MAX;	
	for (i=0;  i<nxy; i++) {
		if (v[i] > lambdaRatioThr*abs(u[i])) {
			bw[i] = v[i];
		} else {
			bw[i] = 0.0;
		}
		if (mmax < bw[i]) mmax = bw[i];
	}		
	//get the threshold using the sparse critiria. 
	gaussianFilter2D(nx, ny, Wxp, Wyp, bw, bw, sigmap);
	thrSparse = getSparseThreshold(nx,ny,bw,sparse);
	//thrSparse = 0.6 * mmax;
	
	for (i=0;i<nxy;i++) {
		if (bw[i] > thrSparse) {
			bw[i] = 1.0;	
		} else {
			bw[i] = 0.0;
		}
	}
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.bwinit.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,bw,filenametemp);
	#endif 	
	
	//fill pixel sized holes. 
	for (i=2; i<nx-3; i++) {
		for (j=2; j<ny-2; j++) {
			i2 = getIndex(i,j,nx,ny);
			if (bw[i2] > 0) continue;
			kk = 1;
			for (ii=i-2; ii<i+3; ii++) {
				jj = getIndex(ii,j-2,nx,ny);
				if (bw[jj] == 0) {
					kk = 0;
					break;
				}
				jj = getIndex(ii,j+2,nx,ny);
				if (bw[jj] == 0) {
					kk = 0;
					break;
				}
			}
			if (kk==0) continue;
			kk = 1;
			for (ii=i-2; ii<i+3; ii++) {
				jj = getIndex(i-2,ii,nx,ny);
				if (bw[jj] == 0) {
					kk = 0;
					break;
				}
				jj = getIndex(i+2,ii,nx,ny);
				if (bw[jj] == 0) {
					kk = 0;
					break;
				}
			}
			if (kk==1) bw[i2] = 1;
		}
	}

	//levelset for boundary and smoothing
	//IMPORTANT:
	//no edge indicator! See vx[i] is set to zero. If edge indicator is needed, remove 0.0. 
	for (i=0; i<nxy; i++) {
		vx[i] = 0.0;
	}
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.edge.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,vx,filenametemp);
	#endif 	
	sparseFieldLevelSet(nx, ny, vx, bw, iter, mu);	

	//remove small areas. 
	printf("Removing small areas ...\n");
	kk = labelObjectsBinaryImage(nx, ny, bw, Is);
	labelCounts = (int *) malloc(kk * sizeof(int));
	for (i=0; i<kk; i++) labelCounts[i] = 0;
	delLabs = (int *) malloc(kk * sizeof(int));
	for (ii=0; ii<nxy; ii++) {
		if (Is[ii] < 1.0) continue;
		for (i=0; i<kk; i++) {
			if ((int)Is[ii] == i+1) labelCounts[i] += 1;
		}
	}
	jj = 0;
	for (i=0; i<kk; i++) {
		if (labelCounts[i] < smallArea) {
			delLabs[jj] = i+1;
			jj += 1;
		}
	}
	for (ii=0; ii<nxy; ii++) {
		if (Is[ii] < 1.0) continue;
		for (i=0; i<jj; i++) {
			if ((int)Is[ii] == delLabs[i]) {
				bw[ii] = 0.0;
				break;
			}
		}
	}
	
	//get rid of boundary effect.
	for (i=0;i<5; i++) {
		for (j=0; j<ny; j++) {
			bw[getIndex(i,j,nx,ny)] = 0.0;
			bw[getIndex(nx-1-i,j,nx,ny)] = 0.0;
		}
	}
	for (j=0;j<5; j++) {
		for (i=0; i<nx; i++) {
			bw[getIndex(i,j,nx,ny)] = 0.0;
			bw[getIndex(i,ny-1-j,nx,ny)] = 0.0;
		}
	}
	#if createMask2D_DIAG
	sprintf(filenametemp,"temp.bwfinale.%d.tif",counter);
	createGrayTiffStackFromFloatArray(nx,ny,1,bw,filenametemp);
	++counter;
	#endif 	
	
	free(Is);
	free(Ix); 
	free(Iy); 
	free(u);
	free(v);
	free(uy);
	free(vx);
	free(Ixy);
	free(labelCounts);
	free(delLabs);
}  

//create a linked point. 
LinkedPoint *CreateLinkedPoint(float x, float y, float z, float r, int ID, int Type)
{
	LinkedPoint *p;
	p = (LinkedPoint *)malloc(sizeof(LinkedPoint));
	p->x = x;
	p->y = y;
	p->z = z;
	p->r = r;
	p->ID = ID;
	p->Type = Type;
	p->conn = (LinkedList *) malloc(sizeof(LinkedList));
	p->conn->val =0;		//the first element in the conn list contains the number of elements in the list.
	p->conn->next = NULL;
	return p;
}


void DeleteLinkedPoint(LinkedPoint *p) 
{
	if (p == NULL) return;
	if (p->conn != NULL) {
		DeleteList(p->conn);
	}
	free(p);
	p = NULL;
}

void AddConnLinkedPoint(LinkedPoint *p, int ID)
{
	LinkedList *iter;
	int flag = 0;
	if (p==NULL) return;
	if (p->ID == ID) return;	//do not connect to itself. 
	//check if this ID is already connected. 
	iter = p->conn->next;
	while (iter != NULL) {
		if (iter->val == ID) {	//already connected. 
			flag = 1;
			break;
		}
		iter = iter->next;
	}
	if (flag == 0) {
		AppendToList(&(p->conn),ID);
		p->conn->val += 1;
	}
}

void DelConnLinkedPoint(LinkedPoint *p, int ID)
{
	if (p == NULL) return;
	LinkedList *cp,*cpp;
	cpp = p->conn;
	cp = p->conn->next;
	while (cp != NULL) {
		if (cp->val == ID) {
			cpp->next = cp->next;
			p->conn->val -= 1;
			free(cp);
			break;
		} else {
			cpp = cp;
			cp = cp->next;
		}
	}	
}

int NumConnLinkedPoint(LinkedPoint *p)
{	
	return p->conn->val;
}


//adjust z of a LinkedPoint. 
#define WSZ 1000	//maximum size of workspace for the z profile.
void adjustZOfLinkedPoint(int x, int y, int z, int r, int nx, int ny, int nz, float *im3d, int *zo) {
	float sw[WSZ],fw[WSZ],smoothfw[WSZ];
	int np,np0;
	int is;
	int zexc;
	int kk;
	int Wx, Wy;
	int sg;
	#if adjustZOfLinkedPoint_DIAG
	FILE *fpt;
	char filename[1000];
	sprintf(filename,"%s","NeuTu.Parameters.dat");
	readParametersFromFile(filename);
	#endif
	(*zo) = z;
	//get the z profile. 
	zexc = (int) fmax(PARAMS.zBoxSize,1.0*r/PARAMS.zfact);	//PARAMETER!
	np = 0;
	np0 = -1;
	for (is = -zexc; is < zexc; is++) {
		if (z + is >= 0 && z + is < nz) {
			sw[np] = is;
			fw[np] = im3d[getIndex3DZfirst(x,y,z+is,nx,ny,nz)];
			if (is == 0) np0 = np;
			np++;
		}
	}
	if (np0 == 0 || np0 == np) return; ////too close to the boundary. 
	if (np0 == -1) return;	//too close to the boundary. 
	//smooth. 
	Wx = 3;
	Wx = Wx/2*2 + 1;
	Wy = Wx;		
	gaussianFilter2D(np,1,Wx,Wy,fw,smoothfw,1);	//PARAMETER!!!
	//find the minimum. 
	sg = 0;
	if (smoothfw[np0] < smoothfw[np0-1]) sg = 1;
	if (smoothfw[np0] < smoothfw[np0+1]) sg = -1;
	kk = np0;
	if (sg !=0 ) {
		while (kk+sg >=0 && kk+sg < np) {
			if (smoothfw[kk+sg] < smoothfw[kk]) {
				kk += sg;
			} else {
				break;
			}
		}
	}
	(*zo) = z + kk - np0;	//adjust z. 

	#if adjustZOfLinkedPoint_DIAG
	//diagonose. copy paste into ipython 
	printf("Writing diagonostic information to temp3.py\n");
	fpt = fopen("temp3.py","w");
	fprintf(fpt,"from pylab import *\n");
	fprintf(fpt,"a=[");
	for (is=0;is<np-1; is++) {
		fprintf(fpt,"%3.2f, ",sw[is]);
	}
	fprintf(fpt,"%3.2f]\n",sw[np-1]);
	fprintf(fpt,"b=[");
	for (is=0;is<np-1; is++) {
		fprintf(fpt,"%3.2f, ",smoothfw[is]);
	}
	fprintf(fpt,"%3.2f]\n",smoothfw[np-1]);		
	fprintf(fpt,"c=[");
	for (is=0;is<np-1; is++) {
		fprintf(fpt,"%3.2f, ",fw[is]);
	}
	fprintf(fpt,"%3.2f]\n",fw[np-1]);		
	
	fprintf(fpt,"figure(20); plot(a,b);\n");
	fprintf(fpt,"figure(20); plot(a,c)\n");
	fclose(fpt);	
	#endif			
}


//check whether a LinkedPoint is valid. 
//return flag, also estimates of new positions and radius and angle if flag = 1.
#define WKSP	10000							//maximum size of workspace sw,fw for lines. 
int checkValidityOfLinkedPoint(int x, int y, int z, int r, int nx, int ny, int nz, float *im3d, int *xo, int *yo, int *ro, float *angle)
{
	int i,j,i1,j1,jj;
	float sw[WKSP],fw[WKSP],smoothfw[WKSP],deriv[WKSP];	//workspace for getting the lines. 
	int np;						//number of points in the line. 
	float ang; 
	int ndiv = 8;				//number of angles to test for the valleyness. 
	int iang;
	int maxs,is;					//maximum of s. 
	float sa, ca;
	int ii, ii0, np0, np1, np2, npm, iam, flag;
	float mmd, sigma, threshold, threshold2, threshold3, mmin, mmax, mmax2, mmax3;
	float fm[ndiv], dm;	//record the valley points. 
	int xm[ndiv],ym[ndiv],rm[ndiv];
	float dd, maxdd;
	int Wx, Wy;
	int itry;
	float maxwithr;
	#if checkValidityOfLinkedPoint_DIAG
	FILE *fpt;
	printf("Writing diagonostic information to temp2.py\n");
	fpt = fopen("temp2.py","w");
	fprintf(fpt,"from pylab import *\n");
	char filename[1000];
	sprintf(filename,"%s","NeuTu.Parameters.dat");
	readParametersFromFile(filename);
	#endif
	
	(*xo) = x; (*yo) = y; (*ro) = r; (*angle) = 0.0;
	
	mmax = FLT_MIN;
	mmax2 = FLT_MIN;
	mmax3 = FLT_MIN;
	maxwithr = FLT_MIN;
	iam = 0;
	for (iang=0; iang < ndiv; iang++) {
		//get sw, fw 
		ang = 3.14156/ndiv*iang;
		sa = sin(ang); ca = cos(ang);
		np = 0;
		maxs = (int) fmax(PARAMS.minRange, 4*r);
		mmin = FLT_MIN; np0 = 0;
		for (is = -maxs; is<maxs; is++) {
			i = round(x + is * ca);
			j = round(y + is * sa);
			if (i<0 || i>=nx || j<0 || j>=ny) continue;
			ii = getIndex(i,j,nx,ny);
			if (np == 0) {
				ii0 = ii;
			} else if (ii == ii0) {
				continue;
			}
			if (abs(is) < mmin) {
				mmin = abs(is);
				np0 = np;
			}
			sw[np] = is;										//distance along the line in this angle
			fw[np] = im3d[getIndex3DZfirst(i,j,z,nx,ny,nz)];	//this is the profile
			np++;
			ii0 = ii;
			if (np >= WKSP) {
				printf("IN checkValidityOfLinkedPoint the WKSP needs to be increased!!!\n");
				break;
			}
		}
		
		Wx = (int)fmax(5,5*PARAMS.sigmaSmoothCurve);
		Wx = Wx/2*2 + 1;
		Wy = Wx;		
		gaussianFilter2D(np,1,Wx,Wy,fw,smoothfw,PARAMS.sigmaSmoothCurve);

		smoothfw[np-1] = smoothfw[np-2]; //get rid of the boundary effects 	
		//get the derivative
		for (is=1; is<np-1; is++) {
			deriv[is] = (smoothfw[is+1] - smoothfw[is-1])/2.0;
		}
		deriv[0] = smoothfw[1] - smoothfw[0];
		deriv[np-1] = smoothfw[np-1] - smoothfw[np-2];		
		gaussianFilter2D(np,1,Wx,Wy,deriv,deriv,PARAMS.sigmaSmoothCurve);			
		deriv[0] = deriv[1];
		deriv[np-1] = deriv[np-2];	//get rid of boundary effect. 
		
		//get threhold
		//get the fluctuation
		mmd = getSparseThreshold(np,1,smoothfw,PARAMS.sparseThreshold);		
		sigma = 0.0;
		for (is=0; is<np; is++) {
			if (smoothfw[is] >= mmd) {	// only include points in the flat part to compute the deviation. Otherwise the peak part produce too much deviation. 
				sigma += (fw[is] - smoothfw[is]) * (fw[is] - smoothfw[is]);
			}
		}
		sigma = sqrt(sigma/np);
		//get threshold
		threshold = mmd - PARAMS.factSigmaThreshold * sigma;	//threshold for judging the validity of the peak. 
		//threshold2 = mmd - sigma;	//threshold for judging if the peak has two flanking flats. 
		
		//get the valley peak. 
		np1 = -1; np2 = -1;
		mmin = smoothfw[np0];
		npm = np0;
		for (is = np0-1; is >=0; is--) { //get the range of the peak
			if (smoothfw[is] < mmin) {
				mmin = smoothfw[is];
				npm = is;
			}
			if (smoothfw[is] >= threshold) {//limit the min to the local minium. 	
				np1 = is;
				break;
			}
		}	
		for (is = np0+1; is < np; is++) { //get the range of the peak
			if (smoothfw[is] < mmin) {
				mmin = smoothfw[is];
				npm = is;
			}
			if (smoothfw[is] >= threshold) {//keep the min to the local minimum.
				np2 = is;
				break;
			}
		}
		xm[iang] = round(x + sw[npm] * ca);
		ym[iang] = round(y + sw[npm] * sa);			
		if (np1 == -1 || np2 == -1 || smoothfw[npm] >= threshold) {	//there is no peak.
			rm[iang] = 0;
			fm[iang] = -1;
		} else {
			//check derivative to get the width. 
			np1 = -1;
			maxdd = 0.0;
			for (is = npm-1; is >=0; is--) {
				dd = -deriv[is];
				if (dd > maxdd) {
					maxdd = dd;
					np1 = is;	
				} else if (maxdd > sigma * PARAMS.factSmallDeriv && dd <= maxdd - sigma * PARAMS.factSmallDeriv) { //PARAMETER!
					break;
				}
			}	
			if (maxdd > mmax) mmax = maxdd;	
			np2 = -1;
			maxdd = 0.0;
			for (is = npm; is < np; is++) {
				dd = deriv[is];
				if (dd > maxdd) {
					maxdd = dd;
					np2 = is;
				} else if (maxdd > sigma * PARAMS.factSmallDeriv && dd <= maxdd - sigma * PARAMS.factSmallDeriv) {	//PARAMETER!
					break;
				}
			}
			if (maxdd > mmax) mmax = maxdd;	
			if (np1 == -1 || np2 == -1) {	//there is no peak.
				rm[iang] = 0;
				fm[iang] = -1;
			} else {
				rm[iang] = abs(sw[np2] - sw[np1])/2.0; 
				fm[iang] = threshold  - smoothfw[npm];
				if ( 1.0/rm[iang] > mmax3) {		//record the angle with smallest radius. 
					mmax3 = 1.0/rm[iang];
					iam = iang;
					maxwithr = fmax(smoothfw[np2],smoothfw[np1]);	//record the maximum value of the intensity. 
				}
			}
		}
		if (fm[iang] > mmax2) {
			mmax2 = fm[iang];
			//iam = iang;
		}
			
		#if checkValidityOfLinkedPoint_DIAG
		//diagonose. copy paste into ipython 
		fprintf(fpt,"angle=%f\n",ang);
		fprintf(fpt,"a=[");
		for (is=0;is<np-1; is++) {
			fprintf(fpt,"%f, ",sw[is]);
		}
		fprintf(fpt,"%f]\n",sw[np-1]);
		fprintf(fpt,"b=[");
		for (is=0;is<np-1; is++) {
			fprintf(fpt,"%f, ",smoothfw[is]);
		}
		fprintf(fpt,"%f]\n",smoothfw[np-1]);		
		fprintf(fpt,"c=[");
		for (is=0;is<np-1; is++) {
			fprintf(fpt,"%f, ",fw[is]);
		}
		fprintf(fpt,"%f]\n",fw[np-1]);		
		fprintf(fpt,"d=[");
		for (is=0;is<np-1; is++) {
			fprintf(fpt,"%f, ",deriv[is]);
		}
		fprintf(fpt,"%f]\n",deriv[np-1]);		
		
		fprintf(fpt,"figure(%d); plot(a,b); title('iang=%f')\n",10+iang+1,ang);
		fprintf(fpt,"figure(%d); plot(a,c)\n",10+iang+1);
		fprintf(fpt,"plot([a[0],a[-1]],[%f,%f])\n",mmd,mmd);
		fprintf(fpt,"plot([a[0],a[-1]],[%f,%f])\n",threshold,threshold);	
		fprintf(fpt,"figure(%d); plot(a,d); title('iang=%f')\n",20+iang+1,ang);
		fprintf(fpt,"plot([a[0],a[-1]],[0,0])\n");
		if (np1 != -1 && np2 != -1 && rm[iang] !=0) 
		fprintf(fpt,"plot([a[%d],a[%d]],[-%f,%f]); plot([a[%d],a[%d]],[-%f,%f])\n",np1,np1,mmax,mmax,np2,np2,mmax,mmax); 
		fprintf(fpt,"figure(5); plot(a,b); title('Smoothed profiles')\n");		
		fprintf(fpt,"figure(6); plot(a,d); title('Derivatives')\n");		
		#endif	
	}	
	flag = 1;
	if (mmax2 < 0) {
		#if checkValidityOfLinkedPoint_DIAG
		printf("Rejecte. No peak found.\n");	
		#endif
		flag = 0;
	} else if (rm[iam] <= PARAMS.minRadius) {
		#if checkValidityOfLinkedPoint_DIAG
		printf("Rejecte. Radius too small.\n");	
		#endif
		flag = 0;
	} else if (rm[iam] >= PARAMS.maxRadius) {
		#if checkValidityOfLinkedPoint_DIAG
		printf("Rejecte. Radius too large.\n");	
		#endif
		flag = 0;
	} else {
		dd = sqrt((xm[iam] - x)*(xm[iam] - x) + (ym[iam] - y)*(ym[iam] - y)); //shift in the estimate.  
		if (dd > rm[iam] * PARAMS.factShift) {
			#if checkValidityOfLinkedPoint_DIAG
			printf("Rejecte. Shift in the center larger than twice of the radius.\n");	
			#endif		
			flag = 0;		//PARAMETER! Tolerance on the shift of the center.
		}
	}
	//check that the profile in the orthogonal direction are dark within the radius. 
	//get the orthogonal profile, recentered to the darkest point. 
	//get sw, fw 
	for (itry =0; itry < 2; itry++) { 
		if (itry == 0) {
			iang = iam;
		} else {
			iang = (int) fmod(iam + ndiv/2,ndiv);	//orthogonal direction.
		} 
		ang = 3.14156/ndiv*iang;
		sa = sin(ang); ca = cos(ang);
		np = 0;
		maxs = (int) fmax(PARAMS.minRange, 4*rm[iam]);
		mmin = FLT_MAX; np0 = 0;
		for (is = -maxs; is<maxs; is++) {
			i = round(xm[iam] + is * ca);
			j = round(ym[iam] + is * sa);
			if (i<0 || i>=nx || j<0 || j>=ny) continue;
			ii = getIndex(i,j,nx,ny);
			if (np == 0) {
				ii0 = ii;
			} else if (ii == ii0) {
				continue;
			}
			if (abs(is) < mmin) {
				mmin = abs(is);
				np0 = np;
			}
			sw[np] = is;										//distance along the line in this angle
			fw[np] = im3d[getIndex3DZfirst(i,j,z,nx,ny,nz)];	//this is the profile
			np++;
			ii0 = ii;
			if (np >= WKSP) {
				printf("IN checkValidityOfLinkedPoint the WKSP needs to be increased!!!\n");
				break;
			}
		}
		gaussianFilter2D(np,1,Wx,Wy,fw,smoothfw,PARAMS.sigmaSmoothCurve);
		smoothfw[np-1] = smoothfw[np-2]; //get rid of the boundary effects 	
		if (itry == 1) {
			//get the maximum value within radius.
			mmax = FLT_MIN;
			for (is = -rm[iam]*0.5; is <= rm[iam]*0.5; is++) {	//PARAMETERS !
				ii = np0 + is; 
				if (ii <0 || ii >= np) continue;
				if (smoothfw[ii] > mmax) mmax = smoothfw[ii];
			}
		}
		#if checkValidityOfLinkedPoint_DIAG
		fprintf(fpt,"aa=[");
		for (is=0; is < np-1; is++) fprintf(fpt,"%f ,",sw[is]);
		fprintf(fpt,"]\n");
		fprintf(fpt,"bb=[");
		for (is=0; is < np-1; is++) fprintf(fpt,"%f ,",smoothfw[is]);
		fprintf(fpt,"]\n");	
		fprintf(fpt,"figure(100); plot(aa,bb)\n");
		fprintf(fpt,"y0,y1=ylim()\n");
		fprintf(fpt,"plot([%f, %f],[y0,y1])\n",-rm[iam]*0.5,-rm[iam]*0.5);
		fprintf(fpt,"plot([%f, %f],[y0,y1])\n",rm[iam]*0.5,rm[iam]*0.5);
		#endif 
	}

	// now judge if the values are dark enough
	threshold3 = maxwithr + sigma;	//PARAMETERS!
	if (mmax > threshold3) {
		//in the orthogonal direction the image does not stay dark with radius, reject. 
		flag = 0;
		#if checkValidityOfLinkedPoint_DIAG
		printf("In orthogonal direction the images does not stay dark enough within the radius. Reject.\n");
		#endif 
	}
	#if checkValidityOfLinkedPoint_DIAG
	fprintf(fpt,"plot([%f,%f],[%f,%f])\n",sw[0],sw[np-1],threshold3,threshold3);
	fprintf(fpt,"title('Orthorgonal profile')\n");
	#endif  
	if (flag == 1) {
		(*xo) = xm[iam];
		(*yo) = ym[iam];
		(*ro) = rm[iam] * PARAMS.factAdjustRadius;			//PARAMETER! ADJUSTMENT of the radius
		(*angle) = 3.14156/ndiv*iam;
	}	
	#if checkValidityOfLinkedPoint_DIAG
	printf("x=%d y=%d r=%d flag=%d fm=%f iam=%d\n",xm[iam],ym[iam],rm[iam],flag,fm[iam],iam);
	fclose(fpt);
	#endif
	return flag;
} 

//mark pixels between P1 and P2 with the IDs of P1 and P2. 
//returns the index of the SWC point closest tp P1 if the pixels between P1 and P2 contains other SWC points.
//useful for connecting ends, espeically in growFromEndPoints.
int markPixelsOccupiedTwoPoints(int nx, int ny, int nz, int *flagOcc, int ii, int x1, int y1, int z1, int r1, int jj, int x2, int y2, int z2, int r2)
{
	int i,j,k;
	float frex, zex, zexc;
	int xss, xee, yss, yee, zss, zee;
	int im;
	float ss, cs;
	long int kk;
	int r11, r22;
	float d, d1, d2;
	int icon;
	float mmin = FLT_MAX;
	float rmm, r1m, r2m;
	
	//mark the points occupied by the SWC points.
	frex = PARAMS.factMarkOccXY; 		
	zex = PARAMS.zOcc;			
		
	r11 = r1 * frex;
	r22 = r2 * frex;
	zexc = fmax(zex,frex*fmax(r1,r2)/PARAMS.zfact);	 
	
	//bounding box of the two points. 
	xss = (int)fmax(0,fmin(x1 - r11, x2 - r22));
	yss = (int)fmax(0,fmin(y1 - r11, y2 - r22));
	zss = (int)fmax(0,fmin(z1 - zexc,   z2 - zexc)); 
	xee = (int)fmin(nx,fmax(x1 + r11 + 1, x2 + r22 +1));
	yee = (int)fmin(ny,fmax(y1 + r11 + 1, y2 + r22 +1));
	zee = (int)fmin(nz,fmax(z1 + zexc +1, z2 + zexc+1)); 

	//go through pixels in the bounding box. 
	icon = -1;
	d  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	for (i=xss; i<xee; i++) {
		for (j=yss; j<yee; j++) {
			im = -1;
			d1 = sqrt((i-x1)*(i-x1) + (j-y1)*(j-y1));
			d2 = sqrt((i-x2)*(i-x2) + (j-y2)*(j-y2));
			if (d1 < r11 && d2 < r22) {	//inside both P1. P2
				if (d1 < d2) {
					im = ii;
				} else {
					im = jj;
				}
			} else if (d1 < r11) {	//inside P1
				im = ii;
			} else if (d2 < r22) {	//inside P2
				im = jj;
			} else if (d > 0) {	//in between, also check for min distance. 
				cs = (d1*d1 + d*d - d2*d2)/(2*d1*d);
				ss = sqrt(1-cs*cs);
				rmm = (d1*cs/d)*(r22-r11)+r11; 
				if ( rmm > d1*ss) { 
					if (d1 < d2) {
						im = ii;	 
					} else {
						im = jj;
					}
					//see if the point is in between P1, P2 and interior. 
					r1m =  (r11/d)*(r22-r11)+r11;
					r2m =  ((d - r11 - r22)/d)*(r22-r11)+r11;
					if (r2m > 0 && r1m > 0) {
						if (((r1m < rmm) && (rmm < r2m)) || (r2m < rmm && (rmm < r1m))) {//interior point
							im = -im - 1; //mark it. 
						}
					}
				}
			}
			if (im != -1) {
				for (k=zss; k<zee; k++) {
					kk = getIndex3DZfirst(i,j,k,nx,ny,nz);
					if (im >= 0) {
						flagOcc[kk] = im;
					} else if (flagOcc[kk] == -1) {
						flagOcc[kk] = -(im + 1);
					} else if (flagOcc[kk] != ii) {	//pixel already marked, do not assgn new label, consider distance to P1. Do not count P1!
						if (d1 < d2 && d1 < mmin) {
							icon = flagOcc[kk];
							mmin = d1;
						}
					}						
				}
			}
		}
	}	
	if (icon == -1) icon = jj;
	return icon;
}

int checkInBound(LinkedPoint *PP, int nx, int ny, int nz)
{
	if (PP->x < 0 || PP->x >=nx || 	//too close to the boundary, do not include as in bound. 
		PP->y < 0 || PP->y >=ny ||
		PP->z < 0 || PP->z >= nz) {
		return 0;
	} else {
		return 1;
	}	 
}


//erase pixels between P1 and P2. 
void erasePixelsOccupiedTwoPoints(int nx, int ny, int nz, float *im3d, int ii, int x1, int y1, int z1, int r1, int jj, int x2, int y2, int z2, int r2)
{
	int i,j,k;
	float frex, zex, zexc;
	int xss, xee, yss, yee, zss, zee;
	int im;
	float ss, cs;
	long int kk;
	int r11, r22;
	float d, d1, d2;
	float mmin = FLT_MAX;
	float rmm, r1m, r2m;
	
	//mark the points occupied by the SWC points.
	frex = PARAMS.factMarkOccXY; 		
	zex = PARAMS.zOcc;			
		
	r11 = r1 * frex;
	r22 = r2 * frex;
	zexc = fmax(zex,1.5*fmax(r1,r2)/PARAMS.zfact);		//Z exclusion is largere for larger radius. PARAMETERS!
	
	//bounding box of the two points. 
	xss = (int)fmax(0,fmin(x1 - r11, x2 - r22));
	yss = (int)fmax(0,fmin(y1 - r11, y2 - r22));
	zss = (int)fmax(0,fmin(z1 - zexc,   z2 - zexc)); 
	xee = (int)fmin(nx,fmax(x1 + r11 + 1, x2 + r22 +1));
	yee = (int)fmin(ny,fmax(y1 + r11 + 1, y2 + r22 +1));
	zee = (int)fmin(nz,fmax(z1 + zexc +1, z2 + zexc+1)); 
	
	//go through pixels in the bounding box. 
	d  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	for (i=xss; i<xee; i++) {
		for (j=yss; j<yee; j++) {
			im = -1;
			d1 = sqrt((i-x1)*(i-x1) + (j-y1)*(j-y1));
			d2 = sqrt((i-x2)*(i-x2) + (j-y2)*(j-y2));
			if (d1 < r11 || d2 < r22) {	//inside both P1. P2
				im = 1;
			} else if (d > 0) {	//in between, also check for min distance. 
				cs = (d1*d1 + d*d - d2*d2)/(2*d1*d);
				ss = sqrt(1-cs*cs);
				rmm = (d1*cs/d)*(r22-r11)+r11; 
				if ( rmm > d1*ss) { 
					if (d1 < d2) {
						im = 1;	 
					} else {
						im = 1;
					}
					//see if the point is in between P1, P2 and interior. 
					r1m =  (r11/d)*(r22-r11)+r11;
					r2m =  ((d - r11 - r22)/d)*(r22-r11)+r11;
					if (r2m > 0 && r1m > 0) {
						if (((r1m < rmm) && (rmm < r2m)) || (r2m < rmm && (rmm < r1m))) {//interior point
							im = 1; //mark it. 
						}
					}
				}
			}
			if (im == 1) {
				for (k=zss; k<zee; k++) {
					kk = getIndex3DZfirst(i,j,k,nx,ny,nz);
					im3d[kk] = 1.0;
				}
			}
		}
	}	
}

//erase image around the SWC points.  
void erasePixelsOccupied(int nx,int ny,int nz,float *im3d, int npoints,LinkedPoint **points)
{
	int i,j,k,ii,jj;
	long int kk, ntot;
	int *iid, *IDInds;
	short *flagDone;
	int mmax;
	int iend, id;
	LinkedPoint *P1, *P2;
	LinkedList *p;
	int ne;
	int x1,y1,z1,r1,x2,y2,z2,r2;

	ntot = nx*ny*nz;
		
	//get the junction points. 
	iid = (int *) malloc((npoints * sizeof(int)));
	flagDone = (short *) malloc(npoints * sizeof(short));

	//setup a table for getting the SWC point from id. 
	mmax = 0;
	for (i=0;i<npoints;i++) {
		if (mmax < points[i]->ID) mmax = points[i]->ID;
		flagDone[i] = 0;
	}
	IDInds = (int *) malloc((mmax+1) * sizeof(int));
	for (i=0; i<npoints; i++) {
		IDInds[points[i]->ID] = i;
	}
	
	//get the end points and junction points.  
	ne = 0;
	for (i=0; i<npoints; i++) {
		if (NumConnLinkedPoint(points[i]) == 2) continue; //interior point. 
		for (j=0; j<NumConnLinkedPoint(points[i]); j++) { //dupilcate the connection points. 
			iid[ne] = i;
			ne++;
		}
	}
	
	//go through the end points and connected points and mark occupied pixels. 
	for (iend = 0; iend < ne; iend++) {
		ii = iid[iend];
		flagDone[ii] = 1;
		P1 = points[ii];
		x1 = P1->x; y1 = P1->y; z1 = P1->z; r1 = P1->r;
		p = P1->conn->next;
		if (checkInBound(P1,nx,ny,nz) == 0) continue;
		while (p != NULL) {
			jj = IDInds[p->val];
			if (flagDone[jj] == 0 && checkInBound(points[jj],nx,ny,nz)==1) break;
			p = p->next;
		}
		if (p == NULL) {
			P2 = NULL;
			x2 = x1; y2 = y1; z2 = z1; r2 = r1;
		} else {
			P2 = points[jj];
			x2 = P2->x; y2 = P2->y; z2 = P2->z; r2 = P2->r;	
		}
		while (1) {
			erasePixelsOccupiedTwoPoints(nx, ny, nz, im3d, ii, x1, y1, z1, r1, jj, x2, y2, z2, r2);
			if (P2 == NULL) break;
			//get the next connected point. 
			P1 = P2; 
			if (NumConnLinkedPoint(P2) != 2) break;	//reached an end or branching point. 
			ii = jj;
			p = P2->conn->next;
			while (p != NULL) {
				jj = IDInds[p->val];
				if (flagDone[jj] == 0 && checkInBound(points[jj],nx,ny,nz)==1) break;
				p = p->next;
			}		
			if (p == NULL) break;
			P2 = points[jj];
			flagDone[jj] = 1;	
			x1 = x2; y1 = y2; z1 = z2; r1 = r2;
			x2 = P2->x; y2 = P2->y; z2 = P2->z; r2 = P2->r;	
		}
	}
	free(iid);
	free(flagDone);
	free(IDInds);
}

//mark pixels in flagOcc with the IDs of neareast SWC points. 
void markPixelsOccupied(int nx,int ny,int nz,int *flagOcc, int npoints,LinkedPoint **points)
{
	int i,j,k,ii,jj;
	int *IDInds;
	int mmax;
	LinkedPoint *P1, *P2;
	LinkedList *p;
	int x1,y1,z1,r1,x2,y2,z2,r2;
	
	//setup a table for getting the SWC point from id. 
	mmax = 0;
	for (i=0;i<npoints;i++) {
		if (mmax < points[i]->ID) mmax = points[i]->ID;
	}
	IDInds = (int *) malloc((mmax+1) * sizeof(int));
	for (i=0; i<npoints; i++) {
		IDInds[points[i]->ID] = i;
	}		
	//go through all points and their connections 
	for (ii=0; ii<npoints; ii++) {
		P1 = points[ii];
		if (checkInBound(P1,nx,ny,nz) == 0) continue;
		x1 = P1->x; y1 = P1->y; z1 = P1->z; r1 = P1->r;
		p = P1->conn->next;
		while (p != NULL) {
			jj = IDInds[p->val];
			P2 = points[jj];
			if (checkInBound(P2,nx,ny,nz)==1) {
				x2 = P2->x; y2 = P2->y; z2 = P2->z; r2 = P2->r;	
				markPixelsOccupiedTwoPoints(nx, ny, nz, flagOcc, ii, x1, y1, z1, r1, jj, x2, y2, z2, r2);
			}
			p = p->next;
		}
	}
	free(IDInds);
}

//get Linked points from the skeleton.
//bw contains the skeleton in the binary image, distMap contains the distance map.  
#define getLinkedPointsFromSkeletonAnd3D_DIAG 0
#define getLinkedPointsFromSkeletonAnd3D_PLOTZMAP 0
#if getLinkedPointsFromSkeletonAnd3D_PLOTZMAP
int c_zmap = 0;	//counter for outputing zmap. Used for illustrating the method of getting z. 
#endif
void getLinkedPointsFromSkeletonAnd3D(int nx,int ny,int nz,float *bw, float *distMap, float *im3d,int z1,int z2,int z1ext,int z2ext,int *npoints,LinkedPoint **points, float zfact)
{	
	float *zmap;
	int i,j,k,i1,j1,i2,j2,j3,ii,jj,kk,k1,k2,dk,nxy,np,ne,nnp,jjp[8],fl,iDone,k3;
	float *flags, mr, *z;
	int *x, *y, *endPoints, maxdp;
	int *iidList, niid, flag, nzz;
	float r1,r2,zJumpFact,*r,dd,xs,ys,xe,ye,dotp,dd2,dd3;
	int nxym;
	int nw,*iidw,ni,iip[8],ii1,jj1,kstart,kend;
	float angle;
	int kk1;
	int xo, yo, ro, zo, f1, f2, f3, dr;
	float ao;
	int *xpp,*ypp,*zpp,*rpp;
	int npp;
	int *flagOcc;
	int ntot;
	int x11, y11, z11, r11, id1, x22, y22, z22, r22, id2;
	int nid;
	int WySmooth;
			
	zJumpFact = PARAMS.zJumpFact;
	nxy = nx * ny;
	nxym = (int)(nxy * 0.2);
	flags = (float *) malloc(nxy * sizeof(float));	
	x = (int *) malloc(nxym * sizeof(int));
	y = (int *) malloc(nxym * sizeof(int));
	z = (float *) malloc(nxym * sizeof(float));
	r = (float *) malloc(nxym * sizeof(float));	
	endPoints = (int *) malloc(nxym * sizeof(int));
	zmap = (float *) malloc(nxym * nz * sizeof(float));
	iidList = (int *) malloc(nxym * sizeof(int));
	iidw = (int *) malloc(nxym * sizeof(int));
	xpp = (int *) malloc(nxym * sizeof(int));
	ypp = (int *) malloc(nxym * sizeof(int));
	zpp = (int *) malloc(nxym * sizeof(int));
	rpp = (int *) malloc(nxym * sizeof(int));	
	ntot = nx * ny * nz;
	flagOcc = (int *) malloc(ntot * sizeof(int));
	for (i=0; i<ntot; i++) flagOcc[i] = -1;	//inilialize. 
	nid = (*npoints);
		
	//mark the points occupied by the SWC points.
	markPixelsOccupied(nx, ny, nz, flagOcc, (*npoints), points);	
	
	//get rid of boundary points to help judging the end points.
	for (i=0; i<nx; i++) {
		bw[getIndex(i,0,nx,ny)] = 0;
		bw[getIndex(i,1,nx,ny)] = 0;
		bw[getIndex(i,ny-1,nx,ny)] = 0;
		bw[getIndex(i,ny-2,nx,ny)] = 0;		
	}	
	for (j=0; j<ny; j++) {
		bw[getIndex(0,j,nx,ny)] = 0;
		bw[getIndex(1,j,nx,ny)] = 0;
		bw[getIndex(nx-1,j,nx,ny)] = 0;
		bw[getIndex(nx-2,j,nx,ny)] = 0;
	}	
	for (ii=0; ii<nxy; ii++) flags[ii] = bw[ii];
	nzz = z2ext - z1ext;

	//get the ids of the white pixels.
	nw = 0; 
	for (ii=0; ii<nxy; ii++) if (flags[ii] > 0) iidw[nw++] = ii;

	//get rid of kinks.
	for (kk=0; kk<nw; kk++) {
		ii = iidw[kk];
		i = ii/ny;
		j = ii - i*ny;
		nnp = 0;		//get white pixels around (i,j)
		for (i1 = i-1; i1 < i+2; i1++) {
		for (j1 = j-1; j1 < j+2; j1++) {
			jj = getIndex(i1,j1,nx,ny);
			if (jj == ii || jj < 0 || jj >nxy-1) continue;
			if (flags[jj] > 0.0) {
				jjp[nnp++] = jj;
			}
		}}
		if (nnp >= 2) {
			//see if these white pixels are all connected.
			iip[0] = jjp[0];
			ni = 1;
			while (1) {
				fl = 0;
				for (ii1=0; ii1<nnp; ii1++) {
					flag = 0;
					for (jj1=0; jj1<ni; jj1++) {
						if (jjp[ii1] == iip[jj1]) {
							flag = 1;
							break;
						} 
					}
					if (flag == 1) continue; //point ii1 is already in the set
					flag = 0;
					for (jj1=0; jj1<ni; jj1++) {
						i1 = jjp[ii1]/ny; j1 = jjp[ii1] - i1 * ny;
						i2 = iip[jj1]/ny; j2 = iip[jj1] - i2 * ny;
						if (sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2)) <= sqrt(2.0)) {
							iip[ni++] = jjp[ii1];	//point ii1 is added to the set. 
							flag = 1;
							break;
						}
					}	
					if (flag == 1) {
						fl = 1;
						break;
					}
				}
				if (fl == 0) break;			
			}
			if (ni == nnp) {	//connected, the center point is redundent.
				flags[ii]= 0;			
			}	
		}
	}

    //break crossing points. 
    for (kk=0; kk<nw; kk++) {
		ii = iidw[kk];
		i = ii/ny;
		j = ii - i*ny;
		nnp = 0;		
		for (i1 = i-1; i1 < i+2; i1++) {
		for (j1 = j-1; j1 < j+2; j1++) {
			jj = getIndex(i1,j1,nx,ny);
			if (jj == ii || jj < 0 || jj >nxy-1) continue;
			if (flags[jj] > 0.0) {
				jjp[nnp++] = jj;
			}
		}}
		if (nnp > 2) { //this is the crossing point, create end points by erasing. 
			//see if these white pixels are all connected.
			iip[0] = jjp[0];
			ni = 1;
			while (1) {
				fl = 0;
				for (ii1=0; ii1<nnp; ii1++) {
					flag = 0;
					for (jj1=0; jj1<ni; jj1++) {
						if (jjp[ii1] == iip[jj1]) {
							flag = 1;
							break;
						} 
					}
					if (flag == 1) continue;
					flag = 0;
					for (jj1=0; jj1<ni; jj1++) {
						i1 = jjp[ii1]/ny; j1 = jjp[ii1] - i1 * ny;
						i2 = iip[jj1]/ny; j2 = iip[jj1] - i2 * ny;
						if (sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2)) <= sqrt(2.0)) {
							iip[ni++] = jjp[ii1];
							flag = 1;
							break;
						}
					}	
					if (flag == 1) {
						fl = 1;
						break;
					}
				}
				if (fl == 0) break;			
			}
			if (ni != nnp) {	//not connected, the center point is a crossing point.
				for (i1 = i-1; i1 < i+2; i1++) {
				for (j1 = j-1; j1 < j+2; j1++) {
					jj = getIndex(i1,j1,nx,ny);
					flags[jj]=0;
				}}					
			}	
		}
	}		
	nw = 0; 
	for (ii=0; ii<nxy; ii++) if (flags[ii] > 0) iidw[nw++] = ii;

	//createGrayTiffStackFromFloatArray(nx,ny,1,flags,"temp.tif");
	//create linked points. 
	//get end points. 
	ne = 0;
	//get end points.
    for (kk=0; kk<nw; kk++) {
		ii = iidw[kk];
		i = ii/ny;
		j = ii - i*ny;
		nnp = 0;		
		for (i1 = i-1; i1 < i+2; i1++) {
		for (j1 = j-1; j1 < j+2; j1++) {
			jj = getIndex(i1,j1,nx,ny);
			if (jj == ii || jj < 0 || jj >nxy-1) continue;
			if (flags[jj] > 0.0) {
				jjp[nnp++] = jj;
			}
		}}
		if (nnp <= 1) {
			endPoints[ne++] = ii;
		}
	}			

	// construct path from the end points to the next
	for (j3=0;j3<ne; j3++) {
		ii = endPoints[j3];
		if (flags[ii] == 0.0) {//this point is already erased. 
			continue;
		}
		np = 0;
		i = ii/ny; j = ii - i*ny;
		x[np] = i; y[np] = j; r[np] = distMap[ii]; np++;
		mr = 0;
		while (1) { //follow the skeleton.
			flags[ii] = 0; 			
			mr += distMap[ii];
			//select next point
			nnp = 0;		
			for (i1 = i-1; i1 < i+2; i1++) {
			for (j1 = j-1; j1 < j+2; j1++) {
				jj = getIndex(i1,j1,nx,ny);
				if (jj == ii || jj < 0 || jj >nxy-1) continue;
				if (flags[jj] > 0.0) {
					jjp[nnp++] = jj;
				}
			}}
			fl = 0;
			if (nnp == 1) {
				kk = jjp[0];	//this is the next point. 
				fl = 1;
			} else if (nnp == 2) {
				i1 = jjp[0]/ny; j1 = jjp[0] - i1 * ny;
				i2 = jjp[1]/ny; j2 = jjp[1] - i2 * ny;
				if (sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2)) == 1.0) {
					fl = 1;
					if (i1 == i || j1 == j) {
						kk = jjp[0];
						flags[jjp[1]] = 0;
					} else if (i2 == i || j2 == j) {
						kk = jjp[1];
						flags[jjp[0]] = 0;
					}
				}
			}
			if (fl == 1) {	//get to the next point. 
				ii = kk;
				i = ii/ny; j = ii - i*ny;
				x[np] = i; y[np] = j; r[np] = distMap[ii]; 
				np++;
			} else {	//end point or cross point, stop. 
				break;
			}
		}			
		if (np > PARAMS.smallLen) { // get LinkedPoints from this branch. 
			//printf("creating branch np=%d ",np);
			//compute z.
			kk = 0;
			for (i=z1ext; i<z2ext; i++) {
				for (j=0; j<np; j++) {
					jj = getIndex3DZfirst(x[j],y[j],i,nx,ny,nz);
					zmap[kk++] = im3d[jj];
				}
			} 
			shortestPathImageLeftToRight(nzz,np,zmap,z,zfact);
			for (i=0; i<np; i++) z[i] += z1ext;
			
			//smooth x, y, z
			WySmooth = (int) fmax(3, 3*PARAMS.sigmaSmooth);
			WySmooth = WySmooth/2*2 + 1;
			
			if (np > WySmooth + 1 ) {
				gaussianFilter2D(1,np,1,WySmooth,(float *)x,(float *)x,PARAMS.sigmaSmooth);
				gaussianFilter2D(1,np,1,WySmooth,(float *)y,(float *)y,PARAMS.sigmaSmooth);
				gaussianFilter2D(1,np,1,WySmooth,(float *)z,(float *)z,PARAMS.sigmaSmooth);
				gaussianFilter2D(1,np,1,WySmooth,r,r,PARAMS.sigmaSmooth);
			}
			
			//save Z map to images. Only useful for writing the paper!
			#if getLinkedPointsFromSkeletonAnd3D_PLOTZMAP
			char filename[1000];
			float *xyimg,*zimg;
			float mmax = FLT_MIN, mmin = FLT_MAX;
			xyimg = (float *) malloc(nx*ny*sizeof(float));
			zimg = (float *) malloc((z2ext-z1ext)*np*sizeof(float));
			for (i=0; i<nxy; i++) xyimg[i] = 0;
			for (j=0; j<np; j++) {
				jj = getIndex(x[j],y[j],nx,ny);
				xyimg[jj] = 1;
			}	
			sprintf(filename,"zmap%d.xy.tif",c_zmap);
			printf("Saving zmap to %s. \nTurn this option off by setting getLinkedPointsFromSkeletonAnd3D_PLOTZMAP to 0 in libNeuTu.c\n",filename);
			createGrayTiffStackFromFloatArray(nx,ny,1,xyimg,filename);		
			sprintf(filename,"zmap%d.zmap.tif",c_zmap);
			createGrayTiffStackFromFloatArray(z2ext-z1ext,np,1,zmap,filename);			
			for (i=z1ext; i<z2ext; i++) {
				for (j=0; j<np; j++) {
					jj = getIndex(i,j,z2ext-z1ext,np);
					zimg[jj] = 0;
				}
			}	 
			for  (j=0; j<np; j++) {
				jj = getIndex(z[j]-z1ext,j,z2ext-z1ext,np);
				zimg[jj] = 1.0;
			}
			sprintf(filename,"zmap%d.z.tif",c_zmap);
			createGrayTiffStackFromFloatArray(z2ext-z1ext,np,1,zimg,filename);			
			++c_zmap;
			free(xyimg);
			free(zimg);
			#endif
			
			//get the list of points to create LinkedPoints. 
			niid = 0;
			//get the start point. start with a point r away from the break point. 
			xs = x[0]; ys = y[0]; 
			for (k3 = 0; k3 <np; k3++) {			
				dd = sqrt((x[k3] - xs)*(x[k3] - xs) + (y[k3] - ys)*(y[k3] - ys));
				if (dd >= r[k3]) break;
			}
			iidList[niid++] = k3; // the first point. 
			xs = x[k3]; ys = y[k3]; r1 = r[k3];
			xe = x[np-1]; ye = y[np-1];
			for (k3 = iidList[0]; k3 <np; k3++) {
				r2 = r[k3];
				dd = sqrt((x[k3] - xs)*(x[k3] - xs) + (y[k3] - ys)*(y[k3] - ys));
				dd3 = sqrt((x[k3] - xe)*(x[k3] - xe) + (y[k3] - ye)*(y[k3] - ye));
				//if (dd3 <= r[k3]) break; //close to the other end. 
				dd2 = (r1+r2)/(1.0+PARAMS.drFactor*abs(r1-r2));	//make the distance smaller if the change in r is large; set upper limit. 
				if (dd >= dd2 || dd3 <= r[k3]) {
					iidList[niid++] = k3;
					flag = 1;
					r1 = r[k3];
					xs = x[k3]; ys = y[k3];
				}
			}
			iidList[niid++] = np-1; // the last point.
						
			//create list of points to be created. 
			npp = 0;
			for (i=0; i<niid; i++) {
				kk = iidList[i];				
				//adjust z of the new point.
				//adjustZOfLinkedPoint(x[kk], y[kk], z[kk], r[kk], nx, ny, nz, im3d, &zo);
				//z[kk] = zo;
				//check validity. 				
				f1 = checkValidityOfLinkedPoint(x[kk],y[kk],z[kk],r[kk],nx,ny,nz,im3d,&xo,&yo,&ro,&ao);		
				f2 = checkValidityOfLinkedPoint(xo,yo,z[kk],ro,nx,ny,nz,im3d,&xo,&yo,&ro,&ao);		
				if (f1 == 1 && f2 == 1 && 
					flagOcc[getIndex3DZfirst(xo,yo,z[kk],nx,ny,nz)] == -1 && 
					z[kk] >=z1 && z[kk] < z2 && 			//restrict to the current slab. 	
					xo > ro && xo < nx-ro-1 &&			//points near the boundary are unreliable. Get rid of.
					yo > ro && yo < ny-ro- 1) {	
					checkValidityOfLinkedPoint(xo,yo,z[kk],ro,nx,ny,nz,im3d,&xo,&yo,&ro,&ao); //further correct radius, position. 						
					xpp[npp] = xo;
					ypp[npp] = yo;
					zpp[npp] = z[kk];
					rpp[npp] = ro;				
					npp++;					
					//mark pixels as occupied. 
					x22 = xo; y22 = yo; z22 = z[kk]; r22 = ro; id2 = nid;
					if (npp == 1) {
						x11 = x22; y11 = y22; z11 = z22; r11 = r22; id1 = nid;
					} 
					markPixelsOccupiedTwoPoints(nx, ny, nz, flagOcc, id1, x11, y11, z11, r11, id2, x22, y22, z22, r22);
					nid++;
					x11 = x22; y11 = y22; z11 = z22; r11 = r22; id1 = id2;
				}
			}
						
			if (npp >=3) {
				//smooth
				//gaussianFilter2D(1,npp,1,3,(float *)xpp,(float *)xpp,1);
				//gaussianFilter2D(1,npp,1,3,(float *)ypp,(float *)ypp,1);
				//gaussianFilter2D(1,npp,1,3,(float *)zpp,(float *)zpp,1);
				gaussianFilter2D(1,npp,1,3,(float *)rpp,(float *)rpp,1);	
			}
			nnp = 0;
			for (i=0; i<npp; i++) {
				//create new point. 
				points[(*npoints)] = CreateLinkedPoint(xpp[i],ypp[i],zpp[i],rpp[i],(*npoints),3);	//use new estimates. 
				if (i==0) {
					(*npoints)++;
					nnp++;
					continue;
				}				
				//check the z jump and angle. 
				f1 = 1;
				if (npp >= 1) {
					xe = points[(*npoints)]->x - points[(*npoints)-1]->x;
					ye = points[(*npoints)]->y - points[(*npoints)-1]->y; 
					dd = fmax(1e-5,sqrt(xe*xe + ye*ye));
					if (dd > (points[(*npoints)-1]->r + points[(*npoints)]->r)* PARAMS.distFactor) { //PARAMETER! Do not connect if the distance is too large. 
						f1 = 0;
					}
					xe /= dd; ye /= dd;
				} 
				dotp = 1.0;
				if (nnp >= 2) {
					xs = points[(*npoints)-1]->x - points[(*npoints)-2]->x; 
					ys = points[(*npoints)-1]->y - points[(*npoints)-2]->y;
					dd = fmax(1e-5,sqrt(xs*xs + ys*ys));
					xs /= dd; ys /= dd; 
					dotp = xs * xe + ys * ye;
				}
				if (f1 == 1
				&& dotp > cos(PARAMS.angle/180.0*3.1415926) //if the angle is more than 90 degrees do not connect. This prevents noisy side branches. 
				&& zfact * abs(points[(*npoints)-1]->z - points[(*npoints)]->z) <= zJumpFact * (points[(*npoints)-1]->r + points[(*npoints)]->r)) {					 
					AddConnLinkedPoint(points[(*npoints)-1],(*npoints));
					AddConnLinkedPoint(points[(*npoints)],(*npoints)-1);
				} else {
					nnp = 0;
				} 
				(*npoints)++;
				nnp++;
			}
		}
	}	 	

	free(flags);
	free(x);
	free(y);
	free(z);
	free(r);
	free(endPoints);	
	free(zmap);
	free(iidList);
	free(iidw);
	free(xpp);
	free(ypp);
	free(zpp);
	free(rpp);
	free(flagOcc);
}

//return 1 if LinkedPoint p is connected with a point with ID.  
int checkConnected(LinkedPoint *p, int ID, int npoints, LinkedPoint **points, int *IDInds, int *flags)
{
	int i, flag;
	LinkedPoint *pnext;
	LinkedList *iter;
	if (p == NULL) return 0;
	flags[p->ID] = 1;
	iter = p->conn->next;
	while (iter != NULL) {
		if (iter->val == ID) {
			return 1;
		}
		iter = iter->next;
	}
	iter = p->conn->next;
	while (iter != NULL) {
		i = IDInds[iter->val];
		pnext = points[i];
		if (pnext == NULL) break;
		if (flags[pnext->ID] == 0) { //not visited yet.
			flag = checkConnected(pnext, ID, npoints, points, IDInds, flags);
			if (flag == 1) {
				return 1;
			} 
		}
		iter = iter->next;
	}		
	return 0;
}

//get think dendrite and soma missed by the previoius method. 
#define getThickDendrites_DAG 0	//set to 1 to out put imFlat and mask for diagnosis
void getThickDendrites(int nx,int ny,int nz,float *im3d, int *npoints, LinkedPoint **points)
{	
	int *flagOcc; 
	float *imFlat, *imFlat2, *bw, minI, minI2;
	int mm, nn;
	float sigma;
	int Wx, Wy;
	int xm, ym;
	float thrSparse;
	int i,j,k,ii,ntot;

	sigma  = (int)fmax(1,PARAMS.somaLengthScale/3);
	Wx = (int) fmax(6, 5*sigma);
	Wy = Wx;

	ntot = nx * ny * nz;
	flagOcc = (int *) malloc(ntot * sizeof(int));
	for (i=0; i<ntot; i++) flagOcc[i] = -1;	//inilialize. 
	imFlat2 = (float *) malloc(nx*ny*sizeof(float));
	imFlat = (float *) malloc(nx*ny*sizeof(float));
	bw = (float *) malloc(nx*ny*sizeof(float));
		
	//expand areas of masking around existing points. THis avoids putting points in shadows of images.  
	PARAMS.zOcc *= 3;						//PARAMETERS!
	PARAMS.factMarkOccXY *= 3;					
	//mark the points occupied by the SWC points.
	markPixelsOccupied(nx, ny, nz, flagOcc, (*npoints), points);	
	//restore parameters
	PARAMS.zOcc /= 3;
	PARAMS.factMarkOccXY /= 3;	
	
	//get new projection. 
	//maximumum intensity projection. 
	mm = 0;
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			minI = FLT_MAX;
			minI2 = FLT_MAX;
			for (k=0; k<nz; k++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				if (flagOcc[ii] == -1) { 
					if (minI > im3d[ii]) minI = im3d[ii];
				} else {
					if (minI2 > im3d[ii]) minI2 = im3d[ii];
				}
			}
			ii = getIndex(i,j,nx,ny);
			imFlat[ii] = minI;
			if (minI2 != FLT_MAX) {
				imFlat2[mm++] = minI2;
			}					
		}
	}
	//check if there is thick dendrite done. 
	if (mm >10) {
		thrSparse = getSparseThreshold(mm,1,imFlat2,0.95); //used occupied to get the threshold.
		nn = 0;
		for (i=0; i<nx*ny; i++) {
			if (imFlat[i] < thrSparse) nn++;
		}
	}
	if (nn < 0.5*PARAMS.somaLengthScale*PARAMS.somaLengthScale) {
		printf("No soma. sqrt(nn)=%d PARAMS.somaLengthScale=%d\n",(int)sqrt(nn),PARAMS.somaLengthScale);
		free(flagOcc);
		free(imFlat2);
		free(imFlat);
		free(bw);	
		return;
	}		
	printf("Trying to get soma...\n");

	//maximumum intensity projection. Do not exlude existing points; otherwise shadow points are created. 
	mm = 0;
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			minI = FLT_MAX;
			for (k=0; k<nz; k++) {
				ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
				if (minI > im3d[ii]) minI = im3d[ii];
			}
			ii = getIndex(i,j,nx,ny);
			imFlat[ii] = minI;
		}
	}			
	//smooth imflat. 
	gaussianFilter2D(nx, ny, Wx, Wy, imFlat, imFlat, sigma);
	//get background
	/*
	sigma *= 5; Wx *= 5; Wy *= 5;
	gaussianFilter2D(nx, ny, Wx, Wy, imFlat, imFlat2, sigma);
	//substract background. 
	for (i=0; i<nx*ny; i++) {
		imFlat[i] -= imFlat2[i];
	}
	* */
	
	#if getThickDendrites_DAG
	char filename[1000];
	sprintf(filename,"temp.Soma.2Dproj.tif");
	printf("Saving the 2D projection to %s...\n",filename);
	createGrayTiffStackFromFloatArray(nx,ny,1,imFlat,filename);
	#endif
	
	//get the threshold. 
	thrSparse = getSparseThreshold(nx*ny,1,imFlat,1-PARAMS.somaSparseThr); //used occupied to get the threshold.
	//get the mask
	for (i=0;i<nx*ny;i++) {
		if (imFlat[i] < thrSparse) {
			bw[i] = 1.0;	
		} else {
			bw[i] = 0.0;
		}
	}
	//get rid of bounary
	for (i=0; i<nx; i++) {
		bw[getIndex(i,0,nx,ny)] = 0;
		bw[getIndex(i,ny-1,nx,ny)] = 0;
	}	
	for (j=0; j<ny; j++) {
		bw[getIndex(0,j,nx,ny)] = 0;
		bw[getIndex(nx-1,j,nx,ny)] = 0;
	}	

	#if getThickDendrites_DAG
	sprintf(filename,"temp.Soma.mask.tif");
	printf("Saving the mask to %s...\n",filename);
	createGrayTiffStackFromFloatArray(nx,ny,1,bw,filename);
	#endif

	//get the distance map, save in imFlat.
	sedt(nx, ny, bw, imFlat);
	#if getThickDendrites_DAG
	sprintf(filename,"temp.Soma.distMap.tif");
	printf("Saving the distance map to %s...\n",filename);
	createGrayTiffStackFromFloatArray(nx,ny,1,imFlat,filename);
	#endif

	//printf("Getting skeleton...\n");
	skeletonizeZS(nx,ny,bw,bw);

	#if getThickDendrites_DAG
	sprintf(filename,"temp.Soma.skeleton.tif");
	printf("Saving the skeleton to %s...\n",filename);
	createGrayTiffStackFromFloatArray(nx,ny,1,bw,filename);
	#endif
		
	//prune small branches. 
	//printf("Pruning small branches...\n");
	//pruneSmallBranches(nx,ny,bw,PARAMS.smallLen); 	
		
	//get Linked points from the skeleton.
	//getLinkedPoints
	printf("Getting linked points for soma...\n");
	getLinkedPointsFromSkeletonAnd3D(nx,ny,nz,bw,imFlat,im3d,0,nz,0,nz,npoints,points,PARAMS.zfact);
	
	free(flagOcc);
	free(imFlat2);
	free(imFlat);
	free(bw);	
}

//connect end points in the linked points. 
void connectEndPoints(int npoints, LinkedPoint **points, float zfact)
{
	int i,j,ne,ii,jj,kk;
	float x1,x2,y1,y2,z1,z2,r1,r2,dd,xs,ys;
	int *eps, nE;
	int *flags,mmax,*IDInds;	
	float *epsVec;
	struct valInd element;
	int *visited; 
	float xv,yv,dpot;
	int flag;
	
	mmax = 0;
	for (i=0;i<npoints;i++) {
		if (mmax < points[i]->ID) mmax = points[i]->ID;
	}
	IDInds = (int *) malloc((mmax+1) * sizeof(int));
	for (i=0; i<npoints; i++) IDInds[points[i]->ID] = i;
	eps = (int *) malloc(npoints * sizeof(int));
	epsVec = (float *) malloc(2*npoints * sizeof(int));
	//get the end points. 
	ne = 0;
	for (i=0; i<npoints; i++) {
		if (NumConnLinkedPoint(points[i]) == 1) {
			eps[ne] = i;
			j = IDInds[points[i]->conn->next->val];
			xs = points[i]->x - points[j]->x;
			ys = points[i]->y - points[j]->y;
			dd = fmax(1e-5,sqrt(xs*xs + ys*ys));
			epsVec[2*ne] = xs/dd;
			epsVec[2*ne+1] = ys/dd;			
			ne++;
		}
	}
	if (ne == 0) {
		free(eps);
		free(IDInds); 
		return;
	}
	//compute distances between the end points. 
	flags = (int *) malloc((mmax+1) * sizeof(int));
	visited = (int *) malloc(ne * sizeof(int));
	for (i=0; i<ne; i++) visited[i] = 0;
	//allocate heat
	if (heap == NULL) free(heap); 
	heap = (struct valInd  *) malloc(ne *ne * sizeof(struct valInd));
	Init(); //heap structure for finding minimum.
	for (i=0; i<ne; i++) {
		ii = eps[i];
		x1 = points[ii]->x;
		y1 = points[ii]->y;
		z1 = points[ii]->z;
		r1 = points[ii]->r;
		for (j=i+1; j<ne; j++) {			
			jj = eps[j];
			for (kk=0;kk<mmax+1;kk++) flags[kk] = 0;
			if (checkConnected(points[ii], points[jj]->ID, npoints, points, IDInds, flags) == 1) continue;
			
			x2 = points[jj]->x;
			y2 = points[jj]->y;
			z2 = points[jj]->z;
			r2 = points[jj]->r;
			
			if (zfact * abs(z2 - z1) > PARAMS.zJumpFact * (r1 + r2)) continue; //z difference is too large. 		
			
			xv = x2 - x1;
			yv = y2 - y1;
			dd = fmax(1e-5,sqrt(xv*xv+yv*yv));
			xv /= dd;
			yv /= dd;
		
			if (dd <= 1.5 * (r1 + r2)) { //these two points are touching. Connect. 
				AddConnLinkedPoint(points[ii],points[jj]->ID);
				AddConnLinkedPoint(points[jj],points[ii]->ID);
				visited[i] = 1;
				visited[j] = 1;	//marke these points as connected. 
				continue;			
			} 
			dpot = epsVec[2*i] * xv + epsVec[2*i+1] * yv;
			if (dpot < cos(PARAMS.angle/180.0*3.14159))  continue;  
			dpot = -epsVec[2*j] * xv - epsVec[2*j+1] * yv;
			if (dpot < cos(PARAMS.angle/180.0*3.14159)) continue;   
						
			if (dd < (r1 + r2) * PARAMS.distFactConn) { //candidate for connection, insert to the heap.  
				dd = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (zfact*zfact)*(z1-z2)*(z1-z2));	//compute the 3D distance.
				element.value = dd;
				element.ind = getIndex(i,j,ne,ne);
				InsertOrUpdate(element);
			}
		}
	}
	while (1) {
		if (heapSize == 0) break;
		element = DeleteMin();
		i = element.ind/ne;
		j = element.ind - i * ne;
		if (visited[i] == 1 || visited[j] == 1) continue;
		ii = eps[i];
		jj = eps[j];
		//check if these points are already connected. 
		for (kk=0;kk<mmax+1;kk++) flags[kk] = 0;
		if (checkConnected(points[ii], points[jj]->ID, npoints, points, IDInds, flags) == 1) continue;
		visited[i] = 1; visited[j] = 1;	
		AddConnLinkedPoint(points[ii],points[jj]->ID);
		AddConnLinkedPoint(points[jj],points[ii]->ID);
	}
	free(IDInds);
	free(eps);
	free(epsVec);
	free(flags);	
	free(heap); heap = NULL;
	free(visited); 
}

//delete isolated points
void deleteIsolatedPoints(int *npoints, LinkedPoint **points)
{
	int i, j, maxID, *IDInds, *iidDelete, nd, ii, jj;

	maxID = 0; 	
	for (i=0; i<(*npoints); i++) {
		if (maxID < points[i]->ID) maxID = points[i]->ID;
	}
	IDInds = (int *)malloc((maxID+1) * sizeof(int));	
	for (i=0; i<(*npoints); i++) {		
		IDInds[points[i]->ID] = i;
	}
		
	//now delete small leaf branches. 
	iidDelete = (int *) malloc((*npoints)*sizeof(int));
	nd = 0;
	for (i=0; i<(*npoints); i++) {
		if (NumConnLinkedPoint(points[i]) == 0) {//isolated point.
			iidDelete[nd++] = i; 
		}
	}		
	//delete connections. 
	for (i=0; i<(*npoints); i++) {
		for (j=0; j<nd; j++) {
			DelConnLinkedPoint(points[i],points[iidDelete[j]]->ID);
		}
	}
	//delete the points
	for (j=0; j<nd; j++) {
		if (points[iidDelete[j]] != NULL) {
			DeleteLinkedPoint(points[iidDelete[j]]);
			points[iidDelete[j]] = NULL;
		}
	}
	//make the remaining points consecutive
	ii = 0; jj = 0;
	while (ii <(*npoints) && jj <(*npoints)) {
		if (points[jj] != NULL) {
			points[ii] = points[jj];
			ii++;
		}
		jj++;
	}
	(*npoints) = ii;
		
	free(IDInds);
	free(iidDelete);
}

//grow linked points from the end points.
void growFromEndPoints(int nx, int ny, int nz, float *im3d, int *npoints, LinkedPoint **points)
{
	int i,j,ii,jj,kk,i3,j3;
	float x1,x2,y1,y2,z1,z2,r1,r2,dd;
	int ne;
	int mmaxID,*IDInds;	
	int xs,*xe,ys,*ye,zs,*ze,*re,i1,j1,k1,k;
	float *epsVec;
	int *iid;
	int ntot;
	LinkedPoint *PP;
	float dot, mmin, angle;
	int x,y,z,r;
	int fs; 
	int *flagOcc;
	int *xg,*yg,*zg,*rg;
	int np;
	float dr,dr2,ddr,drmax,dz,dd2;
	float threshold; 
	int nxym, flag;
	int xo,yo,zo,ro,f1,f2,f3;
	float ao;
	float frex, zex, zexc, frsMax, frsMin;
	float vx,vy;
	long int k3;
	float *im3dS, *dists3d;
	int nb,xc,yc,zc;
	int nxx,nyy,nzz;
	int nid, icon, icon2;
	int xss,yss,zss,xee,yee,zee;
	int x11, y11, z11, r11, id1, x22, y22, z22, r22, id2;
	int ip;
	float d1, d2;
	int flagContinued;
	int boxSize, nptot;
	float rmax;
	int maxTry = 10, ntryext=0;	//maximum number of trials for extending. 
	int i1start,i1end;
	float drs, dr2s, dy, dy2;
	int ii2, niter;
	int js[4];
	//for searching next candidate points with angles. 
	float ang, dang, an, anW;
	int ndiv = 64, iang, iam;	//PARAMETERS!
	int xx[ndiv],yy[ndiv];
	float dds[ndiv];
	
	#if growFromEndPoints_DIAG
	int nout = 151;	//size of patch of image to be plotted, make it an odd number. 
	float imOut[nout*nout];	//image for the extension. 
	#endif 
		
	printf("Growing from the end points...\n");
	
	//delete isolated points. They can disrup growth. 
	deleteIsolatedPoints(npoints, points);
	
	ntot = nx*ny*nz;
	iid = (int *) malloc(2*(*npoints) * sizeof(int));
	flagOcc = (int *) malloc(ntot * sizeof(int));
	for (i=0; i<ntot; i++) flagOcc[i] = -1;	//inilialize. 

	//mark the points occupied by the SWC points.
	markPixelsOccupied(nx, ny, nz, flagOcc, (*npoints), points);	

	//setup a table for getting the SWC point from id. 
	mmaxID = 0;
	for (i=0;i<(*npoints);i++) {
		if (mmaxID < points[i]->ID) mmaxID = points[i]->ID;
	}
	IDInds = (int *) malloc((mmaxID+1) * sizeof(int));
	for (i=0; i<(*npoints); i++) IDInds[points[i]->ID] = i;
	nid = mmaxID + 1;
	//get the end points. 
	ne = 0;
	rmax = FLT_MIN;
	for (i=0; i<(*npoints); i++) {
		//check if the points is in bound.
		if (checkInBound(points[i],nx,ny,nz) == 0) continue;
		if (NumConnLinkedPoint(points[i]) == 1) {
			iid[ne] = i;
			ne++;
		}
		if (points[i]->r > rmax) rmax = r;
	}
	if (ne == 0) {
		free(IDInds); 
		free(iid);
		free(flagOcc);
		return;
	}
	xe = (int *) malloc(ne * sizeof(int));
	ye = (int *) malloc(ne * sizeof(int));
	ze = (int *) malloc(ne * sizeof(int));
	re = (int *) malloc(ne * sizeof(int));
	epsVec = (float *) malloc(2* ne * sizeof(float));
	for (i3=0; i3<ne; i3++) {
		ii = iid[i3];
		xe[i3] = points[ii]->x;
		ye[i3] = points[ii]->y;
		ze[i3] = points[ii]->z;
		re[i3] = points[ii]->r;
		if (points[ii]->conn->next  != NULL) {
			jj = IDInds[points[ii]->conn->next->val];
			xs = points[ii]->x - points[jj]->x;
			ys = points[ii]->y - points[jj]->y;
			dd = fmax(1e-5,sqrt(xs*xs + ys*ys));
			epsVec[2*i3] = xs/dd;		//unit vector pointing from previous point to this end point. 
			epsVec[2*i3+1] = ys/dd;				
		} else {
			epsVec[2*i3] = 0;		//isolated point, no direction. 
			epsVec[2*i3+1] = 0;				
		}	
	}
	nxym = (int)(nx * ny * PARAMS.maxFracTotPoints);
	xg = (int *) malloc(nxym * sizeof(int));
	yg = (int *) malloc(nxym * sizeof(int));
	zg = (int *) malloc(nxym * sizeof(int));
	rg = (int *) malloc(nxym * sizeof(int));	

	frsMax = PARAMS.factSearchMax;		 
	frsMin = PARAMS.factSearchMin;	
	ddr = 2;	//the widths of the shell for seaching.
	nptot = 0;
	
	for (i3=0; i3<ne; i3++) {
		x1 = xe[i3];
		y1 = ye[i3];
		z1 = ze[i3];
		r1 = re[i3];
		ii = iid[i3];
		PP = points[ii];	
		ip = ii;
		np = 0;
		vx = epsVec[2*i3];
		vy = epsVec[2*i3+1];
		x11 = x1; y11 = y1; z11 = z1; r11 = r1; id1 = ii;
		icon = -1;
		while (1) {
			#if growFromEndPoints_DIAG
			printf("\nExtending from x1=%f y1=%f z1=%f r1=%f...\n",x1,y1,z1,r1);
			//copy the image. 
			for (i1=0; i1<nout*nout; i1++) imOut[i1] = 1;
			float mmin2 = FLT_MAX;
			float mmax2 = FLT_MIN;
			for (i1 = fmax(0,x1-nout/2); i1 < fmin(nx,x1+nout/2+1); i1++) {
				for (j1 = fmax(0,y1-nout/2); j1 < fmin(ny,y1+nout/2+1); j1++) {
					k1 = getIndex(nout/2 + i1 - x1, nout/2 + j1 - y1, nout, nout);
					if (k1 <0 || k1 >= nout*nout) continue;
					imOut[k1] = im3d[getIndex3DZfirst(i1,j1,z1,nx,ny,nz)];
					if (imOut[k1] < mmin2) mmin2 = imOut[k1];
					if (imOut[k1] > mmax2) mmax2 = imOut[k1];
				}
			}
			//plot connected SWC points in the view. 
			int xcc, ycc;
			float dd;
			float rr;
			int ii2;
			for (k1 = np-1; k1 >=0; k1--) {
				xcc = xg[k1] - x1 + nout/2;	//coorindates in the image
				ycc = yg[k1] - y1 + nout/2;
				rr = rg[k1];
				if (xcc < 0 || xcc >= nout || ycc <  0 || ycc >= nout) continue; //not in the view
				for (i1=fmax(0,xcc-rr); i1 <=fmin(nout-1,xcc+rr); i1++) {
					for (j1=fmax(0,ycc-rr); j1 <=fmin(nout-1,ycc+rr); j1++) {
						ii2 = getIndex(i1,j1,nout,nout);
						if (ii2 <0 || ii2 >= nout * nout) continue;
						dd = sqrt((i1 - xcc)*(i1 - xcc) + (j1-ycc)*(j1-ycc));
						if (abs(dd - rr) < 0.1 || dd < 0.1) {
							imOut[ii2] = mmin2-0.01;
						}
					}
				}
			}	
			int idPre;
			LinkedPoint *pn;
			idPre = PP->ID;
			pn = PP;
			while (1) {
				xcc = pn->x - x1 + nout/2;	//coorindates in the image
				ycc = pn->y - y1 + nout/2;
				rr = pn->r;
				if (xcc < 0 || xcc >= nout || ycc <  0 || ycc >= nout) {
					printf("Linked point not in view.\n");
					break; //not in the view
				}
				for (i1=fmax(0,xcc-rr); i1 <=fmin(nout-1,xcc+rr); i1++) {
					for (j1=fmax(0,ycc-rr); j1 <=fmin(nout-1,ycc+rr); j1++) {
						ii2 = getIndex(i1,j1,nout,nout);
						if (ii2 <0 || ii2 >= nout * nout) continue;
						dd = sqrt((i1 - xcc)*(i1 - xcc) + (j1-ycc)*(j1-ycc));
						if (abs(dd - rr) < 0.1 || dd < 0.1) {
							imOut[ii2] = mmin2-0.01;
						}
					}
				}
				if (NumConnLinkedPoint(pn) == 1 || pn->conn->next->val == idPre) {
					printf("Linked point reached end or crossing point.\n");
					break;
				}				
				if (pn->conn->next->val != idPre) {
					idPre = pn->ID;
					pn = points[IDInds[pn->conn->next->val]];
				} else if (NumConnLinkedPoint(pn) == 2) {
					idPre = pn->ID;
					pn = points[IDInds[pn->conn->next->next->val]];
				} 		
			}					
			#endif 			
			flagContinued = 0;
			dr2 = frsMin * r1;	 //minium radius	
			drmax = (int) (r1 * frsMax) + 1;
			ddr = (int)fmax(2,(int)drmax*1.0/maxTry+1);	//try maxTry times at most. 
			dr = dr2 + ddr;			//maximum radius
			dz = 0;		//smaller z, PARAMETER.
			//get a small volume and compute shortest path.
			xss = (int)fmax(0, x1-drmax);
			xee = (int)fmin(nx,x1+drmax);
			yss = (int)fmax(0, y1-drmax);
			yee = (int)fmin(ny,y1+drmax);
			zss = (int)fmax(0, z1-dz);
			zee = (int)fmin(nz,z1+dz+1);
			nxx = xee - xss;
			nyy = yee - yss;
			nzz = zee - zss;
			//allocate memory.
			nb = nxx * nyy * nzz; 
			//allocate	
			im3dS = (float *) malloc(nb * sizeof(float));
			dists3d = (float *) malloc(nb * sizeof(float));
			for (i1=0; i1<nxx; i1++) {
				for (j1=0; j1<nyy; j1++) {
					for (k1=0; k1<nzz; k1++) {
						kk = getIndex3D(i1,j1,k1,nxx,nyy,nzz);
						im3dS[kk] = im3d[getIndex3DZfirst(i1+xss,j1+yss,k1+zss,nx,ny,nz)];
						dists3d[kk] = 0.0;
					}
				}
			}
			xc = x1 - xss;
			yc = y1 - yss;
			zc = z1 - zss;
			//smooth.
			for (k1=0; k1<nzz; k1++) {
				gaussianFilter2D(nxx,nyy,7,7,im3dS+k1*(nxx*nyy),im3dS+k1*(nxx*nyy),2.0);	//PARAMETER!
			}
			//get the distance map. 			
			dijstraComputeDists3D(xc, yc, zc, nxx, nyy, nzz, PARAMS.zfact, im3dS, dists3d, 0);
			for (ntryext=0; ntryext < maxTry; ntryext++) {
				//check if the searching radius reaches boundary. If so, stop.
				if (x1 < dr || nx - x1 < dr || y1 < dr || ny - y1 < dr) {
					#if growFromEndPoints_DIAG
					printf("Search region near the boundary. Skip. \n");
					#endif
					break;
				}

				//find the candidate point.  
				ang = atan2(vy,vx);
				anW = PARAMS.angle/180.0*3.1415926;
				dang = 2*anW/ndiv;	//angle increments.
				for (iang = 0; iang < ndiv; iang++) {
					an = ang - anW + dang * iang;
					i1 = (int) (xc + dr2 * cos(an));
					j1 = (int) (yc + dr2 * sin(an));
					dds[iang] = 0;
					for (k1=0; k1<nzz; k1++) {
						jj = getIndex3D(i1,j1,k1,nxx,nyy,nzz);
						dds[iang] += dists3d[jj];
					}
					xx[iang] = i1;
					yy[iang] = j1;
				}
				#if growFromEndPoints_DIAG
				printf("Saving the profile in the ring to temp4.py\n");
				FILE *fpt;
				fpt = fopen("temp4.py","w");
				fprintf(fpt,"from pylab import *\n");
				fprintf(fpt,"a=[");
				for (iang=0;iang<ndiv-1; iang++) {
					an = - anW + dang * iang;
					fprintf(fpt,"%3.2f, ",an);
				}
				an = - anW + dang * iang;
				fprintf(fpt,"%3.2f]\n",an);
				fprintf(fpt,"b=[");
				for (iang=0;iang<ndiv-1; iang++) {
					fprintf(fpt,"%3.2f, ",dds[iang]);
				}
				fprintf(fpt,"%3.2f]\n",dds[iang]);						
				fprintf(fpt,"figure(22); plot(a,b)\n");
				#endif
				
				//smooth, find the minimum.
				gaussianFilter2D(ndiv,1,12,12,dds,dds,4.0);	//PARAMETERS!
				mmin = FLT_MAX;
				for (iang = 0; iang < ndiv; iang++) {
					if (dds[iang] < mmin) {
						mmin = dds[iang];
						iam = iang;
					}
				}
				x = xx[iam] + xss;
				y = yy[iam] + yss;
				z = z1;
				#if growFromEndPoints_DIAG
				fprintf(fpt,"c=[");
				for (iang=0;iang<ndiv-1; iang++) {
					fprintf(fpt,"%3.2f, ",dds[iang]);
				}
				fprintf(fpt,"%3.2f]\n",dds[iang]);						
				fprintf(fpt,"figure(22); plot(a,c)\n");
				fprintf(fpt,"xlabel('angle')\n");
				fprintf(fpt,"ylabel('weighted distance')\n");
				fclose(fpt);	
				#endif
				//end of searching, comes out fs, x, y, z
				
				f1 = 0; f2 = 0; f3 = 0;
				//adjust z of the new point.
				adjustZOfLinkedPoint(x, y, z, r1, nx, ny, nz, im3d, &zo);
				z = zo;
				f1 = checkValidityOfLinkedPoint(x,y,z,r1,nx,ny,nz,im3d,&xo,&yo,&ro,&ao);			//do this twice to correc the position of the point. 
				if (f1 == 1) {
					adjustZOfLinkedPoint(xo, yo, z, ro, nx, ny, nz, im3d, &zo);
					z = zo;
					f2 = checkValidityOfLinkedPoint(xo,yo,z,ro,nx,ny,nz,im3d,&xo,&yo,&ro,&ao);
				}
				if (f2 == 1) {
					adjustZOfLinkedPoint(xo, yo, z, ro, nx, ny, nz, im3d, &zo);
					z = zo;
					f3 = checkValidityOfLinkedPoint(xo,yo,z,ro,nx,ny,nz,im3d,&xo,&yo,&ro,&ao);
				}
				#if growFromEndPoints_DIAG
				xcc = xo - x1 + nout/2;	//coorindates in the image
				ycc = yo - y1 + nout/2;
				rr = ro;
				if (xcc < 0 || xcc >= nout || ycc <  0 || ycc >= nout) {
					break; //not in the view
				}
				for (i1=fmax(0,xcc-rr); i1 <=fmin(nout-1,xcc+rr); i1++) {
					for (j1=fmax(0,ycc-rr); j1 <=fmin(nout-1,ycc+rr); j1++) {
						ii2 = getIndex(i1,j1,nout,nout);
						if (ii2 <0 || ii2 >= nout * nout) continue;
						dd = sqrt((i1 - xcc)*(i1 - xcc) + (j1-ycc)*(j1-ycc));
						if (abs(dd - rr) < 0.1 || dd < 0.1) {
							imOut[ii2] = mmax2+0.01;
						}
					}
				}
				#endif
				
				//check if the new point is overlapping with the existing points. This can happen because the position of the point is adjusted. 
				if (f1 * f2 * f3 == 1) {
					icon = flagOcc[getIndex3DZfirst(xo,yo,z,nx,ny,nz)]; 
				}
				if (icon != -1) {
					#if growFromEndPoints_DIAG
					printf("New point crossing exsisting points. Stop growing.\n");
					#endif 
					break;
				}		
				
				//check the validity of the point. The center of the point can change. Double check. 
				if (f1 * f2 * f3 == 1 &&
					xo > ro && xo < nx-ro-1 && yo > ro && yo < ny-ro-1) {
					//get the radius and adjust the center. 
					//add new point. 
					xg[np] = xo;
					yg[np] = yo;
					zg[np] = z;
					rg[np] = ro;
					np++;

					//update the vector
					xs = (xo - x1); ys = (yo - y1);
					dd = fmax(1e-5,sqrt(xs*xs + ys*ys));
					vx = xs/dd;
					vy = ys/dd;
					//mark pixels as occupied. 
					x22 = xo; y22 = yo; z22 = z; r22 = ro; id2 = (*npoints)+np-1;
					markPixelsOccupiedTwoPoints(nx, ny, nz, flagOcc, id1, x11, y11, z11, r11, id2, x22, y22, z22, r22);
					x11 = x22; y11 = y22; z11 = z22; r11 = r22; id1 = id2;
					//continue from this point.
					x1 = xo;
					y1 = yo; 
					z1 = z;
					r1 = ro;
					nid++;	//increase the maxID. 
					flagContinued = 1;							
					break;		
				} else {				
					#if growFromEndPoints_DIAG
					if (f1 * f2 * f3 ==0) {
						printf("Extension end. Failed validity test.\n");
						printf("f1 = %d f2 = %d f3 = %d\n",f1,f2,f3);
						printf("x = %d y = %d z = %d r = %d \n",(int)x,(int)y,(int)z,(int)r1);
						printf("xo = %d yo = %d zo = %d ro = %d \n",(int)xo,(int)yo,(int)z,(int)ro);
					} else {
						printf("Extension end. Too close to the boundary\n");
						printf("x,y,z,r=%d %d %d %d\n",x,y,z,r);
					}
					#endif				 
					//try again with larger radius
					dr2 = dr;	 //minium radius	
					dr = dr2 + ddr;			//maximum radius
				}
				#if growFromEndPoints_DIAG
				printf("extension failed. try dr=%f \n",dr);		
				printf("Saving local image to temp.growth.tif...");
				createGrayTiffStackFromFloatArray(nout,nout,1,imOut,"temp.growth.tif");
				printf("Press any key..."); getchar();
				#endif
				if (dr > drmax) {
					break;
				}
			}
			//free memory. 
			free(im3dS);
			free(dists3d);
			if (flagContinued == 0) break; //stop growth.
		}	
		if (np > 0) { 
			printf("...added %d points.",np);
			nptot += np;
			//smooth z, r, x, y
			if (np > 5) {
				//gaussianFilter2D(1,np,1,5,(float *)xg,(float *)xg,1);
				//gaussianFilter2D(1,np,1,5,(float *)yg,(float *)yg,1);
				//gaussianFilter2D(1,np,1,5,(float *)zg,(float *)zg,1);
				gaussianFilter2D(1,np,1,3,(float *)rg,(float *)rg,1);
			}
			//create linked points. 
			for (kk=0; kk<np; kk++) {		
				mmaxID += 1;		//this is the new id.	
				points[(*npoints)] = CreateLinkedPoint(xg[kk],yg[kk],zg[kk],rg[kk],mmaxID,3);			
				AddConnLinkedPoint(PP,points[(*npoints)]->ID);
				AddConnLinkedPoint(points[(*npoints)],PP->ID);
				PP = points[(*npoints)];
				ip = (*npoints);
				(*npoints)++;
			}
			//update IDInds
			free(IDInds);
			IDInds = (int *) malloc((mmaxID+1) * sizeof(int));
			for (i=0; i<(*npoints); i++) IDInds[points[i]->ID] = i;			
		}	
		
		//connect the end point to nearby exisiting point.
		if (icon >= 0 && icon < (*npoints) && ip != icon) {
			//mark pixels occupied
			x11 = PP->x; y11 = PP->y; z11 = PP->z; r11 = PP->r;
			x22 = points[icon]->x; y22 = points[icon]->y; z22 = points[icon]->z; r22 = points[icon]->r;
			icon2 = markPixelsOccupiedTwoPoints(nx, ny, nz, flagOcc, ip, x11, y11, z11, r11, icon, x22, y22, z22, r22);
			//icon2 is the closest SWC point to ip. connect to it. 
			d1 = sqrt((x22 - x11)*(x22 - x11) + (y22 - y11)*(y22 - y11) + 
					  (PARAMS.zfact * PARAMS.zfact) * (z22 - z11) * (z22 - z11));
			d2 = sqrt((points[icon2]->x - x11)*(points[icon2]->x - x11) + (points[icon2]->y - y11)*(points[icon2]->y - y11) + 
					  (PARAMS.zfact * PARAMS.zfact) * (points[icon2]->z - z11) * (points[icon2]->z - z11));	
			if (d2 < d1) icon = icon2;
			//check z. 
			if (PARAMS.zfact * abs(z22 -z11) <= PARAMS.zJumpFact * (r11 + r22)) {
				AddConnLinkedPoint(points[icon],PP->ID);
				AddConnLinkedPoint(PP,points[icon]->ID);
			}
		}

	}			
	printf(" Total points added: %d\n",nptot);

	//free memory
	free(xe);
	free(ye);
	free(ze);
	free(re);	
	free(epsVec);
	free(iid);
	free(IDInds);
	free(flagOcc);
	free(xg);
	free(yg);
	free(zg);
	free(rg);
}
 
DLinkedListSWC* GetLastInDListSWC(DLinkedListSWC *dlist)
{
	DLinkedListSWC *p;
	p = dlist;
	if (p == NULL) return NULL;
	while (p->next != NULL) {
		p = p->next;
	}
	return p;	
}

void AppendToDListSWC(DLinkedListSWC **dlist, LinkedPoint *P) 
{
	DLinkedListSWC *new,*last;
	new = (DLinkedListSWC *) malloc(sizeof(DLinkedListSWC));	
	new->P.x = P->x;
	new->P.y = P->y;
	new->P.z = P->z;
	new->P.r = P->r;
	new->P.ID = P->ID;
	new->P.Type = P->Type;
	new->next = NULL;
	new->prev = NULL;
	if ((*dlist) == NULL) {
		new->P.parentID = -1;
		(*dlist) = new;
	} else {
		last = GetLastInDListSWC((*dlist));
		last->next = new;
		new->prev = last;
		new->P.parentID = last->P.ID;
	}
	return;
}
	
void DeleteDListSWC(DLinkedListSWC *dlist)
{
	DLinkedListSWC *p, *next;
	if (dlist == NULL) return;
	p = dlist;
	while (p != NULL) {
		next = p->next;
		free(p);
		p = next;
	}
	dlist = NULL;
}

void DeleteFromListSWC(DLinkedListSWC *dlist, DLinkedListSWC *pdel)
{
	if (pdel->prev == NULL && pdel->next == NULL) {
		dlist = NULL;
	} else {
		if (pdel->prev == NULL) {
			pdel->next->prev = NULL;
		} else if (pdel->next == NULL) {
			pdel->prev->next = NULL;
		} else {
			pdel->prev->next = pdel->next;
			pdel->next->prev = pdel->prev;
		}
	}		
	free(pdel);	
}


// this function gets all branches in the linked points starting from jj.
int getAllChildBranchesLinkedPoints(int jj,int npoints,LinkedPoint **points,int *IDInds,int *flags, int *nbrp, DLinkedListSWC **branches, LinkedList **branchConnIDs)
{
	DLinkedListSWC *br, *biter;
	int i,kk,iid,bID;
	LinkedList *iter;
	DLinkedListSWC *last;
	LinkedList *brIDs;
	int nn=0, njj[1000];	//maximum number of connections is assumed to be 1000.  
		
	br = NULL;
	last = NULL;
	brIDs = NULL;
	biter = NULL;
	iter = NULL;
	while (1) {
		AppendToDListSWC(&br,points[jj]);
		flags[jj] = 1;
		// find the next point.
		nn = 0;
		iter = points[jj]->conn->next;
		while (iter != NULL) {
			kk = iter->val;
			jj = IDInds[kk];
			if (flags[jj] == 0) {
				njj[nn] = jj;
				nn++;
			}				
			iter = iter->next;
		}	
		if (nn ==0) {
			break;
		} else if (nn == 1) {
			jj = njj[0];
		} else {
			for (i=0; i<nn; i++) {
				iid = njj[i];
				if (flags[iid] == 1) continue; //this point could have been used through the recursion. Recheck the flags. 
				bID = getAllChildBranchesLinkedPoints(iid,npoints,points,IDInds,flags,nbrp,branches,branchConnIDs);
				if (bID != -1) {
					AppendToList(&brIDs,bID);
				}
			}
			break;
		}
	}
	branches[*nbrp] = br;
	branchConnIDs[*nbrp] = brIDs;
	*nbrp += 1;			
	return *nbrp-1;
}

//compute the branch tree object size
int compBranchObjSize(int ii,LinkedList **branchConnIDs,int *brLen, int no)
{
	int jj;
	LinkedList *brID;
	
	no += brLen[ii];
	brID = branchConnIDs[ii];
	while (brID != NULL) {
		jj = brID->val;
		no = compBranchObjSize(jj,branchConnIDs,brLen,no);
		brID = brID->next;
	}
	return no;
}

//mark delete branch tree object. 
void markDeleteBranches(int ii,LinkedList **branchConnIDs,int *delMark)
{
	int jj;
	LinkedList *brID;
	
	delMark[ii] = 1;
	brID = branchConnIDs[ii];
	while (brID != NULL) {
		jj = brID->val;
		markDeleteBranches(jj,branchConnIDs,delMark);
		brID = brID->next;
	}
}

//make the connectivity a tree structure, delete small leafs. 
void regularizeAndDeleteSmallBranches(int *npoints, LinkedPoint **points)
{
	float *D, threshold;
	int *E, *IDInds, maxID, i, j, ii, jj;
	LinkedList *iter;
	int nlen, ip, flag, nd, idP;
	int *iidDelete;
	
	D = (float *) malloc((*npoints)*(*npoints) * sizeof(float));
	E = (int *) malloc((2*(*npoints)+1) * sizeof(int));
	threshold = 1;
	for (i=0; i<(*npoints)*(*npoints); i++) D[i] = 2*threshold;
	for (i=0; i<2*(*npoints)+1;i++) E[i] = 0;
	 
	maxID = 0; 	
	for (i=0; i<(*npoints); i++) {
		if (maxID < points[i]->ID) maxID = points[i]->ID;
	}
	IDInds = (int *)malloc((maxID+1) * sizeof(int));	
	for (i=0; i<(*npoints); i++) {		
		IDInds[points[i]->ID] = i;
	}
	
	for (ii=0; ii<(*npoints); ii++) {
		iter = points[ii]->conn->next;
		while (iter != NULL) {
			jj = IDInds[iter->val];
			D[ii*(*npoints)+jj] = 0.5 * threshold;	//set the distance between the connected nodes.
			D[jj*(*npoints)+ii] = 0.5 * threshold;	//make it symmetric		
			iter = iter->next;
		}
	}
	//compute the minium spanning tree.
	kruskal((*npoints), threshold, D, E);
	
	//reset the connections. 
	//delete all connections
	for (i=0; i<(*npoints); i++) {
		DeleteList(points[i]->conn->next);
		points[i]->conn->val = 0;
		points[i]->conn->next = NULL;
	}
	for (i=0; i<E[2*(*npoints)]; i++) {
		ii = E[2*i];
		jj = E[2*i+1];
		AddConnLinkedPoint(points[ii],points[jj]->ID);
		AddConnLinkedPoint(points[jj],points[ii]->ID);	//symmatric connection
	}	
	
	//now delete small leaf branches. 
	iidDelete = (int *) malloc((*npoints)*sizeof(int));
	for (nlen = 1; nlen<PARAMS.minNumPointsBr; nlen++) {
		for (i=0; i<(*npoints); i++) {		
			IDInds[points[i]->ID] = i;
		}
		nd = 0;
		for (i=0; i<(*npoints); i++) {
			if (NumConnLinkedPoint(points[i]) == 0) {//isolated point.
				iidDelete[nd++] = i; 
				continue;
			}		
			if (NumConnLinkedPoint(points[i]) > 1) continue;	//this is not a terminal point. 
			for (ip = 0; ip<2; ip++) {	//two passes.
				if (ip == 0) {
					flag = 0;
				} else {
					if (flag == 0) break;
				}
				idP = points[i]->ID;
				ii = IDInds[points[i]->conn->next->val];
				for (j=0; j<nlen; j++) {	//chase.
					if (NumConnLinkedPoint(points[ii]) > 2)	{//branhcing point.
						if (ip == 0) {
							flag = 1;
						} 
						break;
					}
					if (NumConnLinkedPoint(points[ii]) == 1)	{//end point.
						if (ip == 0) {
							flag = 1;
						} else {
							iidDelete[nd++] = ii;
						}
						break;
					}
					if (ip == 1) iidDelete[nd++] = ii;
					//get to next point. 
					if (points[ii]->conn->next->val == idP) {
						idP = points[ii]->ID;
						ii = IDInds[points[ii]->conn->next->next->val];
					} else {
						idP = points[ii]->ID;
						ii = IDInds[points[ii]->conn->next->val];
					}
				}
			}
			if (flag == 1) iidDelete[nd++] = i;					
		}
		printf("Delete %d points. \n",nd);
		//delete connections. 
		for (i=0; i<(*npoints); i++) {
			for (j=0; j<nd; j++) {
				DelConnLinkedPoint(points[i],points[iidDelete[j]]->ID);
			}
		}
		//delete the points
		for (j=0; j<nd; j++) {
			if (points[iidDelete[j]] != NULL) {
				DeleteLinkedPoint(points[iidDelete[j]]);
				points[iidDelete[j]] = NULL;
			}
		}
		//make the remaining points consecutive
		ii = 0; jj = 0;
		while (ii <(*npoints) && jj <(*npoints)) {
			if (points[jj] != NULL) {
				points[ii] = points[jj];
				ii++;
			}
			jj++;
		}
		(*npoints) = ii;
	}
		
	free(E);
	free(D);
	free(IDInds);
	free(iidDelete);
}

//create SWC from linked points. 
void createSWC(int *npoints, LinkedPoint **points, char *outName)
{
	int *IDInds, *numConn, *flags, i, ii, jj, nbr, mmin;
	DLinkedListSWC **branches, *br, *iter;
	LinkedList **branchConnIDs, *brID;
	int maxID = 0;
	FILE *fid;
	int *brLen, no, *delMark;
	int lastID,flag,nt;

	//regularize and delete small branches. 
	regularizeAndDeleteSmallBranches(npoints, points);

	numConn = (int *)malloc((*npoints) * sizeof(int));
	flags = (int *)malloc((*npoints) * sizeof(int));
	branchConnIDs = (LinkedList **)malloc((*npoints) * sizeof(LinkedList *));
	branches = (DLinkedListSWC **) malloc((*npoints) * sizeof(DLinkedListSWC *));
	for (i=0; i<(*npoints); i++) {
		if (maxID < points[i]->ID) maxID = points[i]->ID;
		branches[i] = NULL;
		branchConnIDs[i] = NULL;
	}
	IDInds = (int *)malloc((maxID+1) * sizeof(int));
	
	for (i=0; i<(*npoints); i++) {		
		IDInds[points[i]->ID] = i;
		numConn[i] = NumConnLinkedPoint(points[i]);
		flags[i] = 0;
	}
	nbr = 0;	
	while (1) {
		jj = -1;
		mmin = 1e5;		
		for (i=0; i<(*npoints); i++) {
			if (flags[i] == 0 && mmin > numConn[i]) {
				mmin = numConn[i];
				jj = i;
			}
		}
		if (jj == -1) break;
		getAllChildBranchesLinkedPoints(jj,(*npoints),points,IDInds,flags,&nbr,branches,branchConnIDs);
	}
	
	// connect branches
	for (ii=0; ii<nbr; ii++) {
		br = branches[ii]; 
		brID = branchConnIDs[ii];
		while (brID != NULL) {
			jj = brID->val;
			if (branches[jj] == NULL) continue;
			branches[jj]->P.parentID = GetLastInDListSWC(br)->P.ID;
			brID = brID->next;
		}
	}
 
	//delete small isolated objects. 
	printf("Delete isolated branches with number of points smaller than %d\n",PARAMS.minNumPointsBrIso);
	brLen = (int *)malloc(nbr * sizeof(int));
	delMark = (int *)malloc(nbr * sizeof(int));	
	for (ii=0; ii<nbr; ii++) {
		brLen[ii] = 0;
		delMark[ii] = 0;
		iter = branches[ii];
		while (iter != NULL) {
			brLen[ii]  += 1;
			lastID = iter->P.ID;
			iter = iter->next;
		}
		if (brLen[ii] < PARAMS.minNumPointsBrIso) {	//PARAMETER!!!
			if (NumConnLinkedPoint(points[IDInds[branches[ii]->P.ID]]) == 1 && 
				NumConnLinkedPoint(points[IDInds[lastID]]) == 1) { //isolated branch. 
				delMark[ii] = 1;
			}	
		}
	}
	//get tree objects from roots and see if the size is small. 
	nt = 0;
	for (ii=0; ii<nbr; ii++) {
		if (branches[ii]->P.parentID == -1) {	//compute size of the tree from this branch, which has a root, and delete all connected branches if the total size is small. 
			//no = compBranchObjSize(ii,branchConnIDs,brLen,0);
			//printf("Root branch id=%d objectsize=%d\n",ii,no);
			nt++;
		}
	}
 
	printf("Saving SWC to %s\n",outName);
	fid = fopen(outName,"w");
	no = 0;
	for (i=0; i<nbr; i++) {
		if (delMark[i] == 1) continue;	//do not record isolated branches.
		no++; 
		iter = branches[i];
		while (iter != NULL) {
			fprintf(fid,"%d %d %d %d %d %d %d\n",iter->P.ID,iter->P.Type,(int)iter->P.y,(int)iter->P.x,(int)iter->P.z,(int)iter->P.r,iter->P.parentID);
			iter = iter->next;
		}
	}
	fclose(fid);
	printf("There %d segments %d trees\n",no,nt);
	
	free(IDInds);
	free(numConn);
	free(flags);
	for (i=0; i<nbr; i++) {
		DeleteDListSWC(branches[i]);
		DeleteList(branchConnIDs[i]);
	}
	free(branchConnIDs);	
	free(branches);
	free(brLen);
}


