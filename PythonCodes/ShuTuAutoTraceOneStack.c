/*
 *  This file is for automatically tracing neurons using 2D project and tracing method. 
 *  Used for tracing neurites in one tiff stack.
 *  Written by Dezhe Z. Jin, dzj2@psu.edu
 *  06/02/2015, (C) Dezhe Jin
 * 
 * 	This program is free software: you can redistribute it and/or modify
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

#define OUTPUTMASK 1	//if set to 1, the mask will be saved. 

//trace neurites in a single tif stack. first argument, tif filename; second argument, parameter filename.   
int main(int argc, char *argv[])
{
	char *filename, *filenameParam, outName[1000];
	int nx, ny, nz, i,j,k,ntot,nxy;
	float *im3d,*imFlat,maxI,minI;
	int islab,z1,z2,z1ext,z2ext;		
	float *bw;
	int tiffinfo[5];
	LinkedPoint **points;
	int npoints;
	int ind, npointsTot;
	int zex;
	int minNumBr;
	long int ii;
	#ifdef OUTPUTMASK
	float *bwCombined;
	#endif 
	
	setbuf(stdout,NULL);	//make print out immeidate, no buffering. 
	
	if (argc > 1) {
		filename = argv[1];
		printf("Using  tiff image filename = %s\n",filename);
	} else {
		printf("ERROR: please sepcify the tiff stack filename as the first commandline argument.\n");
		return 1;
	}	  
	if (argc > 2) {				//read parameters from file
		filenameParam = argv[2];
		readParametersFromFile(filenameParam);
	} else {
		printf("Using default parameters.\n");
		setDefaultParameters();
	}
	printParameters();

	heap = NULL; //initialize the heap pointer. 
	 
	//read image stack. 
	probeFileDimensions(filename,&nx,&ny,&nz,tiffinfo);
	ntot = nx * ny * nz;	
	nxy = nx * ny;
	im3d = (float *) malloc(ntot*sizeof(float));	
	imFlat = (float *) malloc(nxy * sizeof(float)); 
	readImageAndReturnIm3d(filename,nx,ny,nz,im3d,tiffinfo); //this returns normalized im3d

	bw = (float *) malloc(nxy * sizeof(float));
	#ifdef OUTPUTMASK
	bwCombined = (float *) malloc(nxy * sizeof(float));
	#endif

	//pointer arrary of LinkedPoints used in creating swc file. 
	
	npointsTot = (int)(ntot * PARAMS.maxFracTotPoints);		//PARAMETER!
	points = (LinkedPoint **) malloc(npointsTot * sizeof(LinkedPoint *));
	for (i=0; i<npointsTot; i++) points[i] = NULL;
	npoints = 0;
	
	#ifdef OUTPUTMASK
	//initialze bwCombined
	for (i=0; i<nxy; i++) bwCombined[i] = 0.0;
	#endif
		
	zex = (int)PARAMS.zext/PARAMS.zDist;
	for (islab=0;islab <PARAMS.nSplit;islab++) {
		z1 = nz/PARAMS.nSplit * islab;
		z2 = (int) fmin(nz,nz/PARAMS.nSplit * (islab+1));	//work on sub-slab
		z1ext = (int) fmax(0, z1 - zex);				//PARAMETER!! extended sub-slab, to avoid overlap near sub slab boundary.
		z2ext = (int) fmin(nz,z2 + zex);				//PARAMETER!!
		//miniumum intensity projection. 
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				minI = FLT_MAX;
				for (k=z1ext; k<z2ext; k++) {
					ii = getIndex3DZfirst(i,j,k,nx,ny,nz);
					if (minI > im3d[ii]) minI = im3d[ii];
				}
				ii = getIndex(i,j,nx,ny);
				imFlat[ii] = minI;
			}
		}
		sprintf(outName,"%s.imFlat%d.tif",filename,islab);
		printf("Saving maximum projection for slab %d to %s...\n",islab,outName);
		createGrayTiffStackFromFloatArray(nx,ny,1,imFlat,outName);	
				
		//create mask. 
		createMask2D(nx,ny,imFlat,bw);
		sprintf(outName,"%s.imFlat%d.Mask.tif",filename,islab);
		printf("Saving mask for slab %d to %s...\n",islab,outName);
		createGrayTiffStackFromFloatArray(nx,ny,1,bw,outName);	
		
		//get the distance map, save in imFlat.
		sedt(nx, ny, bw, imFlat);
		sprintf(outName,"%s.imFlat%d.distMap.tif",filename,islab);
		printf("Saving distancemap for slab %d to %s...\n",islab,outName);
		createGrayTiffStackFromFloatArray(nx,ny,1,bw,outName);			
				
		//get skeleton. 
		//printf("Getting skeleton...\n");
		skeletonizeZS(nx,ny,bw,bw);
		
		//prune small branches. 
		printf("Pruning small branches...\n");
		pruneSmallBranches(nx,ny,bw,PARAMS.smallLen); 	

		sprintf(outName,"%s.imFlat%d.Skeleton.tif",filename,islab);
		printf("Saving skeleton for slab %d to %s...\n",islab,outName);
		createGrayTiffStackFromFloatArray(nx,ny,1,bw,outName);	

		//get Linked points from the skeleton.
		//getLinkedPoints
		printf("Getting linked points...\n");
		getLinkedPointsFromSkeletonAnd3D(nx,ny,nz,bw,imFlat,im3d,z1,z2,z1ext,z2ext,&npoints,points,PARAMS.zfact);
		#ifdef OUTPUTMASK
		for (i=0; i<nxy; i++) {
			if (bw[i] > 0.0) bwCombined[i] = 1.0;
		}		
		#endif
	}
	printf("npoints=%d\n",npoints);

	#ifdef OUTPUTMASK
	sprintf(outName,"%s.imFlat.Mask.tif",filename);
	printf("Saving the combined mask to %s...\n",outName);
	createGrayTiffStackFromFloatArray(nx,ny,1,bwCombined,outName);
	free(bwCombined);
	#endif	
	
	printf("geting points near soma...\n");
	getThickDendrites(nx,ny,nz,im3d,&npoints,points);	//get think dendrite and soma. 
	printf("npoints=%d\n",npoints);

	//free some memory. 
	free(imFlat);
	free(bw);	

	//grow from the end points. 
	growFromEndPoints(nx, ny, nz, im3d, &npoints, points);	

	//connect end points using distance. 
	printf("Connecting end points...\n");
	connectEndPoints(npoints,points,PARAMS.zfact);
		
	//create swc from the linked points. 
	printf("Creating swc file...\n");
	sprintf(outName,"%s.auto.swc",filename);
	createSWC(&npoints, points, outName);
			
	//final release memeory	
	free(im3d);
	for (i=0; i<npoints; i++) {
		DeleteLinkedPoint(points[i]);
	}
	free(points);	
	return 0;
 }
 
 
 
