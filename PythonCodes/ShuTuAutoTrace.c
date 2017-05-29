/*
 *  This file is for automatically tracing neurons using 2D project and tracing method. 
 *  Used for tracing neuron in stitched tiff stacks. 
 *  The program checks *-StitchCoord.txt for the coordinates of of the tiles. 
 *  The coordinates must have been created by stitchImages in ShuTu.py
 *  The first command line argument is filenameCommon, the tiff images are assumed to be in the format filenameCommon+number+'.tif'
 *  The second command line argument is a text file (such as NeuTuAuto.Params.txt) containing the tracing parameters.  
 * 
 *  Written by Dezhe Z. Jin, dzj2@psu.edu
 *  07/05/2015, (C) Dezhe Jin
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
#include "string.h"
#include "mpi.h"
#include <unistd.h>

#define CONNECT_ENDS 1			//set to 1 if ends are connected based on distance. 
#define FULL_RECONSTRUCTION 1	//set to 0 if only creates SWC file from previous generated linkedPoint files. 

//trace neurites in a single tif stack. first argument, tif filename; second argument, parameter filename.   
int main(int argc, char *argv[])
{
	char *filenameCommon, filename[1000], *filenameParam, outName[1000],*line;
	int tiffinfo[5],nx,ny,nz,i,j,k,ntot,nxy,ii,numTiles;
	float *im3d,*imFlat,maxI,minI;
	int *X, *Y, *Z, *fileIDs, Xmax, Ymax, nxc, nyc, ntotc, kk;
	int islab,z1,z2,z1ext,z2ext;		
	float *bw,sigma,smallLen, x, y, z, r, id;
	int ID, Type;
	unsigned int Wx, Wy;
	LinkedPoint **points, **pointsAll, *lp;
	int npoints;
	FILE *fpt;
	ssize_t read;	
	size_t len;
	int ind;
	int  numProc, rank, mlen, rc, fid; 
	char hostname[MPI_MAX_PROCESSOR_NAME];
	float *linkedP;
	int ns, npointsTot, npointsPart, npointsMax;
	LinkedList *iter, *iterpre;
	int maxID, *IDflags, IDStart,numConn, mmax;
	int *eps, ne, *IDInds, *flags, jj, x1, x2, y1, y2, nD, zex;
	float r1, r2, dd, *epsVec, xs, ys;
	MPI_Status status;
	struct valInd element;
	int *visited; 
	double *DBlock; 
	float xv,yv,dpot;
	int *flagOcc;
	int zOcc; 
	float factMarkOccXY;
	
	setbuf(stdout,NULL);	//make print out immeidate, no buffering. 
	
	//initilize MPI
	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(hostname, &mlen);
	printf ("Number of tasks= %d My rank= %d Running on %s\n", numProc,rank,hostname);

	if (argc > 1) {
		filenameCommon = argv[1];
		if (rank == 0)  printf("Tracing neurons with filenameCommon = %s\n",filenameCommon);
	} else {
		printf("ERROR: please sepcify filenameCommon as the first commandline argument.\n");
		return 1;
	}	
	//check if the coordinate file exists.
	sprintf(outName,"%sStitchCoord.txt",filenameCommon);
	fpt=fopen(outName,"r");
	X = 0; Y = 0; Z = 0; fileIDs = 0;
	if (fpt == NULL) {
		printf("Failed to open %s for reaing the title cooridnates. \nPlease check if filenameCommon is corrected specified and/or the images have been stitched using NeuTu.py. \n",outName);
		return 1;
	} else { //read the coordinates. 
		if (rank == 0) printf("Reading tile coordinates from %s \n",outName);
		if (fscanf(fpt,"%d %d %d %d\n",&numTiles,&nx,&ny,&nz) == EOF) {
			printf("Reading failed.\n");
			return 1;
		}
		if (rank == 0) printf("Number of tiles = %d tiff stack dimensions=(%d, %d, %d)\n",numTiles,nx,ny,nz);
		X = (int *) malloc(numTiles * sizeof(int));
		Y = (int *) malloc(numTiles * sizeof(int));
		Z = (int *) malloc(numTiles * sizeof(int));
		fileIDs = (int *) malloc(numTiles * sizeof(int));
		Xmax = 0; Ymax = 0;
		for (i=0; i<numTiles; i++) {
			if (fscanf(fpt,"%f %f %f %f\n",&id, &x, &y, &z)==EOF) {
				printf("Reading failed.\n");
				return 1;
			}
			fileIDs[i] = (int) id;
			X[i] = (int) x;  if (X[i] > Xmax) Xmax = X[i];
			Y[i] = (int) y;  if (Y[i] > Ymax) Ymax = Y[i];
			Z[i] = (int) z;
			//printf("%d %d %d %d \n",fileIDs[i],X[i],Y[i],Z[i]);
		}
		fclose(fpt);
	}
	if (argc > 2) {				//read parameters from file
		filenameParam = argv[2];
		readParametersFromFile(filenameParam);
	} else {
		printf("Using default parameters.\n");
		setDefaultParameters();
	}
	//print parameters.
	if (rank == 0) { 
		printParameters();
	}
	
	
	#if FULL_RECONSTRUCTION	
	npointsTot = 0;		
	zex = (int)PARAMS.zext/PARAMS.zDist;
	
	for (kk=0; kk<=(numTiles/numProc+1); kk++) {		
		fid = kk * numProc + rank;		
		if (fid >= numTiles) break;
		
		//check whether the .dat file already exists. 
		sprintf(filename,"%s%d.tif.linkedPoints.dat",filenameCommon,fileIDs[fid]);
		if (access( filename, F_OK ) != -1) {
			printf("The file for linked list exists: %s. If need to compute again, delete this file.\n",filename);
			continue;
		}		
		
		sprintf(filename,"%s%d.tif",filenameCommon,fileIDs[fid]);
		printf("Proc rank=%d Loading image from %s\n",rank,filename);		 
		//read image. 
		probeFileDimensions(filename,&nx,&ny,&nz,tiffinfo);
		ntot = nx * ny * nz;	
		nxy = nx * ny;	
		im3d = (float *) malloc(ntot*sizeof(float));	
		imFlat = (float *) malloc(nxy * sizeof(float));
		bw = (float *) malloc(nxy * sizeof(float));

		readImageAndReturnIm3d(filename,nx,ny,nz,im3d,tiffinfo); //this returns normalized im3d
		
		//allocate memory.
		//pointer arrary of LinkedPoints used in creating swc file. 
		npointsMax = (int)(ntot * PARAMS.maxFracTotPoints);		//PARAMETER!	
		points = (LinkedPoint **) malloc(npointsMax * sizeof(LinkedPoint *));
		linkedP = (float *) malloc((int)(npointsMax * (7+5)) * sizeof(float));
		for (i=0; i<npointsMax; i++) points[i] = NULL;
		npoints = 0;

		for (islab=0;islab <PARAMS.nSplit;islab++) {
			z1 = nz/PARAMS.nSplit * islab;
			z2 = (int) fmin(nz,nz/PARAMS.nSplit * (islab+1));	//work on sub-slab
			z1ext = (int) fmax(0, z1 - zex);				//PARAMETER!! extended sub-slab, to avoid overlap near sub slab boundary.
			z2ext = (int) fmin(nz,z2 + zex);				//PARAMETER!!
			//maximum intensity projection. 
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
					
			//create mask. 
			createMask2D(nx,ny,imFlat,bw);
			
			//get the distance map, save in imFlat.
			sedt(nx, ny, bw, imFlat);
					
			//get skeleton. 
			//printf("Getting skeleton...\n");
			skeletonizeZS(nx,ny,bw,bw);
			
			//prune small branches. 
			//printf("Pruning small branches...\n");
			pruneSmallBranches(nx,ny,bw,PARAMS.smallLen); 	
			
			//get Linked points from the skeleton.
			//getLinkedPoints
			printf("Getting linked points...\n");
			getLinkedPointsFromSkeletonAnd3D(nx,ny,nz,bw,imFlat,im3d,z1,z2,z1ext,z2ext,&npoints,points,PARAMS.zfact);
		}
		printf("npoints=%d\n",npoints);

		printf("geting points near soma...\n");
		getThickDendrites(nx,ny,nz,im3d,&npoints,points);	//get think dendrite and soma. 
		printf("npoints=%d\n",npoints);

		//grow from the end points. 
		growFromEndPoints(nx, ny, nz, im3d, &npoints, points);	

		printf("rank=%d npoints=%d\n",rank,npoints);
		//write the linked points to file. 
		ns = 0;
		for (i=0; i<npoints; i++) {
			linkedP[ns++] = points[i]->x + X[fid];
			linkedP[ns++] = points[i]->y + Y[fid];
			linkedP[ns++] = points[i]->z + Z[fid];
			linkedP[ns++] = points[i]->r;
			linkedP[ns++] = points[i]->ID;
			linkedP[ns++] = points[i]->Type;
			iter = points[i]->conn;
			linkedP[ns++] = iter->val;
			iter = iter->next;
			while (iter != NULL) {
				linkedP[ns++] = iter->val;
				iter = iter->next;
			}
		}
		sprintf(outName,"%s.linkedPoints.dat",filename);
		printf("Saving linked points to %s...\n",outName);
		fpt = fopen(outName,"wb");
		fwrite(&npoints,sizeof(int),1,fpt);		
		fwrite(&ns,sizeof(int),1,fpt);
		fwrite(linkedP,sizeof(float),ns,fpt);
		fclose(fpt);
		//delete memory, get ready for next tile. 
		for (i=0; i<npoints; i++) {
			DeleteLinkedPoint(points[i]);
			points[i] = NULL;
		}			
		//release memeory	
		free(points);	
		free(linkedP);
		free(bw);	
		free(im3d);
		free(imFlat);
	}	
	
	MPI_Barrier(MPI_COMM_WORLD);
		
	//construct linked points from the linked point files. 
	if (rank == 0) {
		npointsTot = 0;
		for (kk=0; kk<numTiles; kk++) {
			sprintf(filename,"%s%d.tif.linkedPoints.dat",filenameCommon,fileIDs[kk]);
			//check whether the .dat file already exists. 
			if (access( filename, F_OK ) == -1) {
				printf("The file for linked list does not exists: %s.\n",filename);
				continue;
			}		
						
			fpt = fopen(filename,"rb");
			if (fread(&npointsPart,sizeof(int),1,fpt)==EOF) {
				printf("Failed to read from %s\n",filename);
				return 1;
			}
			fclose(fpt);
			printf("Reading linked points from %s...npoints=%d\n",filename,npointsPart);			
			npointsTot += npointsPart;
		}	
		printf("Rank=%d Total number of points=%d\n",rank,npointsTot);
		npointsTot *= 10;	//allocate more memory to allow for growth. Increate this if the code fails at the final growth point.
		printf("Expanding memory to %d points to allow for growth from the end points.\n",npointsTot);
		pointsAll = (LinkedPoint **) malloc(npointsTot * sizeof(LinkedPoint *));	
		linkedP = (float *) malloc((int)(npointsTot * (7+5)) * sizeof(float));
		for (i=0; i<npointsTot; i++) pointsAll[i] = NULL;
		npoints = 0;
		maxID = -1;

		for (kk=0; kk<numTiles; kk++) {
			if (rank == 0) {
				sprintf(filename,"%s%d.tif.linkedPoints.dat",filenameCommon,fileIDs[kk]);
				//check whether the .dat file already exists. 
				if (access( filename, F_OK ) == -1) {
					printf("The file for linked list does not exists: %s.\n",filename);
					continue;
				}		
				printf("Reading linked points from %s...\n",filename);
				fpt = fopen(filename,"rb");
				if (fread(&npointsPart,sizeof(int),1,fpt) == EOF) {
					printf("Reading failed.\n");
					return 1;
				}
				if (fread(&ns,sizeof(int),1,fpt) == EOF) {
					printf("Reading failed.\n");
					return 1;
				}
				if (fread(linkedP,sizeof(float),ns,fpt) == EOF) {
					printf("Reading failed.\n");
					return 1;
				}
				fclose(fpt);
			}		
			
			sprintf(filename,"%s%d.tif",filenameCommon,fileIDs[kk]);
			//probe image. 
			probeFileDimensions(filename,&nx,&ny,&nz,tiffinfo);
			ntot = nx * ny * nz;	
			flagOcc = (int *) malloc(ntot * sizeof(int));	
			for (i=0; i<ntot; i++) flagOcc[i] = -1;
			zOcc = PARAMS.zOcc;
			factMarkOccXY = PARAMS.factMarkOccXY;
			PARAMS.zOcc *= 5;	//increase zOcc to avoid doubles
			PARAMS.factMarkOccXY *= 2;	//increase the radius of marking
			//shift the positions of the linked points. 
			for (i=0; i<npoints; i++) {
				pointsAll[i]->x -= X[kk];
				pointsAll[i]->y -= Y[kk];
				pointsAll[i]->z -= Z[kk];
			}
			if (npoints > 0) {
				markPixelsOccupied(nx,ny,nz,flagOcc,npoints,pointsAll);
			}
			//shift back. 
			for (i=0; i<npoints; i++) {
				pointsAll[i]->x += X[kk];
				pointsAll[i]->y += Y[kk];
				pointsAll[i]->z += Z[kk];
			}
			PARAMS.zOcc = zOcc;	//restore parameters
			PARAMS.factMarkOccXY = factMarkOccXY;
			//create linked points. 
			IDStart = maxID + 1;
			i = 0;
			while (i < ns) {
				x = linkedP[i++];
				y = linkedP[i++]; 
				z = linkedP[i++];
				r = linkedP[i++];
				ID = (int)linkedP[i++]+IDStart; 				
				if (ID > maxID) maxID = ID;
				Type = (int)linkedP[i++];
				lp = CreateLinkedPoint(x,y,z,r,ID,Type);
				numConn = (int) linkedP[i++];
				for (j=0; j<numConn; j++) {
					ID = (int)linkedP[i++]+IDStart;
					AddConnLinkedPoint(lp,ID);
				}
				//check if this points is already occupied. 
				if (flagOcc[getIndex3DZfirst((int)(x-X[kk]),(int)(y-Y[kk]),(int)(z-Z[kk]),nx,ny,nz)] == -1) {
					pointsAll[npoints++] = lp;
				} else {
					DeleteLinkedPoint(lp);
				}					
			}
			free(flagOcc);	
			//repair linked points, get rid of IDs in linked that do not appear. 
			IDflags = (int *) malloc((maxID+1) * sizeof(int));
			for (i=0; i<maxID+1; i++) IDflags[i] = 0;
			for (i=0; i<npoints; i++) {
				IDflags[pointsAll[i]->ID] = 1;
			}
			for (i=0; i<npoints; i++) {
				iter = pointsAll[i]->conn->next;
				iterpre = pointsAll[i]->conn;
				while (iter != NULL) {
					if (IDflags[iter->val] == 0) {//this id does not exist, delete
						pointsAll[i]->conn->val -= 1;
						iterpre->next = iter->next;
						free(iter);
						iter = iterpre->next;
					} else {
						iterpre = iter;
						iter = iter->next;
					}
				}
			}
			free(IDflags);
		}	
			
		//grow from end points in each tile. 
		printf("Growing from ends in each tile....\n");
		for (kk=0; kk<numTiles; kk++) {
			sprintf(filename,"%s%d.tif",filenameCommon,fileIDs[kk]);
			printf("Loading image from %s\n",filename);		 
			//read image. 
			probeFileDimensions(filename,&nx,&ny,&nz,tiffinfo);
			ntot = nx * ny * nz;	
			nxy = nx * ny;	
			im3d = (float *) malloc(ntot*sizeof(float));	
			imFlat = (float *) malloc(nxy * sizeof(float));			
			readImageAndReturnIm3d(filename,nx,ny,nz,im3d,tiffinfo); 
			
			//shift the positions of the linked points. 
			for (i=0; i<npoints; i++) {
				pointsAll[i]->x -= X[kk];
				pointsAll[i]->y -= Y[kk];
				pointsAll[i]->z -= Z[kk];
			}
			growFromEndPoints(nx, ny, nz, im3d, &npoints, pointsAll);
			//shift back.
			for (i=0; i<npoints; i++) {
				pointsAll[i]->x += X[kk];
				pointsAll[i]->y += Y[kk];
				pointsAll[i]->z += Z[kk];
			}		
			free(im3d);
			free(imFlat);						
		}	
		
		//save linked points to file. 
		//write the linked points to file. 
		ns = 0;
		for (i=0; i<npoints; i++) {
			linkedP[ns++] = pointsAll[i]->x;
			linkedP[ns++] = pointsAll[i]->y;
			linkedP[ns++] = pointsAll[i]->z;
			linkedP[ns++] = pointsAll[i]->r;
			linkedP[ns++] = pointsAll[i]->ID;
			linkedP[ns++] = pointsAll[i]->Type;
			iter = pointsAll[i]->conn;
			linkedP[ns++] = iter->val;
			iter = iter->next;
			while (iter != NULL) {
				linkedP[ns++] = iter->val;
				iter = iter->next;
			}
		}
		sprintf(outName,"%s.linkedPoints.dat",filenameCommon);
		printf("Saving linked points to %s...\n",outName);
		fpt = fopen(outName,"wb");
		fwrite(&npoints,sizeof(int),1,fpt);		
		fwrite(&ns,sizeof(int),1,fpt);
		fwrite(linkedP,sizeof(float),ns,fpt);
		fclose(fpt);
		
		//delete memopry. 
		for (i=0; i<npoints; i++) {
			DeleteLinkedPoint(pointsAll[i]);
			pointsAll[i] = NULL;
		}
		free(pointsAll);	
		free(linkedP);
	}	 
	MPI_Barrier(MPI_COMM_WORLD);
	#endif //if FULL_RECONSTRUCTION

	//read in the linked points. 
	sprintf(filename,"%s.linkedPoints.dat",filenameCommon);
	printf("Reading linked points from %s...\n",filename);
	fpt = fopen(filename,"rb");
	if (fread(&npointsTot,sizeof(int),1,fpt) == EOF) {
		printf("Reading failed.\n");
		return 1;
	}
	if (fread(&ns,sizeof(int),1,fpt) == EOF) {
		printf("Reading failed.\n");
		return 1;
	}
	printf("npointsTot=%d\n",npointsTot);	
	pointsAll = (LinkedPoint **) malloc(npointsTot * sizeof(LinkedPoint *));			
	for (i=0; i<npointsTot; i++) pointsAll[i] = NULL;
	linkedP = (float *) malloc(ns * sizeof(float));
	if (fread(linkedP,sizeof(float),ns,fpt) == EOF) {
		printf("Reading failed.\n");
		free(linkedP);
		for (i=0; i<npoints; i++ ) {
			DeleteLinkedPoint(pointsAll[i]);
		}
		free(pointsAll);
		return 1;
	}
	fclose(fpt);

	//create linked points. 
	i = 0;
	npoints = 0;
	while (i < ns) {
		x = linkedP[i++];
		y = linkedP[i++]; 
		z = linkedP[i++];
		r = linkedP[i++];
		ID = (int)linkedP[i++]; 				
		Type = (int)linkedP[i++];
		lp = CreateLinkedPoint(x,y,z,r,ID,Type);
		numConn = (int) linkedP[i++];
		for (j=0; j<numConn; j++) {
			ID = (int)linkedP[i++];
			AddConnLinkedPoint(lp,ID);
		}
		pointsAll[npoints++] = lp;
	}		

	#if CONNECT_ENDS
	//connect end points. Here is the parallelize version of the function connectEndPoints.
	mmax = 0;
	for (i=0;i<npoints;i++) {
		if (mmax < pointsAll[i]->ID) mmax = pointsAll[i]->ID;
	}
	flags = (int *) malloc((mmax+1) * sizeof(int));
	IDInds = (int *) malloc((mmax+1) * sizeof(int));
	for (i=0; i<npoints; i++) IDInds[pointsAll[i]->ID] = i;

	//get the end points. 
	eps = (int *) malloc(npoints * sizeof(int));
	epsVec = (float *) malloc(2*npoints * sizeof(float));
	//get the end points. 
	
	ne = 0;
	for (i=0; i<npoints; i++) {
		if (NumConnLinkedPoint(pointsAll[i]) == 1) {
			eps[ne] = i;
			j = IDInds[pointsAll[i]->conn->next->val];
			xs = pointsAll[i]->x - pointsAll[j]->x;
			ys = pointsAll[i]->y - pointsAll[j]->y;
			dd = fmax(1e-5,sqrt(xs*xs + ys*ys));
			epsVec[2*ne] = xs/dd;
			epsVec[2*ne+1] = ys/dd;			
			ne++;	
		}
	}		
	printf("rank=%d ne=%d\n",rank,ne);
	MPI_Barrier(MPI_COMM_WORLD);	
		
	if (ne > 0) { //connect end points. 
		if (rank == 0) {
			if (heap != NULL) free(heap); 
			heap = (struct valInd  *) malloc(ne *ne * sizeof(struct valInd));
			Init(); //heap structure for finding minimum.
		}

		//compute distances between the end points. 
		DBlock = (double *) malloc(2 * ne * (int)fmax(1,ne/numProc) * sizeof(double)); 
		//compute the distances between the end points.  

		//compute the distance matrix. 
		nD = 0;
		printf("rank=%d computing the distances between the end points...\n",rank);

		for (i=rank*ne/numProc; i < fmin(ne,(rank+1)*ne/numProc); i++) {			
			ii = eps[i];
			x1 = pointsAll[ii]->x;
			y1 = pointsAll[ii]->y;
			z1 = pointsAll[ii]->z;
			r1 = pointsAll[ii]->r;
			for (j=i+1; j<ne; j++) {			
				jj = eps[j];
				for (kk=0;kk<mmax+1;kk++) flags[kk] = 0;
				if (checkConnected(pointsAll[ii], pointsAll[jj]->ID, npoints, pointsAll, IDInds, flags) == 1) {
					continue;
				} else {	
					x2 = pointsAll[jj]->x;
					y2 = pointsAll[jj]->y;
					z2 = pointsAll[jj]->z;
					r2 = pointsAll[jj]->r;

					if (PARAMS.zfact * abs(z2 - z1) > PARAMS.zJumpFact * (r1 + r2)) continue; //z difference is too large. 		

					xv = x2 - x1;
					yv = y2 - y1;
					dd = fmax(1e-5,sqrt(xv*xv+yv*yv));
					xv /= dd;
					yv /= dd;

					if (dd <= 1.5 * (r1 + r2)) { //these two points are touching. Connect. 
						DBlock[nD++] = dd;	
						DBlock[nD++] = getIndex(i,j,ne,ne);		
						continue;			
					} 

					dpot = epsVec[2*i] * xv + epsVec[2*i+1] * yv;
					if (dpot < cos(PARAMS.angle/180.0*3.14159))  continue;  
					dpot = -epsVec[2*j] * xv - epsVec[2*j+1] * yv;
					if (dpot < cos(PARAMS.angle/180.0*3.14159)) continue;   
						
					if (dd < (r1 + r2) * PARAMS.distFactConn) { //candidate for connection, insert to the heap.  
						dd = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (PARAMS.zfact*PARAMS.zfact)*(z1-z2)*(z1-z2));	//compute the 3D distance.
						DBlock[nD++] = dd;	
						DBlock[nD++] = getIndex(i,j,ne,ne);		
					}
				}
			}				
		}
				
		if (rank == 0) {
			for (i=0; i<nD; i += 2) {
				element.value = DBlock[i];
				element.ind = (long int) DBlock[i+1];
				InsertOrUpdate(element);				
			} 
		} else {
			MPI_Send(DBlock, nD, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);	
		} 
			
		if (rank == 0) {
			for (kk=1; kk<numProc; kk++) {
				// Probe for an incoming message from process kk
				MPI_Probe(kk, 0, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_DOUBLE, &nD);			
				MPI_Recv(DBlock,nD, MPI_DOUBLE, kk, 0, MPI_COMM_WORLD, &status);
				printf("Processing data from rank=%d recieved nD=%d\n",kk,nD);
				for (i=0; i<nD; i += 2) {
					element.value = DBlock[i];
					element.ind = (long int) DBlock[i+1];
					InsertOrUpdate(element);				
				} 
			}

			//connect ends. 
			visited = (int *) malloc(ne * sizeof(int));
			for (i=0; i<ne; i++) visited[i] = 0;
			while (1) {
				if (heapSize == 0) break;
				element = DeleteMin();
				i = element.ind/ne;
				j = element.ind - i * ne;
				ii = eps[i];
				jj = eps[j];

				z1 = pointsAll[ii]->z;
				r1 = pointsAll[ii]->r;
				z2 = pointsAll[jj]->z;
				r2 = pointsAll[jj]->r;
				dd = element.value;

				//for other end points restrict one connection only. 			
				if (visited[i] == 1 || visited[j] == 1) continue;				
				visited[i] = 1; visited[j] = 1;	
				AddConnLinkedPoint(pointsAll[ii],pointsAll[jj]->ID);
				AddConnLinkedPoint(pointsAll[jj],pointsAll[ii]->ID);
			}			
			free(visited);
			free(heap);	heap = NULL;
				
			//create swc file. 			
			printf("Creating swc file...\n");
			sprintf(outName,"%s.auto.swc",filenameCommon);
			createSWC(&npoints, pointsAll, outName);						
		}			
	}	

	//free memeory
	if (IDInds != NULL) free(IDInds);
	if (flags  != NULL) free(flags);
	if (eps  != NULL) free(eps);
	if (epsVec  != NULL) free(epsVec);	
	if (ne > 0) {
		free(DBlock);
	}
	#else
	if (rank == 0) {		//create swc file. 			
		printf("Creating swc file...\n");
		sprintf(outName,"%s.auto.swc",filenameCommon);
		createSWC(&npoints, pointsAll, outName);						
	}
	#endif //	#if CONNECT_ENDS
	 
	for (i=0; i<npoints; i++ ) {
		DeleteLinkedPoint(pointsAll[i]);
	}
	if (pointsAll  != NULL) free(pointsAll);
	if (X != NULL) free(X);
	if (Y != NULL) free(Y);
	if (Z != NULL) free(Z);
	if (fileIDs  != NULL) free(fileIDs);
	if (linkedP  != NULL) free(linkedP);
	MPI_Finalize();
	
	return 0;
}

 
 
