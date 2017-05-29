# This is python code for preprocessing microscopic images, including creating tif stacks, sticthing, and autotrace. 
# The results will be loaded in ShuTu for semi-automatic neuron reconstruction. 
#
# Copyright (C) 2014, Dezhe Z. Jin (dzj2@psu.edu), Ting Zhao

#	This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#------------------Python libraries------------------------------------------------------------------------------

# linux or mac, set PYTHONPATH to the directopry where this file is (thePathToShuTu.py). 
# in ~/.bashrc, add line export PYTHONPATH=$PYTHONPATH:/thePathToShuTu.py/

import os, os.path, glob, re, fnmatch, traceback, numpy, ctypes, PIL, json
from multiprocessing import Pool
from scipy.ndimage import gaussian_filter
from scipy.ndimage.filters import laplace
from pylab import *
import xml.etree.ElementTree as ET
from subprocess import call
from collections import OrderedDict
from llist import dllist	# NOTE: if fails to load, need to install llist, see https://pypi.python.org/pypi/llist/

#------------------ Parameters ----------------------------------------------------------------------------------
#
# Users need to update the following parameters based on their settings. 
#

### The directory of the library libShuTu.so ####

place='NJ'	# O, office, NJ, SC, MC, mac, J, Janelia

if place == 'MC':
	dirCodes = '/Users/dezhejin/projects/ShuTu/'		# mac directory for the c library. 
elif place == 'J':
	dirCodes='/groups/spruston/home/jind/NeuTuRelease/' # linux directory for the c library.
elif place == 'NJ':
	dirCodes = '/home/djin/ShuTu/PythonCodes/'			# ubuntu directory for the c library. 
elif place == 'O' or place == 'SC':
	dirCodes = '/home/dezhe/ShuTu/PythonCodes/'			# ubuntu directory for the c library. 

print ' ' 
print 'Check where these directories are correct:'
print 'dirCodes (libShuTu.so) = ',dirCodes
print ' ' 

### Microscopy settings ######

mag = 100				# microscope magnification
#mag = 63				# microscope magnification
microscopeID = 1		# together with mag, this is used to set xyDist (distance between pixels) and zDist (distance between planes) in micron 
						# 2, Janelia Zeiss
						# 3, Janelia confocal

if microscopeID == 1:
	microscopeName = 'Janelia Zeiss Brightfield Microscope'
	if mag == 100:
		xyDist = 0.065			# micron, pixel distance in x y plane
		zDist = 0.5				# micron, pixel distance in z direction
	elif mag == 63:
		xyDist = 0.102			# micron, pixel distance in x y plane
		zDist = 0.5				# micron, pixel distance in z direction
	else:
		print 'ERROR: Unkown magnification for ',microscopeName
	
print 'Microscope :', microscopeName 
print 'Magnification = ', mag 
print 'Pixel distance (xy, micron):', xyDist
print 'Pixel distance (z, micron):', zDist
print ' '


### Preprocessing paramaters #### 

sigmaBack = mag*100.0	# in pixels, sigma for substracting the background. 
nProc = 2				# number of processors used for parallel computing in images processing 
						# make sure that the tif stack size * nProc does not exceed the internal memory. Otherwise the created tif will be corrupt. 
						
#----------------Set working environment--------------------------------------------------------
#
# Go to the directory where the image files are. 
# Then run setDataSetCurrent()
# This will define the following variables. 
# If the variables are not set correctly, modify setDataSetCurrent() or specify these variables in other ways. 

filenameCommonSlice = '' 	# common string in the tif images. Format: filenameCommon + numbers + .tif
nTiles = 0;					# number of tiles	
filenameCommon = ''			# common string in the tif stacks to be created	
overlapList = ''			# list of pairs of overlapping tiles for stitching. The format is [(1,2),(2,3)] 

def setDataSetCurrent():	# automatically set filenameCommon and load overlapList from file if already exists. 
	global filenameCommon, overlapList
	# define filenameCommon to current directory
	files = glob.glob("*1.tif")
	if len(files) == 0:
		print "There are not tif files."
	else:	
		filenameCommon = re.split('[0-9]+.tif',files[0])[0]
	if len(filenameCommon) == 0 and len(files)>0:
		filenameCommon = raw_input('Set the common prefix for the tif stacks. Input: ')	
	if len(overlapList) == 0:
		try:
			datFile = filenameCommon+'StitchCoord.npz'
			print 'Loading overlapList from ',datFile
			res = numpy.load(datFile)
			overlapList = res['overlapList']
		except:
			print 'overlapList is not loaded from file.'	
	print 'overlapList = ',overlapList			
	print 'filenameCommon is set to ',filenameCommon

setDataSetCurrent()
	
#---------------Create tiff stack using single plane images from the microscpe-------------------------
#
# These functions are used to create tiff stacks from individual tiff images. 
# These processes depend on the setup of imaging process and the output format. 
# They are not needed if the tif stacks are already generated by the microscope's software. 

###
# BUGS: 
# createTifStacksFromSlices hangs when the memory demand in parallel processing exceeds 100%. 
### 

####
# create tiff stacks from plane images in the zeiss microscope at Janelia
# assume filename convention is like slides-09_z035m30.tif, i.e. commonname_z(depth)m(tilenumber).tif
####

def getFilenameCommonSlice():	# tries to find the common filename of the slices.
	global filenameCommonSlice
	ending = '_ORG.tif'			# The ending of the slices are assumed to be _ORG.tif
	print 'Finding filenameCommonSlice by searching files ending with ',ending
	fname = ''
	for root, dirnames, filenames in os.walk('.'):
		for filename in fnmatch.filter(filenames, '*'+ending):	
			fname= os.path.join(root, filename)
			break
		if len(fname) > 0:
			break
	print fname
	iid = fname.find('_z')
	filenameCommonSlice = fname[0:iid]
	print 'Setting filenameCommonSlice to ', filenameCommonSlice
	return filenameCommonSlice

def findNumberOfTiles(filenameCommonSlice):	
	files = glob.glob(filenameCommonSlice+'_*'+"m1_ORG.tif")
	if len(files) == 0:
		files = glob.glob(filenameCommonSlice+'_*'+"m01_ORG.tif")		
		if len(files) == 0:
			files = glob.glob(filenameCommonSlice+'_*'+"m001_ORG.tif")
			if len(files) == 0:
				print 'No slice images found with filenameCommonSlice=',filenameCommonSlice	
				return 0
	NS = len(files)
	files = glob.glob(filenameCommonSlice+'_*_ORG.tif')
	NSAll = len(files)
	nTiles = NSAll/NS
	print 'There are ',nTiles,' tiles, each has ',NS,' slices'
	return nTiles

def createTifStacksFromSlices():
	nTiles = findNumberOfTiles(filenameCommonSlice)	
	if nTiles == 0:
		getFilenameCommonSlice()	# try to automatically get filenameCommonSlice
		nTiles = findNumberOfTiles(filenameCommonSlice)			
	print 'Number of processors used ',nProc
	stacks = range(1,nTiles+1)
	if len(stacks) > 0:		# COMMENT: pool hangs if nProc does not commensurate len(stacks). Need to find a better coding. 
		pool = Pool(processes = min(nProc,nTiles))
		pool.map(createOneTifStackFromSlices,stacks,chunksize=1)
		pool.close()
		pool.join()
	else:
		print 'ERROR in createTifStacksFromSlices, nTiles=0.'
		print 'The common string of slice filenames is set to ',filenameCommonSlice,'. Check if this is correct.'	
	# check the created tif stacks. Some may be corrupt due to too large nProc. Recreate these files. 
	print 'Checking corrupted tif stacks and try to recreate...'
	files = glob.glob(filenameCommon+'[0-9]*.tif')
	filesizes = [os.stat(fn).st_size for fn in files]
	mmax = max(filesizes)
	for ii in range(len(files)):
		if filesizes[ii] < mmax:
			os.remove(files[ii])
			fn = int(files[ii][len(filenameCommon):len(files[ii])-4])
			createOneTifStackFromSlices(fn)
	
def createOneTifStackFromSlices(istack):
	try:
		prefix = filenameCommonSlice+'_z'
		suffix = 'm%02d_ORG.tif' % (istack,)	
		numWidth = 3
		files = sort(glob.glob(prefix+'*'+suffix))
		if len(files) == 0:		# the number of tiles can be larger than 100. 
			suffix = 'm%03d_ORG.tif' % (istack,)	# NOTE if the number of tiles is more than 100, change to 3d. 
			files = sort(glob.glob(prefix+'*'+suffix))
			if len(files) == 0:
				suffix = 'm%01d_ORG.tif' % (istack,)	# NOTE if the number of tiles is smaller than 100, change to 1d. 
				files = sort(glob.glob(prefix+'*'+suffix))
				if len(files) == 0:
					print 'ERROR IN createOneTifStackFromSlices: Cannot locate the tiff slices.'
					return 0					
		filenameStack = filenameCommon+'-'+str(istack)+'.tif'
		# see if the file exists.
		fs = glob.glob(filenameStack)
		if len(fs) > 0:
			print filenameStack,'exists. Skip creating.'
			return 0
		print 'Creating tiff stack ',filenameStack
		# load the C program
		sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
		# set parameter types
		sfmLib.createTiffStackFromSlices.argtypes = [ctypes.c_long,ctypes.POINTER(ctypes.c_char_p),ctypes.c_char_p]
		sfmLib.createTiffStackFromSlices.restype = None
		# call the C function.
		fpointers = (ctypes.c_char_p*len(files))()
		for i in range(len(files)):
			fpointers[i] = files[i]
		sfmLib.createTiffStackFromSlices(ctypes.c_long(len(files)),fpointers,ctypes.c_char_p(filenameStack))
		return len(files)
	except Exception as e:
		traceback.print_exc()
		print()
		raise e
		return 0	

#-----------Processing image stacks for stitching and images used in ShuTu --------------------------------------

# process all images. 
# if the images are dark field, speciy imageType = 1 by calling processAllImages(imageType=1). 
# The dark field images will be convergted into brightfield background. The original tifs will be moved to a directory called originalTiffs. 
def processAllImages(imageType=0,iforce=0):			# if iforce = 0, only process tif stacks yet not processed. if 1, process all. 
	if imageType == 0:
		print "Bright field image."
	elif imageType == 1:
		print "Dark field image."
	else:
		print ("ERROR: specify imageType = 0 or 1.");
	files = glob.glob(filenameCommon+'*[0-9].tif')
	filenames = []
	for fn in files:
		fnsp = fn.split('.')
		if len(fnsp) > 2:	# do not confuse with some of the Flatten.tif file.
			continue
		filenameBase = fnsp[0]
		if iforce==0:
			# see if this file has been processed. Check .Proj.tif file. 
			if os.path.exists(filenameBase+'.Proj.tif'):
				continue
			else:
				filenames.append((filenameBase,imageType))
		else:
			filenames.append((filenameBase,imageType))					
	if len(filenames) == 0:
		print 'IN processAllImages: nothing to process. ' 
		return			
	if len(filenames) > 0:
		pool = Pool(processes = min(nProc,len(filenames)))
		pool.map(processImage,filenames,chunksize = 1)
		pool.close()
		pool.join()

def processImage(params):
	filenameBase,imageType = params
	if imageType == 1 and os.path.isfile(filenameBase+".org.tif"):
		print "The image has been converted to bright field already. If this not the case, delete *.org.tif files."
		imageType = 0
	try:
		# load the C program
		sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
		sfmLib.processImage.argtypes = [ctypes.c_char_p,ctypes.c_long]
		sfmLib.processImage.restype = None
		sfmLib.processImage(ctypes.c_char_p(filenameBase),ctypes.c_long(imageType))		
	except:
		print 'ERROR in processImage: failed.'
		return	
	
def saveImageData(filenameBase,im3d,imFlat,imFlatRGB):
	nx,ny,nz = im3d.shape
	print 'Saving image data to ',filenameBase+'.dat...'
	numpy.savez(filenameBase+'.dat',nx=nx,ny=ny,nz=nz)
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	# set parameter types
	imFlatRGB  = imFlatRGB.astype(float32)
	sfmLib.saveImageData.argtypes = [ctypes.c_char_p,ctypes.c_long,ctypes.c_long,ctypes.c_long, \
		ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float)]
	sfmLib.saveImageData.restype = None
	# call the C function.
	sfmLib.saveImageData(ctypes.c_char_p(filenameBase+'.dat'),  \
			ctypes.c_long(nx), ctypes.c_long(ny), ctypes.c_long(nz), \
			im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
			imFlat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
			imFlatRGB.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))

def readImageData(filenameBase):
	print 'Reading image data from ',filenameBase+'.dat...'
	nx,ny,nz = checkImageSize(filenameBase+'1.tif')
	ntot = nx * ny * nz
	im3d = zeros(ntot).astype(float32)
	imFlat = zeros((nx,ny),dtype='float32')
	imFlatRGB = zeros(nx*ny*3).astype('float32')
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	# set parameter types
	sfmLib.readImageData.argtypes = [ctypes.c_char_p,ctypes.c_long,ctypes.c_long,ctypes.c_long, ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float)]
	sfmLib.readImageData.restype = None
	# call the C function.
	sfmLib.readImageData(ctypes.c_char_p(filenameBase+'.dat'),  \
			ctypes.c_long(nx), ctypes.c_long(ny), ctypes.c_long(nz), \
			im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
			imFlat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
			imFlatRGB.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))			
	im3d = im3d.reshape(nz,nx,ny)
	im3d = swapaxes(im3d,0,1)
	im3d = swapaxes(im3d,1,2)
	imFlatRGB = imFlatRGB.reshape(3,ny,nx)
	imFlatRGB = swapaxes(imFlatRGB,0,1)
	imFlatRGB = swapaxes(imFlatRGB,1,2)
	imFlatRGB = swapaxes(imFlatRGB,0,1)	
	return im3d,imFlat,imFlatRGB

# This file creates a tif stack from im3d and saves to filename.
def createTiffStackFromArray(im3d,filenameBase,iReorder=1):
	filename = filenameBase + '.tif'
	print 'Creating a tiff stack from array to ',filename
	im3d = im3d.astype(float32)
	nx,ny,nz = im3d.shape
	if iReorder == 1:
		# Need to make sure that the memory of the array is in (nz, nx, ny) order. 
		im3d2 = zeros((nz,nx,ny)).astype(float32)
		for k in range(nz):
			im3d2[k,:,:] = im3d[:,:,k]
		im3d = im3d2
	else:
		print 'The order of the memory order is assumed to be nz, nx, ny.'	
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	# set parameter types
	sfmLib.createGrayTiffStackFromFloatArray.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_float),ctypes.c_char_p]
	sfmLib.createGrayTiffStackFromFloatArray.restype = None
	# call the C function.
	sfmLib.createGrayTiffStackFromFloatArray(ctypes.c_int(nx),ctypes.c_int(ny),ctypes.c_int(nz),im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.c_char_p(filename))	

def checkImageSize(tiffFilename):
	# query the size of tiff
	tiffinfo = array([0,0,0,0,0]).astype(int32)
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	# set parameter types
	sfmLib.getTiffAttribute.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int)]
	sfmLib.getTiffAttribute.restype = None
	# call the C function.
	sfmLib.getTiffAttribute(ctypes.c_char_p(tiffFilename),tiffinfo.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))	
	nx,ny,nz,nchan,bitType = tiffinfo
	if bitType == -1:
		print 'Inquariy of file attributions failed. Failed to load ',tiffFilename
		return 0,0,0
	print "Image size (",nx,",",ny,",",nz,")"
	return nx,ny,nz
	
def loadStacks(tiffFilename,imageType=0):	
	print 'Loading stack image from ',tiffFilename
	fn = glob.glob(tiffFilename)
	if len(fn) == 0:
		print 'File ',tiffFilename,'does not exist.'
		return 0,0,0
	try:
		# query the size of tiff
		tiffinfo = array([0,0,0,0,0]).astype(int32)
		# load the C program
		sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
		# set parameter types
		sfmLib.getTiffAttribute.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int)]
		sfmLib.getTiffAttribute.restype = None
		# call the C function.
		sfmLib.getTiffAttribute(ctypes.c_char_p(tiffFilename),tiffinfo.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))	
		nx,ny,nz,nchan,bitType = tiffinfo
		if bitType == -1:
			print 'Inquariy of file attributions failed. Failed to load ',tiffFilename
			return 0,0,0
			
		# now load the tiff file. 
		ntot = nx * ny * nz
		im3dR = zeros(ntot).astype(float32)
		im3dG = zeros(ntot).astype(float32)
		im3dB = zeros(ntot).astype(float32)

		# set parameter types
		sfmLib.readTiffStack.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int),\
						ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float)]
		sfmLib.readTiffStack.restype = None
		# call the C function.
		sfmLib.readTiffStack(ctypes.c_char_p(tiffFilename), \
						tiffinfo.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), \
						im3dR.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
						im3dG.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), im3dB.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))	

		im3dR = im3dR.reshape(nz,nx,ny)
		im3dR = swapaxes(im3dR,0,1)
		im3dR = swapaxes(im3dR,1,2)
		if im3dG.max() == 0 and im3dB.max() == 0:
			nchan = 1	# This is psedu color; there is only 1 channel.
		if nchan > 1:
			im3dG = im3dG.reshape(nz,nx,ny)
			im3dG = swapaxes(im3dG,0,1)
			im3dG = swapaxes(im3dG,1,2)
			im3dB = im3dB.reshape(nz,nx,ny)
			im3dB = swapaxes(im3dB,0,1)
			im3dB = swapaxes(im3dB,1,2)
		
		# smooth the planes.
		sigmaSmooth = 1				# PARAMETER
		for i in range(nz):
			im3dR[:,:,i] = gaussian_filter(im3dR[:,:,i],sigma=sigmaSmooth)
			if nchan > 1:
				im3dG[:,:,i] = gaussian_filter(im3dG[:,:,i],sigma=sigmaSmooth)
				im3dB[:,:,i] = gaussian_filter(im3dB[:,:,i],sigma=sigmaSmooth)

		if bitType == 0:	# 8 bit
			maxVal = 255.0
		elif bitType == 1:  # 16 bit
			maxVal = 65535.0
		else:
			maxVal = im3dR.max()
			if nchan > 1:
				maxVal = max(maxVal,im3dG.max())
				maxVal = max(maxVal,im3dB.max())

		# judge whether the background is bright or dark.
		if imageType != 0:
			print 'Dark background image. Inverting images...'
			im3dR = maxVal - im3dR
			if nchan > 1:
				im3dG = maxVal - im3dG
				im3dB = maxVal - im3dB
		else:
			print 'Bright background image.'

		print 'Projecting red...'
		# nonlinear projection emphasizing brown color (158,85,54).
		#tR = 158.0/255; tG = 85.0/255.0; tB = 54.0/255.0
		#im3dGray = tanh(0.2126 *  abs(im3dR/maxVal - tR) + 0.7152  * abs(im3dG/maxVal - tG) + 0.0722 *  abs(im3dB/maxVal - tB))
		if nchan > 1:
			im3dGray = 0.2126 * im3dR +  0.7152 * im3dG + 0.0722 * im3dB	
		else:
			im3dGray = im3dR
			
		# minimum intensity projection.
		imFlatR = im3dR.min(axis = 2);
		imFlatR = substractBackground(imFlatR,sigmaBackground = 100)	# parameter
		mmin = percentile(imFlatR,0.1)
		mmax = percentile(imFlatR,99.9)
		#mmin = imFlatR.min()
		#mmax = imFlatR.max()
		imFlatR = 255 * (imFlatR - mmin)/(mmax - mmin)
		imFlatR[where(imFlatR <0)] = 0
		imFlatR[where(imFlatR > 255)] = 255	

		if nchan > 1:
			print 'Projecting green and blue...'
			imFlatG = im3dG.min(axis = 2);
			imFlatG = substractBackground(imFlatG,sigmaBackground = 100)
			mmin = percentile(imFlatG,0.1)
			mmax = percentile(imFlatG,99.9)
			#mmin = imFlatG.min()
			#mmax = imFlatG.max()
			imFlatG = 255 * (imFlatG - mmin)/(mmax - mmin)
			imFlatG[where(imFlatG <0)] = 0
			imFlatG[where(imFlatG > 255)] = 255	
					
			imFlatB = im3dB.min(axis = 2);
			imFlatB = substractBackground(imFlatB,sigmaBackground = 100)
			mmin = percentile(imFlatB,0.1)
			mmax = percentile(imFlatB,99.9)
			#mmin = imFlatB.min()
			#mmax = imFlatB.max()
			imFlatB = 255 * (imFlatB - mmin)/(mmax - mmin)
			imFlatB[where(imFlatB <0)] = 0
			imFlatB[where(imFlatB > 255)] = 255	
			
			imFlatRGB = imFlatR.astype('uint8')[...,newaxis]
			imFlatRGB = append(imFlatRGB,imFlatG.astype('uint8')[...,newaxis],2)	
			imFlatRGB = append(imFlatRGB,imFlatB.astype('uint8')[...,newaxis],2)	
			imFlatGray = 0.2126 * imFlatR +  0.7152 * imFlatG  + 0.0722 * imFlatB	
		else:	
			imFlatRGB = imFlatR.astype('uint8')[...,newaxis]
			imFlatRGB = append(imFlatRGB,imFlatR.astype('uint8')[...,newaxis],2)	
			imFlatRGB = append(imFlatRGB,imFlatR.astype('uint8')[...,newaxis],2)
			imFlatGray = imFlatR	
		#imFlatGray = amin(im3dGray,axis=2)
		mmin = imFlatGray.min()
		mmax = imFlatGray.max()
		imFlatGray = (imFlatGray - mmin)/(mmax - mmin)	
		# normalize
		mmin = im3dGray.min()
		mmax = im3dGray.max()
		im3dGray = (im3dGray - mmin)/(mmax - mmin)
		return im3dGray,imFlatGray,imFlatRGB
	except: 
		print 'Failed to load ',tiffFilename
		return 0,0,0

def substractBackground(imFlat,sigmaBackground):
	img = imFlat.copy()
	img = gaussian_filter(img,sigma=sigmaBackground) - img
	#img[where(img < 0)] = 0
	#img = img ** 0.5
	img = 1 - img/img.max()
	img = img ** 2
	return img
	
def convertDarkbackgroundToBrightBackground(tiffFilename):
	# this function read in a darkbackground tif stack and converts to bright background. 
	# the saved tif stack is in 8 bits as required by ShuTu. 
	print 'Loading stack image from ',tiffFilename
	fn = glob.glob(tiffFilename)
	if len(fn) == 0:
		print 'File ',tiffFilename,'does not exist.'
		return
	try:
		# query the size of tiff
		tiffinfo = array([0,0,0,0,0]).astype(int32)
		# load the C program
		sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
		# set parameter types
		sfmLib.getTiffAttribute.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int)]
		sfmLib.getTiffAttribute.restype = None
		# call the C function.
		sfmLib.getTiffAttribute(ctypes.c_char_p(tiffFilename),tiffinfo.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))	
		nx,ny,nz,nchan,bitType = tiffinfo
		if bitType == -1:
			print 'Failed to load ',tiffFilename
			return
		print 'Converting to 8 bit bright field tif stack...'	
		# now load the tiff file. 
		ntot = nx * ny * nz
		im3dR = zeros(ntot).astype(float32)
		im3dG = zeros(ntot).astype(float32)
		im3dB = zeros(ntot).astype(float32)

		# set parameter types
		sfmLib.readTiffStack.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_int),\
						ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float)]
		sfmLib.readTiffStack.restype = None
		# call the C function.
		sfmLib.readTiffStack(ctypes.c_char_p(tiffFilename), \
						tiffinfo.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), \
						im3dR.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
						im3dG.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), im3dB.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))	

		im3dR = im3dR.reshape(nz,nx,ny)
		im3dR = swapaxes(im3dR,0,1)
		im3dR = swapaxes(im3dR,1,2)
		if im3dG.max() == 0 and im3dB.max() == 0:
			nchan = 1	# This is psedu color; there is only 1 channel.
		if nchan > 1:
			im3dG = im3dG.reshape(nz,nx,ny)
			im3dG = swapaxes(im3dG,0,1)
			im3dG = swapaxes(im3dG,1,2)
			im3dB = im3dB.reshape(nz,nx,ny)
			im3dB = swapaxes(im3dB,0,1)
			im3dB = swapaxes(im3dB,1,2)
		print 'Dark background image. Inverting images...'
		if bitType == 0:	# 8 bit
			maxVal = 255.0
		elif bitType == 1:  # 16 bit
			maxVal = 65535.0
		else:
			maxVal = im3dR.max()
			if nchan > 1:
				maxVal = max(maxVal,im3dG.max())
				maxVal = max(maxVal,im3dB.max())
		im3dR = (maxVal - im3dR)/maxVal*255.0
		if nchan > 1:
			im3dG = (maxVal - im3dG)/maxVal*255.0
			im3dB = (maxVal - im3dB)/maxVal*255.0
			im3dGray = 0.2126 * im3dR +  0.7152 * im3dG + 0.0722 * im3dB	
		else:
			im3dGray = im3dR
		
		# save the stack, mono color only. 
		filenameSave = tiffFilename+'.inv'
		createTiffStackFromArray(im3dGray,filenameSave,iReorder=0)	# do not need reorder the memory. 
	except:
		print 'Failed to convert.'

#----------Stitch Images-----------------------------------------------------
#
# There are several options for stiching the image tiles. 
# (1) xmlStitchImages(filenameCommon=filenameCommon)
#		Stitch using the slides-*_info.xml file for positions of the tiles. 
#		The xml file is created by the Zeiss microcope acquisition program ZenBlue
#
# (2) gridStitchImages(nrow,ncol,offsetPerc=80,iplot=0):
# 	Stitch images assuming that the tiles 1,2,3,...,nTiles are arranged in a nrow x ncol grid
# 	such that 
# 		1,2,...,ncol
#		ncol+1,ncol+2,...,2*ncol
#		...
#		..., nrow * ncol
#
# (3) tileSequencesStitch(filename = "tileSequences.txt"):
# 	This function stitches images by reading tile positions from file tileSequences.txt
# 	The format of the file is as follows. Comments starts with #. 
# 	Bright field images are assumed. Use processAllImages(imageType=1) to convert dark field images into bright field images. 
# 	Fisrt line is offset percent.
# 	The subsequence lines are alternations of two types. 
#	first, sequences of adjacent images shifted to the right
# 	then a pair of images shifted down. 
# 	specify sequences of adjacent images in with tile id's (eg filenameCommon20.tif, 20 is the id)
# 	first to the right, then down, then the right, then down, etc.
# 	The sequences should be such that the tile coordinates can be determined uniquely.  
# 	example
# 	1, 2, 3, 4 (2 is right to 1, 3 is right to 2, ect)
# 	1, 5 (5 is down to 1)
# 	5, 6, 7, 8 (6 is right to 5, 7 is right to 6, etc)
# 	6, 9 (9 is down to 6)
# 	9, 10, 11, (10 is right to 9, 11 is right to 10, etc)
#	10,12 (12 is down to 10)
#	12	(the file must end with a right scan, even if there is only one tile). 	
# 
# Several functions are supplied to fix stitching problems. 
#
# (1) replaceElelementInOverlapList(overlapList,(num1,num2),(num1R,num2R)):
# 	this function replaces (num1,num2) in overlpaList to (num1R,num2R)
#
# (2) ManualPairwiseSticting(filenum1,filenum2,offsetDirection,offsetPerc=80):
# 	offsetDirection for filenum1 relative to filenum2, 
# 	offsetDirection, 'up', 'down', 'left', 'right'.	
# 	offsetPerc - percent of movement in that direction. 

def xmlStitchImages(filenameCommon=filenameCommon):
	# This function parses the positions of the tiles from slides-*_info.xml file and produce sticthes of the images. 
	# The xml file is created by the zeiwess microcope aquisition program
	files = glob.glob('*_info.xml')
	if len(files) == 0:
		print "There is no *_info.xml file. Nothing to do."
		return
	print "Parsing the xml file ",files[0] 
	tree = ET.parse(files[0])
	root = tree.getroot()
	tileIDs = []
	tileX = []
	tileY = []
	for child in root:
		filename = child[0].text
		m = re.search('z0+1m[0-9]+_ORG',filename)
		if m:
			nums = re.findall('\d+',m.group(0));
			tilenum = int(nums[-1])			
			x = int(child[1].attrib['StartX'])
			y = int(child[1].attrib['StartY'])	
			nx= int(child[1].attrib['SizeX'])
			ny= int(child[1].attrib['SizeY'])	
			tileIDs.append(tilenum)
			tileX.append(x)
			tileY.append(y)
			print '(tileNum, X,Y,nx,ny):', tilenum, x, y, nx, ny
	return stitchShiftTiles(tileIDs,tileX,tileY,nx,ny)			
	
def gridStitchImages(nrow,ncol,offsetPerc=80):
	# Stitch images assuming that the tiles 1,2,3,...,nTiles are arranged in a nrow x ncol grid
	# such that 
	# 	1,2,...,ncol
	#	ncol+1,ncol+2,...,2*ncol
	#	...
	#	..., nrow * ncol
	nx,ny,nz = checkImageSize(filenameCommon+'1.tif')	
	tileIDs = []
	tileX = []
	tileY = []
	tid = 0
	for i in range(nrow):
		y = i * ny * offsetPerc/100.0 
		for j in range(ncol):
			tid += 1
			x = j * nx * offsetPerc/100.0
			tileIDs.append(tid)
			tileX.append(x)
			tileY.append(y)
	return stitchShiftTiles(tileIDs,tileX,tileY,nx,ny)	

# This function stitches images by reading tile positions from file tileSequences.txt
# The format of the file is as follows. Comments starts with #. 
# Bright field images are assumed. Use processAllImages(imageType=1) to convert dark field images into bright field images. 
# Fisrt line is offset percent.
# The subsequence lines are alternations of two types. 
# first, sequences of adjacent images shifted to the right
# then a pair of images shifted down. 
# specify sequences of adjacent images in with tile id's (eg filenameCommon20.tif, 20 is the id)
# first to the right, then down, then the right, then down, etc.
# The sequences should be such that the tile coordinates can be determined uniquely.  
# example
# 1, 2, 3, 4 (2 is right to 1, 3 is right to 2, ect)
# 1, 5 (5 is down to 1)
# 5, 6, 7, 8 (6 is right to 5, 7 is right to 6, etc)
# 6, 9 (9 is down to 6)
# 9, 10, 11, (10 is right to 9, 11 is right to 10, etc)
def tileSequencesStitch(filename = "tileSequences.txt",iplot=0):
	k = 0	
	scans = []
	downs = []
	with open(filename) as fn:
		for line in fn:
			if len(line) == 0:
				continue
			if line[0] == '#' or line[0] == '\n':
				continue
			if k == 0:				
				offsetPerc = float(line)
				print 'Offset percet = ', offsetPerc
			elif mod(k,2) == 1:	# sequence to the right
				ids = line.split(',')
				iids = []
				print '(right) ',
				for iid in ids:
					iid = int(iid)
					iids.append(iid)
					print iid,
				scans.append(iids)	
			else:	# sequence down
				ids = line.split(',')
				print '(down) ',
				iids = []
				for iid in ids:
					iid = int(iid)
					iids.append(iid)
					print iid,
				downs.append(iids)		
			k += 1
	
	# probe image dimensions, assume *.Proj.tif have been generated by processing images. 
	iid = scans[0][0]
	nx,ny,nz = checkImageSize(filenameCommon+str(iid)+'.Proj.tif')
			
	tileIDs = []
	tileX = []
	tileY = []

	sy0 = []
	sid0 = []
	x = 0
	offsetPerc /= 100.0	# convert to fraction.
	for i in range(len(scans)):
		if len(sy0) > 0:			# record the tile positions and IDs
			tileIDs += sid0
			tileY += sy0
			tileX += [x]*len(sy0)		
		y = -ny * offsetPerc		# make offset fo y to start with 0. 
		sy = []
		sid = []
		for iid in scans[i]:
			sid.append(iid)
			y += ny * offsetPerc
			sy.append(y)
		if len(sy0) == 0:
			sy0 = sy
			sid0 = sid
			continue	
		else:
			k = sid0.index(downs[i-1][0])
			y0 = sy0[k]
			k = sid.index(downs[i-1][1])
			y = sy[k]
			sy0 = [yy - (y - y0) for yy in sy]	# shift y with the information in downs. 
			x += nx * offsetPerc
			sid0 = sid
	if len(sy0) > 0:			# record the tile positions and IDs
		tileIDs += sid0
		tileY += sy0
		tileX += [x]*len(sy0)		
	# visualize the tile positions. 
	# construct a stitched 2D image
	X = array(tileX)
	Y = array(tileY)
	X -= X.min()
	Y -= Y.min()
	imFlatCombined = 255 * ones((int(X.max())+nx,int(Y.max())+ny,3))
	for kk in range(len(tileIDs)-1,-1,-1):
		ii = tileIDs[kk]
		imFlat = array(PIL.Image.open(filenameCommon+str(ii)+'.Proj.tif')).astype('uint8')
		if len(imFlat.shape) == 3:
			imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,:] = imFlat
		else:
			for jj in range(3):
				imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,jj] = imFlat					
	imFlatCombined = imFlatCombined.astype('uint8')	
	if iplot == 1:
		figure(100); clf()
		imshow(imFlatCombined)
		title("Positions of tiles")
		# print the tile IDs. 
		for kk in range(len(tileIDs)-1,-1,-1):
			ii = tileIDs[kk]
			text(Y[kk]+ny/2,X[kk]+nx/2,str(ii))
		show(block=False)	
		raw_input("Press any key to proceed. If this is not good, modify tileSequences.txt and try again.")
	return stitchShiftTiles(tileIDs,tileY,tileX,ny,nx)		# note the switch of X, Y, need to figure out why :)												

def stitchShiftTiles(tileIDs,tileX,tileY,nx,ny):
	# stitch the images using the minimum spanning tree algorithm.
	# nearby tiles are assumed to be shifted to one of the directions: left, right, up, down.  
	tileX = array(tileX)
	tileY = array(tileY)
	overlapList = []
	nt = len(tileIDs)
	D = ones((nt,nt)) * 2.0 	# Distance matrix
	threshold = 1.5  			#threshold beyond which the connection is considered broken in the minimum spanning tree algorithm. 
	# Build up the distance matrix by trying to stitch to nearby tiles. 
	ParamsStitch = []
	IJPairs = []
	for i in range(nt):
		x1 = tileX[i]
		y1 = tileY[i]
		ii = tileIDs[i]
		for j in range(i+1,nt):
			x2 = tileX[j]
			y2 = tileY[j]
			jj = tileIDs[j]
			if abs(x2 - x1) < nx/2 and y2 > y1 and y2 < y1 + ny:	# down
				direct = 'down'
				offsetPerc = (y2 - y1) * 100.0/ny
				flag = 1
			elif abs(x2 - x1) < nx/2 and y2 < y1 and y2 > y1 - ny:	# up
				direct = 'up'
				offsetPerc = (y1 - y2) * 100.0/ny
				flag = 1
			elif abs(y2 - y1) < ny/2 and x2 < x1 and x2 > x1 - nx:	# left
				direct = 'left'
				offsetPerc = (x1 - x2) * 100.0/nx
				flag = 1
			elif abs(y2 - y1) < ny/2 and x2 > x1 and x2 < x1 + nx:	# right
				direct = 'right'
				offsetPerc = (x2 - x1) * 100.0/nx
				flag = 1
			else:
				flag = 0	
			if flag == 1:
				# add safety margin to the offsetPerc, add 10 more percent. 
				offsetPerc = max(offsetPerc - 10, 0)
				ParamsStitch.append((ii,jj,offsetPerc,direct))
				IJPairs.append((i,j))
				
	# compute shifts with multiple processures. 
	if len(ParamsStitch) > 0:
		pool = Pool(processes = min(nProc,nt))
		res=pool.map(DirectionalPairwiseStitch,ParamsStitch,chunksize=1)
		pool.close()
		pool.join()	
	for k in range(len(IJPairs)):
		dx,dy,dz,crr = res[k]
		i,j = IJPairs[k]
		D[i,j] = 1 - crr
		D[j,i] = D[i,j]
				
	# use minimum spanning tree to find the optimal overlapList. 
	D[where(isnan(D))] = 2.0	# TEMPORARY FIX, crr RETURNS NAN SOMETIMES, GETTING RID OF THEM. 

	D = D.astype(float32)
	E = zeros(2*nt+1).astype(int32)			
	NE = zeros(1).astype(int32)
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	# set parameter types
	sfmLib.kruskal.argtypes = [ctypes.c_long,ctypes.c_float, ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_int)]
	sfmLib.kruskal.restype = None
	# call the C function.
	sfmLib.kruskal(ctypes.c_long(nt),ctypes.c_float(threshold), \
					D.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
					E.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
	nEdges = E[-1];	# number of edges in the spanning tree. 
	E = E[0:2*nEdges].reshape(nEdges,2)								
	
	[i,j] = E[0,:]
	ii = tileIDs[i]
	jj = tileIDs[j]
	overlapList = [(ii,jj)]
	inIndR = [ii,jj]
	while 1:
		flag = 0
		for [i,j] in E:
			ii = tileIDs[i]
			jj = tileIDs[j]
			if (ii,jj) in overlapList or (jj,ii) in overlapList:
				continue
			if ii in inIndR:
				overlapList = overlapList + [(ii,jj)]
				inIndR.append(jj)
				flag = 1
				break
			elif jj in inIndR:
				overlapList	= overlapList + [(jj, ii)]
				inIndR.append(ii)
				flag = 1
				break
			else:
				continue
		if flag == 0:
			break
										
	stitchImages(overlapList = overlapList,iplot=0)	
	return overlapList
			 		
# delete old pairwise stitching results to reclauclate. 
def DeleteOldStitchingResults():
	files = glob.glob(filenameCommon+'*.'+filenameCommon+'*.txt')
	for ff in files:
		print 'Delete ',ff
		os.remove(ff)
				
def DirectionalPairwiseStitch(ParamsStitch):
	filenum1,filenum2,offsetPerc,offsetDirection = ParamsStitch
	filenameBase1 = filenameCommon+str(filenum1)
	filenameBase2 = filenameCommon+str(filenum2)
	outfile1 = filenameBase1+'.'+filenameBase2+'.txt'
	outfile2 = filenameBase2+'.'+filenameBase1+'.txt'
	
	# check if this done already. 
	if os.path.isfile(outfile1):
		print outfile1+'exists. Read previously sticthed result from it. Delete the file if need recalculation.'
		fpt = open(outfile1,'r')
		line = fpt.readline()
		fpt.close()
		ret = line.split(',')
		dx = float(ret[0])
		dy = float(ret[1])
		dz = float(ret[2])
		maxcorr = float(ret[3])
		print 'dx=',dx,' dy=',dy,' dz=',dz,' maxcorr=',maxcorr
		return dx,dy,dz,maxcorr
	
	print 'Computing offsets between ',filenameBase1,filenameBase2, ' offset =',offsetPerc, '% direction =', offsetDirection
	if  offsetDirection == 'up':
		offsetD = 1
	elif offsetDirection == 'down':
		offsetD = 2
	elif offsetDirection == 'left':
		offsetD = 3
	elif offsetDirection == 'right':
		offsetD = 4
	filename1 = filenameBase1+'.tif'
	filename2 = filenameBase2+'.tif'	
		
	overlapFract = 1 - offsetPerc/100.0	
	imageType = 0				
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libShuTu.so')
	sfmLib.DirectionalPairwiseStitchC.argtypes = [ctypes.c_long,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_float,\
						ctypes.c_long, ctypes.POINTER(ctypes.c_float)]
	sfmLib.DirectionalPairwiseStitchC.restype = None
	ret = array([0,0,0,0]).astype(float32)
	sfmLib.DirectionalPairwiseStitchC(ctypes.c_long(imageType),ctypes.c_char_p(filename1),ctypes.c_char_p(filename2),\
						ctypes.c_float(overlapFract),ctypes.c_long(offsetD),\
						ret.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))					
	dx,dy,dz,maxcorr = ret
	writeStitchOffsetToFile(outfile1,dx,dy,dz,maxcorr)	
	writeStitchOffsetToFile(outfile2,-dx,-dy,-dz,maxcorr)	
	return (dx,dy,dz,maxcorr)

# This file find stitched one tile to existing tiles. Useful for fixing stitching problem. 
def stitchToNearByTile(i,nx,ny,tileX,tileY,tileIDs,inIDs,tilesNotIn,overlapList):	
	crrs = []
	offs = []
	x = tileX[i]
	y = tileY[i]
	ii= tileIDs[i]
	# find the neighbors of tile i in already stitch tiles
	ind = where((abs(tileX[inIDs] - x) < nx/2) & (tileY[inIDs] < y) & ( tileY[inIDs] > y - ny))
	if len(ind[0]) > 0:
		j = inIDs[ind[0][0]]
		jj= tileIDs[j]
		offsetPerc = (abs(tileY[j] - y) * 100.0)/ny
		dx,dy,dz,crr = DirectionalPairwiseStitch((jj,ii,offsetPerc,'down'))
		crrs.append(crr)
		offs.append((dx,dy,dz,jj))
	ind = where((abs(tileY[inIDs] - y) < ny/2) & (tileX[inIDs] < x) & ( tileX[inIDs] > x - nx))
	if len(ind[0]) > 0:
		j = inIDs[ind[0][0]]
		jj= tileIDs[j]
		offsetPerc = (abs(tileX[j] - x) * 100.0)/nx
		dx,dy,dz,crr = DirectionalPairwiseStitch((jj,ii,offsetPerc,'right'))
		crrs.append(crr)
		offs.append((dx,dy,dz,jj))
	ind = where((abs(tileY[inIDs] - y) < ny/2) & (tileX[inIDs] > x) & ( tileX[inIDs] < x + nx))
	if len(ind[0]) > 0:
		j = inIDs[ind[0][0]]
		jj= tileIDs[j]
		offsetPerc = (abs(tileX[j] - x) * 100.0)/nx
		dx,dy,dz,crr = DirectionalPairwiseStitch((jj,ii,offsetPerc,'left'))
		crrs.append(crr)
		offs.append((dx,dy,dz,jj))
	ind = where((abs(tileX[inIDs] - x) < nx/2) & (tileY[inIDs] > y) & ( tileY[inIDs] < y + ny))
	if len(ind[0]) > 0:
		j = inIDs[ind[0][0]]
		jj= tileIDs[j]
		offsetPerc = (abs(tileY[j] - y) * 100.0)/ny
		dx,dy,dz,crr = DirectionalPairwiseStitch((jj,ii,offsetPerc,'up'))
		crrs.append(crr)
		offs.append((dx,dy,dz,jj))
	# select maximum overlap
	if len(crrs) == 0:
		tilesNotIn.append(i)
		return
	ind = crrs.index(max(crrs))
	dx,dy,dz,jj = offs[ind]
	overlapList.append((jj,ii))
	inIDs.append(i)
	return
	
def locateConnectionInOverlapList(overlapList,num):
	# this function locates the connection to tile num in overlapList
	overlapList = array(overlapList)
	for [n1, n2] in overlapList:
		if n2 == num:
			print n1, n2
			
def manualEditShifts(num1,num2):
	filename = filenameCommon+str(num1)+'.'+filenameCommon+str(num2)+'.txt'
	call(['pico',filename])
				
def writeStitchOffsetToFile(outfile,dx,dy,dz,crr):
	print 'Saving the offsets to ',outfile
	f = open(outfile,'w')
	f.write(str(dx)+','+str(dy)+','+str(dz)+','+str(crr))
	f.close()
	
def stitchImages(overlapList=overlapList,iplot=0):
	# This function sticthes the images together, and creates a summary flattened image of the neuron from the stacks.
	# if the offsets are not computed already, compute them in parallel.
	for (i,j) in overlapList:
		filenameBase1 = filenameCommon+str(i)
		filenameBase2 = filenameCommon+str(j)
		outfile = filenameCommon+str(i)+'.'+filenameCommon+str(j)+'.txt'
	# load the offsets from the computed results.
	sids = []
	for (i,j) in overlapList:
		sids.append(i)
		sids.append(j)
	sids=OrderedDict.fromkeys(sids).keys()
	sidsAppear = [[] for i in range(len(sids))]
	k =0
	offsets = []
	for (i,j) in overlapList:
		sidsAppear[sids.index(i)].append(k); sidsAppear[sids.index(j)].append(k); k += 1
		filenameBase1 = filenameCommon+str(i)
		filenameBase2 = filenameCommon+str(j)
		outfile = filenameCommon+str(i)+'.'+filenameCommon+str(j)+'.txt'
		print 'Loading previous stitching result from '+outfile
		# load the result
		f = open(outfile)
		lines = f.readlines()
		f.close()
		dx,dy,dz,maxcorr = lines[0].split(',')
		dx = float(dx); dy = float(dy); dz = float(dz); maxcorr = float(maxcorr)
		offsets.append(array([dx,dy,dz]))
	for ss in sidsAppear:
		ss = OrderedDict.fromkeys(ss).keys()
	# construct the coordinates of the images. 
	X = zeros(len(sids))
	Y = zeros(len(sids)) 
	Z = zeros(len(sids))
	shapes = []
	flags = zeros(len(sids))
	flags[0] = 1
	findConnectedXYZ(0,sids,sidsAppear,overlapList,offsets,flags,X,Y,Z)
	X = X - X.min()
	Y = Y - Y.min()
	Z = Z - Z.min()

	# get the dimensions of the stacks. 
	shapes = []
	for ii in sids:
		nx,ny,nz = checkImageSize(filenameCommon+str(ii)+'.tif')
		shapes.append((nx,ny,nz))
	
	# construct a stitched 2D image
	imFlatCombined = 255 * ones((int(X.max())+nx,int(Y.max())+ny,3))
	for kk in range(len(sids)-1,-1,-1):
		ii = sids[kk]
		imFlat = array(PIL.Image.open(filenameCommon+str(ii)+'.Proj.tif')).astype('uint8')
		if len(imFlat.shape) == 3:
			imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,:] = imFlat
		else:
			for jj in range(3):
				imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,jj] = imFlat					
	imFlatCombined = imFlatCombined.astype('uint8')	
	
	print 'Saving the combined images to ',filenameCommon+'.flatCombined.tif'
	im = PIL.Image.fromarray(array(imFlatCombined).astype('uint8'))
	im.save(filenameCommon+'.flatCombined.tif')
	
	outFile = filenameCommon+'StitchCoord.npz'
	print 'Saving stitch coordinates to ',outFile
	savez(outFile,X=X,Y=Y,Z=Z,IDs=sids,shapes=shapes,overlapList=overlapList,filenameCommon=filenameCommon)
	
	# save the coordiantes into a text file. 
	outFile = filenameCommon+'StitchCoord.txt'
	print 'Saving cooridnates in text format to ',outFile
	f = open(outFile,'w')
	f.write(str(len(X))+' '+str(nx)+' '+str(ny)+' '+str(nz)+'\n')
	for i in range(len(X)):
		f.write(str(sids[i])+' '+str(X[i])+' '+str(Y[i])+' '+str(Z[i])+'\n')
	f.close()
	
	# modify the project files 
	#getOverlap()
	
	# create the json file to be loaded into ShuTu GUI.
	filenameJS = filenameCommon+'.tiles.json'
	ddr = os.getcwd()
	print 'Writing a JSON format to ',filenameJS

	tileList = [];	
	for kk in range(len(sids)-1,-1,-1):
		ii = sids[kk]
		tile = {"source": filenameCommon+str(ii)+'.tif', \
			"offset": [Y[kk], X[kk], Z[kk]] , \
			"size": [int(ny), int(nx), 1], \
			"image" : filenameCommon+str(ii)+'.Proj.tif'}
		tileList.append(tile)
	tileData = dict({"Tiles": tileList})
	with open(filenameJS, 'w') as f:
		json.dump(tileData, f)

	if iplot == 1:
		figure(25)
		cla(); subplot(1,1,1); imshow(imFlatCombined)
		for ii in range(len(X)):	# write the tile IDs. 
			text(Y[ii]+ny/2.0,X[ii]+nx/2.0,str(sids[ii]))

def getOverlap(filenameCommon=filenameCommon):
	# load the coordinates genrated from stitchImages
	outFile = filenameCommon+'StitchCoord.npz'
	print 'Loading stitch coordinates from ',outFile
	try:
		res=numpy.load(outFile)
		X = res['X']
		Y = res['Y']
		Z = res['Z']
		fileIDs = res['IDs']
		shapes = res['shapes']
	except:
		print 'ERROR: No results from image stitching found. First run stitchImages.'
		return

	# The order of stitching the files SWCs at the end. Assumes reverse ordering. This must agree with stitchSWCs(). 
	fileOrder = range(len(fileIDs)-1,-1,-1)
	nx,ny,nz = shapes[0]
	offset = 20
	for ii in range(len(fileOrder)):
		fid = fileOrder[ii]
		mask = ones((nx,ny))
		x = X[fid]
		y = Y[fid]
		for jj in range(ii+1,len(fileOrder)):
			iid = fileOrder[jj]
			x2 = X[iid]
			y2 = Y[iid]
			if x2 + nx <= x or x2 >= x + nx:
				continue
			if x2 <= x:
				xs = 0; xe = max(1,x2 + nx - x - offset)
			else:
				xs = min(nx-1,x2 - x + offset); xe = x + nx
			if y2 + ny <= y or y2 >= y + ny:
				continue
			if y2 <= y:
				ys = 0; ye = max(1,y2 + ny - y - offset)
			else:
				ys = min(ny-1,y2 - y + offset); ye = y + ny
			mask[xs:xe,ys:ye] = 0
		fn = filenameCommon+str(fileIDs[fid])+'.roi.tif'
		print 'Save mask to ',fn
		im = PIL.Image.fromarray((mask*255).astype('uint8'))
		im.save(fn)

def replaceElelementInOverlapList(overlapList,(num1,num2),(num1R,num2R)):
	# this function replaces (num1,num2) in overlpaList to (num1R,num2R)
	overlapList = array(overlapList)
	for ii in range(len(overlapList)):
		if (overlapList[ii] == [num1,num2]).all():
			print 'Found ',overlapList[ii],' replace with ',[num1R,num2R]
			overlapList[ii] = [num1R,num2R]
			break
	return overlapList

def ManualPairwiseSticting(filenum1,filenum2,offsetDirection,offsetPerc=80):
	# offsetDirection for filenum1 relative to filenum2, 
	# offsetDirection, 'up', 'down', 'left', 'right'.	
	# offsetPerc - percent of movement in that direction. 
	filenameBase1 = filenameCommon+str(filenum1)
	filenameBase2 = filenameCommon+str(filenum2)
	im3d1,imFlat1,imFlatGRB1 = readImageData(filenameBase1)
	im3d2,imFlat2,imFlatGRB2 = readImageData(filenameBase2)	

	nx,ny,nz = im3d1.shape	

	dx =0; dy =0
	if  offsetDirection == 'up':
		dx = -offsetPerc/100.0 * nx
	elif offsetDirection == 'down':
		dx = offsetPerc/100.0 * nx
	elif offsetDirection == 'left':
		dy = -offsetPerc/100.0 * ny
	elif offsetDirection == 'right':
		dy = offsetPerc/100.0 * ny
		
	X = array([0,dx])
	Y = array([0,dy])
	X = X - X.min()
	Y = Y - Y.min()
	# construct a stiched 2D image,
	imFlatCombined = ones((int(X.max()+nx),int(Y.max()+ny)))
	imFlatCombined[int(X[0]):int(X[0]+nx),int(Y[0]):int(Y[0]+ny)] = imFlat1
	imFlatCombined[int(X[1]):int(X[1]+nx),int(Y[1]):int(Y[1]+ny)] = imFlat2
	figure(11)
	cla(); imshow(imFlatCombined,cm.gray)

	# regions of overlap.
	#figure(2); clf()
	nx,ny,nz = im3d1.shape
	if dx < 0:
		x11 = 0; x21 = -dx 
		x12 = int(nx+dx); x22 = nx
	else:
		x21 = 0; x11 = dx 
		x22 = nx- dx; x12 = nx
	if dy < 0:
		y11 = 0; y21 = -dy 
		y12 = ny+dy; y22 = ny
	else:
		y21 = 0; y11 = int(dy); 
		y22 = int(ny - dy); y12 = ny
	figure(12)
	subplot(2,1,1); imshow(imFlat1[x11:x12,y11:y12],cm.gray)
	subplot(2,1,2); imshow(imFlat2[x21:x22,y21:y22],cm.gray)

	# find Z shift. 
	print 'Computing the z shift...',
	dn = 5	# subsample x, y
	z1 = im3d1[x11:x12:dn,y11:y12:dn,:]
	z2 = im3d2[x21:x22:dn,y21:y22:dn,:]
	nx,ny,nz1 = z1.shape
	nx,ny,nz2 = z2.shape
	
	# brutal force computing correlzation for each shift, 
	zscan = range(-nz2+1,nz1-1)
	crrs = []
	for iz in zscan:
		if iz < 0:
			zp = int(min(nz1,nz2 + iz))
			z11 = 0; z12= z11+zp
			z21 = -iz; z22 = z21+zp
		else:
			zp = int(min(nz1 - iz,nz2))
			z11 = iz; z12 = z11 + zp
			z21 = 0; z22 = zp
		zr1 = z1[:,:,z11:z12]
		a,b,c = zr1.shape; zr1 = reshape(zr1,a*b*c)
		zr2 = z2[:,:,z21:z22]
		a,b,c = zr2.shape; zr2 = reshape(zr2,a*b*c)
		crrs.append(corrcoef(zr1,zr2)[0][1])
	crrs = array(crrs)
	figure(13); subplot(2,1,1);cla()
	plot(crrs)	
	dz = zscan[argmax(crrs)]
	print 'Brutal force method returned dz=',dz,' corr=',crrs.max()
	outfile = filenameBase1+'.'+filenameBase2+'.txt'
	print 'Saving the offsets to ',outfile
	f = open(outfile,'w')
	f.write(str(dx)+','+str(dy)+','+str(dz))
	f.close()
	return dx,dy,dz
			 			 
def findConnectedXYZ(kk,sids,sidsAppear,overlapList,offsets,flags,X,Y,Z):
	# find a connected image that has been placed.
	for jj in sidsAppear[kk]:
		(m,n) = overlapList[jj]
		if m == sids[kk]:
			sgn = 1
		else:
			n = m
			sgn = -1
		nind = sids.index(n)
		if flags[nind] == 1:
			continue
		else:
			flags[nind] = 1
			X[nind] = X[kk] + sgn * offsets[jj][0]
			Y[nind] = Y[kk] + sgn * offsets[jj][1]
			Z[nind] = Z[kk] + sgn * offsets[jj][2]
			findConnectedXYZ(nind,sids,sidsAppear,overlapList,offsets,flags,X,Y,Z)
	return

# a useful function for plotting the 2D projection of the entire neuron with tile numbers. 
# this can be used for debugging stitching or autotracing programs. 
def plotProjectionNeuron():
	# prepare for combining the eigenvalues and edges. 
	dataFile = filenameCommon+'StitchCoord.npz'
	print 'Loading stitch coordinates from ',dataFile
	res = numpy.load(dataFile)
	X = res['X'].astype(int)
	Y = res['Y'].astype(int)
	Z = res['Z'].astype(int)
	sids = res['IDs']
	nx,ny,nz=res['shapes'][0]	
	imFlatCombined = 255 * ones((int(X.max())+nx,int(Y.max())+ny,3))
	for kk in range(len(sids)-1,-1,-1):
		ii = sids[kk]
		imFlat = array(PIL.Image.open(filenameCommon+str(ii)+'.Proj.tif')).astype('uint8')
		if len(imFlat.shape) == 3:
			imFlat[0,:,:] = 255; imFlat[nx-1,:,:] = 255	# mark border bright for use in shortest path calculation. 
			imFlat[:,0,:] = 255; imFlat[:,ny-1,:] = 255			
			imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,:] = imFlat
		else:
			imFlat[0,:] = 255; imFlat[nx-1,:] = 255	# mark border bright for use in shortest path calculation. 
			imFlat[:,0] = 255; imFlat[:,ny-1] = 255			
			for jj in range(3):
				imFlatCombined[int(X[kk]):int(X[kk])+nx,int(Y[kk]):int(Y[kk])+ny,jj] = imFlat					
	imFlatCombined = imFlatCombined.astype('uint8')	
	figure()
	cla(); subplot(1,1,1); imshow(imFlatCombined)
	for ii in range(len(X)):	# write the tile IDs. 
		text(Y[ii]+ny/2.0,X[ii]+nx/2.0,str(sids[ii]),color='red')		
	print 'Saving the combined images to ',filenameCommon+'.flatCombined.tif'
	im = PIL.Image.fromarray(array(imFlatCombined).astype('uint8'))
	im.save(filenameCommon+'.flatCombined.tif')

#---SWC Manipulations-------------
#

class SWCPoint:
	def __init__(self,PP):
		self.ID, self.Type, self.x, self.y, self.z, self.r, self.parentID, self.zSig = PP

def saveSWC(filenameSave,branches):
	print 'Saving SWC to file '+filenameSave 
	PPList = []
	IDList = []
	for br in branches:
		for PP in br:
			PPList.append(PP)
			IDList.append(PP.ID)
	inds = argsort(array(IDList))
	f = open(filenameSave,'w')
	for iid in inds:
		PP = PPList[iid]
		ll = str(int(round(PP.ID)))+' ' \
			+str(int(round(PP.Type)))+' '+str(PP.y)+' '+str(PP.x)+' ' \
			+str(PP.z)+' '+str(PP.r)+' '+str(int(round(PP.parentID)))+'\n'
		f.write(ll)
	f.close()

# fixing z mismatch between swc and images. 
# use this function to shift swc in z to match one image tile. 
# when fixing problems with z-mismatch in stitching, the constructed swc can be shifted in z. 	
def zShiftSWCMatch(matchTile):
	# load the coordinates
	fn = filenameCommon+'StitchCoord.npz'
	print 'Loading stitch coordinates from ',fn
	try:
		res=numpy.load(fn)
		X = res['X'].astype(int)
		Y = res['Y'].astype(int)
		Z = res['Z'].astype(int)
		fileIDs = res['IDs']
		shapes = res['shapes']
	except:
		print 'ERROR: No results from image stitching found. First run stitchImages.'
		return

	iid = where(fileIDs == matchTile)[0]
	if len(iid) == 0:
		print 'matchTile=',matchTile,' in the argument does not match any tile ids. '
	x0 = X[iid]
	y0 = Y[iid]
	z0 = Z[iid]
	nx,ny,nz = shapes[iid][0]
	
	swcFilename = filenameCommon+'.swc'
	print 'Shifting z in ',swcFilename,' by matching tile ',matchTile
	LinkedPoints = getLinkedPointsFromSWC(swcFilename)
	
	fn = filenameCommon+str(matchTile)
	im3d,imFlat,imFlatC = readImageData(fn);	
	imFlatYZ = im3d.min(axis = 0)
	ZZ0 = []
	thr = imFlatYZ.min() + (imFlatYZ.max() - imFlatYZ.min()) * 0.2
	for i in range(ny):
		for j in range(nz):
			if imFlatYZ[i,j] < thr:
				ZZ0.append(j)
				
	ZZ = []
	Pinds = []
	for ii in range(len(LinkedPoints)):
		PP = LinkedPoints[ii]
		tind = locateTileInd(PP.x,PP.y,X,Y,nx,ny)
		if tind == iid:
			z = PP.z - z0
			ZZ.append(z[0])
			Pinds.append(ii)
			
	shifts = int(median(ZZ0) - median(ZZ))
	shifts = linspace(shifts-nz/2,shifts+nz/2).astype(int)
	scrs = zeros(len(shifts))
	im3d = im3d.max() -  im3d
	for ii in range(len(shifts)):
		for kk in Pinds:
			PP = LinkedPoints[kk]
			x = int(PP.x - x0)
			y = int(PP.y - y0)
			z = int(PP.z - z0) + shifts[ii]
			if z < 0 or z > nz - 1:
				continue
			scrs[ii] += im3d[x,y,z]
	zshift = shifts[argmax(scrs)]		
	print 'Shifting z in ',swcFilename,' by ',zshift
	for PP in LinkedPoints:
		PP.z += zshift
	branches, branchConnIDs = createBranchesFromLinkedPoints(LinkedPoints) 	
	# save the swc file.
	saveSWC(filenameCommon+'.zshift.swc',branches)					

# add zshitf to all points in swc. 
def shiftSWCZPositions(swcFilename,zshift):
	print 'Shifting z in ',swcFilename,' by ',zshift
	LinkedPoints = getLinkedPointsFromSWC(swcFilename)
	for PP in LinkedPoints:
		PP.z += zshift
	branches, branchConnIDs = createBranchesFromLinkedPoints(LinkedPoints) 	
	saveSWC(swcFilename,branches)

# converts the neuron's swc from pixel to micron using the defaul converstions. 
def scaleSWC():
	print "Loading ",filenameCommon+'.swc'
	f = open(filenameCommon+'.swc')
	lines = f.readlines()
	f.close()
	Points = [];
	for line in lines:
		ll = line.split(' ')
		if ll[0] == '#':
			continue
		nn = int(float(ll[0]))	# label of the point
		tp = int(float(ll[1]))  # point type
		py = float(ll[2])	# note the inversion of x, y.
		px = float(ll[3])	
		z  = float(ll[4])	# z
		r = float(ll[5])	# radius of the sphere. 
		np = int(float(ll[6]))	# parent point id. 			
		# get the length in micron
		py *= xyDist; px *= xyDist; r = r*xyDist; z *= zDist
		Points.append([nn,tp,py,px,z,r,np])
	print 'Saving scaled SWC to file '+filenameCommon+'.scaled.swc'
	f = open(filenameCommon+'.scaled.swc','w')
	for [nn,tp,py,px,z,r,np] in Points:
		ll = str(int(nn))+' '+str(int(tp))+' '+str(py)+' '+str(px)+' '+str(z)+' '+str(r)+' '+str(int(np))+'\n'
		f.write(ll)
	f.close()

class LinkedPoint:
	def __init__(self,x,y,z,r,ID,Type,zSig):
		self.x = x; self.y = y; self.z = z; self.r = r; self.ID = ID; self.Type = Type; self.zSig = zSig; self.conn = set()
	def addConn(self,ID):
		self.conn.add(ID)
	def delConn(self,ID):
		self.conn.discard(ID)	
	def numConn(self):
		return len(self.conn)

def getLinkedPointsFromSWC(swcfilename,nx=0,ny=0,bw=None,bwz=None,bws=None,idOffset=0,xOffset=0,yOffset=0,zOffset=0):
	# reload the swc file
	f = open(swcfilename)
	lines = f.readlines()
	f.close()

	# parse the swc file
	LinkedPoints = []
	pIDs = []
	IDs = []
	for line in lines:
		ll = line.split(' ')
		if ll[0] == '#':
			continue
		nn = float(ll[0]) + idOffset	# label of the point
		tp = float(ll[1])  	# point type
		py = float(ll[2]) + yOffset	# note the inversion of x, y.
		px = float(ll[3]) + xOffset
		r = float(ll[5])	# radius of the sphere. 
		zsig = 1
		z = float(ll[4])
		if bw!= None and bwz != None and bws !=None:
			# find z. 
			zAveSizeMin = 2		# minimum size of pixels for averaging z around each point. 
			zz = []
			zs = []
			zAveSize = min(zAveSizeMin,int(round(r * 0.3)))
			for ii in range(max(0,int(px)-zAveSize),min(nx,int(px)+zAveSize+1)):
				for jj in range(max(0,int(py)-zAveSize),min(ny,int(py)+zAveSize+1)):
					if bw[ii,jj] == False:
						continue;
					zz.append(bwz[ii,jj])
					zs.append(bws[ii,jj])
			if len(zz) > 0:
				z = mean(zz)
				zsig = sum(zs) == len(zs)
		z  += zOffset	# z
		np = float(ll[6]) 
		if np != -1:
			np += idOffset	# parent point id. 			
		pIDs.append(np)
		IDs.append(nn)
		LinkedPoints.append(LinkedPoint(px,py,z,r,nn,tp,zsig))

	for ii in range(len(pIDs)):
		pid = pIDs[ii]
		if pid == -1:
			continue
		iid=IDs.index(pid)
		LinkedPoints[ii].addConn(pid)
		LinkedPoints[iid].addConn(IDs[ii])	
	return LinkedPoints
	
def repairLinkedPoints(LinkedPoints):
	# make sure that all IDs in the connected lists are in the List.
	IDs = [ PP.ID for PP in LinkedPoints ]
	IDs = array(IDs)
	for PP in LinkedPoints:
		for iid in list(PP.conn):
			if len(find(IDs ==iid)) == 0:
				PP.delConn(iid)
	return LinkedPoints			

def getAllChildBranchesLinkedPoints(jj,LinkedPointsAll,IDs,flags,branches,branchConnIDs):
	# this function gets all branches in the linked points starting from jj.
	br = dllist()
	brIDs = []
	while 1:
		PP = LinkedPointsAll[jj]
		if len(br) > 0:
			pid = br.last.value.ID
		else:
			pid = -1
		br.append(SWCPoint(array([PP.ID,PP.Type,PP.x,PP.y,PP.z,PP.r,pid,PP.zSig])))
		flags[jj] = 1
		# find the next point.
		kids = []
		for kk in PP.conn:
			iid = IDs.index(kk)
			if flags[iid] == 0:
				kids.append(kk)
		if len(kids) == 0:
			break
		elif len(kids) == 1:
			jj = IDs.index(kids[0])
		else:
			for kk in kids:
				iid = IDs.index(kk)
				bID  = getAllChildBranchesLinkedPoints(iid,LinkedPointsAll,IDs,flags,branches,branchConnIDs)
				brIDs.append(bID)	
			break		
	branches.append(br)
	branchConnIDs.append(brIDs)
	return len(branches)-1

def createBranchesFromLinkedPoints(LinkedPointsAll):
	IDs = []
	numConn = []
	for PP in LinkedPointsAll:
		IDs.append(PP.ID)
		numConn.append(len(PP.conn))
	flags = zeros(len(LinkedPointsAll))
	numConn = array(numConn)
	
	branches = []		# linked lists of the branches.
	branchConnIDs = []	# IDs of the branches connected to this branch. 
	
	while 1:
		ind = where(flags == 0)[0]
		if len(ind) == 0:
			break
		ii = argmax(-numConn[ind])
		jj = ind[ii]	
		#print ii, jj
		getAllChildBranchesLinkedPoints(jj,LinkedPointsAll,IDs,flags,branches,branchConnIDs)

	# connect branches
	for ii in range(len(branches)):
		br = branches[ii]
		for jj in branchConnIDs[ii]:
			if len(branches[jj]) == 0:
				continue
			branches[jj].first.value.parentID = br.last.value.ID
	# smooth z, radius
	#smoothBranches(branches)
	return branches, branchConnIDs

# poly is a list of (x,y) pairs. 
def point_inside_polygon(x,y,poly):
    n = len(poly)
    inside =False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside

# This function converts the swc into 3point soma format used in NEURON. 
# In the original swc, soma is denoted as Type 1. The outline of the 2D projections of the soma points are assumed to outline
# the shape of the soma in 2D. 	
# IMPORTANT: The swc is assumed to be pixel based. The scaled version is generated. 
def convertTo3PointSomaAndScale(swcFilename):	
	print "Make sure that this file is pixel based: ",swcFilename
	LinkedPoints = getLinkedPointsFromSWC(swcFilename)
	
	indSomaPoints = [i for i in range(len(LinkedPoints)) if LinkedPoints[i].Type == 1]
	IDsSoma = [int(LinkedPoints[i].ID) for i in indSomaPoints]
	indSomaConn = []
	for i in range(len(LinkedPoints)):
		if i in indSomaPoints:
			continue
		PP = LinkedPoints[i]
		for sid in IDsSoma:
			if sid in PP.conn:
				indSomaConn.append(i)
				PP.conn -= set(IDsSoma)
				break
	# project to 2D
	xmin = 1e10; xmax = -1e10; ymin = 1e10; ymax = -1e10
	for i in indSomaPoints:
		PP = LinkedPoints[i]
		xx = PP.x - PP.r 
		if xx < xmin:
			xmin = xx
		xx = PP.x + PP.r 
		if xx > xmax:
			xmax = xx
		yy = PP.y - PP.r 
		if yy < ymin:
			ymin = yy
		yy = PP.y + PP.r 
		if yy > ymax:
			ymax = yy
	xd = int(xmax - xmin) + 10
	yd = int(ymax - ymin) + 10
	bw = zeros((xd,yd))
	doneList = []
	xm = xmin - 5
	ym = ymin - 5
	for i in indSomaPoints:
		PP = LinkedPoints[i]
		for xx in range(int(PP.x - PP.r), int(PP.x + PP.r)+1):
			for yy in range(int(PP.y - PP.r), int(PP.y + PP.r)+1):
				if sqrt((xx - PP.x)**2 + (yy - PP.y)**2) <= PP.r:
					bw[int(xx - xm),int(yy - ym)] = 1
		# find connected soma points. 
		for ids in PP.conn:
			if ids in IDsSoma:
				j = indSomaPoints[IDsSoma.index(ids)]
				if (i,j) in doneList or (j,i) in doneList:
					continue
				doneList.append((i,j))
				PP2 = LinkedPoints[j]
				xmin2 = min(PP.x - PP.r, PP2.x - PP2.r);
				ymin2 = min(PP.y - PP.r, PP2.y - PP2.r);
				xmax2 = max(PP.x + PP.r, PP2.x + PP2.r);
				ymax2 = max(PP.y + PP.r, PP2.y + PP2.r);
				if PP.r < PP2.r:
					r1 = PP.r;  x1 = PP.x;  y1 = PP.y
					r2 = PP2.r; x2 = PP2.x; y2 = PP2.y
				else:
					r1 = PP2.r; x1 = PP2.x; y1 = PP2.y				
					r2 = PP.r;  x2 = PP.x;  y2 = PP.y
				if abs(r2 - r1) < 1e-3:
					theta = pi/2
				else:		
					xc = (r2 * x1 - r1 * x2)/(r2 - r1)				
					yc = (r2 * y1 - r1 * y2)/(r2 - r1)
					d1 = sqrt((x1 - xc)**2 + (y1 - yc)**2)
					theta = arccos(r1/d1)
				beta = arctan2(y1 - y2, x1 - x2)
				cs = cos(beta); ss = sin(beta)
				x0 = r1 * cos(theta); y0 = r1 * sin(theta)
				x11 = x1 + cs * x0 - ss * y0; y11 = y1 + ss * x0 + cs * y0
				x0 = r1 * cos(theta); y0 = -r1 * sin(theta)
				x12 = x1 + cs * x0 - ss * y0; y12 = y1 + ss * x0 + cs * y0			  				
				x0 = r2 * cos(theta); y0 = r2 * sin(theta)
				x21 = x2 + cs * x0 - ss * y0; y21 = y2 + ss * x0 + cs * y0
				x0 = r2 * cos(theta); y0 = -r2 * sin(theta)
				x22 = x2 + cs * x0 - ss * y0; y22 = y2 + ss * x0 + cs * y0			  				
				poly = [(x11,y11),(x21,y21),(x22,y22),(x12,y12)]
				for xx in range(int(xmin2),int(xmax2)+1):
					for yy in range(int(ymin2),int(ymax2)+1):
						if point_inside_polygon(xx,yy,poly):	
							bw[int(xx - xm),int(yy - ym)] = 1				
	# get the list of the boundary points. 
	bpx = []; bpy = []
	for i in range(1,xd-1):
		for j in range(1,yd-1):
			ss = sum(bw[i-1:i+2,j-1:j+2]);
			if ss > 0 and ss < 9:
				bpx.append(i)
				bpy.append(j)
	bpx = array(bpx); bpy = array(bpy)			
	xs = mean(bpx)
	ys = mean(bpy)
	ss = 0
	for i in range(len(bpx)):
		ss += sqrt((bpx[i] - xs)**2 + (bpy[i] - ys)**2)
	rs = ss/len(bpx)
	ss = 0; ssr = 0;
	for i in indSomaPoints:
		ss += LinkedPoints[i].z * LinkedPoints[i].r
		ssr += LinkedPoints[i].r
	zs = ss/ssr
	xs = xm + xs
	ys = ym + ys
	
	# now create 3 point soma. 
	xs *= xyDist; ys *= xyDist; zs *= zDist; rs *= xyDist
	P1 = LinkedPoint(xs,ys,zs,rs,1,1,1)
	P2 = LinkedPoint(xs,ys-rs,zs,rs,2,1,1)
	P3 = LinkedPoint(xs,ys+rs,zs,rs,3,1,1)
	P2.addConn(P1.ID); P1.addConn(P2.ID)
	P3.addConn(P1.ID); P1.addConn(P3.ID)
	LinkedPointsConv = [P1,P2,P3]
	for i in range(len(LinkedPoints)):
		if i in indSomaPoints:
			continue
		PP = LinkedPoints[i]
		PP.ID += 3
		conn = [jj+3 for jj in PP.conn]	
		PP.conn = set(conn)
		PP.x *= xyDist
		PP.y *= xyDist
		PP.z *= zDist
		PP.r *= xyDist
		if i in indSomaConn:
			PP.addConn(P1.ID)
			P1.addConn(PP.ID)
		LinkedPointsConv.append(PP)
	# created the new swc, also scale distances. 
	branches, branchConnIDs = createBranchesFromLinkedPoints(LinkedPointsConv)
	outFilename = swcFilename + '.3PSomaScaled.swc'
	print "Saving the SCALED swc with three point soma into ",outFilename
	saveSWC(outFilename,branches);
	
#---Analyze log file----
#
#	Analyze log file from ShuTu software for number of operation in tracing neuron. 
#
import datetime
import time
def analyzeLogFile():
	
	# parameters
	dt = 0.1 # if time between consecutive same actions smaller than this time the are the same action (this is how the log file is generated.
	iuseStartEndTime = 0	# set to 1 if specifying start and end time.
	if iuseStartEndTime == 1:
		startDate ='2016-06-07'	# start time. keep this format, year-month-day
		startTime ='17:43:00'	# start time. keep this format, hour:minute:second
		endDate ='2016-06-07'	# start time. keep this format, year-month-day
		endTime ='18:08:00'	# start time. keep this format, hour:minute:second
		x = time.strptime(startDate+'T'+startTime,'%Y-%m-%dT%H:%M:%S')
		ts = time.mktime(x)
		x = time.strptime(endDate+'T'+endTime,'%Y-%m-%dT%H:%M:%S')
		tend = time.mktime(x)
	
	#logFilename = '/home/djin/NeuronReconstructionData/DH_10-9-13_100xredo/DH100913X100-.auto.swc.Edit.log.txt'
	#logFilename = '/home/djin/NeuronReconstructionData/DH_7-6-13-2_100x/DH070613C2X100-.auto.swc.Edit.log.txt'
	#logFilename = '/home/djin/NeuronReconstructionData/DH_7-8-13_100x/DH070813100x.auto.swc.Edit.log.txt'
	logFilename = '/home/djin/.neutube.z/log.txt'
	
	editT = []
	editAct = []
	editActUnique = []
	editActC = []
	
	# parse the log file.
	tt = -1
	t0 = -1
	previousAct = ''
	fileana = ''	
	fpt = open(logFilename,'r')
	timeAutoSave = 0;
	fflag = 0
	tt0 = 0
	for line in fpt:
		if line.find('INFO  ') == -1:
			continue

		timestr = line[6:6+19]
		x = time.strptime(timestr,'%Y-%m-%dT%H:%M:%S')
		t = time.mktime(x)

		if line.find('Autosave triggered') != -1 and fflag == 0:
			fflag = 1
			tt0 = t
		if line.find('Autosave triggered') == -1 and fflag == 1:
			timeAutoSave += t - tt0
			fflag = 0	
					
		iid = line.find('"')
		if iid == -1:
			continue			
		if iuseStartEndTime == 1 and (t < ts or t > tend):
			continue
		actstr=line[iid+1:-3]
		if actstr.find('Start reconstruction:') == 0:
			if t0 != -1:
				print 'Another reconstruction started. Stop.'
				break
			t0 = t
			fileana = actstr[30:-3]
			print 'Analyzing reconstruction of neuron in ',fileana
		if t0 == -1:
			continue	
		if t - tt <= dt and previousAct == actstr:
			continue
		else:	
			tt = t
			previousAct = actstr			
		editT.append(t-t0)
		editAct.append(actstr)
	fpt.close()

	# get unique actions and the counts. 
	editActUnique = []
	editActCounts = []
	for i in range(len(editAct)):
		aa = editAct[i]
		flag = 0
		for k in range(len(editActUnique)):
			bb = editActUnique[k]
			if aa == bb:
				editActCounts[k] += 1
				flag = 1
				break;
		if flag == 0:
			editActUnique.append(aa)
			editActCounts.append(1)
	editActCounts = array(editActCounts)		
	sind = argsort(editActCounts)[::-1]					

	# output the results. 		
	print 'Total number of edits: ',len(editAct)
	mm = (editT[-1] - timeAutoSave)/60.0	# time between autosaves are typically idle times. Remove them.
	hh = (editT[-1] - timeAutoSave)/3600.0
	print 'Total time: ',mm,' minutes or ',hh,' hours.'
	print 'List of actions (counts, action):'
	for iid in sind:
		print '  ',editActCounts[iid],'	',editActUnique[iid]
		
	# plot the tally
	figure(11); clf()
	ys = [10,9,8,7,6,5,4,3,2,1]
	ainds = list(sind[0:len(ys)])
	xs = editActCounts[ainds]
	barh(ys,xs,align='center',height=0.3,color='green')
	for i in range(len(ys)):
		text(0,ys[i]+0.25,editActUnique[ainds[i]])
		text(-2,ys[i]-0.15,editActCounts[ainds[i]],horizontalalignment='right') 
	text(-2,0,sum(editActCounts),horizontalalignment='right')
	text(0,0,'= Total number of edits')	
	axis('off')

# number of SWC points added in manual editing phases. 			
def countAdded():	
	swcFilename = 'testImages/DH070613C2X100-.auto.Edit.swc'
	typeAuto = 3		# SWC type of autotraced, all others are manual addition.  
	fpt = open(swcFilename,'r')
	lines = fpt.readlines()
	fpt.close()
	numSWCPoints = 0
	numAdded = 0
	for line in lines:
		if len(line) > 0 and line[0] == '#':
			continue
		numSWCPoints += 1
		if (int(line.split(' ')[1]) != typeAuto):
			numAdded += 1
	print 'Total number of SWC points = ',numSWCPoints
	print 'Total number of SWC points added = ',numAdded
	print 'Fraction of added points = ',numAdded*1.0/numSWCPoints
				
	
