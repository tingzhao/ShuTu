#test programs

import ctypes
from scipy.ndimage import gaussian_filter
from numpy import *
from pylab import *
from NeuTu import *
	
sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')

####test 1d gradient flow program. 
def testGVF1d():
	dt = 0.00001
	mu = 0.2

	x = linspace(0,1,100)
	dx = x[1] - x[0]
	n = len(x)
	maxIter = 1000

	# check the convergence criteria. 

	r = mu * dt /(dx*dx)
	print 'r=',r

	# two Gaussians
	#sigma1 = 0.001
	#sigma2 = 0.005
	#h1 = 0.5
	#h2 = 2
	#c1 = 0.4
	#c2 = 0.6
	#I = h1 * exp(-(x - c1)**2/(2*sigma1)) + h2 * exp(-(x -c2)**2/(2*sigma2)) + 0.1 * rand(n)
	
	# flat top distribution. 
	c1 = 0.3
	c2 = 0.7
	sigma = 0.02
	I = (1 + tanh((x - c1)/sigma)) * (1 + tanh((c2 - x)/sigma)) + 0.1 * rand(n)
	 
	Ix = zeros(n)
	for i in range(1,n-1):
		Ix[i] = (I[i+1] - I[i-1])/(2 * dx)
	Ix[0] = (I[1] - I[0])/dx
	Ix[n-1] = (I[n-1] - I[n-2])/dx
	figure(1); clf()
	subplot(221); plot(x,I); title('I(x)')
	subplot(222); plot(x,Ix); title('Ix')
		
	# load the C program
	# set parameter types
	Ix = Ix.astype(float32)
	v = zeros(n).astype(float32)
	sfmLib.gvf1d.argtypes = [ctypes.c_long, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
									ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_long]
	sfmLib.gvf1d.restype = None
	# call the C function.
	sfmLib.gvf1d(ctypes.c_long(n), v.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Ix.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				ctypes.c_float(mu), ctypes.c_float(dt), ctypes.c_float(dx), ctypes.c_long(maxIter))
				
	subplot(222); plot(x,v)
	# second order derivative
	v2 = zeros(n);
	for i in range(1,n-1):
		v2[i] = (v[i+1] - v[i-1])/(2.0*dx)
	v2[0] = (v[1] - v[0])/dx
	v2[n-1] = (v[n-1] - v[n-2])/dx	
	subplot(223); plot(x,v2); title('Ixx')

###test 2d gradient flow program
def testGVF2d():
	dt = 0.001
	mu = 10
	m = 100
	n = 100
	dx = 1
	dy = 1
	maxIter = 1000
	# check the convergence criteria. 
	r = mu * dt /(dx*dx)
	print 'r=',r

	sigma1 = 0.05		# first ring
	sigma2 = 0.05		# second ring
	sigma3 = 0.1		# a blob at the center
	h1 = 0.1			# height of the outer ring
	h2 = 0.1			# height of the innner ring
	h3 = 0.5			# height of the center blob.
	hrandom = 0.05		# height of noise
	r1 = 0.8
	r2 = 0.6

	I = hrandom * rand(m,n)
	for i in range(m):
		for j in range(n):
			x = (i - m/2.0)/max(m,n)*2
			y = (j - n/2.0)/max(m,n)*2
			r = sqrt(x*x + y*y)
			I[i,j] += h1 * exp(-(r - r1)**2/(2*sigma1**2)) + h2 * exp(-(r -r2)**2/(2*sigma2**2)) \
					 + h3 * exp(-r**2/(2*sigma3**2))
	I = I.max() - I
	for i in range(m):
		for j in range(n):
			if I[i,j] < 0.3:	# mimick satuartion
				I[i,j] = 0.3 + hrandom * rand()		 
	f1=figure(1); clf()
	subplot(421); imshow(I,cm.gray); title('I')
	
	# compute Ix and Iy. 
	sigma = 1	# sigma of the Gaussian for smoothing for taking the derivatives. 
	Is = gaussian_filter(I,sigma)
	Ix = zeros((m,n))
	Iy = zeros((m,n))
	for i in range(1,m-1):
		for j in range(1,n-1):
			Ix[i,j] = (Is[i+1,j] - Is[i-1,j])/(2.0*dx)
			Iy[i,j] = (Is[i,j+1] - Is[i,j-1])/(2.0*dy)
	for j in range(n):
		Ix[0,j] = (Is[1,j] - Is[0,j])/dx
		Ix[m-1,j] = (Is[m-1,j] - Is[m-2,j])/dx
	for i in range(m):
		Iy[i,0] = (Is[i,1] - Is[i,0])/dy
		Iy[i,n-1] = (Is[i,n-1] - Is[i,n-2])/dy
	subplot(422); imshow(Ix,cm.gray); title('Ix')
	subplot(423); imshow(Iy,cm.gray); title('Iy')
				
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# set parameter types
	Ix = Ix.astype(float32)
	Iy = Iy.astype(float32)
	v = zeros((m,n)).astype(float32)
	u = zeros((m,n)).astype(float32)	
	sfmLib.gvf2d.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
							ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_long]
	sfmLib.gvf2d.restype = None
	# call the C function.
	sfmLib.gvf2d(ctypes.c_long(m), ctypes.c_long(n), \
				v.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				u.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Ix.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Iy.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				ctypes.c_float(mu), ctypes.c_float(dt), ctypes.c_float(dx), ctypes.c_float(dy), ctypes.c_long(maxIter))
	v = gaussian_filter(v,sigma)
	u = gaussian_filter(u,sigma)			
	subplot(424); imshow(v,cm.gray); title('v')
	subplot(425); imshow(u,cm.gray); title('u')
	
	# calculate the Hessian matrix
	Ixx = zeros((m,n))
	Ixy = zeros((m,n))
	Iyy = zeros((m,n))
	for i in range(1,m-1):
		for j in range(1,n-1):
			Ixx[i,j] = (v[i+1,j] - v[i-1,j])/(2.0*dx)
			Iyy[i,j] = (u[i,j+1] - u[i,j-1])/(2.0*dy)
			Ixy[i,j] = (v[i,j+1] - v[i,j-1])/(2.0*dy)/2.0 + (u[i+1,j] - u[i-1,j])/(2.0*dx)/2.0
	# compute the eigen values of the Hessian matrix. 	
	lambda1 = 0.5 * (Ixx + Iyy + sqrt(pow(Ixx - Iyy,2) + 4 * Ixy * Ixy))
	lambda2 = 0.5 * (Ixx + Iyy - sqrt(pow(Ixx - Iyy,2) + 4 * Ixy * Ixy))
	lambda1c = lambda1.copy()	# important! 
	ind = where(lambda1 <= 0.0)
	lambda1c[ind] = 0.0
	lambda1c **= 0.2	# compress
	subplot(426); imshow(lambda1c,cm.gray)
		
	# get the valleyness measure. 
	alpha = 10
	lambda1c = lambda1.copy()	# important! 
	mm = lambda1.max() * 0.01
	ind = where(lambda1 <= mm )
	lambda1c[ind] = 0.0
	lambda1c **= 0.01	# compress
	
	# suppress blobs. 
	ind = where(lambda1 >= mm )
	lambda1c[ind] *= exp(-alpha*abs(lambda2[ind])/lambda1[ind])
	# supress edge
	dd = (v*v + u* u)**0.2
	dd = dd/dd.max()
	dd = dd/dd.max()
	lambda1c[ind] = lambda1c[ind] * exp(-dd[ind]*alpha)		
	
	subplot(427)
	imshow(lambda1c,cm.gray); title('lambda1c')
	subplot(428)
	bw = zeros(lambda1c.shape)
	bw[where(lambda1c > lambda1c.max() * 0.1)] = 1
	imshow(bw,cm.gray)
	f1.show()
						
	f2=figure(2); clf()
	subplot(221); plot(I[m/2,:])
	subplot(222); plot(lambda1[m/2,:])
	subplot(223); plot(lambda2[m/2,:])
	subplot(224); plot(lambda1c[m/2,:])
	f2.show()
	
	f3=figure(3); clf()
	l1 = reshape(lambda1,m*n)
	l2 = reshape(lambda2,m*n)
	plot(l1,l2,'.')
	f3.show()
	
###test modified gradient vector flow. 
####test 1d gradient flow program. 
def testModifedGVF1d():
	dt = 0.00001
	mu = 0.5

	x = linspace(0,1,100)
	dx = x[1] - x[0]
	n = len(x)
	maxIter = 1000

	# check the convergence criteria. 

	r = mu * dt /(dx*dx)
	print 'r=',r

	# two Gaussians
	#sigma1 = 0.001
	#sigma2 = 0.005
	#h1 = 0.5
	#h2 = 2
	#c1 = 0.4
	#c2 = 0.6
	#I = h1 * exp(-(x - c1)**2/(2*sigma1)) + h2 * exp(-(x -c2)**2/(2*sigma2)) + 0.1 * rand(n)
	
	# flat top distribution. 
	c1 = 0.3
	c2 = 0.7
	sigma = 0.02
	I = (1 + tanh((x - c1)/sigma)) * (1 + tanh((c2 - x)/sigma)) + 0.1 * rand(n)
	 
	Ix = zeros(n)
	for i in range(1,n-1):
		Ix[i] = (I[i+1] - I[i-1])/(2 * dx)
	Ix[0] = (I[1] - I[0])/dx
	Ix[n-1] = (I[n-1] - I[n-2])/dx
	figure(1); clf()
	subplot(221); plot(x,I); title('I(x)')
	subplot(222); plot(x,Ix); title('Ix')
		
	# load the C program
	# set parameter types
	Ix = Ix.astype(float32)
	v = zeros(n).astype(float32)
	sfmLib.gvmf1d.argtypes = [ctypes.c_long, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
									ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_long]
	sfmLib.gvmf1d.restype = None
	# call the C function.
	sfmLib.gvmf1d(ctypes.c_long(n), v.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Ix.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				ctypes.c_float(mu), ctypes.c_float(dt), ctypes.c_float(dx), ctypes.c_long(maxIter))
				
	subplot(222); plot(x,v)
	# second order derivative
	v2 = zeros(n);
	for i in range(1,n-1):
		v2[i] = (v[i+1] - v[i-1])/(2.0*dx)
	v2[0] = (v[1] - v[0])/dx
	v2[n-1] = (v[n-1] - v[n-2])/dx	
	subplot(223); plot(x,v2); title('Ixx')
	
	
###test squared Euclidean distance transformation. 
def testSEDT():

	m = 100
	n = 100
	sigma1 = 0.1		# first ring
	sigma2 = 0.1		# second ring
	sigma3 = 0.1		# a blob at the center
	h1 = 0.5			# height of the outer ring
	h2 = 0.5			# height of the innner ring
	h3 = 0.5			# height of the center blob.
	hrandom = 0.05		# height of noise
	r1 = 0.8
	r2 = 0.6

	I = hrandom * rand(m,n)
	for i in range(m):
		for j in range(n):
			x = (i - m/2.0)/max(m,n)*2
			y = (j - n/2.0)/max(m,n)*2
			r = sqrt(x*x + y*y)
			I[i,j] += h1 * exp(-(r - r1)**2/(2*sigma1**2)) + h2 * exp(-(r -r2)**2/(2*sigma2**2)) \
					 + h3 * exp(-r**2/(2*sigma3**2))
	I = I.max() - I

	bw = zeros((m,n))
	thr = I.max()*0.5
	bw[where(I < thr)] = 1
	
	figure(1); clf()
	subplot(221); imshow(I,cm.gray); title('I')
	subplot(222); imshow(bw,cm.gray)

	# call the squared distance transformation. 	
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# set parameter types
	bw = bw.astype(float32)
	edt = zeros(bw.shape).astype(float32)
	sfmLib.sedt.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.sedt.restype = None
	# call the C function.
	sfmLib.sedt(ctypes.c_long(m), ctypes.c_long(n), \
				bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				edt.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	subplot(223); imshow(edt,cm.gray)
	subplot(224); plot(edt[m/2,:])
	
# test 3D shortest distance program. 
def testDijstraComputePath3D():
	# create a synthetic 3D image. 
	nx = 50
	ny = 50
	nz = 30
	im3d = ones((nx,ny,nz))
	# create an expanding helix. 
	R = 5
	h = 5
	r = 2
	dr = 0.4
	ic = nx/2
	jc = ny/2
	for k in range(nz):
		RR = R + k * dr		# radius of the circle on this height. 
		ix = ic + int(RR * cos(k * 1.0/h))
		iy = jc + int(RR * sin(k * 1.0/h))
		for ii in range(max(0,ix - r),min(nx,ix+r)):
			for jj in range(max(0,iy - r),min(ny,iy+r)):
				if ((ii-ix)**2 + (jj-iy)**2) < r**2:
					im3d[ii,jj,k] = 0
		if k == 0:
			i1 = ix; j1 = iy; k1 = k	# starting point.
		elif k == nz-1:
			i2 = ix; j2 = iy; k2 = k	# end point. 			
			
	# projections
	imxy = im3d.min(axis = 2)
	imxz = im3d.min(axis = 1)
	imyz = im3d.min(axis = 0)
	figure(22); clf()
	subplot(221)
	imshow(imxy,cm.gray)
	subplot(222)
	imshow(imxz,cm.gray)
	subplot(223)
	imshow(imyz,cm.gray)
		
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	dists3d = zeros((nx,ny,nz)).astype(float32)
	im3d = im3d.astype(float32)
	zfact = 1.0
	sfmLib.dijstraComputeDists3D.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, \
											 ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_float, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.dijstraComputeDists3D.restype = None
	# call the C function.
	sfmLib.dijstraComputeDists3D(ctypes.c_long(i1), ctypes.c_long(j1), ctypes.c_long(k1),\
				ctypes.c_long(nx), ctypes.c_long(ny), ctypes.c_long(nz), ctypes.c_float(zfact),\
				im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				dists3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	# now compute the shortest path. 
	sfmLib.dijstraComputePath3D.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, \
											 ctypes.c_long, ctypes.c_long, ctypes.c_long,\
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.dijstraComputePath3D.restype = None
	# call the C function.
	bw3d = zeros((nx,ny,nz)).astype(float32)
	sfmLib.dijstraComputePath3D(ctypes.c_long(i2), ctypes.c_long(j2), ctypes.c_long(k2),\
				ctypes.c_long(nx), ctypes.c_long(ny), ctypes.c_long(nz), \
				dists3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				bw3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
							
	figure(23); clf()
	subplot(221)
	imshow(bw3d.max(axis=2),cm.gray)
	subplot(222)
	imshow(bw3d.max(axis=1),cm.gray)
	subplot(223)
	imshow(bw3d.max(axis=0),cm.gray)
	
def testGaussianFilter():
	nx = 100
	ny = 100
	Wx = 20
	Wy = 20
	sigma = 1.0
	img = zeros((nx,ny)).astype(float32)
	for i in range(nx):
		for j in range(ny):
			d = sqrt((i-nx/2)**2 + (j-ny/2)**2)
			if d > 20 and d < 30:
				img[i,j] = rand()				
	imgS = rand(nx,ny).astype(float32)	
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	sfmLib.gaussianFilter2D.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_ulong, ctypes.c_ulong, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),ctypes.c_float]
	sfmLib.gaussianFilter2D.restype = None
	# call the C function.
	sfmLib.gaussianFilter2D(ctypes.c_long(nx), ctypes.c_long(ny), ctypes.c_ulong(Wx), ctypes.c_ulong(Wy),\
				img.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				imgS.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), ctypes.c_float(sigma))
	figure(); subplot(211); imshow(img,cm.gray); subplot(212); imshow(imgS,cm.gray)
	
def testGetSparseThreshold():

	imFlat = imread("testImages/DH070613C2X100-40.tif.imFlat0.tif");
	nx, ny = imFlat.shape
	mmax = imFlat.max()
	mmin = imFlat.min()
	imFlat = (mmax - imFlat)*1.0/(mmax - mmin)
	imFlat *= 0.05	#compress
	sparseUp = 0.15

	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	imFlat = imFlat.astype(float32)
	sfmLib.getSparseThreshold.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.c_float]
	sfmLib.getSparseThreshold.restype = ctypes.c_float
	# call the C function.
	thr = sfmLib.getSparseThreshold(ctypes.c_long(nx), ctypes.c_long(ny), \
				imFlat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), ctypes.c_float(sparseUp))
	figure();
	subplot(211); imshow(imFlat,cm.gray)
	bw = imFlat > thr
	subplot(212); imshow(bw,cm.gray)
	
def testLabelObjectsBinaryImage():

	nx = 50; 
	ny = 100; 
	
	bw = zeros((nx,ny))	
	for i in range(nx):
		for j in range(ny):
			ic = nx/2
			m = 4
			for k in range(1,m):
				jc = ny/m*k
				dd = sqrt((i-ic)**2 + (j-jc)**2)
				if (dd < 10):
					bw[i,j] = 1.0
				
	figure(); subplot(211); imshow(bw,cm.gray)

	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	bw = bw.astype(float32)
	labels = zeros((nx,ny)).astype(float32)
	sfmLib.labelObjectsBinaryImage.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.labelObjectsBinaryImage.restype = ctypes.c_long
	# call the C function.
	nlabs = sfmLib.labelObjectsBinaryImage(ctypes.c_long(nx), ctypes.c_long(ny), \
				bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), labels.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	subplot(212); imshow(labels,cm.gray);
			
def testSkeletonizeZS():
	nx = 100
	ny = 100
	bw = zeros((nx,ny))
	ic = nx/2
	jc = ny/2
	r1 = 20
	r2 = 40
	for i in range(nx):
		for j in range(ny):
			dd = sqrt((i-ic)**2 + (j-jc)**2)
			if dd > r1 and dd < r2:
				bw[i,j] = 1.0
	figure(); subplot(211); imshow(bw,cm.gray)

	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	bw = bw.astype(float32)
	sfmLib.skeletonizeZS.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.skeletonizeZS.restype = None
	# call the C function.
	sfmLib.skeletonizeZS(ctypes.c_long(nx), ctypes.c_long(ny), \
				bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))

	subplot(212); imshow(bw,cm.gray);
	
def testPruneSmallBranches():
	nx = 100
	ny = 100
	smallLen = 10
	bw = zeros((nx,ny))
	for i in range(20,80):
		bw[i,ny/2] = 1.0
	for j in range(-10,10):
		bw[nx/2,ny/2+j] = 1.0
	for i in range(nx/2-5,nx/2+5):
		bw[i,ny/3] = 1.0	
	for theta in linspace(0,pi/2,100):
		i = nx/2 + int(20*cos(theta))
		j = ny/2 + int(20*sin(theta))
		bw[i,j] = 1.0
	
	im = imread('testImages/DH070613C2X100-40.tif.imFlat2.Mask.tif')
	bw = im[500:800,900:1200].copy()
	bw = bw/bw.max()
	nx,ny = bw.shape
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	bw = bw.astype(float32)
	sfmLib.skeletonizeZS.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
	sfmLib.skeletonizeZS.restype = None
	# call the C function.
	sfmLib.skeletonizeZS(ctypes.c_long(nx), ctypes.c_long(ny), \
				bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
		
	figure(); subplot(211); imshow(bw,cm.gray)
	
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# first compute shortest distances to all points from i1, j1, k1. 
	# set parameter types
	bw = bw.astype(float32)
	sfmLib.pruneSmallBranches.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.c_float]
	sfmLib.pruneSmallBranches.restype = None
	# call the C function.
	sfmLib.pruneSmallBranches(ctypes.c_long(nx), ctypes.c_long(ny), \
				bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.c_float(smallLen))

	subplot(212); imshow(bw,cm.gray);
	
def testChanVese():
	
	itest = 2	# 2, long thing rod. 
				# 1, circle.
				# 3, neuron image. 
				
	
	if itest == 1:
		# create a test image. 
		nx = 200
		ny = 200
		sigma = 5.0
		l = 100
		rc = 20
		noiseLevel = 0.1

		imFlat = zeros((nx,ny))
		for i in range(nx):
			for j in range(ny):
				dx = i - nx/2
				dy = j - ny/2
				r = sqrt(dx*dx + dy*dy)
				if r < rc:
					imFlat[i,j] = 1.0
				else:
					imFlat[i,j] = 1.0 * exp(-(r-rc)*(r-rc)/(sigma*sigma)) + noiseLevel*rand()
	elif itest == 2:
		# create a test image. 
		nx = 200
		ny = 200
		sigma = 5.0
		l = 100
		rc = 20
		noiseLevel = 0.1
		
		imFlat = zeros((nx,ny))
		for i in range(nx):
			for j in range(ny):
				dx = i - nx/2
				dy = j - ny/2
				if dy > 0:
					if dy > l/2:
						dy = dy - l/2
					else:
						dy = 0
				else:
					if dy <  -l/2:
						dy = dy + l/2
					else:
						dy = 0
				imFlat[i,j] = exp(-dx*dx/(sigma*sigma)) * exp(-dy*dy/(sigma*sigma)) + noiseLevel*rand()	
	elif itest == 3:
		#filename = 'testImages/OneChDH070813100x49.tif.imFlat3.tif'
		filename = 'testImages/OneChannelDH070613C2X100-40.tif.imFlat1.tif'
		imFlat = imread(filename).astype(float32)
		nx,ny = imFlat.shape			
		imFlatBack = gaussian_filter(imFlat,sigma=nx/20.0)
		imFlat = imFlat - imFlatBack
		# compress
		imFlat = pow(imFlat.max() - imFlat,0.5)
	
	# initial bw
	thr = imFlat.max() - (imFlat.max() - imFlat.min())*0.2
	bw = zeros((nx,ny))
	bw[where(imFlat > thr)] = 1.0
	
	# invert to bright field
	imFlat = imFlat.max() - imFlat
			
	f1=figure(1); clf(); 
	subplot(221);imshow(imFlat,cm.gray)
	subplot(222);imshow(bw,cm.gray)

	bw = bw.astype(float32)
	imFlat = imFlat.astype(float32)
	lam = 2
	maxit = 1000
  
	# call the C function.
	sfmLib.sparseFieldChanVese.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_float), \
				ctypes.POINTER(ctypes.c_float), ctypes.c_long, ctypes.c_float]
	sfmLib.sparseFieldChanVese.restype = None
	sfmLib.sparseFieldChanVese(ctypes.c_long(nx), ctypes.c_long(ny), imFlat.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
		bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
		ctypes.c_long(maxit), ctypes.c_float(lam))

	subplot(223); imshow(bw,cm.gray)
	f1.show()
	
def testLevelSet():
	
	# create a test image. 
	nx = 200
	ny = 200
	sigma = 20.0
	l = 100
	rc = 20
	noiseLevel = 0.01
	
	iExam = 1	# 1, circle
				# 2, bar
				# 3, cross
				# 4, neuron image. 
	
	if iExam == 1:
		imFlat = zeros((nx,ny))
		for i in range(nx):
			for j in range(ny):
				dx = i - nx/2
				dy = j - ny/2
				r = sqrt(dx*dx + dy*dy)
				if r < rc:
					imFlat[i,j] = 1.0
				else:
					imFlat[i,j] = 1.0 * exp(-(r-rc)*(r-rc)/(sigma*sigma)) + noiseLevel*rand()
	elif iExam == 2:
		imFlat = zeros((nx,ny))
		for i in range(nx):
			for j in range(ny):
				dx = i - nx/2
				dy = j - ny/2
				if dy > 0:
					if dy > l/2:
						dy = dy - l/2
					else:
						dy = 0
				else:
					if dy <  -l/2:
						dy = dy + l/2
					else:
						dy = 0
				imFlat[i,j] = exp(-dx*dx/(sigma*sigma)) * exp(-dy*dy/(sigma*sigma)) + noiseLevel*rand()	
	elif iExam == 3:
		imFlat = zeros((nx,ny))
		for i in range(nx):
			for j in range(ny):
				dx = i - nx/2
				dy = j - ny/2
				if dy > 0:
					if dy > l/2:
						dy = dy - l/2
					else:
						dy = 0
				else:
					if dy <  -l/2:
						dy = dy + l/2
					else:
						dy = 0
				imFlat[i,j] = exp(-dx*dx/(sigma*sigma)) * exp(-dy*dy/(sigma*sigma))
				dx = i - nx/2
				dy = j - ny/2
				if dx > 0:
					if dx > l/2:
						dx = dx - l/2
					else:
						dx = 0
				else:
					if dx <  -l/2:
						dx = dx + l/2
					else:
						dx = 0
				imFlat[i,j] += exp(-dx*dx/(sigma*sigma)) * exp(-dy*dy/(sigma*sigma))
				if imFlat[i,j] > 1:
					imFlat[i,j] = 1.0
				imFlat[i,j] += noiseLevel * rand()
				
	elif iExam == 4:
		filename = 'testImages/DH100913X100-47.tif.imFlat1.tif'
		#filename = 'testImages/testStack.tif.imFlat3.tif'
		imFlat = imread(filename).astype(float32)
		nx,ny = imFlat.shape
	
	if iExam != 4:					
		imFlat = imFlat.max() - imFlat	# invert.
	
	bw = zeros((nx,ny))
	if iExam != 4:
		for i in range(nx):
			for j in range(ny):
				if sqrt((i-nx/2)**2 + (j-ny/2)**2) < rc*3:
					bw[i,j] = 1
	else:
		bw[where(imFlat < imFlat.max()/4)] = 1
			
	Ix = zeros((nx,ny))
	Iy = zeros((nx,ny))
	for i in range(1,nx-1):
		for j in range(ny):
			Ix[i,j] = (imFlat[i+1,j] - imFlat[i-1,j])/2.0
	for j in range(1,ny-1):
		for i in range(nx):
			Iy[i,j] = (imFlat[i,j+1] - imFlat[i,j-1])/2.0
	for i in range(nx):
		Iy[i,0] = (4*imFlat[i,1] - 3*imFlat[i,0] - imFlat[i,2])/(2.0)
		Iy[i,ny-1] = (4*imFlat[i,ny-2] - 3*imFlat[i,ny-1] - imFlat[i,ny-3])/(-2.0)
	for j in range(ny):
		Ix[0,j] = (4*imFlat[1,j] - 3*imFlat[0,j] - imFlat[2,j])/(2.0)
		Ix[nx-1,j] = (4*imFlat[nx-2,j] - 3*imFlat[nx-1,j] - imFlat[nx-3,j])/(-2.0)
	
	g = (Ix * Ix + Iy * Iy)
	g = (g - g.min())/(g.max() - g.min())		
	g = 1.0/(1 + g**5);
	g = g.astype(float32)
	
	bwini = bw.copy();
	f1=figure(1); clf(); 
	subplot(221);imshow(imFlat,cm.gray)
	subplot(222);imshow(bwini,cm.gray)
	subplot(223);imshow(g,cm.gray)
	f1.show()	
		
	f2=figure(2); clf()
	subplot(321); imshow(imFlat,cm.gray)
	subplot(322); imshow(g,cm.gray)
	subplot(323); imshow(bwini,cm.gray)

	bw = bwini.copy()
	bw = bw.astype(float32)

	# parameters for levelset.
	mu = 1.0			# parameter for regularizing the levelset function
	maxiter = 4000

	sfmLib.sparseFieldLevelSet.argtypes = [ctypes.c_long,ctypes.c_long, \
						ctypes.POINTER(ctypes.c_float), \
						ctypes.POINTER(ctypes.c_float),ctypes.c_long,ctypes.c_float]
	sfmLib.sparseFieldLevelSet.restype = None
	sfmLib.sparseFieldLevelSet(ctypes.c_long(nx),ctypes.c_long(ny),\
						g.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),\
						bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.c_long(maxiter), \
						ctypes.c_float(mu))						
	subplot(324); imshow(bw,cm.gray);
	subplot(325); plot(bw[nx/2,:]); plot(g[nx/2,:])
	subplot(326); plot(bw[:,ny/2]); plot(g[:,ny/2])	
	f2.show()	

#This function creates simple 3D tiff stacks with artificial images	for test purpose. 
def createArtificialImageStacks():
	nx = 100
	ny = 200
	nz = 50
	im3d = zeros((nx,ny,nz))

	# create a tube. 
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				if i > 20 and i < nx - 20:
					r2 = (j - ny/2)**2 + (k - nz/2)**2
					if sqrt(r2) < 5:
						im3d[i,j,k] = 1.0
	
	filenameBase = 'testImages/testStack'
	createTiffStackFromArray(im3d,filenameBase)

# This function tests creation of 2D mask. 	
def testCreateMask2D():

	filename = 'testImages/OneChannelDH070613C2X100-40.tif'
	iBackgrounType = 0	# set to 0 if the image is bright background.			
			
	im3d,I,imRGB = loadStacks(filename,iBackgrounType);
	
	I = I.astype(float32)
	nx,ny = I.shape

	bw = zeros(I.shape).astype(float32)
	sfmLib.createMask2D.argtypes = [ctypes.c_long,ctypes.c_long, \
						ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float)]
	sfmLib.createMask2D.restype = None
	sfmLib.createMask2D(ctypes.c_long(nx),ctypes.c_long(ny),\
						I.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	
	fg = figure(10); 
	subplot(211); imshow(I,cm.gray)
	subplot(212); imshow(bw,cm.gray)
	fg.show()

# Here we explore methods of making 2D mask from projections. 	
def exploreCreateMask2D():

	filename = 'testImages/DH100913X100-47.tif.imFlat1.tif'
	filename = 'testImages/OneChannelDH070613C2X100-40.tif.imFlat1.tif'
	filename = 'testImages/OneChannelDH070613C2X100-40.tif.imFlat2.tif'
	I = imread(filename).astype(float32)
	m,n = I.shape
	sigma = 1	# sigma of the Gaussian for smoothing for taking the derivatives. 
	Is = gaussian_filter(I,sigma)
	
	# get rid of background variation and compress.
	sigmaBack = 200
	Iback = gaussian_filter(Is,sigma = sigmaBack)
	Is = Is - Iback
	Is = (Is.max() - Is)/(Is.max() - Is.min())
	Is = Is.max() - Is	

    # ger 2D gradient flow.  		
	dt = 0.02
	mu = 10
	maxIter = 1000
	dx = 1
	dy = 1
	# check the convergence criteria. 
	r = mu * dt /(dx*dx)
	print 'r=',r
	f1=figure(1); clf()
	subplot(421); imshow(I,cm.gray); title('I')
	
	# compute Ix and Iy. 
	Ix = zeros((m,n))
	Iy = zeros((m,n))
	for i in range(1,m-1):
		for j in range(1,n-1):
			Ix[i,j] = (Is[i+1,j] - Is[i-1,j])/(2.0*dx)
			Iy[i,j] = (Is[i,j+1] - Is[i,j-1])/(2.0*dy)
	for j in range(n):
		Ix[0,j] = (Is[1,j] - Is[0,j])/dx
		Ix[m-1,j] = (Is[m-1,j] - Is[m-2,j])/dx
	for i in range(m):
		Iy[i,0] = (Is[i,1] - Is[i,0])/dy
		Iy[i,n-1] = (Is[i,n-1] - Is[i,n-2])/dy
	subplot(422); imshow(Ix,cm.gray); title('Ix')
	subplot(423); imshow(Iy,cm.gray); title('Iy')
				
	# load the C program
	sfmLib = ctypes.CDLL(dirCodes+'libNeuTu.so')
	# set parameter types
	Ix = Ix.astype(float32)
	Iy = Iy.astype(float32)
	v = zeros((m,n)).astype(float32)
	u = zeros((m,n)).astype(float32)	
	sfmLib.gvf2d.argtypes = [ctypes.c_long, ctypes.c_long, \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
							ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), \
							ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_long]
	sfmLib.gvf2d.restype = None
	# call the C function.
	sfmLib.gvf2d(ctypes.c_long(m), ctypes.c_long(n), \
				v.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				u.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Ix.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				Iy.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
				ctypes.c_float(mu), ctypes.c_float(dt), ctypes.c_float(dx), ctypes.c_float(dy), ctypes.c_long(maxIter))
	v = gaussian_filter(v,sigma)
	u = gaussian_filter(u,sigma)			
	subplot(424); imshow(v,cm.gray); title('v')
	subplot(425); imshow(u,cm.gray); title('u')
	g = 0.5*(v*v + u*u)
	subplot(426); imshow(g**0.1,cm.gray); title('gradient')

	# calculate the Hessian matrix
	Ixx = zeros((m,n))
	Ixy = zeros((m,n))
	Iyy = zeros((m,n))
	for i in range(1,m-1):
		for j in range(1,n-1):
			Ixx[i,j] = (v[i+1,j] - v[i-1,j])/(2.0*dx)
			Iyy[i,j] = (u[i,j+1] - u[i,j-1])/(2.0*dy)
			Ixy[i,j] = (v[i,j+1] - v[i,j-1])/(2.0*dy)/2.0 + (u[i+1,j] - u[i-1,j])/(2.0*dx)/2.0
	# compute the eigen values of the Hessian matrix. 	
	lambda1 = 0.5 * (Ixx + Iyy + sqrt(pow(Ixx - Iyy,2) + 4 * Ixy * Ixy))
	lambda2 = 0.5 * (Ixx + Iyy - sqrt(pow(Ixx - Iyy,2) + 4 * Ixy * Ixy))

	lambdaRatioThr = 2
	expCompr = 0.1
	sparse = 0.15
	lambda1[where(lambda1 < lambdaRatioThr * abs(lambda2))] = 0
	lambda1  =  lambda1**expCompr
	thrSparse = percentile(lambda1.reshape(1,n*m),(1-sparse)*100)
	
	bw2 = zeros((m,n))
	bw2[where(lambda1 > thrSparse)] = 1

	# get mask from image intensitity
	bw = zeros((m,n))
	Is2 = (Is.max()-Is)**0.5
	bw[where(Is2 > 0.8)] = 1
	bw = bw2 + bw 
	bw[where(bw >0)] = 1	
	# smooth the initial mask
	#bw = gaussian_filter(3*(bw-0.5),sigma=3)
	#bw[where(bw > 0)] = 1
	#bw[where(bw <= 0)] = 0
		
	# chan-vese segmentation.
	#Is = Is.astype(float32)
	#bw = bw.astype(float32)
	#lam = 1.0
	#dt = 0.1
	#maxit = 1000
	# call the C function.
	#sfmLib.sparseFieldChanVese.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_float), \
	#			ctypes.POINTER(ctypes.c_float), ctypes.c_long, ctypes.c_float, ctypes.c_float]
	#sfmLib.sparseFieldChanVese.restype = None
	#sfmLib.sparseFieldChanVese(ctypes.c_long(m), ctypes.c_long(n), Is.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
	#	bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
	#	ctypes.c_long(maxit), ctypes.c_float(lam), ctypes.c_float(dt))
	#subplot(427); imshow(bw,cm.gray); title('Chanvese')
	
	# smoothing with level set
	# parameters for levelset.
	mu = 0.2			# parameter for regularizing the levelset function
	maxiter = 2000

	Is = Is.astype(float32)
	bw = bw.astype(float32)
	sfmLib.sparseFieldLevelSet.argtypes = [ctypes.c_long,ctypes.c_long, \
						ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float), \
						ctypes.POINTER(ctypes.c_float),ctypes.c_long,ctypes.c_float]
	sfmLib.sparseFieldLevelSet.restype = None
	sfmLib.sparseFieldLevelSet(ctypes.c_long(m),ctypes.c_long(n),\
						v.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),u.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
						bw.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.c_long(maxiter), \
						ctypes.c_float(mu))						
							
	subplot(428); imshow(bw,cm.gray);

# This function is for testing the validity of SWC points. 	
# IMPORTANT: in libNeuTu.c, set checkValidityOfLinkedPoint_DIAG to 1 to make this work! 
#			Then need to recompile and reload libNeuTu.so
# 			A code for plotting the diagnosis is generated and written to temp2.py	
def testCheckValidityOfLinkedPoint():

	iParSet = 2

	# parameters for the figure in the paper showing checkin validity of a SWC point. 
	if iParSet == 1:
		filename = 'testImages/OneChannelDH070613C2X100-20.tif'
		y = 621
		x = 317
		z = 76
		r = 8
		angle = 0
		iBackgrounType = 1	# set to 0 if the image is bright background.
	elif iParSet == 2:	# other test points. 
		filename = 'testImages/OneChannelDH070613C2X100-40.tif'
		x = 178
		y = 29
		z = 96
		r = 9
		angle = 0
		iBackgrounType = 1	# set to 0 if the image is bright background.			
	elif iParSet == 3: # bright background test points
		#filename='/home/djin/NeuronReconstructionData/DH_7-6-13-2_100x/DH070613C2X100-47.tif'
		filename='/home/djin/NeuronReconstructionData/DH_10-9-13_100xredo/DH100913X100-47.tif'
		x = 486
		y = 1231
		z = 81
		r = 19
		angle = 0
		iBackgrounType = 0	# set to 0 if the image is bright background.			
			
	im3d,imFlat,imRGB = loadStacks(filename,iBackgrounType);
	m,n,nz = im3d.shape
	I = im3d[:,:,z].copy().astype(float32) 
	I = I -  I.min()
	sigma = 1	# sigma of the Gaussian for smoothing for taking the derivatives. 
	#imFlat = gaussian_filter(I,sigma)

	nx,ny = imFlat.shape
	
	sfmLib.checkValidityOfLinkedPoint.argtypes = [ctypes.c_long,ctypes.c_long,ctypes.c_long,ctypes.c_long,\
						ctypes.c_long,ctypes.c_long,ctypes.c_long,\
						ctypes.POINTER(ctypes.c_float), \
						ctypes.POINTER(ctypes.c_long),ctypes.POINTER(ctypes.c_long),ctypes.POINTER(ctypes.c_long), \
						ctypes.POINTER(ctypes.c_float)]
	sfmLib.checkValidityOfLinkedPoint.restype = ctypes.c_long
		
	xo = array([x]).astype(int32); yo = array([y]).astype(int32)
	ro = array([r]).astype(int32); ao = array([angle]).astype(float32)
	
	flag = sfmLib.checkValidityOfLinkedPoint(ctypes.c_long(x),ctypes.c_long(y),ctypes.c_long(z), ctypes.c_long(r), \
						ctypes.c_long(nx),ctypes.c_long(ny),ctypes.c_long(nz),\
						im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
						xo.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), \
						yo.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), \
						ro.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), \
						ao.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	x = xo[0]; y = yo[0]; r = ro[0];						
	figure(1); clf();
	bw = 8 * r
	x0 = max(0,x - bw)
	y0 = max(0,y - bw)
	x1 = min(m,x + bw + 1)
	y1 = min(n,y + bw + 1)
	imshow(I[x0:x1,y0:y1],cm.gray); title('flag='+str(flag))
	for i in range(y-r,y+r+1):						
		for j in range(x-r,x+r+1):
			dd = sqrt((j-x)**2 + (i-y)**2)
			if dd >= r and dd < r+1:
				plot(i-y0,j-x0,'.r')
	rp = 5 * r
	for i in range(8):
		ang = 3.14156/8*i;
		xx1 = rp * cos(ang) 
		yy1 = rp * sin(ang)
		if abs(ang - ao[0]) < 1e-3:
			cl = 'r'
		else:
			cl = 'b'
		plot([y-y0+yy1,y-y0-yy1],[x-x0+xx1,x-x0-xx1],cl)			
	axis('off')
	print 'xo=',xo[0],' yo=',yo[0],' z=',z,' ro=',ro[0],' angle=',ao[0],' flag=',flag

# This function is for testing adjustZOfLinkedPoint. 	
# IMPORTANT: in libNeuTu.c, set adjustZOfLinkedPoint_DAG to 1 to make this work! 
#			Then need to recompile and reload libNeuTu.so
# 			A code for plotting the diagnosis is generated and written to temp3.py	
def testAdjustZOfLinkedPoint():

	iParSet = 3

	# parameters for the figure in the paper showing checkin validity of a SWC point. 
	if iParSet == 1:
		filename = 'testImages/OneChannelDH070613C2X100-20.tif'
		y = 621
		x = 317
		z = 76
		r = 8
		iBackgrounType = 1	# set to 0 if the image is bright background.
	elif iParSet == 2:	# other test points. 
		filename = 'testImages/OneChannelDH070613C2X100-40.tif'
		x = 191
		y = 799
		z = 72
		r = 59
		iBackgrounType = 1	# set to 0 if the image is bright background.			
	elif iParSet == 3: # bright background test points
		filename='/home/djin/NeuronReconstructionData/DH_10-9-13_100xredo/DH100913X100-47.tif'
		x = 486
		y = 1231
		z = 19
		r = 6
		iBackgrounType = 0	# set to 0 if the image is bright background.			
	
	im3d,imFlat,imRGB = loadStacks(filename,iBackgrounType);
	# normalize im3d
	mmin = im3d.min()
	mmax = im3d.max()
	im3d = (im3d - mmin)/(mmax - mmin)
	I = im3d[:,:,z].copy().astype(float32) 
	I = I -  I.min()
	sigma = 1	# sigma of the Gaussian for smoothing for taking the derivatives. 
	#imFlat = gaussian_filter(I,sigma)

	nx,ny,nz = im3d.shape
	im3d = im3d.astype(float32)
	
	sfmLib.adjustZOfLinkedPoint.argtypes = [ctypes.c_long,ctypes.c_long,ctypes.c_long,ctypes.c_long,\
						ctypes.c_long,ctypes.c_long,ctypes.c_long,\
						ctypes.POINTER(ctypes.c_float), \
						ctypes.POINTER(ctypes.c_long)]
	sfmLib.adjustZOfLinkedPoint.restype = None
		
	zo = array([z]).astype(int32)
	
	sfmLib.adjustZOfLinkedPoint(ctypes.c_long(x),ctypes.c_long(y),ctypes.c_long(z),ctypes.c_long(r),\
						ctypes.c_long(nx),ctypes.c_long(ny),ctypes.c_long(nz),\
						im3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), \
						zo.ctypes.data_as(ctypes.POINTER(ctypes.c_long)))
						
	figure(1); clf();
	bw = 5 * r
	x0 = max(0,x - bw)
	y0 = max(0,y - bw)
	x1 = min(nx,x + bw + 1)
	y1 = min(ny,y + bw + 1)
	imshow(I[x0:x1,y0:y1],cm.gray); title('flag='+str(flag))
	for i in range(y-r,y+r+1):						
		for j in range(x-r,x+r+1):
			dd = sqrt((j-x)**2 + (i-y)**2)
			if dd >= r and dd < r+1:
				plot(i-y0,j-x0,'.r')
	axis('off')
	print 'z=',z,' zo=',zo[0]
	import temp3; reload(temp3)

# This functino test the method of getting thrsholds from images. 	
def testThresholding():
	filename = 'testImages/OneChannelDH070613C2X100-40.tif.imFlat1.tif'
	I = imread(filename).astype(float32)
	m,n = I.shape
	ps = 100
	while 1:
		i1 = int((m-2*ps)*rand())
		j1 = int((m-2*ps)*rand())
		IP = I[i1:i1+ps,j1:j1+ps].copy()	
		IP = gaussian_filter(IP,sigma=1)
		figure(1); clf();
		subplot(221); imshow(IP,cm.gray); title('Image Patch')
		I2 = (IP.max() - IP)/(IP.max() - IP.min())
		mp,np = I2.shape
		I2L = I2.reshape(mp*np,1)
		subplot(222);
		hist(I2L,bins=50)
		title('Pixel value distribution (inverted, normalized)')
		bw = zeros((mp,np))
		threshold = median(I2L) + 2.0 * std(I2L)
		bw[where(I2 > threshold)] = 1
		subplot(223)
		imshow(bw,cm.gray)
		title('Threshold='+str(threshold))
		raw_input("Press any key...")
		
	
# This file tests line detector filtering. 
def testValleyDetector():
	sigma = 5.0
	nx = 40
	ny = 40
	
	# plot the second-derivative Gaussian filter. 
	F = zeros((nx,ny))
	s2 = sigma*sigma
	coef = 1.0/(2*pi*s2)
	for i in range(nx):
		for j in range(ny):
			x = i - nx/2
			y = j - ny/2
			F[i,j] = coef * (x*x*1.0/s2 - 1.0)/s2*exp(-(x*x+y*y)/(2.0*s2))

	figure(20,figsize=(5,3)); clf()
	rcParams['pdf.fonttype'] = 42	
	font = {'family' : 'Arial',
			'weight' : 'normal',
			'size'   : 10}
	font2 = {'family' : 'Arial',
			'weight' : 'normal',
			'size'   : 10}
	rec1 = [0.01,0.1,0.25,0.9]
	rec2 = [0.3,0.2,0.3,0.6]
	rec3 = [0.68,0.2,0.3,0.6]		
	ax = axes(rec1)
	imshow(F,cm.gray); axis('off')
	title(r'$\sigma$ = '+str(sigma),font)
	ylabel('$x$',font)
	xlabel('$y$',font)
	ax = axes(rec2)
	plot(F[:,ny/2]*10000)
	locator_params(axis = 'y', nbins = 4)
	yl = ax.get_ylim()
	text(0,yl[1]+(yl[1]-yl[0])*0.02,r'$\times 10^{-4}$',font2)
	x = linspace(0,nx,nx/10+1).astype(int)
	xticks(x,x-nx/2)
	tick_params(axis='both', which='major', labelsize=10, width=1)	
	xlabel('$x$',font)
	ylabel('$f$',font)
	title('$y=0$',font)
	ax = axes(rec3)
	plot(F[nx/2,:]*10000)
	y = linspace(0,ny,ny/10+1).astype(int)
	yl = ax.get_ylim()
	text(0,yl[1]+(yl[1]-yl[0])*0.02,r'$\times 10^{-4}$',font2)
	locator_params(axis = 'y', nbins = 4)
	xticks(y,y-ny/2)
	tick_params(axis='both', which='major', labelsize=10, width=1)	
	xlabel('$y$',font)
	ylabel('$f$',font)	
	title('$x=0$',font)
	show()
	savefig("Figs/Fig.LineDetector.pdf")
		
# Test stitching function. 
def testStitching():

	offsetDirection = 2	# 1, up, 2, down, 3, left, 4, right.
	filename1="/home/djin/NeuronReconstructionData/HP_15_02_02_C1/HP150202C1-26.tif"
	filename2="/home/djin/NeuronReconstructionData/HP_15_02_02_C1/HP150202C1-35.tif"
	overlapFract = 0.2
	
	imageType = 0	# bright field
	
	sfmLib.DirectionalPairwiseStitchC.argtypes = [ctypes.c_long,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_float,\
						ctypes.c_long, ctypes.POINTER(ctypes.c_float)]
	sfmLib.DirectionalPairwiseStitchC.restype = None
	
	ret = array([0,0,0,0]).astype(float32)
		
	sfmLib.DirectionalPairwiseStitchC(ctypes.c_long(imageType),ctypes.c_char_p(filename1),ctypes.c_char_p(filename2),\
						ctypes.c_float(overlapFract),ctypes.c_long(offsetDirection),\
						ret.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
	dx,dy,dz,maxcorr = ret
	print 'dx=',dx,' dy=',dy,' dz=',dz,' maxcorr=',maxcorr
	
# Explore enhancing brown color in converting grayscale. 
import scipy.optimize as optimization
def func(params, xdata, ydata):
    return (ydata - dot(xdata, params))
	
def enhanceBrown():
	brownColors = [[ 0.64  ,  0.58  ,  0.5   ], \
       [ 0.5   ,  0.1647,  0.1647], \
       [ 0.86  ,  0.16  ,  0.16  ], \
       [ 0.53  ,  0.26  ,  0.12  ], \
       [ 0.8706,  0.7216,  0.5294], \
       [ 0.54  ,  0.21  ,  0.06  ], \
       [ 0.54  ,  0.2   ,  0.14  ], \
       [ 0.8235,  0.4118,  0.1176], \
       [ 0.45  ,  0.24  ,  0.1   ], \
       [ 1.    ,  0.49  ,  0.25  ], \
       [ 1.    ,  0.34  ,  0.13  ], \
       [ 0.78  ,  0.47  ,  0.15  ], \
       [ 1.    ,  0.24  ,  0.05  ], \
       [ 0.9412,  0.902 ,  0.549 ], \
       [ 0.7412,  0.7176,  0.4196], \
       [ 0.9608,  0.9608,  0.8627], \
       [ 0.8039,  0.5216,  0.2471], \
       [ 0.7373,  0.5608,  0.5608], \
       [ 0.78  ,  0.38  ,  0.08  ], \
       [ 0.45  ,  0.29  ,  0.07  ], \
       [ 0.37  ,  0.15  ,  0.07  ], \
       [ 0.6275,  0.3216,  0.1765], \
       [ 0.5451,  0.2706,  0.0745], \
       [ 0.9569,  0.6431,  0.3765], \
       [ 0.8235,  0.7059,  0.549 ], \
       [ 0.37  ,  0.15  ,  0.02  ]]
	print 'Going through the named brown colors...'
	img = zeros((100,100,3)).astype(uint8)
	for r,g,b in brownColors:
		img[:,:,0] = r*255
		img[:,:,1] = g*255
		img[:,:,2] = b*255
		imshow(img)
		raw_input("Press any key...")   
	
	# least sqaure fit this to a model b = a * r + b * g 
	brownColors = array(brownColors)
	xdata = transpose(array([brownColors[:,0],brownColors[:,1]]))
	ydata = brownColors[:,2]
	x0 = array([0,0])
	params, dd = optimization.leastsq(func, x0, args=(xdata, ydata))

	# now try colors.
	print 'Sampling rgb from the plane...' 
	while 1:
		r = rand()
		g = rand()
		if g > r:		# for brown, red color should be larger than green color. 
			continue
		b = params[0]*r + params[1]*g
		if b < 0 or b > 1:
			continue
		img[:,:,0] = r*255
		img[:,:,1] = g*255
		img[:,:,2] = b*255
		imshow(img)
		raw_input("Press any key...")   
	
	print 'Testing how brown random colors are...'
	while 1:
		r = rand()
		g = rand()
		b = rand()
		dd1 = abs(params[0]*r + params[1]* g - b)/sqrt(r*r+g*g+b*b)
		dd2 = (1+tanh(-(r-b))/(3.1415/2))/2
		img[:,:,0] = r*255
		img[:,:,1] = g*255
		img[:,:,2] = b*255
		imshow(img)
		title('dd1='+str(dd1)+' dd2='+str(dd2))
		raw_input("Press any key...")   
				
	# load a color image. 
	filename = '/home/dezhe/NeuronReconstructionData/DH_7-8-13_100x/DH070813100x15.Proj.tif'
	imFlat = imread(filename)
	nx,ny,nz = imFlat.shape
	imGray = zeros((nx,ny))
	for i in range(nx):
		for j in range(ny):
			r,g,b=imFlat[i,j,:]
			r /= 255.0
			g /= 255.0
			b /= 255.0
			dd1 = abs(params[0]*r + params[1]* g - b)/sqrt(r*r+b*b+g*g)
			imGray[i,j] = (r+b+g)/3 #* exp(-dd1/0.1)
	imshow(imGray,cm.gray)
			 
def testProcessImage():
	
	#filenameBase = 'testImages/OneChannelDH070613C2X100-20'
	#imageType = 1

	filenameBase = '/home/dezhe/NeuronReconstructionData/DH_7-6-13-2_100x/DH070613C2X100-33'
	imageType = 0
		
	sfmLib.processImage.argtypes = [ctypes.c_char_p,ctypes.c_long]
	sfmLib.processImage.restype = None
	sfmLib.processImage(ctypes.c_char_p(filenameBase),ctypes.c_long(imageType))
		
