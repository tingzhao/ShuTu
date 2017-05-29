# This python script is for creating tiff stacks from the images acquired using  
# Zeiss Zeiss Axio Imager with AxioCam and Zen blue software.
# Assume filename convention for the images is like slides-09_z035m30.tif, i.e. commonname_z(depth)m(tilenumber).tif
# Usage
#
# 	python createTiffStacksZeiss dirData filenameCommon
#
#	dirData 	data directory. The images can be in a subdirectory.
#
# The tiff stacks will be created in dirData. 

import sys, os

# check number of inputs. 
if (len(sys.argv) < 3 or len(sys.argv) > 3):
	print "This script requires two arguments. Usage: python createTiffStacksZeiss dirData filenameCommon"
	sys.exit(0)

dataDir = sys.argv[1]
filenameCommon = sys.argv[2]

print 'The data directory is ', dataDir
print 'filenameCommon is ',filenameCommon

# change directory to dataDir
currentDir = os.getcwd()
os.chdir(dataDir)

import ShuTu
ShuTu.filenameCommon = filenameCommon

ShuTu.getFilenameCommonSlice()
ShuTu.createTifStacksFromSlices()
# process all images. 
ShuTu.processAllImages()

# get back to the original directory
os.chdir(currentDir)




	
