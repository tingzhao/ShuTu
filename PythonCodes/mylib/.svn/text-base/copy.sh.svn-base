#!/bin/sh

# Copy latest mylib code into vendor branch working copy
# Update mylib code before running this script.
# Check in vendor branch after running this script and testing the result.

# Set these variables to the directory locations on your system
MYLIB_SRC_DIR=~/svn/mylib
VENDOR_MYLIB_DEST_DIR=~/svn/vendor_myers_mylib

(cd $MYLIB_SRC_DIR; tar cf - --exclude=.svn .) | (cd $VENDOR_MYLIB_DEST_DIR; tar xvf -)

