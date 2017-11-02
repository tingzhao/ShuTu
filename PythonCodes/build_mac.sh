#!/bin/bash

scriptDir=$(cd $(dirname "$0") && pwd -P)

condaDir=$scriptDir/miniconda
if [ ! -d $condaDir ]
then
  curl -X GET https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh > Miniconda-latest-MacOSX-x86_64.sh
  bash Miniconda-latest-MacOSX-x86_64.sh -b -p $condaDir
fi
source $condaDir/bin/activate root
 
#Checking python packages
echo "Checking python packages ..."
while read package
do
  python -c "import $package" 2> /dev/null
  if [ $? -ne 0 ]
  then
    missing="$missing $package"
  fi
done < $scriptDir/requirements.txt

if [ ! -z "$missing" ]
then
  echo "WARNING: The following python packages are missing:"
  for package in $missing
  do
    echo "*   $package"
    #echo "Please use easy_install or pip to install them and then try agin."
    echo "Installing now ..."
    set -e
    conda install $package -y
    set +e
  done
fi

#command -v mpirun mpicc 2> /dev/null
#if [ $? -ne 0 ]
#then
  #echo "ERROR: OpenMPI is missing. Please install openmpi-bin and libopenmpi-dev."
  #exit 1
set -e
conda install -c mpi4py openmpi=2.0.2 -y
set +e
#fi


echo "Building mylib ..."
cd $scriptDir/mylib
make

cd $scriptDir
make

