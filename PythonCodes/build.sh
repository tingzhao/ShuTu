#!/bin/bash
if [ `uname` == 'Darwin' ]; then
  sh build_mac.sh
  exit
fi

scriptDir=$(dirname `readlink -f "$0"`)


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
    easy_install $package
    set +e
  done
fi

echo "Checking building tools ..."
if [ `uname` == 'Darwin' ]; then
  command -v gcc 2> /dev/null
  if [ $? -ne 0 ]
  then
    #echo "ERROR: gcc is missing. Please install build-essential"
    #exit 1
    set -e
    apt-get install build-essential
    set +e
  fi
fi

command -v mpirun mpicc 2> /dev/null
if [ $? -ne 0 ]
then
  #echo "ERROR: OpenMPI is missing. Please install openmpi-bin and libopenmpi-dev."
  #exit 1
  set -e
  pip install openmpi-bin
  pip install libopenmpi-dev
  set +e
fi


echo "Building mylib ..."
cd $scriptDir/mylib
make

cd $scriptDir
make

