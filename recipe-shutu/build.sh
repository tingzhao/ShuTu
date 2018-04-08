if [ $(uname) == 'Darwin' ]; then
    CC=/usr/bin/cc
    CXX=/usr/bin/clang
fi

if [ $(uname) == 'Darwin' ]; then
    QMAKE_SPEC_PATH=${PREFIX}/mkspecs/macx-clang
else
    QMAKE_SPEC_PATH=${PREFIX}/mkspecs/linux-g++-64
fi

export CONDA_ENV=${PREFIX}

app_name=ShuTu

cd ShuTu_gui
build_dir=neurolabi/build
if [ "$app_name" == 'ShuTu_d' ] 
then
  qtlib_dir=${PREFIX}/lib
  cd neurolabi/shell
  ./fixqtdebug Qt5 $qtlib_dir
  build_flag='-c debug'
  build_dir=neurolabi/build_debug
  cd ../../
fi

edition=biocytin

bash -x -e build.sh ${PREFIX}/bin/qmake ${QMAKE_SPEC_PATH} -e $edition $build_flag

# Install to conda environment
if [ $(uname) == 'Darwin' ]; then
    mv ${build_dir}/${app_name}.app ${PREFIX}/bin/
else
    mv ${build_dir}/${app_name} ${PREFIX}/bin/
    mv ${build_dir}/config.xml ${PREFIX}/bin/
    mv ${build_dir}/doc ${PREFIX}/bin/
fi
