# Generating shared libraries
# (Of course, you can just shove all things in the header file and not do any of below.)
#
#  #How-To:
#  
#  let,
#      icc     = /opt/intel/Compiler/13.1/composer_xe_2013.5.192/bin/intel64/icc 
#      cflags  = -std=c++11 -w1 -xHOST -O3 -ipo -mkl
#  
#  # for gcc;
#      gcc     = /opt/rh/devtoolset-2/root/usr/bin/gcc
#      cflags  = -std=c++11 -Wall -O3
#      # no -fwhole-program option!
#  
#  and let
#      myfile.cpp, myfile.h
#  be in this directory.
#  
#  #Below works for g++ also! Just be careful of the flags.
#  1. make .o file with -fpic option:
#      $icc $cflags -fpic -c myfile.cpp
#  
#  2. make .so file from .o file (so file need to(?) have the prefix lib):
#      $icc $cflags -shared -o libmyfile.so myfile.o
#  
#  3. modify the Makefile:
#      LIBS = -lmyfile
#  
#  4. remove .o file
#  
#  ex)
#      icc -std=c++11 -w1 -xHOST -O3 -ipo -mkl -fpic -c maths.cpp
#      icc -std=c++11 -w1 -xHOST -O3 -ipo -mkl -shared -o libmaths.so maths.o

gcc="/opt/rh/devtoolset-2/root/usr/bin/gcc"
ops="-std=c++11 -Wall -g -O3"

# Add additional file names to the list below without their extensions
#listFileName="consts maths"
listFileName="cnumber"

for fn in $listFileName; do
    echo "Compiling ${fn}.cpp to get ${fn}.o"
    ${gcc} ${ops} -fpic -c ${fn}.cpp

    echo "Generating the shared object file lib${fn}.so"
    ${gcc} ${ops} -shared -o lib${fn}.so ${fn}.o

done

echo "Erasing object files (*.o)"
for fn in $listFileName; do
    rm -f ${fn}.o
done
