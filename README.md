# hpc-final

To run the mpi version of the program, fftw_mpi need to be installed. If this package is not originally installed in your computer, you can download the package from http://fftw.org/ and follow the installation guide of http://micro.stanford.edu/wiki/Install_FFTW3 

After that, when you compile the program, you need to include the 'include' and 'lib' directories as example below:(and thus the make file need to be modified) 
mpicc -I/home/qx344/qx344/fftwmpi/install/include  main_mpi.c  -lfftw3_mpi -lfftw3 -lm -g -fopenmp  -L/home/qx344/qx344/fftwmpi/install/lib -o test3

To run the program, you need to include the path to 'lib' in your shared library path LD_LIBRARY_PATH as below:
LD_LIBRARY_PATH = $LD_LIBRARY_PATH:/home/qx344/qx344/fftwmpi/install/lib