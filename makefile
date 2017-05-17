
EXECS= main main_omp main_mpi main_mpi_omp test_time


all: ${EXECS}

main: main.c
	gcc $^  -lfftw3 -lm -lrt -o test1


main_omp: main_omp.c
	gcc $^  -lfftw3 -lm -lrt -fopenmp -o test2


main_mpi: main_mpi.c
	mpicc -I/home/qx344/qx344/fftwmpi/install/include  $^  -lfftw3_mpi -lfftw3 -lm  -L/home/qx344/qx344/fftwmpi/install/lib -o test3

main_mpi_omp: main_mpi.c
	mpicc -I/home/qx344/qx344/fftwmpi/install/include  $^  -lfftw3_mpi -lfftw3 -lm -fopenmp -L/home/qx344/qx344/fftwmpi/install/lib -o test4


test_time: timetest.c
    mpicc -I/home/qx344/qx344/fftwmpi/install/include  $^  -lfftw3_mpi -lfftw3 -lm -fopenmp -L/home/qx344/qx344/fftwmpi/install/lib -o timetest

clean:
	rm -f ${EXECS}