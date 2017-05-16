#include <stdio.h> 
#include<stdio.h>
#include<stdlib.h>
#include <stddef.h> 
#include <stdint.h>
#include "util.h"
#include <complex.h>
#include <fftw3-mpi.h>

int main(int argc, char **argv)
{

    int iteration=0, maxIter=10, mpirank;
    int Npoints = 4;
    if(argc>=2){sscanf(argv[1], "%d", &Npoints);}  
    if(argc>=3){sscanf(argv[2], "%d", &maxIter);} 
    fftw_plan plan;
    fftw_complex *data, *datahat;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
    

    MPI_Init(&argc, &argv);
    fftw_mpi_init();
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(Npoints, Npoints, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);
    datahat = fftw_alloc_complex(alloc_local);
    /* create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(Npoints, Npoints, data, datahat, MPI_COMM_WORLD,
                                FFTW_FORWARD, FFTW_ESTIMATE);     
    /* initialize data to some function my_function(x,y) */
    fftw_complex iter = 0;

    iter = mpirank*Npoints;
    for (i = 0; i < local_n0; ++i){
        for (j = 0; j < Npoints; ++j){
            data[i*Npoints + j] = iter;
            iter ++;
        }
    } 

    while(iteration<maxIter){
        fftw_execute(plan);
        iteration++;
    }
    

    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputTest%d.txt",mpirank);
    fd = fopen(filename,"w+");
    iter = 0;
    for(i = 0; i < local_n0*Npoints; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%td %11.7f %11.7f\n", i, creal(datahat[i]), cimag(datahat[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
        iter = iter+1;
    }

    fftw_destroy_plan(plan);

    MPI_Finalize();
    return 0;
}