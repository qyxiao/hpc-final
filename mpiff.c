#include <stdio.h> 
#include <stddef.h> 
#include <stdint.h>
#include <complex.h>
#include <fftw3-mpi.h>

int main(int argc, char **argv)
{

    int Npoints = 16, mpirank;
    ptrdiff_t N0 = 4, N1 = 4;
    fftw_plan plan;
    fftw_complex *data;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;

    MPI_Init(&argc, &argv);
    fftw_mpi_init();
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);

    /* create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                                FFTW_FORWARD, FFTW_ESTIMATE);    

    /* initialize data to some function my_function(x,y) */
    fftw_complex iter = mpirank*N1;
    for (i = 0; i < local_n0; ++i) 
        for (j = 0; j < N1; ++j){
            data[i*N1 + j] = iter;
            iter++;
        }
           

    /* compute transforms, in-place, as many times as desired */
    fftw_execute(plan);
    

    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputTest%d.txt",mpirank);
    fd = fopen(filename,"w+");
    iter = 0;
    for(i = 0; i < local_n0*N1; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%f %11.7f %11.7f\n", creal(iter), creal(data[i]), cimag(data[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
        iter = iter+1;
    }

    fftw_destroy_plan(plan);

    MPI_Finalize();
}