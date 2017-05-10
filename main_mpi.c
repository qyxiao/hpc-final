#include<stdio.h>
#include<stdlib.h>
#include <stddef.h> 
#include <stdint.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include <fftw3-mpi.h>
#define pi 3.1415926


void exactVorticity(fftw_complex *vorticity,double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, ptrdiff_t preIndex, ptrdiff_t xLength, ptrdiff_t yLength, double t){
    //fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    ptrdiff_t i,j;
    for(i=0;i<xLength;i++){
    	for(j=0;j<yLength;j++){
    		vorticity[i*yLength+j] = 8*pi/L*exp(-8*pi*pi*miu*t/(L*L))*cos((xMatrix[preIndex+i]-v0*t)*2*pi/L)*cos((yMatrix[j]-v0*t)*2*pi/L);
    	}
    }
}



fftw_complex* CNfun(fftw_complex* vortiCom, fftw_complex* NfluxCom, double miu, fftw_complex* KX, fftw_complex* KY, int xLength, int yLength, double dt){
    fftw_complex* ans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    double factor = -1;
    double timeInv = 1/dt;
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            double temp = (KX[i]*KX[i]+KY[j]*KY[j])*miu*0.5;
            ans[i*yLength+j] = ((timeInv-temp)*vortiCom[i*yLength+j]+factor*NfluxCom[i*yLength+j])/(timeInv+temp);
        }
    }
    return ans;
}



fftw_complex* CNfun_AB(fftw_complex* vortiCom, fftw_complex* NfluxCom, fftw_complex* NfluxComOld, double miu, fftw_complex* KX, fftw_complex* KY, int xLength, int yLength, double dt){
    fftw_complex* ans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    double factor = -1;
    double timeInv = 1/dt;
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            double temp = (KX[i]*KX[i]+KY[j]*KY[j])*miu*0.5;
            ans[i*yLength+j] = ((timeInv-temp)*vortiCom[i*yLength+j]+factor*(1.5*NfluxCom[i*yLength+j]-0.5*NfluxComOld[i*yLength+j]))/(timeInv+temp);
        }
    }
    return ans;
}






int main(int argc, char* argv[])
{
    int iter, maxIter, mpirank;
	int Npoints = 16;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j, xLength=Npoints, yLength=Npoints, preIndex;

	double t = 0, miu = 0.05, v0 = 1, L = 1, dx = L/xLength, dy = L/yLength, dt = 1.0/64;
    double fftDefactor = 1.0/(xLength*yLength);

    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	

    fftw_complex *xMatrix, *yMatrix, *KX, *KY;
    fftw_complex *vorticity, *Uvel, *Vvel, *vortiX, *vortiY, *Nflux, *stream; 
    fftw_complex *vorticityCom, *UvelCom, *VvelCom, *vortiXCom, *vortiYCom, *NfluxCom, *NfluxComOld,*vorticityComNew;
	fftw_plan planVor, planInvVor, planU, planInvU, planV, planInvV, planInvVortiX, planInvVortiY, planN;
    
	   
	xMatrix = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength);
    KX = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength);
    yMatrix = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*yLength);
    KY = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*yLength);
    
    for(i=0; i<(xLength+1)/2; i++){
        KX[i]=i*2*pi/L;
        xMatrix[i]=i*dx;
    }
    for(i=(xLength+1)/2; i<xLength; i++){
        KX[i]=(i-xLength)*2*pi/L;
        xMatrix[i]=i*dx;
    }
    for(i=0; i<(yLength+1)/2; i++){
        KY[i]=i*2*pi/L;
        yMatrix[i]=i*dy;
    }
    for(i=(yLength+1)/2; i<yLength; i++){
        KY[i]=(i-yLength)*2*pi/L;
        yMatrix[i]=i*dy;
    } 

    iter = 0; maxIter = 0;
    
    alloc_local = fftw_mpi_local_size_2d(xLength, yLength, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    preIndex = mpirank*local_n0;

    vorticity = fftw_alloc_complex(alloc_local);

    exactVorticity(vorticity, L, miu, xMatrix, yMatrix, v0, preIndex, local_n0, yLength, t);
    stream = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*yLength);
    Uvel  = fftw_alloc_complex(alloc_local);
    Vvel  = fftw_alloc_complex(alloc_local);
    vortiX = fftw_alloc_complex(alloc_local);
    vortiY = fftw_alloc_complex(alloc_local);
    Nflux = fftw_alloc_complex(alloc_local);

    vorticityCom = fftw_alloc_complex(alloc_local);
    vorticityComNew = fftw_alloc_complex(alloc_local);
    UvelCom = fftw_alloc_complex(alloc_local); 
    VvelCom = fftw_alloc_complex(alloc_local);
    vortiXCom = fftw_alloc_complex(alloc_local);
    vortiYCom = fftw_alloc_complex(alloc_local);
    NfluxCom = fftw_alloc_complex(alloc_local);
    
    planVor = fftw_mpi_plan_dft_2d(xLength, yLength, vorticity, vorticityCom, MPI_COMM_WORLD,FFTW_FORWARD, FFTW_ESTIMATE); 

    planU = fftw_plan_dft_2d(xLength, yLength, Uvel, UvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    planInvU = fftw_plan_dft_2d(xLength, yLength, UvelCom, Uvel, FFTW_BACKWARD, FFTW_ESTIMATE);    
    planV = fftw_plan_dft_2d(xLength, yLength, Vvel, VvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    planInvV = fftw_plan_dft_2d(xLength, yLength, VvelCom, Vvel, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(planVor);

    
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputMPI%d.txt",mpirank);
    fd = fopen(filename,"w+");
    fprintf(fd,"timestep: %f\n",dt);
    for(i = 0; i < local_n0*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%d %11.7f %11.7f\n", i, creal(vorticityCom[i]), cimag(vorticityCom[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
        
    fclose(fd);
    
 //    fftw_free(vorticity); fftw_free(Uvel); fftw_free(Vvel); fftw_free(vortiX); fftw_free(vortiY); fftw_free(Nflux); 
 //    fftw_free(xMatrix); fftw_free(yMatrix); fftw_free(KX); fftw_free(KY); fftw_free(stream);

	// fftw_destroy_plan(planVor); fftw_destroy_plan(planInvVor); fftw_destroy_plan(planU); fftw_destroy_plan(planV);
 //    fftw_destroy_plan(planInvVortiX); fftw_destroy_plan(planInvVortiY); 

	// fftw_free(vorticityCom); fftw_free(UvelCom); fftw_free(VvelCom); fftw_free(vortiXCom); fftw_free(vortiYCom);
 //    fftw_free(NfluxCom);    
    
    MPI_Finalize();  	
	return 0;
}