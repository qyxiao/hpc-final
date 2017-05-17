#include<stdio.h>
#include<stdlib.h>
#include <stddef.h> 
#include <stdint.h>
#include "util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include<math.h>
#include<complex.h>
#include<time.h>
#include <fftw3-mpi.h>
#define pi 3.1415926


void exactVorticity(fftw_complex *vorticity,double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, ptrdiff_t preIndex, ptrdiff_t xLength, ptrdiff_t yLength, double t){
    //fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    ptrdiff_t i,j;
    #pragma omp parallel for default(none) shared(L,miu,xLength,yLength,vorticity,xMatrix,yMatrix,v0,t,preIndex) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
    	for(j=0;j<yLength;j++){
    		vorticity[i*yLength+j] = 8*pi/L*exp(-8*pi*pi*miu*t/(L*L))*cos((xMatrix[preIndex+i]-v0*t)*2*pi/L)*cos((yMatrix[j]-v0*t)*2*pi/L);
    	}
    }
}



void CNfun(fftw_complex* vortiCom, fftw_complex* NfluxCom, double miu, fftw_complex* KX, fftw_complex* KY, ptrdiff_t preIndex, ptrdiff_t xLength, ptrdiff_t yLength, double dt){
    //fftw_complex* ans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    ptrdiff_t i,j;
    double factor = -1;
    double timeInv = 1/dt;
    #pragma omp parallel for default(none) shared(vortiCom,NfluxCom,miu,preIndex,xLength,yLength,KX,KY,dt,factor,timeInv) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            double temp = (KX[preIndex+i]*KX[preIndex+i]+KY[j]*KY[j])*miu*0.5;
            fftw_complex vortiTemp = vortiCom[i*yLength+j];
            vortiCom[i*yLength+j] = vortiTemp; //((timeInv-temp)*vortiTemp+factor*NfluxCom[i*yLength+j])/(timeInv+temp);
        }
    }
}



void CNfun_AB(fftw_complex* vortiCom, fftw_complex* NfluxCom, fftw_complex* NfluxComOld, double miu, fftw_complex* KX, fftw_complex* KY, ptrdiff_t preIndex, ptrdiff_t xLength, ptrdiff_t yLength, double dt){
    //fftw_complex* ans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    ptrdiff_t i,j;
    double factor = -1;
    double timeInv = 1/dt;
    #pragma omp parallel for default(none) shared(vortiCom,NfluxCom,NfluxComOld,miu,preIndex,xLength,yLength,KX,KY,dt,factor,timeInv) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            double temp = (KX[preIndex+i]*KX[preIndex+i]+KY[j]*KY[j])*miu*0.5;
            vortiCom[i*yLength+j] = ((timeInv-temp)*vortiCom[i*yLength+j]+factor*(1.5*NfluxCom[i*yLength+j]-0.5*NfluxComOld[i*yLength+j]))/(timeInv+temp);
        }
    }
}






int main(int argc, char* argv[])
{
    int iter=0, maxIter=1, mpirank;
	int Npoints = 16;
    if(argc>=2){sscanf(argv[1], "%d", &Npoints);}  
    if(argc>=3){sscanf(argv[2], "%d", &maxIter);} 

    ptrdiff_t alloc_local, local_n0, local_0_start, i, j, xLength=Npoints, yLength=Npoints, preIndex;
    printf("number %td ", xLength);
	double t = 0, miu = 0.05, v0 = 1, L = 1, dx = L/xLength, dy = L/yLength, dt = 4.0*dx*dy; //1.0/64;
    double fftDefactor = 1.0/(xLength*yLength);

    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	
    timestamp_type time1, time2;   // time test
    get_timestamp(&time1);

    fftw_complex *xMatrix, *yMatrix, *KX, *KY;
    fftw_complex *vorticity, *Uvel, *Vvel, *vortiX, *vortiY, *Nflux, *stream; 
    fftw_complex *vorticityCom, *UvelCom, *VvelCom, *vortiXCom, *vortiYCom, *NfluxCom, *NfluxComOld,*vorticityComNew;
	fftw_plan planVor, planInvVor, planInvU, planInvV, planInvVortiX, planInvVortiY, planN;
    
	   
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
    planInvU = fftw_mpi_plan_dft_2d(xLength, yLength, UvelCom, Uvel, MPI_COMM_WORLD,FFTW_BACKWARD, FFTW_ESTIMATE);
    planInvV = fftw_mpi_plan_dft_2d(xLength, yLength, VvelCom, Vvel, MPI_COMM_WORLD,FFTW_BACKWARD, FFTW_ESTIMATE);
    planInvVortiX = fftw_mpi_plan_dft_2d(xLength, yLength, vortiXCom, vortiX, MPI_COMM_WORLD,FFTW_BACKWARD, FFTW_ESTIMATE);
    planInvVortiY = fftw_mpi_plan_dft_2d(xLength, yLength, vortiYCom, vortiY, MPI_COMM_WORLD,FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(planVor);

    while(iter<maxIter){
 
        #pragma omp parallel for default(none) shared(preIndex,local_n0,yLength,stream,KX,KY,vorticityCom,UvelCom,VvelCom,vortiXCom,vortiYCom) private(i,j) collapse(2)
        for(i=0;i<local_n0;i++){
            for(j=0;j<yLength;j++){
                stream[i*yLength+j] = vorticityCom[i*yLength+j]/(KX[preIndex+i]*KX[preIndex+i]+KY[j]*KY[j]);
                UvelCom[i*yLength+j] = stream[i*yLength+j]*I*KY[j];
                VvelCom[i*yLength+j] = -1*stream[i*yLength+j]*I*KX[preIndex+i];
                vortiXCom[i*yLength+j] = I*KX[preIndex+i]* vorticityCom[i*yLength+j];
                vortiYCom[i*yLength+j] = I*KY[j]* vorticityCom[i*yLength+j];
            }
        }
        if(mpirank==0){
            stream[0] = vorticityCom[0];
            UvelCom[0] = stream[0]*I*KY[0];
            VvelCom[0] = -1*stream[0]*I*KX[0];
        }


        fftw_execute(planInvU); //Execution of IFFT
        fftw_execute(planInvV); //Execution of IFFT
        fftw_execute(planInvVortiX);
        fftw_execute(planInvVortiY);

        #pragma omp parallel for default(none) shared(local_n0,yLength,Nflux,vortiX,Uvel,fftDefactor,vortiY,Vvel)  private(i,j) collapse(2)
        for(i=0;i<local_n0;i++){
            for(j=0;j<yLength;j++){
                Nflux[i*yLength+j] = vortiX[i*yLength+j]*(Uvel[i*yLength+j]*fftDefactor+1) + vortiY[i*yLength+j]*(Vvel[i*yLength+j]*fftDefactor+1) ;
                Nflux[i*yLength+j] = Nflux[i*yLength+j]*fftDefactor;
            }
        }

        planN = fftw_mpi_plan_dft_2d(xLength, yLength, Nflux, NfluxCom, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); 
        fftw_execute(planN);
        fftw_destroy_plan(planN); 
        if(1){
            CNfun(vorticityCom,NfluxCom,miu,KX,KY,preIndex,local_n0,yLength,dt);
        }else{
            CNfun_AB(vorticityCom,NfluxCom,NfluxComOld,miu,KX,KY,preIndex,local_n0,yLength,dt);
        }

        NfluxComOld = NfluxCom;
        NfluxCom = fftw_alloc_complex(alloc_local);

        //vorticityCom = vorticityComNew;
        iter++;

    }
    planInvVor = fftw_mpi_plan_dft_2d(xLength, yLength, vorticityCom, vorticity, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); 
    fftw_execute(planInvVor);
    
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputMPI%d.txt",mpirank);
    fd = fopen(filename,"w+");
    //fprintf(fd,"timestep: %f\n",dt);
    for(i = 0; i < local_n0*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%11.7f %11.7f\n", fftDefactor*creal(vorticity[i]), fftDefactor*cimag(vorticity[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
    

    fclose(fd);
 //    fftw_free(vorticity); fftw_free(stream); //fftw_free(vorticityCom); 
 //    // fftw_free(Uvel); fftw_free(Vvel); fftw_free(vortiX); fftw_free(vortiY); fftw_free(Nflux); 
 //    // fftw_free(xMatrix); fftw_free(yMatrix); fftw_free(KX); fftw_free(KY); 

	// fftw_destroy_plan(planVor); fftw_destroy_plan(planInvU); fftw_destroy_plan(planInvV); 
 //    fftw_destroy_plan(planInvVortiX); fftw_destroy_plan(planInvVortiY); 

	// fftw_free(UvelCom); fftw_free(VvelCom); fftw_free(vortiXCom); fftw_free(vortiYCom);
 //    fftw_free(NfluxCom);    

    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    if (0 == mpirank) {
        printf("Time elapsed is %f seconds.\n", elapsed);
    }

    MPI_Finalize();  	
	return 0;
}