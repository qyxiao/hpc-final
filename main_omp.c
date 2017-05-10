#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "util.h"
#include<complex.h>
#include<time.h>
#include<fftw3.h>
#define pi 3.1415926


fftw_complex* exactVorticity(double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, int xLength, int yLength, double t){
    fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    #pragma omp parallel for default(none) shared(L,miu,xLength,yLength,ans,xMatrix,yMatrix,v0,t) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
    	for(j=0;j<yLength;j++){
    		ans[i*yLength+j] = 8*pi/L*exp(-8*pi*pi*miu*t/(L*L))*cos((xMatrix[i]-v0*t)*2*pi/L)*cos((yMatrix[j]-v0*t)*2*pi/L);
    	}
    }
    return ans;
}


fftw_complex* exactU (double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, int xLength, int yLength, double t){
    fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    #pragma omp parallel for default(none) shared(L,miu,xLength,yLength,ans,xMatrix,yMatrix,v0,t) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            ans[i*yLength+j] = v0 - 2*exp(-8*pi*pi*miu*t/(L*L))*cos((xMatrix[i]-v0*t)*2*pi/L)*sin((yMatrix[j]-v0*t)*2*pi/L);
        }
    }
    return ans;
} 

fftw_complex* exactV (double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, int xLength, int yLength, double t){
    fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    #pragma omp parallel for default(none) shared(L,miu,xLength,yLength,ans,xMatrix,yMatrix,v0,t) private(i,j) collapse(2)
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            ans[i*yLength+j] = v0 + 2*exp(-8*pi*pi*miu*t/(L*L))*sin((xMatrix[i]-v0*t)*2*pi/L)*cos((yMatrix[j]-v0*t)*2*pi/L);
        }
    }
    return ans;
} 

fftw_complex* CNfun(fftw_complex* vortiCom, fftw_complex* NfluxCom, double miu, fftw_complex* KX, fftw_complex* KY, int xLength, int yLength, double dt){
    fftw_complex* ans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
    double factor = -1;
    double timeInv = 1/dt;
    #pragma omp parallel for default(none) shared(vortiCom,NfluxCom,miu,xLength,yLength,ans,KX,KY,dt,factor,timeInv) private(i,j) collapse(2)
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
    #pragma omp parallel for default(none) shared(vortiCom,NfluxCom,NfluxComOld,miu,xLength,yLength,ans,KX,KY,dt,factor,timeInv) private(i,j) collapse(2)
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
    int i,j, iter, maxIter;
	int Npoints = 32;
	int xLength = Npoints, yLength = Npoints;
	double t = 0, miu = 0.05, v0 = 1, L = 1, dx = L/xLength, dy = L/yLength, dt = 4.0*dx*dy;
    double fftDefactor = 1.0/(xLength*yLength);
	fftw_complex *vorticity, *Uvel, *Vvel, *vortiX, *vortiY, *Nflux, *stream; 
    fftw_complex *xMatrix, *yMatrix, *KX, *KY;
    fftw_complex *vorticityCom, *UvelCom, *VvelCom, *vortiXCom, *vortiYCom, *NfluxCom, *NfluxComOld,*vorticityComNew;
	fftw_plan planVor, planInvVor, planU, planInvU, planV, planInvV, planInvVortiX, planInvVortiY, planN;

    timestamp_type time1, time2;   // time test
    get_timestamp(&time1);
	   
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

    iter = 0; maxIter = 3;
    vorticity = exactVorticity(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
    stream = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
	// Uvel = exactU(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
 //    Vvel = exactV(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
    Uvel = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    Vvel = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    vortiX = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    vortiY = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    Nflux = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    
    vorticityCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength); 
    vorticityComNew = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength); 
    UvelCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength); 
    VvelCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength); 
    vortiXCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength); 
    vortiYCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    NfluxCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);

    // planVor = fftw_plan_dft_2d(xLength, yLength, vorticity, vorticityCom, FFTW_FORWARD, FFTW_ESTIMATE);
    // planInvVor = fftw_plan_dft_2d(xLength, yLength, vorticity, vorticityCom, FFTW_BACKWARD, FFTW_ESTIMATE);
    planVor = fftw_plan_dft_2d(xLength, yLength, vorticity, vorticityCom, FFTW_FORWARD, FFTW_ESTIMATE);
    //planInvVor = fftw_plan_dft_2d(xLength, yLength, vorticityComNew, vorticity, FFTW_BACKWARD, FFTW_ESTIMATE);
    planU = fftw_plan_dft_2d(xLength, yLength, Uvel, UvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    planInvU = fftw_plan_dft_2d(xLength, yLength, UvelCom, Uvel, FFTW_BACKWARD, FFTW_ESTIMATE);    
    planV = fftw_plan_dft_2d(xLength, yLength, Vvel, VvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    planInvV = fftw_plan_dft_2d(xLength, yLength, VvelCom, Vvel, FFTW_BACKWARD, FFTW_ESTIMATE);

    planInvVortiX = fftw_plan_dft_2d(xLength, yLength, vortiXCom, vortiX, FFTW_BACKWARD, FFTW_ESTIMATE);
    planInvVortiY = fftw_plan_dft_2d(xLength, yLength, vortiYCom, vortiY, FFTW_BACKWARD, FFTW_ESTIMATE);

    

    
    fftw_execute(planVor); //Execution of FFT

    while(iter<maxIter){
        #pragma omp parallel for default(none) shared(xLength,yLength,stream,KX,KY,vorticityCom,UvelCom,VvelCom,vortiXCom,vortiYCom) private(i,j) collapse(2)
        for(i=0;i<xLength;i++){
            for(j=0;j<yLength;j++){
                stream[i*yLength+j] = vorticityCom[i*yLength+j]/(KX[i]*KX[i]+KY[j]*KY[j]);
                UvelCom[i*yLength+j] = stream[i*yLength+j]*I*KY[j];
                VvelCom[i*yLength+j] = -1*stream[i*yLength+j]*I*KX[i];
                vortiXCom[i*yLength+j] = I*KX[i]* vorticityCom[i*yLength+j];
                vortiYCom[i*yLength+j] = I*KY[j]* vorticityCom[i*yLength+j];
            }
        }
        stream[0] = vorticityCom[0];
        UvelCom[0] = stream[0]*I*KY[0];
        VvelCom[0] = -1*stream[0]*I*KX[0];
    

        fftw_execute(planInvU); //Execution of IFFT
        fftw_execute(planInvV); //Execution of IFFT
        fftw_execute(planInvVortiX);
        fftw_execute(planInvVortiY);

        #pragma omp parallel for default(none) shared(xLength,yLength,Nflux,vortiX,Uvel,fftDefactor,vortiY,Vvel)  private(i,j) collapse(2)
        for(i=0;i<xLength;i++){
            for(j=0;j<yLength;j++){
                Nflux[i*yLength+j] = vortiX[i*yLength+j]*(Uvel[i*yLength+j]*fftDefactor+1) + vortiY[i*yLength+j]*(Vvel[i*yLength+j]*fftDefactor+1) ;
                Nflux[i*yLength+j] = Nflux[i*yLength+j]*fftDefactor;
            }
        }

        planN = fftw_plan_dft_2d(xLength, yLength, Nflux, NfluxCom, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(planN);
        fftw_destroy_plan(planN); 
        if(iter==0){
            vorticityComNew = CNfun(vorticityCom,NfluxCom,miu,KX,KY,xLength,yLength,dt);
        }else{
            vorticityComNew = CNfun_AB(vorticityCom,NfluxCom,NfluxComOld,miu,KX,KY,xLength,yLength,dt);
        }

        NfluxComOld = NfluxCom;
        NfluxCom = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);

        vorticityCom = vorticityComNew;
        iter++;
    }
    
    planInvVor = fftw_plan_dft_2d(xLength, yLength, vorticityComNew, vorticity, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(planInvVor);

    
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputTestOMP.txt");
    fd = fopen(filename,"w+");
    fprintf(fd,"timestep: %f\n",dt);
    for(i = 0; i < xLength*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%d %11.7f %11.7f\n", i, creal(vorticityComNew[i]), cimag(vorticityComNew[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
    for(i = 0; i < xLength*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%d %11.7f %11.7f\n", i, fftDefactor*creal(vorticity[i]), fftDefactor*cimag(vorticity[i]));
        //fprintf(fd, "%d %f \n", i, Nflux[i]);
    }
        
    fclose(fd);
    
    fftw_free(vorticity); fftw_free(Uvel); fftw_free(Vvel); fftw_free(vortiX); fftw_free(vortiY); fftw_free(Nflux); 
    fftw_free(xMatrix); fftw_free(yMatrix); fftw_free(KX); fftw_free(KY); fftw_free(stream);

	fftw_destroy_plan(planVor); fftw_destroy_plan(planInvVor); fftw_destroy_plan(planU); fftw_destroy_plan(planV);
    fftw_destroy_plan(planInvVortiX); fftw_destroy_plan(planInvVortiY); 

	fftw_free(vorticityCom); fftw_free(UvelCom); fftw_free(VvelCom); fftw_free(vortiXCom); fftw_free(vortiYCom);
    fftw_free(NfluxCom);      	

    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    printf("Time elapsed is %f seconds.\n", elapsed);

	return 0;
}