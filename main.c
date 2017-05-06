#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<fftw3.h>
#define pi 3.1415926


fftw_complex* exactVorticity(double L, double miu, fftw_complex* xMatrix, fftw_complex* yMatrix, double v0, int xLength, int yLength, double t){
    fftw_complex* ans = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);
    int i,j;
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
    double timeInv = 1/dt;
    for(i=0;i<xLength;i++){
        for(j=0;j<yLength;j++){
            double temp = (KX[i]*KX[i]+KY[j]*KY[j])*miu*0.5;
            ans[i*yLength+j] = ((timeInv-temp)*vortiCom[i*yLength+j]+NfluxCom[i*yLength+j])/(timeInv+temp);
        }
    }
    return ans;
}


int main(int argc, char* argv[])
{
    int i,j;
	int Npoints = 16;
	int xLength = Npoints, yLength = Npoints;
	double t = 0, miu = 0.05, v0 = 1, L = 1, dx = L/xLength, dy = L/yLength, dt = 1/64;
	fftw_complex *vorticity, *Uvel, *Vvel, *vortiX, *vortiY, *Nflux; 
    fftw_complex *xMatrix, *yMatrix, *KX, *KY;
    fftw_complex *vorticityCom, *UvelCom, *VvelCom, *vortiXCom, *vortiYCom, *NfluxCom, *vorticityComNew;
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

    vorticity = exactVorticity(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
	Uvel = exactU(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
    Vvel = exactV(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
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
    planInvVor = fftw_plan_dft_2d(xLength, yLength, vorticity, vorticityCom, FFTW_BACKWARD, FFTW_ESTIMATE);
    // planU = fftw_plan_dft_2d(xLength, yLength, Uvel, UvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    // planInvU = fftw_plan_dft_2d(xLength, yLength, Uvel, UvelCom, FFTW_BACKWARD, FFTW_ESTIMATE);    
    // planV = fftw_plan_dft_2d(xLength, yLength, Vvel, VvelCom, FFTW_FORWARD, FFTW_ESTIMATE);
    // planInvV = fftw_plan_dft_2d(xLength, yLength, Vvel, VvelCom, FFTW_BACKWARD, FFTW_ESTIMATE);

    planInvVortiX = fftw_plan_dft_2d(xLength, yLength, vortiX, vortiXCom, FFTW_BACKWARD, FFTW_ESTIMATE);
    planInvVortiY = fftw_plan_dft_2d(xLength, yLength, vortiY, vortiYCom, FFTW_BACKWARD, FFTW_ESTIMATE);

    planN = fftw_plan_dft_2d(xLength, yLength, Nflux, NfluxCom, FFTW_FORWARD, FFTW_ESTIMATE);


    fftw_execute(planVor); //Execution of FFT
    //fftw_execute(planInvVor);
    //fftw_execute(planInvVor); 
    // for(i=0;i<xLength;i++){
    //     for(j=0;j<yLength;j++){
    //         vortiXCom[i*yLength+j] = I*KX[i]* vorticityCom[i*yLength+j];
    //         vortiYCom[i*yLength+j] = I*KY[j]* vorticityCom[i*yLength+j];
    //     }
    // }

    // fftw_execute(planInvVortiX); //Execution of IFFT
    // fftw_execute(planInvVortiY); //Execution of IFFT
    
    // for(i=0;i<xLength;i++){
    //     for(j=0;j<yLength;j++){
    //         Nflux[i*yLength+j] = vortiX[i*yLength+j]*Uvel[i*yLength+j] + vortiY[i*yLength+j]*Vvel[i*yLength+j] ;
    //         //printf("%f \n",Nflux[i*yLength+j]);
    //     }
    // }
   
    // fftw_execute(planN); 

    // vorticityComNew = CNfun(vorticityCom,NfluxCom,miu,KX,KY,xLength,yLength,dt);
   

    
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputTest.txt");
    fd = fopen(filename,"w+");
    for(i = 0; i < xLength*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%d %11.7f %11.7f\n", i, creal(vorticityCom[i]), cimag(vorticityCom[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
    for(i = 0; i < xLength*yLength; i++)
    {
        //fprintf(fd, "%d\n", i);
        //fprintf(fd, "%d %11.7f %11.7f\n", i, creal(vorticityCom[i]), cimag(vorticityCom[i]));
        fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
        
    fclose(fd);
    // for(i = 0; i < xLength; i++)
    //   fprintf(fd, "  %f\n", initialU[i]);
      
    // for(i = 0; i < xLength; i++)
    //   // fprintf(fd, "%d %11.7f %11.7f\n", i, creal(out[i]), cimag(out[i]));
    //   fprintf(fd, "  %f\n", initialV[i*yLength]);
    

    free(vorticity); free(Uvel); free(Vvel); free(vortiX); free(vortiY); free(Nflux); 
    free(xMatrix); free(yMatrix); free(KX); free(KY);

	fftw_destroy_plan(planVor); fftw_destroy_plan(planInvVor); 
    fftw_destroy_plan(planInvVortiX); fftw_destroy_plan(planInvVortiY); fftw_destroy_plan(planN);

	fftw_free(vorticityCom); fftw_free(UvelCom); fftw_free(VvelCom); fftw_free(vortiXCom); fftw_free(vortiYCom);
    fftw_free(NfluxCom);      	
	return 0;
}