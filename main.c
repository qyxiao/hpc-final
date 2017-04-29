#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<fftw3.h>
#define pi 3.1415926


double* exactVorticity(double L, double miu, double* xMatrix, double* yMatrix, double v0, int xLength, int yLength, double t){
    double* ans = (double *) malloc(sizeof(double)*xLength*yLength);
    int i,j;
    for(i=0;i<xLength;i++){
    	for(j=0;j<yLength;j++){
    		ans[i*yLength+j] = 8*pi/L*exp(-8*pi*pi*miu*t/(L*L))*cos((xMatrix[i]-v0*t)*2*pi/L)*cos((yMatrix[j]-v0*t)*2*pi/L);
    	}
    }
    return ans;
}


int main(int argc, char* argv[])
{
    int i,j;
	int Npoints = 16;
	int xLength = Npoints, yLength = Npoints;
	double t = 0, miu = 0.05, v0 = 1, L = 1, dx = L/xLength, dy = L/yLength;
	fftw_complex *in, *out;
	fftw_plan plan;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);			//pay attention
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*xLength*yLength);		//pay attention
    plan = fftw_plan_dft_2d(xLength, yLength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
	// for(i = 0; i < xLength*yLength; i++)
	// {
	// 	in[i] = i;
	// }
	// fftw_execute(plan); //Execution of FFT
    double* xMatrix = (double *) malloc(sizeof(double)*xLength);
    double* yMatrix = (double *) malloc(sizeof(double)*yLength);
    for(i=0; i<xLength; i++){
    	xMatrix[i]=i*dx;
    }
    for(i=0; i<yLength; i++){
    	yMatrix[i]=i*dy;
    } 
   

    double* vorticity = exactVorticity(L, miu, xMatrix, yMatrix, v0, xLength, yLength, t);
	
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputTest.txt");
    fd = fopen(filename,"w+");
    for(i = 0; i < xLength*yLength; i++)
      // fprintf(fd, "%d %11.7f %11.7f\n", i, creal(out[i]), cimag(out[i]));
      fprintf(fd, "  %f\n", vorticity[i]);
    fclose(fd);


	fftw_destroy_plan(plan);	 //Free memory
	fftw_free(in);			 //Free memory
	fftw_free(out);			 //Free memory
	return 0;
}