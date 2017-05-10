#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<fftw3.h>


int main(int argc, char* argv[])
{
    int i=0,j;
	int Npoints = 16;
	int n0 = sqrt(Npoints),n1 = sqrt(Npoints);
	fftw_complex  *in, *out, *LastIn;
	fftw_plan plan,planInv;


    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints); 
	LastIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);			//pay attention
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);		//pay attention
	// plan = fftw_plan_dft_1d(Npoints, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 	//Here we set which kind of transformation we want to perform

    plan = fftw_plan_dft_2d(n0,n1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    planInv = fftw_plan_dft_2d(n0,n1, in, out,  FFTW_BACKWARD, FFTW_ESTIMATE);

	for(i = 0; i < Npoints; i++)
	{
		in[i] = i;
	}

	fftw_execute(plan); //Execution of FFT
    //fftw_execute(planInv); //Execution of FFT
	
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "outputComplex.txt");
    fd = fopen(filename,"w+");
    for(i = 0; i < Npoints; i++)
    {
        //fprintf(fd, "%d\n", i);
        fprintf(fd, "%11.7f %11.7f\n", creal(out[i]), cimag(out[i]));
        //fprintf(fd, "%d %f \n", i, vorticity[i]);
    }
    fclose(fd);


	fftw_destroy_plan(plan);	 //Free memory
	free(in);free(LastIn);			 //Free memory
	fftw_free(out);			 //Free memory
	return 0;
}