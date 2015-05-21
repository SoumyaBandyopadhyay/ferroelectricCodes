#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include"gasdev.c"
#include"ran1.c" 
#define nx 64
#define ny 64
#define nz 64
#define dx 1.0
#define dy 1.0
#define dz 1.0
#define dt 0.25
#define epsilon0 8.854e-12
#define PI acos(-1.0)
#define Re 1
#define Im 0
#define L 1.0
#define num_steps 1000
#define Tolerance 1.0e-08
#define initcount 0 
/****************************************************************************************************/
int main(void){
    fftw_complex *etax, *etay, *etaz, *dfbulkdetax, *dfgraddetax, *dfdipoledetax, *dfappldetax, *newetax; 
    fftw_plan p_up, p_dn;
    double *random_num1, *random_num2, *random_num3;
    double one_by_nxnynz, noise_level = 0.01;
    int nx_half, ny_half, nz_half, count;
    long  SEED = -1000;       
    FILE *fp ;  
     etax          =  (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny *nz);
     etay          =  (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny *nz);
     etaz          =  (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny *nz);
     dfbulkdetax   =  (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny * nz);
     dfgraddetax   =  (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny * nz);
     dfdipoledetax = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx * ny * nz);
     dfappldetax   = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx *ny * nz);
     newetax       = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * nx *ny * nz);
 /***********************************************************************************/ 
     nx_half = nx/2;
     ny_half = ny/2;
     nz_half = nz/2;
   
     one_by_nxnynz = 1.0/(double)(nx * ny * nz);
     
     p_up = fftw_plan_dft_3d(nx, ny, nz, etax, etax, FFTW_FORWARD, FFTW_ESTIMATE);
     p_dn = fftw_plan_dft_3d(nx, ny, nz, etax, etax, FFTW_BACKWARD, FFTW_ESTIMATE);    
  
   /*******************************************************************************/
     //initialization of the non-conserved order parameter and generating gaussian random noise//
        
     float gasdev(long *idum);
     float ran1(long *idum);

     random_num1 =  (double*) malloc(sizeof(double) * nx * ny *nz);

     random_num2 =  (double*) malloc(sizeof(double) * nx * ny *nz);

     random_num3 =  (double*) malloc(sizeof(double) * nx * ny *nz);
   

     for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
           for(int k = 0; k < nz; k++){
              etax[k + nz *(j + i * ny)][Re] = 0.0;
              etay[k + nz *(j + i * ny)][Re] = 0.0;
              etaz[k + nz *(j + i * ny)][Re] = 0.0;
              etax[k + nz *(j + i * ny)][Im] = 0.0;
              etay[k + nz *(j + i * ny)][Im] = 0.0;
              etaz[k + nz *(j + i * ny)][Im] = 0.0;
           }
        }
     }  

     for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
           for(int k = 0; k < nz; k++){
             random_num1[k + nz *(j + i * ny)] = gasdev(&SEED);
             random_num2[k + nz *(j + i * ny)] = gasdev(&SEED);
             random_num3[k + nz *(j + i * ny)] = gasdev(&SEED);
           }
         }    
      } 
 
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
           for(int k = 0; k < nz; k++){
               etax[k + nz *(j + i * ny)][Re] += random_num1[k + nz *(j + i * ny)]* noise_level;
	       etay[k + nz *(j + i * ny)][Re] += random_num2[k + nz *(j + i * ny)]* noise_level;
	       etaz[k + nz *(j + i * ny)][Re] += random_num3[k + nz *(j + i * ny)]* noise_level; 
           }
         }
       }

                fp = fopen("eta.txt", "w");
       for(int i = 0; i < nx; i++){
          for(int j = 0 ; j < ny; j++){
             for(int k = 0; k < nz; k++){
                fprintf(fp,"%d\t%d\t%d\t%le\t%le\t%le\n",i, j, k, etax[k + nz *(j+i*ny)][Re], 
                       etay[k + nz * (j + i * ny)][Re], etaz[k + nz * (j + i * ny)][Re]);
              }
                fprintf(fp,"\n");
           }
       }
                fclose(fp);       

      
   /*************************************************************************************/
    // defining material constants for PbTiO3

        double a1     = 1.339;
        double a11    =-0.3216;
        double a12    = 3.325;
        double a111   = 0.6608;
        double a112   = 1.547;
        double a123   =-9.2804;
        double lambda = 60;
        double zeta   = 2.2;
        double beta   = 2.0;
        double p0     = 0.75;
        double sigma  = 100.0;
        double df     = 7.4e07;
   /*************************************************************************************/
     
      for(int i = 0; i < nx; i++){
         for(int j = 0; j < ny; j++){
            for(int k = 0; k < nz; k++){ 
               dfbulkdetax[k + nz * (j + i * ny)][Re] = 2.0*a1*etax[k + nz *(j + i * ny)][Re] + 4.0*a11*(pow(etax[k + nz *(j + i * ny)][Re], 3.0))
                                                      + a12*(2.0*(etax[k + nz * (j + i * ny)][Re])*(pow(etay[k + nz*(j + i * ny)][Re],2.0) 
                                                      + pow(etaz[k + nz*(j + i * ny)][Re],2.0)));
                                                      + 6.0*a111*pow(etax[k + nz * (j + i * ny)][Re], 5.0) 
                                                      + a112*(2.0 * etax[k + nz*(j + i * ny)][Re] * (pow(etay[k + nz * (j + i * ny)][Re], 4.0) 
                                                      + pow(etaz[k + nz *(j + i * ny)][Re], 4.0))  
                                                      + 4.0 * (pow(etay[k + nz*(j + i * ny)][Re],2.0)*pow(etax[k + nz *(j + i * ny)][Re],3.0))
                                                      + 4.0 * (pow(etaz[k + nz * (j + i * ny)][Re], 2.0)*pow(etax[k + nz * (j + i * ny)][Re], 3.0)))
                                                      + 2.0 * a123*(etax[k + nz * (j + i * ny)][Re])*pow(etay[k + nz *(j + i * ny)][Re],2.0)
                                                       *pow(etaz[k + nz *(j + i * ny)][Re],2.0);
              
              dfbulkdetax[k + nz * (j + i * ny)][Im] = 0.0;
            }
         }
      }    
     

    

  return 0;
 } 




