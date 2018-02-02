#include<iostream>
#include<stdio.h>
#include <stdlib.h>     
#include <math.h>     
#include <stdlib.h>     
#include <omp.h>     

using namespace std;
//using namespace arma;

#define PI 3.14159265358979
#define PI2 9.86960440108936
#define PIHALF 1.57079632679490

//-----------------------------------------------------------
//-----------------------------------------------------------
double* getAlphas(){
  double f;
  double* alphas = (double*) malloc(5*10*sizeof(double));
  FILE * pFile;
  pFile = fopen ("alphas.dat","r");
  for(int i = 0; i < 50; i++){
    fscanf (pFile, "%lf", &alphas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return alphas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getBetas(){
  double f;
  double* betas = (double*) malloc(5*5*10*sizeof(double));
  FILE * pFile;
  pFile = fopen ("betas.dat","r");
  for(int i = 0; i < 5*5*10; i++){
    fscanf (pFile, "%lf", &betas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return betas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getReals(double* x, double* y, int size) {

  double* P = (double*) malloc(size*size*sizeof(double));
  double* M = (double*) malloc(size*size*sizeof(double));
  double* P2 = (double*) malloc(size*size*sizeof(double));
  double* M2 = (double*) malloc(size*size*sizeof(double));
  double* P4 = (double*) malloc(size*size*sizeof(double));
  double* M4 = (double*) malloc(size*size*sizeof(double));
  double* P6 = (double*) malloc(size*size*sizeof(double));
  double* M6 = (double*) malloc(size*size*sizeof(double));
  double* P8 = (double*) malloc(size*size*sizeof(double));
  double* M8 = (double*) malloc(size*size*sizeof(double));
  
  double* Re = (double*) malloc(9*size*size*sizeof(double));
  
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
  
      P[i*size + j] =x[i]*x[j] + y[i]*y[j]; M[i*size+j] =x[i]*y[j] - y[i]*x[j]; 
  
      P2[i*size+j] =P[i*size+j]*P[i*size+j];     M2[i*size+j] =M[i*size+j]*M[i*size+j]; 
      P4[i*size+j] =P2[i*size+j]*P2[i*size+j];   M4[i*size+j] =M2[i*size+j]*M2[i*size+j]; 
      P6[i*size+j] =P4[i*size+j]*P2[i*size+j];   M6[i*size+j] =M4[i*size+j]*M2[i*size+j]; 
      P8[i*size+j] =P4[i*size+j]*P4[i*size+j];   M8[i*size+j] =M4[i*size+j]*M4[i*size+j]; 
    }
  }
  
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
  
  Re[0*size*size + i*size + j] = P[i*size+j];
  
  Re[1*size*size + i*size + j] = P2[i*size+j]
                                 - M2[i*size+j];
  Re[2*size*size + i*size + j] = P[i*size+j]*(P2[i*size+j]
                                 - 3*M2[i*size+j]);
  Re[3*size*size + i*size + j] = P4[i*size+j] + M4[i*size+j]
                                 - 6*P2[i*size+j]*M2[i*size+j];
  Re[4*size*size + i*size + j] = P[i*size+j]*(P4[i*size+j] - 10*P2[i*size+j]*M2[i*size+j]
                                 + 5*M4[i*size+j]);
  Re[5*size*size + i*size + j] = P6[i*size+j] - M6[i*size+j]  - 15*(P4[i*size+j]*M2[i*size+j]
                                 - P2[i*size+j]*M4[i*size+j]);
  Re[6*size*size + i*size + j] = P[i*size+j]*(P6[i*size+j] - 7*M6[i*size+j] - 21*P4[i*size+j]*M2[i*size+j]
                                 + 35*P2[i*size+j]*M4[i*size+j]);
  Re[7*size*size + i*size + j] = P8[i*size+j] + M8[i*size+j] - 28*(P6[i*size+j]*M2[i*size+j] + P2[i*size+j]*M6[i*size+j])
                                 + 70*P4[i*size+j]*M4[i*size+j];
  Re[8*size*size + i*size + j] = P[i*size+j]*(P8[i*size+j] + 9*M8[i*size+j] - 36*P6[i*size+j]*M2[i*size+j]
                                 + 126*P4[i*size+j]*M4[i*size+j] - 84*P2[i*size+j]*M6[i*size+j]);
    }
  }
  free(P);  free(M); free(P2); free(M2); free(P4); free(M4); free(P6); free(M6); free(P8); free(M8);
  
  return Re;

}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getR2(double* x, double* y, double* z, int size){

  double* r2s = (double*) malloc(size*sizeof(double));

  for(int i = 0; i < size; i++){
       r2s[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i]; 
  }

 return r2s; 

}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP0(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize){

  double* P0 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;// = (double*) malloc(Nsize*sizeof(double));
  double oneO1PalphaSqrt;// = (double*) malloc(Nsize*sizeof(double));

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrtCubed = (double*) malloc(Nsize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrtCubed[n] = oneO1PalphaSqrt*oneO1Palpha;
    alphaO1Palpha[n] = alphas[n]*oneO1Palpha;
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += exp(-alphaO1Palpha[k]*r2[i] - alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[n*Nsize + k]*betas[nd*Nsize + kd]*oneO1PalphaSqrtCubed[k]*oneO1PalphaSqrtCubed[kd]*sumsInner;
        }
      }
          P0[Nsize*n + nd] = PI2*0.25*sumsOuter; 
    }
  }

 free(oneO1PalphaSqrtCubed);
 free(alphaO1Palpha);

 return P0; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP1(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P1 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt5 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt5[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[Nsize + n]*oneO1Palpha;
  }
  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = ReX[i*Asize + j] + z[i]*z[j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i] - alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[Nsize*Nsize + n*Nsize + k]*betas[Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt5[k]*oneO1PalphaSqrt5[kd]*sumsInner;
        }
      }
          P1[Nsize*n + nd] = 3.0*PI2*0.25*sumsOuter; 
    }
  }

 free(oneO1PalphaSqrt5);
 free(alphaO1Palpha);

 return P1; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP2(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P2 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt7 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zxy2 = (double*) malloc(Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[2*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt7[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[2*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    zxy2[i] = (2*z[i]*z[i] - x[i]*x[i] - y[i]*y[i]);
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 3*ReX[Asize*Asize + i*Asize + j] + 12*z[i]*z[j]*ReX[i*Asize + j] + zxy2[i]*zxy2[j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[2*Nsize*Nsize + n*Nsize + k]*betas[2*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt7[k]*oneO1PalphaSqrt7[kd]*sumsInner;
        }
      }
          P2[Nsize*n + nd] = 5.0*PI2*0.0625*sumsOuter; 
    }
  }

  free(oneO1PalphaSqrt7);
  free(alphaO1Palpha);

  return P2; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP3(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P3 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt9 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zxy4 = (double*) malloc(Asize*sizeof(double));
  double* zxy233 = (double*) malloc(Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[3*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt9[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[3*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    zxy4[i] = (4*z[i]*z[i] - x[i]*x[i] - y[i]*y[i]);
    zxy233[i] = (2*z[i]*z[i] - 3*x[i]*x[i] - 3*y[i]*y[i]);
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 2.5*ReX[2*Asize*Asize + i*Asize + j] + z[i]*z[j]*(15.0*ReX[Asize*Asize + i*Asize + j] + zxy233[i]*zxy233[j]) + 1.5*zxy4[i]*zxy4[j]*ReX[i*Asize + j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[3*Nsize*Nsize + n*Nsize + k]*betas[3*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt9[k]*oneO1PalphaSqrt9[kd]*sumsInner;
        }
      }
          P3[Nsize*n + nd] = 7.0*PI2*0.0625*sumsOuter; 
    }
  }
    for(int i = 0; i < Nsize; i++){
           for(int j = 0; j < Nsize; j++){
            cout << betas[3*Nsize*Nsize + i*Nsize + j] << " ";
           }
           cout << endl;
    }

  free(oneO1PalphaSqrt9);
  free(alphaO1Palpha);

  return P3; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP4(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P4 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt11 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr7 = (double*) malloc(Asize*sizeof(double));
  double* zr73 = (double*) malloc(Asize*sizeof(double));
  double* zr35 = (double*) malloc(Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[4*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt11[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[4*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i];
    zr7[i] = 7.0*z2 - r2[i];
    zr73[i] = 7.0*z2 - 3.0*r2[i];
    zr35[i] = 35*z2*z2 - 30*z2*r2[i] + 3*r2[i]*r2[i];
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 35.0*ReX[3*Asize*Asize + i*Asize + j] + z[i]*z[j]*(280*ReX[2*Asize*Asize + i*Asize + j] + 40*zr73[i]*zr73[j]*ReX[i*Asize + j]) + 20*zr7[i]*zr7[j]*ReX[Asize*Asize + i*Asize + j] + zr35[i]*zr35[j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[4*Nsize*Nsize + n*Nsize + k]*betas[4*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt11[k]*oneO1PalphaSqrt11[kd]*sumsInner;
        }
      }
          P4[Nsize*n + nd] = 9.0*PI2*0.00390625*sumsOuter; 
    }
  }
    for(int i = 0; i < Nsize; i++){
           for(int j = 0; j < Nsize; j++){
            cout << betas[4*Nsize*Nsize + i*Nsize + j] << " ";
           }
           cout << endl;
    }

  free(oneO1PalphaSqrt11);
  free(alphaO1Palpha);

  return P4; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP5(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P5 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt13 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr9 = (double*) malloc(Asize*sizeof(double));
  double* zr3 = (double*) malloc(Asize*sizeof(double));
  double* zr21 = (double*) malloc(Asize*sizeof(double));
  double* zr63 = (double*) malloc(Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[5*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt13[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[5*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i];
    zr9[i] = 9.0*z2 - r2[i];
    zr3[i] = 3.0*z2 - r2[i];
    zr21[i] = 21*z2*z2 - 14*z2*r2[i] + r2[i]*r2[i];
    zr63[i] = 63*z2*z2 - 70*z2*r2[i] + 15*r2[i]*r2[i];
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 31.5*ReX[4*Asize*Asize + i*Asize + j] + z[i]*z[j]*(315*ReX[3*Asize*Asize + i*Asize + j] + 410*zr3[i]*zr3[j]*ReX[Asize*Asize + i*Asize + j] + zr63[i]*zr63[j]) + 17.5*zr9[i]*zr9[j]*ReX[2*Asize*Asize + i*Asize + j] + 15*zr21[i]*zr21[j]*ReX[i*Asize + j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[5*Nsize*Nsize + n*Nsize + k]*betas[5*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt13[k]*oneO1PalphaSqrt13[kd]*sumsInner;
        }
      }
          P5[Nsize*n + nd] = 11.0*PI2*0.00390625*sumsOuter; 
    }
  }
    for(int i = 0; i < Nsize; i++){
           for(int j = 0; j < Nsize; j++){
            cout << betas[5*Nsize*Nsize + i*Nsize + j] << " ";
           }
           cout << endl;
    }

  free(oneO1PalphaSqrt13);
  free(alphaO1Palpha);

  return P5; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP6(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P6 = (double*) malloc(Nsize*Nsize*sizeof(double));

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt15 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr11 = (double*) malloc(Asize*sizeof(double));
  double* zr113 = (double*) malloc(Asize*sizeof(double));
  double* zr3311 = (double*) malloc(Asize*sizeof(double));
  double* zr3330 = (double*) malloc(Asize*sizeof(double));
  double* zr231 = (double*) malloc(Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[6*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt15[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[6*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i];
    zr11[i] = 11.0*z2 - r2[i];
    zr113[i] = 11.0*z2 - 3.0*r2[i];
    zr3311[i] = 33*z2*z2 - 11*z2*r2[i] + r2[i]*r2[i];
    zr3330[i] = 33*z2*z2 - 30*z2*r2[i] + 5*r2[i]*r2[i];
    zr231[i] = 231*z2*z2*z2 - 315*z2*z2*r2[i] + 105*z2*r2[i]*r2[i] - 5*r2[i]*r2[i]*r2[i];
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 115.5*ReX[5*Asize*Asize + i*Asize + j] + z[i]*z[j]*(1386*ReX[4*Asize*Asize + i*Asize + j] + 210*zr113[i]*zr113[j]*ReX[2*Asize*Asize + i*Asize + j] + 82*zr3330[i]*zr3330[j]*ReX[i*Asize + j]) + 63*zr11[i]*zr11[j]*ReX[3*Asize*Asize + i*Asize + j] + 52.5*zr3311[i]*zr3311[j]*ReX[Asize*Asize + i*Asize + j] + zr231[i]*zr231[j];
    }
  }

  for(int n = 0; n < Nsize; n++){
    for(int nd = 0; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
            }
          }
          sumsOuter +=  betas[6*Nsize*Nsize + n*Nsize + k]*betas[6*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt15[k]*oneO1PalphaSqrt15[kd]*sumsInner;
        }
      }
          P6[Nsize*n + nd] = 13.0*PI2*9.765625e-4*sumsOuter; 
    }
  }
    for(int i = 0; i < Nsize; i++){
           for(int j = 0; j < Nsize; j++){
            cout << betas[6*Nsize*Nsize + i*Nsize + j] << " ";
           }
           cout << endl;
    }

  free(oneO1PalphaSqrt15);
  free(alphaO1Palpha);

  return P6; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
int main(int argc, char* argv[]) {
int Asize = 500;
int Nsize = 5;
double* x = (double*) malloc(Asize*sizeof(double));
double* y = (double*) malloc(Asize*sizeof(double));
double* z = (double*) malloc(Asize*sizeof(double));
double* alphas = getAlphas();
double* betas = getBetas();

x[0] = 0.1;
x[1] = -1.0;
y[0] = 1.0;
y[1] = -1.0;
z[0] = 1.0;
z[1] = -1.0;

double* r2 = getR2(x, y, z, Asize);
double* ReX = getReals(x, y,Asize); 
double* P0;
double* P1;
double* P2;
double* P3;
double* P4;
double* P5;
double* P6;

#pragma omp parallel sections
{
    #pragma omp section
    { 
        P0 = getP0(x,y,z,r2,alphas, betas ,Asize, Nsize);
    }
    #pragma omp section
    { 
        P1 = getP1(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }
    #pragma omp section
    { 
        P2 = getP2(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }
    #pragma omp section
    { 
        P3 = getP3(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }
    #pragma omp section
    { 
        P4 = getP4(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }
    #pragma omp section
    { 
        P5 = getP5(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }
    #pragma omp section
    { 
        P6 = getP6(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
    }

}



  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P2_" << i + 1 << j + 1 <<" : " << P2[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;
  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P3_" << i + 1 << j + 1 <<" : " << P3[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;
  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P4_" << i + 1 << j + 1 <<" : " << P4[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;
  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P5_" << i + 1 << j + 1 <<" : " << P5[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;
  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P6_" << i + 1 << j + 1 <<" : " << P6[i*Nsize + j];
    }
    cout << endl;
  }

//free(x);
//free(y);
//free(f);

}
