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
float* getAlphas(){
  float f;
  float* alphas = (float*) malloc(5*10*sizeof(float));
  FILE * pFile;
  pFile = fopen ("alphas.dat","r");
  for(int i = 0; i < 50; i++){
    fscanf (pFile, "%f", &alphas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return alphas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getBetas(){
  float f;
  float* betas = (float*) malloc(5*5*10*sizeof(float));
  FILE * pFile;
  pFile = fopen ("betas.dat","r");
  for(int i = 0; i < 5*5*10; i++){
    fscanf (pFile, "%f", &betas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return betas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getReals(float* x, float* y, int size) {

  float* P = (float*) malloc(size*size*sizeof(float));
  float* M = (float*) malloc(size*size*sizeof(float));
  float* P2 = (float*) malloc(size*size*sizeof(float));
  float* M2 = (float*) malloc(size*size*sizeof(float));
  float* P4 = (float*) malloc(size*size*sizeof(float));
  float* M4 = (float*) malloc(size*size*sizeof(float));
  float* P6 = (float*) malloc(size*size*sizeof(float));
  float* M6 = (float*) malloc(size*size*sizeof(float));
  float* P8 = (float*) malloc(size*size*sizeof(float));
  float* M8 = (float*) malloc(size*size*sizeof(float));
  
  float* Re = (float*) malloc(9*size*size*sizeof(float));
  
//#pragma omp parallel for schedule(static,chunk) 
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
  
      P[i*size + j] =x[i]*x[j] + y[i]*y[j]; M[i*size+j] =x[i]*y[j] - y[i]*x[j]; 
  
      P2[i*size+j] =P[i*size+j]*P[i*size+j];     M2[i*size+j] =M[i*size+j]*M[i*size+j]; 
      P4[i*size+j] =P2[i*size+j]*P2[i*size+j];   M4[i*size+j] =M2[i*size+j]*M2[i*size+j]; 
      P6[i*size+j] =P4[i*size+j]*P2[i*size+j];   M6[i*size+j] =M4[i*size+j]*M2[i*size+j]; 
      P8[i*size+j] =P4[i*size+j]*P4[i*size+j];   M8[i*size+j] =M4[i*size+j]*M4[i*size+j]; 
    }
  }
  
//#pragma omp parallel for schedule(static,chunk) 
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
float* getR2(float* x, float* y, float* z, int size){

  float* r2s = (float*) malloc(size*sizeof(float));

  for(int i = 0; i < size; i++){
       r2s[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i]; 
  }

 return r2s; 

}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP0(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize){

  float* P0 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;// = (float*) malloc(Nsize*sizeof(float));
  float oneO1PalphaSqrt;// = (float*) malloc(Nsize*sizeof(float));

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrtCubed = (float*) malloc(Nsize*sizeof(float));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrtCubed[n] = oneO1PalphaSqrt*oneO1Palpha;
    alphaO1Palpha[n] = alphas[n]*oneO1Palpha;
  }

//#pragma omp parallel for schedule(static,chunk)
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
float* getP1(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P1 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt5 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
 free(ReXStripe);

 return P1; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP2(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P2 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt7 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));
  float* zxy2 = (float*) malloc(Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
  free(ReXStripe);
  free(zxy2);

  return P2; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP3(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P3 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt9 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));
  float* zxy4 = (float*) malloc(Asize*sizeof(float));
  float* zxy233 = (float*) malloc(Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
//    for(int i = 0; i < Nsize; i++){
//           for(int j = 0; j < Nsize; j++){
//            cout << betas[3*Nsize*Nsize + i*Nsize + j] << " ";
//           }
//           cout << endl;
//    }

  free(oneO1PalphaSqrt9);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zxy4);
  free(zxy233);

  return P3; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP4(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P4 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;
  float z2=0;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt11 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));
  float* zr7 = (float*) malloc(Asize*sizeof(float));
  float* zr73 = (float*) malloc(Asize*sizeof(float));
  float* zr35 = (float*) malloc(Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
//    for(int i = 0; i < Nsize; i++){
//           for(int j = 0; j < Nsize; j++){
//            cout << betas[4*Nsize*Nsize + i*Nsize + j] << " ";
//           }
//           cout << endl;
//    }

  free(oneO1PalphaSqrt11);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr7);
  free(zr73);
  free(zr35);

  return P4; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP5(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P5 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;
  float z2=0;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt13 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));
  float* zr9 = (float*) malloc(Asize*sizeof(float));
  float* zr3 = (float*) malloc(Asize*sizeof(float));
  float* zr21 = (float*) malloc(Asize*sizeof(float));
  float* zr63 = (float*) malloc(Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
//    for(int i = 0; i < Nsize; i++){
//           for(int j = 0; j < Nsize; j++){
//            cout << betas[5*Nsize*Nsize + i*Nsize + j] << " ";
//           }
//           cout << endl;
//    }

  free(oneO1PalphaSqrt13);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr9); 
  free(zr3); 
  free(zr21);
  free(zr63);

  return P5; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
float* getP6(float* x, float* y, float* z,float* r2, float* alphas, float* betas , int Asize, int Nsize, float* ReX){

  float* P6 = (float*) malloc(Nsize*Nsize*sizeof(float));

  float sumsInner = 0;
  float sumsOuter = 0;
  float oneO1Palpha;
  float oneO1PalphaSqrt;
  float z2=0;

  float* alphaO1Palpha = (float*) malloc(Nsize*sizeof(float));
  float* oneO1PalphaSqrt15 = (float*) malloc(Nsize*sizeof(float));
  float* ReXStripe = (float*) malloc(Asize*Asize*sizeof(float));
  float* zr11 = (float*) malloc(Asize*sizeof(float));
  float* zr113 = (float*) malloc(Asize*sizeof(float));
  float* zr3311 = (float*) malloc(Asize*sizeof(float));
  float* zr3330 = (float*) malloc(Asize*sizeof(float));
  float* zr231 = (float*) malloc(Asize*sizeof(float));

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

//#pragma omp parallel for schedule(static,chunk)
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
//    for(int i = 0; i < Nsize; i++){
//           for(int j = 0; j < Nsize; j++){
//            cout << betas[6*Nsize*Nsize + i*Nsize + j] << " ";
//           }
//           cout << endl;
//    }

  free(oneO1PalphaSqrt15);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr11);
  free(zr113);
  free(zr3311);
  free(zr3330);
  free(zr231);

  return P6; 
}
//-----------------------------------------------------------
//-----------------------------------------------------------
int main(int argc, char* argv[]) {
int Asize = 10000;
int Nsize = 5;
float* x = (float*) malloc(Asize*sizeof(float));
float* y = (float*) malloc(Asize*sizeof(float));
float* z = (float*) malloc(Asize*sizeof(float));
float* alphas = getAlphas();
float* betas = getBetas();

x[0] = 0.1;
x[1] = -1.0;
y[0] = 1.0;
y[1] = -1.0;
z[0] = 1.0;
z[1] = -1.0;

float* r2 = getR2(x, y, z, Asize);
float* ReX = getReals(x, y,Asize); 
float* P0;
float* P1;
float* P2;
float* P3;
float* P4;
float* P5;
float* P6;

#pragma omp parallel sections
{
    #pragma omp section
    { 
        P0 = getP0(x,y,z,r2,alphas, betas ,Asize, Nsize);
//        cout << "Done P0" << endl;
    }
    #pragma omp section
    { 
        P1 = getP1(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P1" << endl;
    }
    #pragma omp section
    { 
        P2 = getP2(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P2" << endl;
    }
    #pragma omp section
    { 
        P3 = getP3(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P3" << endl;
    }
    #pragma omp section
    { 
        P4 = getP4(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P4" << endl;
    }
    #pragma omp section
    { 
        P5 = getP5(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P5" << endl;
    }
    #pragma omp section
    { 
        P6 = getP6(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
//        cout << "Done P6" << endl;
    }

}

  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P0_" << i + 1 << j + 1 <<" : " << P0[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;

  for(int i = 0; i <Nsize ; i++){
    for(int j = 0; j <Nsize ; j++){
      cout << " P1_" << i + 1 << j + 1 <<" : " << P1[i*Nsize + j];
    }
    cout << endl;
  }
    cout << endl;
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

free(x);
free(y);
free(z);
free(P0);
free(P1);
free(P2);
free(P3);
free(P4);
free(P5);
free(P6);

}
