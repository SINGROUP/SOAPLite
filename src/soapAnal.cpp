#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

using namespace std;
//using namespace arma;

#define PI     3.14159265358979
#define PI2    9.86960440108936
#define PIHALF 1.57079632679490

//-----------------------------------------------------------
//-----------------------------------------------------------
double* getAlphas(int alphaSize){
  double f;
  double* alphas = (double*) malloc(alphaSize*10*sizeof(double));
  FILE * pFile;
  pFile = fopen ("alphas5.dat","r");
  for(int i = 0; i < alphaSize*10; i++){
    fscanf (pFile, "%lf", &alphas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return alphas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getBetas(int alphaSize){
  double f;
  double* betas = (double*) malloc(alphaSize*alphaSize*10*sizeof(double));
  FILE * pFile;
//  pFile = fopen ("betasCorr.dat","r");
  pFile = fopen ("betas5.dat","r");
  for(int i = 0; i < alphaSize*alphaSize*10; i++){
    fscanf (pFile, "%lf", &betas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return betas;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
int* getInfo(int* totalAN, int* Ntypes){
  FILE* pFile;
  //Getting meta data
  pFile = fopen("metadata.dat","r");
    fscanf(pFile, "%i", totalAN);
    fscanf(pFile, "%i", Ntypes);
  fclose (pFile);
  //Getting atom type counts
  int* typeNs = (int*) malloc(sizeof(int)*Ntypes[0]);
  pFile = fopen("atomtypecount.dat","r");
  for(int i=0; i < Ntypes[0]; i++){
    fscanf(pFile, "%d", &typeNs[i]);
  }
  fclose (pFile);
  return typeNs;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getApos(int* totalAN, int* Ntypes, int* typeNs, int*types){
  FILE* pFile;
  //Getting atom types and positions
  pFile = fopen("type_pos.dat","r");
  int marcher=0;
  double* Apos = (double*) malloc(3*sizeof(double)*totalAN[0]);
  for(int i=0; i < Ntypes[0]; i++){
    fscanf(pFile, "%d", &types[i]);
    for(int j=0; j < typeNs[i] ; j++){
      fscanf(pFile, "%lf", &Apos[3*marcher    ]);
      fscanf(pFile, "%lf", &Apos[3*marcher + 1]);
      fscanf(pFile, "%lf", &Apos[3*marcher + 2]);
      marcher++;
    }
  }
  fclose (pFile);
  return Apos;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getHpos(int Hsize){

  FILE* pFile;
//  pFile = fopen("Hsym.dat","r");
  pFile = fopen("H.dat","r");
//  pFile = fopen("HRot.dat","r");
//  pFile = fopen("PM.dat","r");
  double* Pos = (double*) malloc(3*sizeof(double)*Hsize);
  for(int i=0; i < Hsize; i++){
      fscanf(pFile, "%lf", &Pos[3*i    ]);
      fscanf(pFile, "%lf", &Pos[3*i + 1]);
      fscanf(pFile, "%lf", &Pos[3*i + 2]);
  }
  fclose(pFile);
  return Pos;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
void getPos(double* x, double* y, double* z, double* Apos, double* Hpos, int* typeNs, int Ihpos, int Itype){

  int shiftType = 0;
    for(int i = 0; i < Itype ; i++){
      shiftType += typeNs[i];
    }

    for(int i = 0; i < typeNs[Itype]; i++){
      x[i] =  Apos[3*shiftType + 3*i    ] - Hpos[3*Ihpos    ];
      y[i] =  Apos[3*shiftType + 3*i + 1] - Hpos[3*Ihpos + 1];
      z[i] =  Apos[3*shiftType + 3*i + 2] - Hpos[3*Ihpos + 2];
    }

}
//-----------------------------------------------------------
//-----------------------------------------------------------
int getAllPos(double* x, double* y, double* z, double* Apos, double* Hpos, int* typeNs, int Ihpos,int sizeAll){

  int count = 0;
  int X, Y, Z;
  for(int i = 0; i < sizeAll; i++){
      X =  Apos[3*i    ] - Hpos[3*Ihpos    ];
      Y =  Apos[3*i + 1] - Hpos[3*Ihpos + 1];
      Z =  Apos[3*i + 2] - Hpos[3*Ihpos + 2];
      if( X*X + Y*Y + Z*Z < 50.0 ){
        x[count] = X;
        y[count] = Y;
        z[count] = Z;
        count++;
      }
  }

  return count;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
int getFilteredPos(double* x, double* y, double* z, double* Apos, double* Hpos,
   int* typeNs, int Ihpos, int Itype){

  int shiftType = 0;
  int count = 0;
  double X = 0; double Y = 0; double Z = 0;

    for(int i = 0; i < Itype ; i++){
      shiftType += typeNs[i];
    }

    for(int i = 0; i < typeNs[Itype]; i++){
      X = Apos[3*shiftType + 3*i    ] - Hpos[3*Ihpos    ];
      Y = Apos[3*shiftType + 3*i + 1] - Hpos[3*Ihpos + 1];
      Z = Apos[3*shiftType + 3*i + 2] - Hpos[3*Ihpos + 2];
      if( X*X + Y*Y + Z*Z < 50.0 ){
        x[count] = X;
        y[count] = Y;
        z[count] = Z;
        count++;
      }
    }
  return count;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getReals(double* x, double* y, int size) {

  double* p  = (double*) malloc(size*size*sizeof(double));
  double* m  = (double*) malloc(size*size*sizeof(double));
  double* p2 = (double*) malloc(size*size*sizeof(double));
  double* m2 = (double*) malloc(size*size*sizeof(double));
  double* p4 = (double*) malloc(size*size*sizeof(double));
  double* m4 = (double*) malloc(size*size*sizeof(double));
  double* p6 = (double*) malloc(size*size*sizeof(double));
  double* m6 = (double*) malloc(size*size*sizeof(double));
  double* p8 = (double*) malloc(size*size*sizeof(double));
  double* m8 = (double*) malloc(size*size*sizeof(double));

  double* Re = (double*) malloc(9*size*size*sizeof(double));

//#pragma omp parallel for schedule(static,chunk)
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){

      p[i*size + j] =x[i]*x[j] + y[i]*y[j]; m[i*size+j] =x[i]*y[j] - y[i]*x[j];
      p2[i*size+j] =p[i*size+j]*p[i*size+j];     m2[i*size+j] =m[i*size+j]*m[i*size+j];
      p4[i*size+j] =p2[i*size+j]*p2[i*size+j];   m4[i*size+j] =m2[i*size+j]*m2[i*size+j];
      p6[i*size+j] =p4[i*size+j]*p2[i*size+j];   m6[i*size+j] =m4[i*size+j]*m2[i*size+j];
      p8[i*size+j] =p4[i*size+j]*p4[i*size+j];   m8[i*size+j] =m4[i*size+j]*m4[i*size+j];

    }
  }

//#pragma omp parallel for schedule(static,chunk)
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){

  Re[0*size*size + i*size + j] = p[i*size+j];

  Re[1*size*size + i*size + j] = p2[i*size+j]
                                 - m2[i*size+j];
  Re[2*size*size + i*size + j] = p[i*size+j]*(p2[i*size+j]
                                 - 3*m2[i*size+j]);
  Re[3*size*size + i*size + j] = p4[i*size+j] + m4[i*size+j]
                                 - 6*p2[i*size+j]*m2[i*size+j];
  Re[4*size*size + i*size + j] = p[i*size+j]*(p4[i*size+j] - 10*p2[i*size+j]*m2[i*size+j]
                                 + 5*m4[i*size+j]);
  Re[5*size*size + i*size + j] = p6[i*size+j] - m6[i*size+j]
    - 15*(p4[i*size+j]*m2[i*size+j] - p2[i*size+j]*m4[i*size+j]);
  Re[6*size*size + i*size + j] = p[i*size+j]*(p6[i*size+j] - 7*m6[i*size+j] - 21*p4[i*size+j]*m2[i*size+j]
                                 + 35*p2[i*size+j]*m4[i*size+j]);
  Re[7*size*size + i*size + j] = p8[i*size+j] + m8[i*size+j] - 28*(p6[i*size+j]*m2[i*size+j] + p2[i*size+j]*m6[i*size+j])
                                 + 70*p4[i*size+j]*m4[i*size+j];
  Re[8*size*size + i*size + j] = p[i*size+j]*(p8[i*size+j] + 9*m8[i*size+j] - 36*p6[i*size+j]*m2[i*size+j]
                                 + 126*p4[i*size+j]*m4[i*size+j] - 84*p2[i*size+j]*m6[i*size+j]);
    }
  }

  free(p);  free(m); free(p2); free(m2); free(p4); free(m4); free(p6); free(m6); free(p8); free(m8);
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
double* getP0(double* x, double* y, double* z,double* r2, double* alphas,
   double* betas, int Asize, int Nsize){

  double* P0 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P0[n] = 0.0;
  }
  return P0;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;// = (double*) malloc(Nsize*sizeof(double));
  double oneO1PalphaSqrt;// = (double*) malloc(Nsize*sizeof(double));

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrtCubed = (double*) malloc(Nsize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrtCubed[n] = oneO1PalphaSqrt*oneO1Palpha;
    alphaO1Palpha[n] = alphas[n]*oneO1Palpha;
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
                      = exp(-alphaO1Palpha[k]*r2[i] - alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }

//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[n*Nsize + k]*betas[nd*Nsize + kd]*oneO1PalphaSqrtCubed[k]*oneO1PalphaSqrtCubed[kd]*sumsInner;
        }
      }
          P0[Nsize*n + nd] = PI2*0.25*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P0[Nsize*nd + n] = P0[Nsize*n + nd];
    }
  }

 free(oneO1PalphaSqrtCubed);
 free(alphaO1Palpha);
 free(kkdij);

 return P0;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP1(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P1 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P1[n] = 0.0;
  }
  return P1;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt5 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

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
  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i] - alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }

//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[Nsize*Nsize + n*Nsize + k]*betas[Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt5[k]*oneO1PalphaSqrt5[kd]*sumsInner;
        }
      }
          P1[Nsize*n + nd] = 3.0*PI2*0.25*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P1[Nsize*nd + n] = P1[Nsize*n + nd];
    }
  }

 free(oneO1PalphaSqrt5);
 free(alphaO1Palpha);
 free(ReXStripe);
 free(kkdij);

 return P1;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP2(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P2 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P2[n] = 0.0;
  }
  return P2;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt7 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zxy2 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

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
      ReXStripe[i*Asize + j] = (3*ReX[Asize*Asize + i*Asize + j] + 12*z[i]*z[j]*ReX[i*Asize + j] + zxy2[i]*zxy2[j]);
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[2*Nsize*Nsize + n*Nsize + k]*betas[2*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt7[k]*oneO1PalphaSqrt7[kd]*sumsInner;
        }
      }
          P2[Nsize*n + nd] = 5.0*PI2*0.0625*sumsOuter;///pow(((double) n + 1)*((double) nd + 1),2);
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P2[Nsize*nd + n] = P2[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt7);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zxy2);
  free(kkdij);

  return P2;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP3(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P3 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P3[n] = 0.0;
  }
  return P3;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt9 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zxy4 = (double*) malloc(Asize*sizeof(double));
  double* zxy233 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

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

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[3*Nsize*Nsize + n*Nsize + k]*betas[3*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt9[k]*oneO1PalphaSqrt9[kd]*sumsInner;
        }
      }
          P3[Nsize*n + nd] = 7.0*PI2*0.0625*sumsOuter;///pow(((double) n + 1)*((double) nd + 1),3);
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P3[Nsize*nd + n] = P3[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt9);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zxy4);
  free(zxy233);
  free(kkdij);

  return P3;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP4(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P4 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P4[n] = 0.0;
  }
  return P4;
}

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
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

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

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[4*Nsize*Nsize + n*Nsize + k]*betas[4*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt11[k]*oneO1PalphaSqrt11[kd]*sumsInner;
        }
      }
          P4[Nsize*n + nd] = 9.0*PI2*0.00390625*sumsOuter;///pow(((double) n + 1)*((double) nd + 1),4);
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P4[Nsize*nd + n] = P4[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt11);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr7);
  free(zr73);
  free(zr35);
  free(kkdij);

  return P4;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP5(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P5 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P5[n] = 0.0;
  }
  return P5;
}

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
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

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
      ReXStripe[i*Asize + j] = 31.5*ReX[4*Asize*Asize + i*Asize + j] + z[i]*z[j]*(315*ReX[3*Asize*Asize + i*Asize + j] + 420*zr3[i]*zr3[j]*ReX[Asize*Asize + i*Asize + j] + zr63[i]*zr63[j]) + 17.5*zr9[i]*zr9[j]*ReX[2*Asize*Asize + i*Asize + j] + 15*zr21[i]*zr21[j]*ReX[i*Asize + j];
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[5*Nsize*Nsize + n*Nsize + k]*betas[5*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt13[k]*oneO1PalphaSqrt13[kd]*sumsInner;
        }
      }
          P5[Nsize*n + nd] = 11.0*PI2*0.00390625*sumsOuter;///pow(((double) n + 1)*((double) nd + 1),5);
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P5[Nsize*nd + n] = P5[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt13);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr9);
  free(zr3);
  free(zr21);
  free(zr63);
  free(kkdij);
  return P5;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP6(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P6 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P6[n] = 0.0;
  }
  return P6;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0;
  double r4=0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt15 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr11 = (double*) malloc(Asize*sizeof(double));
  double* zr113 = (double*) malloc(Asize*sizeof(double));
  double* zr3311 = (double*) malloc(Asize*sizeof(double));
  double* zr3330 = (double*) malloc(Asize*sizeof(double));
  double* zr231 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[6*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt15[n] = oneO1PalphaSqrt*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[6*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i];
    r4 = r2[i]*r2[i];
    zr11[i] = 11*z2 - r2[i];
    zr113[i] = 11*z2 - 3*r2[i];
    zr3311[i] = 33*z2*z2 - 18*z2*r2[i] + r4;
    zr3330[i] = 33*z2*z2 - 30*z2*r2[i] + 5*r4;
    zr231[i] = 231*z2*z2*z2 - 315*z2*z2*r2[i] + 105*z2*r4 - 5*r2[i]*r4;
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] = 115.5*ReX[5*Asize*Asize + i*Asize + j]
        + z[i]*z[j]*(
            1386*ReX[4*Asize*Asize + i*Asize + j]
            + 210*zr113[i]*zr113[j]*ReX[2*Asize*Asize + i*Asize + j]
            + 84*zr3330[i]*zr3330[j]*ReX[i*Asize + j]
            )
        + 63*zr11[i]*zr11[j]*ReX[3*Asize*Asize + i*Asize + j]
        + 52.5*zr3311[i]*zr3311[j]*ReX[Asize*Asize + i*Asize + j]
        + zr231[i]*zr231[j];
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[6*Nsize*Nsize + n*Nsize + k]*betas[6*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt15[k]*oneO1PalphaSqrt15[kd]*sumsInner;
        }
      }
          P6[Nsize*n + nd] = 13.0*PI2*9.765625e-4*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P6[Nsize*nd + n] = P6[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt15);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr11);
  free(zr113);
  free(zr3311);
  free(zr3330);
  free(zr231);
  free(kkdij);

  return P6;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP7(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P7 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P7[n] = 0.0;
  }
  return P7;
 }

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt17 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr13 = (double*) malloc(Asize*sizeof(double));
  double* zr133 = (double*) malloc(Asize*sizeof(double));
  double* zr1436 = (double*) malloc(Asize*sizeof(double));
  double* zr1431 = (double*) malloc(Asize*sizeof(double));
  double* zr4294 = (double*) malloc(Asize*sizeof(double));
  double* zr4296 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.00000000+alphas[7*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt17[n] = oneO1PalphaSqrt*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[7*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i];
    zr13[i] = 13.0*z2 - r2[i];
    zr133[i] = 13.0*z2 - 3.0*r2[i];
    zr1436[i] = 143*z2*z2 - 66*z2*r2[i] + 3*r2[i]*r2[i];
    zr1431[i] = 143*z2*z2 - 110*z2*r2[i] + 15*r2[i]*r2[i];
    zr4294[i] = 429*z2*z2*z2 - 495*z2*z2*r2[i] + 135*z2*r2[i]*r2[i] - 5*r2[i]*r2[i]*r2[i];
    zr4296[i] = 429*z2*z2*z2 - 693*z2*z2*r2[i] + 315*z2*r2[i]*r2[i] - 35*r2[i]*r2[i]*r2[i];
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] =
          107.25*ReX[6*Asize*Asize + i*Asize + j]
        + z[i]*z[j]*(1501.5*ReX[5*Asize*Asize + i*Asize + j]
        + 10.5*zr1431[i]*zr1431[j]*ReX[Asize*Asize + i*Asize + j]
        + 231*zr133[i]*zr133[j]*ReX[3*Asize*Asize + i*Asize + j]
        + zr4296[i]*zr4296[j])
        + 57.75*zr13[i]*zr13[j]*ReX[4*Asize*Asize + i*Asize + j]
        + 5.25*zr1436[i]*zr1436[j]*ReX[2*Asize*Asize + i*Asize + j]
        + 1.75*zr4294[i]*zr4294[j]*ReX[i*Asize + j];
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[7*Nsize*Nsize + n*Nsize + k]*betas[7*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt17[k]*oneO1PalphaSqrt17[kd]*sumsInner;
        }
      }
          P7[Nsize*n + nd] = 15.0*PI2*9.765625e-4*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P7[Nsize*nd + n] = P7[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt17);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr13);
  free(zr133);
  free(zr1436);
  free(zr1431);
  free(zr4294);
  free(zr4296);
  free(kkdij);

  return P7;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP8(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P8 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P8[n] = 0.0;
  }
  return P8;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0, z4 = 0, z6 = 0, z8 = 0;
  double r4=0, r6 = 0, r8 = 0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt19 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr15 = (double*) malloc(Asize*sizeof(double));
  double* zr5 = (double*) malloc(Asize*sizeof(double));
  double* zr65 = (double*) malloc(Asize*sizeof(double));
  double* zr39 = (double*) malloc(Asize*sizeof(double));
  double* zr143 = (double*) malloc(Asize*sizeof(double));
  double* zr715 = (double*) malloc(Asize*sizeof(double));
  double* zr6435 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.0 + alphas[8*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt19[n] = oneO1PalphaSqrt*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha;
    alphaO1Palpha[n] = alphas[8*Nsize + n]*oneO1Palpha;
  }


  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i]; z4 = z2*z2; z6 = z2*z4; z8=z4*z4;
    r4 = r2[i]*r2[i]; r6 = r4*r2[i]; r8 = r4*r4;

    zr5[i] = 5.0*z2 - r2[i];
    zr15[i] = 15.0*z2 - r2[i];
    zr65[i] = 65*z4 - 26*z2*r2[i] + r4;
    zr39[i] = 39*z4 - 26*z2*r2[i] + 3*r4;
    zr143[i] = 143*z6 - 143*z4*r2[i] + 33*z2*r4 - r6;
    zr715[i] = 715*z6 - 1001*z4*r2[i] + 385*z2*r4 - 35*r6;
    zr6435[i] = 6435*z8 - 12012*z6*r2[i] + 6930*z4*r4 - 1260*z2*r6 + 35*r8;
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] =
        z[i]*z[j]*(102960*ReX[6*Asize*Asize + i*Asize + j]
        + 144144*zr5[i]*zr5[j]*ReX[4*Asize*Asize + i*Asize + j]
        + 18480*zr39[i]*zr39[j]*ReX[2*Asize*Asize + i*Asize + j]
        + 144*zr715[i]*zr715[j]*ReX[i*Asize + j]
        )
        + 6435*ReX[7*Asize*Asize + i*Asize + j]
        + 3432*zr15[i]*zr15[j]*ReX[5*Asize*Asize + i*Asize + j]
        + 2772*zr65[i]*zr65[j]*ReX[3*Asize*Asize + i*Asize + j]
        + 2520*zr143[i]*zr143[j]*ReX[Asize*Asize + i*Asize + j]
        + zr6435[i]*zr6435[j];
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[8*Nsize*Nsize + n*Nsize + k]*betas[8*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt19[k]*oneO1PalphaSqrt19[kd]*sumsInner;
        }
      }
          P8[Nsize*n + nd] = 17.0*PI2*1.52587890625e-05*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P8[Nsize*nd + n] = P8[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt19);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr15);
  free(zr5);
  free(zr65);
  free(zr39);
  free(zr143);
  free(zr715);
  free(zr6435);
  free(kkdij);

  return P8;
}
//-----------------------------------------------------------
//-----------------------------------------------------------
double* getP9(double* x, double* y, double* z,double* r2, double* alphas, double* betas , int Asize, int Nsize, double* ReX){

  double* P9 = (double*) malloc(Nsize*Nsize*sizeof(double));

  if(Asize == 0){
  for(int n = 0; n < Nsize*Nsize; n++){
      P9[n] = 0.0;
  }
  return P9;
}

  double sumsInner = 0;
  double sumsOuter = 0;
  double oneO1Palpha;
  double oneO1PalphaSqrt;
  double z2=0, z4 = 0, z6 = 0, z8 = 0;
  double r4=0, r6 = 0, r8 = 0;

  double* alphaO1Palpha = (double*) malloc(Nsize*sizeof(double));
  double* oneO1PalphaSqrt21 = (double*) malloc(Nsize*sizeof(double));
  double* ReXStripe = (double*) malloc(Asize*Asize*sizeof(double));
  double* zr17 = (double*) malloc(Asize*sizeof(double));
  double* zr173 = (double*) malloc(Asize*sizeof(double));
  double* zr85 = (double*) malloc(Asize*sizeof(double));
  double* zr1710 = (double*) malloc(Asize*sizeof(double));
  double* zr2211 = (double*) malloc(Asize*sizeof(double));
  double* zr2212 = (double*) malloc(Asize*sizeof(double));
  double* zr2431 = (double*) malloc(Asize*sizeof(double));
  double* zr12155 = (double*) malloc(Asize*sizeof(double));
  double* kkdij = (double*) malloc(Nsize*Nsize*Asize*Asize*sizeof(double));

  for(int n = 0; n < Nsize; n++){
    oneO1Palpha = 1.0/(1.0 + alphas[9*Nsize + n]);
    oneO1PalphaSqrt = sqrt(oneO1Palpha);
    oneO1PalphaSqrt21[n] = oneO1PalphaSqrt*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha
      *oneO1Palpha*oneO1Palpha*oneO1Palpha*oneO1Palpha
      *oneO1Palpha;
    alphaO1Palpha[n] = alphas[9*Nsize + n]*oneO1Palpha;
  }

  for(int i = 0; i < Asize; i++){
    z2 = z[i]*z[i]; z4 = z2*z2; z6 = z2*z4; z8=z4*z4;
    r4 = r2[i]*r2[i]; r6 = r4*r2[i]; r8 = r4*r4;

    zr17[i] = 17.0*z2 - r2[i];
    zr173[i] = 17.0*z2 - 3.0*r2[i];
    zr85[i] = 85*z4 - 30*z2*r2[i] + r4;
    zr1710[i] = 17*z4 - 10*z2*r2[i] + r4;
    zr2211[i] = 221*z6 - 195*z4*r2[i] + 39*z2*r4 - r6;
    zr2212[i] = 221*z6 - 273*z4*r2[i] + 91*z2*r4 - 7*r6;
    zr2431[i] = 2431*z8 - 4004*z6*r2[i] + 2002*z4*r4 - 308*z2*r6 + 7*r8;
    zr12155[i] = 12155*z8 - 25740*z6*r2[i]
               + 18018*z4*r4 - 4620*z2*r6 + 315*r8;
  }

  for(int i = 0; i < Asize; i++){
    for(int j = 0; j < Asize; j++){
      ReXStripe[i*Asize + j] =
        z[i]*z[j]*(109395*ReX[7*Asize*Asize + i*Asize + j]
        + 17160*zr173[i]*zr173[j]*ReX[5*Asize*Asize + i*Asize + j]
        + 180180*zr1710[i]*zr1710[j]*ReX[3*Asize*Asize + i*Asize + j]
        + 3960*zr2212[i]*zr2212[j]*ReX[Asize*Asize + i*Asize + j]
        + zr12155[i]*zr12155[j]
        )
        + 6077.5*ReX[8*Asize*Asize + i*Asize + j]
        + 3217.5*zr17[i]*zr17[j]*ReX[6*Asize*Asize + i*Asize + j]
        + 2574*zr85[i]*zr85[j]*ReX[4*Asize*Asize + i*Asize + j]
        + 2310*zr2211[i]*zr2211[j]*ReX[2*Asize*Asize + i*Asize + j]
        + 45*zr2431[i]*zr2431[j]*ReX[i*Asize + j];
    }
  }

  for(int k = 0; k < Nsize; k++){
    for(int kd = 0; kd < Nsize; kd++){
      for(int i = 0; i < Asize; i++){
        for(int j = 0; j < Asize; j++){
          kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j]
              = ReXStripe[i*Asize + j]*exp(-alphaO1Palpha[k]*r2[i]-alphaO1Palpha[kd]*r2[j]);
        }
      }
    }
  }
//#pragma omp parallel for schedule(static,chunk)
  for(int n = 0; n < Nsize; n++){
    for(int nd = n; nd < Nsize; nd++){
      sumsOuter = 0;
      for(int k = 0; k < Nsize; k++){
        for(int kd = 0; kd < Nsize; kd++){
          sumsInner = 0;
          for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++){
              sumsInner += kkdij[ k*Nsize*Asize*Asize  + kd*Asize*Asize   +  i*Asize + j];
            }
          }
          sumsOuter +=  betas[9*Nsize*Nsize + n*Nsize + k]*betas[9*Nsize*Nsize+nd*Nsize + kd]*oneO1PalphaSqrt21[k]*oneO1PalphaSqrt21[kd]*sumsInner;
        }
      }
          P9[Nsize*n + nd] = 19.0*PI2*1.52587890625e-05*sumsOuter;
    }
  }
  for(int n = 0; n < Nsize; n++){
    for(int nd = n + 1; nd < Nsize; nd++){
          P9[Nsize*nd + n] = P9[Nsize*n + nd];
    }
  }

  free(oneO1PalphaSqrt21);
  free(alphaO1Palpha);
  free(ReXStripe);
  free(zr17);
  free(zr173);
  free(zr85);
  free(zr1710);
  free(zr2211);
  free(zr2212);
  free(zr2431);
  free(zr12155);
  free(kkdij);

  return P9;

}
//-----------------------------------------------------------
//-----------------------------------------------------------
void printP(double* P, int Nsize){
  for(int i = 0; i < Nsize ; i++){
    for(int j = 0; j < Nsize ; j++){
      cout << P[i*Nsize + j] << " ";
    }
  }
}
//-----------------------------------------------------------
//-----------------------------------------------------------
int main(int argc, char* argv[]) {

  std::cout.precision(12);

  int*  totalAN = (int*) malloc(sizeof(int));
  int*  Ntypes = (int*) malloc(sizeof(int));
  int Nsize = 5;
  double* x;
  double* y;
  double* z;

  int l = atoi(argv[1]);

  int* types;

  double* alphas = getAlphas(Nsize); double* betas = getBetas(Nsize);

  int* typeNs = getInfo(totalAN, Ntypes);
  x = (double*) malloc(sizeof(double)*totalAN[0]);
  y = (double*) malloc(sizeof(double)*totalAN[0]);
  z = (double*) malloc(sizeof(double)*totalAN[0]);
  types = (int*) malloc(sizeof(int)*Ntypes[0]);

  int Asize = 0;
  double* Apos = getApos(totalAN, Ntypes, typeNs, types);
  double* Hpos = getHpos(9999);

  double* P0;
  double* P1; double* P2; double* P3;
  double* P4; double* P5; double* P6;
  double* P7; double* P8; double* P9;
  double* r2;
  double* ReX;

  for(int i = 0; i < 9999; i++){
    for(int j = 0; j < Ntypes[0]; j++){

      Asize = getFilteredPos(x, y, z, Apos, Hpos, typeNs, i, j);
      r2 = getR2(x, y, z, Asize);
      ReX = getReals(x, y, Asize);

      P0 = getP0(x,y,z,r2,alphas, betas ,Asize, Nsize);
      printP(P0, Nsize);

     if(l > 0){
       P1 = getP1(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P1, Nsize);
     }

     if(l > 1){
       P2 = getP2(x,y,z,r2,alphas, betas ,Asize, Nsize, ReX);
       printP(P2, Nsize);
     }

     if(l > 2){
       P3 = getP3(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P3, Nsize);
     }

     if(l > 3){
       P4 = getP4(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P4, Nsize);
     }

     if(l > 4){
       P5 = getP5(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P5, Nsize);
     }

     if(l > 5){
       P6 = getP6(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P6, Nsize);
     }

     if(l > 6){
       P7 = getP7(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P7, Nsize);
     }

     if(l > 7){
       P8 = getP8(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P8, Nsize);
     }

     if(l > 8){
       P9 = getP9(x,y,z,r2,alphas, betas ,Asize, Nsize,ReX);
       printP(P9, Nsize);
     }
    }
    cout << endl;
  }
}
