#include<iostream>
#include<stdio.h>
#include <stdlib.h>     /* malloc, free, rand */
#include <math.h>     /* malloc, free, rand */
#include <stdlib.h>     /* malloc, free, rand */

using namespace std;
//using namespace arma;

#define PI 3.14159265358979
#define PI2 6.28318530717959
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
    oneO1Palpha = 1.0/(1.00000000123+alphas[n]);
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
int main(int argc, char* argv[]) {
int Asize = 2;
int Nsize = 3;
double* x = (double*) malloc(Asize*sizeof(double));
double* y = (double*) malloc(Asize*sizeof(double));
double* z = (double*) malloc(Asize*sizeof(double));
double* alphas = getAlphas();
double* betas = getBetas();

x[0] = 1.0;
x[1] = -1.0;
y[0] = 1.0;
y[1] = -1.0;
z[0] = 1.0;
z[1] = -1.0;

//cout << "AA"<< endl;

double* r2 = getR2(x, y, z, Asize);
//cout << "AB"<< endl;
//cout << "AC"<< endl;

double* f = getReals(x, y,Asize); 
//cout << "AD"<< endl;

//cout << f[8*Asize*Asize + 0*Asize + 0] << endl;


//for(int i = 0; i < 50; i++){
//  printf ("I have : %lf \n",alphas[i]);
//}
//for(int i = 0; i < 5*5*10; i++){
//  printf ("I have Beta : %lf \n",betas[i]);
//}

cout << betas[5*5*9 + 5*3 + 4] << endl;
cout << alphas[5*9 + 4] << endl;


double* P0 = getP0(x,y,z,r2,alphas, betas ,Asize, Nsize);


cout << P0[0*Nsize*Asize + 0*Nsize + 0] << endl;
cout << f[8*Asize*Asize + 0*Asize + 0] << endl;

free(x);
free(y);
free(f);

}
