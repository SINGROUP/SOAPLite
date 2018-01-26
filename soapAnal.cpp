#include<iostream>
//#include<armadillo>
//#include"../Header/myArmadillo.h"
//#include"../Header/myMath.h"
//#include"../Header/mySoap.h"

using namespace std;
//using namespace arma;


double* getReals(double* x, double* y, int size) {
double P,M;
double P2,M2;
double P4,M4;
double P8,M8;

for(int i = 0; i < size; i++){
  for(int j = 0; j < size; j++){
    P[i*size + j] =x[i]*x[j] + y[i]*y[j];
    M[i*size+j] =x[i]*y[j] - y[i]*x[j]; 

    P2[i*size+j] =P[i*size+j]*P[i*size+j];     M2[i*size+j] =M[i*size+j]*M[i*size+j]; 
    P4[i*size+j] =P2[i*size+j]*P2[i*size+j];   M4[i*size+j] =M2[i*size+j]*M2[i*size+j]; 
    P6[i*size+j] =P4[i*size+j]*P2[i*size+j];   M6[i*size+j] =M4[i*size+j]*M2[i*size+j]; 
    P8[i*size+j] =P4[i*size+j]*P4[i*size+j];   M8[i*size+j] =M4[i*size+j]*M4[i*size+j]; 
}

Re[0] = P[i*size+j];
Re[1] = P[i*size+j] - M[i*size+j];
Re[2] = P[i*size+j]*(P2[i*size+j] - 3*M2[i*size+j]);
Re[3] = P4[i*size+j] + M4[i*size+j] - 6*P2[i*size+j]*M2[i*size+j];
Re[4] = P[i*size+j]*(P4[i*size+j] - 10*P2[i*size+j]*M2[i*size+j] + 5*M4[i*size+j]);
Re[5] = P6[i*size+j] - M6[i*size+j]  - 15*(P4[i*size+j]*M2[i*size+j] + P2[i*size+j]*M4[i*size+j]);
Re[6] = P[i*size+j]*(P6[i*size+j] - 7*M6[i*size+j] - 21*P4[i*size+j]*M2[i*size+j] + 35*P2[i*size+j]*M4[i*size+j]);
Re[7] = P8[i*size+j] + M8[i*size+j] - 28*(P6[i*size+j]*M2[i*size+j] + P2[i*size+j]*M6[i*size+j]) + 70*P4[i*size+j]*M4[i*size+j];
Re[8] = P[i*size+j]*(P8[i*size+j] + 9*M8[i*size+j] - 36*P6[i*size+j]*M2[i*size+j] + 126*P4[i*size+j]*M4[i*size+j] - 84*P2[i*size+j]*M6[i*size+j]);

return 0;
}
