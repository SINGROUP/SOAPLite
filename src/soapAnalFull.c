#include<stdio.h>
#include<iostream>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define PI2 9.86960440108936
#define PI 3.14159265359
#define PIHalf 1.57079632679490
//===========================================================
double* getReIm2(double* x, double* y, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*x[i]-y[i]*y[i];
    c3[2*i+1] = 2*y[i]*x[i]; 
  }
}
//===========================================================
void getReIm3(double* x, double* y, double* c2, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*c2[2*i] - y[i]*c2[2*i + 1];
    c3[2*i+1] = x[i]*c2[2*i+1] + y[i]*c2[2*i  ]; 
  }
}
//===========================================================
void getMulReIm(double* c1, double* c2, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = c1[2*i  ]*c2[2*i  ] - c1[2*i+1]*c2[2*i+1];
    c3[2*i+1] = c1[2*i  ]*c2[2*i+1] + c1[2*i+1]*c2[2*i  ]; 
  }
}
//===========================================================
void getMulDouble(double* c1, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = c1[2*i]*c1[2*i] - c1[2*i+1]*c1[2*i+1];
    c3[2*i+1] = 2*c1[2*i]*c1[2*i+1]; 
  }
}
//========================================
double* getAlphas(int alphaSize){
  double f;
  double* alphas = (double*) malloc(alphaSize*10*sizeof(double));
  FILE * pFile;
  pFile = fopen ("alphas.dat","r");
  for(int i = 0; i < alphaSize*10; i++){
    fscanf (pFile, "%lf", &alphas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return alphas;
}
//========================================
double* getBetas(int alphaSize){
  double f;
  double* betas = (double*) malloc(alphaSize*alphaSize*10*sizeof(double));
  FILE * pFile;
  pFile = fopen ("betas.dat","r");
  for(int i = 0; i < alphaSize*alphaSize*10; i++){
    fscanf (pFile, "%lf", &betas[i]);
    //  rewind (pFile);
  }
  fclose (pFile);
  return betas;
}
//========================================
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
//========================================
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
//========================================
double* getHpos(int Hsize,char* ch){
  FILE* pFile;
  pFile = fopen(ch,"r");
//  pFile = fopen("HRot.dat","r");
  double* Pos = (double*) malloc(3*sizeof(double)*Hsize);
  for(int i=0; i < Hsize; i++){
      fscanf(pFile, "%lf", &Pos[3*i    ]);
      fscanf(pFile, "%lf", &Pos[3*i + 1]);
      fscanf(pFile, "%lf", &Pos[3*i + 2]);
  }
  fclose(pFile);
  return Pos;
}
//========================================
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
//========================================
int getAllPos(double* x, double* y, double* z, double* Apos, double* Hpos, int* typeNs, int Ihpos,int sizeAll){
  int count = 0;
  int X, Y, Z;
  for(int i = 0; i < sizeAll; i++){
      X =  Apos[3*i    ] - Hpos[3*Ihpos    ];
      Y =  Apos[3*i + 1] - Hpos[3*Ihpos + 1];
      Z =  Apos[3*i + 2] - Hpos[3*Ihpos + 2];
      if( X*X + Y*Y + Z*Z < 100.0 ){
        x[count] = X;
        y[count] = Y;
        z[count] = Z;
        count++;
      }
  }
  return count;
}
//================================================================
int getFilteredPos(double* x, double* y, double* z, double* Apos, double* Hpos,
      int* typeNs, double rCutSqr, int Ihpos, int Itype){

  int shiftType = 0; int count = 0;
  double X = 0; double Y = 0; double Z = 0;

    for(int i = 0; i < Itype ; i++){
      shiftType += typeNs[i];
    }

    for(int i = 0; i < typeNs[Itype]; i++){
      X = Apos[3*shiftType + 3*i    ] - Hpos[3*Ihpos    ];
      Y = Apos[3*shiftType + 3*i + 1] - Hpos[3*Ihpos + 1];
      Z = Apos[3*shiftType + 3*i + 2] - Hpos[3*Ihpos + 2];
      if( X*X + Y*Y + Z*Z < rCutSqr ){
        x[count] = X;
        y[count] = Y;
        z[count] = Z;
        count++;
      }
    }
  return count;
}
//================================================================
double* getRsZs(double* x, double* y, double* z,double* r2,double* r4,double* r6,double* r8,double* z2,double* z4,double* z6,double* z8, int size){
  for(int i = 0; i < size; i++){
    r2[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
    r4[i] = r2[i]*r2[i]; r6[i] = r2[i]*r4[i]; r8[i] = r4[i]*r4[i];
    z2[i] = z[i]*z[i]; z4[i] = z2[i]*z2[i]; z6[i] = z2[i]*z4[i]; z8[i] = z4[i]*z4[i];
  }
}
//================================================================
void getAlphaBeta(double* aOa, double* bOa, double* alphas, double* betas, int Ns){
  int lMax = 9;

  int  NsNs = Ns*Ns;
  double  oneO1alpha;      double  oneO1alpha2; double  oneO1alpha3;
  double  oneO1alpha4; double  oneO1alpha5; double  oneO1alpha6;
  double  oneO1alpha7; double  oneO1alpha8; double  oneO1alpha9;
  double  oneO1alpha10;
  double  oneO1alphaSqrt;// = (double*) malloc(Ns*sizeof(double));
  double  oneO1alphaSqrtX;

  // MY POEWR MISSING (see beggning functions);

  for(int k = 0; k < Ns; k++){
    oneO1alpha = 1.0/(1.0 + alphas[k]);
    oneO1alphaSqrt = sqrt(oneO1alpha);
    aOa[k] = -alphas[k]*oneO1alpha; //got alpha_0k
    oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
    for(int n = 0; n < Ns; n++){ bOa[n*Ns + k] = betas[n*Ns + k]*oneO1alphaSqrtX;} // got beta_0nk
  }
  if(lMax > 0){
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[Ns + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[Ns + k] = -alphas[Ns + k]*oneO1alpha; //got alpha_1k
      oneO1alpha2 = oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha2;
      for(int n = 0; n < Ns; n++){ bOa[NsNs + n*Ns + k] = betas[NsNs + n*Ns + k]*oneO1alphaSqrtX;} // got beta_1nk
    }
  } if(lMax > 1){
    int shift1 = 2*Ns; int shift2 = 2*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_2k
      oneO1alpha3 = oneO1alpha*oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha3;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_2nk
    }
  } if(lMax > 2){
    int shift1 = 3*Ns; int shift2 = 3*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_3k
      oneO1alpha4 = oneO1alpha*oneO1alpha*oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha4;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_3nk
    }
  } if(lMax > 3){
    int shift1 = 4*Ns; int shift2 = 4*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_4k
      oneO1alpha5 = pow(oneO1alpha,5); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha5;
      for(int n = 0; n < Ns; n++){ bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_4nk
    }
  } if(lMax > 4){
    int shift1 = 5*Ns; int shift2 = 5*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_5k
      oneO1alpha6 = pow(oneO1alpha,6); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha6;
      for(int n = 0; n < Ns; n++){ bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_5nk
    }
  } if(lMax > 5){
    int shift1 = 6*Ns; int shift2 = 6*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_6k
      oneO1alpha7 = pow(oneO1alpha,7); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha7;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_6nk
    }
  } if(lMax > 6){
    int shift1 = 7*Ns; int shift2 = 7*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_7k
      oneO1alpha8 = pow(oneO1alpha,8); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha8;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_7nk
    }
  } if(lMax > 7){
    int shift1 = 8*Ns; int shift2 = 8*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_8k
      oneO1alpha9 = pow(oneO1alpha,9); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha9;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_8nk
    }
  }if(lMax > 8){
    int shift1 = 9*Ns; int shift2 = 9*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_9k
      oneO1alpha10 = pow(oneO1alpha,10); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha10;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_9nk
    }
  }
}
//================================================================
void getCfactors(double* preCoef, int Asize, double* x, double* y, double* z, double* z2, double* z4, double* z6, double* z8, double* r2, double* r4, double* r6, double* r8, double* ReIm2, double* ReIm3, double* ReIm4, double* ReIm5, double* ReIm6, double* ReIm7, double* ReIm8, double* ReIm9,int totalAN){
  double c20c;double c30c;double c31c;double c40c;double c41c;double c42c;
  double c50c;double c51c;double c52c;double c53c;double c60c;double c61c;
  double c62c;double c63c;double c64c;double c70c;double c71c;double c72c;
  double c73c;double c74c;double c75c;double c80c;double c81c;double c82c;
  double c83c;double c84c;double c85c;double c86c;double c90c;double c91c;
  double c92c;double c93c;double c94c;double c95c;double c96c;double c97c;

    getReIm2(x, y, ReIm2,Asize);
    getReIm3(x, y, ReIm2, ReIm3, Asize);
    getMulDouble(ReIm2, ReIm4, Asize);
    getMulReIm(ReIm2,ReIm3, ReIm5, Asize);
    getMulDouble(ReIm3, ReIm6, Asize);
    getMulReIm(ReIm3,ReIm4, ReIm7, Asize);
    getMulDouble(ReIm4, ReIm8, Asize);
    getMulReIm(ReIm4,ReIm5, ReIm9, Asize);

  for(int i = 0; i < Asize; i++){
    c20c=3*z2[i]-r2[i];
    c30c=5*z2[i]-3*r2[i];
    c31c=5*z2[i]-r2[i];
    c40c=35*z4[i]-30*z2[i]*r2[i]+3*r4[i];
    c41c=7*z2[i]-3*r2[i];
    c42c=7*z2[i]-r2[i];
    c50c=63*z4[i]-70*z2[i]*r2[i]+15*r4[i];
    c51c=21*z4[i]-14*z2[i]*r2[i]+r4[i];
    c52c=3*z2[i]-r2[i];
    c53c=9*z2[i]-r2[i];
    c60c=231*z6[i] - 315*z4[i]*r2[i] + 105*z2[i]*r4[i] - 5*r6[i];
    c61c=33*z4[i] - 30*z2[i]*r2[i] + 5*r4[i];
    c62c=33*z4[i] - 18*z2[i]*r2[i] + r4[i];
    c63c=11*z2[i] - 3*r2[i];
    c64c=11*z2[i] - r2[i];
    c70c=429*z6[i]-693*z4[i]*r2[i]+315*z2[i]*r4[i]-35*r6[i];
    c71c=429*z6[i]-495*z4[i]*r2[i]+135*z2[i]*r4[i]-5*r6[i];
    c72c=143*z4[i]-110*z2[i]*r2[i]+15*r4[i];
    c73c=143*z4[i]-66*z2[i]*r2[i]+3*r4[i];
    c74c=13*z2[i]-3*r2[i];
    c75c=13*z2[i]-r2[i];
    c80c=6435*z8[i]-12012*z6[i]*r2[i]+6930*z4[i]*r4[i]-1260*z2[i]*r6[i]+35*r8[i];
    c81c=715*z6[i]-1001*z4[i]*r2[i]+385*z2[i]*r4[i]-35*r6[i];
    c82c=143*z6[i]-143*z4[i]*r2[i]+33*z2[i]*r4[i]-r6[i];
    c83c=39*z4[i]-26*z2[i]*r2[i]+3*r4[i];
    c84c=65*z4[i]-26*z2[i]*r2[i]+r4[i];
    c85c=5*z2[i]-r2[i];
    c86c=15*z2[i]-r2[i];
    c90c=12155*z8[i]-25740*z6[i]*r2[i]+18018*z4[i]*r4[i] -4620*z2[i]*r6[i]+315*r8[i];
    c91c=2431*z8[i]-4004*z6[i]*r2[i]+2002*z4[i]*r4[i]-308*z2[i]*r6[i] + 7*r8[i];
    c92c=221*z6[i]-273*z4[i]*r2[i]+91*z2[i]*r4[i]-7*r6[i];
    c93c=221*z6[i]-195*z4[i]*r2[i]+39*z2[i]*r4[i]-r6[i];
    c94c=17*z4[i]-10*z2[i]*r2[i]+r4[i];
    c95c=85*z4[i]-30*z2[i]*r2[i]+r4[i];
    c96c=17*z2[i]-3*r2[i];
    c97c=17*z2[i]-r2[i];
    /*c20  */  preCoef[          +i] = c20c;
    /*c21Re*/  preCoef[totalAN   +i] = z[i]*x[i];
    /*c21Im*/  preCoef[totalAN*2 +i] = z[i]*y[i];
    /*c22Re*/  preCoef[totalAN*3 +i] =      ReIm2[2*i];
    /*c22Im*/  preCoef[totalAN*4 +i] =      ReIm2[2*i+1];
    /*c30  */  preCoef[totalAN*5 +i] = c30c*z[i]; 
    /*c31Re*/  preCoef[totalAN*6 +i] =       x[i]*c31c;
    /*c31Im*/  preCoef[totalAN*7 +i] =       y[i]*c31c;
    /*c32Re*/  preCoef[totalAN*8 +i] = z[i]*ReIm2[2*i];
    /*c32Im*/  preCoef[totalAN*9 +i] = z[i]*ReIm2[2*i+1];
    /*c33Re*/  preCoef[totalAN*10+i] =       ReIm3[2*i  ];
    /*c33Im*/  preCoef[totalAN*11+i] =       ReIm3[2*i+1];
    /*c40  */  preCoef[totalAN*12+i] = c40c;
    /*c41Re*/  preCoef[totalAN*13+i] = z[i]*x[i]*c41c;
    /*c41Im*/  preCoef[totalAN*14+i] = z[i]*y[i]*c41c;
    /*c42Re*/  preCoef[totalAN*15+i] =      ReIm2[2*i  ]*c42c;
    /*c42Im*/  preCoef[totalAN*16+i] =      ReIm2[2*i+1]*c42c;
    /*c43Re*/  preCoef[totalAN*17+i] = z[i]*ReIm3[2*i  ];
    /*c43Im*/  preCoef[totalAN*18+i] = z[i]*ReIm3[2*i+1];
    /*c44Re*/  preCoef[totalAN*19+i] =      ReIm4[2*i  ];
    /*c44Im*/  preCoef[totalAN*20+i] =      ReIm4[2*i+1];
    /*c50  */  preCoef[totalAN*21+i] = c50c*z[i];
    /*c51Re*/  preCoef[totalAN*22+i] =      x[i]*c51c;
    /*c51Im*/  preCoef[totalAN*23+i] =      y[i]*c51c;
    /*c52Re*/  preCoef[totalAN*24+i] = z[i]*ReIm2[2*i  ]*c52c;
    /*c52Im*/  preCoef[totalAN*25+i] = z[i]*ReIm2[2*i+1]*c52c;
    /*c53Re*/  preCoef[totalAN*26+i] =      ReIm3[2*i  ]*c53c;
    /*c53Im*/  preCoef[totalAN*27+i] =      ReIm3[2*i+1]*c53c;
    /*c54Re*/  preCoef[totalAN*28+i] = z[i]*ReIm4[2*i  ];
    /*c54Im*/  preCoef[totalAN*29+i] = z[i]*ReIm4[2*i+1]; 
    /*c55Re*/  preCoef[totalAN*30+i] =      ReIm5[2*i  ];
    /*c55Im*/  preCoef[totalAN*31+i] =      ReIm5[2*i+1];
    /*c60  */  preCoef[totalAN*32+i] = c60c;
    /*c61Re*/  preCoef[totalAN*33+i] = z[i]*x[i]*c61c;
    /*c61Im*/  preCoef[totalAN*34+i] = z[i]*y[i]*c61c;
    /*c62Re*/  preCoef[totalAN*35+i] =      ReIm2[2*i  ]*c62c;
    /*c62Im*/  preCoef[totalAN*36+i] =      ReIm2[2*i+1]*c62c;
    /*c63Re*/  preCoef[totalAN*37+i] = z[i]*ReIm3[2*i  ]*c63c;
    /*c63Im*/  preCoef[totalAN*38+i] = z[i]*ReIm3[2*i+1]*c63c;
    /*c64Re*/  preCoef[totalAN*39+i] =      ReIm4[2*i  ]*c64c;
    /*c64Im*/  preCoef[totalAN*40+i] =      ReIm4[2*i+1]*c64c;
    /*c65Re*/  preCoef[totalAN*41+i] = z[i]*ReIm5[2*i  ];
    /*c65Im*/  preCoef[totalAN*42+i] = z[i]*ReIm5[2*i+1];
    /*c66Re*/  preCoef[totalAN*43+i] =      ReIm6[2*i  ];
    /*c66Im*/  preCoef[totalAN*44+i] =      ReIm6[2*i+1];
    /*c70  */  preCoef[totalAN*45+i] = c70c*z[i];
    /*c71Re*/  preCoef[totalAN*46+i] = x[i]*c71c;
    /*c71Im*/  preCoef[totalAN*47+i] = y[i]*c71c;
    /*c72Re*/  preCoef[totalAN*48+i] = z[i]*ReIm2[2*i  ]*c72c;
    /*c72Im*/  preCoef[totalAN*49+i] = z[i]*ReIm2[2*i+1]*c72c;
    /*c73Re*/  preCoef[totalAN*50+i] =      ReIm3[2*i  ]*c73c;
    /*c73Im*/  preCoef[totalAN*51+i] =      ReIm3[2*i+1]*c73c;
    /*c74Re*/  preCoef[totalAN*52+i] = z[i]*ReIm4[2*i  ]*c74c;
    /*c74Im*/  preCoef[totalAN*53+i] = z[i]*ReIm4[2*i+1]*c74c;
    /*c75Re*/  preCoef[totalAN*54+i] =      ReIm5[2*i  ]*c75c;
    /*c75Im*/  preCoef[totalAN*55+i] =      ReIm5[2*i+1]*c75c;
    /*c76Re*/  preCoef[totalAN*56+i] = z[i]*ReIm6[2*i  ];
    /*c76Im*/  preCoef[totalAN*57+i] = z[i]*ReIm6[2*i+1];
    /*c77Re*/  preCoef[totalAN*58+i] =      ReIm7[2*i  ];
    /*c77Im*/  preCoef[totalAN*59+i] =      ReIm7[2*i+1];
    /*c80  */  preCoef[totalAN*60+i] = c80c;
    /*c81Re*/  preCoef[totalAN*61+i] = z[i]*x[i]*c81c;
    /*c81Im*/  preCoef[totalAN*62+i] = z[i]*y[i]*c81c;
    /*c82Re*/  preCoef[totalAN*63+i] =      ReIm2[2*i  ]*c82c; // ??
    /*c82Im*/  preCoef[totalAN*64+i] =      ReIm2[2*i+1]*c82c; // ??
    /*c83Re*/  preCoef[totalAN*65+i] = z[i]*ReIm3[2*i  ]*c83c; // ??
    /*c83Im*/  preCoef[totalAN*66+i] = z[i]*ReIm3[2*i+1]*c83c; // ??
    /*c84Re*/  preCoef[totalAN*67+i] =      ReIm4[2*i  ]*c84c; // ??
    /*c84Im*/  preCoef[totalAN*68+i] =      ReIm4[2*i+1]*c84c; // ??
    /*c85Re*/  preCoef[totalAN*69+i] = z[i]*ReIm5[2*i  ]*c85c; // ??
    /*c85Im*/  preCoef[totalAN*70+i] = z[i]*ReIm5[2*i+1]*c85c; // ??
    /*c86Re*/  preCoef[totalAN*71+i] =      ReIm6[2*i  ]*c86c; // ??
    /*c86Im*/  preCoef[totalAN*72+i] =      ReIm6[2*i+1]*c86c; // ??
    /*c87Re*/  preCoef[totalAN*73+i] = z[i]*ReIm7[2*i  ];
    /*c87Im*/  preCoef[totalAN*74+i] = z[i]*ReIm7[2*i+1];
    /*c88Re*/  preCoef[totalAN*75+i] =      ReIm8[2*i  ];
    /*c88Im*/  preCoef[totalAN*76+i] =      ReIm8[2*i+1];
    /*c90  */  preCoef[totalAN*77+i] = c90c*z[i];
    /*c91Re*/  preCoef[totalAN*78+i] = x[i]*c91c;
    /*c91Im*/  preCoef[totalAN*79+i] = y[i]*c91c;
    /*c92Re*/  preCoef[totalAN*80+i] = z[i]*ReIm2[2*i  ]*c92c;
    /*c92Im*/  preCoef[totalAN*81+i] = z[i]*ReIm2[2*i+1]*c92c;
    /*c93Re*/  preCoef[totalAN*82+i] =      ReIm3[2*i  ]*c93c;
    /*c93Im*/  preCoef[totalAN*83+i] =      ReIm3[2*i+1]*c93c;
    /*c94Re*/  preCoef[totalAN*84+i] = z[i]*ReIm4[2*i  ]*c94c;
    /*c94Im*/  preCoef[totalAN*85+i] = z[i]*ReIm4[2*i+1]*c94c;
    /*c95Re*/  preCoef[totalAN*86+i] =      ReIm5[2*i  ]*c95c;
    /*c95Im*/  preCoef[totalAN*87+i] =      ReIm5[2*i+1]*c95c;
    /*c96Re*/  preCoef[totalAN*88+i] = z[i]*ReIm6[2*i  ]*c96c;
    /*c96Im*/  preCoef[totalAN*89+i] = z[i]*ReIm6[2*i+1]*c96c;
    /*c97Re*/  preCoef[totalAN*90+i] =      ReIm7[2*i  ]*c97c;
    /*c97Im*/  preCoef[totalAN*91+i] =      ReIm7[2*i+1]*c97c;
    /*c98Re*/  preCoef[totalAN*92+i] = z[i]*ReIm8[2*i  ];
    /*c98Im*/  preCoef[totalAN*93+i] = z[i]*ReIm8[2*i+1];
    /*c99Re*/  preCoef[totalAN*94+i] =      ReIm9[2*i  ];
    /*c99Im*/  preCoef[totalAN*95+i] =      ReIm9[2*i+1];
  }
}
//================================================================
int getC(double* C, double* preCoef, double* x, double* y, double* z,double* r2, double* bOa, double* aOa, double* exes, int Hsize, int totalAN, int Asize, int Ns, int Ntypes, int lMax, int posI, int typeJ){

  if(Asize == 0){return 0;}
  double sumMe = 0; int NsNs = Ns*Ns;  int NsJ = 100*Ns*typeJ; int LNsNs;
  int Nx2 = 2*Ns; int Nx3 = 3*Ns; int Nx4 = 4*Ns; int Nx5 = 5*Ns;
  int Nx6 = 6*Ns; int Nx7 = 7*Ns; int Nx8 = 8*Ns; int Nx9 = 9*Ns;
  int Nx10 = 10*Ns; int Nx11 = 11*Ns; int Nx12 = 12*Ns; int Nx13 = 13*Ns;
  int Nx14 = 14*Ns; int Nx15 = 15*Ns; int Nx16 = 16*Ns; int Nx17 = 17*Ns;
  int Nx18 = 18*Ns; int Nx19 = 19*Ns; int Nx20 = 20*Ns; int Nx21 = 21*Ns;
  int Nx22 = 22*Ns; int Nx23 = 23*Ns; int Nx24 = 24*Ns; int Nx25 = 25*Ns;
  int Nx26 = 26*Ns; int Nx27 = 27*Ns; int Nx28 = 28*Ns; int Nx29 = 29*Ns;
  int Nx30 = 30*Ns; int Nx31 = 31*Ns; int Nx32 = 32*Ns; int Nx33 = 33*Ns;
  int Nx34 = 34*Ns; int Nx35 = 35*Ns; int Nx36 = 36*Ns; int Nx37 = 37*Ns;
  int Nx38 = 38*Ns; int Nx39 = 39*Ns; int Nx40 = 40*Ns; int Nx41 = 41*Ns;
  int Nx42 = 42*Ns; int Nx43 = 43*Ns; int Nx44 = 44*Ns; int Nx45 = 45*Ns;
  int Nx46 = 46*Ns; int Nx47 = 47*Ns; int Nx48 = 48*Ns; int Nx49 = 49*Ns;
  int Nx50 = 50*Ns; int Nx51 = 51*Ns; int Nx52 = 52*Ns; int Nx53 = 53*Ns;
  int Nx54 = 54*Ns; int Nx55 = 55*Ns; int Nx56 = 56*Ns; int Nx57 = 57*Ns;
  int Nx58 = 58*Ns; int Nx59 = 59*Ns; int Nx60 = 60*Ns; int Nx61 = 61*Ns;
  int Nx62 = 62*Ns; int Nx63 = 63*Ns; int Nx64 = 64*Ns; int Nx65 = 65*Ns;
  int Nx66 = 66*Ns; int Nx67 = 67*Ns; int Nx68 = 68*Ns; int Nx69 = 69*Ns;
  int Nx70 = 70*Ns; int Nx71 = 71*Ns; int Nx72 = 72*Ns; int Nx73 = 73*Ns;
  int Nx74 = 74*Ns; int Nx75 = 75*Ns; int Nx76 = 76*Ns; int Nx77 = 77*Ns;
  int Nx78 = 78*Ns; int Nx79 = 79*Ns; int Nx80 = 80*Ns; int Nx81 = 81*Ns;
  int Nx82 = 82*Ns; int Nx83 = 83*Ns; int Nx84 = 84*Ns; int Nx85 = 85*Ns;
  int Nx86 = 86*Ns; int Nx87 = 87*Ns; int Nx88 = 88*Ns; int Nx89 = 89*Ns;
  int Nx90 = 90*Ns; int Nx91 = 91*Ns; int Nx92 = 92*Ns; int Nx93 = 93*Ns;
  int Nx94 = 94*Ns; int Nx95 = 95*Ns; int Nx96 = 96*Ns; int Nx97 = 97*Ns;
  int Nx98 = 98*Ns; int Nx99 = 99*Ns;
  int LNs; int NsTsI = 100*Ns*Ntypes*posI;
  for(int k = 0; k < Ns; k++){
    sumMe = 0; for(int i = 0; i < Asize; i++){ sumMe += exp(aOa[k]*r2[i]);}
    for(int n = 0; n < Ns; n++){ C[NsTsI + NsJ + n] += bOa[n*Ns + k]*sumMe; }
  } if(lMax > 0) { LNsNs=NsNs; LNs=Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c10*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*z[i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns + n] += bOa[LNsNs + n*Ns + k]*sumMe;} 
      sumMe = 0;/*c11Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*x[i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx2 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c11Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*y[i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx3 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }} if(lMax > 1) { LNsNs=2*NsNs; LNs=2*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c20*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx4 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c21Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*1 + i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx5 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c21Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*2 + i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx6 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c22Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*3 + i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx7 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c22Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*4 + i];}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx8 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    }} if(lMax > 2) { LNsNs=3*NsNs; LNs=3*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c30*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*5+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx9 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c31Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*6+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx10 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c31Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*7+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx11 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c32Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*8+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx12 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c32Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*9+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx13 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c33Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*10+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx14 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c33Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*11+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx15 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    }} if(lMax > 3) { LNsNs=4*NsNs; LNs=4*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c40*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*12+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx16 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c41Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*13+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx17 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c41Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*14+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx18 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c42Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*15+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx19 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c42Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*16+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx20 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c43Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*17+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx21 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c43Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*18+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx22 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c44Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*19+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx23 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
      sumMe = 0;/*c44Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*20+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx24 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    }} if(lMax > 4) { LNsNs=5*NsNs; LNs=5*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c50*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*21+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx25 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c51Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*22+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx26 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c51Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*23+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx27 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c52Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*24+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx28 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c52Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*25+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx29 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c53Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*26+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx30 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c53Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*27+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx31 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c54Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*28+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx32 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c54Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*29+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx33 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c55Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*30+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx34 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c55Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*31+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx35 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }} if(lMax > 5) { LNsNs=6*NsNs; LNs=6*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c60*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*32+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx36 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c61Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*33+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx37 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c61Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*34+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx38 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c62Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*35+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx39 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c62Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*36+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx40 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c63Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*37+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx41 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c63Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*38+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx42 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c64Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*39+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx43 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c64Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*40+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx44 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c65Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*41+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx45 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c65Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*42+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx46 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c66Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*43+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx47 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c66Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*44+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx48 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }} if(lMax > 6) { LNsNs=7*NsNs; LNs=7*Ns; 
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c70*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*45+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx49 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c71Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*46+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx50 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c71Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*47+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx51 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c72Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*48+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx52 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c72Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*49+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx53 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c73Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*50+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx54 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c73Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*51+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx55 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c74Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*52+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx56 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c74Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*53+i]);} 
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx57 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c75Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*54+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx58 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c75Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*55+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx59 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c76Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*56+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx60 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c76Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*57+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx61 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c77Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*58+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx62 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c77Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*59+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx63 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }} if(lMax > 7) { LNsNs=8*NsNs; LNs=8*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c80*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*60+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx64 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c81Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*61+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx65 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c81Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*62+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx66 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c82Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*63+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx67 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c82Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*64+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx68 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c83Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*65+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx69 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c83Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*66+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx70 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c84Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*67+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx71 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c84Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*68+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx72 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c85Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*69+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx73 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c85Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*70+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx74 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c86Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*71+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx75 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c86Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*72+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx76 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c87Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*73+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx77 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c87Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*74+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx78 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c88Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*75+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx79 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c88Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*76+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx80 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }} if(lMax > 8) { LNsNs=9*NsNs; LNs=9*Ns;
    for(int k = 0; k < Ns; k++){
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      sumMe = 0;/*c90*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*77+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx81 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c91Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*78+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx82 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c91Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*79+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx83 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c92Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*80+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx84 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c92Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*81+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx85 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c93Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*82+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx86 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c93Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*83+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx87 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c94Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*84+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx88 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c94Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*85+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx89 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c95Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*86+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx90 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c95Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*87+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx91 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c96Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*88+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx92 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c96Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*89+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx93 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c97Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*90+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx94 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c97Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*91+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx95 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c98Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*92+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx96 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c98Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*93+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx97 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c99Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*94+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx98 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
      sumMe = 0;/*c99Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*95+i]);}
      for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx99 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    }}
}
//=======================================================================
double* getP(double* C, int Ns, int Ts, int Hs, int lMax){

  double* P = (double*) malloc(Hs*Ts*Ns*Ns*(lMax+1)*sizeof(double));
  for(int i = 0; i < Hs*Ns*Ns*(lMax+1); i++){P[i] = 0.0;}

  double cs0  = pow(PIHalf,2);
  double cs1  = pow(2.7206990464,2);
  double cs2  = pow(1.9238247452,2); double cs3  = pow(1.7562036828,2); double cs4  = pow(4.3018029072,2);
  double cs5  = pow(2.1509014536,2); double cs6  = pow(2.0779682205,2); double cs7  = pow(1.7995732672,2);
  double cs8  = pow(5.6907503408,2); double cs9  = pow(2.3232390981,2); double cs10 = pow(0.5890486225,2);
  double cs11 = pow(2.6343055241,2); double cs12 = pow(1.8627352998,2); double cs13 = pow(6.9697172942,2);
  double cs14 = pow(2.4641671809,2); double cs15 = pow(0.6512177548,2); double cs16 = pow(1.7834332706,2);
  double cs17 = pow(9.4370418280,2); double cs18 = pow(1.9263280966,2); double cs19 = pow(8.1727179596,2);
  double cs20 = pow(2.5844403427,2); double cs21 = pow(0.3539741687,2); double cs22 = pow(2.2940148014,2);
  double cs23 = pow(1.8135779397,2); double cs24 = pow(3.6271558793,2); double cs25 = pow(1.9866750947,2);
  double cs26 = pow(9.3183321738,2); double cs27 = pow(2.6899707945,2); double cs28 = pow(0.3802292509,2);
  double cs29 = pow(0.3556718963,2); double cs30 = pow(0.8712146618,2); double cs31 = pow(0.6160417952,2);
  double cs32 = pow(4.0863589798,2); double cs33 = pow(2.0431794899,2); double cs34 = pow(10.418212089,2);
  double cs35 = pow(2.7843843014,2); double cs36 = pow(0.0505981185,2); double cs37 = pow(0.4293392727,2);
  double cs38 = pow(1.7960550366,2); double cs39 = pow(4.8637400313,2); double cs40 = pow(1.8837184141,2);
  double cs41 = pow(13.583686661,2); double cs42 = pow(2.0960083567,2); double cs43 = pow(11.480310577,2);
  double cs44 = pow(2.8700776442,2); double cs45 = pow(0.0534917379,2); double cs46 = pow(0.2537335916,2);
  double cs47 = pow(2.3802320735,2); double cs48 = pow(1.8179322747,2); double cs49 = pow(16.055543121,2);
  double cs50 = pow(1.9190044477,2); double cs51 = pow(4.9548481782,2); double cs52 = pow(2.1455121971,2);
  double cs53 = pow(12.510378411,2); double cs54 = pow(2.9487244699,2);

   // SUM M's UP!
  for(int i = 0; i < Hs; i++){
    for(int j = 0; j < Ts; j++){
      for(int k = 0; k < Ns; k++){
        for(int kd = 0; kd < Ns; kd++){
          P[(lMax + 1)*Ns*Ns*Ts*i+ (lMax+1)*Ns*Ns*j+ 0 +Ns*k+kd] =
            cs0*C[100*Ns*Ts*i + 100*Ns*j + 0 + k]*C[100*Ns*Ts*i + 100*Ns*j + 0 + kd];
        }
      }
    }
  }  if(lMax > 0){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ Ns*Ns + Ns*k+kd] =
                cs1*C[100*Ns*Ts*i + 100*Ns*j + 1*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 1*Ns + kd]
             +2*cs2*C[100*Ns*Ts*i + 100*Ns*j + 2*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 2*Ns + kd]
             +2*cs2*C[100*Ns*Ts*i + 100*Ns*j + 3*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 3*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 1){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 2*Ns*Ns + Ns*k+kd] =
                cs3*C[100*Ns*Ts*i + 100*Ns*j + 4*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 4*Ns + kd]
             +2*cs4*C[100*Ns*Ts*i + 100*Ns*j + 5*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 5*Ns + kd]
             +2*cs4*C[100*Ns*Ts*i + 100*Ns*j + 6*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 6*Ns + kd]
             +2*cs5*C[100*Ns*Ts*i + 100*Ns*j + 7*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 7*Ns + kd]
             +2*cs5*C[100*Ns*Ts*i + 100*Ns*j + 8*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 8*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 2){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 3*Ns*Ns + Ns*k+kd] =
                cs6*C[100*Ns*Ts*i + 100*Ns*j + 9*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 9*Ns + kd]
             +2*cs7*C[100*Ns*Ts*i + 100*Ns*j + 10*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 10*Ns + kd]
             +2*cs7*C[100*Ns*Ts*i + 100*Ns*j + 11*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 11*Ns + kd]
             +2*cs8*C[100*Ns*Ts*i + 100*Ns*j + 12*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 12*Ns + kd]
             +2*cs8*C[100*Ns*Ts*i + 100*Ns*j + 13*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 13*Ns + kd]
             +2*cs9*C[100*Ns*Ts*i + 100*Ns*j + 14*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 14*Ns + kd]
             +2*cs9*C[100*Ns*Ts*i + 100*Ns*j + 15*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 15*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 3){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 4*Ns*Ns + Ns*k+kd] =
                cs10*C[100*Ns*Ts*i + 100*Ns*j + 16*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 16*Ns + kd]
             +2*cs11*C[100*Ns*Ts*i + 100*Ns*j + 17*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 17*Ns + kd]
             +2*cs11*C[100*Ns*Ts*i + 100*Ns*j + 18*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 18*Ns + kd]
             +2*cs12*C[100*Ns*Ts*i + 100*Ns*j + 19*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 19*Ns + kd]
             +2*cs12*C[100*Ns*Ts*i + 100*Ns*j + 20*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 20*Ns + kd]
             +2*cs13*C[100*Ns*Ts*i + 100*Ns*j + 21*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 21*Ns + kd]
             +2*cs13*C[100*Ns*Ts*i + 100*Ns*j + 22*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 22*Ns + kd]
             +2*cs14*C[100*Ns*Ts*i + 100*Ns*j + 23*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 23*Ns + kd]
             +2*cs14*C[100*Ns*Ts*i + 100*Ns*j + 24*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 24*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 4){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 5*Ns*Ns + Ns*k+kd] =
                cs15*C[100*Ns*Ts*i + 100*Ns*j + 25*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 25*Ns + kd]
             +2*cs16*C[100*Ns*Ts*i + 100*Ns*j + 26*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 26*Ns + kd]
             +2*cs16*C[100*Ns*Ts*i + 100*Ns*j + 27*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 27*Ns + kd]
             +2*cs17*C[100*Ns*Ts*i + 100*Ns*j + 28*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 28*Ns + kd]
             +2*cs17*C[100*Ns*Ts*i + 100*Ns*j + 29*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 29*Ns + kd]
             +2*cs18*C[100*Ns*Ts*i + 100*Ns*j + 30*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 30*Ns + kd]
             +2*cs18*C[100*Ns*Ts*i + 100*Ns*j + 31*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 31*Ns + kd]
             +2*cs19*C[100*Ns*Ts*i + 100*Ns*j + 32*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 32*Ns + kd]
             +2*cs19*C[100*Ns*Ts*i + 100*Ns*j + 33*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 33*Ns + kd]
             +2*cs20*C[100*Ns*Ts*i + 100*Ns*j + 34*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 34*Ns + kd]
             +2*cs20*C[100*Ns*Ts*i + 100*Ns*j + 35*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 35*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 5){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 6*Ns*Ns + Ns*k+kd] =
                cs21*C[100*Ns*Ts*i + 100*Ns*j + 36*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 36*Ns + kd]
             +2*cs22*C[100*Ns*Ts*i + 100*Ns*j + 37*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 37*Ns + kd]
             +2*cs22*C[100*Ns*Ts*i + 100*Ns*j + 38*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 38*Ns + kd]
             +2*cs23*C[100*Ns*Ts*i + 100*Ns*j + 39*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 39*Ns + kd]
             +2*cs23*C[100*Ns*Ts*i + 100*Ns*j + 40*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 40*Ns + kd]
             +2*cs24*C[100*Ns*Ts*i + 100*Ns*j + 41*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 41*Ns + kd]
             +2*cs24*C[100*Ns*Ts*i + 100*Ns*j + 42*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 42*Ns + kd]
             +2*cs25*C[100*Ns*Ts*i + 100*Ns*j + 43*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 43*Ns + kd]
             +2*cs25*C[100*Ns*Ts*i + 100*Ns*j + 44*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 44*Ns + kd]
             +2*cs26*C[100*Ns*Ts*i + 100*Ns*j + 45*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 45*Ns + kd]
             +2*cs26*C[100*Ns*Ts*i + 100*Ns*j + 46*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 46*Ns + kd]
             +2*cs27*C[100*Ns*Ts*i + 100*Ns*j + 47*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 47*Ns + kd]
             +2*cs27*C[100*Ns*Ts*i + 100*Ns*j + 48*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 48*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 6){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 7*Ns*Ns + Ns*k+kd] =
                cs28*C[100*Ns*Ts*i + 100*Ns*j + 49*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 49*Ns + kd]
             +2*cs29*C[100*Ns*Ts*i + 100*Ns*j + 50*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 50*Ns + kd]
             +2*cs29*C[100*Ns*Ts*i + 100*Ns*j + 51*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 51*Ns + kd]
             +2*cs30*C[100*Ns*Ts*i + 100*Ns*j + 52*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 52*Ns + kd]
             +2*cs30*C[100*Ns*Ts*i + 100*Ns*j + 53*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 53*Ns + kd]
             +2*cs31*C[100*Ns*Ts*i + 100*Ns*j + 54*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 54*Ns + kd]
             +2*cs31*C[100*Ns*Ts*i + 100*Ns*j + 55*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 55*Ns + kd]
             +2*cs32*C[100*Ns*Ts*i + 100*Ns*j + 56*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 56*Ns + kd]
             +2*cs32*C[100*Ns*Ts*i + 100*Ns*j + 57*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 57*Ns + kd]
             +2*cs33*C[100*Ns*Ts*i + 100*Ns*j + 58*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 58*Ns + kd]
             +2*cs33*C[100*Ns*Ts*i + 100*Ns*j + 59*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 59*Ns + kd]
             +2*cs34*C[100*Ns*Ts*i + 100*Ns*j + 60*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 60*Ns + kd]
             +2*cs34*C[100*Ns*Ts*i + 100*Ns*j + 61*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 61*Ns + kd]
             +2*cs35*C[100*Ns*Ts*i + 100*Ns*j + 62*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 62*Ns + kd]
             +2*cs35*C[100*Ns*Ts*i + 100*Ns*j + 63*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 63*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 7){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 8*Ns*Ns + Ns*k+kd] =
                cs36*C[100*Ns*Ts*i + 100*Ns*j + 64*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 64*Ns + kd]
             +2*cs37*C[100*Ns*Ts*i + 100*Ns*j + 65*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 65*Ns + kd]
             +2*cs37*C[100*Ns*Ts*i + 100*Ns*j + 66*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 66*Ns + kd]
             +2*cs38*C[100*Ns*Ts*i + 100*Ns*j + 67*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 67*Ns + kd]
             +2*cs38*C[100*Ns*Ts*i + 100*Ns*j + 68*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 68*Ns + kd]
             +2*cs39*C[100*Ns*Ts*i + 100*Ns*j + 69*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 69*Ns + kd]
             +2*cs39*C[100*Ns*Ts*i + 100*Ns*j + 70*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 70*Ns + kd]
             +2*cs40*C[100*Ns*Ts*i + 100*Ns*j + 71*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 71*Ns + kd]
             +2*cs40*C[100*Ns*Ts*i + 100*Ns*j + 72*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 72*Ns + kd]
             +2*cs41*C[100*Ns*Ts*i + 100*Ns*j + 73*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 73*Ns + kd]
             +2*cs41*C[100*Ns*Ts*i + 100*Ns*j + 74*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 74*Ns + kd]
             +2*cs42*C[100*Ns*Ts*i + 100*Ns*j + 75*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 75*Ns + kd]
             +2*cs42*C[100*Ns*Ts*i + 100*Ns*j + 76*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 76*Ns + kd]
             +2*cs43*C[100*Ns*Ts*i + 100*Ns*j + 77*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 77*Ns + kd]
             +2*cs43*C[100*Ns*Ts*i + 100*Ns*j + 78*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 78*Ns + kd]
             +2*cs44*C[100*Ns*Ts*i + 100*Ns*j + 79*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 79*Ns + kd]
             +2*cs44*C[100*Ns*Ts*i + 100*Ns*j + 80*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 80*Ns + kd];
          }
        }
      }
    }
  }  if(lMax > 8){
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
            P[(lMax + 1)*Ns*Ns*Ts*i+(lMax+1)*Ns*Ns*j+ 9*Ns*Ns + Ns*k+kd] =
                cs45*C[100*Ns*Ts*i + 100*Ns*j + 81*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 81*Ns + kd]
             +2*cs46*C[100*Ns*Ts*i + 100*Ns*j + 82*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 82*Ns + kd]
             +2*cs46*C[100*Ns*Ts*i + 100*Ns*j + 83*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 83*Ns + kd]
             +2*cs47*C[100*Ns*Ts*i + 100*Ns*j + 84*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 84*Ns + kd]
             +2*cs47*C[100*Ns*Ts*i + 100*Ns*j + 85*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 85*Ns + kd]
             +2*cs48*C[100*Ns*Ts*i + 100*Ns*j + 86*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 86*Ns + kd]
             +2*cs48*C[100*Ns*Ts*i + 100*Ns*j + 87*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 87*Ns + kd]
             +2*cs49*C[100*Ns*Ts*i + 100*Ns*j + 88*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 88*Ns + kd]
             +2*cs49*C[100*Ns*Ts*i + 100*Ns*j + 89*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 89*Ns + kd]
             +2*cs50*C[100*Ns*Ts*i + 100*Ns*j + 90*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 90*Ns + kd]
             +2*cs50*C[100*Ns*Ts*i + 100*Ns*j + 91*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 91*Ns + kd]
             +2*cs51*C[100*Ns*Ts*i + 100*Ns*j + 92*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 92*Ns + kd]
             +2*cs51*C[100*Ns*Ts*i + 100*Ns*j + 93*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 93*Ns + kd]
             +2*cs52*C[100*Ns*Ts*i + 100*Ns*j + 94*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 94*Ns + kd]
             +2*cs52*C[100*Ns*Ts*i + 100*Ns*j + 95*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 95*Ns + kd]
             +2*cs53*C[100*Ns*Ts*i + 100*Ns*j + 96*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 96*Ns + kd]
             +2*cs53*C[100*Ns*Ts*i + 100*Ns*j + 97*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 97*Ns + kd]
             +2*cs54*C[100*Ns*Ts*i + 100*Ns*j + 98*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 98*Ns + kd]
             +2*cs54*C[100*Ns*Ts*i + 100*Ns*j + 99*Ns + k]*C[100*Ns*Ts*i + 100*Ns*j + 99*Ns + kd];
          }
        }
      }
    }
   }
  return P;
}
//=======================================================================
int main(int argc, char* argv[]){
  int lMax = 9; int Ns = 5; int Hsize = 9999;
  int*  totalAN = (int*) malloc(sizeof(int));
  int*  Ntypes = (int*) malloc(sizeof(int));
  double* alphas = getAlphas(Ns); double* betas = getBetas(Ns);
  double NsNs = Ns*Ns;
  char* ch =  argv[1];
  int* typeNs = getInfo(totalAN, Ntypes);
  int* types = (int*) malloc(sizeof(int)*Ntypes[0]);
  double* Apos = getApos(totalAN, Ntypes, typeNs, types);
  double* Hpos = getHpos(Hsize,ch);
  double* x  = (double*) malloc(sizeof(double)*totalAN[0]);
  double* y  = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z  = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z2 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z4 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z6 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z8 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r2 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r4 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r6 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r8 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* ReIm2 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm3 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm4 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm5 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm6 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm7 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm8 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* ReIm9 = (double*) malloc(2*sizeof(double)*totalAN[0]);// 2 -> Re + ixIm
  double* exes = (double*) malloc (sizeof(double)*totalAN[0]);
  double* preCoef = (double*) malloc(96*sizeof(double)*totalAN[0]);
  double* bOa = (double*) malloc((lMax+1)*NsNs*sizeof(double));
  double* aOa = (double*) malloc((lMax+1)*Ns*sizeof(double));
  int Asize = 0;
  double* P;
  
  double* cn = (double*) malloc(100*Ntypes[0]*Ns*Hsize*sizeof(double));
  for(int i = 0; i < 100*Ntypes[0]*Ns*Hsize; i++){cn[i] = 0.0;}
  
  //MAKESURE TO NULLIFY THE CNs!!!!!!!
  //Triple Check the implementation, Triple times. Then Triple that again.
  getAlphaBeta(aOa,bOa,alphas,betas,Ns);
  for(int i = 0; i < 9999; i++){
    for(int j = 0; j < Ntypes[0]; j++){
      Asize = getFilteredPos(x, y, z, Apos, Hpos,typeNs, 100.0, i, j);
      getRsZs(x, y, z, r2,r4,r6,r8,z2,z4,z6,z8, Asize);
      getCfactors(preCoef,Asize,x,y,z,z2,z4,z6,z8,r2,r4,r6,r8,ReIm2,ReIm3,ReIm4,ReIm5,ReIm6,ReIm7,ReIm8,ReIm9, totalAN[0]);
      getC(cn,preCoef,x,y,z,r2,bOa,aOa,exes,Hsize,totalAN[0],Asize,Ns,Ntypes[0], lMax, i, j);
    }
  }


  P = getP(cn, Ns, Ntypes[0], Hsize, lMax);


  for(int i = 0; i < Hsize; i++){
    for(int j = 0; j < Ntypes[0]; j++){
      for(int l = 0; l < lMax + 1; l++){
        for(int k = 0; k < Ns; k++){
          for(int kd = 0; kd < Ns; kd++){
             std::cout << P[i*Ntypes[0]*Ns*Ns*(lMax+1) + j*Ns*Ns*(lMax+1) + l*Ns*Ns + k*Ns + kd] << " ";
          }
        }
      }
    }
    std::cout << std::endl;
  }
}











