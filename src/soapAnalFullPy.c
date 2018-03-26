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
//-----------------------------------------------------------
double* getReIm2(double* x, double* y, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*x[i] - y[i]*y[i];
    c3[2*i+1] = 2*y[i]*x[i]; 
  }
}
//-----------------------------------------------------------
double* getReIm3(double* x, double* y, double* c2, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*c2[2*i] - y[i]*c2[2*i + 1];
    c3[2*i+1] = x*c2[2*i+1] + y*c2[2*i  ]; 
  }
}
//-----------------------------------------------------------
double* mulReIm(double* c1, double* c2, double* c3, int Asize){
  
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = c1[2*i  ]*c2[2*i  ] - c1[2*i+1]*c2[2*i+1];
    c3[2*i+1] = c1[2*i  ]*c2[2*i+1] + c1[2*i+1]*c2[2*i  ]; 
  }
}
//-----------------------------------------------------------
void printC(double* C, int Nsize,int Ntypes){
  for(int i = 0; i < 9999 ; i++){
    for(int j = 0; j < Ntypes ; j++){
      for(int k = 0; k < Nsize ; k++){
        std::cout << C[i*Nsize*Ntypes + j*Ntypes + k ] << " ";
      }
    }
    std::cout << std::endl;
  }
}
//-----------------------------------------------------------
double sumMe(double* C, int Nsize,int Ntypes){
  double mySum = 0;
  for(int i = 0; i < 9999 ; i++){
    for(int j = 0; j < Ntypes ; j++){
      for(int k = 0; k < Nsize ; k++){
        mySum += C[i*Nsize*Ntypes + j*Ntypes + k ];
      }
    }
  }
  return mySum;
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
//-----------------------------------------------------------
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
//-----------------------------------------------------------
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
//-----------------------------------------------------------
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
//-----------------------------------------------------------
//========================================
double* getHpos(int Hsize){

  FILE* pFile;
  pFile = fopen("H.dat","r");
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
//========================================
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
double* getRsZs(double* x, double* y, double* z,double* r2,double* r4,double* r6,double* r8,double* z2,double* z4,double* z6,double* z8,, int size){
  for(int i = 0; i < size; i++){
       r2[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
       r4[i] = r2[i]*r2[i]; r6[i] = r2[i]*r4[i]; r8[i] = r4[i]*r4[i];
       z2[i] = z[i]*z[i]; z4[i] = z2[i]*z2[i]; z6[i] = z2[i]*z4[i]; z8[i] = z4[i]*z4[i];
  }

}
//================================================================
double* get_AlphaBeta(double* alphas, double* betas, int Nsize){
  int l = 9;

  int  NsNs = Nsize*Nsize;
  double  oneO1alpha;      double  oneO1alpha2; double  oneO1alpha3;
  double  oneO1alpha4; double  oneO1alpha5; double  oneO1alpha6;
  double  oneO1alpha7; double  oneO1alpha8; double  oneO1alpha9;
  double  oneO1alphaSqrt;// = (double*) malloc(Nsize*sizeof(double));
  double* alphaO1alpha = (double*) malloc((l+1)*Nsize*sizeof(double));
  double* bOa = (double*) malloc((l+1)*NsNs*sizeof(double));

  // MY POEWR MISSING (see beggning functions);

    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[k]);
      oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[k] = -alphas[k]*oneO1alpha; //got alpha_0k

      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[n*Nsize + k] = betas[n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_0nk
    }

  if(l>0){
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[Nsize + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[Nsize + k] = -alphas[Nsize + k]*oneO1alpha; //got alpha_1k
      
      oneO1alpha2 = oneO1alpha*oneO1alpha;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha2;
      for(int n = 0; n < Nsize; n++){
        bOa[NsNs + n*Nsize + k] = betas[NsNs + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_1nk
    }
  } if(l>1){
    int shift1 = 2*Nsize; int shift2 = 2*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_2k

      oneO1alpha3 = oneO1alpha2*oneO1alpha;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_2nk
    }
  } if(l>2){
    int shift1 = 3*Nsize; int shift2 = 3*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_3k

      oneO1alpha4 = oneO1alpha2*oneO1alpha2;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_3nk
    }
  } if(l>3){
    int shift1 = 4*Nsize; int shift2 = 4*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_4k

      oneO1alpha5 = oneO1alpha2*oneO1alpha3;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_4nk
    }
  } if(l>4){
    int shift1 = 5*Nsize; int shift2 = 5*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_5k

      oneO1alpha6 = oneO1alpha3*oneO1alpha3;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_5nk
    }
  } if(l>5){
    int shift1 = 6*Nsize; int shift2 = 6*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_6k

      oneO1alpha7 = oneO1alpha4*oneO1alpha3;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_6nk
    }
  } if(l>6){
    int shift1 = 7*Nsize; int shift2 = 7*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_7k

      oneO1alpha8 = oneO1alpha4*oneO1alpha4;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_7nk
    }
  } if(l>7){
    int shift1 = 8*Nsize; int shift2 = 8*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_8k

      oneO1alpha9 = oneO1alpha4*oneO1alpha5;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_8nk
    }
  }if(l>8){
    int shift1 = 9*Nsize; int shift2 = 9*NsNs;
    for(int k = 0; k < Nsize; k++){
      oneO1alpha = 1.0/(1.0 + alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      alphaO1alpha[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_9k

      oneO1alpha10 = oneO1alpha5*oneO1alpha5;
      oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
      for(int n = 0; n < Nsize; n++){
        bOa[shift2 + n*Nsize + k] = betas[shift2 + n*Nsize + k]*oneO1alphaSqrtX;
      } // got beta_9nk
    }
  }
}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
       for(int i = 0; i < Asize; i++){
      /*c20  */  preCoef[           + i] = (3*z2[i] - r2[i]);

      /*c21Re*/  preCoef[totalAN    + i] = x[i]*z[i];
      /*c21Im*/  preCoef[totalAN*2  + i] = y[i]*z[i];
      /*c22Re*/  preCoef[totalAN*3  + i] = ReIm2[2*i];
      /*c22Im*/  preCoef[totalAN*4  + i] = ReIm2[2*i + 1];

      /*c30  */  preCoef[totalAN*5  + i] = (5*z2[i] - 3*r2[i]);

      /*c31Re*/  preCoef[totalAN*6  + i] = x[i]*(5*z2[i] - r2[i]);
      /*c31Im*/  preCoef[totalAN*7  + i] = y[i]*(5*z2[i] - r2[i]);
      /*c32Re*/  preCoef[totalAN*8  + i] = ReIm2[2*i]*z[i];
      /*c32Im*/  preCoef[totalAN*9  + i] = ReIm2[2*i + 1]*z[i];
      /*c33Re*/  preCoef[totalAN*10 + i] = ReIm3[2*i    ];
      /*c33Im*/  preCoef[totalAN*11 + i] = ReIm3[2*i + 1];

      /*c40  */  preCoef[totalAN*12 + i] = (35*z4[i] - 30*z2[i]*r2[i]* + 3*r4[i]);

      /*c41Re*/  preCoef[totalAN*13 + i] = z[i]*x[i]*(7*z2[i] - 3*r2[i]);
      /*c41Im*/  preCoef[totalAN*14 + i] = z[i]*y[i]*(7*z2[i] - 3*r2[i]);
      /*c42Re*/  preCoef[totalAN*15 + i] =      ReIm2[2*i]*(7*z2[i] - r2[i]);
      /*c42Im*/  preCoef[totalAN*16 + i] =      ReIm2[2*i + 1]*(7*z2[i] - r2[i]);
      /*c43Re*/  preCoef[totalAN*17 + i] = z[i]*ReIm3[2*i    ];
      /*c43Im*/  preCoef[totalAN*18 + i] = z[i]*ReIm3[2*i + 1];
      /*c44Re*/  preCoef[totalAN*19 + i] =      ReIm4[2*i    ];
      /*c44Im*/  preCoef[totalAN*20 + i] =      ReIm4[2*i + 1];

      /*c50  */  preCoef[totalAN*21 + i] = z[i]*(63*z4[i] - 70*z2[i]*r2[i]* + 15*r4[i]);

      /*c51Re*/  preCoef[totalAN*22 + i] =      x[i]*(21*z4[i] - 14*z2[i]*r2[i]* + r4[i]);
      /*c51Im*/  preCoef[totalAN*23 + i] =      y[i]*(21*z4[i] - 14*z2[i]*r2[i]* + r4[i]);
      /*c52Re*/  preCoef[totalAN*24 + i] = z[i]*ReIm2[2*i    ]*(3*z2[i] - r2[i]);
      /*c52Im*/  preCoef[totalAN*25 + i] = z[i]*ReIm2[2*i + 1]*(3*z2[i] - r2[i]);
      /*c53Re*/  preCoef[totalAN*26 + i] =      ReIm3[2*i    ]*(9*z2[i] - r2[i]);
      /*c53Im*/  preCoef[totalAN*27 + i] =      ReIm3[2*i + 1]*(9*z2[i] - r2[i]);
      /*c54Re*/  preCoef[totalAN*28 + i] = z[i]*ReIm4[2*i    ];
      /*c54Im*/  preCoef[totalAN*29 + i] = z[i]*ReIm4[2*i + 1];
      /*c55Re*/  preCoef[totalAN*30 + i] =      ReIm5[2*i    ];
      /*c55Im*/  preCoef[totalAN*31 + i] =      ReIm5[2*i + 1];

      /*c60  */  preCoef[totalAN*32 + i] = (231*z4[i] - 315*z2[i]*r2[i]* + 105*r4[i] - 5*r6[i]);

      /*c61Re*/  preCoef[totalAN*33 + i] = z[i]*x[i]*();
      /*c61Im*/  preCoef[totalAN*34 + i] = z[i]*y[i]*();
      /*c62Re*/  preCoef[totalAN*35 + i] =      ReIm2[2*i    ]*();
      /*c62Im*/  preCoef[totalAN*36 + i] =      ReIm2[2*i + 1]*();
      /*c63Re*/  preCoef[totalAN*37 + i] = z[i]*ReIm3[2*i    ]*();
      /*c63Im*/  preCoef[totalAN*38 + i] = z[i]*ReIm3[2*i + 1]*();
      /*c64Re*/  preCoef[totalAN*39 + i] =      ReIm4[2*i    ]*();
      /*c64Im*/  preCoef[totalAN*40 + i] =      ReIm4[2*i + 1]*();
      /*c65Re*/  preCoef[totalAN*41 + i] = z[i]*ReIm5[2*i    ]*();
      /*c65Im*/  preCoef[totalAN*42 + i] = z[i]*ReIm5[2*i + 1]*();
      /*c66Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c66Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();

      /*c70  */  preCoef[totalAN*32 + i] = (231*z4[i] - 315*z2[i]*r2[i]* + 105*r4[i] - 5*r6[i]);

      /*c71Re*/  preCoef[totalAN*33 + i] = z[i]*x[i]*();
      /*c71Im*/  preCoef[totalAN*34 + i] = z[i]*y[i]*();
      /*c72Re*/  preCoef[totalAN*35 + i] =      ReIm2[2*i    ]*();
      /*c72Im*/  preCoef[totalAN*36 + i] =      ReIm2[2*i + 1]*();
      /*c73Re*/  preCoef[totalAN*37 + i] = z[i]*ReIm3[2*i    ]*();
      /*c73Im*/  preCoef[totalAN*38 + i] = z[i]*ReIm3[2*i + 1]*();
      /*c74Re*/  preCoef[totalAN*39 + i] =      ReIm4[2*i    ]*();
      /*c74Im*/  preCoef[totalAN*40 + i] =      ReIm4[2*i + 1]*();
      /*c75Re*/  preCoef[totalAN*41 + i] = z[i]*ReIm5[2*i    ]*();
      /*c75Im*/  preCoef[totalAN*42 + i] = z[i]*ReIm5[2*i + 1]*();
      /*c76Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c76Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c77Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c77Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();

      /*c80  */  preCoef[totalAN*32 + i] = (231*z4[i] - 315*z2[i]*r2[i]* + 105*r4[i] - 5*r6[i]);

      /*c81Re*/  preCoef[totalAN*33 + i] = z[i]*x[i]*();
      /*c81Im*/  preCoef[totalAN*34 + i] = z[i]*y[i]*();
      /*c82Re*/  preCoef[totalAN*35 + i] =      ReIm2[2*i    ]*();
      /*c82Im*/  preCoef[totalAN*36 + i] =      ReIm2[2*i + 1]*();
      /*c83Re*/  preCoef[totalAN*37 + i] = z[i]*ReIm3[2*i    ]*();
      /*c83Im*/  preCoef[totalAN*38 + i] = z[i]*ReIm3[2*i + 1]*();
      /*c84Re*/  preCoef[totalAN*39 + i] =      ReIm4[2*i    ]*();
      /*c84Im*/  preCoef[totalAN*40 + i] =      ReIm4[2*i + 1]*();
      /*c85Re*/  preCoef[totalAN*41 + i] = z[i]*ReIm5[2*i    ]*();
      /*c85Im*/  preCoef[totalAN*42 + i] = z[i]*ReIm5[2*i + 1]*();
      /*c86Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c86Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c87Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c87Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c88Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c88Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();

      /*c90  */  preCoef[totalAN*32 + i] = (231*z4[i] - 315*z2[i]*r2[i]* + 105*r4[i] - 5*r6[i]);

      /*c91Re*/  preCoef[totalAN*33 + i] = z[i]*x[i]*();
      /*c91Im*/  preCoef[totalAN*34 + i] = z[i]*y[i]*();
      /*c92Re*/  preCoef[totalAN*35 + i] =      ReIm2[2*i    ]*();
      /*c92Im*/  preCoef[totalAN*36 + i] =      ReIm2[2*i + 1]*();
      /*c93Re*/  preCoef[totalAN*37 + i] = z[i]*ReIm3[2*i    ]*();
      /*c93Im*/  preCoef[totalAN*38 + i] = z[i]*ReIm3[2*i + 1]*();
      /*c94Re*/  preCoef[totalAN*39 + i] =      ReIm4[2*i    ]*();
      /*c94Im*/  preCoef[totalAN*40 + i] =      ReIm4[2*i + 1]*();
      /*c95Re*/  preCoef[totalAN*41 + i] = z[i]*ReIm5[2*i    ]*();
      /*c95Im*/  preCoef[totalAN*42 + i] = z[i]*ReIm5[2*i + 1]*();
      /*c96Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c96Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c97Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c97Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c98Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c98Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      /*c99Re*/  preCoef[totalAN*43 + i] =      ReIm6[2*i    ]*();
      /*c99Im*/  preCoef[totalAN*44 + i] =      ReIm6[2*i + 1]*();
      }
//================================================================
int getC(double* C, double* x, double* y, double* z,double* r2, double* ReIm2, double* ReIm3, double* ReIm4, double* ReIm5, double* ReIm6,double* ReIm7, double* ReIm8, double* ReIm9, double* bOa, double* aOa, double* exes, int Hsize, int Asize, int Nsize, int Ntypes,int posI, int typeJ){

  if(Asize == 0){
    return 0;
  }

  double sumMe = 0;
  int NsNs = Nsize*Nsize;
  int NsJ = 55*Nsize*typeJ;
  int NsNtI = Nsize*Ntypes*posI;
  int LNsNs ;
  int Nx2 = 2*Nsize;
  int Nx3 = 3*Nsize;
  int LNs ;
  int LNsx2 = 2*LNs;
  int NsTsI = 2*55*Nsize*Ntypes*posI;

  for(int k = 0; k < Nsize; k++){
    sumMe = 0;
    for(int i = 0; i < Asize; i++){ sumMe += exp(aOa[k]*r2[i]);}
    for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + n] += bOa[n*Nsize + k]*sumMe; }
  }

  if(l > 0){
     LNsNs=NsNs;
     LNs=Nsize;
    for(int k = 0; k < Nsize; k++){
      //exponents
      for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}

      sumMe = 0;//c10
      for(int i = 0; i < Asize; i++){sumMe += exes[i];}
      for(int n = 0; n < Nsize; n++){C[NsTsI + NsJ + Nsize + n] += bOa[LNsNs + n*Nsize + k]*sumMe;} 

      sumMe = 0;//c11Re
      for(int i = 0; i < Asize; i++){sumMe += exes[i]*x[i];}
      for(int n = 0; n < Nsize; n++){C[NsTsI + NsJ + Nx2 + n] += bOa[LNsNs + n*Nsize + k]*sumMe;}

      sumMe = 0;//c11Im
      for(int i = 0; i < Asize; i++){sumMe += exes[i]*y[i];}
      for(int n = 0; n < Nsize; n++){C[NsTsI + NsJ + Nx3 + n] += bOa[LNsNs + n*Nsize + k]*sumMe;}

    }
  }

  if(l > 1){
    LNsNs=2*NsNs;
    LNs=2*Nsize;
    getReIm2(x, y, ReIm2,Asize);
    for(int k = 0; k < Nsize; k++){
      //exponents
      for(int i = 0; i < Asize; i++){ exes[i] = exp(aOa[LNs + k]*r2[i]);}
      sumMe = 0;//c20
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*(3*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx4 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c21Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*x[i]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx5 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c21Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*y[i]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx6 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c22Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx7 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c22Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i + 1];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx8 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
    }
  }
  if(l > 2){
    LNsNs=3*NsNs;
    LNs=3*Nsize;
    getReIm3(x, y, ReIm2, ReIm3, Asize);
    for(int k = 0; k < Nsize; k++){
      //exponents
      for(int i = 0; i < Asize; i++){ exes[i] = exp(aOa[LNs + k]*r2[i]);}
      sumMe = 0;//c30
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*(5*z2[i] - 3*r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx9 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c31Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*x[i]*(5*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx10 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c31Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*y[i]*(5*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx11 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c32Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx12 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c32Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i + 1]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx13 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c33Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i    ];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx14 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c33Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i + 1];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx15 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
    }
  }
  if(l > 3){
    LNsNs=4*NsNs;
    LNs=4*Nsize;
    getMulReIm(ReIm2,ReIm2, ReIm4, Asize);
    for(int k = 0; k < Nsize; k++){
      //exponents
      for(int i = 0; i < Asize; i++){ exes[i] = exp(aOa[LNs + k]*r2[i]);}
      sumMe = 0;//c40
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*(35*z4[i] - 30*z2[i]*r2[i] + 3*r4[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx16 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c41Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*x[i]*z[i]*(7*z2[i] - 3*r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx17 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c41Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*y[i]*z[i]*(7*z2[i] - 3*r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx18 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c42Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i]*(7*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx19 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c42Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i + 1]*(7*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx20 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c43Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i    ]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx21 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c43Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i + 1]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx22 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c44Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i    ];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx23 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c44Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i + 1];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx24 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
    }
  }
  if(l > 4){
    LNsNs=5*NsNs;
    LNs=5*Nsize;
    getMulReIm(ReIm2,ReIm3, ReIm5, Asize);
    for(int k = 0; k < Nsize; k++){
      //exponents
      for(int i = 0; i < Asize; i++){ exes[i] = exp(aOa[LNs + k]*r2[i]);}
      sumMe = 0;//c50
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*(35*z4[i] - 30*z2[i]*r2[i] + 3*r4[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx16 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c51Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*x[i]*z[i]*(7*z2[i] - 3*r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx17 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c51Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*y[i]*z[i]*(7*z2[i] - 3*r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx18 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c52Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i]*(7*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx19 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c52Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm2[2*i + 1]*(7*z2[i] - r2[i]);}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx20 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c53Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i    ]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx21 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c53Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm3[2*i + 1]*z[i];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx22 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c54Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i    ];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx23 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c54Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i + 1];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx24 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c55Re
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i    ];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx23 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
      sumMe = 0;//c55Im
      for(int i = 0; i < Asize; i++){ sumMe += exes[i]*ReIm4[2*i + 1];}
      for(int n = 0; n < Nsize; n++){ C[NsTsI + NsJ + Nx24 + n] += bOa[LNsNs + n*Nsize + k]*sumMe; }
    }
  }
}
//=======================================================================
int main(){
  int l = 9;
  int Nsize = 5;
  int Hsize = 9999;
  int*  totalAN = (int*) malloc(sizeof(int));
  int*  Ntypes = (int*) malloc(sizeof(int));
  int* types;
  double* alphas = getAlphas(Nsize); double* betas = getBetas(Nsize);
  int* typeNs = getInfo(totalAN, Ntypes);
  types = (int*) malloc(sizeof(int)*Ntypes[0]);
  double* Apos = getApos(totalAN, Ntypes, typeNs, types);
  double* Hpos = getHpos(Hsize);
  double* x = (double*) malloc(sizeof(double)*totalAN[0]);
  double* y = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z2 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z4 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z6 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* z8 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r2 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r4 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r6 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* r8 = (double*) malloc(sizeof(double)*totalAN[0]);
  double* ReIm2 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm3 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm4 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm5 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm6 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm7 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm8 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* ReIm9 = (double*) malloc(2*sizeof(double)*totalAN[0])// 2 -> Re + ixIm
  double* exes = (double*) malloc (sizeof(double)*totalAN[0]);
  int Asize = 0;

  double* cn = (double*) malloc(2*55*sizeof(double)*Hsize);

  double cs0  = PIHalf;
  double cs1  = 2.7206990464;
  double cs2  = 1.9238247452; double cs3  = 1.7562036828; double cs4  = 4.3018029072;
  double cs5  = 2.1509014536; double cs6  = 2.0779682205; double cs7  = 1.7995732672;
  double cs8  = 5.6907503408; double cs9  = 2.3232390981; double cs10 = 0.5890486225;
  double cs11 = 2.6343055241; double cs12 = 1.8627352998; double cs13 = 6.9697172942;
  double cs14 = 2.4641671809; double cs15 = 0.6512177548; double cs16 = 1.7834332706;
  double cs17 = 9.4370418280; double cs18 = 1.9263280966; double cs19 = 8.1727179596;
  double cs20 = 2.5844403427; double cs21 = 0.3539741687; double cs22 = 2.2940148014;
  double cs23 = 1.8135779397; double cs24 = 3.6271558793; double cs25 = 1.9866750947;
  double cs26 = 9.3183321738; double cs27 = 2.6899707945; double cs28 = 0.3802292509;
  double cs29 = 0.3556718963; double cs30 = 0.8712146618; double cs31 = 0.6160417952;
  double cs32 = 4.0863589798; double cs33 = 2.0431794899; double cs34 = 10.418212089;
  double cs35 = 2.7843843014; double cs36 = 0.0505981185; double cs37 = 0.4293392727;
  double cs38 = 1.7960550366; double cs39 = 4.8637400313; double cs40 = 1.8837184141;
  double cs41 = 13.583686661; double cs42 = 2.0960083567; double cs43 = 11.480310577;
  double cs44 = 2.8700776442; double cs45 = 0.0534917379; double cs46 = 0.2537335916;
  double cs47 = 2.3802320735; double cs48 = 1.8179322747; double cs49 = 16.055543121;
  double cs50 = 1.9190044477; double cs51 = 4.9548481782; double cs52 = 2.1455121971;
  double cs53 = 12.510378411; double cs54 = 2.9487244699;

  //MAKESURE TO NULLIFY THE CNs!!!!!!!
  //Triple Check the implementation, Triple times. Then Triple that again.
  
  for(int i = 0; i < Hsize; i++){
    for(int j = 0; j < Ntypes[0]; j++){
      Asize = getFilteredPos(x, y, z, Apos, Hpos,typeNs, 50.0, i, j);
      getR2(x, y, z, r2,r4,r6,r8,z2,z4,z6,z8, Asize);
      getC(cn, x, y, z,r2, r4,r8,z2,z4,z6,z8, alphas, betas, Hsize, Asize, Nsize,Ntypes[0],i,j);
    }
  }

  std::cout << sumMe(cn, Nsize, Ntypes[0]) << std::endl;

}
