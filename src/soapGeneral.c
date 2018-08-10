#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

#define tot (double*) malloc(sizeof(double)*totalAN);
#define totrs (double*) malloc(sizeof(double)*totalAN*rsize);
#define sd sizeof(double)
#include "datFile.c"

//=========================================================
double factorY(int l, int m, double* c){ // OK
  return  c[(l*(l+1))/2 + m];//l+1
}
//=========================================================
double* getoOr(double* r, int rsize){
  double* oOr = (double*) malloc(sd*rsize);
  for(int w = 0; w < rsize; w++){
    oOr[w] = 1/r[w];
  }
  return oOr;
}
//=========================================================
double* getrw2(double* r, int rsize){
  double* rw2 = (double*) malloc(sd*rsize);
  for(int w = 0; w < rsize; w++){
    rw2[w] = r[w]*r[w];
  }
  return rw2;
}
//=========================================================
void expMs(double* rExpDiff, double alpha, double* r, double* ri, int isize, int rsize){
  double rDiff;
  for(int i = 0; i < isize; i++){
    for(int w = 0; w < rsize; w++){
      rDiff = r[w] - ri[i];
    if(rDiff > 5.0 ){rExpDiff[rsize*i + w] = 0.0;}
    else {rExpDiff[rsize*i + w] = exp(-alpha*rDiff*rDiff);}
    }
  }
}
//=========================================================
void expPs(double* rExpSum, double alpha, double* r, double* ri, int isize, int rsize){
  double rSum;
  for(int i = 0; i < isize; i++){
    for(int w = 0; w < rsize; w++){
    rSum = r[w] + ri[i];
    if(rSum > 5.0 ){rExpSum[rsize*i + w] = 0.0;}
    else {rExpSum[rsize*i + w] = exp(-alpha*rSum*rSum);}
    }
  }
}
//=========================================================
int getFilteredPos(double* x, double* y, double* z,double* xNow, double* yNow, double* zNow, double* ri, double* rw, double rCut, double* oOri, double* oO4arri, double* minExp, double* pluExp, double alpha, double* Apos, double* Hpos,int* typeNs, int rsize, int Ihpos, int Itype){//OK

  int shiftType = 0;
  int icount = 0;
  double ri2;
  double oOa = 1/alpha;
  double Xi; double Yi; double Zi;
  double* oO4ari = (double*) malloc(sd*typeNs[Itype]);

  for(int i = 0; i < Itype ; i++){
    shiftType += typeNs[i];
  }

  for(int i = 0; i < typeNs[Itype]; i++){
    Xi = Apos[3*shiftType + 3*i    ] - Hpos[3*Ihpos    ];
    Yi = Apos[3*shiftType + 3*i + 1] - Hpos[3*Ihpos + 1];
    Zi = Apos[3*shiftType + 3*i + 2] - Hpos[3*Ihpos + 2];
    ri2 = Xi*Xi + Yi*Yi + Zi*Zi;
    if(ri2 < rCut*rCut){
      xNow[icount] = Xi; yNow[icount] = Yi; zNow[icount] = Zi;
      ri[icount] = sqrt(ri2);
      oOri[icount] = 1/ri[icount];
      oO4ari[icount] = 0.25*oOa*oOri[icount];
      icount++;
    }
  }
  //countMax = isize ----------------------------
  double* oOr = getoOr(rw, rsize);
  for(int i = 0; i < icount; i++){
    for(int w = 0; w < rsize; w++){
      oO4arri[rsize*i + w] = oO4ari[i]*oOr[w];
    }
  }
  expMs(minExp,alpha,rw,ri,icount,rsize);
  expPs(pluExp,alpha,rw,ri,icount,rsize);

//  free(ri);
  free(oO4ari);

  return icount;
}
//=========================================================
double* getFlir(double* oO4arri,double* ri, double* minExp, double* pluExp, int icount, int rsize, int lMax){//OK
  double* Flir = (double*) malloc(sd*(lMax+1)*icount*rsize); 
//  double* rw =   getrw(100, 6);
//  int count = 0;
  //l=0
  for(int i = 0; i < icount; i++){
///    if(ri[i] < 0.01){}
    for(int w = 0; w < rsize; w++){
      Flir[rsize*i + w] = oO4arri[rsize*i + w]*(minExp[rsize*i + w] - pluExp[rsize*i + w]);
 //    exit(1);
    }
  }
  //l=1
  if(lMax>0)
  for(int i = 0; i < icount; i++){
///    if(ri[i] < 0.01){}
    for(int w = 0; w < rsize; w++){
      Flir[rsize*icount + rsize*i + w] = oO4arri[rsize*i + w]*(minExp[rsize*i + w] + pluExp[rsize*i + w] - 2*Flir[rsize*i + w]);
    }
  }
  //l>1
  if(lMax>1)
  for(int l = 2; l < lMax+1; l++){
    for(int i = 0; i < icount; i++){
 ///   if(ri[i] < 0.01){}
      for(int w = 0; w < rsize; w++){
        Flir[l*rsize*icount+rsize*i+w] = Flir[(l-2)*rsize*icount+rsize*i+w] - oO4arri[rsize*i+w]*(4*l-2)*Flir[(l-1)*rsize*icount+rsize*i+w] ;
   if(Flir[l*rsize*icount+rsize*i+w] < 0) Flir[l*rsize*icount+rsize*i+w]=0.0; // Very Important!!!

//    }
      }
    }
  }


  return Flir;

}
//=========================================================
double legendre_poly(int l, int m, double x){ // OK

  double fact,pll,pmm,pmmp1,somx2;
  int ll;

  if (m < 0 || m > l || fabs(x) > 1.0){ printf("ERROR: Bad arguments in routine legendre_poly"); exit(1);}

  pmm = 1.0;

  if(m > 0) {
    somx2=sqrt((1.0 - x)*(1.0 + x));
    fact=1.0;
    for(int i=1; i <= m; i++)
        { 
          pmm *= -fact*somx2;
          fact += 2.0;
           }
   }

  if(l == m) return pmm;

  else{
    pmmp1 = x*(2*m+1)*pmm;
    if(l==(m+1)) return pmmp1;
    else{
      for(ll=m+2; ll<=l; ll++){
        pll=(x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/ (double) (ll-m);
        pmm = pmmp1; 
        pmmp1= pll; 
      
      }
      return pll;
    }
  }

}
//=========================================================
double* getYlmi(double* x, double* y, double* z, double* oOri, double* cf, int icount, int lMax){ // OK
  double* Ylmi = (double*) malloc(2*sd*(lMax+1)*(lMax+1)*icount);
  double* legPol = (double*) malloc(sd*(lMax+1)*(lMax+1)*icount);
  double* ChiCos= (double*) malloc(sd*(lMax+1)*icount);
  double* ChiSin= (double*) malloc(sd*(lMax+1)*icount);
  double myAtan2;

  for(int i = 0; i < icount; i++){
    for(int l = 0; l < lMax + 1; l++){
      for(int m = 0; m < l+1; m++){
        legPol[icount*(lMax+1)*l + icount*m + i] = legendre_poly(l,m,z[i]*oOri[i]);
      }
    }
    
    for(int m = 0; m < lMax+1; m++){
      myAtan2 = m*atan2(y[i],x[i]);
      ChiCos[m*icount + i] = cos(myAtan2);
      ChiSin[m*icount + i] = sin(myAtan2);
//if(y<=0){ ChiSin[m*icount + i] = -sqrt(1 - ChiCos[m*icount + i]*ChiCos[m*icount + i]);}
//else{ChiSin[m*icount + i] = sqrt(1 - ChiCos[m*icount + i]*ChiCos[m*icount + i]);}
    }
  }

  for(int l = 0; l < lMax+1; l++){
    for(int m = 0; m < l+1; m++){//l+1
      for(int i = 0; i < icount; i++){

        Ylmi[2*(lMax+1)*icount*l + 2*icount*m + 2*i]
          =  factorY(l,m,cf)*legPol[icount*(lMax+1)*l + icount*m + i]*ChiCos[m*icount + i];
        Ylmi[2*(lMax+1)*icount*l + 2*icount*m + 2*i + 1]
          = factorY(l,m,cf)*legPol[icount*(lMax+1)*l + icount*m + i]*ChiSin[m*icount + i]; 

      }
    }
  }
  free(legPol); free(ChiCos); free(ChiSin);

  return Ylmi;
}
//=========================================================
double* getIntegrand(double* Flir, double* Ylmi,int rsize, int icount, int lMax){

  double* summed = (double*) malloc(2*sd*(lMax+1)*rsize*(lMax+1));
  double realY;
  double imagY;

  for(int i = 0; i < 2*(lMax+1)*rsize*(lMax+1); i++){summed[i] = 0.0;}

  for(int l = 0; l < lMax+1; l++){
    double summe = 0;
    for(int m = 0; m < l+1; m++){//l+1
     for(int i = 0; i < icount; i++){
      realY = Ylmi[2*(lMax+1)*icount*l + 2*icount*m + 2*i    ];
      imagY = Ylmi[2*(lMax+1)*icount*l + 2*icount*m + 2*i  + 1 ];
      for(int rw = 0; rw < rsize; rw++){
         summed[2*(lMax+1)*l*rsize + 2*m*rsize + 2*rw    ] += Flir[l*rsize*icount + rsize*i + rw] * realY; 
         summed[2*(lMax+1)*l*rsize + 2*m*rsize + 2*rw + 1] += Flir[l*rsize*icount + rsize*i + rw] * imagY; 
      }

    }
  }
  }
  return summed;
}
//=========================================================
void getC(double* Cs, double* ws, double* rw2, double * gns, double* summed, double rCut,int lMax, int rsize, int gnsize){ 

  for(int i = 0; i < 2*(lMax+1)*(lMax+1)*gnsize; i++){ Cs[i] = 0.0;}
  double  theSummedValue = 0;

  for(int n = 0; n < gnsize; n++){

    for(int l = 0; l < lMax+1; l++){
      for(int m = 0; m < l+1; m++){ // l+1
        for(int rw = 0; rw < rsize; rw++){

          Cs[2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m    ] += rw2[rw]*ws[rw]*gns[rsize*n + rw]*summed[2*(lMax+1)*l*rsize + 2*m*rsize + 2*rw    ]; // Re
          Cs[2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m + 1] += rw2[rw]*ws[rw]*gns[rsize*n + rw]*summed[2*(lMax+1)*l*rsize + 2*m*rsize + 2*rw + 1]; //Im


        }
     }
      }
    }
}
//=========================================================
void accumC(double* Cts, double* Cs, int lMax, int gnsize, int typeI){ 

    for(int n = 0; n < gnsize; n++){
      for(int l = 0; l < lMax+1; l++){
        for(int m = 0; m < l+1; m++){//l+1

          Cts[2*typeI*(lMax+1)*(lMax+1)*gnsize +2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m    ] = Cs[2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m    ];
          Cts[2*typeI*(lMax+1)*(lMax+1)*gnsize +2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m + 1] = Cs[2*(lMax+1)*(lMax+1)*n + l*2*(lMax+1) + 2*m + 1];

        }
      }
    }
}
//=========================================================
void getPs(double* Ps, double* Cts,  int Nt, int lMax, int gnsize){ 
  int NN = ((gnsize+1)*gnsize)/2;  int TT = ((Nt+1)*Nt)/2;
  int nshift = 0;
  for(int i = 0; i <TT*(lMax+1)*NN; i++){Ps[i] = 0.0;}
  int tshift = 0;

  for(int t1 = 0; t1 < Nt; t1++){
    for(int t2 = t1; t2 < Nt; t2++){
      for(int l = 0; l < lMax+1; l++){

      nshift = 0;
      for(int n = 0; n < gnsize; n++){
        for(int nd = n; nd < gnsize; nd++){
           for(int m = 0; m < l+1; m++){//l+1
              if(m==0){
                Ps[tshift*(lMax+1)*NN + l*NN + nshift ]
                 +=  Cts[2*t1*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*n  + l*2*(lMax+1)] // m=0
                     *Cts[2*t2*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*nd + l*2*(lMax+1)]; // m=0
              }else{

                Ps[tshift*(lMax+1)*NN + l*NN + nshift] 
                 +=  2*(Cts[2*t1*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*n  + l*2*(lMax+1) + 2*m]
                      *Cts[2*t2*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*nd + l*2*(lMax+1) + 2*m]
		   + Cts[2*t1*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*n  + l*2*(lMax+1) + 2*m + 1]
                      *Cts[2*t2*(lMax+1)*(lMax+1)*gnsize + 2*(lMax+1)*(lMax+1)*nd + l*2*(lMax+1) + 2*m + 1]);
              }
            }

            nshift++;
          }
        }
      }
              tshift++;
    }
  }

}
//=========================================================
void accumP(double* Phs, double* Ps, int Nt, int lMax, int gnsize, double rCut2, int Ihpos){ 
  int tshift=0;
  int NN = ((gnsize+1)*gnsize)/2;
  int TT = ((Nt+1)*Nt)/2;
  for(int t1 = 0; t1 < Nt; t1++){
    for(int t2 = t1; t2 < Nt; t2++){
      for(int l = 0; l < lMax+1; l++){
       int nshift=0;
        for(int n = 0; n < gnsize; n++){
          for(int nd = n; nd < gnsize; nd++){
            Phs[Ihpos*TT*(lMax+1)*NN + tshift*(lMax+1)*NN + l*NN + nshift] = 39.478417604*rCut2*Ps[tshift*(lMax+1)*NN + l*NN + nshift];// 16*9.869604401089358*Ps[tshift*(lMax+1)*NN + l*NN + nshift];
            nshift++;
          }
        }
      }
      tshift++;
    }
  }
}
//=========================================================
//=========================================================
//=========================================================
double* soap(double* c, double* Apos,double* Hpos, double* alphas,double* betas, int* typeNs, double rCut, int totalAN,int Nt,int gnsize, int lMax, int Hs, double alpha, double* rw, double* gss);
double* soap(double* c, double* Apos,double* Hpos, double* alphas,double* betas, int* typeNs, double rCut, int totalAN,int Nt,int gnsize, int lMax, int Hs, double alpha, double* rw, double* gss){

  double* cf = factorListSet();

  int rsize = 100; // constant
  double rCut2 = rCut*rCut;


  double* x    = tot  double* y    = tot  double* z    = tot double* xNow    = tot double* yNow    = tot double* zNow    = tot
  double* ris  = tot double* oOri = tot

  double* ws  = getws();
  double* oOr = getoOr(rw, rsize);  double* rw2 = getrw2(rw, rsize);

  double* oO4arri = totrs  double* minExp = totrs double* pluExp = totrs

  int Asize = 0;
  double* Cs = (double*) malloc(2*sd*(lMax+1)*(lMax+1)*gnsize);
  double* Cts = (double*) malloc(2*sd*(lMax+1)*(lMax+1)*gnsize*Nt);
  double* Ps = (double*) malloc((Nt*(Nt+1))/2*sd*(lMax+1)*((gnsize+1)*gnsize)/2);
  double* Phs = (double*) malloc(Hs*(Nt*(Nt + 1))/2*sd*(lMax+1)*((gnsize+1)*gnsize)/2);
  int icount;

  for(int Ihpos = 0; Ihpos < Hs; Ihpos++){
    for(int Itype = 0; Itype < Nt; Itype++){

        double* Ylmi; double* Flir; double* summed;

        icount = getFilteredPos(x, y, z,xNow,yNow,zNow,ris,rw,rCut, oOri, oO4arri, minExp, pluExp, alpha, Apos, Hpos,typeNs, rsize, Ihpos, Itype);

        Flir   = getFlir(oO4arri, ris, minExp, pluExp, icount, rsize, lMax);
        Ylmi   = getYlmi(xNow, yNow, zNow, oOri,cf,icount, lMax);
        summed = getIntegrand(Flir, Ylmi, rsize, icount, lMax);

        getC(Cs, ws, rw2, gss, summed, rCut,lMax, rsize, gnsize);
        accumC(Cts, Cs, lMax, gnsize, Itype);

        free(Flir); free(Ylmi); free(summed);

    }

    getPs(Ps, Cts,  Nt, lMax, gnsize); 
    accumP(Phs, Ps, Nt, lMax, gnsize,rCut2, Ihpos);
  }

//  int NN =((gnsize+1)*gnsize)/2;
//  int TT =((Nt+1)*Nt)/2;
//
//  for(int Ihpos = 0; Ihpos < Hs; Ihpos++){
//  int tshift=0;
//    for(int t1 = 0; t1 < Nt; t1++){
//      for(int t2 = t1; t2 < Nt; t2++){
//        for(int l = 0; l < lMax+1; l++){
//          int nshift=0;
//          for(int n = 0; n < gnsize; n++){
//            for(int nd = n; nd < gnsize; nd++){
//               cout <<  Phs[Ihpos*TT*(lMax+1)*NN + tshift*(lMax+1)*NN + l*NN + nshift] << " " ;
//               nshift++;
//          }
//        }
//      }
//      tshift++;
//    }
//  }
//  cout << endl;
//  }
  return Phs;
}
//=========================================================
//=========================================================
//=========================================================
