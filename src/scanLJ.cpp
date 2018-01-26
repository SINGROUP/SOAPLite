#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>

#include "../Header/fileoperation.h"
#include "../Header/scanSurface.h"

using namespace std;
using namespace arma;
//-----------------------------------------------------------------------
int main(int argc,char *argv[]){
  double r;
  int gridSize = atoi(argv[2]);
  double mypi = 3.14159265359;
  double mypi2 =2*3.14159265359;
  vec thetaScan = linspace<vec>(0,mypi,gridSize);
  vec phiScan = linspace<vec>(0,mypi2,gridSize);
  vec theta = zeros<vec>(phiScan.n_rows*thetaScan.n_rows);
  vec phi = zeros<vec>(phiScan.n_rows*thetaScan.n_rows);
  mat sphPos(100000,3);
  int n = 0;

  mat dome(theta.n_rows,3);
  vec buffvec(3)  ;
  double sparseNess =atof(argv[3]);
  double sparseNessSquared = sparseNess*sparseNess;
  double innDist = atof(argv[4]);
  double innDistSquared = innDist*innDist;
  double maxDist = atof(argv[6]);
  double maxDistSquared = maxDist*maxDist;
  double minDist = atof(argv[5]);
  double minDistSquared = minDist*minDist;

  double buffDist;
  double a;
  double b;
  double c;
  double buffDistNS;
  double d;
  double e;
  double f;

//-----------------------------------------------------------------------
  for(int i = 0; i < gridSize; i++ ){for(int j = 0; j < gridSize; j++ ){

   theta(n) = thetaScan(i); 
   phi(n) = phiScan(j); 

   n = n + 1;

  }}
//-----------------------------------------------------------------------
  mat B = getPos(argv[1]);  
  mat A(B.n_rows + 1, B.n_cols);
    for(int i = 0; i < B.n_rows; i++){
           for(int j = 0; j < B.n_cols; j++){
               A(i,j)=B(i,j); 
           }
    }
    
  vec aveAll = sum(B,0).t()/B.n_rows;
  A.row(A.n_rows - 1) = aveAll.t();


  
  string* mytype = getType(argv[1]);


  uvec q =  getSurf(A,innDist);
  uvec p = getNonSurf(A,q);
  mat surfX = A.rows(q);

  mat nonSurfY = A.rows(p);

  vec aveX = sum(surfX,0).t()/q.n_rows;
//-----------------------------------------------------------------------
for(r=0.1; r < 30 ; r = r + sparseNess){
//-----------------------------------------------------------------------
//cout << r << endl;
    for(int j = 0; j < theta.n_rows;j++){
      buffvec.row(0) = sin(theta(j))*cos(phi(j)) ;
      buffvec.row(1) =  sin(theta(j))*sin(phi(j)) ;
      buffvec.row(2) =  cos(theta(j));
      dome.row(j) = r*buffvec.t() + aveX.t();
    }
//----------------------------------------------------------------------
for(int i = 0; i < theta.n_rows - 1  ;i++){
  for(int j = i + 1; j < theta.n_rows ;j++){
      
        if(dome(i,0) <  1e3 ){
          a = dome(i,0) - dome(j,0) ;
          b = dome(i,1) - dome(j,1) ;
          c = dome(i,2) - dome(j,2) ;

          buffDist = pow(a,2) + pow(b,2) + pow(c,2);

          if(buffDist < sparseNessSquared )
            dome(j,0) = 1e7;
      }
   }}

  uvec domeEvenedUvec = find(dome.col(0) < 1e5);
  mat domeEvened = dome.rows(domeEvenedUvec);
  int evenedSize = domeEvened.n_rows;
//----------------------------------------------------------------------
    int trigger;
    int trigger2;
  for(int j = 0; j < evenedSize ;j++){
    trigger = 0;
    trigger2 = 1;
    for(int i = 0; i < surfX.n_rows  ;i++){
      
          a = surfX(i,0) - domeEvened(j,0) ;
          b = surfX(i,1) - domeEvened(j,1) ;
          c = surfX(i,2) - domeEvened(j,2) ;

          buffDist = pow(a,2) + pow(b,2) + pow(c,2);

          if(maxDistSquared  > buffDist && buffDist >minDistSquared){trigger = 1;} 
               }

     for(int k = 0; k < p.n_rows  ;k++){

        d = nonSurfY(k,0) - domeEvened(j,0) ;
        e = nonSurfY(k,1) - domeEvened(j,1) ;
        f = nonSurfY(k,2) - domeEvened(j,2) ;

        buffDistNS = pow(d,2) + pow(e,2) + pow(f,2);

        if(buffDistNS < innDistSquared ){

          trigger2 =0; }

                             
                       } 

      
             if(trigger == 0 || trigger2 == 0) domeEvened(j,0) = 1e8;
                    }

  uvec domeEvenedElimUvec = find(domeEvened.col(0) < 1e5);
  mat domeEvenedElim = domeEvened.rows(domeEvenedElimUvec);
  int evenedElimSize = domeEvenedElim.n_rows;
//------------------------------------------------------------------------

  for(int j = 0; j < evenedElimSize ;j++){

    trigger = 0;
    for(int i = 0; i < A.n_rows  ;i++){
      
          a = A(i,0) - domeEvenedElim(j,0) ;
          b = A(i,1) - domeEvenedElim(j,1) ;
          c = A(i,2) - domeEvenedElim(j,2) ;

          buffDist = pow(a,2) + pow(b,2) + pow(c,2);

          if(buffDist < minDistSquared){trigger = 1;} 
               }

      
             if(trigger == 1) domeEvenedElim(j,0) = 1e8;
                    }
  uvec domeEvenedElimElimUvec = find(domeEvenedElim.col(0) < 1e5);
  mat domeEvenedElimElim = domeEvenedElim.rows(domeEvenedElimElimUvec);
  int evenedElimSizeElim = domeEvenedElimElim.n_rows;

  vec Elj = zeros(domeEvenedElimElimUvec.n_elem);
  for(int filtI = 0; filtI < Elj.n_elem; filtI++){
    for(int filtJ = 0; filtJ < A.n_rows; filtJ++){
        a = A(filtJ,0);
        b = A(filtJ,1);
        c = A(filtJ,2);

        d = domeEvenedElimElim(filtI,0);
        e = domeEvenedElimElim(filtI,1);
        f = domeEvenedElimElim(filtI,2);

     } 
   }

///cout << "aaa" << endl;
//------------------------------------------------------------------------
  for(int j = 0; j < evenedElimSizeElim ;j++){
             cout <<  domeEvenedElimElim(j,0) << " " <<  domeEvenedElimElim(j,1) << " " <<
             domeEvenedElimElim(j,2) <<  endl; }
//----------------------------------------------------------------------
      }
//----------------------------------------------------------------------

}




