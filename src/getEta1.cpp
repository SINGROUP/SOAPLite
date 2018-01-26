#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"
#include "../Header/scanSurface.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

// Non-Flexible variables:
  double pi = 3.14159265358979324;
  double halfPi = 3.14159265358979324*0.5;

// Getting R, Theta and Phi rescaled for the Gaull-Legendre quadrature
  vec R = linspace<vec>(1.1 ,2.1,30);
  vec The =linspace<vec>(0,pi,90);
  vec Phi =linspace<vec>(0,2*pi,180);
  double dist1 = atof(argv[4]);
  double dist2 = atof(argv[5]);
  double epsilon=0.01;
  double neighEpsi=0.01;
  double potDiff = 0.00;
  double bubble = atof(argv[6]);

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec bufval(1);
  vec newPot = zeros<vec>(1);
  vec bestPot = zeros<vec>(1);
  newPot(0) = 0;
  bestPot(0) = 1e100;
 
  vec allDist=zeros<vec>(1);

  uvec q =  getSurf(coord,bubble);


  vec typeA = zeros<vec>(coord.n_rows);
  vec typeB = zeros<vec>(coord.n_rows);
//  coord = posAve(coord); 

  rowvec buffvec1 = zeros<rowvec>(3);
  rowvec buffvec2 = zeros<rowvec>(3);
  rowvec buffvecR = zeros<rowvec>(3);

  rowvec DVecNew = zeros<rowvec>(3);
  rowvec DVecOld = zeros<rowvec>(3);

  rowvec BestVec = zeros<rowvec>(3);


//Finding where atom A and atom B are in .xyz.
  for(int i=0; i < q.n_elem; i++)  { 
     if(type[q(i)] == argv[2] ){typeA(q(i)) =1;}
   }
  for(int i=0; i < q.n_elem; i++)  { 
     if(type[q(i)] == argv[3] ){typeB(q(i)) =1;}
   }

// New coordinates -> type A + Hydrogen and type B + Hydrogen
  mat coord_a = zeros<mat>(sum(typeA) ,3);
  mat coord_b = zeros<mat>(sum(typeB) ,3);

  double newJ = 0;
  for(int i=0; i < q.n_rows; i++)  { 
     if(type[q(i)] == argv[2] ){coord_a.row(newJ) = coord.row(q(i)); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < q.n_rows; i++)  { 
     if(type[q(i)] == argv[3] ){coord_b.row(newJ) = coord.row(q(i)); newJ++;}
   }
  newJ = 0;


//------------------------------------------------------------------------------------------------------------
  for(int i=0; i < coord_a.n_rows; i++)  { 
      buffvec1 = coord_a.row(i);
        for(int t=0; t < The.n_elem; t++)  { 
          for(int p=0; p < Phi.n_elem; p++)  { 

            buffvecR(0)= sin(The(t))*cos(Phi(p)); 
            buffvecR(1)= sin(The(t))*sin(Phi(p)); 
            buffvecR(2)= cos(The(t)); 
            DVecNew = buffvec1 + dist1*buffvecR;
             
            for(int x=0; x < coord.n_rows; x++)  { 
              allDist = sqrt((DVecNew- coord.row(x))*(DVecNew - coord.row(x)).t());

              newPot = newPot + 1/(allDist%allDist);
            }

            if(newPot(0) < bestPot(0)){
              bestPot = newPot;
              BestVec = DVecNew;
            }
           newPot(0)= 0;
          }
        }
           bestPot(0)= 1e100;
           BestVec.print();
    }
//------------------------------------------------------------------------------------------------------------
  for(int i=0; i < coord_b.n_rows; i++)  { 
      buffvec1 = coord_b.row(i);
        for(int t=0; t < The.n_elem; t++)  { 
          for(int p=0; p < Phi.n_elem; p++)  { 

            buffvecR(0)= sin(The(t))*cos(Phi(p)); 
            buffvecR(1)= sin(The(t))*sin(Phi(p)); 
            buffvecR(2)= cos(The(t)); 
            DVecNew = buffvec1 + dist2*buffvecR;
             
            for(int x=0; x < coord.n_rows; x++)  { 
              allDist = sqrt((DVecNew- coord.row(x))*(DVecNew - coord.row(x)).t());
              newPot = newPot + 1/(allDist%allDist);
            }

            if(newPot(0) < bestPot(0)){
              bestPot = newPot;
              BestVec = DVecNew;
            }
           newPot(0)= 0;
          }
        }
           bestPot(0)= 1e100;
           BestVec.print();
    }
//------------------------------------------------------------------------------------------------------------

return 0;
}
















