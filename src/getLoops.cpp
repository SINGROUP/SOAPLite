#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

// Non-Flexible variables:
  double pi = 3.14159265358979324;
  double halfPi = 3.14159265358979324*0.5;

// Getting R, Theta and Phi rescaled for the Gaull-Legendre quadrature
  vec R = linspace<vec>(1.1 ,2.1,30);
  vec The =linspace<vec>(0,pi,180);
  vec Phi =linspace<vec>(0,2*pi,360);
  double dist1 = 1.8;
  double dist2 = 1.6;
  double epsilon=0.01;
  double neighEpsi=0.01;
  double potDiff = 0.00;

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec bufval(1);
  vec newPot = zeros<vec>(1);
  vec bestPot = zeros<vec>(1);


  vec typeA = zeros<vec>(coord.n_rows);
  vec typeB = zeros<vec>(coord.n_rows);
//  coord = posAve(coord); 

  rowvec buffvec1 = zeros<rowvec>(3);
  rowvec buffvec2 = zeros<rowvec>(3);
  rowvec buffvecR = zeros<rowvec>(3);

  rowvec DVecNew = zeros<rowvec>(3);
  rowvec DVecOld = zeros<rowvec>(3);


//Finding where atom A and atom B are in .xyz.
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){typeA(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){typeB(i) =1;}
   }

// New coordinates -> type A + Hydrogen and type B + Hydrogen
  mat coord_a = zeros<mat>(sum(typeA) ,3);
  mat coord_b = zeros<mat>(sum(typeB) ,3);

  double newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){coord_a.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){coord_b.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;


  for(int i=0; i < coord_a.n_rows; i++)  { 
    for(int j=i+1; j < coord_a.n_rows; j++)  { 
      buffvec1 = coord_a.row(i);
      buffvec2 = coord_a.row(j);
      bufval = sqrt((buffvec1- buffvec2)*(buffvec1 - buffvec2).t());
//      cout << bufval(0) << endl;
      if(bufval(0) <3){
        for(int t=0; t < The.n_elem; t++)  { 
          for(int p=0; p < Phi.n_elem; p++)  { 

            buffvecR(0)= sin(The(t))*cos(Phi(p)); 
            buffvecR(1)= sin(The(t))*sin(Phi(p)); 
            buffvecR(2)= cos(The(t)); 
 //           DVecOld = DVecNew;
            DVecNew = buffvec1 + dist1*buffvecR;
            bufval = sqrt((DVecNew- buffvec2)*(DVecNew - buffvec2).t());
             
            if(bufval(0) < (dist1 + epsilon) && bufval(0) > (dist1 - epsilon)){
              DVecNew.print(); 
            }
          }
        }
      }
           newPot= zeros<vec>(1);
           bestPot= zeros<vec>(1);
    }
  }


return 0;
}
















