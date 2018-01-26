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
  vec The =linspace<vec>(0,pi,450);
  vec Phi =linspace<vec>(0,2*pi,900);
  double dist1 = atof(argv[2]);
  double epsilon=0.01;
  double neighEpsi=0.01;
  double potDiff = 0.00;
  double bubble = atof(argv[3]);

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec bufval1(1);
  vec bufval2(1);
  vec bufval3(1);
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
  rowvec buffvec3 = zeros<rowvec>(3);

  rowvec buffvecR = zeros<rowvec>(3);

  rowvec DVec = zeros<rowvec>(3);

  rowvec BestVec = zeros<rowvec>(3);

//------------------------------------------------------------------------------------------------------------
//GOLD-GOLD-GOLD
  for(int i=0; i < q.n_rows; i++)  { 
    for(int j=i+1; j < q.n_rows; j++)  { 
      for(int k=j+1; k < q.n_rows; k++)  { 

        buffvec1 = coord.row(q(i));
        buffvec2 = coord.row(q(j));
        buffvec3 = coord.row(q(k));

        bufval1 = sqrt((buffvec1- buffvec2)*(buffvec1 - buffvec2).t());
        bufval2 = sqrt((buffvec1- buffvec3)*(buffvec1 - buffvec3).t());
        bufval3 = sqrt((buffvec2- buffvec3)*(buffvec2 - buffvec3).t());

        if(bufval1(0) < 3 && bufval2(0) < 3 && bufval3(0) < 3){
          for(int t=0; t < The.n_elem; t++)  { 
            for(int p=0; p < Phi.n_elem; p++)  { 

              buffvecR(0)= sin(The(t))*cos(Phi(p)); 
              buffvecR(1)= sin(The(t))*sin(Phi(p)); 
              buffvecR(2)= cos(The(t)); 

              DVec = buffvec1 + dist1*buffvecR;

              bufval1 = sqrt((DVec - buffvec2)*(DVec - buffvec2).t());
              bufval2 = sqrt((DVec - buffvec3)*(DVec - buffvec3).t());
//              bufval3 = sqrt((DVec - buffvec2)*(DVec - buffvec2).t());
             
              for(int x=0; x < coord.n_rows; x++)  { 

                allDist = sqrt((DVec- coord.row(x))*(DVec - coord.row(x)).t());
                newPot = newPot + 1/(allDist%allDist);

              }

              if(bufval1(0) < (dist1 + epsilon) && bufval1(0) > (dist1 - epsilon)
                 && bufval2(0) < (dist1 + epsilon) && bufval2(0) > (dist1 - epsilon)
//                 && bufval3(0) < (dist1 + epsilon) && bufval3(0) > (dist1 - epsilon)
                 && newPot(0) < bestPot(0)){
                    bestPot = newPot;
                    BestVec = DVec;
              }
             newPot(0)= 0;
            }
          }
             bestPot(0)= 1e100;
             BestVec.print();
        }
      }
    }
  }
//------------------------------------------------------------------------------------------------------------

return 0;
}
















