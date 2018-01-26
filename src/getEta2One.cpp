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
  vec The =linspace<vec>(0,pi,180);
  vec Phi =linspace<vec>(0,2*pi,360);
  double dist1 = atof(argv[2]);
  double epsilon=0.01;
  double neighEpsi=0.01;
  double potDiff = 0.00;
  double bubble = atof(argv[3]);

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec bufval(1);
  vec newPot = zeros<vec>(1);
  vec bestPot = zeros<vec>(1);
  newPot(0) = 0;
  bestPot(0) = 1e100;
 
  vec allDist=zeros<vec>(1);

  uvec q =  getSurf(coord,bubble);

//  coord = posAve(coord); 

  rowvec buffvec1 = zeros<rowvec>(3);
  rowvec buffvec2 = zeros<rowvec>(3);
  rowvec buffvecR = zeros<rowvec>(3);

  rowvec DVecNew = zeros<rowvec>(3);
  rowvec DVecOld = zeros<rowvec>(3);

  rowvec BestVec = zeros<rowvec>(3);

  for(int i=0; i < q.n_rows; i++)  { 
    for(int j=i+1; j < q.n_rows; j++)  { 
      buffvec1 = coord.row(q(i));
      buffvec2 = coord.row(q(j));
      bufval = sqrt((buffvec1- buffvec2)*(buffvec1 - buffvec2).t());
      if(bufval(0) <3){
        for(int t=0; t < The.n_elem; t++)  { 
          for(int p=0; p < Phi.n_elem; p++)  { 

            buffvecR(0)= sin(The(t))*cos(Phi(p)); 
            buffvecR(1)= sin(The(t))*sin(Phi(p)); 
            buffvecR(2)= cos(The(t)); 
            DVecNew = buffvec1 + dist1*buffvecR;
            bufval = sqrt((DVecNew- buffvec2)*(DVecNew - buffvec2).t());
             
            for(int x=0; x < coord.n_rows; x++)  { 
              allDist = sqrt((DVecNew- coord.row(x))*(DVecNew - coord.row(x)).t());

              newPot = newPot + 1/(allDist%allDist);
            }

            if(bufval(0) < (dist1 + epsilon) && bufval(0) > (dist1 - epsilon) && newPot(0) < bestPot(0)){
              bestPot = newPot;
              BestVec = DVecNew;
            }
           newPot(0)= 0;
          }
        }
           bestPot(0)= 1e100;
           BestVec.print();
      }
    }
  }

//------------------------------------------------------------------------------------------------------------

return 0;
}
















