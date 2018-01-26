#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

  uvec q ;
  q.load(argv[1]); 
  mat finalHs = zeros<mat>(q.n_rows, 3);

  int x = 0;

  double coulmbNew = 0;
  double coulmbOld = 100000000000000;
  vec distSqr = zeros<mat>(1);

  mat coord = getPos(argv[2]);
  string* type = getType(argv[2]);

  rowvec bufvec = zeros<rowvec>(3);
  rowvec bufvec2 = zeros<rowvec>(3);
  rowvec bufvecBest = zeros<rowvec>(3);

  vec the = 3.1416*linspace<vec>(0,3.1416/2.0,45);
  vec phi = 3.1416*linspace<vec>(0,3.1416,90);


  double radius = atof(argv[3]);

  for(int p = 0; p < q.n_elem; p++){
    coulmbNew = 0;
    coulmbOld = 100000000000000;
    for(int i = 0; i < q.n_elem; i++){
      for(int j = 0; j < the.n_elem ; j++){
        for(int k = 0; k <  phi.n_elem; k++){

          bufvec(0) = sin(the.at(j))*cos(phi.at(k));
          bufvec(1) = sin(the.at(j))*sin(phi.at(k));
          bufvec(2) = cos(the.at(j));
          bufvec = radius*bufvec + coord.row(q(p) - 1);
          coulmbNew = 0;

          for(int x = 0; x < coord.n_rows; x++){
             bufvec2 = bufvec - coord.row(x); 
             distSqr = bufvec2*bufvec2.t();
             coulmbNew += 1/distSqr.at(0);
          }

          if(coulmbNew < coulmbOld){
             coulmbOld = coulmbNew;
             bufvecBest = bufvec;
          }

        }

      }

    }

    bufvecBest.print("");

  }

        
return 0;
}


