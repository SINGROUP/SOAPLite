#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"

using namespace std;
using namespace arma;


int main(int argc, char** argv) {

  mat soapMatCol;
  soapMatCol.load(argv[1]); 
  mat soapMat = soapMatCol; 
  vec stillExists = ones<mat>(soapMat.n_rows);
  rowvec multMe = zeros<rowvec>(soapMat.n_cols);
  vec buffer(1);
//  vec uniqueList = -ones<mat>(soapMat.n_rows);
  int x = 0;

  double Ecut = atof(argv[2]); 
  vec soapDiffs(soapMat.n_rows);
  for(int j=0; j < soapMat.n_rows ; j++){
    if(stillExists(j) > 0.5){
      cout << j+1 << endl;
      for(int i=0; i < soapMat.n_rows ; i++){
          multMe = soapMat.row(j) - soapMat.row(i);
          buffer = sqrt(multMe*multMe.t())/soapMat.n_cols;
          soapDiffs(i) = buffer(0);

      }
    }
  uvec q = find(soapDiffs <= Ecut);

  for(int i=0; i < q.n_elem ; i++){
    stillExists(q(i)) = 0.0;
  }

  }

return 0;
}


