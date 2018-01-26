#include<iostream>
#include<armadillo>

#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"
#include "../Header/scanSurface.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {
  int selectSize = 100000;
  mat hydrogenSOAP;
  mat hydrogenPos;

  hydrogenSOAP.load(argv[1]);
  hydrogenPos.load(argv[2]);

  double Ecut = atof(argv[3]); 
  int Nhyd = atof(argv[4]); 

  int hSSiz =  hydrogenSOAP.n_rows;
  int randInt = 0;
// cout << "CCA" << endl;
  
  rowvec buffRow(Nhyd*hydrogenSOAP.n_cols);
  vec buffVal(1);
  rowvec buffSOAP(hydrogenSOAP.n_cols);
  imat hydConfigs(selectSize,Nhyd);
  int uppTensorSizeAcc = 0; 

  irowvec indxVec = linspace<irowvec>(0, hSSiz - 1, hSSiz);
//  indxVec.print("ind");
  imat randNumVec = randi<imat>(selectSize, Nhyd);

//    for(int i = 0; i < selectSize; i++){
//       for(int n = 0; n < Nhyd; n++){
//         randNumVec(i,n) = randNumVec(i,n)%hydrogenSOAP.n_rows; 
//       }
//    }
// cout << "CAA" << endl;

  for(int i = 0; i < selectSize; i++){
    for(int s = 0; s < Nhyd; s++){
      indxVec = shuffle(indxVec);
         randNumVec.at(i,s) = indxVec[s]; 
//         cout << indxVec(s) << endl;
    }
  }

//  randNumVec.print("");
// cout << "AAA" << endl;

//  mat combinedAverage(selectSize, Nhyd*hydrogenSOAP.n_cols);
  mat combinedAverage = zeros<mat>(selectSize, hydrogenSOAP.n_cols);
  
//getting upper tensor size
  for(int k = 0; k < selectSize; k++){
    for(int n = 0; n < Nhyd; n++){
        randInt = randNumVec(k,n);
        hydConfigs(k,n) = randInt;
//        buffSOAP = hydrogenSOAP.row(randInt);
        combinedAverage.row(k) = combinedAverage.row(k) + hydrogenSOAP.row(randInt);
//      for(int x = 0; x < hydrogenSOAP.n_cols; x++){
        
//        combinedAverage(k,n*Nhyd + x) = buffSOAP(x);
//        k++;
//      }
    }
  }

        combinedAverage = combinedAverage/Nhyd;
// cout << "AAB" << endl;

  vec stillExists = ones<vec>(combinedAverage.n_rows);
  vec soapDiffs(combinedAverage.n_rows);

//  k = 0;

//  for(int l=0; l < hydrogenSOAP.n_rows; l++){
//    for(int m=l + 1; m < hydrogenSOAP.n_rows; m++){
//      for(int n=m + 1; n < hydrogenSOAP.n_rows; n++){
    for(int k=0; k < selectSize; k++){
//      cout << "k: " << k << endl;
    
        if(stillExists(k) > 0.5){

          for(int l=0; l < Nhyd; l++){

            cout << hydrogenPos.row(hydConfigs(k,l)) << " ";

          }
// cout << "AAC" << endl;
          cout << endl;

          for(int i=0; i < combinedAverage.n_rows ; i++){

            buffRow = combinedAverage.row(k) - combinedAverage.row(i);
            buffVal = sqrt(buffRow*buffRow.t())/hydrogenSOAP.n_cols;
            soapDiffs(i) = buffVal(0);

          }
        }
// cout << "AAD" << endl;
        uvec q = find(soapDiffs <= Ecut);
  
// cout << "AAE" << endl;
        for(int i=0; i < q.n_elem ; i++){
          stillExists(q(i)) = 0.0;
        }
// cout << "AAF" << endl;

      }
//    }

//  }

return 0;
}

