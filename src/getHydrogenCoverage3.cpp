#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"
#include "../Header/scanSurface.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

  mat hydrogenSOAP;
  mat hydrogenPos;

  hydrogenSOAP.load(argv[1]);
  hydrogenPos.load(argv[2]);

  double Ecut = atof(argv[3]); 

  int hSSiz =  hydrogenSOAP.n_rows;
  
  rowvec buffRow(2*hydrogenSOAP.n_cols);
  vec buffVal(1);
  int uppTensorSizeAcc = 0; 
  int k = 0;

  for(int l = 0; l < hydrogenSOAP.n_rows; l++){
    for(int m = l + 1; m < hydrogenSOAP.n_rows; m++){
      for(int n = m + 1; n < hydrogenSOAP.n_rows; n++){
        uppTensorSizeAcc++;
      }
    }
  }

  mat combinedH1(uppTensorSizeAcc, 3*hydrogenSOAP.n_cols);

//getting upper tensor size
  for(int l = 0; l < hydrogenSOAP.n_rows; l++){
    for(int m = l + 1; m < hydrogenSOAP.n_rows; m++){
      for(int n = m + 1; n < hydrogenSOAP.n_rows; n++){
        combinedH1.row(k)=  join_horiz(join_horiz(hydrogenSOAP.row(l),hydrogenSOAP.row(m)),hydrogenSOAP.row(n));
        k++;
      }
    }
  }


  vec stillExists = ones<vec>(combinedH1.n_rows);
  vec soapDiffs(combinedH1.n_rows);

  k = 0;

  for(int l=0; l < hydrogenSOAP.n_rows; l++){
    for(int m=l + 1; m < hydrogenSOAP.n_rows; m++){
      for(int n=m + 1; n < hydrogenSOAP.n_rows; n++){
    
        if(stillExists(k) > 0.5){
          cout << hydrogenPos.row(l) << "" << hydrogenPos.row(m) << "" << hydrogenPos.row(n) <<  endl;
  
          for(int i=0; i < combinedH1.n_rows ; i++){
            buffRow = combinedH1.row(k) - combinedH1.row(i);
            buffVal = sqrt(buffRow*buffRow.t())/hydrogenSOAP.n_cols;
            soapDiffs(i) = buffVal(0);
          }
        }
        uvec q = find(soapDiffs <= Ecut);
  
        for(int i=0; i < q.n_elem ; i++){
          stillExists(q(i)) = 0.0;
        }
        k++;
      }
    }
//    cout << l << endl;
  }

return 0;
}
















