#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"
#include "../Header/scanSurface.h"

using namespace std;
using namespace arma;


int main(int argc, char** argv) {
  //----------------------------------------------------------------------------------------------------------------
  // Part 0) Setting Parameters
  //----------------------------------------------------------------------------------------------------------------
  
  // Flexible variables:
  double ao = atof(argv[4]); // [0.1 - 5.0] EXP. BEST = 0.5 5 3 5
  double rcut = atof(argv[5]); // [7.0 - 100.0] -> make sure to use fine grid with large cut
  int radialN = atoi(argv[6]); // [1-4]
  int lMax = atoi(argv[7]) ; // [1-9]
  int grid = atoi(argv[8]) ; // [1-9]

  // Non-Flexible variables:
  double pi = 3.14159265358979324;
  double halfPi = 3.14159265358979324*0.5;

  double rsc = pi*pi*0.5*0.5*rcut; // rescaleing the integration for gauss-legendre quaduature.

  double sig = 1.0; //Ghost variable, not used for now.

  double z = 1.0;
  double norm = pow(sqrt(z/ao),3);
  int typesN = 3;
  cube C = 100*ones<cube>(typesN,radialN,100); // (Type, n, all the coeffs) -- exactly 100 for l = 9
  cube intMea;
  cube intMeb;
  cube intAll;

//cout << "AAA" << endl;  
  //----------------------------------------------------------------------------------------------------------------
  // Part 1) Retrieving Data -> W, Gauss-Legendre, XYZ-Smeared
  //----------------------------------------------------------------------------------------------------------------
  
  mat GL; // [http://keisan.casio.com/exec/system/1329114617 (June 5th 2017)] , produced by Octave. W(:,0) -> GL coord. pos. W(:,1) -> GL weights.
  cube X;
  cube Y;
  cube Z;
  cube GLC(GL.n_rows,GL.n_rows,GL.n_rows); // GL weights in 3D which is just an outer product of GL;
//  if(rcut <= 2.5){GL.load("LegendreWeights/parameters20.txt");GLC.load("BI/GLC20.bi");} //
//  else if (rcut >2.5 && rcut <= 3.5 ){GL.load("LegendreWeights/parameters30.txt");GLC.load("BI/GLC30.bi");} //
//  else if (rcut >3.5 && rcut <= 4.5 ){GL.load("LegendreWeights/parameters40.txt");GLC.load("BI/GLC40.bi");} //
//  else if (rcut >4.5 && rcut <= 5.5 ){GL.load("LegendreWeights/parameters50.txt");GLC.load("BI/GLC50.bi");} //
//  else if (rcut >5.5 && rcut <= 6.5 ){GL.load("LegendreWeights/parameters60.txt");  GLC.load("BI/GLC60.bi");}//
//  else if (rcut >6.5 && rcut <= 7.5 ){GL.load("LegendreWeights/parameters70.txt");  GLC.load("BI/GLC70.bi");}//
//  else if (rcut >7.5 && rcut <= 8.5 ){GL.load("LegendreWeights/parameters80.txt");GLC.load("BI/GLC80.bi");} //
//  else if (rcut >8.5 && rcut <= 9.5 ){GL.load("LegendreWeights/parameters90.txt");GLC.load("BI/GLC90.bi");} //
//  else if (rcut >9.5 && rcut <= 10.5 ){ GL.load("LegendreWeights/parameters100.txt");GLC.load("BI/GLC100.bi");} //
////  else if (rcut >9.5 && rcut <= 10.5 ){GL.load("parameters100.txt");GLC.load("BI/GLC100.bi");}//
//  else    {cout << "Error::rcut too large..\n Exiting.." << cout; exit(0);}
 
  if(grid == 0 ){GL.load("LegendreWeights/parameters20.txt");GLC.load("BI/GLC20.bi");} //
  else if (grid == 1){GL.load("LegendreWeights/parameters30.txt");GLC.load("BI/GLC30.bi");} //
  else if (grid == 2){GL.load("LegendreWeights/parameters40.txt");GLC.load("BI/GLC40.bi");} //
  else if (grid == 3){GL.load("LegendreWeights/parameters50.txt");GLC.load("BI/GLC50.bi");} //
  else if (grid == 4){GL.load("LegendreWeights/parameters60.txt");  GLC.load("BI/GLC60.bi");}//
  else if (grid == 5){GL.load("LegendreWeights/parameters70.txt");  GLC.load("BI/GLC70.bi");}//
  else if (grid == 6){GL.load("LegendreWeights/parameters80.txt");GLC.load("BI/GLC80.bi");} //
  else if (grid == 7){GL.load("LegendreWeights/parameters90.txt");GLC.load("BI/GLC90.bi");} //
  else if (grid == 8){ GL.load("LegendreWeights/parameters100.txt");GLC.load("BI/GLC100.bi");} //
  else    {cout << "Error::grid too large..\n Exiting.." << cout; exit(0);}

// Getting R, Theta and Phi rescaled for the Gaull-Legendre quadrature
  vec R = rcut*0.5*GL.col(0) + rcut*0.5 ;
  vec The = pi*GL.col(0)*0.5 + pi*0.5;
  vec Phi = pi*GL.col(0) + pi;       

//Setting Cartesian coordinates from the R, Theta and Phi. 
  X = getSphericalToCartCubeX( R, The, Phi);
  Y = getSphericalToCartCubeY( R, The, Phi);
  Z = getSphericalToCartCubeZ( R, The, Phi);

//----Holding memory places.----
  double sumMe = 0;
  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec typeA = zeros<vec>(coord.n_rows);
  vec typeB = zeros<vec>(coord.n_rows);
  coord = posAve(coord); 

  mat coordReload = coord;
//cout << "CAA" << endl;  

  double newJ = 0;
  mat gn(R.n_rows,3);// Radial Basis Functions [APB eq.25]. g(*,:) -> n's. g(:,*) -> r's of GL coord. pos. Slater used.
  mat  g;
  mat Y0 = zeros<mat>(GL.n_rows, GL.n_rows);// Tesseral Spherical Harmonics at GL coord. pos. first lMax is l, second lMax is m but patted with 0's. 
  cube rho_a;
  cube rho_b;
  cube rhoAll;
  int globalI = 0;
  int incrementN = 0; 
  double P[typesN][radialN][radialN][lMax]; // Power Spectrum P[A-type][n1][n2][l]
// New coordinates -> type A + Hydrogen and type B + Hydrogen

  gn.col(0)= hydrogenRDF(1, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(1)= hydrogenRDF(2, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(2)= hydrogenRDF(3, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
//  gn.col(3)= hydrogenRDF(4, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)

  g = gn.t();

  for(int t=0; t < GL.n_rows; t++){ 
    for(int p=0; p < GL.n_rows; p++){ 
       Y0.at(t,p) = tesseral_spherical_harm(0,0,The.at(t),Phi.at(p)); // rescaled
      }
    }

  cube Y1 = getY(1,The, Phi);
  cube Y2 = getY(2,The, Phi);
  cube Y3 = getY(3,The, Phi);
  cube Y4 = getY(4,The, Phi);
  cube Y5 = getY(5,The, Phi);
  cube Y6 = getY(6,The, Phi);
  cube Y7 = getY(7,The, Phi);
  cube Y8 = getY(8,The, Phi);
  cube Y9 = getY(9,The, Phi);

  vec Cbuff(3); // for type a and type b

//Finding where atom A and atom B are in .xyz.
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){typeA(i) =1; }
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){typeB(i) =1;}
   }

  mat coord_a = zeros<mat>(sum(typeA) + 1,3);
  mat coord_b = zeros<mat>(sum(typeB) + 1,3);
  mat coordH = zeros<mat>(coord.n_rows + 1, 3);


// FOR For
for(int bigI = 0; bigI < coord.n_rows; bigI++){
// if(bigI != 0 && bigI != 6 && bigI != 161) continue;
   coord = coordReload;
//  coord.print("cord");

//  coord(bigI,0) = -1000000;   coord(bigI,1) = -1000000;   coord(bigI,2) = -1000000; 

  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){coord_a.row(newJ) = coord.row(i); newJ++;}
   }

  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){coord_b.row(newJ) = coord.row(i); newJ++;}
   }

  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     coordH.row(i) = coord.row(i); newJ++;
   }

// Not adding Hydrogen at the end, but adding the loop position of atoms.
  coord_a.row(coord_a.n_rows - 1) = coordReload.row(bigI);
  coord_b.row(coord_b.n_rows - 1) = coordReload.row(bigI);
  coordH.row(coordH.n_rows - 1) = coordReload.row(bigI);

// Gaussian Smearing at Gauss-Lgedendre quadratue points
  rho_a  = getGaussDistr(coord_a,R, The, Phi, X, Y, Z, sig);
  rho_b  = getGaussDistr(coord_b,R, The, Phi, X, Y, Z, sig);
  rhoAll = getGaussDistr(coordH ,R, The, Phi, X, Y, Z, sig);


  //----------------------------------------------------------------------------------------------------------------
  // Part 3) get coefs c_anlm by Integration ->  where C(a,n,I) where a is the type, n is the radial basis function
  //         and I is (l,m) sequenially. intMa and intMb multiplied with getT are integrands.
  //----------------------------------------------------------------------------------------------------------------

  intMea = rho_a%GLC;
  intMeb = rho_b%GLC;
  intAll = rhoAll%GLC;


  globalI = 0;

  for(int n=0; n < radialN; n++) {
    Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getTMat(n,0,0,g,Y0,R,The,Phi));
    for(int i=0; i < typesN; i++){
       C(i,n,0) = Cbuff(i);
    }
  }


for(int l = 1; l <= lMax; l++){
  for(int m=-l; m <= l; m++) {
    globalI++;
    for(int n=0; n < radialN; n++) {
        if(l == 1){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y1,R,The,Phi));}
        else if(l == 2){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y2,R,The,Phi));}
        else if(l == 3){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y3,R,The,Phi));}
        else if(l == 4){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y4,R,The,Phi));}
        else if(l == 5){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y5,R,The,Phi));}
        else if(l == 6){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y6,R,The,Phi));}
        else if(l == 7){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y7,R,The,Phi));}
        else if(l == 8){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y8,R,The,Phi));}
        else if(l == 9){Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getT(n,l,m,g,Y9,R,The,Phi));}
        else{cout << "ERROR, l too large" << endl;}
      for(int i=0; i<typesN;i++){
      C(i,n,globalI) = Cbuff(i);
      }
    }
  }
}

  //----------------------------------------------------------------------------------------------------------------
  // Part 4) get Power Spectrum
  //----------------------------------------------------------------------------------------------------------------
  
  incrementN = 0; 

  for(int a=0; a < typesN; a++){ // Types + All
    for(int n1=0; n1 < radialN; n1++){  
      for(int n2=0; n2 < radialN; n2++){ 

          incrementN = 0;

        for(int l=0; l <= lMax; l++){ 
          sumMe = 0;
          for(int m=-l; m <= l; m++){ 
            sumMe += C(a,n1,incrementN)*C(a,n2,incrementN);
            incrementN++;
          }
            cout << sumMe << " ";
        }
      }
    }
  }
cout << endl;
//cout << bigI << endl;
}
return 0;
}
