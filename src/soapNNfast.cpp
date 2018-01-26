#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"

using namespace std;
using namespace arma;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a SOAP algorithm with real spherical harmonics with the radial basis function of polynomials larger than //
// order of 2. The algorithm is based on [On representing chemical environments by Albert P. Batrok et al.] and     // 
// [Comparing molecules and solids across structural and alchemical space by Sandip De et al.]. In the cose, the    //
// equations will be reffered as [APB eq.#] and [SD eq.#]. The armadillo package [http://arma.sourceforge.net/]  is //
// heavily used due to its Octave/Matlab/Numpy like syntax which accelerates the developement speed.                //
//   The aim is to attempt to use the power spectrum for Neural Network instead for KRR.                            //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // The algorithm is devided in 5 parts:                                                                          //
  // 1) Retrieving data. -> W[APB eq.25], Gauss-Legendre values for integration,                                   //
  // 2) Constructing Basis functions -> g_n(r)[APB eq.25], Y_lm(The,Phi)[Tesseral Harmonics], where The=Theta      //
  // 3) Get coeffs c_nlm[APD eq.24] = Integrate Rho(r,The,Phi) T dV, where Rho[SD eq.14] is the Gaussian smeared   //
  //                                    xyz atomic positions.                                                      //
  // 4) Get power spectrums by P_b1b2lab = sum_m (c_b1lma  c_b2lmb) [SD eq.17]                                     //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //----------------------------------------------------------------------------------------------------------------
  // Part 0) Setting Parameters
  //----------------------------------------------------------------------------------------------------------------
  
  // Flexible variables:
  double ao = atof(argv[4]); // [0.1 - 5.0] EXP. BEST = 0.5 5 3 5
  double rcut = atof(argv[5]); // [7.0 - 100.0] -> make sure to use fine grid with large cut
  int radialN = atoi(argv[6]); // [1-4]
  int lMax = atoi(argv[7]) ; // [1-9]
  double sumMe = 0;

  // Non-Flexible variables:
  double pi = 3.14159265358979324;
  double halfPi = 3.14159265358979324*0.5;

  double rsc = pi*pi*0.5*0.5*rcut; // rescaleing the integration for gauss-legendre quaduature.

  double sig = 1.0; //Ghost variable, not used for now.

  double z = 1.0;
  double norm = pow(sqrt(z/ao),3);
  int typesN = 3;
  
  //----------------------------------------------------------------------------------------------------------------
  // Part 1) Retrieving Data -> W, Gauss-Legendre, XYZ-Smeared
  //----------------------------------------------------------------------------------------------------------------
  
  mat GL; // [http://keisan.casio.com/exec/system/1329114617 (June 5th 2017)] , produced by Octave. W(:,0) -> GL coord. pos. W(:,1) -> GL weights.
  cube X;
  cube Y;
  cube Z;
  cube GLC(GL.n_rows,GL.n_rows,GL.n_rows); // GL weights in 3D which is just an outer product of GL;
//  if(rcut <= 1.5){GL.load("LegendreWeights/parameters20.txt");GLC.load("BI/GLC20.bi");} //
//  else if (rcut >1.5 && rcut <= 2.5 ){GL.load("LegendreWeights/parameters30.txt");GLC.load("BI/GLC30.bi");} //
//  else if (rcut >2.5 && rcut <= 3.5 ){GL.load("LegendreWeights/parameters40.txt");GLC.load("BI/GLC40.bi");} //
//  else if (rcut >3.5 && rcut <= 4.5 ){GL.load("LegendreWeights/parameters50.txt");GLC.load("BI/GLC50.bi");} //
//  else if (rcut >4.5 && rcut <= 5.5 ){GL.load("LegendreWeights/parameters60.txt");  GLC.load("BI/GLC60.bi");}//
//  else if (rcut >5.5 && rcut <= 6.5 ){GL.load("LegendreWeights/parameters70.txt");  GLC.load("BI/GLC70.bi");}//
//  else if (rcut >6.5 && rcut <= 7.5 ){GL.load("LegendreWeights/parameters80.txt");GLC.load("BI/GLC80.bi");} //
//  else if (rcut >7.5 && rcut <= 8.5 ){GL.load("LegendreWeights/parameters90.txt");GLC.load("BI/GLC90.bi");} //
//  else if (rcut >8.5 && rcut <= 9.5 ){ GL.load("LegendreWeights/parameters100.txt");GLC.load("BI/GLC100.bi");} //
//  else if (rcut >9.5 && rcut <= 10.5 ){GL.load("LegendreWeights/parameters100.txt");GLC.load("BI/GLC100.bi");}//
//  else    {cout << "Error::rcut too large..\n Exiting.." << cout; exit(0);}

  GL.load("LegendreWeights/parameters50.txt");GLC.load("BI/GLC50.bi");//

// Getting R, Theta and Phi rescaled for the Gaull-Legendre quadrature
  vec R = rcut*0.5*GL.col(0) + rcut*0.5 ;
  vec The = pi*GL.col(0)*0.5 + pi*0.5;
  vec Phi = pi*GL.col(0) + pi;       

//Setting Cartesian coordinates from the R, Theta and Phi. 
  X = getSphericalToCartCubeX( R, The, Phi);
  Y = getSphericalToCartCubeY( R, The, Phi);
  Z = getSphericalToCartCubeZ( R, The, Phi);

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);
  vec typeA = zeros<vec>(coord.n_rows);
  vec typeB = zeros<vec>(coord.n_rows);
  coord = posAve(coord); 

// Checking the Soap Kernel by rotating it and changing the shape.
  //  coord.print("Before:");
//    coord = rotate3d(coord,1,1,1);
  //  coord.print("After:");
//    coord(coord.n_rows - 1,2) += 0.1;

//Finding where atom A and atom B are in .xyz.
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){typeA(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){typeB(i) =1;}
   }
// New coordinates -> type A + Hydrogen and type B + Hydrogen
  mat coord_a = zeros<mat>(sum(typeA) + 1,3);
  mat coord_b = zeros<mat>(sum(typeB) + 1,3);

  double newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[2] ){coord_a.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == argv[3] ){coord_b.row(newJ) = coord.row(i); newJ++;}
   }

// Adding Hydrogen at the end
  coord(coord.n_rows - 1,0) = coord(coord.n_rows - 1,0) + 1000.0;
  coord(coord.n_rows - 1,1) = coord(coord.n_rows - 1,1) + 1000.0;
  coord(coord.n_rows - 1,2) = coord(coord.n_rows - 1,2) + 1000.0;
  coord_a.row(coord_a.n_rows - 1) = coord.row(coord.n_rows - 1);
  coord_b.row(coord_b.n_rows - 1) = coord.row(coord.n_rows - 1);

// Gaussian Smearing at Gauss-Lgedendre quadratue points
  cube rho_a = getGaussDistr(coord_a,R, The, Phi, X, Y, Z, sig);
  cube rho_b = getGaussDistr(coord_b,R, The, Phi, X, Y, Z, sig);
  cube rhoAll = getGaussDistr(coord,R, The, Phi, X, Y, Z, sig);


  //----------------------------------------------------------------------------------------------------------------
  // Part 2) Constructing Basis Functions -> g_n(r),  and Yn(Theta, Phi) 
  //----------------------------------------------------------------------------------------------------------------
  
  mat gn(R.n_rows,4);// Radial Basis Functions [APB eq.25]. g(*,:) -> n's. g(:,*) -> r's of GL coord. pos. Slater used.
  gn.col(0)= hydrogenRDF(1, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(1)= hydrogenRDF(2, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(2)= hydrogenRDF(3, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(3)= hydrogenRDF(4, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)

// Printing out gn for debugging.
//  gn.save("gn.dat",raw_ascii);
//  R.save("R.dat",raw_ascii);
//  The.save("The.dat",raw_ascii);
//  Phi.save("Phi.dat",raw_ascii);
// Transpose for gn
  mat    g = gn.t();

  mat Y0 = zeros<mat>(GL.n_rows, GL.n_rows);// Tesseral Spherical Harmonics at GL coord. pos. first lMax is l, second lMax is m but patted with 0's. 
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

  //----------------------------------------------------------------------------------------------------------------
  // Part 3) get coefs c_anlm by Integration ->  where C(a,n,I) where a is the type, n is the radial basis function
  //         and I is (l,m) sequenially. intMa and intMb multiplied with getT are integrands.
  //----------------------------------------------------------------------------------------------------------------

  cube intMea = rho_a%GLC;
  cube intMeb = rho_b%GLC;
  cube intAll = rhoAll%GLC;
  cube C = 100*ones<cube>(typesN,radialN,100); // (Type, n, all the coeffs) -- exactly 100 for l = 9
  int globalI = 0;



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
//C.slice(0).save("C.save",raw_ascii);
//rho_b.save("rho_b.save",raw_ascii);
//rhoAll.save("rhoAll.save",raw_ascii);

  //----------------------------------------------------------------------------------------------------------------
  // Part 4) get Power Spectrum
  //----------------------------------------------------------------------------------------------------------------
  
 int incrementN = 0; 

  double P[typesN][radialN][radialN][lMax]; // Power Spectrum P[A-type][n1][n2][l]
  memset(P, 0.0, sizeof P);

  for(int a=0; a < typesN; a++){ // Types + All
    for(int n1=0; n1 < radialN; n1++){  
      for(int n2=0; n2 < radialN; n2++){ 

          incrementN = 0;

        for(int l=0; l <= lMax; l++){ 
          sumMe = 0;
          for(int m=-l; m <= l; m++){ 
//            P[a][n1][n2][l] += C(a,n1,incrementN)*C(a,n2,incrementN);
            sumMe += C(a,n1,incrementN)*C(a,n2,incrementN);
            incrementN++;
          }
//            cout << a  << " " <<n1 << " "  << n2 << " " << l << " " << P[a][n1][n2][l] << endl;
//            cout << P[a][n1][n2][l] << " \n ";
            cout << sumMe << " ";
        }
      }
    }
  }
cout << endl;
return 0;
}
















