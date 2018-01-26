#include<iostream>
#include<armadillo>
#include"../Header/myArmadillo.h"
#include"../Header/myMath.h"
#include"../Header/mySoap.h"

using namespace std;
using namespace arma;


int main(int argc, char** argv) {
  //----------------------------------------------------------------------------------------------------------------
  // Part 0) Setting Parameters
  //----------------------------------------------------------------------------------------------------------------
  
  // Flexible variables:
  mat hydrogenPos;
  hydrogenPos.load(argv[2]); // [0.1 - 5.0] EXP. BEST = 0.5 5 3 5
  double ao = atof(argv[3]); // [0.1 - 5.0] EXP. BEST = 0.5 5 3 5
  double rcut = atof(argv[4]); // [7.0 - 100.0] -> make sure to use fine grid with large cut
  int radialN = 3;// atoi(argv[5]); // [1-4]
  int lMax = 9 ;//atoi(argv[6]) ; // [1-9]
  int grid = atoi(argv[5]);//atoi(argv[7]) ; // [1-9]

  // Non-Flexible variables:
  double pi = 3.14159265358979324;
  double halfPi = 3.14159265358979324*0.5;

  double rsc = pi*pi*0.5*0.5*rcut; // rescaleing the integration for gauss-legendre quaduature.

  double sig = 1.0; //Ghost variable, not used for now.

  double z = 1.0;
  double norm = pow(sqrt(z/ao),3);
  int typesN = 95;

//cout << "AAA" << endl;  
  //----------------------------------------------------------------------------------------------------------------
  // Part 1) Retrieving Data -> W, Gauss-Legendre, XYZ-Smeared
  //----------------------------------------------------------------------------------------------------------------
  
  mat GL; // [http://keisan.casio.com/exec/system/1329114617 (June 5th 2017)] , produced by Octave. W(:,0) -> GL coord. pos. W(:,1) -> GL weights.
  cube X;
  cube Y;
  cube Z;
  cube GLC(GL.n_rows,GL.n_rows,GL.n_rows); 

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
//  vec typeA = zeros<vec>(coord.n_rows);
//  vec typeB = zeros<vec>(coord.n_rows);
    vec type_H = zeros<vec>(coord.n_rows); 
    vec type_He = zeros<vec>(coord.n_rows); 
    vec type_Li = zeros<vec>(coord.n_rows); 
    vec type_Be = zeros<vec>(coord.n_rows); 
    vec type_B = zeros<vec>(coord.n_rows); 
    vec type_C = zeros<vec>(coord.n_rows); 
    vec type_N = zeros<vec>(coord.n_rows); 
    vec type_O = zeros<vec>(coord.n_rows); 
    vec type_F = zeros<vec>(coord.n_rows); 
    vec type_Ne = zeros<vec>(coord.n_rows); 
    vec type_Na = zeros<vec>(coord.n_rows); 
    vec type_Mg = zeros<vec>(coord.n_rows); 
    vec type_Al = zeros<vec>(coord.n_rows); 
    vec type_Si = zeros<vec>(coord.n_rows); 
    vec type_P = zeros<vec>(coord.n_rows); 
    vec type_S = zeros<vec>(coord.n_rows); 
    vec type_Cl = zeros<vec>(coord.n_rows); 
    vec type_Ar = zeros<vec>(coord.n_rows); 
    vec type_K = zeros<vec>(coord.n_rows); 
    vec type_Ca = zeros<vec>(coord.n_rows); 
    vec type_Sc = zeros<vec>(coord.n_rows); 
    vec type_Ti = zeros<vec>(coord.n_rows); 
    vec type_V = zeros<vec>(coord.n_rows); 
    vec type_Cr = zeros<vec>(coord.n_rows); 
    vec type_Mn = zeros<vec>(coord.n_rows); 
    vec type_Fe = zeros<vec>(coord.n_rows); 
    vec type_Co = zeros<vec>(coord.n_rows); 
    vec type_Ni = zeros<vec>(coord.n_rows); 
    vec type_Cu = zeros<vec>(coord.n_rows); 
    vec type_Zn = zeros<vec>(coord.n_rows); 
    vec type_Ga = zeros<vec>(coord.n_rows); 
    vec type_Ge = zeros<vec>(coord.n_rows); 
    vec type_As = zeros<vec>(coord.n_rows); 
    vec type_Se = zeros<vec>(coord.n_rows); 
    vec type_Br = zeros<vec>(coord.n_rows); 
    vec type_Kr = zeros<vec>(coord.n_rows); 
    vec type_Rb = zeros<vec>(coord.n_rows); 
    vec type_Sr = zeros<vec>(coord.n_rows); 
    vec type_Y = zeros<vec>(coord.n_rows); 
    vec type_Zr = zeros<vec>(coord.n_rows); 
    vec type_Nb = zeros<vec>(coord.n_rows); 
    vec type_Mo = zeros<vec>(coord.n_rows); 
    vec type_Tc = zeros<vec>(coord.n_rows); 
    vec type_Ru = zeros<vec>(coord.n_rows); 
    vec type_Rh = zeros<vec>(coord.n_rows); 
    vec type_Pd = zeros<vec>(coord.n_rows); 
    vec type_Ag = zeros<vec>(coord.n_rows); 
    vec type_Cd = zeros<vec>(coord.n_rows); 
    vec type_In = zeros<vec>(coord.n_rows); 
    vec type_Sn = zeros<vec>(coord.n_rows); 
    vec type_Sb = zeros<vec>(coord.n_rows); 
    vec type_Te = zeros<vec>(coord.n_rows); 
    vec type_I = zeros<vec>(coord.n_rows); 
    vec type_Xe = zeros<vec>(coord.n_rows); 
    vec type_Cs = zeros<vec>(coord.n_rows); 
    vec type_Ba = zeros<vec>(coord.n_rows); 
    vec type_La = zeros<vec>(coord.n_rows); 
    vec type_Ce = zeros<vec>(coord.n_rows); 
    vec type_Pr = zeros<vec>(coord.n_rows); 
    vec type_Nd = zeros<vec>(coord.n_rows); 
    vec type_Pm = zeros<vec>(coord.n_rows); 
    vec type_Sm = zeros<vec>(coord.n_rows); 
    vec type_Eu = zeros<vec>(coord.n_rows); 
    vec type_Gd = zeros<vec>(coord.n_rows); 
    vec type_Tb = zeros<vec>(coord.n_rows); 
    vec type_Dy = zeros<vec>(coord.n_rows); 
    vec type_Ho = zeros<vec>(coord.n_rows); 
    vec type_Er = zeros<vec>(coord.n_rows); 
    vec type_Tm = zeros<vec>(coord.n_rows); 
    vec type_Yb = zeros<vec>(coord.n_rows); 
    vec type_Lu = zeros<vec>(coord.n_rows); 
    vec type_Hf = zeros<vec>(coord.n_rows); 
    vec type_Ta = zeros<vec>(coord.n_rows); 
    vec type_W = zeros<vec>(coord.n_rows); 
    vec type_Re = zeros<vec>(coord.n_rows); 
    vec type_Os = zeros<vec>(coord.n_rows); 
    vec type_Ir = zeros<vec>(coord.n_rows); 
    vec type_Pt = zeros<vec>(coord.n_rows); 
    vec type_Au = zeros<vec>(coord.n_rows); 
    vec type_Hg = zeros<vec>(coord.n_rows); 
    vec type_Tl = zeros<vec>(coord.n_rows); 
    vec type_Pb = zeros<vec>(coord.n_rows); 
    vec type_Bi = zeros<vec>(coord.n_rows); 
    vec type_Po = zeros<vec>(coord.n_rows); 
    vec type_At = zeros<vec>(coord.n_rows); 
    vec type_Rn = zeros<vec>(coord.n_rows); 
    vec type_Fr = zeros<vec>(coord.n_rows); 
    vec type_Ra = zeros<vec>(coord.n_rows); 
    vec type_Ac = zeros<vec>(coord.n_rows); 
    vec type_Th = zeros<vec>(coord.n_rows); 
    vec type_Pa = zeros<vec>(coord.n_rows); 
    vec type_U = zeros<vec>(coord.n_rows); 
    vec type_Np = zeros<vec>(coord.n_rows); 
    vec type_Pu = zeros<vec>(coord.n_rows); 
  hydrogenPos = posAveExt(coord,hydrogenPos); 
  coord = posAve(coord); 
  mat coordReload = coord;
//cout << "CAA" << endl;  
  cube C = 100*ones<cube>(typesN,radialN,100); // (Type, n, all the coeffs) -- exactly 100 for l = 9
//  cube intMea;
//  cube intMeb;
    cube intMeH; 
    cube intMeHe; 
    cube intMeLi; 
    cube intMeBe; 
    cube intMeB; 
    cube intMeC; 
    cube intMeN; 
    cube intMeO; 
    cube intMeF; 
    cube intMeNe; 
    cube intMeNa; 
    cube intMeMg; 
    cube intMeAl; 
    cube intMeSi; 
    cube intMeP; 
    cube intMeS; 
    cube intMeCl; 
    cube intMeAr; 
    cube intMeK; 
    cube intMeCa; 
    cube intMeSc; 
    cube intMeTi; 
    cube intMeV; 
    cube intMeCr; 
    cube intMeMn; 
    cube intMeFe; 
    cube intMeCo; 
    cube intMeNi; 
    cube intMeCu; 
    cube intMeZn; 
    cube intMeGa; 
    cube intMeGe; 
    cube intMeAs; 
    cube intMeSe; 
    cube intMeBr; 
    cube intMeKr; 
    cube intMeRb; 
    cube intMeSr; 
    cube intMeY; 
    cube intMeZr; 
    cube intMeNb; 
    cube intMeMo; 
    cube intMeTc; 
    cube intMeRu; 
    cube intMeRh; 
    cube intMePd; 
    cube intMeAg; 
    cube intMeCd; 
    cube intMeIn; 
    cube intMeSn; 
    cube intMeSb; 
    cube intMeTe; 
    cube intMeI; 
    cube intMeXe; 
    cube intMeCs; 
    cube intMeBa; 
    cube intMeLa; 
    cube intMeCe; 
    cube intMePr; 
    cube intMeNd; 
    cube intMePm; 
    cube intMeSm; 
    cube intMeEu; 
    cube intMeGd; 
    cube intMeTb; 
    cube intMeDy; 
    cube intMeHo; 
    cube intMeEr; 
    cube intMeTm; 
    cube intMeYb; 
    cube intMeLu; 
    cube intMeHf; 
    cube intMeTa; 
    cube intMeW; 
    cube intMeRe; 
    cube intMeOs; 
    cube intMeIr; 
    cube intMePt; 
    cube intMeAu; 
    cube intMeHg; 
    cube intMeTl; 
    cube intMePb; 
    cube intMeBi; 
    cube intMePo; 
    cube intMeAt; 
    cube intMeRn; 
    cube intMeFr; 
    cube intMeRa; 
    cube intMeAc; 
    cube intMeTh; 
    cube intMePa; 
    cube intMeU; 
    cube intMeNp; 
    cube intMePu; 



    cube intAll;

  double newJ = 0;
  mat gn(R.n_rows,4);// Radial Basis Functions [APB eq.25]. g(*,:) -> n's. g(:,*) -> r's of GL coord. pos. Slater used.
  mat  g;
  mat Y0 = zeros<mat>(GL.n_rows, GL.n_rows);// Tesseral Spherical Harmonics at GL coord. pos. first lMax is l, second lMax is m but patted with 0's. 
//  cube rho_a;
//  cube rho_b;
    cube rho_H; 
    cube rho_He; 
    cube rho_Li; 
    cube rho_Be; 
    cube rho_B; 
    cube rho_C; 
    cube rho_N; 
    cube rho_O; 
    cube rho_F; 
    cube rho_Ne; 
    cube rho_Na; 
    cube rho_Mg; 
    cube rho_Al; 
    cube rho_Si; 
    cube rho_P; 
    cube rho_S; 
    cube rho_Cl; 
    cube rho_Ar; 
    cube rho_K; 
    cube rho_Ca; 
    cube rho_Sc; 
    cube rho_Ti; 
    cube rho_V; 
    cube rho_Cr; 
    cube rho_Mn; 
    cube rho_Fe; 
    cube rho_Co; 
    cube rho_Ni; 
    cube rho_Cu; 
    cube rho_Zn; 
    cube rho_Ga; 
    cube rho_Ge; 
    cube rho_As; 
    cube rho_Se; 
    cube rho_Br; 
    cube rho_Kr; 
    cube rho_Rb; 
    cube rho_Sr; 
    cube rho_Y; 
    cube rho_Zr; 
    cube rho_Nb; 
    cube rho_Mo; 
    cube rho_Tc; 
    cube rho_Ru; 
    cube rho_Rh; 
    cube rho_Pd; 
    cube rho_Ag; 
    cube rho_Cd; 
    cube rho_In; 
    cube rho_Sn; 
    cube rho_Sb; 
    cube rho_Te; 
    cube rho_I; 
    cube rho_Xe; 
    cube rho_Cs; 
    cube rho_Ba; 
    cube rho_La; 
    cube rho_Ce; 
    cube rho_Pr; 
    cube rho_Nd; 
    cube rho_Pm; 
    cube rho_Sm; 
    cube rho_Eu; 
    cube rho_Gd; 
    cube rho_Tb; 
    cube rho_Dy; 
    cube rho_Ho; 
    cube rho_Er; 
    cube rho_Tm; 
    cube rho_Yb; 
    cube rho_Lu; 
    cube rho_Hf; 
    cube rho_Ta; 
    cube rho_W; 
    cube rho_Re; 
    cube rho_Os; 
    cube rho_Ir; 
    cube rho_Pt; 
    cube rho_Au; 
    cube rho_Hg; 
    cube rho_Tl; 
    cube rho_Pb; 
    cube rho_Bi; 
    cube rho_Po; 
    cube rho_At; 
    cube rho_Rn; 
    cube rho_Fr; 
    cube rho_Ra; 
    cube rho_Ac; 
    cube rho_Th; 
    cube rho_Pa; 
    cube rho_U; 
    cube rho_Np; 
    cube rho_Pu; 
    cube rhoAll;
  int globalI = 0;
  int incrementN = 0; 
  double P[typesN][radialN][radialN][lMax]; // Power Spectrum P[A-type][n1][n2][l]
// New coordinates -> type A + Hydrogen and type B + Hydrogen

//cout << "CAB" << endl;  
  gn.col(0)= hydrogenRDF(1, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(1)= hydrogenRDF(2, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(2)= hydrogenRDF(3, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)
  gn.col(3)= hydrogenRDF(4, z,ao,norm,R); // a0 = 0.5, Norm. Const. = 2^(3/2)

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

  vec Cbuff(95); // for type a and type b

//Finding where atom A and atom B are in .xyz.
//  for(int i=0; i < coord.n_rows; i++)  { 
//     if(type[i] == argv[2] ){typeA(i) =1;}
//   }
//  for(int i=0; i < coord.n_rows; i++)  { 
//     if(type[i] == argv[3] ){typeB(i) =1;}
//   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "H" ){type_H(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "He" ){type_He(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Li" ){type_Li(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Be" ){type_Be(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "B" ){type_B(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "C" ){type_C(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "N" ){type_N(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "O" ){type_O(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "F" ){type_F(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ne" ){type_Ne(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Na" ){type_Na(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mg" ){type_Mg(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Al" ){type_Al(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Si" ){type_Si(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "P" ){type_P(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "S" ){type_S(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cl" ){type_Cl(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ar" ){type_Ar(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "K" ){type_K(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ca" ){type_Ca(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sc" ){type_Sc(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ti" ){type_Ti(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "V" ){type_V(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cr" ){type_Cr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mn" ){type_Mn(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Fe" ){type_Fe(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Co" ){type_Co(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ni" ){type_Ni(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cu" ){type_Cu(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Zn" ){type_Zn(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ga" ){type_Ga(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ge" ){type_Ge(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "As" ){type_As(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Se" ){type_Se(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Br" ){type_Br(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Kr" ){type_Kr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rb" ){type_Rb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sr" ){type_Sr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Y" ){type_Y(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Zr" ){type_Zr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Nb" ){type_Nb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mo" ){type_Mo(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tc" ){type_Tc(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ru" ){type_Ru(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rh" ){type_Rh(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pd" ){type_Pd(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ag" ){type_Ag(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cd" ){type_Cd(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "In" ){type_In(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sn" ){type_Sn(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sb" ){type_Sb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Te" ){type_Te(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "I" ){type_I(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Xe" ){type_Xe(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cs" ){type_Cs(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ba" ){type_Ba(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "La" ){type_La(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ce" ){type_Ce(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pr" ){type_Pr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Nd" ){type_Nd(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pm" ){type_Pm(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sm" ){type_Sm(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Eu" ){type_Eu(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Gd" ){type_Gd(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tb" ){type_Tb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Dy" ){type_Dy(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ho" ){type_Ho(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Er" ){type_Er(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tm" ){type_Tm(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Yb" ){type_Yb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Lu" ){type_Lu(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Hf" ){type_Hf(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ta" ){type_Ta(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "W" ){type_W(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Re" ){type_Re(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Os" ){type_Os(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ir" ){type_Ir(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pt" ){type_Pt(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Au" ){type_Au(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Hg" ){type_Hg(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tl" ){type_Tl(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pb" ){type_Pb(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Bi" ){type_Bi(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Po" ){type_Po(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "At" ){type_At(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rn" ){type_Rn(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Fr" ){type_Fr(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ra" ){type_Ra(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ac" ){type_Ac(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Th" ){type_Th(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pa" ){type_Pa(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "U" ){type_U(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Np" ){type_Np(i) =1;}
   }
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pu" ){type_Pu(i) =1;}
   }
//  mat coord_a = zeros<mat>(sum(typeB) + 1,3);
//  mat coord_b = zeros<mat>(sum(typeB) + 1,3);
  mat coordH = zeros<mat>(coord.n_rows + 1,3);
//-----
  mat coord_H = zeros<mat>(sum(type_H ) + 1,3);
  mat coord_He = zeros<mat>(sum(type_He ) + 1,3);
  mat coord_Li = zeros<mat>(sum(type_Li ) + 1,3);
  mat coord_Be = zeros<mat>(sum(type_Be ) + 1,3);
  mat coord_B = zeros<mat>(sum(type_B ) + 1,3);
  mat coord_C = zeros<mat>(sum(type_C ) + 1,3);
  mat coord_N = zeros<mat>(sum(type_N ) + 1,3);
  mat coord_O = zeros<mat>(sum(type_O ) + 1,3);
  mat coord_F = zeros<mat>(sum(type_F ) + 1,3);
  mat coord_Ne = zeros<mat>(sum(type_Ne ) + 1,3);
  mat coord_Na = zeros<mat>(sum(type_Na ) + 1,3);
  mat coord_Mg = zeros<mat>(sum(type_Mg ) + 1,3);
  mat coord_Al = zeros<mat>(sum(type_Al ) + 1,3);
  mat coord_Si = zeros<mat>(sum(type_Si ) + 1,3);
  mat coord_P = zeros<mat>(sum(type_P ) + 1,3);
  mat coord_S = zeros<mat>(sum(type_S ) + 1,3);
  mat coord_Cl = zeros<mat>(sum(type_Cl ) + 1,3);
  mat coord_Ar = zeros<mat>(sum(type_Ar ) + 1,3);
  mat coord_K = zeros<mat>(sum(type_K ) + 1,3);
  mat coord_Ca = zeros<mat>(sum(type_Ca ) + 1,3);
  mat coord_Sc = zeros<mat>(sum(type_Sc ) + 1,3);
  mat coord_Ti = zeros<mat>(sum(type_Ti ) + 1,3);
  mat coord_V = zeros<mat>(sum(type_V ) + 1,3);
  mat coord_Cr = zeros<mat>(sum(type_Cr ) + 1,3);
  mat coord_Mn = zeros<mat>(sum(type_Mn ) + 1,3);
  mat coord_Fe = zeros<mat>(sum(type_Fe ) + 1,3);
  mat coord_Co = zeros<mat>(sum(type_Co ) + 1,3);
  mat coord_Ni = zeros<mat>(sum(type_Ni ) + 1,3);
  mat coord_Cu = zeros<mat>(sum(type_Cu ) + 1,3);
  mat coord_Zn = zeros<mat>(sum(type_Zn ) + 1,3);
  mat coord_Ga = zeros<mat>(sum(type_Ga ) + 1,3);
  mat coord_Ge = zeros<mat>(sum(type_Ge ) + 1,3);
  mat coord_As = zeros<mat>(sum(type_As ) + 1,3);
  mat coord_Se = zeros<mat>(sum(type_Se ) + 1,3);
  mat coord_Br = zeros<mat>(sum(type_Br ) + 1,3);
  mat coord_Kr = zeros<mat>(sum(type_Kr ) + 1,3);
  mat coord_Rb = zeros<mat>(sum(type_Rb ) + 1,3);
  mat coord_Sr = zeros<mat>(sum(type_Sr ) + 1,3);
  mat coord_Y = zeros<mat>(sum(type_Y ) + 1,3);
  mat coord_Zr = zeros<mat>(sum(type_Zr ) + 1,3);
  mat coord_Nb = zeros<mat>(sum(type_Nb ) + 1,3);
  mat coord_Mo = zeros<mat>(sum(type_Mo ) + 1,3);
  mat coord_Tc = zeros<mat>(sum(type_Tc ) + 1,3);
  mat coord_Ru = zeros<mat>(sum(type_Ru ) + 1,3);
  mat coord_Rh = zeros<mat>(sum(type_Rh ) + 1,3);
  mat coord_Pd = zeros<mat>(sum(type_Pd ) + 1,3);
  mat coord_Ag = zeros<mat>(sum(type_Ag ) + 1,3);
  mat coord_Cd = zeros<mat>(sum(type_Cd ) + 1,3);
  mat coord_In = zeros<mat>(sum(type_In ) + 1,3);
  mat coord_Sn = zeros<mat>(sum(type_Sn ) + 1,3);
  mat coord_Sb = zeros<mat>(sum(type_Sb ) + 1,3);
  mat coord_Te = zeros<mat>(sum(type_Te ) + 1,3);
  mat coord_I = zeros<mat>(sum(type_I ) + 1,3);
  mat coord_Xe = zeros<mat>(sum(type_Xe ) + 1,3);
  mat coord_Cs = zeros<mat>(sum(type_Cs ) + 1,3);
  mat coord_Ba = zeros<mat>(sum(type_Ba ) + 1,3);
  mat coord_La = zeros<mat>(sum(type_La ) + 1,3);
  mat coord_Ce = zeros<mat>(sum(type_Ce ) + 1,3);
  mat coord_Pr = zeros<mat>(sum(type_Pr ) + 1,3);
  mat coord_Nd = zeros<mat>(sum(type_Nd ) + 1,3);
  mat coord_Pm = zeros<mat>(sum(type_Pm ) + 1,3);
  mat coord_Sm = zeros<mat>(sum(type_Sm ) + 1,3);
  mat coord_Eu = zeros<mat>(sum(type_Eu ) + 1,3);
  mat coord_Gd = zeros<mat>(sum(type_Gd ) + 1,3);
  mat coord_Tb = zeros<mat>(sum(type_Tb ) + 1,3);
  mat coord_Dy = zeros<mat>(sum(type_Dy ) + 1,3);
  mat coord_Ho = zeros<mat>(sum(type_Ho ) + 1,3);
  mat coord_Er = zeros<mat>(sum(type_Er ) + 1,3);
  mat coord_Tm = zeros<mat>(sum(type_Tm ) + 1,3);
  mat coord_Yb = zeros<mat>(sum(type_Yb ) + 1,3);
  mat coord_Lu = zeros<mat>(sum(type_Lu ) + 1,3);
  mat coord_Hf = zeros<mat>(sum(type_Hf ) + 1,3);
  mat coord_Ta = zeros<mat>(sum(type_Ta ) + 1,3);
  mat coord_W = zeros<mat>(sum(type_W ) + 1,3);
  mat coord_Re = zeros<mat>(sum(type_Re ) + 1,3);
  mat coord_Os = zeros<mat>(sum(type_Os ) + 1,3);
  mat coord_Ir = zeros<mat>(sum(type_Ir ) + 1,3);
  mat coord_Pt = zeros<mat>(sum(type_Pt ) + 1,3);
  mat coord_Au = zeros<mat>(sum(type_Au ) + 1,3);
  mat coord_Hg = zeros<mat>(sum(type_Hg ) + 1,3);
  mat coord_Tl = zeros<mat>(sum(type_Tl ) + 1,3);
  mat coord_Pb = zeros<mat>(sum(type_Pb ) + 1,3);
  mat coord_Bi = zeros<mat>(sum(type_Bi ) + 1,3);
  mat coord_Po = zeros<mat>(sum(type_Po ) + 1,3);
  mat coord_At = zeros<mat>(sum(type_At ) + 1,3);
  mat coord_Rn = zeros<mat>(sum(type_Rn ) + 1,3);
  mat coord_Fr = zeros<mat>(sum(type_Fr ) + 1,3);
  mat coord_Ra = zeros<mat>(sum(type_Ra ) + 1,3);
  mat coord_Ac = zeros<mat>(sum(type_Ac ) + 1,3);
  mat coord_Th = zeros<mat>(sum(type_Th ) + 1,3);
  mat coord_Pa = zeros<mat>(sum(type_Pa ) + 1,3);
  mat coord_U = zeros<mat>(sum(type_U ) + 1,3);
  mat coord_Np = zeros<mat>(sum(type_Np ) + 1,3);
  mat coord_Pu = zeros<mat>(sum(type_Pu ) + 1,3);
//cout << "CAD" << endl;  
//-----
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "H" ){coord_H.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "He" ){coord_He.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Li" ){coord_Li.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Be" ){coord_Be.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "B" ){coord_B.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "C" ){coord_C.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "N" ){coord_N.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "O" ){coord_O.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "F" ){coord_F.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ne" ){coord_Ne.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Na" ){coord_Na.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mg" ){coord_Mg.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Al" ){coord_Al.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Si" ){coord_Si.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "P" ){coord_P.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "S" ){coord_S.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cl" ){coord_Cl.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ar" ){coord_Ar.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "K" ){coord_K.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ca" ){coord_Ca.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sc" ){coord_Sc.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ti" ){coord_Ti.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "V" ){coord_V.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cr" ){coord_Cr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mn" ){coord_Mn.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Fe" ){coord_Fe.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Co" ){coord_Co.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ni" ){coord_Ni.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cu" ){coord_Cu.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Zn" ){coord_Zn.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ga" ){coord_Ga.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ge" ){coord_Ge.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "As" ){coord_As.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Se" ){coord_Se.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Br" ){coord_Br.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Kr" ){coord_Kr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rb" ){coord_Rb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sr" ){coord_Sr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Y" ){coord_Y.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Zr" ){coord_Zr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Nb" ){coord_Nb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Mo" ){coord_Mo.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tc" ){coord_Tc.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ru" ){coord_Ru.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rh" ){coord_Rh.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pd" ){coord_Pd.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ag" ){coord_Ag.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cd" ){coord_Cd.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "In" ){coord_In.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sn" ){coord_Sn.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sb" ){coord_Sb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Te" ){coord_Te.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "I" ){coord_I.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Xe" ){coord_Xe.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Cs" ){coord_Cs.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ba" ){coord_Ba.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "La" ){coord_La.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ce" ){coord_Ce.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pr" ){coord_Pr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Nd" ){coord_Nd.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pm" ){coord_Pm.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Sm" ){coord_Sm.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Eu" ){coord_Eu.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Gd" ){coord_Gd.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tb" ){coord_Tb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Dy" ){coord_Dy.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ho" ){coord_Ho.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Er" ){coord_Er.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tm" ){coord_Tm.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Yb" ){coord_Yb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Lu" ){coord_Lu.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Hf" ){coord_Hf.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ta" ){coord_Ta.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "W" ){coord_W.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Re" ){coord_Re.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Os" ){coord_Os.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ir" ){coord_Ir.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pt" ){coord_Pt.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Au" ){coord_Au.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Hg" ){coord_Hg.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Tl" ){coord_Tl.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pb" ){coord_Pb.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Bi" ){coord_Bi.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Po" ){coord_Po.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "At" ){coord_At.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Rn" ){coord_Rn.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Fr" ){coord_Fr.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ra" ){coord_Ra.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Ac" ){coord_Ac.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Th" ){coord_Th.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pa" ){coord_Pa.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "U" ){coord_U.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Np" ){coord_Np.row(newJ) = coord.row(i); newJ++;}
   }
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     if(type[i] == "Pu" ){coord_Pu.row(newJ) = coord.row(i); newJ++;}
   }
//-----

//  newJ = 0;
//  for(int i=0; i < coord.n_rows; i++)  { 
//     if(type[i] == argv[2] ){coord_a.row(newJ) = coord.row(i); newJ++;}
//   }
//
//  newJ = 0;
//  for(int i=0; i < coord.n_rows; i++)  { 
//     if(type[i] == argv[3] ){coord_b.row(newJ) = coord.row(i); newJ++;}
//   }

//cout << "CAC" << endl;  
  newJ = 0;
  for(int i=0; i < coord.n_rows; i++)  { 
     coordH.row(newJ) = coord.row(i); newJ++;
   }

//cout << "CBC" << endl;  

// FOR For
for(int bigI = 0; bigI < hydrogenPos.n_rows; bigI++){

// Adding Hydrogen at the end

  coord_H.row(coord_H.n_rows - 1) = hydrogenPos.row(bigI);
  coord_He.row(coord_He.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Li.row(coord_Li.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Be.row(coord_Be.n_rows - 1) = hydrogenPos.row(bigI);
  coord_B.row(coord_B.n_rows - 1) = hydrogenPos.row(bigI);
  coord_C.row(coord_C.n_rows - 1) = hydrogenPos.row(bigI);
  coord_N.row(coord_N.n_rows - 1) = hydrogenPos.row(bigI);
  coord_O.row(coord_O.n_rows - 1) = hydrogenPos.row(bigI);
  coord_F.row(coord_F.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ne.row(coord_Ne.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Na.row(coord_Na.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Mg.row(coord_Mg.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Al.row(coord_Al.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Si.row(coord_Si.n_rows - 1) = hydrogenPos.row(bigI);
  coord_P.row(coord_P.n_rows - 1) = hydrogenPos.row(bigI);
  coord_S.row(coord_S.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Cl.row(coord_Cl.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ar.row(coord_Ar.n_rows - 1) = hydrogenPos.row(bigI);
  coord_K.row(coord_K.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ca.row(coord_Ca.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Sc.row(coord_Sc.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ti.row(coord_Ti.n_rows - 1) = hydrogenPos.row(bigI);
  coord_V.row(coord_V.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Cr.row(coord_Cr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Mn.row(coord_Mn.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Fe.row(coord_Fe.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Co.row(coord_Co.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ni.row(coord_Ni.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Cu.row(coord_Cu.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Zn.row(coord_Zn.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ga.row(coord_Ga.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ge.row(coord_Ge.n_rows - 1) = hydrogenPos.row(bigI);
  coord_As.row(coord_As.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Se.row(coord_Se.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Br.row(coord_Br.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Kr.row(coord_Kr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Rb.row(coord_Rb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Sr.row(coord_Sr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Y.row(coord_Y.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Zr.row(coord_Zr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Nb.row(coord_Nb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Mo.row(coord_Mo.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Tc.row(coord_Tc.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ru.row(coord_Ru.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Rh.row(coord_Rh.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pd.row(coord_Pd.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ag.row(coord_Ag.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Cd.row(coord_Cd.n_rows - 1) = hydrogenPos.row(bigI);
  coord_In.row(coord_In.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Sn.row(coord_Sn.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Sb.row(coord_Sb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Te.row(coord_Te.n_rows - 1) = hydrogenPos.row(bigI);
  coord_I.row(coord_I.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Xe.row(coord_Xe.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Cs.row(coord_Cs.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ba.row(coord_Ba.n_rows - 1) = hydrogenPos.row(bigI);
  coord_La.row(coord_La.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ce.row(coord_Ce.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pr.row(coord_Pr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Nd.row(coord_Nd.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pm.row(coord_Pm.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Sm.row(coord_Sm.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Eu.row(coord_Eu.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Gd.row(coord_Gd.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Tb.row(coord_Tb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Dy.row(coord_Dy.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ho.row(coord_Ho.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Er.row(coord_Er.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Tm.row(coord_Tm.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Yb.row(coord_Yb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Lu.row(coord_Lu.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Hf.row(coord_Hf.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ta.row(coord_Ta.n_rows - 1) = hydrogenPos.row(bigI);
  coord_W.row(coord_W.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Re.row(coord_Re.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Os.row(coord_Os.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ir.row(coord_Ir.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pt.row(coord_Pt.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Au.row(coord_Au.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Hg.row(coord_Hg.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Tl.row(coord_Tl.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pb.row(coord_Pb.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Bi.row(coord_Bi.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Po.row(coord_Po.n_rows - 1) = hydrogenPos.row(bigI);
  coord_At.row(coord_At.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Rn.row(coord_Rn.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Fr.row(coord_Fr.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ra.row(coord_Ra.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Ac.row(coord_Ac.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Th.row(coord_Th.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pa.row(coord_Pa.n_rows - 1) = hydrogenPos.row(bigI);
  coord_U.row(coord_U.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Np.row(coord_Np.n_rows - 1) = hydrogenPos.row(bigI);
  coord_Pu.row(coord_Pu.n_rows - 1) = hydrogenPos.row(bigI);

//  coord_a.row(coord_a.n_rows - 1) = hydrogenPos.row(bigI);
//  coord_b.row(coord_b.n_rows - 1) = hydrogenPos.row(bigI);
  coordH.row(coordH.n_rows - 1) = hydrogenPos.row(bigI);

//cout << "CCC" << endl;  
//coord_a.print("a");
//coord_b.print("b");
//coordH.print("H");

// Gaussian Smearing at Gauss-Lgedendre quadratue points
  rho_H = getGaussDistr(coord_H ,R,The,Phi,X,Y,Z,sig);
  rho_He = getGaussDistr(coord_He ,R,The,Phi,X,Y,Z,sig);
  rho_Li = getGaussDistr(coord_Li ,R,The,Phi,X,Y,Z,sig);
  rho_Be = getGaussDistr(coord_Be ,R,The,Phi,X,Y,Z,sig);
  rho_B = getGaussDistr(coord_B ,R,The,Phi,X,Y,Z,sig);
  rho_C = getGaussDistr(coord_C ,R,The,Phi,X,Y,Z,sig);
  rho_N = getGaussDistr(coord_N ,R,The,Phi,X,Y,Z,sig);
  rho_O = getGaussDistr(coord_O ,R,The,Phi,X,Y,Z,sig);
  rho_F = getGaussDistr(coord_F ,R,The,Phi,X,Y,Z,sig);
  rho_Ne = getGaussDistr(coord_Ne ,R,The,Phi,X,Y,Z,sig);
  rho_Na = getGaussDistr(coord_Na ,R,The,Phi,X,Y,Z,sig);
  rho_Mg = getGaussDistr(coord_Mg ,R,The,Phi,X,Y,Z,sig);
  rho_Al = getGaussDistr(coord_Al ,R,The,Phi,X,Y,Z,sig);
  rho_Si = getGaussDistr(coord_Si ,R,The,Phi,X,Y,Z,sig);
  rho_P = getGaussDistr(coord_P ,R,The,Phi,X,Y,Z,sig);
  rho_S = getGaussDistr(coord_S ,R,The,Phi,X,Y,Z,sig);
  rho_Cl = getGaussDistr(coord_Cl ,R,The,Phi,X,Y,Z,sig);
  rho_Ar = getGaussDistr(coord_Ar ,R,The,Phi,X,Y,Z,sig);
  rho_K = getGaussDistr(coord_K ,R,The,Phi,X,Y,Z,sig);
  rho_Ca = getGaussDistr(coord_Ca ,R,The,Phi,X,Y,Z,sig);
  rho_Sc = getGaussDistr(coord_Sc ,R,The,Phi,X,Y,Z,sig);
  rho_Ti = getGaussDistr(coord_Ti ,R,The,Phi,X,Y,Z,sig);
  rho_V = getGaussDistr(coord_V ,R,The,Phi,X,Y,Z,sig);
  rho_Cr = getGaussDistr(coord_Cr ,R,The,Phi,X,Y,Z,sig);
  rho_Mn = getGaussDistr(coord_Mn ,R,The,Phi,X,Y,Z,sig);
  rho_Fe = getGaussDistr(coord_Fe ,R,The,Phi,X,Y,Z,sig);
  rho_Co = getGaussDistr(coord_Co ,R,The,Phi,X,Y,Z,sig);
  rho_Ni = getGaussDistr(coord_Ni ,R,The,Phi,X,Y,Z,sig);
  rho_Cu = getGaussDistr(coord_Cu ,R,The,Phi,X,Y,Z,sig);
  rho_Zn = getGaussDistr(coord_Zn ,R,The,Phi,X,Y,Z,sig);
  rho_Ga = getGaussDistr(coord_Ga ,R,The,Phi,X,Y,Z,sig);
  rho_Ge = getGaussDistr(coord_Ge ,R,The,Phi,X,Y,Z,sig);
  rho_As = getGaussDistr(coord_As ,R,The,Phi,X,Y,Z,sig);
  rho_Se = getGaussDistr(coord_Se ,R,The,Phi,X,Y,Z,sig);
  rho_Br = getGaussDistr(coord_Br ,R,The,Phi,X,Y,Z,sig);
  rho_Kr = getGaussDistr(coord_Kr ,R,The,Phi,X,Y,Z,sig);
  rho_Rb = getGaussDistr(coord_Rb ,R,The,Phi,X,Y,Z,sig);
  rho_Sr = getGaussDistr(coord_Sr ,R,The,Phi,X,Y,Z,sig);
  rho_Y = getGaussDistr(coord_Y ,R,The,Phi,X,Y,Z,sig);
  rho_Zr = getGaussDistr(coord_Zr ,R,The,Phi,X,Y,Z,sig);
  rho_Nb = getGaussDistr(coord_Nb ,R,The,Phi,X,Y,Z,sig);
  rho_Mo = getGaussDistr(coord_Mo ,R,The,Phi,X,Y,Z,sig);
  rho_Tc = getGaussDistr(coord_Tc ,R,The,Phi,X,Y,Z,sig);
  rho_Ru = getGaussDistr(coord_Ru ,R,The,Phi,X,Y,Z,sig);
  rho_Rh = getGaussDistr(coord_Rh ,R,The,Phi,X,Y,Z,sig);
  rho_Pd = getGaussDistr(coord_Pd ,R,The,Phi,X,Y,Z,sig);
  rho_Ag = getGaussDistr(coord_Ag ,R,The,Phi,X,Y,Z,sig);
  rho_Cd = getGaussDistr(coord_Cd ,R,The,Phi,X,Y,Z,sig);
  rho_In = getGaussDistr(coord_In ,R,The,Phi,X,Y,Z,sig);
  rho_Sn = getGaussDistr(coord_Sn ,R,The,Phi,X,Y,Z,sig);
  rho_Sb = getGaussDistr(coord_Sb ,R,The,Phi,X,Y,Z,sig);
  rho_Te = getGaussDistr(coord_Te ,R,The,Phi,X,Y,Z,sig);
  rho_I = getGaussDistr(coord_I ,R,The,Phi,X,Y,Z,sig);
  rho_Xe = getGaussDistr(coord_Xe ,R,The,Phi,X,Y,Z,sig);
  rho_Cs = getGaussDistr(coord_Cs ,R,The,Phi,X,Y,Z,sig);
  rho_Ba = getGaussDistr(coord_Ba ,R,The,Phi,X,Y,Z,sig);
  rho_La = getGaussDistr(coord_La ,R,The,Phi,X,Y,Z,sig);
  rho_Ce = getGaussDistr(coord_Ce ,R,The,Phi,X,Y,Z,sig);
  rho_Pr = getGaussDistr(coord_Pr ,R,The,Phi,X,Y,Z,sig);
  rho_Nd = getGaussDistr(coord_Nd ,R,The,Phi,X,Y,Z,sig);
  rho_Pm = getGaussDistr(coord_Pm ,R,The,Phi,X,Y,Z,sig);
  rho_Sm = getGaussDistr(coord_Sm ,R,The,Phi,X,Y,Z,sig);
  rho_Eu = getGaussDistr(coord_Eu ,R,The,Phi,X,Y,Z,sig);
  rho_Gd = getGaussDistr(coord_Gd ,R,The,Phi,X,Y,Z,sig);
  rho_Tb = getGaussDistr(coord_Tb ,R,The,Phi,X,Y,Z,sig);
  rho_Dy = getGaussDistr(coord_Dy ,R,The,Phi,X,Y,Z,sig);
  rho_Ho = getGaussDistr(coord_Ho ,R,The,Phi,X,Y,Z,sig);
  rho_Er = getGaussDistr(coord_Er ,R,The,Phi,X,Y,Z,sig);
  rho_Tm = getGaussDistr(coord_Tm ,R,The,Phi,X,Y,Z,sig);
  rho_Yb = getGaussDistr(coord_Yb ,R,The,Phi,X,Y,Z,sig);
  rho_Lu = getGaussDistr(coord_Lu ,R,The,Phi,X,Y,Z,sig);
  rho_Hf = getGaussDistr(coord_Hf ,R,The,Phi,X,Y,Z,sig);
  rho_Ta = getGaussDistr(coord_Ta ,R,The,Phi,X,Y,Z,sig);
  rho_W = getGaussDistr(coord_W ,R,The,Phi,X,Y,Z,sig);
  rho_Re = getGaussDistr(coord_Re ,R,The,Phi,X,Y,Z,sig);
  rho_Os = getGaussDistr(coord_Os ,R,The,Phi,X,Y,Z,sig);
  rho_Ir = getGaussDistr(coord_Ir ,R,The,Phi,X,Y,Z,sig);
  rho_Pt = getGaussDistr(coord_Pt ,R,The,Phi,X,Y,Z,sig);
  rho_Au = getGaussDistr(coord_Au ,R,The,Phi,X,Y,Z,sig);
  rho_Hg = getGaussDistr(coord_Hg ,R,The,Phi,X,Y,Z,sig);
  rho_Tl = getGaussDistr(coord_Tl ,R,The,Phi,X,Y,Z,sig);
  rho_Pb = getGaussDistr(coord_Pb ,R,The,Phi,X,Y,Z,sig);
  rho_Bi = getGaussDistr(coord_Bi ,R,The,Phi,X,Y,Z,sig);
  rho_Po = getGaussDistr(coord_Po ,R,The,Phi,X,Y,Z,sig);
  rho_At = getGaussDistr(coord_At ,R,The,Phi,X,Y,Z,sig);
  rho_Rn = getGaussDistr(coord_Rn ,R,The,Phi,X,Y,Z,sig);
  rho_Fr = getGaussDistr(coord_Fr ,R,The,Phi,X,Y,Z,sig);
  rho_Ra = getGaussDistr(coord_Ra ,R,The,Phi,X,Y,Z,sig);
  rho_Ac = getGaussDistr(coord_Ac ,R,The,Phi,X,Y,Z,sig);
  rho_Th = getGaussDistr(coord_Th ,R,The,Phi,X,Y,Z,sig);
  rho_Pa = getGaussDistr(coord_Pa ,R,The,Phi,X,Y,Z,sig);
  rho_U = getGaussDistr(coord_U ,R,The,Phi,X,Y,Z,sig);
  rho_Np = getGaussDistr(coord_Np ,R,The,Phi,X,Y,Z,sig);
  rho_Pu = getGaussDistr(coord_Pu ,R,The,Phi,X,Y,Z,sig);

//  rho_a = getGaussDistr(coord_a,R, The, Phi, X, Y, Z, sig);
//  rho_b = getGaussDistr(coord_b,R, The, Phi, X, Y, Z, sig);
  rhoAll = getGaussDistr(coordH,R, The, Phi, X, Y, Z, sig);
//cout << "CDC" << endl;  

  //----------------------------------------------------------------------------------------------------------------
  // Part 3) get coefs c_anlm by Integration ->  where C(a,n,I) where a is the type, n is the radial basis function
  //         and I is (l,m) sequenially. intMa and intMb multiplied with getT are integrands.
  //----------------------------------------------------------------------------------------------------------------

  intMeH = rho_H %GLC;
  intMeHe = rho_He %GLC;
  intMeLi = rho_Li %GLC;
  intMeBe = rho_Be %GLC;
  intMeB = rho_B %GLC;
  intMeC = rho_C %GLC;
  intMeN = rho_N %GLC;
  intMeO = rho_O %GLC;
  intMeF = rho_F %GLC;
  intMeNe = rho_Ne %GLC;
  intMeNa = rho_Na %GLC;
  intMeMg = rho_Mg %GLC;
  intMeAl = rho_Al %GLC;
  intMeSi = rho_Si %GLC;
  intMeP = rho_P %GLC;
  intMeS = rho_S %GLC;
  intMeCl = rho_Cl %GLC;
  intMeAr = rho_Ar %GLC;
  intMeK = rho_K %GLC;
  intMeCa = rho_Ca %GLC;
  intMeSc = rho_Sc %GLC;
  intMeTi = rho_Ti %GLC;
  intMeV = rho_V %GLC;
  intMeCr = rho_Cr %GLC;
  intMeMn = rho_Mn %GLC;
  intMeFe = rho_Fe %GLC;
  intMeCo = rho_Co %GLC;
  intMeNi = rho_Ni %GLC;
  intMeCu = rho_Cu %GLC;
  intMeZn = rho_Zn %GLC;
  intMeGa = rho_Ga %GLC;
  intMeGe = rho_Ge %GLC;
  intMeAs = rho_As %GLC;
  intMeSe = rho_Se %GLC;
  intMeBr = rho_Br %GLC;
  intMeKr = rho_Kr %GLC;
  intMeRb = rho_Rb %GLC;
  intMeSr = rho_Sr %GLC;
  intMeY = rho_Y %GLC;
  intMeZr = rho_Zr %GLC;
  intMeNb = rho_Nb %GLC;
  intMeMo = rho_Mo %GLC;
  intMeTc = rho_Tc %GLC;
  intMeRu = rho_Ru %GLC;
  intMeRh = rho_Rh %GLC;
  intMePd = rho_Pd %GLC;
  intMeAg = rho_Ag %GLC;
  intMeCd = rho_Cd %GLC;
  intMeIn = rho_In %GLC;
  intMeSn = rho_Sn %GLC;
  intMeSb = rho_Sb %GLC;
  intMeTe = rho_Te %GLC;
  intMeI = rho_I %GLC;
  intMeXe = rho_Xe %GLC;
  intMeCs = rho_Cs %GLC;
  intMeBa = rho_Ba %GLC;
  intMeLa = rho_La %GLC;
  intMeCe = rho_Ce %GLC;
  intMePr = rho_Pr %GLC;
  intMeNd = rho_Nd %GLC;
  intMePm = rho_Pm %GLC;
  intMeSm = rho_Sm %GLC;
  intMeEu = rho_Eu %GLC;
  intMeGd = rho_Gd %GLC;
  intMeTb = rho_Tb %GLC;
  intMeDy = rho_Dy %GLC;
  intMeHo = rho_Ho %GLC;
  intMeEr = rho_Er %GLC;
  intMeTm = rho_Tm %GLC;
  intMeYb = rho_Yb %GLC;
  intMeLu = rho_Lu %GLC;
  intMeHf = rho_Hf %GLC;
  intMeTa = rho_Ta %GLC;
  intMeW = rho_W %GLC;
  intMeRe = rho_Re %GLC;
  intMeOs = rho_Os %GLC;
  intMeIr = rho_Ir %GLC;
  intMePt = rho_Pt %GLC;
  intMeAu = rho_Au %GLC;
  intMeHg = rho_Hg %GLC;
  intMeTl = rho_Tl %GLC;
  intMePb = rho_Pb %GLC;
  intMeBi = rho_Bi %GLC;
  intMePo = rho_Po %GLC;
  intMeAt = rho_At %GLC;
  intMeRn = rho_Rn %GLC;
  intMeFr = rho_Fr %GLC;
  intMeRa = rho_Ra %GLC;
  intMeAc = rho_Ac %GLC;
  intMeTh = rho_Th %GLC;
  intMePa = rho_Pa %GLC;
  intMeU = rho_U %GLC;
  intMeNp = rho_Np %GLC;
  intMePu = rho_Pu %GLC;

//  intMea = rho_a%GLC;
//  intMeb = rho_b%GLC;
  intAll = rhoAll%GLC;

  globalI = 0;
//cout << "CEC" << endl;  

  for(int n=0; n < radialN; n++) {
//    Cbuff=rsc*integ3Dvec(intMea,intMeb,intAll, getTMat(n,0,0,g,Y0,R,The,Phi));
    Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getTMat(n,0,0,g,Y0,R,The,Phi));
    for(int i=0; i < typesN; i++){
       C(i,n,0) = Cbuff(i);
    }
  }
//cout << "CFC" << endl;  

  for(int l = 1; l <= lMax; l++){
    for(int m=-l; m <= l; m++) {
      globalI++;
      for(int n=0; n < radialN; n++) {
          if(l == 1){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
              intAll, getT(n,l,m,g,Y1,R,The,Phi));}
          else if(l == 2){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y2,R,The,Phi));}
          else if(l == 3){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y3,R,The,Phi));}
          else if(l == 4){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y4,R,The,Phi));}
          else if(l == 5){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y5,R,The,Phi));}
          else if(l == 6){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y6,R,The,Phi));}
          else if(l == 7){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y7,R,The,Phi));}
          else if(l == 8){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y8,R,The,Phi));}
          else if(l == 9){Cbuff=rsc*integFull(
  intMeH,
  intMeHe,
  intMeLi,
  intMeBe,
  intMeB,
  intMeC,
  intMeN,
  intMeO,
  intMeF,
  intMeNe,
  intMeNa,
  intMeMg,
  intMeAl,
  intMeSi,
  intMeP,
  intMeS,
  intMeCl,
  intMeAr,
  intMeK,
  intMeCa,
  intMeSc,
  intMeTi,
  intMeV,
  intMeCr,
  intMeMn,
  intMeFe,
  intMeCo,
  intMeNi,
  intMeCu,
  intMeZn,
  intMeGa,
  intMeGe,
  intMeAs,
  intMeSe,
  intMeBr,
  intMeKr,
  intMeRb,
  intMeSr,
  intMeY,
  intMeZr,
  intMeNb,
  intMeMo,
  intMeTc,
  intMeRu,
  intMeRh,
  intMePd,
  intMeAg,
  intMeCd,
  intMeIn,
  intMeSn,
  intMeSb,
  intMeTe,
  intMeI,
  intMeXe,
  intMeCs,
  intMeBa,
  intMeLa,
  intMeCe,
  intMePr,
  intMeNd,
  intMePm,
  intMeSm,
  intMeEu,
  intMeGd,
  intMeTb,
  intMeDy,
  intMeHo,
  intMeEr,
  intMeTm,
  intMeYb,
  intMeLu,
  intMeHf,
  intMeTa,
  intMeW,
  intMeRe,
  intMeOs,
  intMeIr,
  intMePt,
  intMeAu,
  intMeHg,
  intMeTl,
  intMePb,
  intMeBi,
  intMePo,
  intMeAt,
  intMeRn,
  intMeFr,
  intMeRa,
  intMeAc,
  intMeTh,
  intMePa,
  intMeU,
  intMeNp,
  intMePu,
  intAll, getT(n,l,m,g,Y9,R,The,Phi));}
          else{cout << "ERROR, l too large" << endl;}
        for(int i=0; i<typesN;i++){
        C(i,n,globalI) = Cbuff(i);
        }
      }
    }
  }

//cout << "CAD" << endl;  
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
//            cout << a  << " " <<n1 << " "  << n2 << " " << l << " " << P[a][n1][n2][l] << endl;
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
