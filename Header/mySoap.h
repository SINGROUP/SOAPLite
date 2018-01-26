#ifndef MYSOAP/* Include guard */
#define SOAP

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>
#include <iomanip>
#include "myArmadillo.h"


using namespace std;
using namespace arma;
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getY(int l, vec Theta, vec Phi){
 cube Yn = zeros<cube>(2*l + 1,Theta.n_elem, Phi.n_elem); 
    for(int m=-l; m <= l; m++){ 
     for(int t=0; t < Theta.n_elem; t++){ 
       for(int p=0; p < Phi.n_elem; p++){ 
        Yn.at(l+m,t,p) = tesseral_spherical_harm(l,m,Theta.at(t),Phi.at(p)); 
        }
      }
    }
    return Yn;
};
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getT(int n, int l,int m, mat g, cube Yl, vec R, vec Theta, vec Phi){
 cube T =  zeros<cube>(R.n_elem,Theta.n_elem,Phi.n_elem); // (r, Theta, Phi)
 vec RR = R%R;
 vec sinT = sin(Theta);
  for(int r=0; r < R.n_elem; r++){ 
    for(int t=0; t < Theta.n_elem; t++){ 
      for(int p=0; p < Phi.n_elem; p++){ 
        T.at(r,t,p) = g.at(n,r)*RR.at(r)*sinT.at(t)*Yl.at(l+m,t,p);  // R*R*sin(Theta) is the Jacobian for the integration.
      }
    }
  }
  return T;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
//cube getT(int n, int l,int m, mat g, cube Yl, vec R, vec Theta, vec Phi){
// cube T =  zeros<cube>(R.n_elem,Theta.n_elem,Phi.n_elem); // (r, Theta, Phi)
//  for(int r=0; r < R.n_elem; r++){ 
//    for(int t=0; t < Theta.n_elem; t++){ 
//      for(int p=0; p < Phi.n_elem; p++){ 
//        T.at(r,t,p) = g.at(n,r)*R.at(r)* R.at(r)*sin(Theta.at(t))*Yl.at(l+m,t,p);  // R*R*sin(Theta) is the Jacobian for the integration.
//      }
//    }
//  }
//
//  return T;
//}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getTMat(int n, int l,int m, mat g, mat Y00, vec R, vec Theta, vec Phi){
 cube T =  zeros<cube>(R.n_elem,Theta.n_elem,Phi.n_elem); // (r, Theta, Phi)

  for(int r=0; r < R.n_elem; r++){ 
    for(int t=0; t < Theta.n_elem; t++){ 
      for(int p=0; p < Phi.n_elem; p++){ 
        T.at(r,t,p) = R.at(r)* R.at(r)*sin(Theta.at(t))*g.at(n,r)*Y00.at(t,p); 
      }
    }
  }
  return T;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getGaussDistr(mat coord,vec R, vec The, vec Phi, cube X, cube Y , cube Z, double sig){
//!!!! DOUNT FORGET IN MAIN TO DO:
//  sig = 1/sig;
//  sig = sig*sig;
//  sig = 0.5*sig;
//  DONT FORGET TO RESCALE ORIGIN FOR GAUSS-LEGENDRE QUADUATURE!
  
  cube G = zeros<cube>(R.n_elem,The.n_elem,Phi.n_elem);
  double x1, x2, x3; // X1 = R*sin(Theta)*cos(Phi), X2 = R*sin(Theta)*sin(Phi), X3 = R*cos(Theta) 
  rowvec origin = coord.row(coord.n_rows - 1); // last xyz position of the .xyz file.

  for(int i=0; i < coord.n_rows; i++){
  coord.row(i) = coord.row(i) - origin; // shifting all xyz position so that H is on the origin. 
  }
if(coord.n_rows > 0){
  for(int p=0; p < coord.n_rows - 1; p++)
    for(int i=0; i < R.n_rows; i++)
      for(int j=0; j < The.n_rows; j++)
        for(int k=0; k < Phi.n_rows; k++){ 
          x1 = X.at(i,j,k);
          x2 = Y.at(i,j,k);
          x3 = Z.at(i,j,k);
//          G(i,j,k) = G(i,j,k) + exp(-abs((coord(p,0) - x1)) - abs((coord(p,1) - x2)) - abs((coord(p,2) - x3)));
          G.at(i,j,k) = G.at(i,j,k) + exp(-pow(((coord(p,0) - x1)),2) - pow((coord(p,1) - x2),2) - pow((coord(p,2) - x3),2));
         }
 }

 return G;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double integ3D(cube intMea,cube Tnlm){

  // DANGER!!! NOT RESCALED -> MUST BE RESCALED BY 0.5*0.5*0.5*pi*pi*rcut*integ3D() in main.
  // This is to increase the computation speed.
  double z;

  mat x1 = sum(intMea%Tnlm);
  rowvec y1 = sum(x1);

  z = sum(y1);

return z;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec integ3Dvec(cube intMea,cube intMeb, cube intAll ,cube Tnlm){

  // DANGER!!! NOT RESCALED -> MUST BE RESCALED BY 0.5*0.5*0.5*pi*pi*rcut*integ3D() in main.
  // This is to increase the computation time.
  vec z(3);

  mat x1 = sum(intMea%Tnlm);
  mat x2 = sum(intMeb%Tnlm);
  mat x3 = sum(intAll%Tnlm);

  rowvec y1 = sum(x1);
  rowvec y2 = sum(x2);
  rowvec y3 = sum(x3);

  z(0) = sum(y1);
  z(1) = sum(y2);
  z(2) = sum(y3);
//z.print("z");
return z;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec integ3Dvec3Atom(cube intMea,cube intMeb, cube intMec,  cube intAll ,cube Tnlm){

  // DANGER!!! NOT RESCALED -> MUST BE RESCALED BY 0.5*0.5*0.5*pi*pi*rcut*integ3D() in main.
  // This is to increase the computation time.
  vec z(4);

  mat x1 = sum(intMea%Tnlm);
  mat x2 = sum(intMeb%Tnlm);
  mat x3 = sum(intMec%Tnlm);
  mat x4 = sum(intAll%Tnlm);

  rowvec y1 = sum(x1);
  rowvec y2 = sum(x2);
  rowvec y3 = sum(x3);
  rowvec y4 = sum(x4);

  z(0) = sum(y1);
  z(1) = sum(y2);
  z(2) = sum(y3);
  z(3) = sum(y4);
//z.print("z");
return z;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec integFull(
    cube intMe_1,
    cube intMe_2,
    cube intMe_3,
    cube intMe_4,
    cube intMe_5,
    cube intMe_6,
    cube intMe_7,
    cube intMe_8,
    cube intMe_9,
    cube intMe_10,
    cube intMe_11,
    cube intMe_12,
    cube intMe_13,
    cube intMe_14,
    cube intMe_15,
    cube intMe_16,
    cube intMe_17,
    cube intMe_18,
    cube intMe_19,
    cube intMe_20,
    cube intMe_21,
    cube intMe_22,
    cube intMe_23,
    cube intMe_24,
    cube intMe_25,
    cube intMe_26,
    cube intMe_27,
    cube intMe_28,
    cube intMe_29,
    cube intMe_30,
    cube intMe_31,
    cube intMe_32,
    cube intMe_33,
    cube intMe_34,
    cube intMe_35,
    cube intMe_36,
    cube intMe_37,
    cube intMe_38,
    cube intMe_39,
    cube intMe_40,
    cube intMe_41,
    cube intMe_42,
    cube intMe_43,
    cube intMe_44,
    cube intMe_45,
    cube intMe_46,
    cube intMe_47,
    cube intMe_48,
    cube intMe_49,
    cube intMe_50,
    cube intMe_51,
    cube intMe_52,
    cube intMe_53,
    cube intMe_54,
    cube intMe_55,
    cube intMe_56,
    cube intMe_57,
    cube intMe_58,
    cube intMe_59,
    cube intMe_60,
    cube intMe_61,
    cube intMe_62,
    cube intMe_63,
    cube intMe_64,
    cube intMe_65,
    cube intMe_66,
    cube intMe_67,
    cube intMe_68,
    cube intMe_69,
    cube intMe_70,
    cube intMe_71,
    cube intMe_72,
    cube intMe_73,
    cube intMe_74,
    cube intMe_75,
    cube intMe_76,
    cube intMe_77,
    cube intMe_78,
    cube intMe_79,
    cube intMe_80,
    cube intMe_81,
    cube intMe_82,
    cube intMe_83,
    cube intMe_84,
    cube intMe_85,
    cube intMe_86,
    cube intMe_87,
    cube intMe_88,
    cube intMe_89,
    cube intMe_90,
    cube intMe_91,
    cube intMe_92,
    cube intMe_93,
    cube intMe_94,
    cube intAll ,cube Tnlm){

  // DANGER!!! NOT RESCALED -> MUST BE RESCALED BY 0.5*0.5*0.5*pi*pi*rcut*integ3D() in main.
  // This is to increase the computation time.
  vec z(95);
  //cout << "DAA" << endl;
 mat x1 = sum(intMe_1%Tnlm);
 mat x2 = sum(intMe_2%Tnlm);
 mat x3 = sum(intMe_3%Tnlm);
 mat x4 = sum(intMe_4%Tnlm);
 mat x5 = sum(intMe_5%Tnlm);
 mat x6 = sum(intMe_6%Tnlm);
 mat x7 = sum(intMe_7%Tnlm);
 mat x8 = sum(intMe_8%Tnlm);
 mat x9 = sum(intMe_9%Tnlm);
 mat x10 = sum(intMe_10%Tnlm);
 mat x11 = sum(intMe_11%Tnlm);
 mat x12 = sum(intMe_12%Tnlm);
 mat x13 = sum(intMe_13%Tnlm);
 mat x14 = sum(intMe_14%Tnlm);
 mat x15 = sum(intMe_15%Tnlm);
 mat x16 = sum(intMe_16%Tnlm);
 mat x17 = sum(intMe_17%Tnlm);
 mat x18 = sum(intMe_18%Tnlm);
 mat x19 = sum(intMe_19%Tnlm);
 mat x20 = sum(intMe_20%Tnlm);
 mat x21 = sum(intMe_21%Tnlm);
 mat x22 = sum(intMe_22%Tnlm);
 mat x23 = sum(intMe_23%Tnlm);
 mat x24 = sum(intMe_24%Tnlm);
 mat x25 = sum(intMe_25%Tnlm);
 mat x26 = sum(intMe_26%Tnlm);
 mat x27 = sum(intMe_27%Tnlm);
 mat x28 = sum(intMe_28%Tnlm);
 mat x29 = sum(intMe_29%Tnlm);
 mat x30 = sum(intMe_30%Tnlm);
 mat x31 = sum(intMe_31%Tnlm);
 mat x32 = sum(intMe_32%Tnlm);
 mat x33 = sum(intMe_33%Tnlm);
 mat x34 = sum(intMe_34%Tnlm);
 mat x35 = sum(intMe_35%Tnlm);
 mat x36 = sum(intMe_36%Tnlm);
 mat x37 = sum(intMe_37%Tnlm);
 mat x38 = sum(intMe_38%Tnlm);
 mat x39 = sum(intMe_39%Tnlm);
 mat x40 = sum(intMe_40%Tnlm);
 mat x41 = sum(intMe_41%Tnlm);
 mat x42 = sum(intMe_42%Tnlm);
 mat x43 = sum(intMe_43%Tnlm);
 mat x44 = sum(intMe_44%Tnlm);
 mat x45 = sum(intMe_45%Tnlm);
 mat x46 = sum(intMe_46%Tnlm);
 mat x47 = sum(intMe_47%Tnlm);
 mat x48 = sum(intMe_48%Tnlm);
 mat x49 = sum(intMe_49%Tnlm);
 mat x50 = sum(intMe_50%Tnlm);
 mat x51 = sum(intMe_51%Tnlm);
 mat x52 = sum(intMe_52%Tnlm);
 mat x53 = sum(intMe_53%Tnlm);
 mat x54 = sum(intMe_54%Tnlm);
 mat x55 = sum(intMe_55%Tnlm);
 mat x56 = sum(intMe_56%Tnlm);
 mat x57 = sum(intMe_57%Tnlm);
 mat x58 = sum(intMe_58%Tnlm);
 mat x59 = sum(intMe_59%Tnlm);
 mat x60 = sum(intMe_60%Tnlm);
 mat x61 = sum(intMe_61%Tnlm);
 mat x62 = sum(intMe_62%Tnlm);
 mat x63 = sum(intMe_63%Tnlm);
 mat x64 = sum(intMe_64%Tnlm);
 mat x65 = sum(intMe_65%Tnlm);
 mat x66 = sum(intMe_66%Tnlm);
 mat x67 = sum(intMe_67%Tnlm);
 mat x68 = sum(intMe_68%Tnlm);
 mat x69 = sum(intMe_69%Tnlm);
 mat x70 = sum(intMe_70%Tnlm);
 mat x71 = sum(intMe_71%Tnlm);
 mat x72 = sum(intMe_72%Tnlm);
 mat x73 = sum(intMe_73%Tnlm);
 mat x74 = sum(intMe_74%Tnlm);
 mat x75 = sum(intMe_75%Tnlm);
 mat x76 = sum(intMe_76%Tnlm);
 mat x77 = sum(intMe_77%Tnlm);
 mat x78 = sum(intMe_78%Tnlm);
 mat x79 = sum(intMe_79%Tnlm);
 mat x80 = sum(intMe_80%Tnlm);
 mat x81 = sum(intMe_81%Tnlm);
 mat x82 = sum(intMe_82%Tnlm);
 mat x83 = sum(intMe_83%Tnlm);
 mat x84 = sum(intMe_84%Tnlm);
 mat x85 = sum(intMe_85%Tnlm);
 mat x86 = sum(intMe_86%Tnlm);
 mat x87 = sum(intMe_87%Tnlm);
 mat x88 = sum(intMe_88%Tnlm);
 mat x89 = sum(intMe_89%Tnlm);
 mat x90 = sum(intMe_90%Tnlm);
 mat x91 = sum(intMe_91%Tnlm);
 mat x92 = sum(intMe_92%Tnlm);
 mat x93 = sum(intMe_93%Tnlm);
 mat x94 = sum(intMe_94%Tnlm);
 mat x95 = sum(intAll%Tnlm);

  //cout << "DBB" << endl;
  //cout << "DAB" << endl;
  //

rowvec y1 = sum(x1);
rowvec y2 = sum(x2);
rowvec y3 = sum(x3);
rowvec y4 = sum(x4);
rowvec y5 = sum(x5);
rowvec y6 = sum(x6);
rowvec y7 = sum(x7);
rowvec y8 = sum(x8);
rowvec y9 = sum(x9);
rowvec y10 = sum(x10);
rowvec y11 = sum(x11);
rowvec y12 = sum(x12);
rowvec y13 = sum(x13);
rowvec y14 = sum(x14);
rowvec y15 = sum(x15);
rowvec y16 = sum(x16);
rowvec y17 = sum(x17);
rowvec y18 = sum(x18);
rowvec y19 = sum(x19);
rowvec y20 = sum(x20);
rowvec y21 = sum(x21);
rowvec y22 = sum(x22);
rowvec y23 = sum(x23);
rowvec y24 = sum(x24);
rowvec y25 = sum(x25);
rowvec y26 = sum(x26);
rowvec y27 = sum(x27);
rowvec y28 = sum(x28);
rowvec y29 = sum(x29);
rowvec y30 = sum(x30);
rowvec y31 = sum(x31);
rowvec y32 = sum(x32);
rowvec y33 = sum(x33);
rowvec y34 = sum(x34);
rowvec y35 = sum(x35);
rowvec y36 = sum(x36);
rowvec y37 = sum(x37);
rowvec y38 = sum(x38);
rowvec y39 = sum(x39);
rowvec y40 = sum(x40);
rowvec y41 = sum(x41);
rowvec y42 = sum(x42);
rowvec y43 = sum(x43);
rowvec y44 = sum(x44);
rowvec y45 = sum(x45);
rowvec y46 = sum(x46);
rowvec y47 = sum(x47);
rowvec y48 = sum(x48);
rowvec y49 = sum(x49);
rowvec y50 = sum(x50);
rowvec y51 = sum(x51);
rowvec y52 = sum(x52);
rowvec y53 = sum(x53);
rowvec y54 = sum(x54);
rowvec y55 = sum(x55);
rowvec y56 = sum(x56);
rowvec y57 = sum(x57);
rowvec y58 = sum(x58);
rowvec y59 = sum(x59);
rowvec y60 = sum(x60);
rowvec y61 = sum(x61);
rowvec y62 = sum(x62);
rowvec y63 = sum(x63);
rowvec y64 = sum(x64);
rowvec y65 = sum(x65);
rowvec y66 = sum(x66);
rowvec y67 = sum(x67);
rowvec y68 = sum(x68);
rowvec y69 = sum(x69);
rowvec y70 = sum(x70);
rowvec y71 = sum(x71);
rowvec y72 = sum(x72);
rowvec y73 = sum(x73);
rowvec y74 = sum(x74);
rowvec y75 = sum(x75);
rowvec y76 = sum(x76);
rowvec y77 = sum(x77);
rowvec y78 = sum(x78);
rowvec y79 = sum(x79);
rowvec y80 = sum(x80);
rowvec y81 = sum(x81);
rowvec y82 = sum(x82);
rowvec y83 = sum(x83);
rowvec y84 = sum(x84);
rowvec y85 = sum(x85);
rowvec y86 = sum(x86);
rowvec y87 = sum(x87);
rowvec y88 = sum(x88);
rowvec y89 = sum(x89);
rowvec y90 = sum(x90);
rowvec y91 = sum(x91);
rowvec y92 = sum(x92);
rowvec y93 = sum(x93);
rowvec y94 = sum(x94);

rowvec y95 = sum(x95);
  //cout << "DAC" << endl;

 z(0) = sum(y1);
 z(1) = sum(y2);
 z(2) = sum(y3);
 z(3) = sum(y4);
 z(4) = sum(y5);
 z(5) = sum(y6);
 z(6) = sum(y7);
 z(7) = sum(y8);
 z(8) = sum(y9);
 z(9) = sum(y10);
 z(10) = sum(y11);
 z(11) = sum(y12);
 z(12) = sum(y13);
 z(13) = sum(y14);
 z(14) = sum(y15);
 z(15) = sum(y16);
 z(16) = sum(y17);
 z(17) = sum(y18);
 z(18) = sum(y19);
 z(19) = sum(y20);
 z(20) = sum(y21);
 z(21) = sum(y22);
 z(22) = sum(y23);
 z(23) = sum(y24);
 z(24) = sum(y25);
 z(25) = sum(y26);
 z(26) = sum(y27);
 z(27) = sum(y28);
 z(28) = sum(y29);
 z(29) = sum(y30);
 z(30) = sum(y31);
 z(31) = sum(y32);
 z(32) = sum(y33);
 z(33) = sum(y34);
 z(34) = sum(y35);
 z(35) = sum(y36);
 z(36) = sum(y37);
 z(37) = sum(y38);
 z(38) = sum(y39);
 z(39) = sum(y40);
 z(40) = sum(y41);
 z(41) = sum(y42);
 z(42) = sum(y43);
 z(43) = sum(y44);
 z(44) = sum(y45);
 z(45) = sum(y46);
 z(46) = sum(y47);
 z(47) = sum(y48);
 z(48) = sum(y49);
 z(49) = sum(y50);
 z(50) = sum(y51);
 z(51) = sum(y52);
 z(52) = sum(y53);
 z(53) = sum(y54);
 z(54) = sum(y55);
 z(55) = sum(y56);
 z(56) = sum(y57);
 z(57) = sum(y58);
 z(58) = sum(y59);
 z(59) = sum(y60);
 z(60) = sum(y61);
 z(61) = sum(y62);
 z(62) = sum(y63);
 z(63) = sum(y64);
 z(64) = sum(y65);
 z(65) = sum(y66);
 z(66) = sum(y67);
 z(67) = sum(y68);
 z(68) = sum(y69);
 z(69) = sum(y70);
 z(70) = sum(y71);
 z(71) = sum(y72);
 z(72) = sum(y73);
 z(73) = sum(y74);
 z(74) = sum(y75);
 z(75) = sum(y76);
 z(76) = sum(y77);
 z(77) = sum(y78);
 z(78) = sum(y79);
 z(79) = sum(y80);
 z(80) = sum(y81);
 z(81) = sum(y82);
 z(82) = sum(y83);
 z(83) = sum(y84);
 z(84) = sum(y85);
 z(85) = sum(y86);
 z(86) = sum(y87);
 z(87) = sum(y88);
 z(88) = sum(y89);
 z(89) = sum(y90);
 z(90) = sum(y91);
 z(91) = sum(y92);
 z(92) = sum(y93);
 z(93) = sum(y94);
 z(94) = sum(y95);
//z.print("z");
return z;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------

#endif 

