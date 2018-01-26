#ifndef MYMATH/* Include guard */
#define MYMATH

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>
#include <iomanip>
#include "myArmadillo.h"
#include "fileoperation.h"


using namespace std;
using namespace arma;
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double factorialR(int b){
  if(b<0) {cout << "ERROR: Negative Input..." << endl; exit(1);}
  if(b>20) {cout << "ERROR: Number Too Large..." << endl; exit(1);}
  int c = b;
long  int  a = 1;
  if(b==0) return 1; 

for(int i=1; i <= b; i++)
    { 
      a *= c;
      c = c - 1;
       }

return (double) a;

};
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double legendre_poly(int l, int m, double x){

  double fact,pll,pmm,pmmp1,somx2;
  int ll;

  if (m < 0 || m > l || fabs(x) > 1.0){ cout << "ERROR: Bad arguments in routine legendre_poly" << endl; exit(1);}

  pmm = 1.0;

  if(m > 0) {
    somx2=sqrt((1.0 - x)*(1.0 + x));
    fact=1.0;
    for(int i=1; i <= m; i++)
        { 
          pmm *= -fact*somx2;
          fact += 2.0;
           }
   }

  if(l == m) return pmm;

  else{
    pmmp1 = x*(2*m+1)*pmm;
    if(l==(m+1)) return pmmp1;
    else{
      for(ll=m+2; ll<=l; ll++){
        pll=(x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/ (double) (ll-m);
        pmm = pmmp1; 
        pmmp1= pll; 
      
      }
      return pll;
    }
  
  }

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseral_spherical_harm(int l, int m, double theta, double phi){
if(m==0) return sqrt((2.0*l +1)/(12.5663706144))*legendre_poly(l,m,cos(theta));
if(m > 0) return 1.414213562373095 * sqrt((2.0*l +1)/(12.5663706144)*factorialR(l-m)/factorialR(l+m))*legendre_poly(l,m,cos(theta))*cos(m*phi);
if(m < 0){int mfabs = fabs(m);
 return 1.414213562373095 * sqrt((2.0*l +1)/(12.5663706144)*factorialR(l-mfabs)/factorialR(l+mfabs)) *legendre_poly(l,mfabs,cos(theta))*sin(mfabs*phi);}
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseral_spherical_harmNoN(int l, int m, double theta, double phi){
if(m==0) return legendre_poly(l,m,cos(theta));
if(m > 0) return legendre_poly(l,m,cos(theta))*cos(m*phi);
if(m < 0){int mfabs = fabs(m);
 return legendre_poly(l,mfabs,cos(theta))*sin(mfabs*phi);}
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double int2D(double (*func)(double,double), double a, double b, double c , double d){

   double shift1m = (b - a) / 2;
   double shift1p = (b + a) / 2;
   double shift2m = (d - c) / 2;
   double shift2p = (d + b) / 2;

    mat A;
    A.load("parameters100.txt");

    

vec w = A.col(2);
vec x = A.col(1);
double mysum = 0;

 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
      mysum = mysum + shift2m*shift1m*w(i)*w(j)*func(shift1m*x(i) + shift1p,shift2m*x(j) + shift2p);

}}
  return mysum;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double getTesseralNorm(int l, int m){
double a = 0;
double b = 3.14159265358979;
double c = 0;
double d = 2*b;
   double shift1m = (b - a) / 2;
   double shift1p = (b + a) / 2;
   double shift2m = (d - c) / 2;
   double shift2p = (d + b) / 2;

    mat A;
    A.load("parameters100.txt");
vec ww = A.col(2);
vec xx = A.col(1);
double newX;
double newY;
double mysum = 0;
 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
newX = shift1m*xx(i) + shift1p;
newY = shift2m*xx(j) + shift2p;
      mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*pow(tesseral_spherical_harmNoN(l,m,newX,newY),2)*sin(newX);

}}

  return 1/sqrt(mysum);
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dBig(int l, int m, double r,vec x ,vec y,vec z, double boxSize ,cube V){
double a = 0;
double b = 3.14159265358979;
double c = 0;
double d = 2*b;
   double shift1m = (b - a) / 2;
   double shift1p = (b + a) / 2;
   double shift2m = (d - c) / 2;
   double shift2p = (d + b) / 2;

    mat A;
    A.load("parameters100.txt");
vec ww = A.col(2);
vec xx = A.col(1);
double newX;
double newY;
double mysumNorm = 0;
double mysum = 0;
 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
newX = shift1m*xx(i) + shift1p;
newY = shift2m*xx(j) + shift2p;
      mysumNorm = mysumNorm + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harm(l,m,newX,newY)*tesseral_spherical_harm(l,m,newX,newY)*sin(newX);

}}
cout << mysumNorm;

 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
newX = shift1m*xx(i) + shift1p;
newY = shift2m*xx(j) + shift2p;
      mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harm(l,m,newX,newY)*RtoX3d(r,newX,newY,x,y,z, boxSize, V)*sin(newX);

}}
  return mysum/mysumNorm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dTable(int l, int m, double r,vec x ,vec y,vec z, double boxSize, double norm, cube V){
double a = 0;
double b = 3.14159265358979;
double c = 0;
double d = 2*b;
   double shift1m = (b - a) / 2;
   double shift1p = (b + a) / 2;
   double shift2m = (d - c) / 2;
   double shift2p = (d + b) / 2;

    mat A;
    A.load("parameters100.txt");
vec ww = A.col(2);
vec xx = A.col(1);
double newX;
double newY;
double mysum = 0;

 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
newX = shift1m*xx(i) + shift1p;
newY = shift2m*xx(j) + shift2p;
      mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harmNoN(l,m,newX,newY)*RtoX3d(r,newX,newY,x,y,z, boxSize, V)*sin(newX);
//      mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harmNoN(l,m,newX,newY)*sin(newX);

}}
  return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2D(int l, int m, double r,vec x ,vec y,vec z, double boxSize ,cube V){
double a = 0;
double b = 3.14159265358979;
double c = 0;
double d = 2*b;
   double shift1m = (b - a) / 2;
   double shift1p = (b + a) / 2;
   double shift2m = (d - c) / 2;
   double shift2p = (d + b) / 2;

    mat A;
    A.load("parameters100.txt");
vec ww = A.col(2);
vec xx = A.col(1);
double newX;
double newY;
double mysum = 0;

 for(int i = 0;i < 100 ; i++){for(int j = 0;j < 100 ; j++){
newX = shift1m*xx(i) + shift1p;
newY = shift2m*xx(j) + shift2p;
      mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harm(l,m,newX,newY)*RtoX3d(r,newX,newY,x,y,z, boxSize, V)*sin(newX);

}}
  return mysum;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------ANALYTICL---------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec getGaussAnal(mat coord, double q1, double q2, double q3, double sig){

  vec x= zeros<vec>(2);
  double buff = 0;

 for(int i=0; i < coord.n_rows; i++)
     { 
        buff = coord(i,0) * q1 + coord(i,1)*q2 + coord(i,2)*q3;


         x(0) =x(0) + cos(buff);
         x(1) =x(1) - sin(buff);

        }

// cout << "x(0)=" <<x(0)<< endl;
// cout << "x(1)=" <<x(1)<< endl;

 return 5.56832799683171*sqrt(sig*sig*sig)*x*exp(-(q1*q1 + q2*q2 + q3*q3)*sig*0.5);

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getSphericalToCartCubeX( vec R, vec Theta, vec Phi){
  
  cube X(R.n_elem,Theta.n_elem,Phi.n_elem);

  for(int i=0; i < R.n_rows; i++)
    for(int j=0; j < Theta.n_rows; j++)
      for(int k=0; k < Phi.n_rows; k++){ 
        X(i,j,k) = R(i)*sin(Theta(j))*cos(Phi(k));
        }

 return X;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getSphericalToCartCubeY( vec R, vec Theta, vec Phi){
  
  cube Y(R.n_elem,Theta.n_elem,Phi.n_elem);

  for(int i=0; i < R.n_rows; i++)
    for(int j=0; j < Theta.n_rows; j++)
      for(int k=0; k < Phi.n_rows; k++){ 
        Y(i,j,k) = R(i)*sin(Theta(j))*sin(Phi(k));

        }

 return Y;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube getSphericalToCartCubeZ( vec R, vec Theta, vec Phi){
  
  cube Z(R.n_elem,Theta.n_elem,Phi.n_elem);

  for(int i=0; i < R.n_rows; i++)
    for(int j=0; j < Theta.n_rows; j++)
      for(int k=0; k < Phi.n_rows; k++){ 
        Z(i,j,k) = R(i)*cos(Theta(j));

        }

 return Z;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec getGaussAnalType(mat coord, double q1, double q2, double q3, double sig, vec type){

  sig = 1/sig;
  sig = sig*sig;

  vec x= zeros<vec>(2);
  double buff = 0;
  double buffType = 0;

 for(int i=0; i < coord.n_rows; i++)
     { 
        buff = coord(i,0) * q1 + coord(i,1)*q2 + coord(i,2)*q3;
        buffType = type(i);


         x(0) =x(0) + buffType*cos(buff);
         x(1) =x(1) - buffType*sin(buff);

        }

 return 5.56832799683171*sqrt(sig*sig*sig)*x*exp(-(q1*q1 + q2*q2 + q3*q3)*sig*0.5);

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double vecAbs(vec x){
 
// cout << "x(0)=" <<1.0*x(0)*x(0)<< endl;
// cout << "x(1)=" <<1.0*x(1)*x(1)<< endl;
// cout << "x=" <<1.0*x(0)*x(0) + 1.0*x(1)*x(1)<< endl;

 return sqrt(x(0)*x(0) + x(1)*x(1));

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double vecRe(vec x){
 
// cout << "x(0)=" <<1.0*x(0)*x(0)<< endl;
// cout << "x(1)=" <<1.0*x(1)*x(1)<< endl;
// cout << "x=" <<1.0*x(0)*x(0) + 1.0*x(1)*x(1)<< endl;

 return x(0);

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double vecIm(vec x){
 
// cout << "x(0)=" <<1.0*x(0)*x(0)<< endl;
// cout << "x(1)=" <<1.0*x(1)*x(1)<< endl;
// cout << "x=" <<1.0*x(0)*x(0) + 1.0*x(1)*x(1)<< endl;

 return x(1);

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dTable(int l, int m, double r, double norm, mat coord){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

	mat A;
	A.load("parameters100.txt");
	vec ww = A.col(2);
	vec xx = A.col(1);
	double newX;
	double newY;
	double mysum = 0;
	double mySin = 0;
	// newX = Theta, newY = Phi...
	for(int i = 0;i < 100 ; i++){
		for(int j = 0;j < 100 ; j++){
			newX = shift1m*xx(i) + shift1p;
			newY = shift2m*xx(j) + shift2p;
			mySin = sin(newX);
			mysum = mysum + shift2m*shift1m*ww(i)*ww(j)*tesseral_spherical_harmNoN(l,m,newX,newY)*vecAbs(getGaussAnal(coord, r*mySin*cos(newY),r*mySin*sin(newY), r*cos(newX), 2))*mySin;
		}
	}
	return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dTableType(int l, int m, double r, double norm, mat coord, vec type, double sig, mat A){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

//	mat A(100,2);
//	A.load("parameters100.txt");
//	A.load("P200_both.txt");

	vec ww = A.col(1);
	vec xx = A.col(0);
//	vec ww = A.col(1);
//	vec xx = A.col(0);

	double newX;
	double newY;
	double mysum = 0;
	double mySin = 0;
	// newX = Theta, newY = Phi...
	for(int i = 0;i < ww.n_elem ; i++){
		for(int j = 0;j < ww.n_elem ; j++){
			newX = shift1m*xx[i] + shift1p;
			newY = shift2m*xx[j] + shift2p;
			mySin = sin(newX);
			mysum = mysum + shift2m*shift1m*ww[i]*ww[j]*tesseral_spherical_harmNoN(l,m,newX,newY)*vecAbs(getGaussAnalType(coord, r*mySin*cos(newY),r*mySin*sin(newY), r*cos(newX), sig, type))*mySin;
		}
	}
	return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dRe(int l, int m, double r, double norm, mat coord, vec type, double sig, mat A){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

//	mat A(100,2);
//	A.load("parameters100.txt");
//	A.load("P200_both.txt");

	vec ww = A.col(1);
	vec xx = A.col(0);
//	vec ww = A.col(1);
//	vec xx = A.col(0);

	double newX;
	double newY;
	double mysum = 0;
	double mySin = 0;
	// newX = Theta, newY = Phi...
	for(int i = 0;i < ww.n_elem ; i++){
		for(int j = 0;j < ww.n_elem ; j++){
			newX = shift1m*xx[i] + shift1p;
			newY = shift2m*xx[j] + shift2p;
			mySin = sin(newX);
			mysum = mysum + shift2m*shift1m*ww[i]*ww[j]*tesseral_spherical_harmNoN(l,m,newX,newY)*vecRe(getGaussAnalType(coord, r*mySin*cos(newY),r*mySin*sin(newY), r*cos(newX), sig, type))*mySin;
		}
	}
	return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dIm(int l, int m, double r, double norm, mat coord, vec type, double sig, mat A){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

//	mat A(100,2);
//	A.load("parameters100.txt");
//	A.load("P200_both.txt");

	vec ww = A.col(1);
	vec xx = A.col(0);
//	vec ww = A.col(1);
//	vec xx = A.col(0);

	double newX;
	double newY;
	double mysum = 0;
	double mySin = 0;
	// newX = Theta, newY = Phi...
	for(int i = 0;i < ww.n_elem ; i++){
		for(int j = 0;j < ww.n_elem ; j++){
			newX = shift1m*xx[i] + shift1p;
			newY = shift2m*xx[j] + shift2p;
			mySin = sin(newX);
			mysum = mysum + shift2m*shift1m*ww[i]*ww[j]*tesseral_spherical_harmNoN(l,m,newX,newY)*vecIm(getGaussAnalType(coord, r*mySin*cos(newY),r*mySin*sin(newY), r*cos(newX), sig, type))*mySin;
		}
	}
	return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
int TabulateTesseralPos(int l, int m, mat A){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

	vec ww = A.col(1);
	vec xx = A.col(0);

	double theta;
	double phi;

	for(int i = 0;i < ww.n_elem ; i++){
		for(int j = 0;j < ww.n_elem ; j++){
			theta = shift1m*xx[i] + shift1p;
			phi = shift2m*xx[j] + shift2p;
			cout << setprecision(12) <<  tesseral_spherical_harmNoN(l,m,theta,phi) << " ";
		}
			cout <<  endl;
	}
	return 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double tesseralInt2dTableTypeTable(int l, int m, double r, double norm, mat coord, vec type, double sig, mat A, mat tessTable){
	double a = 0;
	double b = 3.14159265358979;
	double c = 0;
	double d = 2*b;
	double shift1m = (b - a) / 2;
	double shift1p = (b + a) / 2;
	double shift2m = (d - c) / 2;
	double shift2p = (d + b) / 2;

	vec ww = A.col(1);
	vec xx = A.col(0);

	double newX;
	double newY;
	double mysum = 0;
	double mySin = 0;
	// newX = Theta, newY = Phi...
	for(int i = 0;i < ww.n_elem ; i++){
		for(int j = 0;j < ww.n_elem ; j++){
			newX = shift1m*xx[i] + shift1p;
			newY = shift2m*xx[j] + shift2p;
			mySin = sin(newX);
			mysum = mysum + shift2m*shift1m*ww[i]*ww[j]*tessTable(i,j)*vecAbs(getGaussAnalType(coord, r*mySin*cos(newY),r*mySin*sin(newY), r*cos(newX), sig, type))*mySin;
		}
	}
	return mysum*norm;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec hydrogenRDF(int n, double z,double ao,double norm,vec r){
  vec rhon = 2*z*r/ao/(double) n; 
  vec exprhon = exp(-rhon*0.5);

  if(n==1){
	  return (norm*exprhon*2);
  }
  if(n==2){
	  return (norm*exprhon%(2 - rhon)/sqrt(2)*0.5);
  }
  if(n==3){
   	return (norm*exprhon%(6-6*rhon+rhon%rhon)/9/sqrt(3)); 
  }
  if(n==4){
   	return (norm*exprhon/96%(24 - 36*rhon+12*rhon%rhon - rhon%rhon%rhon)); 
  }

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double integ1D(vec R, vec F, mat GL, double rcut){

  double result =	0.5*rcut*sum(GL.col(1)%F);

  return result;


}
//----------------------------------------------------------------------------------------------------------------------------------------------
#endif 

