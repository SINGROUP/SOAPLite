#ifndef MYARMADILLO/* Include guard */
#define MYARMADILLO

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>


using namespace std;
using namespace arma;
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube meshgrid3D(vec a, vec b, vec c, int d){

 int na=a.n_rows;
 int nb=b.n_rows;
 int nc=c.n_rows;
 cube X(na,nb,nc);

 for(int i=0; i < na; i++){
   for(int j=0; j < nb; j++){
     for(int k=0; k < nc; k++){
         if(d==0){X(i,j,k) = a(i) ;}
	 else if(d==1){X(i,j,k) = b(j) ;}
	 else{X(i,j,k) = c(k) ;}

     }}}


return X;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
mat meshgrid2D(vec a, vec b, int d){

 int na=a.n_rows;
 int nb=b.n_rows;
 mat X(na,nb);

 for(int i=0; i < na; i++){
   for(int j=0; j < nb; j++){
         if(d==0){X(i,j) = a(i) ;}
	 else if(d==1){X(i,j) = b(j) ;}
     }}


return X;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double getValue(cube A, vec a, vec b, vec c, double x, double y, double z){

    vec xx = abs(a-x-0.00001);
    vec yy = abs(b-y-0.00001);
    vec zz = abs(c-z-0.00001);

    int aInx=index_min(xx);
    int bInx=index_min(yy);
    int cInx=index_min(zz);

    return A(aInx,bInx,cInx);


}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
mat  myGauss2(double x0,double y0,mat X,mat Y, double sig)
{
  double oversigma = 0.5/sig/sig;
 mat A(X.n_rows,X.n_cols);
for(int i=0; i < A.n_rows; i++)

    { 
 for(int j=0; j < A.n_cols; j++)
     { 
            A(i,j)=exp( -(pow((X(i,j) - x0),2) +  pow((Y(i,j) - y0),2))*oversigma) ;
          }
       }
 return A;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
mat  fftshift2d(mat X)
{
  int nrows = X.n_rows;
  int ncols = X.n_cols;

mat center = zeros<mat>(nrows,nrows);
for(int i=0; i < center.n_rows; i++)
    { 
     for(int j=0; j < center.n_cols; j++)
         { 
           if(i < nrows/2 && j < nrows/2 ) center(i+nrows/2,j+nrows/2) = X(i,j);
            if(i > nrows/2 && j < nrows/2 ) center(i-nrows/2,j+nrows/2) = X(i,j);
            }
       }
for(int i=0; i < center.n_rows; i++)
    { 
     for(int j=0; j < center.n_cols/2; j++)
         { 
            center(i,j) = center(i,nrows - j - 1);
            }
       }

return center;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube  fftshift3d(cube X)
{
  int nrows = X.slice(0).n_rows;
  int ncols = X.slice(0).n_cols;
  int nslices = X.n_slices;
//  int nslicesDouble = 2*(X.n_slices - 1);

cube center = zeros<cube>(nrows,nrows,nslices);
for(int i=0; i < nrows; i++)
     for(int j=0; j < ncols; j++)
        for(int k=0; k < nslices; k++)
           { 
             if(i < nrows/2 && j < nrows/2  && k < nrows/2) center(i+nrows/2,j+nrows/2,k+nrows/2) = X(i,j,k);
             if(i > nrows/2 && j < nrows/2  && k < nrows/2) center(i-nrows/2,j+nrows/2,k+nrows/2) = X(i,j,k);
             if(i < nrows/2 && j > nrows/2  && k < nrows/2) center(i+nrows/2,j-nrows/2,k+nrows/2) = X(i,j,k);
             if(i > nrows/2 && j > nrows/2  && k < nrows/2) center(i-nrows/2,j-nrows/2,k+nrows/2) = X(i,j,k);
              }
for(int i=0; i < nrows; i++) for(int j=0; j < ncols; j++) for(int k=0; k < nslices/2; k++)
            center(i,j,k) = center(i,j,nslices - k - 1);

return center;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube  myGauss3(double x0,double y0,double z0,cube X,cube Y,cube Z, double sig)
{
 cube A(X.slice(0).n_rows,X.slice(0).n_cols, X.n_slices);

for(int i=0; i < A.n_rows; i++)
    { 
 for(int j=0; j < A.n_cols; j++)
     { 
   for(int k=0; k < A.n_cols; k++)
       { 
            A(i,j,k)=exp( -(pow((X(i,j,k) - x0),2) +  pow((Y(i,j,k) - y0),2) + pow((Z(i,j,k) - z0),2)) /2/sig/sig) ;

          }
       }
   }
 return A;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double  RtoX3d(double r,double Theta,double Phi,vec x,vec y,vec z, double  boxSize, cube V)
{
double xVal = r*sin(Theta)*cos(Phi);
double yVal = r*sin(Theta)*sin(Phi);
double zVal = r*cos(Theta);

double ddx = (x(1) - x(0));
double ddy = (y(1) - y(0));
double ddz = (z(1) - z(0));

int Ix =floor((xVal + boxSize)/ddx) ;
int Iy =floor((yVal + boxSize)/ddy) ;
int Iz =floor((zVal + boxSize)/ddz) ;

cout << "Ix: " << Ix << " Iy: " << Iy << " Iz:" << Iz << endl;

double dx = (xVal - x(Ix))/ddx;
double dy = (yVal - y(Iy))/ddy;
double dz = (zVal - z(Iz))/ddz;

double V00 = V(Ix,Iy,Iz)*(1 - dx) + V(Ix + 1,Iy,Iz)*dx;
double V01 = V(Ix,Iy,Iz+1)*(1 - dx) + V(Ix+1,Iy,Iz+1)*dx;
double V10 = V(Ix,Iy+1,Iz)*(1 - dx) + V(Ix+1,Iy+1,Iz)*dx;
double V11 = V(Ix,Iy+1,Iz+1)*(1 - dx) + V(Ix+1,Iy+1,Iz+1)*dx;

double V0=V00*(1 - dy) + V10*dy ;
double V1=V01*(1 - dy) + V11*dy ;

double C=V0*(1 - dz) + V1*dz ;

return C;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double  RtoX3dSimple(double r,double Theta,double Phi,vec x,vec y,vec z, double  boxSize, cube V)
{
double xVal = r*sin(Theta)*cos(Phi);
double yVal = r*sin(Theta)*sin(Phi);
double zVal = r*cos(Theta);

double ddx = (x(1) - x(0));
double ddy = (y(1) - y(0));
double ddz = (z(1) - z(0));

int Ix =round((xVal + boxSize)/ddx) ;
int Iy =round((yVal + boxSize)/ddy) ;
int Iz =round((zVal + boxSize)/ddz) ;


double C=V(Ix,Iy,Iz);

return C;
}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube  getR3d(double r,double Theta,double Phi,vec x,vec y,vec z, double  boxSize, cube V)
{
int nrows =  V.slice(0).n_rows;
int ncols =   V.slice(0).n_cols;
int nslices =  V.n_slices;
  cube C = zeros<cube>(nrows,ncols, nslices);
for(int i=0; i < nrows; i++) {cout << nrows - i << endl; for(int j=0; j < ncols; j++)for(int k=0; k < nslices; k++) {
  
        C(i,j,k) = RtoX3d(r/nrows * (double) i,Theta/ncols * (double) j, Phi/nslices * (double) k ,x,y,z,boxSize,V);

       }}

return C;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube  getR3dSimple(double r,double Theta,double Phi,vec x,vec y,vec z, double  boxSize, cube V)
{
int nrows =  V.slice(0).n_rows;
int ncols =   V.slice(0).n_cols;
int nslices =  V.n_slices;
  cube C = zeros<cube>(nrows,ncols, nslices);
for(int i=0; i < nrows; i++) {cout << nrows - i << endl; for(int j=0; j < ncols; j++)for(int k=0; k < nslices; k++) {
  
        C(i,j,k) = RtoX3dSimple(r/nrows * (double) i,Theta/ncols * (double) j, Phi/nslices * (double) k ,x,y,z,boxSize,V);

       }}

return C;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
mat  getR2d(double r,double Theta,double Phi,vec x,vec y,vec z, double  boxSize, cube V)
{
int  nrows =  V.slice(0).n_rows ;
int ncols =   V.slice(0).n_cols ;
//int  nrows =  pow(2,3); // DOESN'T WORK, WHY !?!?
//int ncols =  pow(2,3) ; // DOESN'T WORK, WHY !?!?
  mat C = zeros<mat>(nrows,ncols);
for(int j=0; j < nrows; j++) {cout << nrows - j << endl; for(int k=0; k < ncols; k++){
  
        C(j,k) = RtoX3d(r,Theta/nrows * (double) j, Phi/ncols * (double) k ,x,y,z,boxSize,V);
//cout << C(j,k)<< endl;
       }}

return C;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double sum3d(cube A){
  double buffDiff3 = 0;
for(int i=0; i < A.slice(0).n_rows; i++)  for(int j=0; j < A.slice(0).n_cols; j++) for(int k=0; k < A.n_slices; k++)
           buffDiff3 = buffDiff3 + abs(A(i,j,k));
return buffDiff3;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
double sum2d(mat A){
  double buffDiff2 = 0;
for(int i=0; i < A.n_rows; i++)  for(int j=0; j < A.n_cols; j++)
           buffDiff2 = buffDiff2 + abs(A(i,j));
return buffDiff2;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
cube gauss3d(mat A, cube X, cube Y, cube Z){
  cube B = zeros<cube>(X.slice(0).n_rows,X.slice(0).n_cols, X.n_slices);
  double buffDiff2 = 0;
for(int i=0; i < A.n_rows; i++)
           B = B +  myGauss3(A(i,0), A(i,1), A(i,2),X,Y,Z,1.5);
return B;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
mat posAve(mat A){

  rowvec v = zeros<rowvec>(A.n_cols);
  mat B = A;
for(int i=0; i < A.n_rows; i++)
    { 
         v = v + A.row(i);
       }

v = v / A.n_rows;

for(int i=0; i < A.n_rows; i++)
    { 
     B.row(i) = A.row(i) - v; 
       }

return B;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
mat posAveExt(mat A, mat B){

  rowvec v = zeros<rowvec>(A.n_cols);

for(int i=0; i < A.n_rows; i++)
    { 
         v = v + A.row(i);
       }

v = v / A.n_rows;

for(int i=0; i < B.n_rows; i++)
    { 
     B.row(i) = B.row(i) - v; 
       }

return B;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
mat rotate3d(mat A,double myalpha, double mybeta, double mygamma){
  mat Rx = zeros<mat>(3,3) ;
  mat Ry = zeros<mat>(3,3) ;
  mat Rz = zeros<mat>(3,3) ;
  mat rotationMat = zeros<mat>(3,3) ;
  

  mat B = zeros<mat>(A.n_rows,A.n_cols);

  Rx(0,0)= 1; Rx(0,1)= 0;            Rx(0,2)= 0;
  Rx(1,0)= 0; Rx(1,1)= cos(myalpha); Rx(1,2)= -sin(myalpha);
  Rx(2,0)= 0; Rx(2,1)= sin(myalpha); Rx(2,2)= cos(myalpha);


  Ry(0,0)= cos(mybeta);   Ry(0,1)= 0;  Ry(0,2)=sin(mybeta);
  Ry(1,0)= 0;             Ry(1,1)= 1;  Ry(1,2)= 0;
  Ry(2,0)= -sin(mybeta);  Ry(2,1)= 0;  Ry(2,2)= cos(mybeta);

  Rz(0,0)= cos(mygamma);     Rz(0,1)= -sin(mygamma); Rz(0,2)=0;
  Rz(1,0)= sin(mygamma);    Rz(1,1)= cos(mygamma);  Rz(1,2)= 0;
  Rz(2,0)= 0;               Rz(2,1)= 0;             Rz(2,2)= 1;

 rotationMat = Rx*Ry*Rz;

  for(int i = 0; i < A.n_rows; i++){ 
  B.row(i) = A.row(i)*rotationMat; 
  } 

return B;

}
//----------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------
vec getTypeVal(string* A, mat coord){
  int coordN = coord.n_rows;
  vec x = zeros<vec>(coordN);


for(int i=0; i < coordN; i++)
    { // cout << A[i] << endl;
           if(A[i] == "H") x(i) = 1;
           else if(A[i] == "C") x(i) = 6;
           else if(A[i] == "N") x(i) = 7;
           else if(A[i] == "O") x(i) = 8;
           else if(A[i] == "F") x(i) = 9;
           else if(A[i] == "Mo") x(i) = 42;
           else if(A[i] == "S") x(i) = 16;
           else if(A[i] == "Cu") x(i) = 29;
           else if(A[i] == "Au") x(i) = 79;
           else{cout << "New Atom Introduced, must define atom number!" << endl; exit(1);}
       }


return x;

}
//----------------------------------------------------------------------------------------------------------------------------------------------


#endif 

// (b - a)^3/2 int(-1,1)int(-1,1)int(-1,1)dxdydz f ( x-beta)/alfa, (x - beta)/alfa, (x-beta)/alfa) where alpha = 2/(b - a), beta = -(a+b)/(b-a)
