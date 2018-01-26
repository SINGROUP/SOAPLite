#ifndef CLUS_GEO   /* Include guard */
#define CLUS_GEO


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

//--------------------------------------------------------------------------

uvec getSurf(mat A,double surfDefN){

    srand(time(NULL));

    double extra = 10.0; //10.0 for MoS2, 5.0 for AuCu
//    double surfDefN = 3.0;//5.0 for MoS2, 2.0 for AuCu
    double surfDefN2 = pow(surfDefN,2);
    bool q;

    double xmax = A.col(0).max(); double xmaxInx = A.col(0).index_max();
    double ymax = A.col(1).max(); double ymaxInx = A.col(1).index_max();
    double zmax = A.col(2).max(); double zmaxInx = A.col(2).index_max();

    double xmin = A.col(0).min(); double xminInx = A.col(0).index_min();
    double ymin = A.col(1).min(); double yminInx = A.col(1).index_min();
    double zmin = A.col(2).min(); double zminInx = A.col(2).index_min();

    xmax = xmax + extra; xmin = xmin - extra;
    ymax = ymax + extra; ymin = ymin - extra;
    zmax = zmax + extra; zmin = zmin - extra;


    mat randPos = randu<mat>(1000000,3);

    randPos.col(0) = (xmax - xmin) * randPos.col(0) + xmin ;
    randPos.col(1) = (ymax - ymin) * randPos.col(1) + ymin ;
    randPos.col(2) = (zmax - zmin) * randPos.col(2) + zmin ;

    vec V = zeros<vec>(A.n_rows);

    uvec a(1000000) ;
    a.ones();

    a = 100000 * a;

    for(int j=0; j < randPos.n_rows; j++){
        for(int i=0; i < A.n_rows; i++){
          V(i) = sum(pow(randPos.row(j) - A.row(i),2)); 
        }

        if(all(V > surfDefN2)) {
            a(j) = V.index_min(); 
        } 
     }

    a = unique(sort(a));

    uvec b = a.rows(0,a.n_rows - 2);

  return b;

}

uvec getNonSurf(mat x, uvec q){

    uvec p = find(x.col(1) > -10000000 );


    uvec r(p.n_rows - q.n_rows);

    int n=0;
    int trigger = 0;

    for(int i =0; i < p.n_rows; i++){
        trigger = 0;
        for(int j =0; j < q.n_rows; j++){

       if(p(i) == q(j)) {trigger = 1;}

      }

      if(trigger == 0){ r(n) = p(i); n++;}

      }

  return r;

}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#endif // 
