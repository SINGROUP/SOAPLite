#ifndef MY_HEADER_H
# define MY_HEADER_H

#include<armadillo>
#include<mpi.h>

using namespace arma;

//----------------------------------------------------------------------
void MPI_SendMat(mat &A, int tocore, int tag){
        int kk = 0;
        int mysize = A.size();
        double buff[mysize];
        for(int i=0; i <A.n_rows; i++){ for(int j=0; j <A.n_cols; j++){
                
                buff[kk] = A(i,j);
                kk = kk + 1; 
                
        }}
        MPI_Send(buff, mysize, MPI_DOUBLE, tocore, tag, MPI_COMM_WORLD);
    }
//----------------------------------------------------------------------
void MPI_RecvMat (mat &A, int fromcore, int tag){
        MPI_Status status;
        int kk = 0;
        int mysize = A.size();
        double buff[mysize];
        MPI_Recv(buff, mysize ,MPI_DOUBLE,fromcore, tag,MPI_COMM_WORLD, &status);
        for(int i=0; i <A.n_rows; i++){ for(int j=0; j <A.n_cols; j++){
                A(i,j) = buff[kk];
                kk = kk + 1; 
        }}
    }
//----------------------------------------------------------------------
void MPI_SendUvec(uvec &q, int tocore, int tag){
        int mysize = q.n_rows;
        double buff[mysize];
        for(int i=0; i < mysize ; i++){
                
                buff[i] = q(i);
                
        }
        MPI_Send(buff, mysize, MPI_DOUBLE, tocore, tag, MPI_COMM_WORLD);
    }
//----------------------------------------------------------------------   
void MPI_RecvUvec (uvec &q, int fromcore, int tag){
        MPI_Status status;
        int qsize = q.n_rows;
        double buff[qsize];
        MPI_Recv(buff, qsize, MPI_DOUBLE,fromcore, tag,MPI_COMM_WORLD, &status);
        for(int i=0; i < qsize ; i++){

                q(i) = buff[i];

        }
    }
//---------------------------------------------------------------------- 
//----------------------------------------------------------------------
void MPI_SendVec(vec &q, int tocore, int tag){
        int mysize = q.n_rows;
        double buff[mysize];
        for(int i=0; i < mysize ; i++){
                
                buff[i] = q(i);
                
        }
        MPI_Send(buff, mysize, MPI_DOUBLE, tocore, tag, MPI_COMM_WORLD);
    }
//----------------------------------------------------------------------   
//----------------------------------------------------------------------   
void MPI_RecvVec (vec &q, int fromcore, int tag){
        MPI_Status status;
        int qsize = q.n_rows;
        double buff[qsize];
        MPI_Recv(buff, qsize, MPI_DOUBLE,fromcore, tag,MPI_COMM_WORLD, &status);
        for(int i=0; i < qsize ; i++){

                q(i) = buff[i];

        }
    }
//---------------------------------------------------------------------- 

#endif
