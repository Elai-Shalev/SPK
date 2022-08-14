#define PY_SSIZE_T_CLEAN
// include "Python.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* Global Variables */
int MAX_ITER;
int K;
double const EPSILON;
int dim;
int num_of_vectors;

/* Macros */
#define SQR(x) ((x)*(x))

/* Structs */
typedef struct {
    double* coordinate;
    int cluster;
} Vector;

//Read the File
//Matricies are flattened arrays 

double* calc_weighted_matrix(double* points){
    int i;
    int j;
    double* weighted_adj_matrix = (double*)malloc(num_of_vectors, sizeof(double)*pow(num_of_vectors,2));
    for(i = 0; i<n; i++){
        for(j=i; j<n; j++){
            if(i==j){
                weighted_adj_matrix[n*i+j] = 0;
            }
            else{
                weighted_adj_matrix[n*i+j] = exp(-0.5*norm2(points,i,j));

            }
        }
    }

}
double norm2(double* matrix, int i, int j){

}

double* calc_diagonal_deg_matrix(double* mat){
}
double* calc_lnorm_matrix(double* mat){
}
double* calc_eigen(double* mat){
}


int main(int argc, char* args[]){
    return 0;
}
