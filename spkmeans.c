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

double* calc_weighted_matrix(double* points, int n){
}
double* calc_diagonal_deg_matrix(double* mat, int n){
}
double* calc_lnorm_matrix(double* mat, int n){
}
double* calc_eigen(double* mat, int n){
}


int main(int argc, char* args[]){
    return 0;
}
