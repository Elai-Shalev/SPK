#ifndef SPKMEANS_HEADER_FILE_H
#define SPKMEANS_HEADER_FILE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Global Variables */
int MAX_ITER;
int K;
double const EPSILON;
const char* filename;
int dim;
int num_of_vectors;

/* Macros */
#define SQR(x) ((x)*(x))
#define ABS(x) ((x<0)?-x:x)
#define NULL_ERROR_CHECK(x) {\
    if(x == NULL){\
        printf("An Error Has Occurred");\
        exit(1);\
    }}


/* Structs */
typedef struct {
    double* coordinate;
    int cluster;
} Vector;

/* Dimension Reduction Functions */
/* Calculates T according to the SPK algorithm */
double* dimension_reduction_spk(double* points);
/* Transposes inplace a double* matrix */
void transpose_square_matrix_inplace(double** matrix, int size);
/* Calls for calc_norma and normalizes the matrix inplace */
void normalize_first_k_vectors(double** eigen_vectors_matrix, int k);
/* Calculates the norm of a vector */
double calc_norma(double* vector, int size);
/* Implementation of the Eigen-gap Heuristic */
int determine_k(double* sorted_eigen_values);
/* Converts to a pointer-array representation of a matrix */
double** convert_double_array_to_matrix(double* array, 
                                        int num_rows,
                                        int num_cols);
/* Sorts Eigenvalues and Eigenvectors  in decreasing order */
void sort_eigen_v(double* eigen_values, double** eigen_vectors);
/* Calculates the Weighted Adjacency Matrix */
double* calc_weighted_matrix(double* points);
/* Returnes the Euclidian Distance between two vectors */
double euclidian_dist(double* points, int i, int j);
/* Calculates the Diagonal Degree Matrix of weighted adjacency matrix */
double* calc_diagonal_deg_matrix(double* mat);
/* Calculates the Normalized Graph Laplacian of the input matrix*/
double* calc_lnorm_matrix(double* points);
/* Returns the maximum off-diagonal entry of the matrix */
int* max_abs_off_diagonal_entry(double* mat, int size);
/* Calculates the sum of squares for all off-diagonal entries */
double sum_squares_off_diagonal(double* mat, int size);
/* Creates initial rotation matrix for Jacobi algorithm  */
double* create_initial_p_matrix(int i, int j, double c, double s);
/* Calculates the Eigenvalues and Eigenvectors of a Normalized Laplacian */
double** calc_eigen(double* A);
/* Pivot step implementation according to the Jacobi algorithm */
double * pivot_jacobi(double * A, int max_i, int max_j);
/* Implementation for rotation-matrix multiplication */
void rotation_matrix_multiply_simplified(double * mat, 
                                         int a, int b, 
                                         double c, double s);
/* Read Input File */
double* read_file(char* file_in);
/* Prints Matrix */
void print_matrix(double* matrix, char deli, int num_rows, int num_cols);
/* Prints double array */
void print_double_array(double* arr, char deli, int length);
/* Converts pointer-array matrix representation to single array */
double* convert_matrix_to_double_array(double** matrix, 
                                        int new_num_rows, 
                                        int new_num_cols);

/* Clustering Functions */
/* Returns square distance of two Vectors */
double square_distance(Vector* vec1, Vector* vec2);
/* Assigns to nearest cluster */
void assign_to_nearest_cluster(Vector* vec, Vector** centroids);
/* Initilizes new zero Vector */
Vector* new_zero_vector();
/* Initilizes new zero Vector */
Vector* new_deep_copy_vector(Vector* vec);
/* Calculates inplace the sum of a vector */
void vectoric_sum(Vector* vec1, Vector* vec2);
/* Calculates inplace the quotient between two vectors */
void divide_vector(Vector* vec, int divisor);
/* Updates centroids according to kmeans algorithm */
int update_centroids(Vector** vector_list, Vector** centroids);
/* Runs kmeans clustering algorithm */
void run_kmeans(Vector** vectors, Vector** centroids);
/* Convets double array to vector array */
Vector** double_array_to_vector_array(double* vector_array);

/* C-API Functions */
/* Runs kmeans++ implemenation*/
Vector** fit_c(double* vector_array, double* centroid_array);
/* Runs C implementation of WAM */
void wam_c(double* points);
/* Runs C implementation of DDG */
void ddg_c(double* points);
/* Runs C implementation of LNORM */
void lnorm_c(double* points);
/* Runs C implementation of JACOBI */
void jacobi_c(double* points);

#endif
