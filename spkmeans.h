#ifndef SPKMEANS_HEADER_FILE_H
#define SPKMEANS_HEADER_FILE_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
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
double* dimension_reduction_spk(double* points);
double* transpose_matrix_copy(double** matrix, int rows, int cols);
void transpose_square_matrix_inplace(double** matrix, int size);
void normalize_first_k_vectors(double** eigen_vectors_matrix, int k);
double calc_norma(double* vector, int size);
int determine_k(double* sorted_eigen_values);
double** convert_double_array_to_matrix(double* array, 
                                        int num_rows,
                                        int num_cols);
void sort_eigen_v(double* eigen_values, double** eigen_vectors);
double* calc_weighted_matrix(double* points);
double euclidian_dist(double* points, int i, int j);
double* calc_diagonal_deg_matrix(double* mat);
double* calc_lnorm_matrix(double* points);
int* max_abs_off_diagonal_entry(double* mat, int size);
double sum_squares_off_diagonal(double* mat, int size);
double* create_initial_p_matrix(int i, int j, double c, double s);
double** calc_eigen(double* A);
double * pivot_jacobi(double * A, int max_i, int max_j);
void rotation_matrix_multiply_simplified(double * mat, 
                                         int a, int b, 
                                         double c, double s);
double* read_file(char* file_in);
void print_matrix(double* matrix, char deli, int num_rows, int num_cols);
void print_double_array(double* arr, char deli, int length);
void print_matrix_doublestar(double** matrix, char deli, 
                             int num_rows, int num_cols);
double* convert_matrix_to_double_array(double** matrix, 
                                        int new_num_rows, 
                                        int new_num_cols);

/* Clustering Functions */
double square_distance(Vector* vec1, Vector* vec2);
void assign_to_nearest_cluster(Vector* vec, Vector** centroids);
Vector* new_zero_vector();
Vector* new_deep_copy_vector(Vector* vec);
void vectoric_sum(Vector* vec1, Vector* vec2);
void divide_vector(Vector* vec, int divisor);
int update_centroids(Vector** vector_list, Vector** centroids);
void run_kmeans(Vector** vectors, Vector** centroids);
Vector** double_array_to_vector_array(double* vector_array);
double* python_list_to_c_array(PyObject* float_list);
PyObject* c_array_to_python_list(double* float_list);
Vector** fit_c(double* vector_array, double* centroid_array);
void wam_c(double* points);

#endif