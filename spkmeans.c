#define PY_SSIZE_T_CLEAN
// include "Python.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"


/* Global Variables */
int MAX_ITER;
int K;
double const EPSILON;
double const JACOBIAN_EPSILON = 0.00001;
int const JACOBIAN_MAX_ITER = 100;
int dim;
int num_of_vectors;

/* Macros */
#define SQR(x) ((x)*(x))
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

double* dimension_reduction_spk(double* points){
    double* l_norm; 
    double** data;
    double* eigen_values;
    double* eigen_vectors;
    double** eigen_vectors_matrix;
    double* normalized_eigen_vectors;
    int k;
    l_norm = calc_lnorm_matrix(points);
    data = calc_eigen(l_norm);
    eigen_values = data[0];
    eigen_vectors = data[1];
    eigen_vectors_matrix = convert_double_array_to_matrix(
        eigen_vectors,num_of_vectors,num_of_vectors);
    transpose_square_matrix_inplace(eigen_vectors_matrix, num_of_vectors);
    sort_eigen_v(eigen_values, eigen_vectors_matrix);
    k = determine_k(eigen_values);
    normalize_first_k_vectors(eigen_vectors_matrix, k);
    normalized_eigen_vectors = transpose_matrix_copy(eigen_vectors_matrix, 
                                                     k, num_of_vectors);
    return normalized_eigen_vectors; /* T in the algorithm */
}

double* transpose_matrix_copy(double** matrix, int rows, int cols)
{
    double* transposed_matrix;
    int i, j;
    transposed_matrix = (double*)malloc(sizeof(double)*rows*cols);
    NULL_ERROR_CHECK(transposed_matrix);

    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            transposed_matrix[i*rows+j] = matrix[j][i];
        }
    }
    return transposed_matrix;
}

void transpose_square_matrix_inplace(double** matrix, int size){
    int i, j;
    double temp;

    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

void normalize_first_k_vectors(double** eigen_vectors_matrix, int k){
    int i,j;
    double curr_norma;
    for(i = 0; i<k; i++){
        curr_norma = calc_norma(eigen_vectors_matrix[i],num_of_vectors);
        for(j=0; j<k; j++){
            eigen_vectors_matrix[i][j] /= curr_norma;
        }
    }
}

double calc_norma(double* vector, int size){
    double sum;
    int i;
    for(i=0; i<num_of_vectors; i++){
        sum += SQR(vector[i]);
    }
    return sqrt(sum);
}

int determine_k(double* sorted_eigen_values){
    int i, k;
    double gap, curr_gap;
    gap = INFINITY*-1;
    for(i=0; i<num_of_vectors-1; i++){
        curr_gap = abs(sorted_eigen_values[i] - sorted_eigen_values[i-1]);
        if(curr_gap > gap){
            gap = curr_gap;
            k=i;
        }
    }
    return k;
}

double** convert_double_array_to_matrix(double* array, 
                                        int size_row, 
                                        int num_rows){
    int i;
    double** mat = (double**)malloc(sizeof(double*)*num_rows);
    NULL_ERROR_CHECK(mat);
    for(i=0; i<num_rows; i++){
        mat[i] = &array[i*size_row];
    }
    return mat;
}

void sort_eigen_v(double* eigen_values, double** eigen_vectors){
    int i, j;
    double temp_eigenvalue;
    double* temp_eigenvector;

    for(i = 0; i < num_of_vectors; i++){
        for (j = i; j < num_of_vectors; j++){
            if (eigen_values[j] < eigen_values[i]){
                temp_eigenvalue = eigen_values[i];
                temp_eigenvector = eigen_vectors[i];
                eigen_values[i] = eigen_values[j];
                eigen_vectors[i] = eigen_vectors[j];
                eigen_values[j] = temp_eigenvalue;
                eigen_vectors[j] = temp_eigenvector;
            }
        }
    }
}

double* calc_weighted_matrix(double* points){
    int i;
    int j;
    double norm_i_j;
    double* weighted_adj_matrix = 
    (double*)malloc(sizeof(double)*SQR(num_of_vectors));
    NULL_ERROR_CHECK(weighted_adj_matrix);
    for(i = 0; i<num_of_vectors; i++){
        for(j=i; j<num_of_vectors; j++){
            if(i==j){
                weighted_adj_matrix[num_of_vectors*i+j] = 1;
            }
            else{
                double norm_i_j = euclidian_dist(points,i,j);
                weighted_adj_matrix[num_of_vectors*i+j] = exp(-0.5*norm_i_j);
                weighted_adj_matrix[num_of_vectors*j+i] = exp(-0.5*norm_i_j);
            }
        }
    }
    return weighted_adj_matrix;

}
double euclidian_dist(double* points, int i, int j){
    int k;
    double sum_diff = 0;
    for(k=0; k<dim; k++){
        sum_diff += pow((points[dim*i+k] - points[dim*j+k]),2);
    }
    return sqrt(sum_diff);
}

double* calc_diagonal_deg_matrix(double* mat){
    int z,i;
    double sum_row=0;
    double* diagonal_deg_matrix = 
    (double*)malloc(sizeof(double)*num_of_vectors);
    NULL_ERROR_CHECK(diagonal_deg_matrix);
    for(i=0; i<num_of_vectors; i++){
        sum_row = 0;
        for(z = 0; z<num_of_vectors; z++){
            sum_row += mat[num_of_vectors*i+z];
        }
        diagonal_deg_matrix[i] = sum_row;
    }
    return diagonal_deg_matrix;
}

double* calc_lnorm_matrix(double* points){
    double* lnorm_matrix;
    double* weighted_matrix = calc_weighted_matrix(points);
    double* diagonal_deg_matrix = calc_diagonal_deg_matrix(weighted_matrix);
    int i,j;

    lnorm_matrix = (double*)malloc(sizeof(double)*SQR(num_of_vectors));
    NULL_ERROR_CHECK(lnorm_matrix);
    
    for (i = 0; i < num_of_vectors; i++){
        for (j = 0; j < num_of_vectors; j++){
            lnorm_matrix[i*num_of_vectors+j] = 
                                    -1* pow(diagonal_deg_matrix[i], -0.5)*
                                    weighted_matrix[i*num_of_vectors+j]*
                                    pow(diagonal_deg_matrix[j], -0.5);
            
            if (i==j){
                lnorm_matrix[i*num_of_vectors+j] += 1.0000;
            }
        }
    }

    return lnorm_matrix;
}


int* max_abs_off_diagonal_entry(double* mat, int size){
    int* data = (int*)malloc(sizeof(int)*2);
    NULL_ERROR_CHECK(data);
    int i, j;

    data[0] = 0;
    data[1] = 1;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            if ((i != j) &&
                (abs(mat[(data[0])*size+(data[1])]) < abs(mat[i*size+j]))){
                data[0] = i;
                data[1] = j;
            }
        }
    } 

    return data;
}


double sum_squares_off_diagonal(double* mat, int size){
    double sum = 0;
    int i,j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            if (i != j){
                sum += mat[i*size+j];
            }
        }
    } 

    return sum;
}

double* create_initial_p_matrix(int i, int j, double c, double s){
    double* P = (double*)calloc(SQR(num_of_vectors), sizeof(double));
    NULL_ERROR_CHECK(P);
    int k;

    for (k = 0; k < num_of_vectors; k++){
        P[k*num_of_vectors + k] = 1;
    }

    P[i*num_of_vectors + i] = c;
    P[j*num_of_vectors + j] = c;
    P[i*num_of_vectors + j] = s;
    P[j*num_of_vectors + i] = -1*s;

    return P;
}

double** calc_eigen(double* A){
    double** result;
    double* P;
    double* eigenvalues;
    double* V;
    int* data;
    double c, s, sum_A_squared, sum_A_next_squared;
    double temp_ii, temp_ij, temp_jj, temp_ri, temp_rj;
    int max_i, max_j, r, l, i;
    int iteration = 1;

    result = (double**)malloc(sizeof(double*)*2);
    NULL_ERROR_CHECK(result);
    eigenvalues = (double*)malloc(sizeof(double)*num_of_vectors);
    NULL_ERROR_CHECK(eigenvalues);

    do{
        data = max_abs_off_diagonal_entry(A, num_of_vectors);
        max_i = data[0];
        max_j = data[1];
        P = pivot_jacobi(A, max_i, max_j);
        c = P[0];
        s = P[1];
        free(data);
        free(P);

        sum_A_squared = sum_squares_off_diagonal(A, num_of_vectors);

        temp_ii = A[max_i*num_of_vectors+max_i];
        temp_ij = A[max_i*num_of_vectors+max_j];
        temp_jj = A[max_j*num_of_vectors+max_j];
        
        for (r = 0; r < num_of_vectors; r++){
            if ((r != max_i) && (r != max_j)){
                temp_ri = A[r*num_of_vectors+max_i];
                temp_rj = A[r*num_of_vectors+max_j];

                A[r*num_of_vectors+max_i] = c*temp_ri - s*temp_rj;

                A[r*num_of_vectors+max_j] = c*temp_rj + s*temp_ri;
            }
        }

        A[max_i*num_of_vectors+max_i] = 
                SQR(c)*temp_ii + SQR(s)*temp_jj - 2*s*c*temp_ij;

        A[max_j*num_of_vectors+max_j] = 
                SQR(s)*temp_ii + SQR(c)*temp_jj + 2*s*c*temp_ij;

        A[max_i*num_of_vectors+max_j] = 0;
        A[max_j*num_of_vectors+max_i] = 0;

        if (iteration == 1){
            V = create_initial_p_matrix(max_i, max_j, c, s);
        }
        else{
            rotation_matrix_multiply_simplified(V, max_i, max_j, c, s);
        }

        sum_A_next_squared = sum_squares_off_diagonal(A, num_of_vectors);
        iteration++;
    }
    while((sum_A_squared - sum_A_next_squared > JACOBIAN_EPSILON) && 
    (iteration < JACOBIAN_MAX_ITER));

    for (i = 0; i < num_of_vectors; i++){
        eigenvalues[i] = A[num_of_vectors*i+i];
    }
    free(A);

    result[0] = eigenvalues;
    result[1] = V;
    return result;
}


double * pivot_jacobi(double * A, int max_i, int max_j){
    int sign;
    double t, c, s;
    double theta = (A[max_j*num_of_vectors+max_j] 
                    - A[max_i*num_of_vectors+max_i]) /
                    (2*A[num_of_vectors*max_i + max_j]);
    if(theta < 0){
        sign = -1;
    }
    else{
        sign = 1;
    }
    t = sign / (abs(theta)+ sqrt(SQR(theta)+1));
    c = 1 / (sqrt(SQR(t)+1));
    s = t*c;
    double * return_vals = (double*)malloc(2*sizeof(double));
    NULL_ERROR_CHECK(return_vals);
    return_vals[0] = c;
    return_vals[1] = s;
    return return_vals;
}

void rotation_matrix_multiply_simplified(double * mat, 
                                         int a, int b, 
                                         double c, double s){
    int i;
    double v_ia, v_ib;
    for (i=0; i<num_of_vectors; i++){
        v_ia = mat[num_of_vectors*i+a];
        v_ib = mat[num_of_vectors*i+b];
        mat[num_of_vectors*i+a] = c*v_ia - s*v_ib;
        mat[num_of_vectors*i+b] = s*v_ia + c*v_ib;
    }
}


double* read_file(char* file_in){
    FILE *ifp;
    int vec_size = 1;
    int flag =0;
    int vector_count=0;
    char c = 'a';
    int nch =0; 
    int i = 0;
    int size = 1000;
    char *buf = malloc(size);
    NULL_ERROR_CHECK(buf);
    int startplace;
    int vec_idx;
    char* g;
    double* points;

    ifp = NULL;
    ifp = fopen(file_in, "r");
    NULL_ERROR_CHECK(ifp);
    while((c = fgetc(ifp)) != EOF){
        if(nch >= size-1){
            size *= 2;
            buf = realloc(buf, size);
            NULL_ERROR_CHECK(buf);
        }
        if(flag==0){
            if (c == ','){
                vec_size++;
            }
        }
        if(c == '\n'){
            flag = 1;
            vector_count++;
        }
        buf[nch++] = c;
    }

    buf[nch++] = '\000';
    num_of_vectors = vector_count;
    dim = vec_size;
    fclose(ifp);
    
    startplace = 0;
    vec_idx = 0;
    points = (double*)malloc(num_of_vectors*dim*sizeof(double));
    NULL_ERROR_CHECK(points);
    while(buf[i] != '\000'){
        startplace = i;
        while(buf[i] != ',' && buf[i] != '\n' && buf[i] != '\000'){
            i++;
        }
        points[vec_idx] = (double)strtod(&buf[startplace], &g);
        if (buf[i] != '\000'){
            i++;
        }
        vec_idx++;
    }
    free(buf);
    return points;
}

void print_matrix(double* matrix, char deli, int num_rows, int num_cols){
    int i,j;
    for(i = 0; i < num_rows; i++){
        for(j = 0; j < num_cols; j++){
            printf("%.4f", matrix[num_rows*i + j]);
            if(j!=(num_cols-1)){
                printf("%c", deli);
            }
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]){
    int i, j;
    double* points;
    double** data;
    double* weighted_matrix;
    double* diag_deg_matrix;
    double* lnorm;
    if(argc != 3){
        printf("Invalid Input!");
        exit(1);
    }

    if(!(strcmp(argv[1],"jacobi") != 0 || strcmp(argv[1], "wam") != 0 || 
       strcmp(argv[1],"ddg") != 0 || strcmp(argv[1], "lnorm") != 0)){
        printf("Invalid Input!");
        exit(1);
    }
    points = read_file(argv[2]);
    if(strcmp(argv[1],"jacobi") == 0){
        data = calc_eigen(calc_lnorm_matrix(points));
        for(i = 0; i<num_of_vectors; i++){
            printf("%.4f", data[0][i]);
            if(i!=(num_of_vectors-1)){
                    printf("%c", ',');
                }
        }
        printf("\n");
        print_matrix(data[1], ' ', num_of_vectors, num_of_vectors);
        free(data);
    }

    if(strcmp(argv[1], "wam") == 0){
        weighted_matrix = calc_weighted_matrix(points);
        print_matrix(weighted_matrix, ',', num_of_vectors, num_of_vectors);
    }

    if(strcmp(argv[1],"ddg") == 0){
        diag_deg_matrix = 
        calc_diagonal_deg_matrix(calc_weighted_matrix(points));
        for(i=0; i<num_of_vectors; i++){
            for(j=0; j<num_of_vectors; j++){
                if(i==j){
                    printf("%.4f", diag_deg_matrix[i]);
                }
                else{
                    printf("0.0000");
                }
                if(j!=(num_of_vectors-1)){
                    printf("%c", ',');
                }
            }
            printf("\n");
        }
    }

    if(strcmp(argv[1],"lnorm") == 0){
        lnorm = calc_lnorm_matrix(points);
        print_matrix(lnorm, ',', num_of_vectors, num_of_vectors);
    }
}


