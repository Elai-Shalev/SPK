#define PY_SSIZE_T_CLEAN
// include "Python.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


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
    double norm_i_j;
    double* weighted_adj_matrix = (double*)malloc(sizeof(double)*SQR(num_of_vectors));
    for(i = 0; i<num_of_vectors; i++){
        for(j=i; j<num_of_vectors; j++){
            if(i==j){
                weighted_adj_matrix[num_of_vectors*i+j] = 0;
            }
            else{
                double norm_i_j = norm2(points,i,j);
                weighted_adj_matrix[num_of_vectors*i+j] = exp(-0.5*norm_i_j);
                weighted_adj_matrix[num_of_vectors*j+i] = exp(-0.5*norm_i_j);
            }
        }
    }
    return weighted_adj_matrix;

}
double norm2(double* points, int i, int j){
    int k;
    int sum_diff;
    for(k=0; k<dim; k++){
        sum_diff += pow((points[dim*i+k] - points[dim*j+k]),2);
    }
    return sqrt(sum_diff);
}

double* calc_diagonal_deg_matrix(double* mat){
    int z,i;
    int sum_row;
    double* diagonal_deg_matrix = (double*)malloc(sizeof(double)*num_of_vectors);
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
    
    for (i = 0; i < num_of_vectors; i++){
        for (j = 0; j < num_of_vectors; j++){
            lnorm_matrix[i*num_of_vectors+j] = -1* pow(diagonal_deg_matrix[i], -0.5)*
                                    weighted_matrix[i*num_of_vectors+j]*
                                    pow(diagonal_deg_matrix[j], -0.5);
            
            if (i==j){
                lnorm_matrix[i*num_of_vectors+j] += 1;
            }
        }
    }

    return lnorm_matrix;
}


int* max_abs_off_diagonal_entry(double* mat, int size){
    int* data = (int*)malloc(sizeof(int)*2);
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

double* create_initial_p_matrix(i, j, c, s){
    double* P = (double*)calloc(sizeof(double)*SQR(num_of_vectors));
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
    double* V;
    int* data;
    double c, s, sum_A_squared, sum_A_next_squared;
    int max_i, max_j, r, l;
    int iteration = 1;

    result = (double**)malloc(sizeof(double*)*2);
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
        
        for (r = 0; r < num_of_vectors; r++){
            if ((r != max_i) && (r != max_j)){
                A[r*num_of_vectors+max_i] = c*A[r*num_of_vectors+max_i]-
                                                s*A[r*num_of_vectors+max_j];

                A[r*num_of_vectors+max_j] = c*A[r*num_of_vectors+max_j]+
                                                s*A[r*num_of_vectors+max_i];
            }
        }
        A[max_i*num_of_vectors+max_i] = 
                SQR(c)*A[max_i*num_of_vectors+max_i]+
                SQR(s)*A[max_j*num_of_vectors+max_j]-
                2*s*c*A[max_i*num_of_vectors+max_j];

        A[max_j*num_of_vectors+max_j] = 
                SQR(s)*A[max_i*num_of_vectors+max_i]+
                SQR(c)*A[max_j*num_of_vectors+max_j]+
                2*s*c*A[max_i*num_of_vectors+max_j];

        A[max_i*num_of_vectors+max_j] = 0;

        if (iteration == 1){
            V = create_initial_p_matrix(max_i, max_j, c, s);
        }
        else{
            rotation_matrix_multiply_simplified(V, max_i, max_j, c, s);
        }

        sum_A_next_squared = sum_squares_off_diagonal(A, num_of_vectors);
    }
    while((sum_A_squared - sum_A_next_squared > JACOBIAN_EPSILON) || 
    (iteration < JACOBIAN_MAX_ITER));

    result[0] = A;
    result[1] = V;
}


double * pivot_jacobi(double * A, int max_i, int max_j){
    int sign;
    double t, c, s;
    double theta = (A[max_j*num_of_vectors+max_j] - A[max_i*num_of_vectors+max_i]) /
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
    return_vals[0] = c;
    return_vals[1] = s;
    return return_vals;
}

void rotation_matrix_multiply_simplified(double * mat, int a, int b, double c, double s){
    int i;
    double v_ia, v_ib;
    for (i=0; i<num_of_vectors; i++){
        v_ia = mat[num_of_vectors*i+a];
        v_ib = mat[num_of_vectors*i+b];
        mat[num_of_vectors*i+a] = c*v_ia - s*v_ib;
        mat[num_of_vectors*i+b] = s*v_ia + c*v_ib;
    }
}

int main(int argc, char* argv[]){
    int i, j;
    if(argc != 2){
        printf("Invalid Input!");
        exit(1);
    }

    enum command{
        wam, ddg, lnorm, jacobi,
    };
    enum command choice;
    if(argv[1] == "wam"){
        choice = wam;
    }
    if(argv[1] == "ddg"){
        choice = ddg;
    }
    if(argv[1] == "lnorm"){
        choice = lnorm;
    }
    if(argv[1] == "jacobi"){
        choice = jacobi;
    }
    else{
        printf("Invalid Input!");
        exit(1);
    }
    
    //File Read
    FILE* fp;
    fp = fopen(argv[2], "r"); 
    assert(fp);

    if(choice != jacobi){
        double* input_matrix = (double*)malloc(num_of_vectors*dim*sizeof(double));
        assert(input_matrix);
        int res = fread(input_matrix, sizeof(double), num_of_vectors*dim, fp);
        if(res != num_of_vectors*dim){
            exit(1);
        }
    }
    if(choice == jacobi){
        double* symmetrical_matrix = (double*)malloc(SQR(num_of_vectors)*sizeof(double));
        assert(symmetrical_matrix);
        int res = fread(symmetrical_matrix, sizeof(double), SQR(num_of_vectors), fp);
        if(res != SQR(num_of_vectors)){
            exit(1);
        }

        double** eigen_result = calc_eigen(symmetrical_matrix);
        for(i = 0; i<num_of_vectors; i++){
            printf("%d, ", eigen_result[0][num_of_vectors*i+i]);
        }
        printf("\n");
        for(i = 0; i < num_of_vectors; i++){
            for(j = 0; j < num_of_vectors; j++){
                printf("%d ", eigen_result[1][num_of_vectors*i + j]);
            }
            printf("\n");
    }
    
    fclose(fp);
    return 0;
    }
}

