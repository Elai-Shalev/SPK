#include "spkmeans.h"


/* Global Variables */
double const JACOBIAN_EPSILON = 0.00001;
int const JACOBIAN_MAX_ITER = 100;

double* dimension_reduction_spk(double* points){
    double* l_norm; 
    double** data;
    double* eigen_values;
    double* eigen_vectors;
    double** eigen_vectors_matrix;
    double* normalized_eigen_vectors;
    l_norm = calc_lnorm_matrix(points);
    data = calc_eigen(l_norm);
    eigen_values = data[0];
    eigen_vectors = data[1];
    free(data);

    eigen_vectors_matrix = convert_double_array_to_matrix(
        eigen_vectors,num_of_vectors,num_of_vectors);

    transpose_square_matrix_inplace(eigen_vectors_matrix, num_of_vectors);
    sort_eigen_v(eigen_values, eigen_vectors_matrix);
    transpose_square_matrix_inplace(eigen_vectors_matrix, num_of_vectors);

    if (K == 0){
        K = determine_k(eigen_values);
    }
        
    normalize_first_k_vectors(eigen_vectors_matrix, K);

    normalized_eigen_vectors = convert_matrix_to_double_array(eigen_vectors_matrix, num_of_vectors, K);

    free(l_norm);
    free(eigen_values);
    free(eigen_vectors);
    free(eigen_vectors_matrix);

    return normalized_eigen_vectors; /* T in the algorithm */
}

void transpose_square_matrix_inplace(double** matrix, int size){
    int i, j;
    double temp;

    for(i = 0; i < size; i++){
        for(j = i; j < size; j++){
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

void normalize_first_k_vectors(double** eigen_vectors_matrix, int k){
    int i,j;
    double curr_norma = 0;
    for(i = 0; i < num_of_vectors; i++){
        curr_norma = calc_norma(eigen_vectors_matrix[i], k);
        for(j=0; j < k; j++){
            if (curr_norma != 0){
                eigen_vectors_matrix[i][j] /= curr_norma;
            }
        }
    }
}

double calc_norma(double* vector, int size){
    double sum = 0;
    int i;
    for(i=0; i<size; i++){
        sum += SQR(vector[i]);
    }
    return pow(sum, 0.5);
}

int determine_k(double* sorted_eigen_values){
    int i; 
    int k = 0;
    double gap, curr_gap;
    gap = 0;
    for(i=0; i<(num_of_vectors/2); i++){
        curr_gap = ABS(sorted_eigen_values[i] - sorted_eigen_values[i+1]);
        if(curr_gap > gap){
            gap = curr_gap;
            k=i;
        }
    }
    return k+1;
}

double** convert_double_array_to_matrix(double* array, 
                                        int num_rows,
                                        int num_cols){
    int i;
    double** mat = (double**)malloc(sizeof(double*)*num_rows);
    NULL_ERROR_CHECK(mat);
    for(i=0; i<num_rows; i++){
        mat[i] = &array[i*num_cols];
    }
    return mat;
}

double* convert_matrix_to_double_array(double** matrix, 
                                        int new_num_rows, 
                                        int new_num_cols){
    double* array;
    int i, j;
    array = (double*)malloc(sizeof(double)*new_num_cols*new_num_rows);
    NULL_ERROR_CHECK(array);

    for(i = 0; i < new_num_rows; i++){
        for(j = 0; j < new_num_cols; j++){
            array[i*new_num_cols+j] = matrix[i][j];
        }
    }
    return array;
}

void sort_eigen_v(double* eigen_values, double** eigen_vectors){
    int i, j;
    double temp_eigenvalue;
    double* temp_eigenvector;

    for(i = 0; i < num_of_vectors; i++){
        for (j = i; j < num_of_vectors; j++){
            if (eigen_values[j] > eigen_values[i]){
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
                weighted_adj_matrix[num_of_vectors*i+j] = 0;
            }
            else{
                norm_i_j = euclidian_dist(points,i,j);
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

    free(weighted_matrix);
    free(diagonal_deg_matrix);
    return lnorm_matrix;
}


int* max_abs_off_diagonal_entry(double* mat, int size){
    int i=0;
    int j=0;
    int* data = (int*)malloc(sizeof(int)*2);
    double currmax, currnum;
    NULL_ERROR_CHECK(data);
    

    data[0] = 0;
    data[1] = 1;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            currmax = ABS(mat[(data[0])*size+(data[1])]);
            currnum = ABS(mat[i*size+j]);
            if ((i != j) && (currmax < currnum)){
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
                sum += SQR(mat[i*size+j]);
            }
        }
    } 

    return sum;
}

double* create_initial_p_matrix(int i, int j, double c, double s){
    int k=0;
    double* P = (double*)calloc(SQR(num_of_vectors), sizeof(double));
    NULL_ERROR_CHECK(P);
    

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
    int max_i, max_j;
    int r=0;
    int i=0;
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
                A[max_i*num_of_vectors+r] = c*temp_ri - s*temp_rj;

                A[r*num_of_vectors+max_j] = c*temp_rj + s*temp_ri;
                A[max_j*num_of_vectors+r] = c*temp_rj + s*temp_ri;
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
    while((ABS(sum_A_squared - sum_A_next_squared) > JACOBIAN_EPSILON) && 
    (iteration < JACOBIAN_MAX_ITER));

    for (i = 0; i < num_of_vectors; i++){
        eigenvalues[i] = A[num_of_vectors*i+i];
    }

    result[0] = eigenvalues;
    result[1] = V;
    return result;
}


double * pivot_jacobi(double * A, int max_i, int max_j){
    int sign;
    double * return_vals;
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
    t = sign / (ABS(theta)+ sqrt(SQR(theta)+1));
    c = 1 / (sqrt(SQR(t)+1));
    s = t*c;
    return_vals = (double*)malloc(2*sizeof(double));
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
    int startplace=0;
    int flag =0;
    int vector_count=0;
    char c = 'a';
    int nch =0; 
    int i = 0;
    int vec_idx;
    char* g;
    double* points;
    int size = 1000;
    char *buf = malloc(size);
    NULL_ERROR_CHECK(buf);
    

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
            printf("%.4f", matrix[num_cols*i + j]);
            if(j!=(num_cols-1)){
                printf("%c", deli);
            }
        }
        printf("\n");
    }
}

void print_matrix_doublestar(double** matrix, char deli, int num_rows, int num_cols){
    int i,j;
    for(i = 0; i < num_rows; i++){
        for(j = 0; j < num_cols; j++){
            printf("%.4f", matrix[i][j]);
            if(j!=(num_cols-1)){
                printf("%c", deli);
            }
        }
        printf("\n");
    }
}

void print_double_array(double* arr, char deli, int length){
    int i;
    for(i = 0; i < length; i++){
        printf("%.4f", arr[i]);
        if(i!=(length-1)){
            printf("%c", deli);
        }
        printf("\n");
    }
}

void wam_c(double* points){
    double* weighted_matrix;
    weighted_matrix = calc_weighted_matrix(points);
    print_matrix(weighted_matrix, ',', num_of_vectors, num_of_vectors);
    free(weighted_matrix);
}

void ddg_c(double* points){
    double* diag_deg_matrix;
    double* weighted_matrix;
    int i,j;
    weighted_matrix = calc_weighted_matrix(points);
    diag_deg_matrix = calc_diagonal_deg_matrix(weighted_matrix);
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
    free(weighted_matrix);
    free(diag_deg_matrix);
}

void lnorm_c(double* points){
    double* lnorm;
    lnorm = calc_lnorm_matrix(points);
    print_matrix(lnorm, ',', num_of_vectors, num_of_vectors);
    free(lnorm);
}

void jacobi_c(double* lnorm_input){
    double** data;
    int i;
    data = calc_eigen(lnorm_input);
    for(i = 0; i<num_of_vectors; i++){
        printf("%.4f", data[0][i]);
        if(i!=(num_of_vectors-1)){
                printf("%c", ',');
            }
    }
    printf("\n");
    print_matrix(data[1], ' ', num_of_vectors, num_of_vectors);
    free(data[0]);
    free(data[1]);
    free(data);
}

int main(int argc, char* argv[]){
    double* points;
    double* lnorm;
    K = 0;
    if(argc != 3){
        printf("Invalid Input!");
        exit(1);
    }

    /* DELETE SPK */
    if(!(strcmp(argv[1],"jacobi") != 0 || strcmp(argv[1], "wam") != 0 || 
        strcmp(argv[1],"ddg") != 0 || strcmp(argv[1], "lnorm") != 0)){
        printf("Invalid Input!");
        exit(1);
    }
    points = read_file(argv[2]);
    if(strcmp(argv[1],"jacobi") == 0){
        jacobi_c(points);
    }
    else if(strcmp(argv[1], "wam") == 0){
        wam_c(points);
    }
    else if(strcmp(argv[1],"ddg") == 0){
        ddg_c(points);
    }
    else if(strcmp(argv[1],"lnorm") == 0){
       lnorm_c(points);
    }
    else{
        printf("Invalid Input!");
        exit(1);
    }
    free(points);
    return 0;
}


double square_distance(Vector* vec1, Vector* vec2){
    double sum_of_squares = 0;
    int i;

    for (i = 0; i < dim; i++){
        sum_of_squares += SQR(((vec1->coordinate)[i]) - ((vec2->coordinate)[i]));  
    }
    
    return sum_of_squares;
}


void assign_to_nearest_cluster(Vector* vec, Vector** centroids){
    double min_distance;
    int min_cluster = -1;
    double distance_from_cluster;
    int i;
    
    for(i = 0; i < K; i++){
        distance_from_cluster = square_distance(vec, centroids[i]);
        if (i == 0 || (distance_from_cluster < min_distance)){
            min_distance = distance_from_cluster;
            min_cluster = i;
        }
    }

    (vec->cluster) = min_cluster;
    return;
}


Vector* new_zero_vector(){
    Vector* v = (Vector*)malloc(sizeof(Vector));
    double* vec_coordinates;

    if(v==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    vec_coordinates = (double*)calloc(dim, sizeof(double));
    if(vec_coordinates == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    v->coordinate = vec_coordinates;
    v->cluster = -1;
    return v;
}


Vector* new_deep_copy_vector(Vector* vec){
    int i;
    Vector* new_vec = new_zero_vector();

    for(i = 0; i < dim; i++){
        (new_vec->coordinate[i]) = (vec->coordinate[i]);
    }
    new_vec->cluster = vec->cluster;

    return new_vec;
}


void vectoric_sum(Vector* vec1, Vector* vec2){
    int i;

    for(i = 0; i < dim; i++){
        (vec1->coordinate[i]) += (vec2->coordinate[i]);
    }

    return;
}


void divide_vector(Vector* vec, int divisor){
    int i;

    for(i = 0; i < dim; i++){
        (vec->coordinate[i]) /= ((float)divisor);
    }

    return;
}


int update_centroids(Vector** vector_list, Vector** centroids){
    int i, j;
    int returnValue = 1;
    Vector** cluster_sums;
    int* cluster_sizes;

    cluster_sums = (Vector**)malloc(K*sizeof(Vector*));
    if(cluster_sums == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    cluster_sizes = (int*)malloc(K*sizeof(int));
    if(cluster_sizes == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    for(i = 0; i < K; i++){
        cluster_sums[i] = new_zero_vector();
        cluster_sizes[i] = 0;
    }

    for(i = 0; i < num_of_vectors; i++){
        vectoric_sum(cluster_sums[vector_list[i]->cluster], vector_list[i]);
        cluster_sizes[vector_list[i]->cluster]++;
    }

    for(i = 0; i < K; i++){
        divide_vector(cluster_sums[i], cluster_sizes[i]);
        if (sqrt(square_distance(centroids[i], cluster_sums[i])) >= EPSILON){
            returnValue = 0;
        }
        for(j = 0; j < dim; j++){
            (centroids[i]->coordinate[j]) = (cluster_sums[i]->coordinate[j]);
        }
    }

    free(cluster_sums);
    free(cluster_sizes);
    return returnValue;
}


void run_kmeans(Vector** vectors, Vector** centroids){
    int i;
    int is_converged = 0;
    int iter_count = 0;

    while((is_converged == 0) && (iter_count < MAX_ITER)){
        for(i = 0; i < num_of_vectors; i++){
            assign_to_nearest_cluster(vectors[i], centroids);
        }

        is_converged = update_centroids(vectors, centroids);
        iter_count++;
    }

    return;
}


Vector** double_array_to_vector_array(double* vector_array){
    int vec_idx;
    Vector* v;
    Vector** vectors;
    double* vec_coordinates;
    int vector_array_idx = 0;
    int dim_idx;

    vectors = (Vector**)malloc(num_of_vectors*(sizeof(Vector*)));
    if(vectors == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    for(vec_idx = 0; vec_idx < num_of_vectors; vec_idx++){
        v = (Vector*)malloc(sizeof(Vector));
        if(v == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }

        vec_coordinates = (double*)malloc(dim*sizeof(double));
        if(vec_coordinates == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }

        for(dim_idx = 0; dim_idx < dim; dim_idx++){
            vec_coordinates[dim_idx] = vector_array[vector_array_idx];
            vector_array_idx++;
        }
        
        v -> coordinate = vec_coordinates;
        v -> cluster = -1;
        vectors[vec_idx] = v;
    }

    free(vector_array);
    return vectors;
}


Vector** fit_c(double* vector_array, double* centroid_array){
    Vector** vector_list;
    Vector** centroid_list;

    vector_list = double_array_to_vector_array(vector_array);
    centroid_list = double_array_to_vector_array(centroid_array);

    run_kmeans(vector_list, centroid_list);
    free(vector_list);

    return centroid_list;
}
