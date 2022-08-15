double* dimension_reduction_spk(double* points);
double* transpose_matrix_copy(double** matrix, int rows, int cols);
void transpose_square_matrix_inplace(double** matrix, int size);
void normalize_first_k_vectors(double** eigen_vectors_matrix, int k);
double calc_norma(double* vector, int size);
int determine_k(double* sorted_eigen_values);
double** convert_double_array_to_matrix(double* array, 
                                        int size_row, 
                                        int num_rows);
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
