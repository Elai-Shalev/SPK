#include "spkmeans.h"

static PyObject* kmpp_capi(PyObject *self, PyObject *args);
static PyObject* wam_capi(PyObject *self, PyObject *args);
static PyObject* ddg_capi(PyObject *self, PyObject *args);
static PyObject* lnorm_capi(PyObject *self, PyObject *args);
static PyObject* jacobi_capi(PyObject *self, PyObject *args);
static PyObject* spk_capi(PyObject *self, PyObject *args);

double* python_list_to_c_array(PyObject* float_list){
    double* double_arr;
    int pr_length;
    int index;
    PyObject *item;

    pr_length = PyObject_Length(float_list);
    if (pr_length < 0)
        return NULL;

    double_arr = (double *) malloc(sizeof(double *) * pr_length);
    if(double_arr == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    if (double_arr == NULL)
        return NULL;
    for (index = 0; index < pr_length; index++) {
        item = PyList_GetItem(float_list, index);
        if (!PyFloat_Check(item))
            double_arr[index] = 0.0;
        double_arr[index] = PyFloat_AsDouble(item);
    }

    return double_arr;
}


PyObject* c_array_to_python_list(double* float_list){
    PyObject* python_list = PyList_New(K);
    PyObject* temp_list;
    int i, j;
    PyObject* temp_float;

    for(i = 0; i < K; i++){
        temp_list = PyList_New(dim);
        for(j = 0; j < dim; j++){
            temp_float = PyFloat_FromDouble(float_list[i*dim + j]);
            PyList_SetItem(temp_list, j, temp_float);
        }
        PyList_SetItem(python_list, i, temp_list);
    }
    free(float_list);
    return python_list;
}


static PyObject* kmpp_capi(PyObject *self, PyObject *args)
{
    PyObject* vector_float_list;
    PyObject* centroid_float_list;
    PyObject* python_list_result;
    double* vector_list;
    double* centroid_list;
    double* centroid_flattened_list;
    Vector** k_means_result;
    int i, j, idx;

    if (!PyArg_ParseTuple(args, "OOiiiiO", &vector_float_list, &centroid_float_list, &num_of_vectors, &dim, &K, &MAX_ITER, &EPSILON)){
        return NULL;
    }

    vector_list = python_list_to_c_array(vector_float_list);
    centroid_list = python_list_to_c_array(centroid_float_list);
    k_means_result = fit_c(vector_list, centroid_list);

    centroid_flattened_list = (double *) malloc(sizeof(double *) * K * dim);
    if(centroid_float_list == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    idx = 0;
    for(i = 0; i < K; i++){
        for(j = 0; j < dim; j++){
            centroid_flattened_list[idx] = (k_means_result[i] -> coordinate)[j];
            idx++;
        }
    }

    python_list_result = c_array_to_python_list(centroid_flattened_list);

    free(k_means_result);
    
    return Py_BuildValue("O", python_list_result);
}

static PyObject* wam_capi(PyObject *self, PyObject *args)
{
    PyObject* temp;
    PyObject* python_filename;
    double* points;
    
    if (!PyArg_ParseTuple(args, "si", &python_filename, &K)){
        return NULL;
    }

    points = read_file(python_filename);
    wam_c(points);
    return temp;
}

static PyObject* ddg_capi(PyObject *self, PyObject *args)
{
    PyObject* temp;
    PyObject* python_filename;
    double* points;
    
    if (!PyArg_ParseTuple(args, "si", &python_filename, &K)){
        return NULL;
    }

    points = read_file(python_filename);
    ddg_c(points);
    return temp;
}

static PyObject* lnorm_capi(PyObject *self, PyObject *args)
{
    PyObject* temp;
    PyObject* python_filename;
    double* points;
    
    if (!PyArg_ParseTuple(args, "si", &python_filename, &K)){
        return NULL;
    }

    points = read_file(python_filename);
    lnorm_c(points);
    return temp;
}

static PyObject* jacobi_capi(PyObject *self, PyObject *args)
{
    PyObject* temp;
    PyObject* python_filename;
    double* points;
    
    if (!PyArg_ParseTuple(args, "si", &python_filename, &K)){
        return NULL;
    }

    points = read_file(python_filename);
    jacobi_c(points);
    return temp;
}

static PyObject* spk_capi(PyObject *self, PyObject *args)
{
    PyObject* temp;
    return temp;
}

static PyMethodDef capiMethods[] = {
    {"fit", 
        (PyCFunction) kmpp_capi,
        METH_VARARGS,
        PyDoc_STR("A list of K centroid objects")
    },
    {"wam", 
        (PyCFunction) wam_capi,
        METH_VARARGS,
        PyDoc_STR("Runs WAM")
    },
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    capiMethods
};


PyMODINIT_FUNC PyInit_spkmeansmodule(void){
    PyObject* m;
    m = PyModule_Create(&moduledef);
    
    if (!m){
        return NULL;
    }
    return m;
}
