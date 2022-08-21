#include "spkmeans.h"

static PyMethodDef capiMethods[] = {
    {"fit", 
        (PyCFunction) fit_capi,
        METH_VARARGS,
        PyDoc_STR("A list of K centroid objects")
    },
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods
};


PyMODINIT_FUNC PyInit_mykmeanssp(void){
    PyObject* m;
    m = PyModule_Create(&moduledef);
    
    if (!m){
        return NULL;
    }
    return m;
}


int main(int argc, char** argv){
    return 0;
}
