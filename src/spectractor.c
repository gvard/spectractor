#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <math.h>

static PyObject *SpectractorError;

static PyObject *
write_fds(PyObject *self, PyObject *args)
{
    PyArrayObject *py_fds;
    const char *fdsfile;
    double *fds;
    FILE *ffil;
    static unsigned int nx, ny, iwl;
    unsigned short i, j;
    static double wl;

    if(!PyArg_ParseTuple(args, "O!s", &PyArray_Type, &py_fds, &fdsfile))
        return NULL;
    if(NULL == py_fds) return NULL;

    fds = (double *)PyArray_DATA(py_fds);
    nx = PyArray_DIM(py_fds, 0);
    ny = PyArray_DIM(py_fds, 1);

    if ((ffil = fopen(fdsfile, "w")) == NULL){
        printf("Error for open file %s", fdsfile);
        return NULL;
    }
    for (i=0; i < ny; i++){
        for (j=0; j < nx; j++){
            wl = round(fds[i*nx+j] * 10000.);
            iwl = (int) wl;
            fwrite(&iwl, sizeof(int), 1, ffil);
        }
    }
    fclose(ffil);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
write_bin(PyObject *self, PyObject *args)
{
    PyArrayObject *py_bin;
    const char *fdsfile, *objname;
    double *dat;
    FILE *ffil;
    static unsigned int nx, ny;
    short iwl;
    unsigned short i, j;
    static double wl;

    if(!PyArg_ParseTuple(args, "O!ss", &PyArray_Type, &py_bin,
        &fdsfile, &objname))
        return NULL;
    if(NULL == py_bin) return NULL;

    dat = (double *)PyArray_DATA(py_bin);
    nx = PyArray_DIM(py_bin, 0);
    ny = PyArray_DIM(py_bin, 1);

    if ((ffil = fopen(fdsfile, "w")) == NULL){
        printf("Error for open file %s\n", fdsfile);
        return NULL;
    }
    fwrite(objname, 10, 1, ffil);
    fwrite(&nx, sizeof(short), 1, ffil);
    fwrite(&ny, sizeof(short), 1, ffil);
    for (i=0; i < ny; i++){
        for (j=0; j < nx; j++){
            wl = round(dat[i*nx+j]);
            iwl = (short) wl;
            fwrite(&iwl, sizeof(short), 1, ffil);
        }
    }
    fclose(ffil);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
read_fds(PyObject *self, PyObject *args)
{
    const char *binfile;
    const unsigned short nx, ny;
    FILE *ffil;
    unsigned short i, j;
    PyArrayObject *py_bin;
    npy_intp dims[2];
    double leng;
    double *dat;
    unsigned int iwl;

    if(!PyArg_ParseTuple(args, "sHH", &binfile, &nx, &ny))
        return NULL;

    dims[0] = nx;
    dims[1] = ny;

    if ((ffil = fopen(binfile, "r")) == NULL){
        printf("Error for open file %s\n", binfile);
        return NULL;
    }
    py_bin = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    dat = (double *)PyArray_DATA(py_bin);
    for (i=0; i < ny; i++){
        for (j=0; j < nx; j++){
            fread(&iwl, sizeof(int), 1, ffil);
            leng = (double) iwl;
            dat[i*nx+j] = leng/10000;
        }
    }
    fclose(ffil);
    return PyArray_Return(py_bin);
}


static PyMethodDef SpectractorMethods[] = {
    {"write_fds", write_fds, METH_VARARGS,
        "Create fds file from numpy array of floats."},
    {"write_bin", write_bin, METH_VARARGS,
        "Create bin file from numpy array of floats."},
    {"read_fds", read_fds, METH_VARARGS,
        "Read fds file to numpy array."},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_spectractor(void)
{
    PyObject *m;
    m = Py_InitModule("_spectractor", SpectractorMethods);
    if (m == NULL)
        return;
    SpectractorError = PyErr_NewException("_spectractor.error", NULL, NULL);
    Py_INCREF(SpectractorError);
    PyModule_AddObject(m, "error", SpectractorError);
    import_array();
}
