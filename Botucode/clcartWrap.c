#include <Python.h>
//#include <cstdio>
extern void cluster_c();

PyObject* clcart_cluster_c(PyObject* self, PyObject* args)
{
        cluster_c();
	return Py_BuildValue("");
}


static PyMethodDef clcartmethods[] = {
	{"cluster_c", clcart_cluster_c, METH_VARARGS},
	{NULL},
};

static struct PyModuleDef clcartmodule =
{
    PyModuleDef_HEAD_INIT,
    "clcartmodule",
    "",
    -1,
    clcartmethods
};

PyMODINIT_FUNC PyInit_clcartmodule(void)
{
	return PyModule_Create(&clcartmodule);
}
