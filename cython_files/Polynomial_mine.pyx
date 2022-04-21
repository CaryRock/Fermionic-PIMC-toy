import numpy as np
#cimport numpy as np

cpdef double GenPolynomialValue_c(int order, double[:] vals, double x):
    cdef int i = 0
    cdef double summmation 
    summation = 0.0
    for i in range(0, order):
        summation += vals[i] * x**i

    return summation
