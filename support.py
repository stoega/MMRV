#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:51:28 2024

@author: Kurt
"""

# In[0]

import numpy as numpy
import control as control

# In[1]

def n_I(n):
    """.

    Parameters
    ----------
    n : int, Anzahl der Integratoren

    Returns
    -------
    TYPE control.ss. ABCD Darstellung
    """
    A = numpy.zeros((n, n), dtype=float)
    b = numpy.zeros((n, 1), dtype=float)
    c = numpy.zeros((1, n), dtype=float)
    d = numpy.zeros((1, 1), dtype=float)
    for i in range(n-1):
        A[i, i+1] = 1.0
    b[n-1, 0] = 1.0
    c[0, n-1] = 1.0
    return control.ss(A, b, c, d)


def Equ(A, b, N):
    """.

    Parameters
    ----------
    A : n x n float array
    b : n x n float array
    N : int

    Returns x x N float array [A^(N-1)b Ab ... b]
    """
    Res = b
    for i in range(1, N):
        Res = numpy.hstack((numpy.matmul(A, Res[:, 0:1]), Res))
    return Res


def Sol(A, b, U, x0):
    """.

    Parameters
    ----------
    A : n x n float array
    b : n x 1 float array
    U : 1 x N float array
    x0 : n x 1 float array

    Returns
    -------
    Res : n x (N+1) float array [x_0 x_1 ... X_N]
    """
    Res = x0
    for i, u in enumerate(U):
        Res = numpy.hstack((Res, numpy.matmul(A, Res[:, i:(i+1)]) + b*u))
    return Res

