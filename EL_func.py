from numpy import *
from scipy import linalg

def EL_func (t, x):
    M = zeros((2, 2))
    R = zeros(2)
    M[0, 0] = 0.0517835709769978
    M[0, 1] = 0.00134192639324236
    M[1, 1] = 0.00134192639324236
    R[0] = -1.7720667778888*cos(x[0])
    R[1] = 0
    return hstack((x[2:], linalg.solve(M, R, assume_a='pos')))
