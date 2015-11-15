import numpy as np
from itertools import product
from scipy.linalg import solve_banded

u = 500
l = 500
M = 500000 # A is a MxM  matrix with 1 non-zero upper diagonal and 2 non-zero lower diagonals

def getA(ab, i,j):
    if i - j <= l and j - i <= u:
        return ab[u+i-j,j]
    else:
        return 0

def to_dense(ab):
    A = np.empty((M,M))
    for i,j in product(range(M), range(M)):
        A[i,j] = getA(ab, i,j)
    return A

print('creating matrices')
ab = np.random.randn(u+l+1, M)
b = np.random.randn(M)

print('solving system')
x = solve_banded((l,u), ab, b)
print(x)
