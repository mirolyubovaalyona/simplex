import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve

import SimplexAlgorithm


t = Symbol('t')

#ab=[-float("inf"),  float("inf")]
ab=[0,  float("inf")] #промежуток t

z=[3-6*t, 2-2*t, 5+5*t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 2, 1, 1]] #коэффициенты из ограничений-равенств 
b_eq=[40] 
a_ub=[[3, 0, 2, 0], [1, 4, 0, 0]] #коэффициенты из ограничений-неравенств
b_ub=[60, 30]


SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
#print(Res)

#var('t')
z=[7*t+10, 2*t-9]
for i in range(len(z)):
    z[i]=z[i]+t>=0
print(z)
r=solve(z, t)
print(str(r), str(r)[1])

import random

z = 4
x = random.sample(range(50), z)
y = random.sample(range(50), z)

for i in x:
    for u in y:        
        print('{:6d}'.format(i*u), end='| ')   
    print('')