import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve


import SimplexAlgorithm


t = Symbol('t')

print("пример 1")
#ab=[-float("inf"),  float("inf")]
ab=[0,  float("inf")] #промежуток t

z=[3-6*t, 2-2*t, 5+5*t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 2, 1, 1]] #коэффициенты из ограничений-равенств 
b_eq=[40] 
a_ub=[[3, 0, 2, 0], [1, 4, 0, 0]] #коэффициенты из ограничений-неравенств
b_ub=[60, 30]


SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
 
print("пример 2")
ab=[0,  float("inf")]
z=[3, 2, 5] #тут хрвняться выражения перед иксами в функции
a_eq=[] #коэффициенты из ограничений-равенств 
b_eq=[] 
a_ub=[[1, 2, 1], [3, 0, 2], [1, 4, 0]] #коэффициенты из ограничений-неравенств
b_ub=[40-t, 60+2*t, 30-7*t]

SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
