
import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve


import SimplexAlgorithm


t = Symbol('t')


print("пример 2.60")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[6, 0, 4+t, 12-t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, -1, 1, 0, 0],
     [-1, 3, 0, 1, 0],
     [-1/2, 2, 0, 0, 1]]#коэффициенты из ограничений-равенств 
b_eq=[28, 20, 24] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)




print("пример 2.69")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[-(1-t), 4-t, -(2-t), 2-t, -(3-2*t)] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 1, 1, 1, 0],
     [-2, 1, -1, 0, 1]]#коэффициенты из ограничений-равенств 
b_eq=[2, 1] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)



print("пример 2.70")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[6, 0, -(4+t), 12-t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, -1, 1, 0, 0],
     [-1, 3, 0, 1, 0],
     [-1/2, 2, 0, 0, 1]]#коэффициенты из ограничений-равенств 
b_eq=[28, 20, 24] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)

print("пример 2.71")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[2, 5, 0, 0, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, -2, 1, 0, 0],
     [-2, 1, 0, 1, 0],
     [3, -1, 0, 0, 1]]#коэффициенты из ограничений-равенств 
b_eq=[4+2*t, 6+t, 8-3*t] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)
 
