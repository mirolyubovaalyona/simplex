#починить minus_b

import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve


import SimplexAlgorithm


t = Symbol('t')


print("пример 2.59")
ab=[0,  10] #промежуток t
z=[2+t, 13-t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[] #коэффициенты из ограничений-равенств 
b_eq=[] 
a_ub=[[4, 1, 0], [2, 2, 0], [6, 0, 3]] #коэффициенты из ограничений-неравенств
b_ub=[16, 22, 36]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)
 
print("пример 2.64")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[2, 3+4*t, 0, 0, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 1, 1, 0, 0], 
      [1, -1, 0, 1, 0],
      [-1, 1, 0, 0, 1]] #коэффициенты из ограничений-равенств 
b_eq=[12, 10, 6] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

# Работвет
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)


print("пример 2.66")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[3, -2, 5, 0, -4] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 1, 1, 0, 0], 
      [2, -1, 0, 1, 0],
      [-2, 2, 0, 0, 1]] #коэффициенты из ограничений-равенств 
b_eq=[12+t, 8+4*t, 10-6*t] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#ротаеьт
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)

print("пример 2.68")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[8-5*t, 9-3*t, -3+5*t, -2-4*t] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, -1, 1, 0], 
      [-1, 2, 0, 1]] #коэффициенты из ограничений-равенств 
b_eq=[24-12*t, -18+10*t] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#ротаеьт
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)


print("пример 2.73")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[2-3*t, -4*t, 0, (2+6*t), 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, -1, 1, 0, 0], 
      [2, 1, 0, 1, 0],
      [-1, -2, 0, 0, 1]] #коэффициенты из ограничений-равенств 
b_eq=[12, 6+8*t, -10+12*t] 
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#ротаеьт
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)

print("пример 2.74")
ab=[-float("inf"),  float("inf")] #промежуток t
z=[2+t, -(3-t), 0, 0, 3*(2+4*t)] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 2, 1, 0, 0],
      [-1, 3, 0, 1, 0],
      [1, -2, 0, 0, 1]] #коэффициенты из ограничений-равенств 
b_eq=[4+6*t, 9-12*t, 8+9*t]
a_ub=[] #коэффициенты из ограничений-неравенств
b_ub=[]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)

print('таха')
print("пример 7.6.1")
#ab=[-float("inf"),  float("inf")]
ab=[0,  float("inf")] #промежуток t

z=[3-6*t, 2-2*t, 5+5*t, 0] #тут хрвняться выражения перед иксами в функции
a_eq=[[1, 2, 1, 1]] #коэффициенты из ограничений-равенств 
b_eq=[40] 
a_ub=[[3, 0, 2, 0], [1, 4, 0, 0]] #коэффициенты из ограничений-неравенств
b_ub=[60, 30]

#Работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)


print("пример 7.6.2")
ab=[0,  float("inf")]
z=[3, 2, 5] #тут хрвняться выражения перед иксами в функции
a_eq=[] #коэффициенты из ограничений-равенств 
b_eq=[] 
a_ub=[[1, 2, 1], [3, 0, 2], [1, 4, 0]] #коэффициенты из ограничений-неравенств
b_ub=[40-t, 60+2*t, 30-7*t]

#работает
Simplex_Res=SimplexAlgorithm.ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub)
SimplexAlgorithm.print_Res(Simplex_Res)
print('Результат :')
SimplexAlgorithm.print_short_Res(Simplex_Res, ab)



###############################################################################3
class Table():
    b=[]
    x=[]
    z=[]
    bazis=[]
    min=[]
    F=0

s=Table()
s.bazis=[3, 4, 5]
s.b=[-33, -23, -12]
s.x=[[-4, -3, -1],
     [-3, -2, -1],
     [-2, -1, -2],
     [1, 0, 0],
     [0, 1,0],
     [0, 0, 1]]
s.z=[-20, -20, -10, 0, 0, 0]
s.min=[-0, -0, -0]

#table, y=SimplexAlgorithm.minus_b(s, 0)

#ss=[]
#ss.append(table)
#SimplexAlgorithm.print_table(ss)