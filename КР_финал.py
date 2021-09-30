import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve
import copy


import SimplexAlgorithm


t = Symbol('t')


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

#симплекс метод без параметра
ab=[0,  10] #промежуток t
z=[2, 13] #тут хрвняться выражения перед иксами в функции
a_eq=[] #коэффициенты из ограничений-равенств 
b_eq=[] 
a_ub=[[4, 1], [2, 2], [6, 3]] #коэффициенты из ограничений-неравенств
b_ub=[16, 22, 36]


z, a_eq, b_eq, a_ub, b_ub=SimplexAlgorithm.init(z, a_eq, b_eq, a_ub, b_ub)
simplex_table=[]
table=SimplexAlgorithm.build_simplex_table(z, a_eq, b_eq, a_ub, b_ub)
simplex_table.append(table)
ti=0
s, y=copy.deepcopy(SimplexAlgorithm.simplex(simplex_table,ti))
SimplexAlgorithm. print_table(s)

 