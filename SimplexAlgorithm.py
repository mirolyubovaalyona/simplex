import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve

t = Symbol('t')

#класс для хранения результатов вычислений
class Result():
    F= []  #целевая функция на промежутках
    X = [] #вектора иксов X = (x_1, ... , x_n)
    T = [] # ...


class T_Res():
    a1=None
    a2=None
    b1=''
    b2=''

def build_simplex_table(t0, ab, z, a_eq, b_eq, a_ub, b_ub):

    #построение симплекс таблицы
    b=b_eq+b_ub
    bazis=[0]*len(b)
    min=[]
    n=len(b_ub)+len(a_eq[0])
    x=[ [0 for j in range(len(b))] for i in range(n)]
    z_str=[0]*(len(z)+1)
    for i in range(len(z)):
        z_str[i]=-z[i]
    for i in range(len(b_eq)):
        for j in range(len(a_eq[i])):
            x[j][i]=a_eq[i][j]
    for i in range(len(b_ub)):
        for j in range(len(a_eq[0])):
            x[j][len(b_eq)+i]=a_ub[i][j]
    for i in range(len(a_eq[0]), len(x)):
        x[i][i-len(a_eq[0])+len(b_eq)]=1
    for i in range(n):
        if x[i].count(1)==1 and x[i].count(0)==len(b)-1:
            bazis[x[i].index(1)]=i
    #коец построения симплекс таблицы   
    print(b)
    print(x)
    print(z)
    print(bazis)


   
    return

def simplex():
    return

def ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub):
    z, a_eq, b_eq, a_ub, b_ub=init(z, a_eq, b_eq, a_ub, b_ub)
    t0=0
    build_simplex_table(t0, ab, z, a_eq, b_eq, a_ub, b_ub)
    return

def print_table():
    return

def init(z, a_eq, b_eq, a_ub, b_ub):
    for i in range(len(z)):
        z[i]+=0*t
    for i in range(len(b_eq)):
        b_eq[i]+=0*t
    for i in range(len(b_ub)):
        b_ub[i]+=0*t
    for i in range(len(a_eq)):
        for j in range(len(a_eq[0])):
            a_eq[i][j]+=0*t
    for i in range(len(a_ub)):
        for j in range(len(a_ub[0])):
            a_ub[i][j]+=0*t
    return z, a_eq, b_eq, a_ub, b_ub