import scipy
from scipy.optimize import linprog
import numpy as np
import copy
import sympy
from sympy import *
from sympy.solvers import solve

t = Symbol('t')
e=10**(-10)
#класс для хранения результатов вычислений
class Result():
    F= []  #целевая функция на промежутках
    X = [] #вектора иксов X = (x_1, ... , x_n)
    T = [] # ...T_Res


class T_Res():
    a1=None
    a2=None
    b1=''
    b2=''

class Table():
    b=[]
    x=[]
    z=[]
    bazis=[]
    min=[]
    F=0

def find_ti(ab):
    if ab==[-float("inf"),  float("inf")]:
        ti=0
    if ab[0]==-float("inf"):
        ti=ab[1]
    if ab[1]==float("inf"):
        ti=ab[0]
    return ti

def build_simplex_table(ab, z, a_eq, b_eq, a_ub, b_ub):
        #построение симплекс таблицы
    b=b_eq+b_ub
    bazis=[0]*len(b)
    n=len(b_ub)+len(a_eq[0])
    x=[ [0 for j in range(len(b))] for i in range(n)]
    z_str=[0*t]*(n)
    for i in range(len(z)):
        z_str[i]=-z[i]
    for i in range(len(b_eq)):
        for j in range(len(a_eq[i])):
            x[j][i]=a_eq[i][j]
    for i in range(len(b_ub)):
        for j in range(len(a_eq[0])):
            x[j][len(b_eq)+i]=a_ub[i][j]
    for i in range(len(a_eq[0]), len(x)):
        x[i][i-len(a_eq[0])+len(b_eq)]=1+0*t
    for i in range(n):
        if x[i].count(1)==1 and x[i].count(0)==len(b)-1:
            bazis[x[i].index(1)]=i
    #коец построения симплекс таблицы  
    
    simplex_table=Table()
    simplex_table.b=b
    simplex_table.x=x
    simplex_table.z=z_str
    simplex_table.bazis=bazis
    simplex_table.min=[0*t]*len(b)
    simplex_table.F=0
    return simplex_table


def simplex(simplex_t,ti, ab):
    simplex_table=simplex_t[len(simplex_t)-1]
    table=copy.deepcopy(simplex_table)
    table.z[0]=simplex_table.z[0].subs(t, ti)
    min=table.z[0]
    i_min=0
    for i in range(1, len(simplex_table.z)):
        table.z[i]=simplex_table.z[i].subs(t, ti)
        if min>table.z[i]:
            min=table.z[i]
            i_min=i
    if min>=0:
        return simplex_t
    for i in range(len(simplex_table.b)):
        if simplex_table.b[i]==0 or simplex_table.x[i_min][i]==0:
            simplex_table.min[i]=float("inf")+0*t
        else:
            simplex_table.min[i]=simplex_table.b[i]/simplex_table.x[i_min][i]
    minmin=simplex_table.min[0]
    i_minmin=0
    for i in range(1, len(simplex_table.min)):
        if minmin>simplex_table.min[i]:
            minmin=simplex_table.min[i]
            i_minmin=i
    print('hhhhhhhhhhhhhh')
    #постоегие новой таблицы
    new_table=copy.deepcopy(simplex_table)
    new_table.bazis[i_minmin]=i_min
    print(new_table.bazis)
    for i in range(len(simplex_table.z)):
        new_table.x[i][i_minmin]=simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]
        print(new_table.x[i][i_minmin])
    #new b
    new_table.b[i_minmin]=simplex_table.b[i_minmin]/simplex_table.x[i_min][i_minmin]
    for i in range(i_minmin):
        new_table.b[i]=simplex_table.b[i]-simplex_table.b[i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.x[i_min][i]
    for i in range(i_minmin+1, len(new_table.b)):
        new_table.b[i]=simplex_table.b[i]-simplex_table.b[i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.x[i_min][i]
    #new x
    for i in range(len(simplex_table.z)):
        new_table.x[i][i_minmin]=simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]
    for j in range(i_minmin):
        for i in range(len(simplex_table.z)):
            new_table.x[i][j]=simplex_table.x[i][j]-simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.x[i_min][j]
    for j in range(i_minmin+1, len(simplex_table.b)):
        for i in range(len(simplex_table.z)):
            new_table.x[i][j]=simplex_table.x[i][j]-simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.x[i_min][j]
    print(new_table.x)
    #new z
    for i in range(len(simplex_table.z)):
        new_table.z[i]=simplex_table.z[i]-simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.z[i_min]
    print(new_table.z)
    #new F
    new_table.F=simplex_table.F-simplex_table.b[i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.z[i_min]
    print(new_table.F)

    simplex_t[len(simplex_t)-1]=simplex_table
    simplex_t.append(new_table)

    return  simplex(simplex_t,ti, ab)

def ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub):
    z, a_eq, b_eq, a_ub, b_ub=init(z, a_eq, b_eq, a_ub, b_ub)
    simplex_table=[]
    table=build_simplex_table( ab, z, a_eq, b_eq, a_ub, b_ub)
    simplex_table.append(table)
    ti=find_ti(ab)
    print(ti)
    Res=simplex(simplex_table,ti, ab )
    return Res

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