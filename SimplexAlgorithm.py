import scipy
from scipy.optimize import linprog
import numpy as np
import copy
import sympy
from sympy import *
from sympy.solvers import solve
from fractions import Fraction
from prettytable import PrettyTable

t = Symbol('t')
e=10**(-10)

#класс для хранения результатов вычислений
class Result():
    simplex_table=[] # таблица на прмежутке
    T = [] # ...T_Res


class T_Res():
    a=''
    b=''
    z1=''
    z2=''

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

def build_simplex_table(z, a_eq, b_eq, a_ub, b_ub):
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


def simplex(simplex_t,ti):
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

    # нахождение элементов
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

    #постоегие новой таблицы
    new_table=copy.deepcopy(simplex_table)
    new_table.bazis[i_minmin]=i_min
    print('bazis', new_table.bazis)
    for i in range(len(simplex_table.z)):
        new_table.x[i][i_minmin]=simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]
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
    print('x', new_table.x)
    #new z
    for i in range(len(simplex_table.z)):
        new_table.z[i]=simplex_table.z[i]-simplex_table.x[i][i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.z[i_min]
    print('z', new_table.z)
    #new F
    new_table.F=simplex_table.F-simplex_table.b[i_minmin]/simplex_table.x[i_min][i_minmin]*simplex_table.z[i_min]
    print('f', new_table.F) 

    simplex_t[len(simplex_t)-1]=simplex_table
    simplex_t.append(new_table)

    return  simplex(simplex_t,ti)

def ParameterInSimplex(ab, z, a_eq, b_eq, a_ub, b_ub):
    z, a_eq, b_eq, a_ub, b_ub=init(z, a_eq, b_eq, a_ub, b_ub)
    simplex_table=[]
    table=build_simplex_table(z, a_eq, b_eq, a_ub, b_ub)
    simplex_table.append(table)
    ti=find_ti(ab)
    print('simplex')
    s=simplex(simplex_table,ti)
    print('interval')
    t1=interval_t(s)
    r1=Result()
    r1.simplex_table=s
    r1.T=t1
    Res=[]
    Res.append(r1)
    print('loop')
    Simplex_Res=loop(ab, Res, table)
    print('Prrrrrrrrrrriiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiint')
    print_Res(Res)
    return 

def print_table(simplex_table): 
    table= PrettyTable()
    x=[]
    for i in range(len(simplex_table[0].x)):
        x.append('x'+str(i))
    row=['Базис ', 'Решение']+x+['min ']
    for i in range(len(row)):
        print(row[i], end='\t|') 
    table.field_names=row
    for i in range(len(simplex_table)):
        print()
        row=[' ']
        row.append('F='+str(simplex_table[i].F))
        main_row=[]
        for j in range(len(simplex_table[0].x)):
            row.append(str(simplex_table[i].z[j]))
        for j in range(len(simplex_table[0].b)):
            main_row.append('x'+str(simplex_table[i].bazis[j]))
            print('x'+str(simplex_table[i].bazis[j]), end='\t|')
            main_row.append(simplex_table[i].b[j])
            print(simplex_table[i].b[j], end='\t|') 
            for q in range(len(simplex_table[0].x)):
                main_row.append(simplex_table[i].x[q][j])
                print(simplex_table[i].x[q][j], end='\t|') 
            main_row.append(simplex_table[i].min[j])
            print(simplex_table[i].min[j], end='\t|') 
            print('len', len(main_row))
            table.add_row(main_row)
            main_row=[]
            print()
        row.append(' ')
        print(row)
        print(simplex_table[i].z)
        table.add_row(row)
        table.add_row([' ']*len(row))
        for j in range(len(row)):
            print(row[j], end='\t|') 
        print()
        print(table)
        
        

def print_Res(Res):
    for i in range(len(Res)):
        print('i=', i)
        print_table(Res[i].simplex_table)
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

def interval_t(simplex_table):
    print(simplex_table[len(simplex_table)-1])
    z=copy.deepcopy(simplex_table[len(simplex_table)-1].z)
    print('z', z)
    t_result=T_Res()
    for i in range(len(z)):
        z[i]=z[i]+0*t>=0
    print(z)
    r=str(solve(z, t))
    print(str(r))
    r=r.replace('&', '')
    if r.find('(t < oo)')!=-1:
        r=r.replace('(t < oo)', '')
        t_result.z2='<'
        t_result.b=float("inf")
    if r.find('(-oo < t)')!=-1:
        r=r.replace('(-oo < t)', '')
        t_result.z1='<'
        t_result.a=-float("inf")
    k=0
    while r.find('(') != -1 and k!=2:
        k+=1
        a=r.find('(')
        b=r.find(')')
        s=''
        for i in range(a+1, b):
            s+=r[i]
        r=r[:a] + r[b+1:]
        if s[len(s)-1]=='t':
            if s.find('<=')!=-1:
                t_result.z1='<='
            else:
                t_result.z1='<'
            s=s.replace('t', '')
            s=s.replace('<=', '')
            s=s.replace('<', '')
            s=s.replace(' ', '')
            t_result.a=float(Fraction(s))
        if s[0]=='t':
            if s.find('<=')!=-1:
                t_result.z2='<='
            else:
                t_result.z2='<'
            s=s.replace('t', '')
            s=s.replace('<=', '')
            s=s.replace('<', '')
            s=s.replace(' ', '')
            t_result.b=float(Fraction(s))
    return t_result


def loop(ab, Res, table):
    for i in range(len(Res)):
        print(i, Res[i].T.a, Res[i].T.b, Res[i].T.z1, Res[i].T.z2)
    if (ab[0]==Res[0].T.a and Res[0].T.z1=='<=' and ab[1]==Res[len(Res)-1].T.b and Res[len(Res)-1].T.z2=='<=') or (Res[0].T.a==-float("inf") and Res[len(Res)-1].T.b==float("inf")) or (ab[0]>Res[0].T.a and ab[1]<Res[len(Res)-1].T.b) or (ab[0]>Res[0].T.a  and Res[len(Res)-1].T.b==float("inf")) or (Res[0].T.a==-float("inf") and ab[1]<Res[len(Res)-1].T.b):
        return Res
    simplex_table=[]
    simplex_table.append(table)
    if (ab[0]==Res[0].T.a and Res[0].T.z1=='<=') or (Res[0].T.a==-float("inf")) or (ab[0]>Res[0].T.a):
        if Res[len(Res)-1].T.z2=='<=':
            t1=Res[len(Res)-1].T.b+e
        else:
            t1=Res[len(Res)-1].T.b
        simplex2=simplex(simplex_table,t1)
        print('loop2')
        t2=interval_t(simplex2)
        r2=Result()
        r2.simplex_table=simplex2
        r2.T=t2
        Res.append(r2)
        return loop(ab, Res, table)
    if (ab[1]==Res[len(Res)-1].T.b and Res[len(Res)-1].T.z2=='<=') or (Res[len(Res)-1].T.b==float("inf")) or (ab[1]<Res[len(Res)-1].T.b):
        if Res[0].T.z1=='<=':
             t0=Res[0].T.a-e
        else:
            t0=Res[0].T.a
        simplex1=simplex(simplex_table,t0)
        print('loop3')
        t1=interval_t(simplex1)
        r1=Result()
        r1.simplex_table=simplex1
        r1.T=t1
        Res.insert(0, r1)
        return loop(ab, Res, table)

    if Res[0].T.z1=='<=':
        t0=Res[0].T.a-e
    else:
        t0=Res[0].T.a
    simplex1=simplex(simplex_table,t0)
    if Res[len(Res)-1].T.z2=='<=':
        t1=Res[len(Res)-1].T.b+e
    else:
        t1=Res[len(Res)-1].T.b
    simplex2=simplex(simplex_table,t1)

    print('loop1')

    t1=interval_t(simplex1)
    t2=interval_t(simplex2)
    r1=Result()
    r1.simplex_table=simplex1
    r1.T=t1
    r2=Result()
    r2.simplex_table=simplex2
    r2.T=t2
    Res.insert(0, r1)
    Res.append(r2)
    return loop(ab, Res, table)

