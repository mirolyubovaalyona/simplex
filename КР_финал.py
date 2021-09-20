import scipy
from scipy.optimize import linprog
import numpy as np
import sympy
from sympy import *
from sympy.solvers import solve


import SimplexAlgorithm



t = Symbol('t')
class T_Res():
    a=''
    b=''
    z1=''
    z2=''

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
    z[i]=z[i]+t<=0
#print(z)
r=solve([ t<=0], t)
#print(str(r), str(r)[1])

import random

z = 4
x = random.sample(range(50), z)
y = random.sample(range(50), z)

#for i in x:
  #  for u in y:        
      #  print('{:6d}'.format(i*u), end='| ')   
    #print('')

r=str(r)
t_result=T_Res()
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
        t_result.a=float(s)
    if s[0]=='t':
        if s.find('<=')!=-1:
            t_result.z2='<='
        else:
            t_result.z2='<'
        s=s.replace('t', '')
        s=s.replace('<=', '')
        s=s.replace('<', '')
        s=s.replace(' ', '')
        t_result.b=float(s)

#print(t_result.a, t_result.b, t_result.z1, t_result.z2)