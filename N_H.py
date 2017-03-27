# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 19:43:19 2017

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as un
from astropy import constants as con

e_t = e_tau.clip(min=1e-32)
x   = v
y   = -np.log(e_t)
T_s = 50
S   = 0
for i in y:
    S = S + i*1.5
N_HI = 1.9e18*S*T_s

print(N_HI)

x = v_co
y = T_on_co
S   = 0
for i in y:
    S = S + i*0.21
N_H2 = 1.8e20*S*30

print(N_H2)

d=14*un.kpc
theta1=0.03/180*np.pi
theta2=0.02/180*np.pi
M_H2 = N_H2*theta1*theta2*d.to('cm').value**2*2*con.m_p.value/con.M_sun.value
                              
print(M_H2)

r=12*un.pc
n=5.1e22/r.to('cm').value
           
print(n)

t=4.5e15/n
def f(E):
    return (59e-15/7**-2.88)*E**-2.88

S   = 0
delta = (7-0.5)/1000
for i in range(1001):
    E = 0.5 + i*delta
    S = S + f(E)*delta*E
S=S*un.TeV
print(S.to('erg').value*4*np.pi*d.to('cm').value**2)
print(S.to('erg').value*4*np.pi*d.to('cm').value**2*t)
