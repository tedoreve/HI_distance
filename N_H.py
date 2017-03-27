# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 19:43:19 2017

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as un
from astropy import constants as con

e_t = e_tau.clip(min=1e-2)
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

r=4*un.pc
n=1.7e22/r.to('cm').value
           
print(n)

