# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 19:43:19 2017

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as un
from astropy import constants as con
from astropy.modeling import models, fitting
e_t = e_tau.clip(min=1e-8)
x   = v
y   = -np.log(e_t)
T_s = 100
S   = 0
for i in y:
    S = S + i*1.5
N_HI = 1.9e18*S*T_s

print(N_HI)

x = v_co
y = T_on_co
S   = 0
for i in y[0:600]:
    S = S + i*0.21
N_H2 = 1.8e20*S*30

print(N_H2)

d=14*un.kpc
theta1=0.03/180*np.pi
theta2=0.02/180*np.pi
M_H2 = N_H2*theta1*theta2*d.to('cm').value**2*2*con.m_p.value/con.M_sun.value
                              
print(M_H2)

S=6*un.pc
n=N_H2/S.to('cm').value
           
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

#T=np.sum(spec_on_co[120:125],axis=0)
#plt.imshow(T,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
#cbar = plt.colorbar()
#cbar.set_label('T(K)')
#plt.contour(cont_reg,levels=levels,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
#plt.xlabel('l(deg)')
#plt.ylabel('b(deg)')
#onoff = [15.91,0.09,0.01]
#an = np.linspace(0, 2*np.pi, 100)
#plt.plot(onoff[2]*np.cos(an)+onoff[0], onoff[2]*np.sin(an)+onoff[1],'r')
#onoff = [15.91,0.09,0.2]
#plt.plot(onoff[2]*np.cos(an)+onoff[0], onoff[2]*np.sin(an)+onoff[1],'pink')
#plt.xlim(l2,l1)
#plt.ylim(b1,b2)
