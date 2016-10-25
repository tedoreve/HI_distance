# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 08:58:41 2016

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
#matplotlib.use('Agg')
import matplotlib.path as mpath
from astropy.io import fits
from astropy import wcs
from astropy import units as un
from astropy import constants as con

file    = ['../data/VGPS_cont_MOS049.fits','../data/MOS_049.Tb.fits','../data/rotation_model.txt']
region  = [48.5,49.5,-1,0]      #region l1,l2,b1,b2
on      = [49,49.3,-0.76,-0.65]
off     = [49,49.3,-0.9,-0.65]
analyze  = ''               # cont,spec,both
spec_v  = 100
xlim    = [-100000,100000]
ylim    = [-0.5,2]
spectrum= True
dist    = True
model   = 'constant'            #constant, model
V       = 220                   #km/s
d       = np.linspace(1,40,100)
l       = 49.2
b       = 0
y2lim   = [-1,120]      
#=============================assistant code===================================
def plot_origin(data,head,contrast,name):
    w = wcs.WCS(head)
    if head['NAXIS'] == 2:
        l1,b1 = w.wcs_pix2world(0,0,0)
        l2,b2 = w.wcs_pix2world(data.shape[1],data.shape[0],0)
    if head['NAXIS'] == 3:
        l1,b1,v = w.wcs_pix2world(0,0,0,0)
        l2,b2,v = w.wcs_pix2world(data.shape[1],data.shape[0],0,0)    
    if head['NAXIS'] == 4:
        l1,b1,v,s = w.wcs_pix2world(0,0,0,0,0)
        l2,b2,v,s = w.wcs_pix2world(data.shape[1],data.shape[0],0,0,0)

    result = data
    result = (result-np.mean(result))+np.mean(result)*contrast
    
    fig, ax = plt.subplots() 
    plt.title(name)
    plt.imshow(np.log(result),origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    ax.grid()
    

    
    
def plot(data,head,contrast,name,region,*args):
    l1,l2,b1,b2 = region
    x1,y1,x2,y2 = coo_box(head,region)
    
    if data.ndim == 2:
        result = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result = data[args[0],y1:y2,x1:x2] 
    result = (result-np.mean(result))+np.mean(result)*contrast
    
    fig, ax = plt.subplots() 
    plt.title(name)
    plt.imshow(np.log(result),origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    ax.grid()
    
def coo_box(head,region):
    l1,l2,b1,b2 = region
    w = wcs.WCS(head)
    if head['NAXIS'] == 2:
        x1,y1 = w.wcs_world2pix(l2,b1,0)
        x2,y2 = w.wcs_world2pix(l1,b2,0)   
    if head['NAXIS'] == 3:
        x1,y1,v = w.wcs_world2pix(l2,b1,0,0)
        x2,y2,v = w.wcs_world2pix(l1,b2,0,0)
    if head['NAXIS'] == 4:
        x1,y1,v,s = w.wcs_world2pix(l2,b1,0,0,0)
        x2,y2,v,s = w.wcs_world2pix(l1,b2,0,0,0)
    return x1,y1,x2,y2

def box(data,head,contrast,name,region,onoff,*args):
    l1,l2,b1,b2 = region  
    x1,y1,x2,y2 = coo_box(head,region)
    if data.ndim == 2:
        result = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result = data[args[0],y1:y2,x1:x2] 
    result = (result-np.mean(result))+np.mean(result)*contrast
    
    fig, ax = plt.subplots() 
    plt.title(name)
    plt.imshow(np.log(result),origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    
    Path = mpath.Path
    path_data = [
        (Path.MOVETO, (onoff[0], onoff[2])),
        (Path.CURVE4, (onoff[0], onoff[3])),
        (Path.CURVE4, (onoff[1], onoff[3])),
        (Path.CURVE4, (onoff[1], onoff[2])),
        (Path.CLOSEPOLY, (onoff[0], onoff[2])),   
        ]
    codes, verts = zip(*path_data)
    path = Path(verts, codes)
    x, y = zip(*path.vertices)
    ax.plot(x, y, 'go-')
    ax.grid()
    
    x1,y1,x2,y2 = coo_box(head,onoff)  
    return data[y1:y2,x1:x2]

def circle(head,l,b,r):
    pass

def spec_box(data,head,onoff):
    x1,y1,x2,y2 = coo_box(head,onoff)
    return data[:,y1:y2,x1:x2]
    
def velocity(data,head):
    w   = wcs.WCS(head)
    pix = np.linspace(1,data.shape[0],data.shape[0])
    x,y,v,s   = w.wcs_pix2world(0,0,pix,0,0)
    return v
    
def mod():
    V_R = np.loadtxt(file[2])
    fig, ax = plt.subplots() 
    plt.title('rotation curve')
    plt.plot(V_R[:,0],V_R[:,1]+30)
    return V_R

def v_d(model,l,b,d,V = 220,v_sun = 220,r_sun = 8.5):
    if model == 'constant':
        v_mod = V
    if model == 'model':
        v_mod = V
        mod()            
    b = np.deg2rad(b)
    l = np.deg2rad(l)
    r = (r_sun**2+(d*np.cos(b))**2-2*r_sun*d*np.cos(b)*np.cos(l))**0.5
    v = v_mod*r_sun*np.sin(l)*np.cos(b)/r-v_sun*np.sin(l)*np.cos(b) 
    fig, ax = plt.subplots() 
    plt.title('distance-velocity')
    plt.plot(d,v)
    return v,r   
#=============================continuum========================================
cont = fits.open(file[0])
cont_head = cont[0].header
cont_data = cont[0].data[0,0,:,:]
cont_head['CUNIT3'] = 'm/s'
plot_origin(cont_data,cont_head,0,'origin')

if analyze == 'cont' or analyze == 'both':
    plot(cont_data,cont_head,0,'continuum',region)
    if on != []:
        cont_on  = box(cont_data,cont_head,0,'cont_on',region,on)
    if off != []:
        cont_off = box(cont_data,cont_head,0,'cont_off',region,off)

#=============================spectrum=========================================
spec = fits.open(file[1])
spec_head = spec[0].header
spec_data = spec[0].data[0,:,:,:]
spec_head['CUNIT3'] = 'm/s'
    
if analyze == 'spec' or analyze == 'both':
    plot(spec_data,spec_head,0,'spectrum',region,spec_v)
    if on != []:
        spec_on0  = box(spec_data,spec_head,0,'cont_on',region,on,spec_v)
    if off != []:    
        spec_off0 = box(spec_data,spec_head,0,'cont_off',region,off,spec_v)

#=============================optical depth====================================
if spectrum:
    spec_on  = spec_box(spec_data,spec_head,on)
    spec_off = spec_box(spec_data,spec_head,off)
    e_tau    = 1+(np.mean(np.mean(spec_on,axis=1),axis=1)                     \
                -np.mean(np.mean(spec_off,axis=1),axis=1))                    \
                /(np.mean(cont_on)-np.mean(cont_off))               
    
    v        = velocity(spec_data,spec_head) 
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    l1 = ax1.plot(v, np.mean(np.mean(spec_on,axis=1),axis=1),label='on')
    l2 = ax1.plot(v, np.mean(np.mean(spec_off,axis=1),axis=1),label='off')
    props = font_manager.FontProperties(size=10)
    ax1.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
    ax1.plot()
    ax1.set_title('absorption spectrum')
    ax1.set_ylim(y2lim[0],y2lim[1])
    l3 = ax2.plot(v, e_tau)
    ax2.set_xlim(xlim[0],xlim[1])
    ax2.set_ylim(ylim[0],ylim[1])
    fig.subplots_adjust(hspace=0)
    plt.legend()
    plt.show()
    

#=============================distance=========================================    
if dist:
    v_model = v_d(model,l,b,d)

