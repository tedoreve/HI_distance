# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 14:20:42 2016

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.tests import zmf as z
from scipy.optimize import leastsq

#=============================assistant code===================================
def box(data,head,contrast,name,region,onoff,*args):
    '''
    plot box pixel coordinates
    '''
    l1,l2,b1,b2 = region  
    x1,y1,x2,y2 = z.coo_box(head,region)
    if data.ndim == 2:
        result0 = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result0 = data[args[0],y1:y2,x1:x2] 
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast
    
    plt.subplots() 
    plt.title(name)
    plt.imshow(result0,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    plt.colorbar()
    plt.xlabel(r'$l (deg)$')
    plt.ylabel(r'$b (deg)$')    
    
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
    plt.plot(x, y, 'go-')
    plt.xlim(l2,l1)
    plt.ylim(b1,b2)
    plt.grid()
    
    x1,y1,x2,y2 = z.coo_box(head,onoff)  
    if data.ndim == 2:
        result = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result = data[:,y1:y2,x1:x2]
    return result


def circle(data,head,contrast,name,region,onoff,*args):
    '''
    plot circle pixel coordinates
    '''
    l1,l2,b1,b2 = region  
    x1,y1,x2,y2 = z.coo_box(head,region)
    if data.ndim == 2:
        result0 = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result0 = data[args[0],y1:y2,x1:x2] 
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast
    
    plt.subplots()
    plt.title(name)
    plt.imshow(result0,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    plt.colorbar()
    plt.xlabel(r'$l (deg)$')
    plt.ylabel(r'$b (deg)$')  
    
    x,y,r =z.coo_circle(head,onoff)
    
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(onoff[2]*np.cos(an)+onoff[0], onoff[2]*np.sin(an)+onoff[1])
    plt.xlim(l2,l1)
    plt.ylim(b1,b2)
    plt.grid()
    
    if data.ndim == 2:
        result = cp.copy(data[int(y-r):int(y+r),int(x-r):int(x+r)])
        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                if (i-result.shape[0]/2)**2+(j-result.shape[1]/2)**2 > r**2:
                    result[i,j]=0
    if data.ndim == 3:
        result = cp.copy(data[:,int(y-r):int(y+r),int(x-r):int(x+r)])
        for i in range(result.shape[1]):
            for j in range(result.shape[2]):
                if (i-result.shape[1]/2)**2+(j-result.shape[2]/2)**2 > r**2:
                    result[:,i,j]=0

    return result,result0

def read(file):
    head = []
    data = []
    for i in range(len(file)):
        cont = fits.open(file[i])
        head.append(cont[0].header)
        data.append(cont[0].data[0,:,:])#*z.conversion(cont_head['BMAJ'],cont_head['BMIN'])    
        cont.close()
 
    return head,data
    
def fun(v,p):
    a,b = p
    return a*v + b
    
def residuals(p, flux, v):
    """
    实验数据x, y和拟合函数之间的差，p为拟合需要找到的系数
    """
    return (flux - fun(v,p))

#=============================continuum========================================
def continuum(data,head,analyze,region,on,contrast):
#    head['CUNIT3'] = 'Hz'
    
    z.plot_fits(data[2],head[2],1,'origin')
    flux_on  = []
    flux_reg = []
    for i in range(len(head)):
        if analyze == 'box':
            if on != []:
                cont_on  = box(data,head,1,'cont_on',region,on)
        elif analyze == 'circle':
                cont_on,cont_reg  = circle(data[i],head[i],1,str(i)+' continuum (K)',region,on)
        flux_on.append(cont_on)
        flux_reg.append(cont_reg)
    return flux_on,flux_reg  
    
#===============================main===========================================
if __name__=='__main__':
    file    = ['../data/THOR_cont_1060MHz_L49.25deg_25arcsec_image.fits',     \
               '../data/THOR_cont_1310MHz_L49.25deg_25arcsec_image.fits',     \
               '../data/THOR_cont_1440MHz_L49.25deg_25arcsec_image.fits',     \
               '../data/THOR_cont_1690MHz_L49.25deg_25arcsec_image.fits',     \
               '../data/THOR_cont_1820MHz_L49.25deg_25arcsec_image.fits',     \
               '../data/THOR_cont_1950MHz_L49.25deg_25arcsec_image.fits']
    region  = [49.4,49.6,-0.5,-0.3]      #region l1,l2,b1,b2
    on      = [49.5,-0.4,0.005]
    contrast = 1
    analyze  = 'circle'               # box,circle
    spec_v   = 85
    head, data = read(file)
    flux_on,flux_reg = continuum(data,head,analyze,region,on,contrast)
    
    v = [1060,1310,1440,1690,1820,1950]
    
    flux = []
    for i in range(len(head)):
        flux.append(np.mean(flux_reg[i]))
    
    v = np.log(v)
    flux = np.log(flux)
#==============================analysis========================================
    p0=[0.1,0.1]
    plsq, pcov, infodict, errmsg, success= leastsq(residuals,p0,args=(flux,v),full_output=1, epsfcn=0.0001)
    plt.scatter(v,flux)
    plt.plot(v,fun(v,plsq))


