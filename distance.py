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
import copy as cp
from astropy.tests import zmf as z
#from astropy import constants as con

#=============================assistant code===================================
#--------------------------------plot------------------------------------------
    
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
    cbar = plt.colorbar()
    cbar.set_label('S (Jy/beam)')
    
    plt.xlabel('l (deg)')
    plt.ylabel('b (deg)')  
    
    x,y,r =z.coo_circle(head,onoff)
    
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(onoff[2]*np.cos(an)+onoff[0], onoff[2]*np.sin(an)+onoff[1],'r')
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
    
#------------------about spectra,rotation curve,Jy/beam->K---------------------    
def velocity(data,head):
    '''
    return spec velocity axis
    '''
    w   = wcs.WCS(head)
    pix = np.linspace(1,data.shape[0],data.shape[0])
    if head['NAXIS'] == 3:

        x,y,v   = w.wcs_pix2world(0,0,pix,0)
    if head['NAXIS'] == 4:

        x,y,v,s   = w.wcs_pix2world(0,0,pix,0,0)
    return v
    
def mod(file):
    '''
    get Milky Way rotation curve
    '''
    V_R = np.loadtxt(file)
    fig, ax = plt.subplots() 
    plt.title('rotation curve')
    plt.plot(V_R[:,0],V_R[:,1]+30)
    return V_R 

#=============================continuum========================================
def continuum(file,analyze,region,on,off,contrast):
    cont = fits.open(file)
    cont_head = cont[0].header
    cont_data = cont[0].data[0,:,:]*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
    cont.close()
    cont_head['CUNIT3'] = 'Hz'
    
    z.plot_fits(cont_data,cont_head,1,'origin')
    
    if analyze == 'box':
        if on != []:
            cont_on  = box(cont_data,cont_head,1,'cont_on',region,on)
        if off != []:
            cont_off = box(cont_data,cont_head,1,'cont_off',region,off)
    elif analyze == 'circle':
            cont_on,cont_reg  = circle(cont_data,cont_head,1,'1420MHz continuum (K)',region,on)
            cont_off,cont_reg = circle(cont_data,cont_head,1,'cont_off',region,off)
    return cont_on,cont_off

#=============================spectra==========================================
def spectra(file,analyze,region,on,off,contrast,spec_v):
    spec = fits.open(file)
    spec_head = spec[0].header
    spec_data = spec[0].data[:,:,:]*z.conversion(1.4,spec_head['BMAJ'],spec_head['BMIN'])
    spec.close()
    spec_head['CUNIT3'] = 'm/s'
        
    v = velocity(spec_data,spec_head)
        
    if analyze == 'box':
        if on != []:
            spec_on  = box(spec_data,spec_head,1,'cont_on',region,on,spec_v)
        if off != []:    
            spec_off = box(spec_data,spec_head,1,'cont_off',region,off,spec_v)
    elif analyze == 'circle':
        if on != []:
            spec_on,spec_reg  = circle(spec_data,spec_head,1,'1720MHz spectrum map (K) at '+str(int(v[spec_v]))+' m/s',region,on,spec_v)
        if off != []:    
            spec_off,spec_reg = circle(spec_data,spec_head,1,'cont_off',region,off,spec_v)
    return spec_on,spec_off,v
#=============================absorption=======================================
def absorption_spec(spec_on,spec_off,v,cont_on,cont_off,on,off,analyze,method): 
    if analyze   == 'box':
        if method == 'tww':
            T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
            T_off    = (np.sum(np.sum(spec_off,axis=1),axis=1)                    \
                        -(np.sum(np.sum(spec_on,axis=1),axis=1)))                 \
                        /(spec_off.shape[1]*spec_off.shape[2]-spec_on.shape[1]*spec_on.shape[2])
            T_con    = np.mean(cont_on)
            T_coff   = (np.sum(cont_off)-np.sum(cont_on))                         \
                        /(cont_off.shape[0]*cont_off.shape[1]-cont_on.shape[0]*cont_on.shape[1])
            e_tau    = (T_on-T_off)/(T_con-T_coff)       
        else:
            T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
            T_off    = np.mean(np.mean(spec_off,axis=1),axis=1)
            T_con    = np.mean(cont_on)
            T_coff   = np.mean(cont_off)
            e_tau    = (T_on-T_off)/(T_con-T_coff)       
#            e_tau    = (np.mean(np.mean(spec_on,axis=1),axis=1)                 \
#                        -np.mean(np.mean(spec_off,axis=1),axis=1))                \
#                        /(np.mean(cont_on)-np.mean(cont_off))
    elif analyze == 'circle':
        T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
        T_off    = (np.sum(np.sum(spec_off,axis=1),axis=1)                    \
                    -(np.sum(np.sum(spec_on,axis=1),axis=1)))                 \
                    /(spec_off.shape[1]*spec_off.shape[2]-spec_on.shape[1]*spec_on.shape[2])
        T_con    = np.mean(cont_on)
        T_coff   = (np.sum(cont_off)-np.sum(cont_on))                         \
                    /(cont_off.shape[0]*cont_off.shape[1]-cont_on.shape[0]*cont_on.shape[1])
        e_tau    = 1+(T_on-T_off)/(T_con-T_coff)            
 
    v = v/1000
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(v[70:], T_on[70:],label='on')
    ax1.plot(v[70:], T_off[70:],label='off')
    ax1.set_ylabel('T(K)')
    props = font_manager.FontProperties(size=10)
    ax1.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
    ax1.set_title('absorption spectrum')
#    ax1.set_ylim(y2lim[0],y2lim[1])
    ax2.plot(v[70:], e_tau[70:])
    ax2.set_ylabel(r'$e^{-\tau}$',fontsize=20)
    ax2.set_xlabel('velocity(km/s)')
#    ax2.set_xlim(xlim[0],xlim[1])
#    ax2.set_ylim(ylim[0],ylim[1])
    fig.subplots_adjust(hspace=0.05)
    plt.legend()
    plt.show()   

#=============================distance=========================================    
def dist(model,file,l,b,d,V = 220,v_sun = 220,r_sun = 8.5):
    if model == 'constant':
        v_mod = V
    if model == 'model':
        v_mod = V
        mod(file)            
    b = np.deg2rad(b)
    l = np.deg2rad(l)
    r = (r_sun**2+(d*np.cos(b))**2-2*r_sun*d*np.cos(b)*np.cos(l))**0.5
    v = v_mod*r_sun*np.sin(l)*np.cos(b)/r-v_sun*np.sin(l)*np.cos(b) 
    fig, ax = plt.subplots() 
    plt.title('distance-velocity')
    plt.plot(d,v)
    plt.xlabel('distance (kpc)')
    plt.ylabel('velocity (km/s)')
    return v,d   
#===============================main===========================================
if __name__=='__main__':
    file1   = '../data/CONT_16deg_1400mhz_25arc_thor_vgps.fits'
    file2   = '../data/THOR_HI_with_continuum_L16.25_image.fits'
    file3   = '../data/rotation_model.txt'
    region  = [15.8,15.95,0.1,0.25]      #region l1,l2,b1,b2
    on      = [15.9,15.92,0.17,0.187] 
    off     = [15.9,15.94,0.17,0.187]
    contrast = 1
    analyze  = 'box'               # box,circle
    spec_v   = 85
    model   = 'constant'            #constant, model
    V       = 220                   #km/s
    d       = np.linspace(1,40,100)
    l       = 15.9
    b       = 0.2
    method  = 'classic' #获得吸收谱的方法，tww或者classic
    cont_on,cont_off    = continuum(file1,analyze,region,on,off,contrast)
    spec_on,spec_off,v  = spectra(file2,analyze,region,on,off,contrast,spec_v)
    absorption_spec(spec_on,spec_off,v,cont_on,cont_off,on,off,analyze,method)
    v,d = dist(model,file3,l,b,d,V = 220,v_sun = 220,r_sun = 8.5)
    


