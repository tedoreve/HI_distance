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
from astropy import units as u
import copy as cp
#from astropy import constants as con
    
#=============================assistant code===================================

#------------------------------coordinates-------------------------------------
def coo_box(head,region):
    '''
    get box pixel coordinates 
    '''
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
        
def coo_circle(head,region):
    '''
    get circle pixel coordinates
    '''
    l,b,r = region
    r = r/np.abs(head['CDELT1'])
    w = wcs.WCS(head)
    if head['NAXIS'] == 2:
        x,y = w.wcs_world2pix(l,b,0)
    if head['NAXIS'] == 3:
        x,y,v = w.wcs_world2pix(l,b,0,0)
    if head['NAXIS'] == 4:
        x,y,v,s = w.wcs_world2pix(l,b,0,0,0)
    return x,y,r

#--------------------------------plot------------------------------------------
def plot_origin(data,head,contrast,name):
    '''
    plot original continuum figure
    '''
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

    result0 = data
    result0 = np.nan_to_num(result0)
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast

    plt.subplots() 
    plt.title(name)
    plt.imshow(np.log(result0),origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    plt.grid()
    

def box(data,head,contrast,name,region,onoff,*args):
    '''
    plot box pixel coordinates
    '''
    l1,l2,b1,b2 = region  
    x1,y1,x2,y2 = coo_box(head,region)
    if data.ndim == 2:
        result0 = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result0 = data[args[0],y1:y2,x1:x2] 
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast
    
    plt.subplots() 
    plt.title(name)
    plt.imshow(np.log(result0),origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    
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
    
    x1,y1,x2,y2 = coo_box(head,onoff)  
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
    x1,y1,x2,y2 = coo_box(head,region)
    if data.ndim == 2:
        result0 = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result0 = data[args[0],y1:y2,x1:x2] 
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast
    
    plt.subplots()
    plt.title(name)
    plt.imshow(result0,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    
    x,y,r = coo_circle(head,onoff)
    
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
    return result
    
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
    
def conversion(bmaj,bmin):
    '''
    Jy/beam -> K
    '''
    bmaj = bmaj*u.deg
    bmin = bmin*u.deg
    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
    freq = 1.4*u.GHz
    equiv = u.brightness_temperature(beam_area, freq)
    return u.Jy.to(u.K, equivalencies=equiv)  

#=============================continuum========================================
def continuum(file,analyze,region,on,off,contrast):
    cont = fits.open(file)
    cont_head = cont[0].header
    cont_data = cont[0].data[0,:,:]*conversion(cont_head['BMAJ'],cont_head['BMIN'])
    cont.close()
    cont_head['CUNIT3'] = 'Hz'
    
    plot_origin(cont_data,cont_head,0,'origin')
    
    
    if analyze == 'box':
        if on != []:
            cont_on  = box(cont_data,cont_head,0,'cont_on',region,on)
        if off != []:
            cont_off = box(cont_data,cont_head,0,'cont_off',region,off)
    elif analyze == 'circle':
            cont_on  = circle(cont_data,cont_head,0,'cont_on',region,on)
            cont_off = circle(cont_data,cont_head,0,'cont_off',region,off)
    return cont_on,cont_off

#=============================spectra==========================================
def spectra(file,analyze,region,on,off,contrast,spec_v):
    spec = fits.open(file)
    spec_head = spec[0].header
    spec_data = spec[0].data[:,:,:]*conversion(spec_head['BMAJ'],spec_head['BMIN'])
    spec.close()
    spec_head['CUNIT3'] = 'm/s'
        
    v = velocity(spec_data,spec_head)
        
    if analyze == 'box':
        if on != []:
            spec_on  = box(spec_data,spec_head,0,'cont_on',region,on,spec_v)
        if off != []:    
            spec_off = box(spec_data,spec_head,0,'cont_off',region,off,spec_v)
    elif analyze == 'circle':
        if on != []:
            spec_on  = circle(spec_data,spec_head,0,'cont_on',region,on,spec_v)
        if off != []:    
            spec_off = circle(spec_data,spec_head,0,'cont_off',region,off,spec_v)
    return spec_on,spec_off,v

#=============================absorption=======================================
def absorption_spec(spec_on,spec_off,v,cont_on,cont_off,on,off,analyze): 
    if analyze   == 'box':
        e_tau    = 1+(np.mean(np.mean(spec_on,axis=1),axis=1)                 \
                    -np.mean(np.mean(spec_off,axis=1),axis=1))                \
                    /(np.mean(cont_on)-np.mean(cont_off))     
    elif analyze == 'circle':
        T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
        T_off    = (np.sum(np.sum(spec_off,axis=1),axis=1)                    \
                    -(np.sum(np.sum(spec_on,axis=1),axis=1)))                 \
                    /(spec_off.shape[1]*spec_off.shape[2]-spec_on.shape[1]*spec_on.shape[2])
        T_con    = np.mean(cont_on)
        T_coff   = (np.sum(cont_off)-np.sum(cont_on))                         \
                    /(cont_off.shape[0]*cont_off.shape[1]-cont_on.shape[0]*cont_on.shape[1])
        e_tau    = 1+(T_on-T_off)#/(T_con-T_coff)            
 
    v = v/1000
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(v[70:], T_on[70:],label='on')
    ax1.plot(v[70:], T_off[70:],label='off')
    ax1.set_ylabel('T(K)')
    props = font_manager.FontProperties(size=10)
    ax1.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
    ax1.set_title('absorption spectra')
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
    return v,r   
#=============================OH 1720 1665 1612================================
def OH(file):
    v     = []
    T_on  = []
    T_off = []
    for i in range(len(file)):
        spec = fits.open(file[i])
        spec_head = spec[0].header
        spec_data = spec[0].data[:,:,:]*conversion(spec_head['BMAJ'],spec_head['BMIN'])
        spec.close()
        spec_head['CUNIT3'] = 'm/s'
            
        v.append(velocity(spec_data,spec_head))
        T_on.append(np.mean(np.mean(spec_on,axis=1),axis=1))
        T_off.append(np.mean(np.mean(spec_off,axis=1),axis=1))
        
    v = v/1000
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(v[0], T_on[0],label='on')
    ax1.plot(v[0], T_off[0],label='off')
    ax1.set_ylabel('T(K)')
    props = font_manager.FontProperties(size=10)
    ax1.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
    ax1.set_title('absorption spectra')
    ax2.plot(v[0], T_on[0],label='on')
    ax2.plot(v[0], T_off[0],label='off')
    ax2.set_ylabel('T(K)')
    ax2.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
    ax3.plot(v[0], T_on[0],label='on')
    ax3.plot(v[0], T_off[0],label='off')
    ax3.set_ylabel('T(K)')
    ax3.set_xlabel('velocity(km/s)')
    ax3.legend(loc='upper left', shadow=True, fancybox=True, prop=props)    
#    ax1.set_ylim(y2lim[0],y2lim[1])
#    ax2.set_xlim(xlim[0],xlim[1])
#    ax2.set_ylim(ylim[0],ylim[1])
    fig.subplots_adjust(hspace=0.05)
    plt.legend()
    plt.show()   

#===============================main===========================================
if __name__=='__main__':
    file1   = '../data/CONT_49deg_1400mhz_25arc_thor_vgps.fits'
    file2   = '../data/OH_1665mhz_L49.25_deg.smooth20sec.fits'
    file3   = '../data/rotation_model.txt'
    file4   = ['../data/OH_1720mhz_L49.25_deg.smooth20sec.fits','../data/OH_1665mhz_L49.25_deg.smooth20sec.fits','../data/OH_1612mhz_L49.25_deg.smooth20sec.fits']
    region  = [49.1,49.3,-0.4,-0.2]      #region l1,l2,b1,b2
    on      = [49.210,-0.340,0.015]
    off     = [49.210,-0.340,0.025]
#    on      = [49.176,-0.326,0.006]
#    off     = [49.160,-0.310,0.006] OH1720
    contrast = 1
    analyze  = 'circle'               # box,circle
    spec_v   = 84
    #xlim    = [-100000,100000]
    #ylim    = [-0.5,2]
    model   = 'constant'            #constant, model
    V       = 220                   #km/s
    d       = np.linspace(1,40,100)
    l       = 49.2
    b       = 0
    #y2lim   = [-1,120]  
    cont_on,cont_off    = continuum(file1,analyze,region,on,off,contrast)
    spec_on,spec_off,v  = spectra(file2,analyze,region,on,off,contrast,spec_v)
    absorption_spec(spec_on,spec_off,v,cont_on,cont_off,on,off,analyze)
    v,r = dist(model,file3,l,b,d,V = 220,v_sun = 220,r_sun = 8.5)
    


