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
        result0 = data[:,y1:y2,x1:x2] 
    result0 = (result0-np.mean(result0))+np.mean(result0)*contrast
    
#    plt.subplots() 
#    plt.title(name)
#    plt.imshow(result0,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
#    plt.colorbar()
#    plt.xlabel(r'$l (deg)$')
#    plt.ylabel(r'$b (deg)$')    
#    
#    Path = mpath.Path
#    path_data = [
#        (Path.MOVETO, (onoff[0], onoff[2])),
#        (Path.CURVE4, (onoff[0], onoff[3])),
#        (Path.CURVE4, (onoff[1], onoff[3])),
#        (Path.CURVE4, (onoff[1], onoff[2])),
#        (Path.CLOSEPOLY, (onoff[0], onoff[2])),   
#        ]
#    codes, verts = zip(*path_data)
#    path = Path(verts, codes)
#    x, y = zip(*path.vertices)
#    plt.plot(x, y, 'go-')
#    plt.xlim(l2,l1)
#    plt.ylim(b1,b2)
#    plt.grid()
    
    x1,y1,x2,y2 = z.coo_box(head,onoff)  
    if data.ndim == 2:
        result = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result = data[:,y1:y2,x1:x2]
    return result,result0


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

def con(head,cont_on,region,levels):
    '''
    plot box pixel coordinates
    '''
    l1,l2,b1,b2 = region  
    x1,y1,x2,y2 = z.coo_box(head,region)
    
    plt.subplots() 
    plt.contour(cont_on,levels=levels,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    plt.colorbar()
    plt.xlabel(r'$l (deg)$')
    plt.ylabel(r'$b (deg)$')    
    
    plt.xlim(l2,l1)
    plt.ylim(b1,b2)
    plt.grid()
    
    
    
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
    if cont_head['BUNIT'] == 'K' and cont_head['NAXIS'] == 4:
        cont_data = cont[0].data[0,0,:,:]
    elif cont_head['BUNIT'] == 'K' and cont_head['NAXIS'] == 3:
        cont_data = cont[0].data[0,:,:]
    else:
        cont_data = cont[0].data[0,:,:]*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
    cont.close()
    cont_head['CUNIT3'] = 'm/s'
    
    z.plot_fits(cont_data,cont_head,1,'origin')
    
    if analyze == 'box':
        if on != []:
            cont_on,cont_reg  = box(cont_data,cont_head,1,'cont_on',region,on)
        if off != []:
            cont_off,cont_reg = box(cont_data,cont_head,1,'cont_off',region,off)
    elif analyze == 'circle':
            cont_on,cont_reg  = circle(cont_data,cont_head,1,'1420MHz continuum (K)',region,on)
            cont_off,cont_reg = circle(cont_data,cont_head,1,'cont_off',region,off)
    return cont_on,cont_off

#=============================spectra==========================================
def spectra(file,analyze,region,on,off,contrast,spec_v):
    spec = fits.open(file)
    spec_head = spec[0].header
    if spec_head['BUNIT'] == 'K' and spec_head['NAXIS'] == 4:
        spec_data = spec[0].data[0,:,:,:]
    elif spec_head['BUNIT'] == 'K' and spec_head['NAXIS'] == 3:
        spec_data = spec[0].data[:,:,:]
    else:
        spec_data = spec[0].data[:,:,:]*z.conversion(1.4,spec_head['BMAJ'],spec_head['BMIN'])
    spec.close()
    spec_head['CUNIT3'] = 'm/s'
        
    v = velocity(spec_data,spec_head)
        
    if analyze == 'box':
        if on != []:
            spec_on,spec_reg  = box(spec_data,spec_head,1,'spec_on',region,on,spec_v)
        if off != []:    
            spec_off,spec_reg = box(spec_data,spec_head,1,'spec_off',region,off,spec_v)
    elif analyze == 'circle':
        if on != []:
            spec_on,spec_reg  = circle(spec_data,spec_head,1,'1720MHz spectrum map (K) at '+str(int(v[spec_v]))+' m/s',region,on,spec_v)
        if off != []:    
            spec_off,spec_reg = circle(spec_data,spec_head,1,'spec_off',region,off,spec_v)
    
    return spec_on,spec_off,v
#=============================absorption=======================================
def absorption_spec(spec_on,spec_off,v,spec_on_co,spec_off_co,v_co,v0,d0,cont_on,cont_off,on,off,analyze,method): 
    if analyze   == 'box':
        if method == 'tww':
            T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
            T_off    = (np.sum(np.sum(spec_off,axis=1),axis=1)                    \
                        -(np.sum(np.sum(spec_on,axis=1),axis=1)))                 \
                        /(spec_off.shape[1]*spec_off.shape[2]-spec_on.shape[1]*spec_on.shape[2])
            
            T_con    = np.mean(cont_on)
            T_coff   = (np.sum(cont_off)-np.sum(cont_on))                         \
                        /(cont_off.shape[0]*cont_off.shape[1]-cont_on.shape[0]*cont_on.shape[1])
                        
            T_on_co  = np.mean(np.mean(spec_on_co,axis=1),axis=1)
#            T_on_vgps= np.mean(np.mean(spec_on_vgps-cont_on_vgps,axis=1),axis=1)
#            T_off_co = (np.sum(np.sum(spec_off_co,axis=1),axis=1)                    \
#                        -(np.sum(np.sum(spec_on_co,axis=1),axis=1)))                 \
#                        /(spec_off_co.shape[1]*spec_off_co.shape[2]-spec_on_co.shape[1]*spec_on_co.shape[2])

            e_tau    = 1+(T_on-T_off)/(T_con-T_coff)       
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
 
    v     = v/1000
    v_co  = v_co/1000
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    x1 = ax1.plot(v , T_on )
    x2 = ax1.plot(v , T_off )
    ax1.plot([138]*len(list(range(-80,20))),list(range(-80,20)),'--',color='purple')
#    ax4= ax1.twinx()
#    x4 = ax4.plot(v_vgps,T_on_vgps,color='r')
   
    xx = x1 + x2
    labs = ['HI_on','HI_off']
    props = font_manager.FontProperties(size=10)
    ax1.legend(xx, labs, loc='lower right', shadow=True, prop=props)
    
    ax1.set_ylabel('T(K)')
#    ax4.set_ylabel('T(K)')
    ax1.set_title('Spectrum of G15.9+0.2')
#    ax1.set_ylim(y2lim[0],y2lim[1])

    x1 = ax2.plot(v , e_tau )
    ax2.plot(v ,[1]*len(v ),'--',color='purple')
    ax2.plot([138]*len(np.arange(-0.25,1.5,0.05)),np.arange(-0.25,1.5,0.05),'--',color='purple')
    ax22  = ax2.twinx()
    x2 = ax22.plot(v_co, T_on_co,color='r')
    
    xx = x1 + x2
    labs = ['HI','CO']
    ax2.legend(xx, labs, loc='lower right', shadow=True, prop=props)
    
    ax2.set_ylabel(r'$e^{-\tau}$',fontsize=15) 
    ax22.set_ylabel('T(K)')
    
    ax3.plot(v0,d0)
    ax3.plot(v[0:177] ,[7.5]*len(v[0:177]),'--',color='purple')
    ax3.plot([138]*len(list(range(0,20))),list(range(0,20)),'--',color='purple')
    ax3.set_xlabel('velocity(km/s)')
    ax3.set_ylabel('distance(kpc)')
    ax3.set_xlim(v[75],v[199])
    ax3.set_ylim(0,20)
#    ax3.tick_params('y')
    
#    ax2.plot(v_co[70:], T_off_co[70:],label='off')

#    ax2.set_ylim(ylim[0],ylim[1])
    fig.subplots_adjust(hspace=0.05)
#    plt.legend()
    plt.show()
    return v,e_tau

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
    file1   = '../data/THOR_cont_1440MHz_L16.25deg_25arcsec_image.fits'
    file2   = '../data/THOR_HI_without_continuum_L16.25.fits'
    file3   = '../data/grs-16-cube.fits'
    file4   = '../data/MOS_017.Tb.fits'
    file5   = '../data/VGPS_cont_MOS017.fits'
    file6   = '../data/rotation_model.txt'
    region  = [15.7,16.0,0.0,0.3]      #region l1,l2,b1,b2
    on      = [15.9,15.92,0.16,0.19] 
    off     = [15.9,15.94,0.16,0.19]
    on_co   = [15.88,15.94,0.14,0.2]
    contrast = 1
    analyze  = 'box'               # box,circle
    spec_v   = 87
    model   = 'constant'            #constant, model
    V       = 220                   #km/s
    d       = np.linspace(1,40,100)
    l       = 15.9
    b       = 0.2
    method  = 'tww' #获得吸收谱的方法，tww或者classic
    levels=[5,10,20,60,100,140]
    cont_on,cont_off    = continuum(file1,analyze,region,on,off,contrast)
    spec_on,spec_off,v  = spectra(file2,analyze,region,on,off,contrast,spec_v)
    spec_on_co,spec_off_co,v_co  = spectra(file3,analyze,region,on_co,off,contrast,122)
    #v_co[120:125] T=np.sum(spec_on_co[120:125],axis=0)
    v0,d0 = dist(model,file6,l,b,d,V = 220,v_sun = 220,r_sun = 8.5)
#    spec_on_vgps,spec_off_vgps,v_vgps  = spectra(file4,analyze,region,on,off,contrast,spec_v)
#    cont_on_vgps,cont_off_vgps    = continuum(file5,analyze,region,on,off,contrast)
#    v,e_tau = absorption_spec(spec_on,spec_off,v,spec_on_co,spec_off_co,v_co,v0,d0,cont_on,cont_off,on,off,analyze,method)
    
    


