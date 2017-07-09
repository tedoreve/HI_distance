# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 08:58:41 2016

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
#matplotlib.use('Agg')
from astropy.io import fits
import zmf as z
 
#============================ continuum =======================================
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
    
    z.plot_2Dfits(cont_data,cont_head,1)
    
    if analyze == 'box':
        if on != []:
            cont_on,cont_reg  = z.box(cont_data,cont_head,1,region,on)
        if off != []:
            cont_off,cont_reg = z.box(cont_data,cont_head,1,region,off)
    elif analyze == 'circle':
            cont_on,cont_reg  = z.circle(cont_data,cont_head,1,region,on)
            cont_off,cont_reg = z.circle(cont_data,cont_head,1,region,off)
#    z.con(cont_reg,cont_head,region,levels=[5,10,20,60,100,140])
    return cont_on,cont_off,cont_reg

#============================ spectra =========================================
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
        
    v = z.velocity(spec_data,spec_head)
        
    if analyze == 'box':
        if on != []:
            spec_on,spec_reg  = z.box(spec_data,spec_head,1,region,on,spec_v)
        if off != []:    
            spec_off,spec_reg = z.box(spec_data,spec_head,1,region,off,spec_v)
    elif analyze == 'circle':
        if on != []:
            spec_on,spec_reg  = z.circle(spec_data,spec_head,1,region,on,spec_v)
        if off != []:    
            spec_off,spec_reg = z.circle(spec_data,spec_head,1,region,off,spec_v)
    
    return spec_on,spec_off,v
#============================ absorption ======================================
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
            e_tau    = 1+(T_on-T_off)/(T_con-T_coff)       
        else:
            T_on     = np.mean(np.mean(spec_on,axis=1),axis=1)
            T_off    = np.mean(np.mean(spec_off,axis=1),axis=1)
            T_con    = np.mean(cont_on)
            T_coff   = np.mean(cont_off)
            e_tau    = (T_on-T_off)/(T_con-T_coff)       
            
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
    
    #开始画图
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    x1 = ax1.plot(v , T_on )
    x2 = ax1.plot(v , T_off )
    ax1.plot([138]*len(list(range(-80,20))),list(range(-80,20)),'--',color='purple')
    ax1.plot([21]*len(list(range(-80,20))),list(range(-80,20)),'--',color='purple')
#    ax4= ax1.twinx()
#    x4 = ax4.plot(v_vgps,T_on_vgps,color='r')   
    xx = x1 + x2
    labs = ['HI_on','HI_off']
    props = font_manager.FontProperties(size=10)
    ax1.legend(xx, labs, loc='lower right', shadow=True, prop=props)    
    ax1.set_ylabel('T(K)')
#    ax4.set_ylabel('T(K)')
    ax1.set_title('Spectra of G15.9+0.2')
#    ax1.set_ylim(y2lim[0],y2lim[1])
    x1 = ax2.plot(v , e_tau )
    ax2.plot(v ,[1]*len(v ),'--',color='purple')
    ax2.plot([138]*len(np.arange(-0.25,1.5,0.05)),np.arange(-0.25,1.5,0.05),'--',color='purple')
    ax2.plot([21]*len(np.arange(-0.25,1.5,0.05)),np.arange(-0.25,1.5,0.05),'--',color='purple')
    ax22  = ax2.twinx()
    x2 = ax22.plot(v_co, T_on_co,color='r')
    
    xx = x1 + x2
    labs = ['HI','CO']
    ax2.legend(xx, labs, loc='lower right', shadow=True, prop=props)    
    ax2.set_ylabel(r'$e^{-\tau}$',fontsize=15) 
    ax22.set_ylabel('T(K)')    
    
    ax3.plot(v0,d0)
    ax3.plot(v[0:177] ,[7.5]*len(v[0:177]),'--',color='purple')
    ax3.plot(v[0:99] ,[14]*len(v[0:99]),'--',color='purple')
    ax3.plot([138]*len(list(range(0,20))),list(range(0,20)),'--',color='purple')
    ax3.plot([21]*len(list(range(0,20))),list(range(0,20)),'--',color='purple')
    ax3.set_xlabel('velocity(km/s)')
    ax3.set_ylabel('distance(kpc)')
    ax3.set_xlim(0,170)
    ax3.set_ylim(0,20)
    
    
    fig.subplots_adjust(hspace=0.05)
    plt.show()
    
    return T_on_co,e_tau
#===============================main===========================================
if __name__=='__main__':
    file1   = '../data/THOR_cont_1440MHz_L16.25deg_25arcsec_image.fits'
    file2   = '../data/THOR_HI_without_continuum_L16.25.fits'
    file3   = '../../data/grs-16-cube.fits'
    file4   = '../../data/MOS_017.Tb.fits'
    file5   = '../../data/VGPS_cont_MOS017.fits'
    file6   = '../data/rotation_model.txt'
    region  = [15.8,16.1,0.0,0.3]      #region l1,l2,b1,b2
    on      = [15.83,15.95,0.1,0.25] 
    off     = [15.83,16.1,0.1,0.25]
    on_co   = [15.88,15.94,0.14,0.2]
    contrast = 1
    analyze  = 'box'               # box,circle
    spec_v   = 87
    model   = 'constant'            #constant, model
    V       = 220                   #km/s
    d       = np.linspace(1,40,100)
    l       = 15.9
    b       = 0.2
    method  = 'tww' #the way to get absorption spectrum，'tww' or 'classic'
    levels=[5,10,20,60,100,140]
    cont_on,cont_off,cont_reg    = continuum(file5,analyze,region,on,off,contrast)
    spec_on,spec_off,v  = spectra(file4,analyze,region,on,off,contrast,spec_v)
    spec_on_co,spec_off_co,v_co  = spectra(file3,analyze,region,on,off,contrast,122)
    #v_co[120:125] T=np.sum(spec_on_co[120:125],axis=0)
    v0,d0 = z.dist(l,b,d,V = 240,v_sun = 240,r_sun = 8.34)
#    spec_on_vgps,spec_off_vgps,v_vgps  = spectra(file4,analyze,region,on,off,contrast,spec_v)
#    cont_on_vgps,cont_off_vgps    = continuum(file5,analyze,region,on,off,contrast)
    T_on_co,e_tau = absorption_spec(spec_on,spec_off,v,spec_on_co,spec_off_co,v_co,v0,d0,cont_on,cont_off,on,off,analyze,method)
    
