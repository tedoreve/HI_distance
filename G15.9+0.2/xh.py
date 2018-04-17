# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:41:31 2018

@author: tedoreve
"""
import numpy as np
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
from astropy.io import fits
import zmf as z
#============================ continuum =======================================
def continuum(file,region,on,contrast):
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
    
    cont_on,cont_reg  = z.box(cont_data,cont_head,1,region,on)

    return cont_reg,cont_on
#============================ spectra =========================================
def spectra(file,region,on,contrast,spec_v):
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

    spec_on,spec_reg  = z.box(spec_data,spec_head,1,region,on,spec_v)

    
    return spec_reg,spec_on,v

if __name__=='__main__':
    file1   = '../../data/THOR_cont_1440MHz_L16.25deg_25arcsec_image.fits'
    file2   = '../../data/grs-16-cube.fits'
    region  = [15.7,16.0,0.0,0.3]      #region l1,l2,b1,b2
    on      = [15.9,15.92,0.16,0.19] 
    contrast = 1
    spec_v   = 87
    cont_reg,cont_on   = continuum(file1,region,on,contrast)
    spec_reg_co,spec_on_co,v_co  = spectra(file2,region,on,contrast,122)
    
    l1,l2,b1,b2 = region
    T=np.sum(spec_reg_co[120:126],axis=0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    #neg = ax.imshow(T,cmap='gray',origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    neg = ax.imshow(T,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
    cbar=fig.colorbar(neg,ax=ax)
    
    fig.subplots_adjust(top=0.98,bottom=0.1,left=0.07,right=1.0)
    cbar.set_label('T(K)')
    ax.set_xlabel('l(deg)')
    ax.set_ylabel('b(deg)')
    levels=[10,20,30,40,60]
    plt.contour(cont_reg,colors='white',levels=levels,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
