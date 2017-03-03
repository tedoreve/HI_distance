# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 08:58:41 2016

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.tests import zmf as z

#=============================assistant code===================================
#--------------------------------plot------------------------------------------
    
def box(data,head,region,on,*args):
    '''
    plot box pixel coordinates
    '''  
    x1,y1,x2,y2 = z.coo_box(head,on)  
    if data.ndim == 2:
        result = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result = data[:,y1:y2,x1:x2]
    if data.ndim == 4:
        result = data[0,:,y1:y2,x1:x2]
        
    x1,y1,x2,y2 = z.coo_box(head,region)
    if data.ndim == 2:
        result0 = data[y1:y2,x1:x2]
    if data.ndim == 3:
        result0 = data[args[0],y1:y2,x1:x2]
    if data.ndim == 4:
        result0 = data[0,args[0],y1:y2,x1:x2]
        
    return result,result0
#------------------about spectra,rotation curve--------------------------------    
def velocity(data,head):
    '''
    return spec velocity axis
    '''
    w   = wcs.WCS(head)
    
    if head['NAXIS'] == 3:
        pix = np.linspace(1,data.shape[0],data.shape[0])
        x,y,v   = w.wcs_pix2world(0,0,pix,0)
    if head['NAXIS'] == 4:
        pix = np.linspace(1,data.shape[0],data.shape[0])
        x,y,v,s = w.wcs_pix2world(0,0,pix,0,0)
    return v

#=============================mapp=============================================
def mapp(file,region,res,index,v_seq,m,name,on):
    #read continuum data
    cont = fits.open(file[0])
    cont_head = cont[0].header
    if cont_head['NAXIS'] == 2:
        cont_data = cont[0].data[:,:]#*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
    if cont_head['NAXIS'] == 3:
        cont_data = cont[0].data[0,:,:]#*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
    if cont_head['NAXIS'] == 4:
        cont_data = cont[0].data[0,0,:,:]#*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
  
    cont.close()
    cont_head['CUNIT3'] = 'm/s'
    
    #read spectral data
    spec = fits.open(file[1])
    spec_head = spec[0].header
    if spec_head['NAXIS'] == 3:
        spec_data = spec[0].data[:,:,:]#*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
    if spec_head['NAXIS'] == 4:
        spec_data = spec[0].data[0,:,:,:]#*z.conversion(1.4,cont_head['BMAJ'],cont_head['BMIN'])
  
    spec.close()
    spec_head['CUNIT3'] = 'm/s'
    
    #get velocity dimension in 3D spectral tube
    v = velocity(spec_data,spec_head)
    if v_seq == 'H':
        v = v[::-1]/1000 #from m/s to km/s
    elif v_seq =='OH':
        v = v/1000
    
    
    #get chosen region continuum and spectral data and plot it

    
    cont_on,cont_reg  = box(cont_data,cont_head,region,on)
    spec_on,spec_reg  = box(spec_data,spec_head,region,on,87)
    plt.imshow(cont_reg,origin='lower')
    
    
    plt.show()   
    return cont_on,cont_reg

#===============================main===========================================
if __name__=='__main__':
    file     = ['../../../VGPS/VGPS_cont_MOS049.fits',    \
                '../../../VGPS/MOS_049.Tb.fits']
    region   = [49.15,49.25,-0.4,-0.3]      #region l1,l2,b1,b2
    on       = [49.20,49.21,-0.35,-0.34]    #on 
    res      = 8
    index    = 48
    v_seq    = 'H'  # H OH
    m        = 0.05  # control the scale of spectra
    name     = '1720 MHz spectra on 1440 MHz continuum'
    cont_on,cont_reg  = mapp(file,region,res,index,v_seq,m,name,on)

    


