#!/usr/bin/python
#Filename: estra_spec; this is a script to extra spectra automately with given source size

import  sys
import  math
import  pyfits

from    sys         import  stdout
import  numpy       as      np
from    astropy.io  import  fits
from    astropy.io  import  ascii

from    astropy     import  wcs
from    astropy.wcs import  WCS
import  wcsaxes

from spectral_cube import SpectralCube
from    astropy.coordinates import SkyCoord
from    astropy             import units as u
from    astropy.table       import Table, Column
from    matplotlib.patches  import Rectangle, Circle, Ellipse
from    astropy.io.fits     import getdata
from    matplotlib.lines    import Line2D
from    astropy             import units as u
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

from wcsaxes import datasets
from astropy.stats import sigma_clipped_stats
from astropy import coordinates

import matplotlib.patches as patches

import  matplotlib.pyplot   as plt
import  matplotlib.image    as mpimg
import  pylab               as py
import  aplpy
from astropy.nddata import Cutout2D

from matplotlib.ticker import FormatStrFormatter
import montage_wrapper as montage


CRED = '\033[91m'  ##print word in red
CEND = '\033[0m'

def extra_data(name, infits, l_deg, b_deg, source_size,outfile):
    '''
        l_deg: longititude, degree;
        b_deg: latitude, degree;
        infits: input fits name; 
        v_lsr: local standard velocity - 20; km/s;
        source_size: the extract spectral region; arcseconds;
        #vlsr: local stardard of rest velocity from literature, km/s;
        '''
    
    #print 'begin cut the spectral cube '
    stdout.write("\t name:%20s In_fits:%20s l:%15.8f b:%15.8f Source_size:%15.8f Outfile:%20s\n" %(name, infits, l_deg, b_deg, source_size, outfile))
    cube        = SpectralCube.read(infits)
    vlsr        = cube.spectral_axis.value  #return quantity array without unit and just value
    print 'the max and min velocity',np.max(vlsr),np.min(vlsr)
    channels = cube.shape[0]
    moment0   = cube.moment(order=0)  #obtain the 2d data header information.
    header    = moment0.header
    w        = wcs.WCS(header)
    if source_size < 19: #beam size
        source_size = 19
    else:
        source_size = source_size
    pixelsize   = header['CDELT2']  #read the pixel size in unit of degree
    pixelscale  = ((source_size)/ (60*60*pixelsize)) #cut the pixel of 1.0arcmin, image size 2arcmin*2arcmin
    print 'pixel size in deg',pixelsize, 'source size in arcsec', source_size,'cutpixel size number',pixelscale, 'veloctiy channels',channels

    ########obtain the pixel coordinates for the extract region#########
    pixcrd      = w.wcs_world2pix([[l_deg, b_deg]],1)
    # convert world coordinate to pixel coordinate, if deal with 3d data, then should be [l,b,0], but sometimes will come out wrong if you deal with 3d data, so better change to 2d;
    #print 'pixel coord',pixcrd
    mypixa      = pixcrd[0,0]
    mypixb      = pixcrd[0,1]
#crdpix      = w.wcs_pix2world([[mypixa,mypixb]],1) #for double check the world coordinates;
    print 'pixel coord',mypixa, mypixb
#print 'the check world coor', crdpix
    
    xlow        = (mypixa - 0.5*pixelscale)
    xupper      = (mypixa + 0.5*pixelscale)
    #print xlow, xupper
    if xlow < 0:
        xlow = 0
    else:
        xlow = xlow
    
    ylow        = (mypixb - 0.5*pixelscale)
    yupper      = (mypixb + 0.5*pixelscale)
#print ylow, yupper
    if ylow < 0:
        ylow = 0
    else:
        ylow = ylow
    print 'extra min and max pixel coor respect to x and y',xlow, xupper,ylow,yupper
##############set the temperature#############
    data = open(outfile,'w+')
    print >> data, '#V_LSR','T_a' # 'T_a2', 'T_a3'
    for i in xrange(0, channels):
        #subcube1 = subcube[i,211:217,226:232]  #for check the result with ds9
        subcube = cube[i,ylow:yupper,xlow:xupper] # this is the right way after check with ds9
        #subcube2 = Cutout2D(subcube[i,:,:], center, size=[source_size,source_size]*u.arcsec,wcs=w)
        #print subcube2, subcube.spectral_axis[-1]
        #print((subcube2.input_position_original, subcube2.input_position_cutout))
        #Cutout2D result is not quite right comppare to find the pixel coordinates.
        print >> data, vlsr[i], np.mean(subcube) #np.mean(subcube2) #np.mean(subcube3.data)
    print(CRED+'extra the source '+name+' with a size of '+'%.1f'%source_size+' arcsec successfullt'+CEND)


