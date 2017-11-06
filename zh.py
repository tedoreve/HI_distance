from astropy.io import fits
import copy
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
import zmf as z
import numpy as np

hdu1 = fits.open('../data/Herschel500um.fits')[1]


h = copy.deepcopy(hdu1.header)
c = SkyCoord(h['CRVAL1'],h['CRVAL2'], frame='icrs', unit='deg')
h['CRVAL1'],h['CRVAL2'] = c.galactic.l.value, c.galactic.b.value
h['CTYPE1'],h['CTYPE2'] = 'GLON-CAR', 'GLAT-CAR'

array, footprint = reproject_interp(hdu1, h)
region  = [206.5,208,-2.5,-1.3]
x1,y1,x2,y2 = z.coo_box(h,region)

ax2 = plt.subplot(1,1,1, projection=WCS(h))
ax2.imshow(array, origin='lower')
ax2.coords.grid(color='white')
ax2.coords['glon'].set_axislabel('Galactic Longitude')
ax2.coords['glat'].set_axislabel('Galactic Latitude')
ax2.set_xlim(x1,x2)
ax2.set_ylim(y1,y2)
ax2.set_title('l-b')
#

star = np.loadtxt('../data/rosStar.txt')
w = WCS(h)
x = []
y = []
for i in range(30):
    x1,y1 = w.wcs_world2pix(star[i][2],star[i][3],0)
    x.append(x1)
    y.append(y1)
plt.plot(x,y,'o',color='r')