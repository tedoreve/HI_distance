# -*- coding: utf-8 -*-
"""
Created on Wed May  3 22:31:54 2017

@author: tedoreve
"""

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

l1,l2,b1,b2 = region

#T=np.sum(spec_reg_co[630:],axis=0)

fig = plt.figure(figsize=(7,8.3))

ax = fig.add_subplot(111, aspect='equal')
#neg = ax.imshow(T,cmap='gray',origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
neg = ax.imshow(cont_reg,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar=fig.colorbar(neg,cax=cax)
#plt.colorbar(neg, cax=cax)

fig.subplots_adjust(top=0.99,bottom=0.1,left=0.0,right=1.0)
#fig.subplots_adjust(top=0.98,bottom=0.1,left=0.2,right=0.0)
cbar.set_label('T(K)',fontsize=14)
cbar.ax.tick_params(labelsize=14)
ax.set_xlabel('l(deg)',fontsize=14)
ax.set_ylabel('b(deg)',fontsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
#deltal = (l2-l1)/4
#xtick  = [l2,l2+deltal,l2+deltal*2,l2+deltal*3,l2+deltal*4]
#ax.set_xticklabels(xtick,fontsize=14)
#deltab = (b2-b1)/4
#ytick  = [b1,b1+deltab,b1+deltab*2,b1+deltab*3,b1+deltab*4]
#ax.set_yticklabels(ytick,fontsize=14)

#
onoff = off
ax.add_patch(patches.Rectangle((onoff[0], onoff[2]),onoff[1]-onoff[0],    \
                               onoff[3]-onoff[2],color='r',fill=False))
onoff = on
ax.add_patch(patches.Rectangle((onoff[0], onoff[2]),onoff[1]-onoff[0],    \
                               onoff[3]-onoff[2],color='r',fill=False))
#
xmajorLocator   = MultipleLocator(0.03)
xmajorFormatter = FormatStrFormatter('%2.2f')
ax.xaxis.set_major_locator(xmajorLocator)  
ax.xaxis.set_major_formatter(xmajorFormatter)
#
levels=[10,20,30,40,60]
cs=ax.contour(cont_reg,colors='white',levels=levels,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
cbar.add_lines(cs)
ax.add_patch(patches.Circle((16.79, 0.0),radius=20/3600/2,color='white',fill=True))
#ax.add_patch(patches.Circle((15.95, 0.14),radius=20/3600/2,color='white',fill=True))