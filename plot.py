# -*- coding: utf-8 -*-
"""
Created on Wed May  3 22:31:54 2017

@author: tedoreve
"""

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as patches

l1,l2,b1,b2 = region

T=np.sum(spec_reg_co[140:146],axis=0)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
#neg = ax.imshow(T,cmap='gray',origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
neg = ax.imshow(cont_reg,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])
cbar=fig.colorbar(neg,ax=ax)

fig.subplots_adjust(top=0.98,bottom=0.1,left=0.07,right=1.0)
cbar.set_label('T(K)')
ax.set_xlabel('l(deg)')
ax.set_ylabel('b(deg)')
#deltal = (l2-l1)/4
#xtick  = [l2,l2+deltal,l2+deltal*2,l2+deltal*3,l2+deltal*4]
#ax.set_xticklabels(xtick,fontsize=20)
#deltab = (b2-b1)/4
#ytick  = [b1,b1+deltab,b1+deltab*2,b1+deltab*3,b1+deltab*4]
#ax.set_yticklabels(ytick,fontsize=20)


onoff = on
ax.add_patch(patches.Rectangle((onoff[0], onoff[2]),onoff[1]-onoff[0],    \
                               onoff[3]-onoff[2],color='r',fill=False))

xmajorLocator   = MultipleLocator(0.03)
xmajorFormatter = FormatStrFormatter('%2.2f')
ax.xaxis.set_major_locator(xmajorLocator)  
ax.xaxis.set_major_formatter(xmajorFormatter)

#levels=[10,20,30,40,60]
#plt.contour(cont_reg,colors='white',levels=levels,origin='lower',interpolation='nearest',extent=[l2,l1,b1,b2])