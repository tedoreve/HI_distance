# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 18:05:07 2018

@author: tedoreve
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
#fig.subplots_adjust(top=0.98,bottom=0.02,left=0.02,right=0.98)
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

#circle1 = plt.Circle((0, 0), 8.5, fill=False , label='Solar Circle')
#circle2 = plt.Circle((0, 8.5), 0.5,color='r', fill=False )
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
        
r = 8.5
x = np.linspace(-r,r,1000)
y = (r**2-x**2)**0.5
ax.plot(x,y,color='black',label='Solar Circle')
ax.plot(x,-y,color='black')

x = [1.0,3.0]
y = [2.0,-7.9]
ax.scatter(0, 8.5, s=200, c='r', alpha=0.5,label='Sun')
ax.scatter(0, 0, s=200, c='orange', alpha=0.5,label='Galactic Center')
ax.scatter(x, y, marker='x', s=100, c='b', alpha=0.5,label='G15.9+0.2')
#ax.scatter(3.0, -7.9,marker='x', s=100, c='b', alpha=0.5,label='G15.9+0.2')
ax.scatter(3.5, -5.5, marker='*',s=100, c='g', alpha=0.5,label='G16.7+0.1')

ax.arrow(1.0, 2.0, 1.8, -8.9, head_width=0.5, head_length=1,color='black')
ax.arrow(2.8, -6.9, -1.6, 7.9, head_width=0.5, head_length=1,color='black')
ax.text(0.5, 6.5, 'l=16 deg',rotation=-74 , fontdict=font)
ax.text(0.6, -1.0, 'Distance range',rotation=-76 , fontdict=font)


x1, y1 = np.array([[0.0, 0.0], [-8.5, 8.5]])
line1 = mlines.Line2D(x1, y1 , lw=2, ls='--', color='r')
x2, y2 = np.array([[0.0, 3.4], [8.5, -7.9]])
line2 = mlines.Line2D(x2, y2 , lw=2, ls='--', color='g',label='Line of Sight')
x3, y3 = np.array([[0.0, 7.83], [0.0, 3.43]])
line3 = mlines.Line2D(x3, y3 , lw=2, ls='--', color='b',label='Radius to Tangent Point')

ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)




#colors = np.linspace(0, 1, len(patches))

ax.add_line(line1)
ax.add_line(line2)
ax.add_line(line3)
#ax.add_patch(circle1)
#ax.add_patch(circle2)
ax.set_xlim((-10, 10))
ax.set_ylim((-10, 10))
#ax.legend([circle1], ['Stress State'])

plt.legend(loc='lower left', shadow=True, scatterpoints=1)
plt.show()