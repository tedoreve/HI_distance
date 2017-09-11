# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:49:43 2017

@author: tedoreve
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

# Generate fake data
np.random.seed(42)
g1 = models.Gaussian1D(0.4, 63, 1)
g2 = models.Gaussian1D(0.4, 43, 1)
g3 = models.Gaussian1D(0.2, 28, 1)
x = np.linspace(-50, 150, 1000)
y = g1(x) + g2(x) + g3(x) + np.random.normal(0., 0.01, x.shape)
x = v_co/1000
y = T_on_co
# Now to fit the data create a new superposition with initial
# guesses for the parameters:
gg_init = models.Gaussian1D(0.4, 63, 1) + models.Gaussian1D(0.4, 43, 1) + models.Gaussian1D(0.2, 28, 1)
#gg_init = models.Gaussian1D(1, 0, 0.1) + models.Gaussian1D(2, 0.5, 0.1)
fitter = fitting.SLSQPLSQFitter()
gg_fit = fitter(gg_init, x, y)

# Plot the data with the best-fit model
#plt.figure(figsize=(8,5))
plt.plot(x, y, 'ko')
plt.plot(x, gg_fit(x))
plt.xlabel('Position')
plt.ylabel('Flux')

