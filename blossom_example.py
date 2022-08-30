# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 15:57:11 2021

Example run of BLOSSOM

@author: Judith Boekee
"""


import os
os.chdir('...')

import blossom as blossom
import matplotlib.pyplot as plt


""" 
Create model_input and run case
"""
run1input                = blossom.model_input()

run1input.model          = 'RV'

run1input.dt             = 1       # time step [s]
run1input.runtime        = 300     # total run time [s]

run1input.Tair           = 278      # initial mixed-layer potential temperature [K]
run1input.Tground        = 280      # initial ground surface temperature [K]
run1input.Tsky           = 260      # temperature clear sky [K]
run1input.u              = 4        # initial mixed-layer u-wind speed [m s-1] 
run1input.Ps             = 101300   # initial surface pressure [Pa]
run1input.Twater         = 280

run1input.plantpart      = 'leaf'
run1input.heatcap_leaf   = 1000e3     # heat capacity plant  [J m-3 K-1]
run1input.d_leaf         = 0.0007     # thickness plant [m]
run1input.r_leaf         = 0.00657    # radius plant organ [m]
run1input.Tleaf          = 275        # initial plant temperature [K]
run1input.SVF            = 0.8

run1input.Tair_series    = None     # observations air temperature [K]
run1input.Tleaf_series   = None     # observations leaf temperature [K]
run1input.u_series       = None     # observations wind speed [m/s]
run1input.Tground_series = None     # observations soil temperature [K]
run1input.Tsky_series    = None     # observations sky temperature [K]


r1 = blossom.model(run1input)
r1.run()

""" 
Create plots to show model outcomes
"""

fig, ax = plt.subplots(1,1)
ax.plot(r1.out.t, r1.out.Tair-273.16, label = 'air temperature', ls = '--')
ax.plot(r1.out.t, r1.out.Tleaf-273.16, label = 'leaf temperature')

ax.set_xlabel('time [s]')
ax.set_ylabel('T [degC]')

ax.legend()

