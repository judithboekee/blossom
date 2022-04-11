# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 14:32:28 2021

BLOSSOM

Non-Steady State Plant Energy Balance model

@author: Judith Boekee
"""

#%% Packages

import numpy as np
import copy as cp


#%% Class

# Initialize the bud model
class model:
    def __init__(self, model_input):
        # initialize the different components of the model
        self.input = cp.deepcopy(model_input)

  
    def run(self):
        # initialize model variables
        self.init()
  
        # time integrate model 
        for self.t in range(self.tsteps):
         
            # time integrate components
            self.timestep()
  
        # delete unnecessary variables from memory
        self.exitmodel()


    def init(self):
        # assign variables from input data
        
        # initialize constants air
        self.cp         = 1005.                 # specific heat of dry air [J kg-1 K-1]
        self.kappa      = 0.4                   # Von Karman constant [-]
        self.g          = 9.81                  # gravity acceleration [m s-2]
        self.Rd         = 287.                  # gas constant for dry air [J kg-1 K-1]
        self.bolz       = 5.67e-8               # Bolzman constant [W  m -2  K -4 ]
        self.dyn_vis    = 1.718e-5              # dynamical viscocity air [N s m-2]

        # initialize non-dimension numbers
        self.Nu                = None                       # Nusselt number []
        self.Re                = None                       # Reynolds number []
        self.Gr                = None                       # Grashof number []
        self.Ri                = None                       # Richards number (Gr/Re2) []
        
        # initialize leaf characteristics 
        self.heatcap_leaf  = self.input.heatcap_leaf         # heat capacity plant [J kg-1 K-1]
        self.r_leaf        = self.input.r_leaf               # diameter plant [m]
        self.d_leaf        = self.input.d_leaf               # thickness plant [m]
        self.E_leaf        = 0.96                            # emissivity plant
        self.Tleaf         = self.input.Tleaf                # plant temperature [K]
        self.SVF           = self.input.SVF                  # skyviewfactor plant (fraction visible to sky) []
        
        self.plantpart    = self.input.plantpart                 #'leaf', 'flower' or 'branch'
        if self.plantpart == 'branch':
            self.A_leaf        = 4* np.pi * (self.r_leaf)**2     # surface area leaf [m2] 
            self.V_leaf        = 4/3 * np.pi * (self.r_leaf)**3  # volume leaf [m2]
        elif self.plantpart == 'flower':
            self.A_leaf        = 2 * np.pi * (self.r_leaf)**2    # surface area leaf [m2] 
            self.V_leaf        = self.A_leaf * self.d_leaf       # volume leaf [m2]
        else:
            self.A_leaf        = 2 * np.pi * (self.r_leaf)**2    # surface area leaf [m2] 
            self.V_leaf        = self.A_leaf * self.d_leaf       # volume leaf [m2]     
       
        self.adv           = self.input.adv
        
        # initialize mixed-layer
        self.Swin       = 0        # Incoming shortwave radiation [W m-2]
        self.Swout      = 0        # Outgoing shortwave radiation [W m-2]
        self.Eair       = 0.8      # Emmisivity air
        self.Eground    = 0.98     # Emissivity ground

        
        self.Ps         = self.input.Ps         # surface pressure [Pa]
        self.Tair       = self.input.Tair       # air temperature [K]
        self.u          = self.input.u          # initial mixed-layer u-wind speed [m s-1]
        
        self.rho        = self.Ps / (self.Rd * self.Tair) # air density [kg m-3]
        self.kin_vis    = self.dyn_vis / self.rho         # kinematische viscosity [m2 s-1]
        self.k          = 0.024                           # thermal conductivity of air [W m-1 K-1]      
        
        self.Tsky       = self.input.Tsky       # temperature clear sky
        self.Twater       = self.input.Twater   # ground water temperature [K]
        self.Tground    = self.input.Tground    # ground surface temperature [K]
        
        # initialize observations
        self.Tair_series    = self.input.Tair_series        # observations air temperature [K]
        self.u_series       = self.input.u_series           # observation wind speed [m/s]
        self.Tleaf_series   = self.input.Tleaf_series       # observations leaf temperature [K]
        self.Tground_series = self.input.Tground_series     # observations of ground temperature [K]
        self.Tsky_series    = self.input.Tsky_series        # observations of upper air temperature [K]

        # intialize time variables
        self.tsteps     = int(np.floor(self.input.runtime / self.input.dt))
        self.dt         = self.input.dt
        self.t          = 0
        
        # initialize tendency
        self.Tleaftend    = 0.                   # tendency of leaf temperature [K s-1]
     
        # initialize output
        self.out        = model_output(self.tsteps)

    def timestep(self):
        # store output before time integration
        self.T = self.Tleaf
        
        self.calc_conductance_plant()
        self.run_radiation()
        self.run_plant()
        
        self.store()
        self.integrate_plant()
        
    def calc_conductance_plant(self):        
        # calculates the radiative and convective conductence of the bud and leaf
       
        self.Gr = abs(1 / self.Tair * self.g * (2 * self.r_leaf)**3 * (self.T - self.Tair) / (self.kin_vis**2))

        if self.u > 0:
            self.Re = abs(self.u) * (self.r_leaf *2) / self.kin_vis
            
        self.Ri = abs(self.Gr / self.Re**2)
        
        
        if self.plantpart == 'branch': 
        
            if self.Ri > 10:
                print("free")
                if self.Tleaf < self.Tair:
                    self.Nu = 0.48 * self.Gr**0.25
                elif self.Tleaf > self.Tair:
                    self.Nu = 0.09 * self.Gr**0.33
            
            # Forced convection
            elif self.Ri < 0.1:
                if self.Re <= 10e3:
                    self.Nu = (0.32 + 0.51*self.Re**0.52) 
                    
                elif self.Re > 10e3:
                    self.Nu = 0.23 * self.Re**0.6 
                    
            
            # Mixed convection
            else: 
                self.Nu = ((0.23 * self.Re**0.6 )**3.5 + (0.48 * self.Gr**0.25)**3.5)**-3.5
               
        
        if self.plantpart == 'leaf' or self.plantpart == 'flower':      
        # Free convection
            if self.Ri > 10:
               
                if self.Tleaf < self.Tair:
                    self.Nu = 0.13 * self.Gr**0.33
                elif self.Tleaf > self.Tair:
                    self.Nu = 0.5 * self.Gr**0.25
            
            # Forced convection
            elif self.Ri < 0.1:
                if self.Re <= 2e4:
                    self.Nu = 0.6 * self.Re**0.5
                   
                elif self.Re > 2e4:
                    self.Nu = 0.032 * self.Re**0.8 
                   
            
            # Mixed convection
            else: 
                self.Nu = ((0.032 * self.Re**0.8 )**3.5 + (0.13 * self.Gr**0.33)**3.5)**-3.5
                
        
        
    def run_radiation(self):  
        # Calculate the radiation balance of a bud
                                   
        self.Lwin  = 0.5 * self.SVF * self.Eair * self.bolz * self.Tsky ** 4 +  0.5 * self.SVF * self.Eground * self.bolz * self.Tground ** 4. \
            + (1 - self.SVF) *  self.E_leaf * self.bolz * self.Tleaf ** 4 # we assume there is half sky , half ground 
            
        self.Lwout = self.E_leaf * self.bolz * self.Tleaf ** 4. 
          
        self.Q     = (self.Swin - self.Swout + self.Lwin - self.Lwout) 
    
    
    def run_plant(self):

        self.rC =  (self.r_leaf*2) / (self.Nu * self.k) * (self.cp * self.rho)
        
        self.H =  self.cp * self.rho * self.rC**-1 * (self.Tleaf - self.Tair)
        self.W = self.adv * (self.Tleaf - 273.16 - self.Twater)
        
        self.Tleaftend =  (self.Q - self.H - self.W) * self.A_leaf  * self.V_leaf**-1 * self.heatcap_leaf**-1 
        
        
    def integrate_plant(self):
        
        #Air temperature timeseries
        if self.Tair_series is None:
            # set values previous time step
            self.Tair    = self.Tair          
        else:
            # set values from observations
            self.Tair     = self.Tair_series[self.t]
        
        # Ground timeseries
        if self.Tground_series is None:
            # set values previous time step
            self.Tground  = self.Tground    
        else:
            # set values from observations
            self.Tground     = self.Tground_series[self.t]        
            
        # Sky timeseries
        if self.Tsky_series is None:
            # set values previous time step
            self.Tsky  = self.Tsky    
        else:
            # set values from observations
            self.Tsky = self.Tsky_series[self.t]            
            
        # Wind timeseries
        if self.u_series is None:
            # set values previous time step
            self.u     = self.u
        else:
            # set values from observations
            self.u        = self.u_series[self.t]
            

        # Leaf timeseries    
        if self.Tleaf_series is None: 
            # set values previous time step
            Tleaf0  = self.Tleaf
            # integrate equations
            self.Tleaf   = Tleaf0 + self.dt * self.Tleaftend 

        else:
            # set values from observations
            self.Tleaf    = self.Tleaf_series[self.t]



    def store(self):
            t                         = self.t
            self.out.t[t]             = t * self.dt
            self.out.Tleaf[t]         = self.Tleaf
            self.out.Tair[t]          = self.Tair

            self.out.u[t]             = self.u
            self.out.Tsky[t]          = self.Tsky
            self.out.Tground[t]       = self.Tground
            self.out.Twater[t]        = self.Twater
            
            self.out.Ri[t]            = self.Ri
            self.out.Nu[t]            = self.Nu
            self.out.Gr[t]            = self.Gr
            
            self.out.Lwin[t]          = self.Lwin
            self.out.Lwout[t]         = self.Lwout

            self.out.Q[t]             = self.Q 
            self.out.H[t]             = self.H 
            self.out.W[t]             = self.W
            
    def exitmodel(self):    
        del(self.cp)

        del(self.kappa)
        del(self.g)
        del(self.Rd)
        del(self.bolz)
        del(self.dyn_vis)
        
        del(self.Nu)
        del(self.Ri)
        del(self.Re)
        del(self.Gr)
        
        del(self.E_leaf)
        del(self.heatcap_leaf)
        del(self.d_leaf)
        del(self.A_leaf)
        del(self.V_leaf)
        del(self.r_leaf)

        del(self.SVF)
        del(self.plantpart)
    
        del(self.Swin)
        del(self.Swout)
        del(self.Lwin)
        del(self.Lwout)
        
        del(self.rC)
        
        del(self.H)
        del(self.Q)
        del(self.W)
        
        del(self.Eair)
        del(self.Eground)
        
        del(self.Ps)
        del(self.u)
        
        del(self.rho)
        del(self.kin_vis)
        del(self.k)
        
        del(self.Tleaf)
        del(self.Tair)
        del(self.Tsky)
        del(self.Tground)
        del(self.Twater)
        
        del(self.Tair_series)
        del(self.Tleaf_series)
        del(self.Tground_series)
        del(self.Tsky_series)
        del(self.u_series)
        
        del(self.t)
        del(self.dt)
        del(self.tsteps)
        
        del(self.Tleaftend)

        
   
# Class for storing bud model output data
class model_output:
    def __init__(self, tsteps):
       
        
        self.Nu           = np.zeros(tsteps)      # Nusselt number []
        self.Re           = np.zeros(tsteps)      # Reynolds number []
        self.Gr           = np.zeros(tsteps)      # Grashof number []
        self.Ri           = np.zeros(tsteps)      # Richards number (Gr/Re2) []

        self.Swin       = np.zeros(tsteps)        # incoming shortwave radiation [W m-2]
        self.Swout      = np.zeros(tsteps)        # outgoing shortwave radiation [W m-2]
        self.Lwin       = np.zeros(tsteps)        # incoming longwave radiation [W m-2]
        self.Lwout      = np.zeros(tsteps)        # outgoing longwave radiation [W m-2]
        self.Q          = np.zeros(tsteps)        # net radiation [W m-2]
        
        self.H          = np.zeros(tsteps)        # convection [W m-2]
        self.W          = np.zeros(tsteps)        # sap flow [W m-2]
        
        self.Tair       = np.zeros(tsteps)        # air temperature [K]
        self.Tground    = np.zeros(tsteps)        # ground surface temperature [K]
        self.Tsky       = np.zeros(tsteps)        # sky temperature [K]
        self.Twater     = np.zeros(tsteps)        # ground water temperature [K]
        self.u          = np.zeros(tsteps)        # initial mixed-layer u-wind speed [m s-1]
        
        self.tsteps     = np.zeros(tsteps)        # number of timestep [s]
        self.dt         = np.zeros(tsteps)        # length timestep [s]
        self.t          = np.zeros(tsteps)        # time [s]
        
        self.Tleaf      = np.zeros(tsteps)        # leaf temperature [s]

        
# Class for storing bud model input data
class model_input:
    def __init__(self):
        # general model variables
       
        self.runtime        = None  # duration of model run [s]
        self.dt             = None  # time step [s]

        self.Tleaf          = None  # initial leaf temperature [K]
        self.heatcap_leaf   = None  # heat capacity leaf [J K-1 m-2]
        self.r_leaf         = None  # radius leaf [m]
        self.d_leaf         = None  # thickness leaf [m]
        self.SVF            = 1     # skyviewfactor plant (fraction visible to sky) []
        self.plantpart      = None  # 'leaf', 'flower' or 'branch'
        
        self.Ps             = None  # initial surface pressure [Pa]
        self.Tair           = None  # initial air temperature [K]
        self.Tground        = None  # ground surface temperature [K]
        self.Tsky           = None  # temperature clear sky [K]
        self.u              = None  # initial mixed-layer u-wind speed [m s-1]
        
        self.Tair_series    = None  # observations air temperature [K]
        self.u_series       = None  # observation wind speed [m/s]
        self.Tleaf_seres    = None  # observations leaf temperature [K]
        self.Tground_series = None  # observations of ground temperature [K]
        self.Tsky_series    = None  # observations of upper air temperature [K]
        self.adv            = 0     # advection parameter alpha []

        