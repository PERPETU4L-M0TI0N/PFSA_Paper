# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 20:12:52 2018

@author: Anderw
"""

import scipy as sp
import numpy as np
from numpy import exp, sqrt, pi, log
from scipy.constants import e, N_A, epsilon_0

F = sp.constants.physical_constants['Faraday constant'][0]
R = sp.constants.R
T = 298.15         #K
M0 = 18.01528/1000      #kg/mol
V0 = M0*1000/0.997           #cm3/mol
rho0 = 0.997     #kg/cm^3
nu = 0.5
A = 1.1777  #(mol/cm^3)^(1/2)
B = 0.3291E8    #(mol/cm^3)^(1/2)/cm
eta0 = 8.90E-4   #Pa*s    viscosity of water


color_dict = {}
color_dict['H'] = (1, 1, 1)
color_dict['Li'] = (0.71,0.098,0.263)
color_dict['Na'] = (0.784,0.62,0.11)
color_dict['K'] = (0.576,0.078,0.455)
color_dict['Cs'] = (0.078,0.408,0.486)

color_dict['Fe'] = (0.702,0.098,0.278)
color_dict['Cu'] =  (0.098,0.333,0.502)
color_dict['Ni'] = (0.341,0.69,0.098)
color_dict['V(III)'] = (0.078,0.553,0.298)
color_dict['V(IV)'] = (0.133,0.188,0.537)
color_dict['V(V)'] = (183/255, 178/255, 37/225)
color_dict['Ca'] = (0.9556, 0.493, 0.1422)

symbol_dict = {}
symbol_dict['H'] = 'o'
symbol_dict['Li'] = '^'
symbol_dict['Na'] = 's'
symbol_dict['K'] = 'p'
symbol_dict['Cs'] = 'D'
symbol_dict['Fe'] = 'P'
symbol_dict['Cu'] = '<'
symbol_dict['Ni'] = 'v'
symbol_dict['V(III)'] = 'p'
symbol_dict['V(IV)'] = '8'
symbol_dict['V(V)'] = 'd'
symbol_dict['Ca'] = '>'

#line_dict = {}
#line_dict['H'] = 'o'
#line_dict['Li'] = '^'
#line_dict['Na'] = 's'
#line_dict['K'] = 'p'
#line_dict['Cs'] = 'D'
#line_dict['Fe3'] = '8'
#line_dict['Cu2'] = '>'
#line_dict['Ni2'] = '<'
#line_dict['VIII'] = '+'
#line_dict['VIV'] = 'x'

          
