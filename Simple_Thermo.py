
import scipy as sp  #scientific computer package
from scipy import optimize
import numpy as np   #math and array package
import pandas as pd
import matplotlib.pyplot as plt
from numpy import exp, sqrt, pi, log  #special functions and constants

from constants import F, R, T, M0, rho0, V0, A, B #from the constants.py file in the same folder

""" universal constants """
R = 8.314 #gas constant, SI unit


""" operating conditions """
T = 298 #temperature, K
def Global_Temp_set(exp_T):
    global T 
    T= exp_T

""" water material properties """
M0 = 18.01528/1000      #kg/mol, molar mass of water
rho0 = 0.997     #kg/cm^3, density of water
V0 = 18e-6 #m^3/mol, partial molar volume of water


""" polymer material properties """
Vp = (1100/2.1)*10**-6 #m^3/mol, partial molar volume of polymer


""" Charge numbers of ions """
zM = -1.0
zH = 1.0
zs = np.array([zH,zM])


"""parameters: swelling chemical potential"""
# NON-FITTED
phi_domain0 = 40.94/(1100/2.1) #volume fraction of hydrophilic polymer (sulfonic acid sites) 
d0 = 2.7E-7 # dry ionomer domain spacing
n_scaling = 2 # ionomer geometric parameter that is 1, 2, and 3 for lamellar, cylindrical, and spherical domain
m_scaling = 1.33 # ionomer morphology scaling parameter, between 1 and 1.5

"""parameters: charging, Debye-Huckel chemical potential"""
#NON-FITTED PARAMETERS    
a = 3E-8 #Table II main paper, same as a_i
b = 0.47E-9 #Table II main paper
A = 1.1777 #List of symbols in main paper
B = 3.291E9 #List of symbols in main paper
#RELEVANT SUB_FUNCTIONS
def sigma_fun(x):  #debye Huckel function
    """sigma function"""
    return 3*(x - 2*log(1 + x) - 1/(1 + x) + 1)/x**3

   
"""parameters: solvation energy chemical potential"""
kH = 1.13839062 #supplementary material table C1
kM = 0.0 #from transport conductivity code
ks = np.array([kH, kM])

ai = 3E-8

#From supplementary material appendix A
NH = 5*zH
NM = 4*zM
Ns = np.array([NH, NM])

"""
exp_data = pd.read_excel('./sorption.xlsx', sheet_name="H-M")  #importing experimental data
RH_exp = exp_data.RH.values
lam_exp = exp_data.lam.values
"""

""" External Phase Chemical Potential """
def mu_external(xs_ext):
    ## Ideal
    mu0_id = R*T*log(xs_ext)
    
    ## Total external phase chemical potential
    mu0_ext = mu0_id
    return mu0_ext
""" Membrane Phase Chemical Potential """

def mu_membrane(lam,E0,beta_HM):
    ## Ideal
    nu = 2 #total stoichiometric coefficient in the membrane
    # mole fractions of species
    x0 = lam/(lam+nu)
    xM = (1-x0)/2
    xH = xM #equal to membrane because of electroneutrality
    
    xs = np.array([xH,xM]) #array with mole fractions except water 
    
    mu0_id = R*T*log(x0)   

    ##########################################################################

    ###########################   Swelling   #################################
    phi0 = lam/(lam+Vp/V0) #volume fraction of water
    doverd0 = (1-phi0)**(-m_scaling) #ratio of hydrated to dry domain spacing
    phi_domain = phi0 + (1-phi0)*phi_domain0 #hydrophilic hydrated ionomer volume fraction (water volume fraction + hydropholic ionomer fraction)
    mu0_swe = V0*E0*(1-doverd0*((1-phi_domain**(1/n_scaling))/(1-phi_domain0**(1/n_scaling))))
    ##########################################################################
   
    
     ###########################   Physical   #################################
    # molalities of species
    mM = 1/(lam*M0) 
    mH = mM    
    ms = np.array([mH, mM])   
    betas_mem = np.array([[0, beta_HM],[beta_HM, 0]])
    mu0_phy = - M0*np.sum(betas_mem.dot(ms)*ms)*R*T    
    ##########################################################################
    
    
    
    ##############################    Total   ################################
    mu0_mem =  mu0_id + mu0_phy +mu0_swe
    ##########################################################################
    
    return mu0_mem

##############################################################################
##############################################################################

""" Condition of Phase Equilibrium for Water """
def phase_equib(membrane_phase_composition,RH,E0,beta_HM):
    lam = membrane_phase_composition    
    mu_ext = mu_external(RH)
    # print(mu_ext)
    mu_mem = mu_membrane(lam,E0,beta_HM)
    res = mu_ext - mu_mem
    return res

##############################################################################
##############################################################################

def Res_PhaseEqb(params,RH_exp): 
    lam_pred = [] # to store lam for each xs_ext in xs_ext_list
    for RH in RH_exp:
        lam_0 = 2
        #solution = sp.optimize.least_squares(phase_equib, lam_0, args=(RH,params[0],params[1]), bounds = [0,20], method ='trf',loss='soft_l1',f_scale=1)
        solution = sp.optimize.root(phase_equib, lam_0, args=(RH,params[0],params[1]), method="broyden1") #solving for internal condition of phase equilibrium!!!!
        lam = solution.x  
        lam_pred = lam_pred + [lam]     
        # print(lam_pred)
    return lam_pred


# def Res_Params(params):
#     lam_pred = np.array(Res_PhaseEqb(params,RH_exp))  
#     res = np.reshape(lam_exp,(lam_exp.shape[0],1)) - lam_pred
#     return res.flatten()

def Res_Params_CurveFit1(RH_exp, E0, beta_HM):
    lam_pred = np.array(Res_PhaseEqb((E0,beta_HM),RH_exp))  
    # res = np.reshape(lam_exp,(lam_exp.shape[0],1)) - lam_pred
    return lam_pred


### DATA FITTING

# Initial guesses for fitting parameters
param00 = 210e6 #E0
param10 = 0.133 #beta_HM

init_guess = np.array([param00,param10])


# Parameter Bounds
lb_0 = 150e6 #for param[0] (E0)
ub_0 = 400e6

lb_1 = 0.04 #for param[1] (beta_HM)
ub_1 = 0.2

bounds_array = np.array([[lb_0,lb_1],[ub_0,ub_1]])

"""
# Solver
# sol = sp.optimize.least_squares(Res_Params, init_guess, bounds = bounds_array, method ='trf')
popt,pocv = sp.optimize.curve_fit(Res_Params_CurveFit, RH_exp, lam_exp, bounds = bounds_array)

# Results
E0 = popt[0]
beta_HM = popt[1]


### PRINTING FITTED PARAMETERS
print('E0=',E0)
print('beta_HM=',beta_HM)

RH_plot = np.linspace(5e-7,1,200)
lam_pred = Res_PhaseEqb((E0,beta_HM),RH_plot)
   
plt.plot(RH_exp, lam_exp, 'o', c='k')
plt.plot(RH_plot, lam_pred, c='g')
plt.xlabel('RH')
plt.ylabel('lambda')  
plt.show()

"""
    
    