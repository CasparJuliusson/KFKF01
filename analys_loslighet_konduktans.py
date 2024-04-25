from math import *
from numpy import *
from linregmc import *
import matplotlib.pyplot as plt

# Script for analysis of conductivity experiment.
# Below, everything spelled NN should be numbers and everything spelled 
# EXPRESSION should be some kind of Python code. 

nmc = 1000       # Number of Monte-Carlo points
T = 293         # Temperature

pconf = 0.95  # the desired confidence level, for example 0.95  (or 0.683 which gives the standard error)

# Conductivites, molar conductivities and dielectricity constant
# Order: butanol etanol metanol vatten
kappa = array([7.56*10**-4, 236.3*10**-4, 0.95, 24.03])     # S/m
Lam = array([0.13, 0.66, 1.6, 3.2])       # S/m/M
er = array([17.10, 24.30, 32.63, 78.54])

# Concentrations, their logarithms and dG0 (in J/mol)
# Remember that ln([Na][Cl]) = 2*ln([salt])
R = 8.3145
c = kappa/Lam            # M
lnc = log(c)
dG0 = -2*R*T*lnc     # J/mol

## Generate MC concentrations and errors
# Relative standard deviation for conductivity. Use value below or values 
# given by supervisor
sigrel = 0.3

#To avoid negative mc-concentrations, we use the lognormal
#distribution. 
c_mc = createmcdata(c,sigrel*c,nmc,'lognorm')
lnc_mc = log(c_mc) 
dG0_mc = -2*R*T*lnc_mc 
dG0_sig = std(dG0_mc,0)  #Standard errors of the dG0 values, to be used for error bars. The second argument 0 means take standard deviation over the 0th dimension (rows).


## Perform linear fit with MC. 

# The function linregmc takes two arguments:
# xvalues - a 1d-array with length N (the number of data points)
# yvalues - a 2d-array with dimensions (nmc x N)

# Remember that pp contains the coefficients of the polynomial that you fit to the
# data (compare to POLYFIT, which you might be familiar with).

x_val = 1/er -1

(pp,psig,pchi2,pmc) = linregmc(x_val, dG0_mc)
(perr, plow, phigh) = mcerrconf(pmc, 0.95)


# Plot data and result of fit. This is the main result graph of the
# lab. It hould contain: ∆G0 with errorbars as a function of (1/eps_r-1). 
# It should also contain the line that resulted from the
# linear regression above. 
# Use plt.errorbar(xvalues, yvalues, yerr=yerrors) instead of plt.plot to plot error bars.
# Then add your fitted line by using polyval function on your polynomial pp.
# For nicer y axis, please convert to kJ/mol

plt.plot(x_val, polyval(pp/1000, x_val))
plt.errorbar(x_val, dG0/1000, yerr=dG0_sig/1000, fmt = 'none', capsize = 4, ecolor= 'black')
plt.xlabel("1/er - 1")
plt.ylabel("kJ/mol")
plt.show()

## Calculate radius and do error propagation to get error in ion radius

# constant before (1/a+-) in expression
k = 96485**2/(4*pi*8.854e-12*6.022e23)   

# Calculate a from constant k and determined slope pp[0] from linregmc.
# Calculate gitter energy from intercept pp[1] from linregmc.
# Also convert to nicer units: Å and kJ/mol
a = k/pp[0] * pow(10, 10)
dGgitter = pp[1]/1000

# Calculate the confidence intervals for the radius and gitter energy by first computing 
# their values for all the MC fits.
a_mc = k/pmc[:,0] * pow(10, 10)          #Use same formula as you did to calculate a, now using pmc[:,0] instead of pp[0]
dGgitter_mc = pmc[:,1]/1000


#Perform statistical analysis of these values to get the error and coverage at the pconf level 
(aerr,alow,ahigh) = mcerrconf(a_mc, pconf)
(dGerr,dGlow,dGhigh) = mcerrconf(dGgitter_mc, pconf)    


print('a = %.2f +/- %.2f Å (pconf = %.2f)\n'%(a,aerr,pconf))
print('dG0gitter = %.1f +/- %.1f kJ/mol (pconf = %.2f)\n'%(dGgitter,dGerr,pconf))

#Calculate a+- from literature values of ionic radii
a_plus = 1.16 # for Na+
a_minus = 1.67 # for Cl-
a_lit = 2/(1/a_plus + 1/a_minus)
print("Literature: %.2f Å "%a_lit)

print(pchi2)