import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values
import math
from scipy.optimize import curve_fit
from pint import UnitRegistry
import latex as l
r = l.Latexdocument('results.tex')
u = UnitRegistry()
Q_ = u.Quantity
import os
from scipy.integrate import quad
x = 2

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

def noise(y):  
    error_y = np.random.normal(0, 0.01, size=len(y))
    return y + error_y

import matplotlib as mlp
mlp.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = '[lmroman10-regular]:+tlig;'
rcParams['text.usetex'] = True
rcParams['pgf.texsystem'] = 'lualatex'
rcParams['font.size'] = 10
rcParams['mathtext.fontset'] = 'custom'
rcParams['pgf.preamble'] = r'\usepackage[locale=DE]{siunitx}'#\usepackage{tgpagella}
rcParams['text.latex.preamble'] = r'\usepackage[math-style=ISO,bold-style=ISO,sans-style=italic,nabla=upright,partial=upright,]{unicode-math}'
#rcParams['text.latex.preamble'] = r'\usepackage[locale=DE]{siunitx}'
#rcParams['text.latex.preamble'] = r'\DeclareSIUnit{\counts}{Counts}'
rcParams['axes.formatter.use_mathtext'] = True
rcParams['legend.fontsize'] = 10
rcParams['savefig.dpi'] = 300
from siunitx import *


#Naturkonstanten


#Fits für Ausdehnungskoeffizient und T-R-Charakteristik
T, alpha = np.genfromtxt(abs_path('data/alpha.txt'), unpack = True)


def func_alpha(T, a, b, c, d):
    return a * T**3 + b * T**2 + c * T + d

params_alpha, cov_alpha = curve_fit(func_alpha, T, alpha)

correlated_params = correlated_values(params_alpha, cov_alpha)

params_alpha_units = []
for i in range(len(correlated_params)):
    params_alpha_units.append(Q_(correlated_params[i], f'1 / kelvin**{4 - i}'))
    r.add_result(f'param_alpha{i}', params_alpha_units[i])

def model_alpha(T):
    return func_alpha(T, *correlated_params) * 1e-6



rcParams['figure.figsize'] = 5.906, 3
fig, ax = plt.subplots(1, 2)

T_plot = np.linspace(60, 310, 1000)
ax[1].plot(T, alpha, '.', label = 'Daten')
ax[1].plot(T_plot, func_alpha(T_plot, *params_alpha), 'r-', label = 'Fit')
ax[1].set_xlabel('$T / \si{\kelvin}$')
ax[1].set_ylabel(r'$\alpha / \si{10^{-6} \kelvin}^{-1}$')
ax[1].set_xlim(T_plot[0], T_plot[-1])
ax[1].legend(loc = 'lower right')
ax[1].text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax[1].transAxes)



T, R = np.genfromtxt(abs_path('data/T_R.txt'), unpack = True)
T += 273.15 # to kelvin

def T_R(R, a, b, c):
    return a * R**2 + b * R + c

params_T, cov_T = curve_fit(T_R, R, T)

correlated_params = correlated_values(params_T, cov_T)
params_T_units = []
for i in range(len(correlated_params)):
    params_T_units.append(Q_(correlated_params[i], f'kelvin / ohm**{2 - i}'))
    r.add_result(f'param_T{i}', params_T_units[i])

def model_T(R):
    return T_R(R, *params_T_units) 


R_plot = np.linspace(10, 120, 1000)
ax[0].plot(R, T, '.', label = 'Daten')
ax[0].plot(R_plot, T_R(R_plot, *params_T), 'r-', label = 'Fit')
ax[0].set_xlabel('$R / \si{\ohm}$')
ax[0].legend(loc = 'lower right')
ax[0].set_ylabel('$T / \si{\kelvin}$')
ax[0].set_xlim(R_plot[0], R_plot[-1])
ax[0].text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax[0].transAxes)

fig.tight_layout()
fig.savefig(abs_path('results/alpha_T_R.pdf'), bbox_inches='tight', pad_inches = 0)












#Bestimmung der Debye-Temperatur

def c_v_debye(T, T_D):
    return [(T / T_D)**3 * (quad(lambda x: x**4 * np.exp(x) / (np.exp(x) - 1)**2, 0, T_D / T))[0] for T in T]

T = np.linspace(80, 170, 100)
C_V = c_v_debye(T, 300)
C_V = noise(C_V)
params, cov = curve_fit(c_v_debye, T, C_V)

fig, ax = plt.subplots(1, 1)
ax.plot(T, C_V, '.') 
ax.plot(T, c_v_debye(T, *params), 'r-')
ax.set_xlabel('T / \si{\kelvin}')
#fig.savefig('test.pdf', bbox_inches='tight', pad_inches = 0)   
