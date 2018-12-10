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
from pandas import DataFrame

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
import scipy.constants as const


#Konstanten und Funktionen

#Molmasse
M_cu = Q_(63.546 * const.N_A, 'dalton / mol').to('kg / mol')
# Probenmasse Versuchsanleitung
m_probe = Q_(342, 'g')

#Kompressionsmodul
kappa = Q_(140, 'gigapascal')

#Molvolumen
V_0 = Q_(7.11 / 10**(6), 'm**3 / mol')


#Gaskonstanstante 
R_gas = Q_(const.gas_constant, 'joule / kelvin / mol')

def C_p(U, I, dt, dT):
    return M_cu / m_probe * U * I * dt / dT

def C_v(C_p, alpha, T):
    return C_p - 9 * alpha**2 * kappa * V_0 * T    


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
    return func_alpha(T, *params_alpha_units) * 1e-6



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

def R_T(T, a, b, c):
    return -b / (2 * a) + np.sqrt((b / (2 * a))**2 - c / a + T / a)  

params_T, cov_T = curve_fit(T_R, R, T)

correlated_params = correlated_values(params_T, cov_T)
params_T_units = []
for i in range(len(correlated_params)):
    params_T_units.append(Q_(correlated_params[i], f'kelvin / ohm**{2 - i}'))
    r.add_result(f'param_T{i}', params_T_units[i])

def model_T(R):
    return T_R(R, *params_T_units)

def model_R(T):
    return R_T(T, *params_T) 


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


# Berechnungen mit den Messdaten 
R_zyl, U, I, R_prob = np.genfromtxt(abs_path('data/data.txt'), unpack = True)
R_zyl = Q_(R_zyl, 'kiloohm')
R_prob = Q_(R_prob, 'kiloohm')
U = Q_(U, 'volt')
I = Q_(I, 'milliampere')
t = np.array([300 * i for i in range(9900 // 300)])

T_zyl = model_T(R_zyl).to('kelvin')
T_prob = model_T(R_prob).to('kelvin')

mean_diff = np.mean(T_prob - T_zyl)
r.add_result('mean_T_diff', mean_diff)

fig, ax = plt.subplots(1, 2)

ax[0].plot(t, noms(T_prob), '.', label = r'$T_{\text{Prob}}$')
ax[0].plot(t, noms(T_zyl), '.', label = r'$T_{\text{Zyl}}$')
ax[0].set_ylabel('$T / \si{\kelvin}$')
ax[0].set_xlabel('$t / \si{\second}$')
ax[0].legend()

ax[1].plot(t, noms(T_prob - T_zyl), '.', color = 'r')
ax[1].set_ylabel(r'$(T_{\text{Prob}} - T_{\text{Zyl}}) / \si{\kelvin}$')
ax[1].set_xlabel('$t / s$')
ax[1].set_ylim(-6, 6)


fig.tight_layout()
fig.savefig(abs_path('results/temperaturen.pdf'), bbox_inches='tight', pad_inches = 0)










#Bestimmung der Debye-Temperatur

C_p = C_p(U[:-1], I[:-1], Q_(np.diff(t), 'second'), Q_(np.diff(T_prob), 'kelvin'))
C_v = C_v(C_p, model_alpha(T_prob[1:]), T_prob[1:]).to('joule / kelvin / mol')

rcParams['figure.figsize'] = 5.906, 4.5
fig, ax = plt.subplots(1, 1)
ax.plot(noms(T_prob[1:]), noms(C_v), 'o', label = 'Daten')




def c_v_debye(T, T_D):
    return [9 * R_gas.magnitude * (T / T_D)**3 * (quad(lambda x: x**4 * np.exp(x) / (np.exp(x) - 1)**2, 0, T_D / T))[0] for T in T]

params, cov = curve_fit(c_v_debye, noms(T_prob[1:]), noms(C_v), p0 = [300])
T_plot = np.linspace(60, 350, 1000)

ax.set_xlim(T_plot[0], T_plot[-1])
ax.plot(T_plot, c_v_debye(T_plot, *params), 'r-', label = 'Fit')
ax.axhline(y = 3 * R_gas.magnitude, color = 'k', linestyle = '--', label = '$3 R$')
ax.plot(T_plot, c_v_debye(T_plot, 345), 'g-', label = 'Theorie')
ax.set_xlabel(r'$T / \si{\kelvin}$')
ax.set_ylabel(r'$C_V / \si{\joule  \per \kelvin \per \mol}$')
ax.legend()

fig.tight_layout()
fig.savefig(abs_path('results/C_V.pdf'), bbox_inches='tight', pad_inches = 0)
