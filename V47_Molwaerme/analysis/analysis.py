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
from scipy.interpolate import interp1d
from tab2tex import make_table

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
rcParams['axes.formatter.use_mathtext'] = True
rcParams['legend.fontsize'] = 10
rcParams['savefig.dpi'] = 300
import scipy.constants as const


#Konstanten und Funktionen

#Molmasse
M_cu = Q_(63.546 * const.N_A, 'dalton / mol').to('kg / mol') #https://www.webqc.org/molecular-weight-of-Cu.html
# Probenmasse Versuchsanleitung
m_probe = Q_(342, 'g')

#Kompressionsmodul
kappa = Q_(14e10, 'pascal') #gross

#Molvolumen
V_0 = Q_(7.0922e-6, 'm**3 / mol') #http://periodictable.com/Elements/029/data.html

#theoriewert
#T_debye_lit = Q_(343, 'kelvin')

N_L = m_probe / M_cu * const.N_A * u('1 / mol')
V = V_0 * m_probe / M_cu
v_long = Q_(4.7, 'km / s').to('m / s')
v_trans = Q_(2.26, 'km / s').to('m / s')

omega_debye = (18 * np.pi**2 * N_L / V / (1 / v_long**3 + 2 / v_trans**3))**(1 / 3)

r.add_result(name = 'omega_debye', value = omega_debye)

T_debye_theorie = Q_(const.hbar, 'J * s') * omega_debye / Q_(const.k, 'J / kelvin')
r.add_result(name = 'T_debye_theorie', value = T_debye_theorie)


#Gaskonstanstante 
R_gas = Q_(const.gas_constant, 'joule / kelvin / mol')

def C_p(U, I, dt, dT):
    return M_cu / m_probe * U * I * dt / dT

def C_v(C_p, alpha, T):
    return C_p - 9 * alpha**2 * kappa * V_0 * T     


#Fits für Ausdehnungskoeffizient und T-R-Charakteristik
T, alpha = np.genfromtxt(abs_path('data/alpha.txt'), unpack = True)
model_alpha = interp1d(T, alpha  * 1e-6)

#def model_alpha(Temp):
#    T, alpha = np.genfromtxt(abs_path('data/alpha.txt'), unpack = True)
#    return np.interp(Temp, T, alpha)

#def func_alpha(T, a, b, c, d):
#    return a * T**3 + b * T**2 + c * T + d
#
#params_alpha, cov_alpha = curve_fit(func_alpha, T, alpha)
#
#correlated_params = correlated_values(params_alpha, cov_alpha)
#
#params_alpha_units = []
#for i in range(len(correlated_params)):
#    params_alpha_units.append(Q_(correlated_params[i], f'1 / kelvin**{4 - i}'))
#    r.add_result(f'param_alpha{i}', params_alpha_units[i])


rcParams['figure.figsize'] = 5.906, 3
fig, ax = plt.subplots(1, 2)

T_plot = np.linspace(T[0], T[-1], 1000)
ax[1].plot(T, alpha, '.', label = 'Daten')
ax[1].plot(T_plot, model_alpha(T_plot) * 1e6, 'r-', label = 'Interpolation')
ax[1].set_xlabel('$T / \si{\kelvin}$')
ax[1].set_ylabel(r'$\alpha / \si{10^{-6} \kelvin}^{-1}$')
ax[1].set_xlim(T_plot[0] - 10, T_plot[-1] + 10)
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
R_zyl = Q_(unp.uarray(R_zyl, np.repeat(0.0001, len(R_zyl))), 'kiloohm')
R_prob = Q_(unp.uarray(R_prob, np.repeat(0.0001, len(R_prob))), 'kiloohm')
U = Q_(unp.uarray(U, np.repeat(0.01, len(U))), 'volt')
I = Q_(unp.uarray(I, np.repeat(0.1, len(I))), 'milliampere')
t = Q_(unp.uarray([300 * i for i in range(9900 // 300)], np.repeat(1, 9900 // 300)), 'second')


T_zyl = model_T(R_zyl).to('kelvin')


T_prob = model_T(R_prob).to('kelvin')

make_table(filename = abs_path('tabs/data.tex'),
    data = [t.magnitude, U.magnitude, I.magnitude, R_zyl.magnitude, T_zyl.magnitude, R_prob.magnitude, T_prob.magnitude], 
    header = [r'$t$ / \second', r'$U$ / \volt', r'$I$ / \milli\ampere', r'$R_{\text{zyl}}$ / \kilo\ohm', r'$T_{\text{zyl}}$ / \kelvin', r'$R_{\text{prob}}$ / \kilo\ohm', r'$T_{\text{prob}}$ / \kelvin'],
    places = [(4.0, 1.0), (2.2, 1.2), (3.1, 1.1), (1.4, 1.4), (3.2, 1.2), (1.4, 1.4), (3.2, 1.2)],
    caption = r'Aufgenommene Messdaten zur Bestimmung der Debye Temperatur von Kupfer. Teit $t$, anliegende Spannung $U$ an der Probe, Strom $I$ durch die Probe, Widerstand $R_{\text{zyl}}$ des Zylinder-Pt-Elements und Widerstand $R_{\text{probe}}$ des Proben-Pt-Elements. Aus den Widerständen werden die Temperaturen $T_{\text{zyl}}$ und $T_{\text{probe}}$ berechnet.',
    label = 'tab: data'
) 

T_mean = (T_prob[1:] + T_prob[:-1]) / 2


fig, ax = plt.subplots(1, 2)

ax[0].plot(noms(t), noms(T_prob), '.', label = r'$T_{\text{Prob}}$')
ax[0].plot(noms(t), noms(T_zyl), '.', label = r'$T_{\text{Zyl}}$')
ax[0].set_ylabel('$T / \si{\kelvin}$')
ax[0].set_xlabel('$t / \si{\second}$')
ax[0].legend()

ax[1].plot(noms(t), noms(T_prob - T_zyl), '.', color = 'b')
ax[1].plot(noms(t)[1:3], noms(T_prob - T_zyl)[1:3], '.', color = 'r')
ax[1].set_ylabel(r'$(T_{\text{Prob}} - T_{\text{Zyl}}) / \si{\kelvin}$')
ax[1].set_xlabel('$t / s$')
ax[1].set_ylim(-6, 6)


fig.tight_layout()
fig.savefig(abs_path('results/temperaturen.pdf'), bbox_inches='tight', pad_inches = 0)








#Bestimmung der Debye-Temperatur

C_p = C_p(U[:-1], I[:-1], Q_(np.diff(t), 'second'), Q_(np.diff(T_prob), 'kelvin'))
alpha = Q_(model_alpha(noms(T_mean[:-1])), '1 / kelvin')

C_v = C_v(C_p[:-1], alpha, T_mean[:-1]).to('joule / kelvin / mol')






rcParams['figure.figsize'] = 5.906, 6.2
fig, ax = plt.subplots(1, 1)
ax.errorbar(x = noms(T_mean[:-1]), y = noms(C_v), xerr = stds(T_mean[:-1]), yerr = stds(C_v), fmt = '.', label = 'Daten', linestyle = None)
ax.errorbar(x = noms(T_mean[1]), y = noms(C_v)[1], xerr = stds(T_mean[1]), yerr = stds(C_v)[1], fmt = '.', label = 'Ausreißer', linestyle = None, color = 'r')




def c_v_debye(T, T_D):
    return [9 * R_gas.magnitude * (T / T_D)**3 * (quad(lambda x: x**4 * np.exp(x) / (np.exp(x) - 1)**2, 0, T_D / T))[0] for T in T]

#params, cov = curve_fit(c_v_debye, noms((T_prob[1:] + T_prob[:-1]) / 2), noms(C_v), p0 = [300])
#print(correlated_values(params, cov))
T_plot = np.linspace(60, 350, 1000)

ax.set_xlim(T_plot[0], T_plot[-1])
#ax.plot(T_plot, c_v_debye(T_plot, *params), 'r-', label = 'Fit')
ax.axhline(y = 3 * R_gas.magnitude, color = 'k', linestyle = '--', label = '$3 R$')
ax.plot(T_plot, c_v_debye(T_plot, T_debye_theorie.magnitude), 'g-', label = 'Theorie')
ax.set_xlabel(r'$\overline{T} / \si{\kelvin}$')
ax.set_ylabel(r'$C_V / \si{\joule  \per \kelvin \per \mol}$')



def weight_mean(X):
    n = noms(X)
    s = stds(X)
    mean = sum(n * 1 / s**2) / sum(1 / s**2)
    std = np.sqrt(1 / sum(1 / s**2))
    return ufloat(mean, std)

# Nochmal mit anderer Methode

data = np.genfromtxt(abs_path('data/tab.txt'), unpack = True) 
data = data.T
best_theta_D = []
best_theta_std = []
mask_25 = noms(C_v) < 40
for C_v_i in C_v[mask_25].magnitude:
    index = np.unravel_index((abs(data - C_v_i.n)).argmin(), data.shape)
    theta = index[0] + 0.1 * index[1]
    best_theta_D.append(theta)
    index_upper = np.unravel_index((abs(data - C_v_i.n - C_v_i.s)).argmin(), data.shape)
    index_lower = np.unravel_index((abs(data - C_v_i.n + C_v_i.s)).argmin(), data.shape)
    theta_upper = index_upper[0] + 0.1 * index_upper[1]
    theta_lower = index_lower[0] + 0.1 * index_lower[1]
    best_theta_std.append(max(abs(theta - theta_lower), abs(theta - theta_upper)))

best_theta_D = unp.uarray(best_theta_D, best_theta_std)
theta_D = best_theta_D * (T_mean[:-1])[mask_25]


make_table(filename = abs_path('tabs/results.tex'),
    data = [np.diff(t)[:-1], np.diff(T_prob)[:-1], C_p.magnitude[:-1], T_mean.magnitude[:-1], alpha.magnitude * 1e6, C_v.magnitude, best_theta_D, theta_D.magnitude], 
    header = [r'$\Delta t$ / \second', r'$\Delta T$ / \kelvin', r'$C_p$ / \joule\per\kelvin\per\mol', r'$\overline{T}$ / \kelvin', r'$\alpha(\overline{T})$ / 10^{-6}\per\kelvin', r'$C_V$ / \joule\per\kelvin\per\mol', r'$\frac{\Theta}{\overline{T}}$', r'$\Theta_D$ / \kelvin'],
    places = [(3.1, 1.1), (1.2, 1.2), (2.1, 1.1), (3.2, 1.2), 2.2, (2.1, 1.1), (1.1, 1.1), (3.0, 3.0)],
    caption = r'Ergebnisse für die spezifische Wärmekapazität bei konstantem Druck und konstantem Volumen, sowie der Debye-Temperatur.',
    label = 'tab: results'
) 

theta_D_mean = weight_mean(theta_D)   
r.add_result(name = 'T_debye_exp', value = Q_(theta_D_mean, 'kelvin'))  

ax.plot(T_plot, c_v_debye(T_plot, theta_D_mean.n), 'r-', label = 'Experiment')
ax.fill_between(T_plot, c_v_debye(T_plot, theta_D_mean.n + theta_D_mean.s), c_v_debye(T_plot, theta_D_mean.n - theta_D_mean.s), color='r', alpha=0.3)
from matplotlib.patches import Rectangle
rect = Rectangle((0,0), 0,0, color='r', alpha=0.3, label = '$1\sigma$')
ax.add_patch(rect)

ax.legend()


fig.tight_layout()
fig.savefig(abs_path('results/C_V.pdf'), bbox_inches='tight', pad_inches = 0)





