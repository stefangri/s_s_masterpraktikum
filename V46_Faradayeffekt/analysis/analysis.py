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

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

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

#Naturkonstanten
from scipy.constants import physical_constants as const
c = const['speed of light in vacuum']
c =  Q_(c[0], c[1])

e0 = const['elementary charge']
e0 =  Q_(e0[0], e0[1])

epsilon0 = const['electric constant']
epsilon0 = Q_(epsilon0[0], epsilon0[1])

me = const['electron mass']
me = Q_(me[0], me[1])

n = 3.4 #Brechungsindex

def mass_from_slope(slope, N, B):
    return ((e0**3 * N * B / (8 * np.pi**2  * epsilon0 * c**3 * slope * n)).to('kilogram**2'))**0.5
 


# Fit des Magnetfeldes
def gauss(x, x0, A, sig):
    return A / np.sqrt(2 * np.pi * sig**2) * np.exp(- (x - x0)**2 / (2 * sig**2))

def spulenfunc(x, x0, A, R):
    return A * R**2 / (R**2 + (x - x0)**2)**(3 / 2)

def bogenminuten_to_grad(theta):
    return theta / 60    

def lin_model(x, a, b):
    return a * x + b    




z, B = np.genfromtxt(abs_path('data/magnetfeld.txt'), unpack = True)

max_B = Q_(max(B), 'millitesla')

params_poly, cov_poly = np.polyfit(z, B, deg = 4, cov = True)
poly = np.poly1d(params_poly)
params, cov = curve_fit(gauss, z, B, p0 = [110, max(B), 20])
#params, cov = curve_fit(spulenfunc, z, B, p0 = [110, max(B) / np.sqrt(20), 20])
rcParams['figure.figsize'] = 5.906, 3

fig, ax = plt.subplots(1, 1)
ax.plot(z, B, 'o', color = (128/255, 186/255, 38/255), label = 'Messwerte')
z_plot = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)

ax.plot(z_plot, gauss(z_plot, *params), 'r-')
ax.plot(z_plot, poly(z_plot), 'r-')
#ax.plot(z_plot, spulenfunc(z_plot, *params), 'r-')
ax.set_xlabel(r'$z / \si{\milli\meter}$')
ax.set_ylabel(r'$B / \si{\milli\tesla}$')
ax.legend()
fig.savefig(abs_path('results/magnetfeld.pdf'), bbox_inches='tight', pad_inches = 0)


#GaAs rein, L = 5.11mm
L_undot = Q_(5.11, 'millimeter')

lam = np.genfromtxt(abs_path('data/wellenlängen.txt'), unpack = True)
lam = Q_(lam, 'micrometer')

theta_pB1, theta_pB2, theta_nB1, theta_nB2 = np.genfromtxt(abs_path('data/ga_as_undotiert.txt'), unpack = True)

theta_pB = Q_(theta_pB1 + bogenminuten_to_grad(theta_pB2), 'degree')
theta_nB = Q_(theta_nB1 + bogenminuten_to_grad(theta_nB2), 'degree')

dif_theta_undot = 0.5 * (theta_nB - theta_pB)
dif_theta_undot_normed = (dif_theta_undot / L_undot).to('radian / millimeter')


fig, ax = plt.subplots(1, 1)
ax.plot(lam**2, dif_theta_undot_normed, 'o', label = 'Messwerte')
ax.set_xlim(ax.get_xlim())
lam_plot = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)

ax.set_xlabel(r'$\lambda^2 / \si{\micro\meter^2}$')
ax.set_ylabel(r'Normierter Winkel $\theta / \si{\radian / \milli\meter}$')
ax.legend()
fig.savefig(abs_path('results/fit_undot.pdf'), bbox_inches='tight', pad_inches = 0)




# GaAs n dotiert, N = 2.8e18/cm^^3, L = 1.296mm
L_dot_duenn = Q_(1.296, 'millimeter')
N_duenn = Q_(2.8e18, '1/cm**3')


theta_pB1, theta_pB2, theta_nB1, theta_nB2 = np.genfromtxt(abs_path('data/ga_as_dotiert_duenn.txt'), unpack = True)

theta_pB = Q_(theta_pB1 + bogenminuten_to_grad(theta_pB2), 'degree')
theta_nB = Q_(theta_nB1 + bogenminuten_to_grad(theta_nB2), 'degree')

dif_theta_dot_duenn = 0.5 * (theta_nB - theta_pB)
dif_theta_dot_duenn_normed = (dif_theta_dot_duenn / L_dot_duenn).to('radian / millimeter')

theta_duenn = dif_theta_dot_duenn_normed - dif_theta_undot_normed 

params, cov = curve_fit(lin_model, (lam**2).magnitude, theta_duenn.magnitude)
params = correlated_values(params, cov)

slope_duenn = Q_(params[0], 'radian / millimeter / micrometer**2').to('radian / millimeter**3')
offset_duenn = Q_(params[1], 'radian / millimeter')
r.add_result(name = 'slope_duenn', value = slope_duenn)
r.add_result(name = 'offset_duenn', value = offset_duenn)

mass_duenn = mass_from_slope(slope_duenn, N_duenn, max_B)
r.add_result(name = 'mass_duenn', value = mass_duenn)
r.add_result(name = 'mass_ratio_duenn', value = mass_duenn / me)



fig, ax = plt.subplots(1, 1)
ax.plot(lam**2, theta_duenn, 'o', color = (128/255, 186/255, 38/255), label = 'Messwerte')
ax.set_xlim(ax.get_xlim())
lam_plot = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
ax.plot(lam_plot, lin_model(lam_plot, *noms(params)), 'r-', label = 'Fit')

ax.set_xlabel(r'$\lambda^2 / \si{\micro\meter^2}$')
ax.set_ylabel(r'Normierter Winkel $\theta / \si{\radian / \milli\meter}$')
ax.legend()
fig.savefig(abs_path('results/fit_dot_duenn.pdf'), bbox_inches='tight', pad_inches = 0)








# GaAs n dotiert, N = 1.2e18/cm^^3, L = 1.36mm
L_dot_dick = Q_(1.36, 'millimeter')
N_dick = Q_(1.2e18, '1/cm**3')


theta_pB1, theta_pB2, theta_nB1, theta_nB2 = np.genfromtxt(abs_path('data/ga_as_dotiert_dick.txt'), unpack = True)

theta_pB = Q_(theta_pB1 + bogenminuten_to_grad(theta_pB2), 'degree')
theta_nB = Q_(theta_nB1 + bogenminuten_to_grad(theta_nB2), 'degree')

dif_theta_dot_dick = 0.5 * (theta_nB - theta_pB)
dif_theta_dot_dick_normed = (dif_theta_dot_dick / L_dot_dick).to('radian / millimeter')

theta_dick = dif_theta_dot_dick_normed - dif_theta_undot_normed 

#maske um den ausreißer zu entfernen
mask = np.arange(len(lam)) != len(lam) - 2


params, cov = curve_fit(lin_model, (lam[mask]**2).magnitude, theta_dick[mask].magnitude)
params = correlated_values(params, cov)

slope_dick = Q_(params[0], 'radian / millimeter / micrometer**2').to('radian / millimeter**3')
offset_dick = Q_(params[1], 'radian / millimeter')
r.add_result(name = 'slope_dick', value = slope_dick)
r.add_result(name = 'offset_dick', value = offset_dick)

mass_dick = mass_from_slope(slope_dick, N_dick, max_B)
r.add_result(name = 'mass_dick', value = mass_dick)
r.add_result(name = 'mass_ratio_dick', value = mass_dick / me)


fig, ax = plt.subplots(1, 1)
ax.plot(lam**2, theta_dick, 'o', color = (128/255, 186/255, 38/255), label = 'Messwerte')
ax.plot(lam[mask == False]**2, theta_dick[mask == False], 'o', label = 'Ignorierter Datenpunkt')
ax.set_xlim(ax.get_xlim())
lam_plot = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 1000)
ax.plot(lam_plot, lin_model(lam_plot, *noms(params)), 'r-', label = 'Fit')

ax.set_xlabel(r'$\lambda^2 / \si{\micro\meter^2}$')
ax.set_ylabel(r'Normierter Winkel $\theta / \si{\radian / \milli\meter}$')
ax.legend()
fig.savefig(abs_path('results/fit_dot_dick.pdf'), bbox_inches='tight', pad_inches = 0)


r.makeresults()




