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



#KONTRAST MESSUNG


def contrast_model(phi, a, phi_0, const):
    return a * abs(np.sin(2 * phi + phi_0)) + const

angle, U_max, U_min = np.genfromtxt(abs_path('data/kontrast.txt'), unpack = True)


K = (U_max - U_min) / (U_max + U_min)
angle = angle / 180 * np.pi #convert to rad

params, cov = curve_fit(contrast_model, angle, K, p0 = [0.5, -np.pi / 8, 0.5])

rcParams['figure.figsize'] =  5.906, 3
fig, ax = plt.subplots(1, 1)
ang_plot = np.linspace(-1, 4, 1000)
ax.set_xlim(-0.1, np.pi + 0.1)
ax.plot(ang_plot, contrast_model(ang_plot, *params), 'r-', label = 'Fit')
#ax.plot(ang_plot, abs(np.sin(2 * ang_plot)), 'g-', label = r'$\left|\sin(2\vartheta)\right|$')
ax.plot(angle, K, 'o', label = 'Messdaten')
ax.set_xlabel(r'$\vartheta / \si{\radian}$')
ax.set_ylabel(r'Kontrast $K$')
ax.legend()
fig.tight_layout()
fig.savefig(abs_path('results/kontrast_fit.pdf'))


params = correlated_values(params, cov)
r.add_result(name = 'A', value = Q_(params[0]))
r.add_result(name = '\\vartheta_0', value = Q_(params[1] / np.pi * 180, 'degree'))
print(params[1])
r.add_result(name = r'C', value = Q_(params[2]))



# DRUCKABHÄNGIGKEIT
p_1, M_1 = np.genfromtxt(abs_path('data/n_gas1.txt'), unpack = True)
p_2, M_2 = np.genfromtxt(abs_path('data/n_gas2.txt'), unpack = True)
p_3, M_3 = np.genfromtxt(abs_path('data/n_gas3.txt'), unpack = True)

lam = Q_(633, 'nanometer')
L = Q_(ufloat(100, 0.1), 'mm')
n_1 = lam / L / 2 * M_1 + 1
n_2 = lam / L / 2 * M_2 + 1
n_3 = lam / L / 2 * M_3 + 1

print(noms(n_1), n_2, n_3)





fig.savefig(abs_path('results/druck_fits.pdf'))   
