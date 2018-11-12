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
 


# Fit des Magnetfeldes
def gauss(x, x0, A, sig):
    return A / np.sqrt(2 * np.pi * sig**2) * np.exp(- (x - x0)**2 / (2 * sig**2))

def spulenfunc(x, x0, A, R):
    return A * R**2 / (R**2 + (x - x0)**2)**(3 / 2)




z, B = np.genfromtxt(abs_path('data/magnetfeld.txt'), unpack = True)

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