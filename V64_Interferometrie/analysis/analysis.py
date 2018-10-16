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
def lin_model(x, A, B):
    return A * x + B


p_1, counts_1 = np.genfromtxt(abs_path('data/n_gas1.txt'), unpack = True)
p_2, counts_2 = np.genfromtxt(abs_path('data/n_gas2.txt'), unpack = True)
p_3, counts_3 = np.genfromtxt(abs_path('data/n_gas3.txt'), unpack = True)

p = [p_1, p_2, p_3]
#p = unp.uarray(p, np.full(np.shape(p), 1))
#print(p)
counts = [counts_1, counts_2, counts_3]
unp.uarray(counts, np.full(np.shape(counts), 1))

lam = Q_(632.99, 'nanometer') #wavelength
L = Q_(ufloat(100, 0.1), 'mm') #length of the gas chamber

n = lam / 2 / L * counts + 1 #index of refraction

plt.clf()
all_params = []

p_plot = np.linspace(-50, 1050, 10000)
for i in range(3):
    fig, ax = plt.subplots(1, 1)
    params, cov = curve_fit(lin_model, xdata = p[i], ydata = noms(n[i]),
                            #sigma = stds(n[i]),
                            p0 = [0.0000001, 0])
    all_params.append(correlated_values(params, cov))

    ax.plot(p[i], noms(n[i]), 'o', label = 'Messwerte')
    ax.plot(p_plot, lin_model(p_plot, *params), 'r-', label = 'Regression')
    ax.set_xlim(-50, 1050)
    ax.set_xlabel(r'Druck $p / \si{\milli \bar}$')
    ax.set_ylabel(r'Brechungsindex $n$')
    ax.legend()
    fig.tight_layout()
    fig.savefig(abs_path(f'results/druck_fit_{i + 1}.pdf'))

all_params = np.array(all_params).T
n_1000 = all_params[0] * 1000 + all_params[1]
print(n_1000[0])
for i in [1, 2, 3]:
    r.add_result(name = f'n_{i}', value = Q_(n_1000[i - 1]))



#p = [p_1, p_2, p_3]
#M = [M_1, M_2, M_3]
#p_mid = unp.uarray(np.mean(p, axis = 0), np.std(p, axis = 0))
#M_mid = unp.uarray(np.mean(M, axis = 0), np.std(M, axis = 0) / np.sqrt(3))
#plt.clf()
#plt.errorbar(x = noms(p_mid), y = noms(M_mid), xerr = stds(p_mid), yerr = stds(M_mid),
#            fmt = '.')

plt.savefig(abs_path('results/beispiel_druck_fit.pdf'))


#print(noms(n_1), n_2, n_3)
