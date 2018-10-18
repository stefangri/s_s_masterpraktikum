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


l.Latexdocument(filename = abs_path('tabs/kontrast_tabelle.tex')).tabular(
    data = [angle * 180 / np.pi, U_max, U_min, K * 100], 
    header = [r'\text{Winkel} / \degree', r'U_{max} / \volt', r'U_{min} / \volt', r'\text{Kontrast} / \percent'],
    places = [0, 2, 2, 1],
    caption = 'Winkel des Polarisationsfilter, sowie maximal und minimal gemessene Spannung zur Bestimmung des Kontrast des Interferometers.',
    label = 'messwerte_kontrast'
)   



# DRUCKABHÄNGIGKEIT
def lin_model(x, A, B):
    return A * x + B


p_1, counts_1 = np.genfromtxt(abs_path('data/n_gas1.txt'), unpack = True)
p_2, counts_2 = np.genfromtxt(abs_path('data/n_gas2.txt'), unpack = True)
p_3, counts_3 = np.genfromtxt(abs_path('data/n_gas3.txt'), unpack = True)

p = [p_1, p_2, p_3]
counts = [counts_1, counts_2, counts_3]
unp.uarray(counts, np.full(np.shape(counts), 1))

lam = Q_(632.99, 'nanometer') #wavelength
L = Q_(ufloat(100, 0.1), 'mm') #length of the gas chamber

n = lam / 2 / L * counts + 1 #index of refraction
n = unp.uarray(noms(n), stds(n))
#print(n)

l.Latexdocument(filename = abs_path('tabs/gas_tabelle.tex')).tabular(
    data = [p_1, counts_1, p_2, counts_2, p_3, counts_3], 
    header = [r'p_1 / \milli\bar', r'\text{Counts}_1 / ', r'p_2 / \milli\bar', 
            r'\text{Counts}_2 / ', r'p_3 / \milli\bar', r'\text{Counts}_3 / ', 
            ],
    places = [0, 0, 0, 0, 0, 0],
    caption = r'Messwerte der drei Reihen zur Bestimmung des Brechungsindex von Luft. Hierhi bezeichnet $p_i$ den gemessenen Druck und $\text{Counts}_i$ die gezählte Anzahl an $2\pi$ Phasenverschiebungen.',
    label = 'messwerte_n_gas'
) 

l.Latexdocument(filename = abs_path('tabs/n_gas_tabelle.tex')).tabular(
    data = [n[0], n[1], n[2]], 
    header = [ r'n_1 / ', r'n_2 / ', r'n_3 / '],
    places = [(1.8, 1.8), (1.8, 1.8), (1.8, 1.8)],
    caption = r'Berechnete Brechungsindices aus den Werten der Tabelle~\ref{tab: messwerte_n_gas}.',
    label = 'berechnete_n_gas'
) 

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
    fig.savefig(abs_path(f'results/druck_fit_{i + 1}.pdf'), bbox_inches = 'tight', pad_inches = 0)

all_params = np.array(all_params).T
all_params = unp.uarray(noms(all_params), stds(all_params))
n_1000 = all_params[0] * 1000 + all_params[1]
n_1000 = unp.uarray(noms(n_1000), stds(n_1000))
n_1000_mean = ufloat(np.mean(noms(n_1000)), np.std(noms(n_1000)))
r.add_result(name = r'n', value = Q_(n_1000_mean))

#for i in [1, 2, 3]:
#    r.add_result(name = f'n_{i}', value = Q_(n_1000[i - 1]))

l.Latexdocument(filename = abs_path('tabs/fitparams_druck.tex')).tabular(
    data = [[1, 2, 3], unp.uarray(noms(all_params[0]*1e7), stds(all_params[0]*1e7)), all_params[1], n_1000], 
    header = ['Messung / ', r'A / 10^{-7}\milli\bar^-1', 'B / ', r'n_\mathup{norm} / '],
    places = [0, (1.2, 1.2), (1.7, 1.7), (1.7, 1.7)],
    caption = 'Ermittelte Regressionsparameter $A$ und $B$ der Messung zur Bestimmung der Abhängigkeit zwischen Brechungsindex $n$ und Druck p. Zudem ist jeweils der Wert $n_{\mathup{norm}}$ zum Vergleich mit der Literatur angegeben.',
    label = 'fitparams_druck'
)    




#Brechungsindex Glas
T = Q_(1, 'millimeter')

def n_glas(counts, theta):
    alpha = counts * lam / 4 / T 
    return (((alpha**2 + 2 *  (1 - np.cos(theta)) * (1 - alpha)) / ( 2 * (1 - np.cos(theta) - alpha))).to('dimensionless')).magnitude

def n_glas2(counts, theta1, theta2):
    return ((1 / (1 - counts * lam /  T / (theta1**2 - theta2**2))).to('dimensionless')).magnitude


def n_glas3(counts, deltheta, theta0 = 10 / 180 * np.pi):
    return ((1 / (1 - counts * lam  /T / (theta0**2 + 2 * theta0 * deltheta))).to('dimensionless')).magnitude   

#def n_glas3(counts, theta):
#    return 

counts_glas = np.genfromtxt(abs_path('data/n_glas.txt')) 

#print((n_glas(counts_glas, 10 / 180 * np.pi)))
#print(n_glas2(counts_glas, 20 / 180 * np.pi, 10 / 180 * np.pi))
n_glas = n_glas3(counts = counts_glas, deltheta = 10 / 180 * np.pi)
print(n_glas)

l.Latexdocument(filename = abs_path('tabs/counts_glas.tex')).tabular(
    data = [counts_glas, n_glas], 
    header = [r'\text{Counts} / ', r'n / '],
    places = [0, 2],
    caption = 'Gemessene Anzahl der $2\pi$ Phasenverschiebungen (Counts) unter Drehung der Glasplatten um $\SI{10}{\degree}$, sowie daraus berechnete Brechungsindices $n$.',
    label = 'counts_glas'
)  


n_glas = ufloat(np.mean(n_glas), np.std(n_glas, ddof = 1))
print(n_glas)
r.add_result(name = 'n_glas', value = Q_(n_glas))




