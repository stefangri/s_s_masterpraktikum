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
import pandas as pd
from scipy.signal import find_peaks


def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

#import matplotlib as mlp
#mlp.use('pgf')
#import matplotlib.pyplot as plt
#from matplotlib import rcParams

#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = '[lmroman10-regular]:+tlig;'
#rcParams['text.usetex'] = True
#rcParams['pgf.texsystem'] = 'lualatex'
#rcParams['font.size'] = 10
#rcParams['mathtext.fontset'] = 'custom'
#rcParams['pgf.preamble'] = r'\usepackage[locale=DE]{siunitx}'#\usepackage{tgpagella}
#rcParams['text.latex.preamble'] = r'\usepackage[math-style=ISO,bold-style=ISO,sans-style=italic,nabla=upright,partial=upright,]{unicode-math}'
##rcParams['text.latex.preamble'] = r'\usepackage[locale=DE]{siunitx}'
##rcParams['text.latex.preamble'] = r'\DeclareSIUnit{\counts}{Counts}'
#rcParams['axes.formatter.use_mathtext'] = True
#rcParams['legend.fontsize'] = 10
#rcParams['savefig.dpi'] = 300

import matplotlib.pyplot as plt

# #############################################################################
# ## Getting the data of the spectrum measured with only one photodiode #######
# #############################################################################

df_spectrum = pd.read_csv('./data/spectrum_angled.csv')

peak_indexes = find_peaks(-1*df_spectrum['2'].values[1:].astype('float'),
                          height=-5.3, distance=30)

peak_labels = ['87b', '85b', '85a', '87a']

plt.plot(df_spectrum['x-axis'].values[1:].astype('float'),
         df_spectrum['1'].values[1:].astype('float'), label=r'$U_{piezo}$')

plt.plot(df_spectrum['x-axis'].values[1:].astype('float'),
         df_spectrum['2'].values[1:].astype('float'), label=r'$U_{pd}$')

#plt.plot(df_spectrum['x-axis'].values[1:].astype('float')[peak_indexes[0][:-2]],
#         df_spectrum['2'].values[1:].astype('float')[peak_indexes[0][:-2]],'.',
#         label=r'$U_{pd}$')

for i in range(len(peak_indexes[0])-2):
    plt.text(df_spectrum['x-axis'].values[1:].astype('float')[peak_indexes[0][i]]- 0.0039,
             df_spectrum['2'].values[1:].astype('float')[peak_indexes[0][i]]+0.05,
             s=peak_labels[i], rotation=-20., label='Transition')

plt.xlabel(r'$t \, / \, \mathrm{s}$')
plt.ylabel(r'$U \, / \, \mathrm{V}$')
plt.legend()
#plt.show()
plt.savefig('./plots/spectrum_angled.pdf')


# #############################################################################
# ## Getting the data of the straight spectrum  #######
# #############################################################################


df_spectrum_straight = pd.read_csv('./data/spectrum_straight.csv')

peak_indexes = find_peaks(-1*df_spectrum_straight['2'].values[1:].astype('float'),
                          height=-5, distance=100)

plt.clf()
plt.plot(df_spectrum_straight['x-axis'].values[1:].astype('float'),
         df_spectrum_straight['2'].values[1:].astype('float'),
         label=r'$\Delta U_{pd}$')

#plt.plot(df_spectrum_straight['x-axis'].values[1:].astype('float')[peak_indexes[0]],
#         df_spectrum_straight['2'].values[1:].astype('float')[peak_indexes[0]],'.' ,
#         label=r'$U_{pd}$')

for i in range(len(peak_indexes[0])):
    plt.text(df_spectrum_straight['x-axis'].values[1:].astype('float')[peak_indexes[0][i]]- 0.0039,
             df_spectrum_straight['2'].values[1:].astype('float')[peak_indexes[0][i]]+0.05,
             s=peak_labels[i], rotation=-20., label='Transition')

plt.xlabel(r'$t \, / \, \mathrm{s}$')
plt.ylabel(r'$U \, / \, \mathrm{V}$')
plt.legend()
plt.savefig('./plots/spectrum_straight.pdf')


# #############################################################################
# ## Getting the data of the straight spectrum  #######
# #############################################################################

plt.clf()

df_full = pd.read_csv('./data/spectrum_with_mode_shift.csv')

peak_indexes = find_peaks(df_full['2'].values[1:].astype('float'),
                          height=(-5.1, -3), distance=30)


plt.plot(df_full['x-axis'].values[1:].astype('float'),
         df_full['2'].values[1:].astype('float'),
         label=r'$U_{pd}$')

plt.plot(df_full['x-axis'].values[1:].astype('float'),
         df_full['1'].values[1:].astype('float'),
         label=r'$U_{piezo}$')

plt.vlines(df_full['x-axis'].values[1:].astype('float')[peak_indexes[0]][1:],
           ymin=-6,
           ymax=df_full['2'].values[1:].astype('float')[peak_indexes[0]][1:]+0.2,
           color='C3',
           linestyle='dashed',
           label='Mode hope')
plt.xlabel(r'$t \, / \, \mathrm{s}$')
plt.ylabel(r'$U \, / \, \mathrm{V}$')
plt.legend()

plt.savefig('./plots/modeshift.pdf')
