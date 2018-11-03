import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pint import UnitRegistry
import scipy.constants as const
import latex as l
r = l.Latexdocument('plots.tex')
u = UnitRegistry()
Q_ = u.Quantity
import pandas as pd
from pandas import Series, DataFrame
from scipy.signal import find_peaks

# Import measurment datas
channel_content_unkown = np.genfromtxt('./2018-10-29_becker_grisard/unknwon.txt',
                                   unpack=True)


# --- Import params to transform binnumber to energy ---#
params, errors = np.genfromtxt('./umrechnungen_bins_to_energy.txt',
                               unpack=True)
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])


def g(x, m, b):
    '''Define linear function for Transformation purpose'''
    return m * x + b

# ---  Import params for efficiency calculations --- #


params_exp, errors_exp = np.genfromtxt('params_efficency.txt', unpack=True)

a = ufloat(params_exp[0], errors_exp[0])
b = ufloat(params_exp[1], errors_exp[1])
c = ufloat(params_exp[2], errors_exp[2])


def exp_efficency(energy, a, b, c):
    return a * unp.exp(-b * (energy)+c


def decay_rate(area, angle_distribution, prohability, efficiency, time):
    print(area, prohability, efficiency)
    return area/(angle_distribution * prohability * efficiency * time)



plt.clf()
plt.xlim(0, 1200)
plt.hist(range(0, len(channel_content_unknown), 1),
         bins=np.linspace(0, len(channel_content_unknown),
         len(channel_content_unknown)),
         weights=channel_content_unknown, label='Spektrum')

plt.plot(peak_indexes[0], channel_content_unknown[peak_indexes[0]], '.')
plt.xlabel(r'$\mathrm{Binnummer}$')
plt.ylabel(r'$\mathrm{Counts}$')
plt.legend()
plt.show()
#plt.savefig(f'./plots/unkwon/spectrum.pdf')
