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
channel_content_sb_ba = np.genfromtxt('./2018-10-29_becker_grisard/sb_or_ba.txt',
                                   unpack=True)

# --- Import params to transform binnumber to energy ---#
params, errors = np.genfromtxt('./umrechnungen_bins_to_energy.txt',
                               unpack=True)
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])


def g(x, m, b):
    '''Define linear function for Transformation purpose'''
    return m * x + b


# --- Find index of peak
peak_indexes = find_peaks(channel_content_sb_ba, height=120, distance=15)

# Wirte a fuction for automatic peak analysis


def gaus(x, amplitude, sigma, offset):
    return amplitude * np.exp(-1/2 * (x - offset)**2/sigma**2)


def automatic_spectrum_anaysis(channel_content, index_of_peak, fit_function):

    # use as fit_intervall \pm 15, create for fit numpy array with index values
    index = np.arange(index_of_peak-15, index_of_peak+15, 1)

    params_gaus, cov_gaus = curve_fit(fit_function, index,
                                      channel_content[index],
                                      p0=[channel_content[index_of_peak], 1,
                                          index_of_peak])

    error_gaus = np.sqrt(np.diag(cov_gaus))

    amplitude = ufloat(params_gaus[0], error_gaus[0])
    sigma = ufloat(params_gaus[1], error_gaus[1])
    offset = ufloat(params_gaus[2], error_gaus[2])

    # --- Plot function ---  #

    index_fit_plot = np.linspace(index_of_peak-15, index_of_peak+15, 1e4)

    plt.clf()
    plt.xlim(index_of_peak-20, index_of_peak+20)
    plt.ylim(0, channel_content[index_of_peak] * 1.2)
    plt.hist(range(0, len(channel_content_sb_ba), 1),
             bins=np.linspace(0, len(channel_content_sb_ba),
             len(channel_content_sb_ba)),
             weights=channel_content_sb_ba, label='Spektrum')

    plt.plot(index_fit_plot, fit_function(index_fit_plot, *params_gaus),
             label='Fit')
    plt.xlabel(r'$\mathrm{Binnummer}$')
    plt.ylabel(r'$\mathrm{Counts}$')
    plt.legend()
    plt.savefig(f'./plots/sb_or_ba/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

    # --- Return values --- #

    return amplitude, sigma, offset

# ########################################################################### #
# ########################## Using the function ############################# #
# ########################################################################### #


amplitude_of_peaks = []

sigma_of_peaks = []
sigma_of_peaks_energy = []

offset_of_peak = []
offset_of_peak_in_energy = []

print(type(m), type(b))
for index in peak_indexes[0]:
    amplitude, sigma, offset = automatic_spectrum_anaysis(channel_content_sb_ba,
                                                          index, gaus)

    amplitude_of_peaks.append(amplitude)

    sigma_of_peaks.append(sigma)
    sigma_of_peaks_energy.append(g(sigma, m, b))

    offset_of_peak.append(offset)
    offset_of_peak_in_energy.append(g(offset, m, b))

print(offset_of_peak)
print('\n')
print(offset_of_peak_in_energy)
# ########################################################################### #
# ########################################################################### #
# ########################################################################### #


plt.clf()
plt.xlim(0, 1200)
plt.hist(range(0, len(channel_content_sb_ba), 1),
         bins=np.linspace(0, len(channel_content_sb_ba),
         len(channel_content_sb_ba)),
         weights=channel_content_sb_ba, label='Spektrum')

plt.plot(peak_indexes[0], channel_content_sb_ba[peak_indexes[0]], '.')
plt.xlabel(r'$\mathrm{Binnummer}$')
plt.ylabel(r'$\mathrm{Counts}$')
plt.legend()
plt.savefig(f'./plots/sb_or_ba/spectrum.pdf')

# ### ToDo:
# -calculating the sum
# -
#
#
