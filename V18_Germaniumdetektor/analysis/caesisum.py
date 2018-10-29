import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pint import UnitRegistry
import latex as l
r = l.Latexdocument('plots.tex')
u = UnitRegistry()
Q_ = u.Quantity
import pandas as pd
from pandas import Series, DataFrame
from scipy.signal import find_peaks

# Import measurment datas
channel_content_ca = np.genfromtxt('./2018-10-29_becker_grisard/caesium.txt',
                                   unpack=True)

# --- Import params to transform binnumber to energy ---#
params, errors = np.genfromtxt('./umrechnungen_bins_to_energy.txt',
                               unpack=True)
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])


def g(x, m, b):
    return m * x + b


# --- Find index of peak
peak_indexes = find_peaks(channel_content_ca, height=200, distance=10)

# ########################################################################### #
# ###################### --- Untersuche den Photopeak --- ################### #
# ########################################################################### #

# --- Definiere Gauß-Funktion um Photopeak zu fitten --- #


def gaus(x, amplitude, sigma, offset):
    return amplitude * np.exp(-1/2 * (x - offset)**2/sigma**2)


index = np.arange(peak_indexes[0]-15, peak_indexes[0]+15)

params_gaus, cov_gaus = curve_fit(gaus, index, channel_content_ca[peak_indexes[0][0]-15:peak_indexes[0][0]+15], p0=[channel_content_ca[peak_indexes[0]][0], 10, 1648])
errors_gaus = np.sqrt(np.diag(cov_gaus))

print('\n\n------------------------------------------------------------------')
print('----------- Fit Parameter Photpeakvermessung -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Amplitude/E_gamma:', ufloat(params_gaus[0], errors_gaus[0]))
amplitude = ufloat(params_gaus[0], errors_gaus[0])

print('Halbertsbreite:', ufloat(params_gaus[1], errors_gaus[1]))
sigma = ufloat(params_gaus[1], errors_gaus[1])

print('Offset/bins:', ufloat(params_gaus[2], errors_gaus[2]))
offset = ufloat(params_gaus[2], errors_gaus[2])
offset_energy = g(offset, m, b)
print('Offset/energy [keV]:', offset_energy)
print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# --- Bestimme den Flächeninhalt des Photopeaks --- #

area_photo_peak = sum(channel_content_ca[peak_indexes[0][0]-15:
                                         peak_indexes[0][0]+15])

print('\n\n------------------------------------------------------------------')
print('----------- Fläche unter dem Photopeak -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Flächephotopeak:', area_photo_peak)
print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')


index = np.linspace(peak_indexes[0]-20, peak_indexes[0]+20,1e3)

plt.xlim( peak_indexes[0]-30, peak_indexes[0]+30)
plt.hist(range(0, len(channel_content_ca), 1),
          bins=np.linspace(0, len(channel_content_ca), len(channel_content_ca)),
          weights=channel_content_ca, label='Spektrum')
plt.plot(index, gaus(index, *params_gaus), label='Fit')

plt.xlabel(r'$\mathrm{Binnummer}$')
plt.ylabel(r'$\mathrm{Counts}$')
plt.legend()
#plt.show()
plt.savefig('./plots/caesium/photopeak.pdf')
