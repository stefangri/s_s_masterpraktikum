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

# ---  Import params for efficiency calculations --- #


params_exp, errors_exp = np.genfromtxt('params_efficency.txt', unpack=True)

a = ufloat(params_exp[0], errors_exp[0])
b = ufloat(params_exp[1], errors_exp[1])
c = ufloat(params_exp[2], errors_exp[2])


def exp_efficency(energy, a, b, c):
    return a * unp.exp(-b * (energy))+c


def decay_rate(area, angle_distribution, prohability, efficiency, time):
    print(area, prohability, efficiency)
    return area / (angle_distribution * prohability * efficiency * time)


# --- Find index of peak
peak_indexes = find_peaks(channel_content_sb_ba, height=120, distance=15)

# Wirte a fuction for automatic peak analysis


def gaus(x, amplitude, sigma, offset):
    return amplitude * np.exp(-1/2 * (x - offset)**2/sigma**2)


def automatic_spectrum_anaysis(channel_content, index_of_peak, fit_function,):

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

    # --- Calculate area under peak --- #

    area_under_peak = sum(channel_content[index_of_peak-5:index_of_peak+5])

    # --- Plot function ---  #

    index_fit_plot = np.linspace(index_of_peak-15, index_of_peak+15, 1e4)

    #plt.clf()
    #plt.xlim(index_of_peak-20, index_of_peak+20)
    #plt.ylim(0, channel_content[index_of_peak] * 1.2)
    #plt.hist(range(0, len(channel_content_sb_ba), 1),
    #         bins=np.linspace(0, len(channel_content_sb_ba),
    #         len(channel_content_sb_ba)),
    #         weights=channel_content_sb_ba, label='Spektrum')

    #plt.plot(index_fit_plot, fit_function(index_fit_plot, *params_gaus),
    #         label='Fit')
    #plt.xlabel(r'$\mathrm{Binnummer}$')
    #plt.ylabel(r'$\mathrm{Counts}$')
    #plt.legend()
    #plt.savefig(f'./plots/sb_or_ba/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

    # --- Return values --- #

    return amplitude, sigma, offset, area_under_peak

# ########################################################################### #
# ########################## Using the function ############################# #
# ########################################################################### #


amplitude_of_peaks = []

sigma_of_peaks = []
sigma_of_peaks_energy = []

offset_of_peak = []
offset_of_peak_in_energy = []

area_under_peak = []

for index in peak_indexes[0]:
    amplitude, sigma, offset, area = automatic_spectrum_anaysis(channel_content_sb_ba,
                                                          index, gaus)

    amplitude_of_peaks.append(amplitude)

    sigma_of_peaks.append(sigma)
    sigma_of_peaks_energy.append(g(sigma, m, b))

    offset_of_peak.append(offset)
    offset_of_peak_in_energy.append(g(offset, m, b))

    area_under_peak.append(area)


print(offset_of_peak_in_energy)
print(area)


# --- Calculate the normal_acitvity for all --- #

measurment_time = 2941  # in seconds

# --- Get the angle_distribution --- #

angle_distribution = np.genfromtxt('winkelverteilung.txt', unpack=True)

# source: https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html
prohability = [0.33, 0.07, 0.18, 0.62, 0.09]


efficency_calculatet = exp_efficency(np.array(offset_of_peak_in_energy),
                                     a, b, c)

decay_rate_calculated = decay_rate(area_under_peak, angle_distribution,
                                   prohability, efficency_calculatet,
                                   measurment_time)


#print(decay_rate_calculated, decay_rate_calculated.mean())

# ########################################################################### #
# #### --- Speicherergebnisse Peak Eigenschaften in eine Tabelle --- ######## #
# ########################################################################### #



#print(len(unp.uarray(noms(sigma_of_peaks_energy), stds(sigma_of_peaks_energy))), type(unp.uarray(noms(sigma_of_peaks_energy), stds(sigma_of_peaks_energy))))
#print(len(unp.uarray(noms(offset_of_peak_in_energy), stds(offset_of_peak_in_energy))), type(unp.uarray(noms(offset_of_peak_in_energy), stds(offset_of_peak_in_energy))))
#print(len(area_under_peak), type(area_under_peak))
#print(len(unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))), type(unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))))


l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/sb_or_ba/peak_charakteristiken.tex').tabular(
        data=[peak_indexes[0], prohability, unp.uarray(noms(amplitude_of_peaks), stds(amplitude_of_peaks)),
          unp.uarray(noms(sigma_of_peaks_energy), stds(sigma_of_peaks_energy)),
          unp.uarray(noms(offset_of_peak_in_energy), stds(offset_of_peak_in_energy)),
          area_under_peak,
          unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))],
    header=['Binnummer / ', r'Übergangswahrscheinlichkeit / ', r'Amplitude / ',
            r'\sigma / \kilo\eV', r'\mu / \kilo\eV', r'Fläche / ',
            r'Aktivität / \Bq'],
    places=[0, 2, (1.4, 1.4), (1.4, 1.4), (1.4, 1.4), 2, (1.2, 1.2)],
    caption='Bestimmte Eigenschaften der Peaks.',
    label='results_peaks'
)


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
