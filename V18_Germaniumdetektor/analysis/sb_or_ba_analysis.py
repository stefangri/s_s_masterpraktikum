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

peak_index_energies = noms(g(peak_indexes[0], m, b))
# Wirte a fuction for automatic peak analysis


def gaus(x, amplitude, sigma, offset):
    return amplitude * np.exp(-1/2 * (x - offset)**2/sigma**2)


def automatic_spectrum_anaysis(channel_content, index_of_peak,
                               fit_function,):

    # use as fit_intervall \pm 15, create for fit numpy array with index values
    index_channel = np.arange(index_of_peak-15, index_of_peak+15, 1)
    index = noms(g(index_channel, m, b))

    params_gaus, cov_gaus = curve_fit(fit_function, index,
                                      channel_content[index_channel],
                                      p0=[channel_content[index_of_peak], 1,
                                          noms(g(index_of_peak, m, b))])

    error_gaus = np.sqrt(np.diag(cov_gaus))

    amplitude = ufloat(params_gaus[0], error_gaus[0])
    sigma = ufloat(params_gaus[1], error_gaus[1])
    offset = ufloat(params_gaus[2], error_gaus[2])

    # --- Calculate area under peak --- #

    area_under_peak = sum(channel_content[index_of_peak-5:index_of_peak+5])

    # --- Plot function ---  #

    index_fit_plot = np.linspace(index_of_peak-15, index_of_peak+15, 1e4)

    #plt.clf()
    #plt.xlim(noms(g(index_of_peak, m, b))-1, noms(g(index_of_peak, m, b))+1)
    #plt.ylim(0, channel_content[index_of_peak] * 1.2)
    #plt.hist(noms(g(np.arange(0, len(channel_content_sb_ba), 1), m, b)),
    #         bins=noms(g(np.linspace(0, len(channel_content_sb_ba),
    #         len(channel_content_sb_ba)), m, b)),
    #         weights=channel_content_sb_ba, label='Spektrum')

    #plt.plot(index_fit_plot, fit_function(index_fit_plot, *params_gaus),
    #         label='Fit')
    #plt.xlabel(r'$\mathrm{Channel}$')
    #plt.ylabel(r'$\mathrm{Count}$')
    #plt.legend()
    #plt.savefig(f'./plots/sb_or_ba/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

    # --- Return values --- #

    return amplitude, sigma, offset, int(area_under_peak)

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


print(amplitude_of_peaks)
print(sigma_of_peaks)
print(offset_of_peak)

# --- Calculate the normal_acitvity for all --- #

measurment_time = 2941  # in seconds

# --- Get the angle_distribution --- #

angle_distribution = np.genfromtxt('winkelverteilung.txt', unpack=True)

# source: https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html
prohability = [0.33, 0.07, 0.18, 0.62, 0.09]


efficency_calculated = exp_efficency(np.array(offset_of_peak_in_energy),
                                     a, b, c)

decay_rate_calculated = decay_rate(area_under_peak, angle_distribution,
                                   prohability, efficency_calculated,
                                   measurment_time)


print('Gemittelte Aktivität', decay_rate_calculated.mean())

# ########################################################################### #
# #### --- Speicherergebnisse Peak Eigenschaften in eine Tabelle --- ######## #
# ########################################################################### #

l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/sb_or_ba/peak_fit_parameter.tex').tabular(
        data=[peak_indexes[0], unp.uarray(noms(amplitude_of_peaks), stds(amplitude_of_peaks)),
          unp.uarray(noms(sigma_of_peaks), stds(sigma_of_peaks)),
          unp.uarray(noms(offset_of_peak), stds(offset_of_peak))],
    header=['Kanal / ', r'Amplitude / None ',
            r'\sigma /\kilo\eV', r'\mu / \kilo\eV'],
    places=[0, (1.2, 1.2), (1.2, 1.2), (2.2, 1.2)],
    caption='Regressionsparameter der Peak-Anpassung.',
    label='results_peaks'
)

l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/sb_or_ba/peak_charakteristiken.tex').tabular(
        data=[peak_indexes[0],
          unp.uarray(noms(offset_of_peak), stds(offset_of_peak)),
           prohability,
           area_under_peak,
          unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))],
    header=['Kanal / ', r'\mu / \kilo\eV', r'P\ua{über} / ', r'Fläche / ',
            r'A\ua{k} / \becquerel'],
    places=[0,  (2.2, 1.2), 2, 0, (1.2, 3.2)],
    caption='Bestimmte Aktivität für jeden Peak der $^{133}\ce{Ba}$ Quelle.',
    label='decay_rate_peak'
)


l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/sb_or_ba/calculated_efficencies.tex').tabular(
        data=[peak_indexes[0],
              unp.uarray(noms(peak_index_energies), stds(peak_index_energies)),
              unp.uarray(noms(efficency_calculated),
              stds(efficency_calculated))],
    header=['Kanal / ', r'Energie / \kilo\eV', r'Q / '],
    places=[0, (2.2, 1.2), (1.2, 1.2)],
    caption='Berchente Vollenergienachweiseffizienz $^{133}\ce{Ba}$.',
    label='effizienz'
)


# ########################################################################### #
# ########################################################################### #
# ########################################################################### #


plt.clf()
plt.xlim(0, noms(g(1200, m, b)))
plt.yscale('log')
plt.hist(noms(g(np.arange(0, len(channel_content_sb_ba), 1), m, b)),
         bins=noms(g(np.linspace(0, len(channel_content_sb_ba),
                                 len(channel_content_sb_ba)), m, b)),
         weights=channel_content_sb_ba, label='Spektrum')

plt.plot(noms(g(peak_indexes[0],m ,b)), channel_content_sb_ba[peak_indexes[0]],
         '.', label='Peak')

plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.ylabel(r'$\mathrm{Zählung}$')
plt.legend()
plt.savefig(f'./plots/sb_or_ba/spectrum.pdf')

# ### ToDo:
# -calculating the sum
# -
#
#
