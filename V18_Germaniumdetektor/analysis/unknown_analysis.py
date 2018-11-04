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
channel_content_unknown= np.genfromtxt('./2018-10-29_becker_grisard/unknown.txt',
                                       unpack=True)


peak_index_unknown = find_peaks(channel_content_unknown, height=40,
                                distance=15)
# --- Import params to transform binnumber to energy ---#
params, errors = np.genfromtxt('./umrechnungen_bins_to_energy.txt',
                               unpack=True)
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])


def g(x, m, b):
    '''Define linear function for Transformation purpose'''
    return m * x + b

# Transform index to energy


peak_energy_unknown = g(peak_index_unknown[0], m, b)
# ---  Import params for efficiency calculations --- #


params_exp, errors_exp = np.genfromtxt('params_efficency.txt', unpack=True)

a = ufloat(params_exp[0], errors_exp[0])
b = ufloat(params_exp[1], errors_exp[1])
c = ufloat(params_exp[2], errors_exp[2])


def exp_efficency(energy, a, b, c):
    return a * unp.exp(-b * energy)+c


def decay_rate(area, angle_distribution, prohability, efficiency, time):
    return area/(angle_distribution * prohability * efficiency * time)


# --- Automatic Spectrum anaysis --- #

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

    area_under_peak = sum(channel_content[index_of_peak-10:index_of_peak+10])

    # --- Plot function ---  #

    index_fit_plot = np.linspace(index_of_peak-15, index_of_peak+15, 1e4)

    # plt.clf()
    # plt.xlim(index_of_peak-20, index_of_peak+20)
    # plt.ylim(0, channel_content[index_of_peak] * 1.2)
    # plt.hist(range(0, len(channel_content), 1),
             # bins=np.linspace(0, len(channel_content),
             # len(channel_content)),
             # weights=channel_content, label='Spektrum')

    # plt.plot(index_fit_plot, fit_function(index_fit_plot, *params_gaus),
             # label='Fit')
    # plt.xlabel(r'$\mathrm{Binnummer}$')
    # plt.ylabel(r'$\mathrm{Counts}$')
    # plt.legend()
    # plt.savefig(f'./plots/unknown/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

    # --- Return values --- #

    return amplitude, sigma, offset, area_under_peak

# ########################################################################
# ################# Make all the necessary anaylisies ###################


amplitude_of_peaks = []

sigma_of_peaks = []
sigma_of_peaks_energy = []

offset_of_peak = []
offset_of_peak_in_energy = []

area_under_peak = []

for index in peak_index_unknown[0]:
    amplitude, sigma, offset, area = automatic_spectrum_anaysis(channel_content_unknown,
                                                          index, gaus)

    amplitude_of_peaks.append(amplitude)

    sigma_of_peaks.append(sigma)
    sigma_of_peaks_energy.append(g(sigma, m, b))

    offset_of_peak.append(offset)
    offset_of_peak_in_energy.append(g(offset, m, b))

    area_under_peak.append(area)


print(offset_of_peak_in_energy)
# --- Calculate the normal_acitvity for all --- #

measurment_time = 3412  # in seconds

# --- Get the angle_distribution --- #

angle_distribution = np.genfromtxt('winkelverteilung.txt', unpack=True)

# source: https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html
prohability = [0.9958, 0.9998]


efficency_calculated = exp_efficency(np.array(offset_of_peak_in_energy),
                                     a, b, c)


decay_rate_calculated = decay_rate(area_under_peak, angle_distribution,
                                   prohability, efficency_calculated,
                                   measurment_time)

### --- Save data --- ###

l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/unknown/peak_charakteristiken_unknown.tex').tabular(
        data=[peak_index_unknown[0], prohability, unp.uarray(noms(amplitude_of_peaks), stds(amplitude_of_peaks)),
          unp.uarray(noms(sigma_of_peaks_energy), stds(sigma_of_peaks_energy)),
          unp.uarray(noms(offset_of_peak_in_energy), stds(offset_of_peak_in_energy)),
          area_under_peak,
          unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))],
    header=['Binnummer / ', r'Übergangswahrscheinlichkeit / ', r'Amplitude / ',
            r'\sigma / \kilo\eV', r'\mu / \kilo\eV', r'Fläche / ',
            r'Aktivität / \Bq'],
    places=[0, 2, (1.4, 1.4), (1.4, 1.4), (1.4, 1.4), 2, (1.2, 1.2)],
    caption='Bestimmte Eigenschaften der Peaks von $^{60}\ce{Co}$.',
    label='results_peaks_unknown'
)


l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/unknown/calculated_efficencies.tex').tabular(
        data=[peak_index_unknown[0],
              unp.uarray(noms(peak_energy_unknown), stds(peak_energy_unknown)),
              unp.uarray(noms(efficency_calculated),
              stds(efficency_calculated))],
    header=['Binnummer / ', r'Energie / \kilo\eV', r'Effizienz / '],
    places=[0, (1.2, 1.2), (1.2, 1.2)],
    caption='Berchente Vollenergienachweiseffizienz von $^{60}\ce{Co}$.',
    label='effizienz'
)
# ########################################################################
# ########################################################################
# plt.clf()
# plt.xlim(0, 4000)
# plt.hist(range(0, len(channel_content_unknown), 1),
         # bins=np.linspace(0, len(channel_content_unknown),
         # len(channel_content_unknown)),
         # weights=channel_content_unknown, label='Spektrum')
#
# plt.plot(peak_index_unknown[0],
         # channel_content_unknown[peak_index_unknown[0]], '.')
# plt.xlabel(r'$\mathrm{Binnummer}$')
# plt.ylabel(r'$\mathrm{Counts}$')
# plt.legend()
# plt.show()
# plt.savefig(f'./plots/unknown/spectrum.pdf')
