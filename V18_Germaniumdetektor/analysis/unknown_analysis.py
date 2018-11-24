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
import os

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)


# Import measurment datas
channel_content_unknown= np.genfromtxt('./2018-10-29_becker_grisard/unknown.txt',
                                       unpack=True)


peak_index_unknown = find_peaks(channel_content_unknown, height=40,
                                distance=15)
peak_indexes = peak_index_unknown
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
    #plt.hist(noms(g(np.arange(0, len(channel_content_unknown), 1), m, b)),
    #         bins=noms(g(np.linspace(0, len(channel_content_unknown),
    #         len(channel_content_unknown)), m, b)),
    #         weights=channel_content_unknown, label='Spektrum')

    #plt.plot(index_fit_plot, fit_function(index_fit_plot, *params_gaus),
    #         label='Fit')
    #plt.xlabel(r'$\mathrm{Channel}$')
    #plt.ylabel(r'$\mathrm{Count}$')
    #plt.legend()
    #plt.savefig(f'./plots/sb_or_ba/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

    # --- Return values --- #

    return amplitude, sigma, offset, int(area_under_peak)

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


efficency_calculated = exp_efficency(np.array(offset_of_peak),
                                     a, b, c)


decay_rate_calculated = decay_rate(area_under_peak, angle_distribution,
                                   prohability, efficency_calculated,
                                   measurment_time)

print('Gemitelte Aktivität', decay_rate_calculated.mean())
### --- Save data --- ###


l.Latexdocument(filename =abs_path('tabs/unknown/peak_fit_parameter.tex')).tabular(
        data=[peak_indexes[0], unp.uarray(noms(amplitude_of_peaks), stds(amplitude_of_peaks)),
          unp.uarray(noms(sigma_of_peaks), stds(sigma_of_peaks)),
          unp.uarray(noms(offset_of_peak), stds(offset_of_peak))],
    header=['Kanal / ', r'Amplitude / None ',
            r'\sigma /\kilo\eV', r'\mu / \kilo\eV'],
    places=[0, (1.2, 1.2), (1.2, 1.2), (2.2, 1.2)],
    caption='Regressionsparameter der Peak-Anpassung.',
    label='results_peaks'
)
l.Latexdocument(filename =abs_path('tabs/unknown/peak_charakteristiken.tex')).tabular(
        data=[peak_indexes[0],
          unp.uarray(noms(offset_of_peak), stds(offset_of_peak)),
           prohability,
           area_under_peak,
          unp.uarray(noms(decay_rate_calculated), stds(decay_rate_calculated))],
    header=['Kanal / ', r'\mu / \kilo\eV', r'P\ua{über} / ', r'Fläche / ',
            r'A\ua{k} / \becquerel'],
    places=[0,  (2.2, 1.2), 4, 0, (1.2, 3.2)],
    caption='Bestimmte Aktivität für jeden Peak der $^{133}\ce{Ba}$ Quelle.',
    label='decay_rate_peak'
)


peak_index_energies = g(peak_indexes[0], m, b)
E_lit = [1173.22, 1332.49]

l.Latexdocument(filename =abs_path('tabs/unknown/calculated_efficencies.tex')).tabular(
        data=[peak_indexes[0],
              unp.uarray(noms(peak_index_energies), stds(peak_index_energies)),
              E_lit,
              unp.uarray(noms(efficency_calculated),
              stds(efficency_calculated))],
    header=['Kanal / ', r'Energie / \kilo\eV',  r'E\ua{\gamma,lit}/ \kilo\eV', r'Q / '],
    places=[0, (2.2, 1.2),2, (1.2, 1.2)],
    caption='Berchente Vollenergienachweiseffizienz $^{133}\ce{Ba}$.',
    label='effizienz'
)
# ########################################################################
# ########################################################################

plt.clf()
plt.xlim(0, noms(g(4000, m, b)))
plt.yscale('log')
plt.hist(noms(g(np.arange(0, len(channel_content_unknown), 1), m, b)),
         bins=noms(g(np.linspace(0, len(channel_content_unknown),
                                 len(channel_content_unknown)), m, b)),
         weights=channel_content_unknown, label='Spektrum')

plt.plot(noms(g(peak_index_unknown[0], m, b)),
         channel_content_unknown[peak_index_unknown[0]], '.', label='Peak')
plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.ylabel(r'$\mathrm{Zählung}$')
plt.legend()
#plt.show()
plt.savefig(f'./plots/unknown/spectrum_unknown.pdf')
