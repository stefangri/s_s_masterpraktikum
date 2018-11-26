import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values
import math
import datetime
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
import os


def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)


# ########################################################################### #
# ###################### --- Bestimmung der Peaks --- ####################### #
# ########################################################################### #


channel_content_eu = np.genfromtxt('./2018-10-29_becker_grisard/europium.txt',
                                   unpack=True)

# Peaksearch for different intervals
indexes_lower = find_peaks(channel_content_eu[:1000], height=290)
indexes_upper = find_peaks(channel_content_eu[1000:], height=50, distance=10)

# The Upper Index has to be shiftet
index_upper_shiftet = indexes_upper[0]+1000

# Create on array with all peak indices
index = np.append(indexes_lower[0], index_upper_shiftet)

# ########################################################################### #
# ## --- Transformation der Skala mit Hilfe einer Linearen Regression --- ## #
# ########################################################################### #

# energies from the manual
energies = [121.78, 244.70, 344.3, 411.96, 443.96, 778.90, 867.37, 946.08,
            1085.90, 1112.10, 1408.00]

# ########################################################################### #
# ## --- Use Automatic automatic_spectrum_anaysis to get peak_index--- ## #
# ########################################################################### #


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

    # --- Calculate area under peak --- #

    area_under_peak = sum(channel_content[index_of_peak-5:index_of_peak+5])

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
    # plt.xlabel(r'$\mathrm{Channel}$')
    # plt.ylabel(r'$\mathrm{Count}$')
    # plt.legend()
    # plt.savefig(f'./plots/europium/spectrum_fit_at_index_{str(index_of_peak)}.pdf')

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

for i in index:
    amplitude, sigma, offset, area = automatic_spectrum_anaysis(channel_content_eu,
                                                                i, gaus)

    amplitude_of_peaks.append(amplitude)

    sigma_of_peaks.append(sigma)

    offset_of_peak.append(offset)

    area_under_peak.append(area)


# Save data in table
print(index)

l.Latexdocument(filename = abs_path('tabs/europium/peak_charakteristiken_eu.tex')).tabular(
        data=[index, unp.uarray(noms(amplitude_of_peaks), stds(amplitude_of_peaks)),
          unp.uarray(noms(sigma_of_peaks), stds(sigma_of_peaks)),
          unp.uarray(noms(offset_of_peak), stds(offset_of_peak))],
    header=['Kanal / ', r'A / ',
            r'\sigma / ', r'\mu / '],
    places=[0, (2.2, 2.2), (1.2, 1.2), (4.2, 1.2)],
    caption='Bestimmte Eigenschaften der Peaks von $^{152}\ce{Eu}$.',
    label='results_peaks_eu'
)
# ########################################################################### #
# ## --- Linear Regression to get transformations parameters--- ## #
# ########################################################################### #


def g(x, m, b):
    '''Define linear function for Transformation purpose'''
    return m * x + b


# Get Regressionparameters with curve_fit

params, cov = curve_fit(g, noms(offset_of_peak), energies)

print('Kovarianzmatrix Energiefit', cov)
errors = np.sqrt(np.diag(cov))

index_intervall = np.linspace(-20, 3650, 10000)



# Print FitParameters
print('\n\n------------------------------------------------------------------')
print('--------------- Fit Parameter Skalen Trafo ---------------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Steigung:', ufloat(params[0], errors[0]))
m = ufloat(params[0], errors[0])
print('y-Abschnitt:', ufloat(params[1], errors[1]))
b = ufloat(params[1], errors[1])

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# Speichere die Parameter in eine .txt damit ich sie für die spätere Auswertung
# verwenden kann.

np.savetxt('./umrechnungen_bins_to_energy.txt',
           np.column_stack([params, errors]), header='param error')

# Plotte Fit and values
plt.xlim(-10, 3600)
plt.plot(np.append(noms(offset_of_peak), 0), np.append(energies, 0), '.',
         label='Fitpunkte')
plt.plot(index_intervall, g(index_intervall, *params), label='Fit')

plt.xlabel('Kanal')
plt.ylabel(r'$\mathrm{Energie} \, / \, \mathrm{keV}$')
plt.legend()
plt.savefig('./plots/europium/skalen_trafo_fit.pdf')

# Save the  transformend channelnumbers into a table

offset_of_peak_in_energy = g(unp.uarray(noms(offset_of_peak),
                             stds(offset_of_peak)), m, b)

l.Latexdocument(filename =abs_path('tabs/europium/peak_in_energy_eu.tex')).tabular(
        data=[index, energies,
          unp.uarray(noms(offset_of_peak_in_energy), stds(offset_of_peak_in_energy))],
    header=['Kanal / ', r'E_{\gamma,theo} / \kilo\eV ',
            r'E_{\gamma} / \kilo\eV '],
    places=[0 , 2, (3.2, 1.2)],
    caption='Energiewerte der Peaks von $^{152}\ce{Eu}$.',
    label='energy__peaks_eu'
)

# ########################################################################### #
# ################# --- Bestimmung des Effizienz --- ######################## #
# ########################################################################### #

# --- Bestimmung der Winkelverteilung --- #

distance_probe_detektor = 73.1e-3 + 1.5e-2  # Abstand Probe Alu + Dicke Alu
radiant_of_detector = 45e-3 / 2


def winkelverteilung(distance_probe_detektor, radiant_of_detector):
    ''' Function to calculate angle distribution OMEGA '''

    return 1/2 * (1 - (distance_probe_detektor) /
                  (np.sqrt(distance_probe_detektor**2
                   + radiant_of_detector**2)))


angle_distribution = winkelverteilung(distance_probe_detektor,
                                      radiant_of_detector)

print('\n\n----------------------------------------------------------------------')
print('--------------- Winkelverteilung ---------------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Winkelverteilung:', angle_distribution)

print('\n----------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('----------------------------------------------------------------------\n\n')

# --- Save angle_distribution to use it in later projects

np.savetxt('./winkelverteilung.txt',
           np.column_stack([angle_distribution]), header='Winkelverteilung')

# --- Bestimme die heutige Aktivität --- #

halbwertszeit = ufloat(4943, 5)  # days
A_0 = ufloat(4130, 60)  # Bq


def decay_rate(test_time, half_life, start_decay_rate):
    ''' Function to calculate the decay_rate of europium at any given date '''
    ''' half_life has to be in units of days '''

    test_time = test_time.split('.')
    test_date = datetime.date(int(test_time[0]), int(test_time[1]), int(test_time[2]))  # test_time should be someting like: '2018.10.29'
    # Source: user_manual
    creation_date = datetime.date(2000, 10, 1)
    # How many days past since the creation
    past_time = (test_date-creation_date).days
    # Dif has to be bigger then zero
    assert past_time >= 0
    # Calculate decay_rate of test_date
    return start_decay_rate * unp.exp(- (past_time *np.log(2) ) / half_life)


decay_rate_29_10_18 = decay_rate('2018.10.29', halbwertszeit, A_0)

print('\n\n----------------------------------------------------------------------')
print('--------------- Aktivität am Messtag ---------------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Aktivität am 29.10.18:', decay_rate_29_10_18)

print('\n----------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('----------------------------------------------------------------------\n\n')


# --- Definiere Funktion um die Effizienz zu bestimmen --- #


def full_energy_efficiency(measured_decay_rate, angle_distribution,
                           soll_decay_rate, transition_probaility,
                           measurment_time):
    ''' Function to calculate the full_energy_efficiency '''
    return measured_decay_rate / (angle_distribution * soll_decay_rate
                                  * transition_probaility * measurment_time)


# --- Bestimmte die gemessene Aktivität um jeden Peak herum --- ###
'''  Um die Aktivität um jeden Peak herum zu bestimmen, gehe ich vom Peaks
    4 Bins nach rechts und links summiere dann über den Bininhalt  '''

transition_probaility = [0.286, 0.076, 0.265, 0.022, 0.031, 0.129, 0.042,
                         0.146, 0.102, 0.136, 0.210]

efficiency = []
summierte_akti = []
measurment_time = 3380
for i in range(len(index)):
    summierte_aktivität = sum(channel_content_eu[index[i]-4:index[i]+4])
    print(i, summierte_aktivität)

    summierte_akti.append(summierte_aktivität)
    print(type(angle_distribution), angle_distribution)
    efficiency.append(full_energy_efficiency(summierte_aktivität,
                                             angle_distribution,
                                             decay_rate_29_10_18,
                                             transition_probaility[i],
                                             measurment_time))

print('\n\n----------------------------------------------------------------------')
print('--------------- Bestimmten Effizienz ---------------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Energie-Effizienz-Wertepaar:')

for i in range(len(energies)):
    print('E/keV:', energies[i],'Q:', efficiency[i])

print('\n----------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('----------------------------------------------------------------------\n\n')


# --- Fitte an die Effizienz eine exp() --- #


def exp(x, a, b, c):
    return a * np.exp(-b * (x))+c


def efficeny_function(E, a, b):
    return a * E**b


'Beachte das ich hier den ersten Wert wie in der Anleitung gefordert wegwerfe!'
params_exp, cov_exp = curve_fit(exp, energies[1:], noms(efficiency[1:]),
                                p0=[2, 0.01, 0])
errors_exp = np.sqrt(np.diag(cov_exp))

params_eff, cov_eff = curve_fit(f=efficeny_function, xdata=energies[1:],
                                ydata=noms(efficiency[1:]))

errors_efficency = np.sqrt(np.diag(cov_eff))

print('\n\n------------------------------------------------------------------')
print('----------- Fit Parameter Effizienzbestimmung -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('a:', ufloat(params_exp[0], errors_exp[0]))
a = ufloat(params_exp[0], errors_exp[0])
print('b:', ufloat(params_exp[1], errors_exp[1]))
b = ufloat(params_exp[1], errors_exp[1])
print('c:', ufloat(params_exp[2], errors_exp[2]))
c = ufloat(params_exp[2], errors_exp[2])

print('--------------------- Effizienzfunktion -------------------------------')

a_efficency = ufloat(params_eff[0], errors_efficency[0])
b_efficency = ufloat(params_eff[1], errors_efficency[1])

print('a:', a_efficency, a_efficency/62)
print('b:', b_efficency)

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# Save the exp parameters in file to reuse in a later project

np.savetxt('./params_efficency.txt',
           np.column_stack([params_exp, errors_exp]), header='param error')

#  Generating some values to plot the fit

e = np.linspace(energies[1]-100, energies[-1]+100, 1000)

plt.clf()
plt.errorbar(energies, noms(efficiency), yerr=stds(efficiency), fmt='.',
             label='Datenpunkte')
plt.plot(e, exp(e, *params_exp), label=r'$Q_2(E)$')
plt.plot(e, efficeny_function(e, *params_eff), label=r'$Q_1(E)$')
plt.xlim(energies[1]-50, energies[-1]+50)
plt.xlabel(r'$\mathrm{Energie} \, / \, \mathrm{keV}$')
plt.ylabel(r'$\mathrm{Effizienz}$')
plt.legend()

plt.savefig('./plots/europium/effizienz.pdf')

# ########################################################################### #
# ############ --- Speicherergebnisse in eine Tabelle --- ################### #
# ########################################################################### #

#print(summierte_akti, transition_probaility)
#l.Latexdocument(filename =abs_path('tabs/europium/results_europium.tex')).tabular(
#    data = [index, energies, summierte_akti, transition_probaility, unp.uarray(noms(efficiency), stds(efficiency))],
#    header = ['Kanal / ', r'Energie \, / \kilo\eV', r'Z / \becquerel', r'W /',
#              r'Effizienz /'],
#    places = [0, 2, 4, 1.3, (1.4, 1.4)],
#    caption = 'Bestimmten Energie und Effizienzwerte.',
#    label = 'results_europium'
#)

l.Latexdocument(filename =abs_path('tabs/europium/results_europium.tex')).tabular(
    data = [index, energies, summierte_akti,  transition_probaility,
            unp.uarray(noms(efficiency), stds(efficiency))],
    header = ['Kanal / ', r'Energie \, / \kilo\eV', r'Z / \becquerel', r'W /',
              r'Effizienz /'],
    places = [0, 2, 4, 3,  (1.4, 1.4)],
    caption = 'Bestimmten Energie und Effizienzwerte.',
    label = 'results_europium'
)

# ########################################################################### #
# ############ --- Plotte das Spektrum mit Peak-Detection --- ############### #
# ########################################################################### #

#Plot with Binnumbers as x axis
plt.clf()
plt.hist(range(0, len(channel_content_eu), 1),
         bins=np.linspace(0, len(channel_content_eu), len(channel_content_eu)),
         weights=channel_content_eu, label='Spektrum')

plt.plot(indexes_lower[0], channel_content_eu[indexes_lower[0]], '.',
         markersize=2, label='Peaks', color='C1', alpha=0.8)
plt.plot(index_upper_shiftet, channel_content_eu[index_upper_shiftet], '.',
         markersize=2, color='C1', alpha=0.8)
# plt.ylim(0, )
plt.xlim(0, 4000)

plt.ylabel(r'$\mathrm{Zählungen}$')
plt.xlabel(r'$\mathrm{Kanal}$')
plt.legend()
plt.yscale('log')
plt.savefig('./plots/europium/spektrum_index.pdf')

# Plot with Energies as x axis

plt.clf()

index_to_energie = g(range(0, len(channel_content_eu), 1), *params)

plt.hist(np.linspace(0., max(index_to_energie), len(channel_content_eu)),
         bins=np.linspace(0., max(index_to_energie), len(channel_content_eu)),
         weights=channel_content_eu, label='Spektrum')

plt.xlim(0, max(index_to_energie))
plt.ylabel(r'$\mathrm{Zählungen}$')
plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.legend()
plt.savefig('./plots/europium/spektrum_energie.pdf')
# plt.show()
