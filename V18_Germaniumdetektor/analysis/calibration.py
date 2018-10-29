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


def g(x, m, b):
    '''Define linear function for Transformation purpose'''
    return m * x + b


# Get Regressionparameters with curve_fit

params, cov = curve_fit(g, index, energies)
errors = np.sqrt(np.diag(cov))

index_intervall = np.linspace(min(index)-10, max(index)+10, 10000)

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
plt.plot(index, energies, '.', label='Datenpunkte')
plt.plot(index_intervall, g(index_intervall, *params), label='Fit')

plt.xlabel('Binnummer')
plt.ylabel(r'$\mathrm{Energie} \, / \, \mathrm{keV}$')
plt.legend()
plt.savefig('./plots/europium/skalen_trafo_fit.pdf')


# ########################################################################### #
# ################# --- Bestimmung des Effizienz --- ######################## #
# ########################################################################### #

# --- Bestimmung der Winkelverteilung --- #

distance_probe_detektor = 73.1e-3 + 1.5e-2  # Abstand Probe Alu + Dicke Alu
radiant_of_detector = 45e-3 / 2


def winkelverteilung(distance_probe_detektor, radiant_of_detector):
    ''' Function to calculate angle distribution OMEGA '''

    return 1/2 * (1 - (distance_probe_detektor /
                  (np.sqrt(distance_probe_detektor**2
                   + radiant_of_detector**2))))


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
    return start_decay_rate * unp.exp(- past_time / half_life)


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
measurment_time = 3380
for i in range(len(index)):
    summierte_aktivität = sum(channel_content_eu[index[i]-4:index[i]+4])
    print(i, summierte_aktivität)

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
    return a * np.exp(-b * x) + c


'Beachte das ich hier den ersten Wert wie in der Anleitung gefordert wegwerfe!'
params_exp, cov_exp = curve_fit(exp, energies[1:], noms(efficiency[1:]),
                                p0=[0.5, 0.01, 0.1])
errors_exp = np.sqrt(np.diag(cov_exp))

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

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

#  Generating some values to plot the fit

e = np.linspace(energies[1]-100, energies[-1]+100, 1000)

plt.clf()
plt.errorbar(energies, noms(efficiency), yerr=stds(efficiency), fmt='.',
             label='Datenpunkte')
plt.plot(e, exp(e, *params_exp), label='Fit')

plt.xlim(energies[1]-50, energies[-1]+50)
plt.xlabel(r'$\mathrm{Energie} \, / \, \mathrm{keV}$')
plt.ylabel(r'$\mathrm{Effizienz}$')
plt.legend()

plt.savefig('./plots/europium/effizienz.pdf')

# ########################################################################### #
# ############ --- Speicherergebnisse in eine Tabelle --- ################### #
# ########################################################################### #


l.Latexdocument(filename ='/home/beckstev/Documents/s_s_masterpraktikum/V18_Germaniumdetektor/analysis/tabs/europium/results_europium.tex').tabular(
    data = [index, energies, unp.uarray(noms(efficiency), stds(efficiency))],
    header = ['Binnummer / ', r'Energie \, / \kilo\eV', r'Effizienz /'],
    places = [0, 2, (1.4, 1.4)],
    caption = 'Bestimmten Energie und Effizienzwerte.',
    label = 'results_europium'
)


# ########################################################################### #
# ############ --- Plotte das Spektrum mit Peak-Detection --- ############### #
# ########################################################################### #

# Plot with Binnumbers as x axis
# plt.clf()
# plt.hist(range(0, len(channel_content_eu), 1),
#          bins=np.linspace(0, len(channel_content_eu), len(channel_content_eu)),
#          weights=channel_content_eu, label='Spektrum')
#
# plt.plot(indexes_lower[0], channel_content_eu[indexes_lower[0]], 'x',
#          markersize=1, label='Peaks', color='C1', alpha=0.6)
# plt.plot(index_upper_shiftet, channel_content_eu[index_upper_shiftet], 'x',
#          markersize=1, color='C1', alpha=0.6)
# # plt.ylim(0, )
# plt.xlim(0, 4000)
#
# plt.ylabel(r'$\mathrm{Counts}$')
# plt.xlabel(r'$\mathrm{Binnummer}$')
# plt.legend()
# plt.savefig('./plots/europium/spektrum_index.pdf')

# Plot with Energies as x axis

# plt.clf()
#
# index_to_energie = g(range(0, len(channel_content_eu), 1), *params)
#
# plt.hist(np.linspace(0., max(index_to_energie), len(channel_content_eu)),
#          bins=np.linspace(0., max(index_to_energie), len(channel_content_eu)),
#          weights=channel_content_eu, label='Spektrum')
#
# plt.xlim(0, max(index_to_energie))
# plt.ylabel(r'$\mathrm{Counts}$')
# plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
# plt.legend()
# plt.savefig('./plots/europium/spektrum_energie.pdf')
# ült.show()
