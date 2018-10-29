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
r = l.Latexdocument('results.tex')
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
print('\n\n----------------------------------------------------------------------')
print('--------------- Fit Parameter Skalen Trafo ---------------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Steigung:', ufloat(params[0], errors[0]))
print('y-Abschnitt:', ufloat(params[0], errors[0]))

print('\n----------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('----------------------------------------------------------------------\n\n')

# Plotte Fit and values
plt.plot(index, energies, '.', label='Datenpunkte')
plt.plot(index_intervall, g(index_intervall, *params), label='Fit')

plt.xlabel('Binnummer')
plt.ylabel(r'$\mathrm{Energie} \, / \, \mathrm{eV}$')
plt.legend()
plt.savefig('./results/europium/skalen_trafo_fit.pdf')


# ########################################################################### #
# ################# --- Bestimmung des Effizienz --- ######################## #
# ########################################################################### #

# --- Bestimmung der Winkelverteilung --- #

distance_probe_detektor = 73.1e-3 + 1.5e-2  # Abstand Probe Alu + Dicke Alu
radiant_of_detector = 45e-3 / 2


def winkelverteilung(distance_probe_detektor, radiant_of_detector):
    ''' Function to calculate angle distribution OMEGA '''

    return 1/2 * (1 - (distance_probe_detektor / (np.sqrt(distance_probe_detektor**2 + radiant_of_detector**2))))


angle_distribution = winkelverteilung(distance_probe_detektor,
                                      radiant_of_detector)


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

# --- Definiere Funktion um die Effizienz zu bestimmen --- #


def full_energy_efficiency(measured_decay_rate, angle_distribution,
                           soll_decay_rate, transition_probaility):
    ''' Function to calculate the full_energy_efficiency '''
    return measured_decay_rate / (angle_distribution * soll_decay_rate
                                  * transition_probaility)


# --- Bestimmte die gemessene Aktivität um jeden Peak herum --- ###
'''  Um die Aktivität um jeden Peak herum zu bestimmen, gehe ich vom Peaks
    4 Bins nach rechts und links summiere dann über den Bininhalt  '''


plt.clf()
plt.hist(range(0, len(channel_content_eu), 1),
    bins=np.linspace(0,len(channel_content_eu),len(channel_content_eu)),
    weights=channel_content_eu)

plt.plot(indexes_lower[0], channel_content_eu[indexes_lower[0]], 'x', markersize=0.8)
plt.plot(index_upper_shiftet, channel_content_eu[index_upper_shiftet], 'x', markersize=0.8)
plt.ylim(0,500)
plt.xlim(250, 350)
# plt.show()
plt.savefig('./results/europium/test.pdf')
