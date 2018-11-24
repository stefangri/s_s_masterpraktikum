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

#peak_energy = g(peak_indexes, m, b)
# ########################################################################### #
# ###################### --- Untersuche den Photopeak --- ################### #
# ########################################################################### #

# --- Definiere Gauß-Funktion um Photopeak zu fitten --- #


def gaus(x, amplitude, sigma, offset):
    return amplitude * np.exp(-1/2 * (x - offset)**2/sigma**2)


index = np.arange(peak_indexes[0]-15, peak_indexes[0]+15)


photo_peak_energy = g(index, m, b)

print()
params_gaus, cov_gaus = curve_fit(gaus, noms(photo_peak_energy),
                                  channel_content_ca[peak_indexes[0][0] - 15:
                                                     peak_indexes[0][0] + 15],
                                  p0=[channel_content_ca[peak_indexes[0]][0],
                                      0.9, 661])

errors_gaus = np.sqrt(np.diag(cov_gaus))

print('\n\n------------------------------------------------------------------')
print('----------- Fit Parameter Photpeakvermessung -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Amplitude/E_gamma:', ufloat(params_gaus[0], errors_gaus[0]))
amplitude = ufloat(params_gaus[0], errors_gaus[0])

sigma = ufloat(params_gaus[1], errors_gaus[1])
print('Halbertsbreite:', ufloat(params_gaus[1], errors_gaus[1]))

print('Offset/bins:', ufloat(params_gaus[2], errors_gaus[2]))
offset = ufloat(params_gaus[2], errors_gaus[2])

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

# --- Bestimmung der Zehntel Breite --- #

index = np.arange(0, len(channel_content_ca), 1)
index_zehntel = index[np.abs(channel_content_ca - 1 / 10 * params_gaus[0])
                      < 100]
index_zehntel = g(index_zehntel, m, b)

index_halb = index[np.abs(channel_content_ca - 1 / 2 * params_gaus[0]) < 200]
index_halb = g(index_halb, m, b)

halbebreite_berechnet = 2 * np.sqrt(2 * np.log(2)) * sigma
zehntelbreite_berechnet = 1.823 * halbebreite_berechnet


print('\n\n------------------------------------------------------------------')
print('----------- Halbwerts- und zehntelbreite -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Halbwertsbreite gemessen: Energie[keV]',
      2 * np.abs(offset-index_halb), '\n')

print('Halbwertsbreite berechnet: Channel, Energie[keV]',
      halbebreite_berechnet, '\n')

print('Zehntelbreite gemessen: Channel, Energie[keV]',
      2 * np.abs(offset-index_zehntel), '\n')

print('Zehntelbreite berechnet aus Sigma: Channel, Energie[keV]',
      zehntelbreite_berechnet, '\n')

print('Zehntelbreite berechnet aus Messwert: Channel, Energie[keV]',
     2* np.abs(offset-index_halb) * 1.823)

print('Verhältnis Halbwerts- zu Zehntelbreite gemessen:',
       np.abs(offset-index_zehntel) / np.abs(offset-index_halb) , '\n')

print('Verhältnis Halbwerts- zu Zehntelbreite berechnet:',
        zehntelbreite_berechnet / halbebreite_berechnet, '\n')

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

print('\n\n------------------------------------------------------------------')
print('--------- Vergleiche Halbwertsbreite mit Delta E Formel 20 -----------')
print('----------------------------------------------------------------------',
      '\n')

print('Halbwertsbreite gemessen: Channel, Energie[keV]',
      np.abs(offset-index_halb), g(np.abs(offset-index_halb),  m, b))

print('Halbwertsbreite berechnet: Channel, Energie[keV]',
      halbebreite_berechnet, g(halbebreite_berechnet,  m, b), '\n')

print('Auflösungsvermögen nach Formel 20: Channel, Energie[keV]',
      2.35 * unp.sqrt(0.1 * offset * 1e3 * 2.9) * 1e-3)

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# --- Bestimmte die Absoprtionswahrcheinlichkeit --- #

# mu aus Abbilung bei einer Energie von  661 kev abgelsen
mu_compton = ufloat(0.35, 0.05)  # in 1/cm
mu_photo = ufloat(0.06, 0.005)
length_of_detector = 3.9  # in cm

print('\n\n------------------------------------------------------------------')
print('--------- Absorptionswahrcheinlichkeit -----------')
print('----------------------------------------------------------------------',
      '\n')

print('Wahrscheinlichkeit Compton:',
      (1 - unp.exp(-mu_compton * length_of_detector)) * 100)

print('Wahrscheinlichkeit Photo:',
      (1 - unp.exp(-mu_photo * length_of_detector)) * 100)

print('Verhältnis Compton/Photo',
     (1 - unp.exp(-mu_compton * length_of_detector) * 100)/(1 - unp.exp(-mu_photo * length_of_detector) * 100))


print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')



index = np.linspace(noms(photo_peak_energy[0]), noms(photo_peak_energy[-1])+5,
                    1e3)

plt.xlim(noms(photo_peak_energy[0])-5, noms(photo_peak_energy[-1])+5)

plt.hist(noms(g(np.arange(0, len(channel_content_ca), 1), m, b)),
         bins=noms(g(np.linspace(0, len(channel_content_ca),
          len(channel_content_ca)), m, b)),
         weights=channel_content_ca, label='Spektrum')

plt.plot(index, gaus(index, *params_gaus),
         label='Fit')

plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.ylabel(r'$\mathrm{Zählung}$')
# plt.plot(index_halb, channel_content_ca[index_halb], '.',
         # label='Abgelesene Höhe\nHalbwertsbreite')
# plt.plot(index_zehntel, channel_content_ca[index_zehntel], '.',
         # label='Abgelesene Höhe\nZehntelbreite')

plt.legend()
plt.savefig('./plots/caesium/photopeak.pdf')


# ########################################################################### #
# ################# --- Untersuchung des Comptpon-Spektrums --- ############ #
# ########################################################################### #

# --- Suche den Rückstreu und "Compton-Kanten-Peak" --- #

# Habe den Mittelpunktes des Peaks erst mit PeakDetection abgeschätzt und dann
# aus dem Histogramm abgelesen, hierbei wollte ich den Mittelpunkt des Peaks
# treffen. Bin analog für die Compton Kante vorgenagen

index_peak_rückstreu = 480
index_compton_kante = 1182

# --- Transformiere die obgien Indicies in full_energy_efficiency --- #

energy_peak_rückstreuung = g(ufloat(index_peak_rückstreu, 5),  m, b)
energy_compton_kante = g(ufloat(index_compton_kante, 5),  m, b)

print('\n\n------------------------------------------------------------------')
print('----------- Lage Rückstreupeak und Comptonkante ---------------------')
print('----------------------------------------------------------------------',
      '\n')

print('RückstreuPeak: Index, Energie[keV]:', index_peak_rückstreu,
      energy_peak_rückstreuung)

print('Compton_peak: Index, Energie[keV]:', index_compton_kante,
      energy_compton_kante)

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# --- Vergleiche Comptonkante mit Formel 9 und ?? --- #


def compton_edge_energy(E_gamma):
    # Using mass of the electron directly, beaucse this is the transformed
    # Ruheenergie
    m_e = ufloat(const.physical_constants["electron mass energy equivalent in MeV"][0],
                 const.physical_constants["electron mass energy equivalent in MeV"][2])
    m_e *= 1e3
    epsilon = E_gamma / (m_e)
    return E_gamma * ((2 * epsilon) / (1 + 2 * epsilon))

def rueckstreu_energy(E_gamma):
    # Using mass of the electron directly, beaucse this is the transformed
    # Ruheenergie
    m_e = ufloat(const.physical_constants["electron mass energy equivalent in MeV"][0],
                 const.physical_constants["electron mass energy equivalent in MeV"][2])
    m_e *= 1e3
    epsilon = E_gamma / (m_e)
    return E_gamma * (1 / (1 + 2 * epsilon))

print('\n\n------------------------------------------------------------------')
print('----------- Berchneten Comptonedge energie ---------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Compton_edge berechnet: Energie[keV]:',
      compton_edge_energy(offset))

print('Rückstreupeak berechnet: Energie[keV]:',
      rueckstreu_energy(offset))

print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')




# --- Bestimme den Flächeninhalt --- #
# Indices wurden aus der Grafik abgelsen

lower_index_compton = 51
inhalt_compton_spektrum = sum(channel_content_ca[lower_index_compton
                              :index_compton_kante])

print('\n\n------------------------------------------------------------------')
print('------------------------ Inhalt Compton Fläche -----------------------')
print('----------------------------------------------------------------------',
      '\n')

print('Inhalt Compton Spektrum:', inhalt_compton_spektrum)

print('Verhältnis Photo/Compton', inhalt_compton_spektrum/area_photo_peak)
print('\n--------------------------------------------------------------------')
print('----------------------------------------------------------------------')
print('------------------------------------------------------------------\n\n')

# --- Plotte den Compton Bereich --- #
plt.clf()
# index = np.linspace(peak_indexes[0]-20, peak_indexes[0]+20, 1e3)

plt.xlim(0, noms(g(1400, m, b)))
# plt.ylim(0, 400)
plt.hist(noms(g(np.arange(0, len(channel_content_ca), 1), m, b)),
         bins=noms(g(np.linspace(0, len(channel_content_ca),
                                 len(channel_content_ca)), m, b)),
         weights=channel_content_ca, label='Spektrum')

plt.plot(noms(g(index_peak_rückstreu, m, b)), channel_content_ca[index_peak_rückstreu], '.',
         label='Rückstreupeak')
plt.plot(noms(g(index_compton_kante, m, b)), channel_content_ca[index_compton_kante], '.',
         label='Compton Kante')
# plt.plot(lower_index_compton, channel_content_ca[lower_index_compton], '.',
# label='Beginn Compton Spektrum')

plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.ylabel(r'$\mathrm{Zählung}$')
plt.legend()
plt.savefig('./plots/caesium/compton.pdf')

# ##########################################
# ## --- Plotte das gesamte Spektrum --- ###
# ##########################################
plt.clf()
# index = np.linspace(peak_indexes[0]-20, peak_indexes[0]+20, 1e3)

plt.yscale('log')
plt.xlim(0, noms(g(1800, m, b)))
# plt.ylim(0, 400)
plt.hist(noms(g(np.arange(0, len(channel_content_ca), 1), m, b)),
         bins=noms(g(np.linspace(0, len(channel_content_ca),
                                 len(channel_content_ca)), m, b)),
         weights=channel_content_ca,
         label='Spektrum')

plt.xlabel(r'$\mathrm{Energie}\, / \, keV$')
plt.ylabel(r'$\mathrm{Zählung}$')
plt.legend()
plt.savefig('./plots/caesium/caesium_spektrum.pdf')
