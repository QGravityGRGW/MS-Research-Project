#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 15:22:40 2020

@author: jarvis-astro
"""

import numpy as np
import math
import os
#import astropy.units as u
from astropy import units as u
from astropy import constants as const
#from astropy.constants import c, m_e
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib import font_manager as fm, rcParams
import mpdaf
from mpdaf import obj
import scipy.integrate as integrate


# Defining variables

c = 299792458.0 # velocity of light in m/sec
#pc = 3.0856775814671916e+16 # parsec in m
H0 = 0.0    # Hubble Constant
Om = 0.0    # Omega_matter
Ol = 0.0    # Omega_lambda(Omega_vacuum)
Ok = 0.0    
Or = 0.0    # Omega_radiation
tz = 0.5    # Time from z to now in units of 1/H0
t0 = 0.5    # Age of Universe in units of 1/H0
Gyr = (1.0 / (365 * 24 * 60 * 60 * 1e9) )   # Transformation from s to Gyr
yr = (1 / (365 * 24 * 60 * 60) )    # Transformation from s to yr
age_z = 0.1 # Age of Universe at redshift z in units of 1/H0
dCMR = 0.0     # Comoving radial distance in units of c/H0
dCMT = 0.0     # Comoving transverse distance in units of c/H0
dA = 0.0       # Angular size distance
dA_Mpc = 0.0
dA_kpc = 0.0
dA_Gyr = 0.0
dL = 0.0       # Luminosity distance
dL_Mpc = 0.0
dL_kpc = 0.0
dL_Gyr = 0.0
VCM = 0.0      # Comoving volume
DM = 0.0       # Distance modulus

sname = list()    # source name
sname1 = list()
linename = list()
redshift = list()
fl = list()     # calculated flux in units of 1e-20 erg/s/cm^2


# Importing redshift and flux from the file

f = open('/home/jarvis-astro/PROJECT/spectra/SFR_z.txt','r')
for line in f:    
    line = line.strip('"')
    columns = line.split('"')
    sname_val = str(columns[0])
    sname.append(sname_val)
    redshift_val = float(columns[1])
    redshift.append(redshift_val)
#print('Source Name: ', sname)
#print('Redshift: ', redshift)
f.close()

g = open('/home/jarvis-astro/PROJECT/spectra/SFR_flux.txt','r')
for line in g:    
    line = line.strip('"')
    columns = line.split('"')
    sname1_val = str(columns[0])
    sname1.append(sname1_val)
    linename_val = str(columns[1])
    linename.append(linename_val)
    fl_val = float(columns[2])
    fl.append(fl_val)
#print('Source Name: ', sname1)
#print('Line name: ', linename)
#print('Flux: ', fl)
g.close()


# Changing the list to numpy array

sname = np.asarray(sname,dtype=str)
redshift = np.asarray(redshift,dtype=float)
sname1 = np.asarray(sname1,dtype=str)
linename = np.asarray(linename,dtype=str)
fl = np.asarray(fl,dtype=float)
#print('Source Name: ', sname)
#print('Redshift: ', redshift)
#print('Source Name: ', sname1)
#print('Line name: ', linename)
#print('Flux: ', fl)

n = len(redshift)


# Importing spectrum data from the file

a = np.genfromtxt('1-G1_176-145.txt')
#a = np.genfromtxt('4-G3_116-224.txt')
#a = np.genfromtxt('6-G5_180-264.txt')
wav = a[:,0]
wav_a = wav
#print(wav_a)
flux = a[:,1]


# Converting air wavelengths to vacuum wavelengths

for j in range(len(wav)):
    sig = np.empty(len(wav))
    fact = np.empty(len(wav))
    #wav_v = np.empty(len(wav))
    sig[j] = (1e4 / wav[j])**2.
    fact[j] = 1.0 + (6.4328e-5 + (2.94981e-2 / (146. - sig[j])) + (2.5540e-4 / (41. - sig[j])))
    #if not isinstance(wav[j], np.ndarray):
    #    if wav[j] < 2000:
    #        fact = 1.0
    #else:
    #    ignore = np.where(wav[j] < 2000)
    #    fact[ignore] = 1.0
    #wav_v = np.empty(len(wav))
    #wav[j] = wav[j] * (1 + 6.4328e-5 + (2.94981e-2 / (146 - sig[j])) + (2.5540e-4 / (41 - sig[j])))
    wav[j] = wav[j] * fact[j]
    wav_v = wav
#print(wav_v)


# Cosmology Calculator

print("Enter the present value of the Hubble constant (in km/(s.Mpc)) (H0) :")
H0 = float(input(""))
print("Enter the value of the ratio of Matter density of universe to critical densiy (Omega_M) :")
Om = float(input(""))

print("Enter OPEN / FLAT / GENERAL ")
a = str(input(""))

h = H0 / 100
Or = 4.165e-5 / (h * h)

if (a == "GENERAL"):
 	print("Enter the value of the ratio of Vacuum density of universe to critical density (Omega_V) :")
 	Ol = float(input(""))
 	Ok = 1 - Ol - Om - Or
elif (a == "OPEN"):
 	Ol = 0
 	Ok = 1 - Ol - Om - Or
elif (a == "FLAT"):
 	Ok = 0
 	Ol = 1 - Om
else:
 	print("INVALID VALUE ENTERED")
    
def E(z):
    return (Om * ((1 + z) ** 3) + Or * ((1 + z) ** 4) + Ok * ((1 + z) ** 2) + Ol)

dH = c / (H0 * 1e3)               # Hubble distance (in MPC)

def d(z):
    return (dH / (E(z) ** 0.5))

dCMR = np.zeros(n, dtype=float)
for p in range(n):
    dCMR[p] = integrate.quad(d, 0, redshift[p])[0]

dCMT = np.zeros(n, dtype=float)
for q in range(n):
    if (Ok > 0):
        dCMT[q] = (dH / (Ok ** 0.5)) * np.sinh((Ok ** 0.5) * dCMR[q] / dH)
    elif (Ok == 0):
        dCMT[q] = dCMR[q]
    elif (Ok < 0):
        dCMT[q] = (dH / ((-Ok) ** 0.5)) * np.sin(((-Ok) ** 0.5) * dCMR[q] / dH)

# Calculating Angular Diameter(Size) Distance    
dA = np.zeros(n, dtype=float)
for r in range(n):
    dA[r] = dCMT[r] / (1 + redshift[r]) # Mpc
    #print('The angular size distance dA of ' +sname[r] + ' is : ' + '\n' + str(dA[r]) + ' Mpc or ' + str(dA[r]*(3.08568 / (9.461*1e2))) + ' Gly ')    
 
# Calculating Luminosity Distance
dL = np.zeros(n, dtype=float)
for s in range(n):
    dL[s] = ((1 + redshift[s]) ** 2) * dA[s]    # Mpc
    print('The luminosity distance dL of ' +sname[s] + ' is : ' + '\n' + str(dL[s]) + ' Mpc or ' + str(dL[s]*(3.08568 / (9.461*1e2))) + ' Gly ')


# Calculating Luminosity

fl_OII = fl[0:3] * 1e-20
#print(fl_OII)
fl_Halpha = fl[3:6] * 1e-20
#print(fl_Halpha)

linename_OII = linename[0:3]
#print(linename_OII)
linename_Halpha = linename[3:6]
#print(linename_Halpha)

pc = 3.0856775814671916e+18     # in cm
L_Halpha = np.zeros(n, dtype=float)    # Luminosity in erg/s
L_OII = np.zeros(n, dtype=float)
for t in range(n):
    L_Halpha[t] = fl_Halpha[t] * 4 * np.pi * (dL[t] * 1e6 * pc)**2
    print('The Luminosity of ' +linename_Halpha[t] + 'for ' +sname[t] + ' is : ' + '\n' + str(L_Halpha[t]) + ' erg/sec ')
    L_OII[t] = fl_OII[t] * 4 * np.pi * (dL[t] * 1e6 * pc)**2
    print('The Luminosity of ' +linename_OII[t] + 'for ' +sname[t] + ' is : ' + '\n' + str(L_OII[t]) + ' erg/sec ')


# Calculating SFR
SFR_Halpha = np.zeros(n, dtype=float)    # Star Formation Rate in M_sun/year
SFR_OII = np.zeros(n, dtype=float)
err_SFR_OII = np.zeros(n, dtype=float)
L_OII_I = np.zeros(n, dtype=float)
SFR_OII_I = np.zeros(n, dtype=float)
err_SFR_OII_I = np.zeros(n, dtype=float)
for i in range(n):
    SFR_Halpha[i] = 7.9e-42 * L_Halpha[i]
    print('The SFR from' +linename_Halpha[i] + 'luminosity for ' +sname[i] + ' is : ' + '\n' + str(SFR_Halpha[i]) + ' $M_{\odot}$/year ')
    SFR_OII[i] = 1.4e-41 * L_OII[i]
    err_SFR_OII[i] = 0.4e-41 * L_OII[i]
    print('The SFR from' +linename_OII[i] + 'luminosity for ' +sname[i] + ' is : ' + '\n' + str(SFR_OII[i]) + ' $M_{\odot}$/year ')
    print('Error in the SFR from' +linename_OII[i] + 'luminosity for ' +sname[i] + ' is : ' + '\n' + str(err_SFR_OII[i]) + ' $M_{\odot}$/year ')
    
    L_OII_I[i] = 3.11e-20 * ((L_OII[i])**1.495)
    SFR_OII_I[i] = 6.58e-42 * L_OII_I[i]
    err_SFR_OII_I[i] = 1.65e-42 * L_OII_I[i]
    print('The SFR from' +linename_OII[i] + 'luminosity (reddening corrected) for ' +sname[i] + ' is : ' + '\n' + str(SFR_OII_I[i]) + ' $M_{\odot}$/year ')
    print('Error in the SFR from' +linename_OII[i] + 'luminosity (reddening corrected) for ' +sname[i] + ' is : ' + '\n' + str(err_SFR_OII_I[i]) + ' $M_{\odot}$/year ')    



