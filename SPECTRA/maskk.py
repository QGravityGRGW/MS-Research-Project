#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 21:31:09 2020

@author: jarvis-astro
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os
from matplotlib import font_manager as fm, rcParams
import astropy.units as u
import mpdaf
from mpdaf import obj
from mpdaf.obj import Spectrum, WaveCoord, airtovac


#Defining variables

wav = list()    #Air Wavelength
#wav = list()    #Vacuum Wavelength
flux = list()   #Flux
err = list()   #1-sigma error
cont = list()  #continuum

#Importing spectrum data from the file

#f = open('/home/jarvis-astro/spectra/1-G1_176-145.txt', 'r')
#for line in f: 
#    line = line.strip()
#    columns = line.split()
#    wav_val = float(columns[0])
#    wav.append(wav_val)
#    flux_val = float(columns[1])
#    flux.append(flux_val)
#    err_val = float(columns[2])
#    err.append(err_val)
#    print("wavelength: ", wav)
#    print('flux: ', flux)
#f.close()

a = np.genfromtxt('1-G1_176-145.txt')
#a = np.genfromtxt('2-G2_379-175.txt')
#a = np.genfromtxt('4-G3_116-224.txt')
#a = np.genfromtxt('5-G4_83-249.txt')
#a = np.genfromtxt('6-G5_180-264.txt')
#a = np.genfromtxt('7-G6_117-317.txt')
#a = np.genfromtxt('8-G7_265-319.txt')

wav = a[:,0]
wav1 = a[:,0]
flux = a[:,1]
#print("wavelength: ", wav)
#print('flux: ', flux)

spec = Spectrum(filename='/home/jarvis-astro/PROJECT/MUSE Data/1-G1_176-145.fits', ext=0)
spec1 = Spectrum(filename='/home/jarvis-astro/PROJECT/MUSE Data/1-G1_176-145.fits', ext=[0, 1])


#Converting air wavelengths to vacuum wavelengths

# for i in range(len(wav)):
#     #wav = np.empty(len(wav))
#     wav[i] = airtovac(wav[i])
# print(wav)
#print(airtovac(4749.7553710938))
#for i in range(len(wav)):
#    sig = np.empty(len(wav))
#    #wav_v = np.empty(len(wav))
#    sig[i] = 1e4 / wav[i]
#    wav[i] = wav[i] * (1 + ((5.792105 * 1e-2) / (238.0185 - sig[i]**2)) + ((1.67917 * 1e-3) / (57.362 - sig[i]**2)))
#    #wav[i] = wav[i] * (1 + (6.4328 * 1e-5) + ((2.94981 * 1e-2) / (146 - sig[i]**2)) + ((2.5540 * 1e-4) / (41 - sig[i]**2)))
#print(wav)

for j in range(len(wav1)):
    sig = np.empty(len(wav1))
    fact = np.empty(len(wav1))
    sig[j] = (1e4 / wav1[j])**2.
    fact[j] = 1.0 + (6.4328e-5 + (2.94981e-2 / (146. - sig[j])) + (2.5540e-4 / (41. - sig[j])))
    #if not isinstance(wav[j], np.ndarray):
    #    if wav[j] < 2000:
    #        fact = 1.0
    #else:
    #    ignore = np.where(wav[j] < 2000)
    #    fact[ignore] = 1.0
    #wav_v = np.empty(len(wav))
    #wav[j] = wav[j] * (1 + 6.4328e-5 + (2.94981e-2 / (146 - sig[j])) + (2.5540e-4 / (41 - sig[j])))
    wav1[j] = wav1[j] * fact[j]
print(wav1)


#Masking the Spectrum

spec1.mask_region(lmin=5575, lmax=5590, unit=u.angstrom)
spec1.mask_region(lmin=5890, lmax=5915, unit=u.angstrom)
spec1.mask_region(lmin=6298, lmax=6311, unit=u.angstrom)
spec1.mask_region(lmin=7243, lmax=7259, unit=u.angstrom)
spec1.mask_region(lmin=8340, lmax=8350, unit=u.angstrom)
spec1.mask_region(lmin=7274, lmax=7280, unit=u.angstrom)
spec1.mask_region(lmin=7313, lmax=7619, unit=u.angstrom)
spec1.mask_region(lmin=7337, lmax=7343, unit=u.angstrom)
spec1.mask_region(lmin=7365, lmax=7371, unit=u.angstrom)
spec1.mask_region(lmin=7595, lmax=7601, unit=u.angstrom)
spec1.mask_region(lmin=7791, lmax=7797, unit=u.angstrom)
spec1.mask_region(lmin=7819, lmax=7825, unit=u.angstrom)
spec1.mask_region(lmin=7851, lmax=7857, unit=u.angstrom)
spec1.mask_region(lmin=7909, lmax=7915, unit=u.angstrom)
spec1.mask_region(lmin=7960, lmax=7966, unit=u.angstrom)
spec1.mask_region(lmin=7989, lmax=7995, unit=u.angstrom)
# 7277,7316,7340,7368,7598,7750,7794,7822,7854,7912,7963,7992
speccut = spec1.subspec(lmin=math.floor(wav[0]), lmax=math.ceil(wav[len(wav)-1]), unit=u.angstrom) 
speccut.interp_mask()
#speccut.interp_mask(spline=True)


#Plotting the normalized spectra

#fig, ax  = plt.subplots()
#fig=plt.figure(1,figsize=(20,7))

plt.subplots_adjust(bottom=0.12, top=0.78)
f1 = plt.figure(1, figsize=(20,7), dpi=600, frameon=True)
ax = f1.add_subplot(111)
#f2 = plt.figure(2, figsize=(20,7), dpi=600, frameon=True)
#ax2 = f2.add_subplot(111)

plt.rcParams['axes.unicode_minus']=False

#RomanT.tt is for Hershey fonts
fpath = os.path.join(plt.rcParams["datapath"], "/home/jarvis-astro/spectra/RomanT.ttf")
prop = fm.FontProperties(fname=fpath,weight='bold')
#prop1 = fm.FontProperties(fname=fpath,weight='bold',size=16)
fname = os.path.split(fpath)[1]

#plt.rcParams['font.serif'] = "Times New Roman"
#plt.rcParams['font.family'] = "serif"
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
#ax.set_xlim(math.floor(wav[0]), math.ceil(wav[len(wav)-1]))
#ax.set_xlim(math.floor(wav[0]), 7000)
#ax.set_xlim(8600,math.ceil(wav[len(wav)-1]))
ax.set_xlim(8000, 9000)
##plt.gca().set_xlim(right=6400)
#ax.set_ylim(-0.2,2.25)
#plt.title('Spectrum of G16 (RA = 17.734525$\mathregular{^{o}}$, DEC = -2.200093$\mathregular{^{o}}$)',fontsize=20)
ax.set_xlabel('Observed Wavelength ($\mathregular{\AA}$)',fontproperties=prop,fontsize=20,fontweight='heavy')
ax.set_ylabel('Flux ($\mathregular{10^{-20}}$ erg s$\mathregular{^{-1}}$ cm$\mathregular{^{-2}}$ $\mathregular{\AA}\mathregular{^{-1}}$)',fontproperties=prop,fontsize=20,fontweight='heavy')
#ax.set_ylabel('Flux (1e-20 erg Ang$^{-1}$ cm$^{-2}$ s$^{-1}$)')
#ax.set_xlabel('Observed Wavelength ($\mathregular{\AA}$)')
ax.axhline(y=0,color='grey',ls='--')
#ax.set_title('Spectral Analysis before interpolation')

#ax2.set_xlim(math.floor(wav[0]), math.ceil(wav[len(wav)-1]))
#ax2.set_xlabel('Observed Wavelength ($\mathregular{\AA}$)',fontproperties=prop,fontsize=20,fontweight='heavy')
#ax2.set_ylabel('Flux ($\mathregular{10^{-20}}$ erg s$\mathregular{^{-1}}$ cm$\mathregular{^{-2}}$ $\mathregular{\AA}\mathregular{^{-1}}$)',fontproperties=prop,fontsize=20,fontweight='heavy')
#ax2.axhline(y=0,color='grey',ls='--')
#ax2.set_title('Spectral Analysis after interpolation')

#ax.step(wav, flux)
##ax.bar(l, y, align='center')
#ax.plot(wav, flux, drawstyle='steps-mid', c='b')  # If plotting from already given galaxy spectrum
#spec.plot() #If lotting from MUSE cube
#speccut.plot()
#ax.plot(wav,err*10**17,drawstyle='steps-mid',c='b')

#spec1.unmask()
plt.figure(2)
#spec1.plot(lmin=math.floor(wav[0]), lmax=math.ceil(wav[len(wav)-1]), title='Spectrum before interpolation', unit=u.angstrom)
spec.plot(title='Spectrum before interpolation', unit=u.angstrom)
plt.figure(3)
speccut.plot(title='Spectrum after interpolation', unit=u.angstrom)

plt.show()

