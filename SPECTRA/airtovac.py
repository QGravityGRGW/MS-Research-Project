#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 20:00:39 2021

@author: jarvis-astro
"""

import numpy as np
import mpdaf
from mpdaf import obj

#Defining variables

wav_emm = list()  #emission line wavelength


# Importing emission line data from the file

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

a = np.genfromtxt('emission_lines_air.txt')

wav = a[:,3]
wav_a = a[:,3]
#print("wavelength: ", wav_a)



# Converting air wavelengths to vacuum wavelengths

for j in range(len(wav)):
    sig = np.empty(len(wav))
    fact = np.empty(len(wav))
    #wav_v = np.empty(len(wav))
    if wav[j]>2000:
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
#print('Calculated vacuum wavelength: ', wav_v)

#for i in range(len(wav_a)):
#    wav_v1 = np.empty(len(wav_a))
#    if wav_a[i]>2000:
#        wav_v1[i] = obj.airtovac(wav_a[i])
#print('Vacuum wavelength from mpdaf: ', wav_v1)

#print('Error: ', wav_v1 - wav_v)

with open("file1.txt", "w") as output:
    output.write(str(wav_v))
    
#with open("file2.txt", "w") as output:
#    output.write(str(wav_v1))