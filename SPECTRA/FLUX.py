#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:11:37 2020

@author: jarvis-astro
"""

import numpy as np
import math
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib import font_manager as fm, rcParams
import mpdaf
from mpdaf import obj
import scipy.integrate as integrate


# Extracting Cube and subcube

cube = obj.Cube('/home/jarvis-astro/PROJECT/MUSE Data/ADP.2017-03-24T13:24:09.696.fits')
ra = obj.hms2deg('02:09:31.298')  #Enter RA in hours, minutes and seconds
dec = obj.dms2deg('-04:38:15.80') #Enter DEC in degress, minutes and seconds
size = 10 #Enter circular radius
sub = cube.subcube((dec,ra), size)
img = sub.sum(axis=0)
spec = sub.sum(axis=(1,2))

#G1	02:09:31.351	-04:38:39.80   9
#G3	02:09:32.154	-04:38:24.00   10
#G5 02:09:31.298	-04:38:15.80   10


# Defining variables

wav = list()    #Wavelength
flux = list()   #Flux
err = list()   #1-sigma error
cont = list()  #continuum
#z = 0.39045 #Redshift of the galaxy (Change here!!!!!!)

# Importing spectrum data from the file

#a = np.genfromtxt('1-G1_176-145.txt')
#a = np.genfromtxt('4-G3_116-224.txt')
a = np.genfromtxt('6-G5_180-264.txt')

wav = a[:,0]
wav_a = a[:,0]
#print(wav_a)
flux = a[:,1]
#print("wavelength: ", wav)
#print('flux: ', flux)


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


# Plotting the normalized spectra

#fig, ax  = plt.subplots()
#fig=plt.figure(1,figsize=(20,7))

plt.subplots_adjust(bottom=0.2)
f1 = plt.figure(1, figsize=(20,7))
ax = f1.add_subplot(111)
#ax.step(wav, flux)
##ax.bar(l, y, align='center')
#ax.set_xlim(4750, 5020)
##plt.gca().set_xlim(right=6400)
#ax.set_ylim(500, 7000)
#ax.set_ylabel('Flux (1e-20 erg Ang$^{-1}$ cm$^{-2}$ s$^{-1}$)')
#ax.set_xlabel('Observed Wavelength ($\mathregular{\AA}$)')
#ax.set_title('Spectral Analysis')
#ax.set_xlim(1140, 1180)

plt.rcParams['axes.unicode_minus']=False

# RomanT.tt is for Hershey fonts
fpath = os.path.join(plt.rcParams["datapath"], "/home/jarvis-astro/PROJECT/spectra/RomanT.ttf")
prop = fm.FontProperties(fname=fpath, weight='bold')
prop1 = fm.FontProperties(fname=fpath, weight='bold', size=16)
fname = os.path.split(fpath)[1]

#plt.rcParams['font.serif'] = "Times New Roman"
#plt.rcParams['font.family'] = "serif"
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
#ax1=plt.subplot(111)
#ax.set_title('Spectrum of G1 (RA = 17.734525$\mathregular{^{o}}$, DEC = -2.200093$\mathregular{^{o}}$)', fontsize=20)
ax.set_xlabel('Observed Wavelength ($\mathregular{\AA}$)', fontproperties=prop, fontsize=20, fontweight='heavy')
ax.set_ylabel('Flux ($\mathregular{10^{-20}}$ erg s$\mathregular{^{-1}}$ cm$\mathregular{^{-2}}$ $\mathregular{\AA}\mathregular{^{-1}}$)', fontproperties=prop, fontsize=20, fontweight='heavy')
ax.axhline(y=0, color='grey', ls='--')
ax.tick_params(direction='in', which='major', length=8)
ax.tick_params(which='minor', direction='in', length=4)

ax.step(wav_v, flux)
ax.set_xlim(7285, 7310)
#ax.set_xlim(math.floor(wav_v[0]), math.ceil(wav_v[len(wav_v)-1]))
#ax.set_xlim(math.floor(wav_v[0]), 5000)
#ax.set_xlim(9000, math.ceil(wav_v[len(wav_v)-1]))
#ax.set_xlim(math.floor(wav_v[0]), math.floor(wav_v[0]) + 500)
##plt.gca().set_xlim(right=6400)
#plt.gca().set_ylim(top=math.ceil(flux[len(flux)-1]))
#ax.set_ylim(-0.2, 2.25)
#ax.set_ylim(math.floor(flux[0]), math.ceil(flux[len(flux)-1]))
#ax.plot(wav, flux, drawstyle='steps-mid', c='b')  # If plotting from already given galaxy spectrum
#spec.plot() #If lotting from MUSE cube
#ax.plot(wav_v, err*10**17, drawstyle='steps-mid', c='b')


# Calculating Flux by clicking to get λmin, λmax and dλ

b = int(input("Enter width: "))

def mean(p, q, r):
    s = 0
    for i in range(q, (r + 1)):
        s += p[i]
    return s / (r + 1 - q)

def FLUX(point, wav, flux):
    i = s = 0
    while (point[0][0] > wav[i]):
        i += 1
    e = i
    while (point[1][0] > wav[e]):
        e += 1
#    x3 = point[2][0]
#    x4 = point[3][0]
    #print(l[e], l[i])
    s += mean(flux, i, e) * (wav[e] - wav[i])
#    s += mean(flux, i, e) * (x4 - x3)
    return s

def next(event):
    x = ax.get_xlim()
    ax.set_xlim(x[0]+b,x[1]+b)
#   print('next')
    plt.draw()
    
def prev(event):
    x = ax.get_xlim()
    ax.set_xlim(x[0]-b,x[1]-b)
#   print('prev')
    plt.draw()
    
def onclick(event):
#    global ix, iy
#    ix, iy=event.xdata,event.ydata
#    print(ix,iy)
    Flux = 0.0
    Flux_mpdaf = 0.0
    Flux_int = 0.0
    x_ip = plt.ginput(2)
#    x1 = [x1[0] for x1 in x_ip]
#    x2 = [x2[0] for x2 in x_ip]
    x1 = x_ip[0][0]
    x2 = x_ip[1][0]
#    x3 = x_ip[2][0]
#    x4 = x_ip[3][0]
    print("The coordinates are: ", x_ip)
    print("The first x-point is: ", x1)
    print("The second x-point is: ", x2)
#    print("The third x-point is: ", x3)
#    print("The fourth x-point is: ", x4)
#    a = np.genfromtxt('1-G1_176-145.txt')
#    a = np.genfromtxt('4-G3_116-224.txt')
#    wav = a[:,0]
#    wav_a = wav
#    flux = a[:,1]
#    Flux = FLUX(x_ip, wav, flux)
#    print('Flux is: ' + str(Flux) + ' ' + 'e-20' + ' ' +'erg/s/cm\u00b2')
    Flux_mpdaf = spec.integrate(x1, x2)
    print('Flux_mpdaf is: ' + str(Flux_mpdaf))
    #Flux_int = integrate.quad(flux, x1, x2)
#    print('Flux is: ' + str(Flux_int) + ' ' + 'e-20' + ' ' +'erg/s/cm\u00b2')
#    ('The luminosity distance dL of ' +sname[s] + ' is : ' + '\n' + str(dL[s]) + ' Mpc or ' + str(dL[s]*(3.08568 / (9.461*1e2))) + ' Gly ')
    #with open("fluxG1.txt", "w") as output:
    #    output.write(str(Flux_mpdaf))

#Flux_mpdaf = 0.0
#Flux_mpdaf = spec.integrate(4982.391961280496,4994.897182025895)
#print('1. Flux_mpdaf is: ' + str(Flux_mpdaf))

    
#cid = f1.canvas.mpl_connect('button_press_event', onclick)
axprev = plt.axes([0.7, 0.01, 0.08, 0.05])
axnext = plt.axes([0.8, 0.01, 0.08, 0.05])
axonclick = plt.axes([0.07, 0.02, 0.25, 0.05])
bnext = Button(axnext, 'Next')
bnext.on_clicked(next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(prev)
clck = Button(axonclick, 'Click here for selection of the points ')
clck.on_clicked(onclick)




plt.show()