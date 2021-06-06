#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 15:22:40 2020

@author: jarvis-astro
"""

import numpy as np
import math
#import astropy.units as u
from astropy import units as u
#from astropy.constants import c, m_e
from astropy.coordinates import Angle
import mpdaf
from mpdaf import obj
import scipy.integrate as integrate


# Defining variables

c = 299792458.0 #velocity of light in m/sec
H0 = 0.0    #Hubble Constant
Om = 0.0    #Omega_matter
Ol = 0.0    #Omega_lambda(Omega_vacuum)
Ok = 0.0    
Or = 0.0    #Omega_radiation
tz = 0.5    #Time from z to now in units of 1/H0
t0 = 0.5    #Age of Universe in units of 1/H0
Gyr = (1.0 / (365 * 24 * 60 * 60 * 1e9) )   #transformation from s to Gyr
yr = (1 / (365 * 24 * 60 * 60) )    #Transformation from s to yr
age_z = 0.1 #Age of Universe at redshift z in units of 1/H0
dCMR = 0.0     #Comoving radial distance in units of c/H0
dCMT = 0.0     #Comoving transverse distance in units of c/H0
dA = 0.0       #Angular size distance
dA_Mpc = 0.0
dA_kpc = 0.0
dA_Gyr = 0.0
dL = 0.0       #Luminosity distance
dL_Mpc = 0.0
dL_kpc = 0.0
dL_Gyr = 0.0
VCM = 0.0      #Comoving volume
DM = 0.0       #Distance modulus

sname = list()    #source name
redshift = list()
RA = list() #in hours, minutes and seconds
DEC = list()    #degress, minutes and seconds
#Quasar 02:09:30.780	-04:38:26.70
#G1 02:09:31.351	-04:38:39.80
#G3 02:09:32.154	-04:38:24.00


# Importing redshift, RA and DEC data from the file

#f = open('/home/jarvis-astro/PROJECT/spectra/vac/final/z_RA_Dec.txt','r')
f = open('/home/jarvis-astro/PROJECT/spectra/vac/final/z_RA_Dec_G0.txt','r')
for line in f:    
    line = line.strip('"')
    columns = line.split()
    sname_val = str(columns[0])
    sname.append(sname_val)
    redshift_val = float(columns[1])
    redshift.append(redshift_val)
    RA_val = str(columns[2])
    RA.append(RA_val)
    DEC_val = str(columns[3])
    DEC.append(DEC_val)
#print('Source Name: ', sname)
#print('Redshift: ', redshift)
#print('RA: ', RA)
#print('DEC: ', DEC)
f.close()


# Changing the list to numpy array

sname = np.asarray(sname,dtype=str)
redshift = np.asarray(redshift,dtype=float)
RA = np.asarray(RA,dtype=str)
DEC = np.asarray(DEC,dtype=str)
print('Source Name: ', sname)
print('Redshift: ', redshift)
print('RA: ', RA)
print('DEC: ', DEC)


#Converting RA and DEC from string to float

#for i in range(len(RA)):
#    a = RA[i]
#    print(a)
#for j in range(len(DEC)):
#    b = DEC[j]
#    print(b)

#for i in range(len(RA)):
#    ra = float(a[i])
    #dec = float(b[i])
#print('RA: ', ra)
#print('DEC: ', dec)

# Converting RA and DEC in degree and radian

n = len(redshift)
ra = np.zeros(n, dtype=float)
ra_rad = np.zeros(n, dtype=float)
for i in range(n):
    ra[i] = obj.hms2deg(RA[i])  #in degree
    ra_rad[i] = ra[i] * (np.pi / 180)  #in radian 
    print('RA in degree of ' +sname[i], ra[i])
    print('RA in radian of ' +sname[i], ra_rad[i])

dec = np.zeros(n, dtype=float)
dec_rad = np.zeros(n, dtype=float)
for j in range(n):
    dec[j] = obj.dms2deg(DEC[j])    #in degree
    dec_rad[j] = dec[j] * (np.pi / 180)    #in radian 
    print('DEC in degree of ' +sname[j], dec[j])
    print('DEC in radian of ' +sname[j], dec_rad[j])


# Velocity difference between systems at two redshifts along the same LOS

def veldif(z1, z2):
    delv = c * ((1 + z1)**2 - (1 + z2)**2) / ((1 + z1)**2 + (1 + z2)**2)
    return delv * 1e-3      #in km/s


# Angular Separation between two point sources

def angsep(RA1, RA2, DEC1, DEC2):
    rad = math.acos((math.cos(np.pi/2 - DEC1) * math.cos(np.pi/2 - DEC2)) + (math.sin(np.pi/2 - DEC1) * math.sin(np.pi/2 - DEC2) * math.cos(RA1 - RA2)))    #in radian
    deg = rad * (180 / np.pi)   #in degree
    return deg


# Calculation
delv = np.zeros(n, dtype=float)
for k in range(1, n):
    delv[k] = veldif(redshift[k], redshift[0])
    print('Velocity difference between Quasar and Galaxy ' +sname[k] + ' is : ' + '\n' + str(delv[k]) + ' km/s ')
    delv1 = np.trim_zeros(delv)

delth = np.zeros(n, dtype=float) 
delth_arcsec = np.zeros(n, dtype=float)   
for l in range(1, n):
    delth[l] = angsep(ra_rad[l], ra_rad[0], dec_rad[l], dec_rad[0]) #in degree
    delth_arcsec[l] = delth[l] * 60 * 60    #in arcsec
    print('Angular Separation between Quasar and Galaxy ' +sname[l] + ' is : ' + '\n' + str(delth[l]) + ' degree ')
    print('Angular Separation between Quasar and Galaxy ' +sname[l] + ' is : ' + '\n' + str(delth_arcsec[l]) + ' arcsec ')
    delth_deg = np.trim_zeros(delth)
    delth_arcsec1 = np.trim_zeros(delth_arcsec)

m = len(delth_deg)
    
#rad1 = math.acos((math.cos(np.pi/2 - (-0.08105989724669943)) * math.cos(np.pi/2 - (-0.08099636726192684))) + (math.sin(np.pi/2 - (-0.08105989724669943)) * math.sin(np.pi/2 - (-0.08099636726192684)) * math.cos((0.5651480837712759) - (0.5651061958692281))))
#deg1 = rad1 * (180 / np.pi)    
##print(rad1)
#print(deg1)
#rad2 = math.acos((math.cos(np.pi/2 - (-0.08098327729253689)) * math.cos(np.pi/2 - (-0.08099636726192684))) + (math.sin(np.pi/2 - (-0.08098327729253689)) * math.sin(np.pi/2 - (-0.08099636726192684)) * math.cos((0.5652074249658436) - (0.5651061958692281))))
#deg2 = rad2 * (180 / np.pi)    
##print(rad2)
#print(deg2)


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

dH = c / (H0 * 1e3)               #Hubble distance (in MPC)

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
dA_Mpc = np.zeros(n, dtype=float)
dA_kpc = np.zeros(n, dtype=float)
for r in range(n):
    dA[r] = dCMT[r] / (1 + redshift[r])
    #dA_Mpc[r] = (c / H0) * dA[r]
    dA_kpc[r] = dA[r] / 206.264806
    print('The angular size distance dA of ' +sname[r] + ' is : ' + '\n' + str(dA[r]) + ' Mpc or ' + str(dA[r]*(3.08568 / (9.461*1e2))) + ' Gly ')
    #print('This gives a scale for ' +sname[r] + '%.2f' % dA_kpc[r] + ' kpc/" ')
    print('This gives a scale for ' +sname[r] + ' of : ' + '\n' + str(dA_kpc[r]) + ' kpc/" ')
    #dA = np.round(dA,6)
    #print('The angular size distance dA of ' +sname[r] + str(np.round(dA[r],3)) + ' Mpc or ' + str(np.round(dA[r]*(3.08568 / (9.461*1e2)),3)) + ' Gly ')
    
 
# Calculating Luminosity Distance
dL = np.zeros(n, dtype=float)
#dL_Mpc = np.zeros(n, dtype=float)
#dL_kpc = np.zeros(n, dtype=float)
for t in range(n):
    dL[t] = ((1 + redshift[t]) ** 2) * dA[t]
    #dA_Mpc[r] = (c / H0) * dA[r]
    #dL_kpc[t] = dL[t] / 206.264806
    print('The luminosity distance dL is ' +sname[t] + ' is : ' + '\n' + str(dL[t]) + ' Mpc or ' + str(dL[t]*(3.08568 / (9.461*1e2))) + ' Gly ')


# Calculating Physical Separation by Cosine Rule

#def physep(a, b, th):
#    c = np.sqrt(a**2 + b**2 - (2 * a * b * math.cos(th)))
#    return c

# delth_rad = delth_deg * (np.pi / 180)
#delth_rad = delth * (np.pi / 180)
physep_d = np.zeros(n, dtype=float)
#for  s in range(1, n):
#    physep_d[s] = physep(dA[s], dA[0], delth_rad[s])
#    print('The physical separation distance between Quasar and Galaxy ' +sname[s] + ' is : ' + '\n' + str(physep_d[s]) + ' Mpc or ' + str(physep_d[s]*(3.08568 / (9.461*1e2))) + ' Gly ')
#    physep_d1 = np.trim_zeros(physep_d)


for  s in range(1, n):
    physep_d[s] = dA_kpc[0] * delth_arcsec[s]
    print('The physical separation distance in absorber redshift between Quasar and Galaxy ' +sname[s] + ' is : ' + '\n' + str(physep_d[s]) + ' kpc ')
    physep_d1 = np.trim_zeros(physep_d)

