#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 13:07:33 2021

@author: jarvis-astro
"""

import numpy as np
import math

c = 299792458.0 #velocity of light in m/sec
z2 = 0.39047    #absorber redshift
delv1 = 1.7 * 1e3
delv2 = -121 * 1e3

def redshift(delv, z2):
    z1 = np.sqrt(((delv + c) * (1 + z2)**2) / (c - delv))
    return z1 - 1

print('Redshift 1: ', redshift(1.7 * 1e3, z2))
print('Redshift 2: ', redshift(-121 * 1e3, z2))
print('Redshift 3: ', redshift(-147 * 1e3, z2))
