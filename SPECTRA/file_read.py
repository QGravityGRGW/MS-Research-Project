#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 00:22:09 2021

@author: jarvis-astro
"""

import numpy as np

ionname = list()
data_col = list()
sigma_col = list()
f = open('/home/jarvis-astro/cloudy_run/ion list.txt','r')
for line in f:
    line = line.strip('"')
    columns = line.split('"')
    ionname_val = str(columns[0])
    ionname.append(ionname_val)
    data_col_val = float(columns[1])
    data_col.append(data_col_val)
    sigma_col_val = float(columns[2])
    sigma_col.append(sigma_col_val)
    #print('Ion Name: ', ionname)
    #print('Column Density: ', data_col)
    #print('Error in Column Density: ', sigma_col)
f.close() 
#ionname = np.asarray(ionname,dtype=str)
#data_col = np.asarray(data_col,dtype=float)
#sigma_col = np.asarray(sigma_col,dtype=float)
print('Ion Name: ', ionname)
print('Column Density: ', data_col)
print('Error in Column Density: ', sigma_col)