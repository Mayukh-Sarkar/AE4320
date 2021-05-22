# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:16:20 2021

@author: Mayukh
"""

import scipy.io
import scipy as sc
import numpy as np
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
from findpeaks import findpeaks


data = scipy.io.loadmat('E:\MSc Aero-20200916T124909Z-001\MSc Aero\Q3\system identification\AE4320\AE4320_dataset1.mat')
e = scipy.io.loadmat('E:\MSc Aero-20200916T124909Z-001\MSc Aero\Q3\system identification\AE4320\error.mat')
u = scipy.io.loadmat('E:\MSc Aero-20200916T124909Z-001\MSc Aero\Q3\system identification\AE4320\input.mat')


f_t = np.transpose(data["ft"])
f_d = data["fd"]
t = np.transpose(data["t"])
e_id = e["error"]
u_id = u["input"]

plt.figure(1)

ft = findpeaks(lookahead=1)
results = ft.fit(f_t)
ft.plot()

