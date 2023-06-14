#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys
from scipy.optimize import curve_fit
import cmath
import pickle

def mean(angles, deg=True):
    '''Circular mean of angle data
    '''
    a = np.deg2rad(angles) if deg else np.array(angles)
    angles_complex = np.frompyfunc(cmath.exp, 1, 1)(a * 1j)
    mean = cmath.phase(angles_complex.sum()) % (2 * np.pi)
    return round(np.rad2deg(mean) if deg else mean, 4) 


def corrcoef(x, y, deg=True):
    '''Circular correlation coefficient of two angle data
    '''    
    convert = np.pi / 180.0 if deg else 1
    sx = np.frompyfunc(np.sin, 1, 1)((x - mean(x, deg)) * convert)
    sy = np.frompyfunc(np.sin, 1, 1)((y - mean(y, deg)) * convert)
    r = (sx * sy).sum() / np.sqrt((sx ** 2).sum() * (sy ** 2).sum())
    
    return round(r, 4)

if __name__ == "__main__":
    with open("result_table2.pickle",mode="rb") as fi: data = pickle.load(fi)

    c_arr = np.array(data["c"])
    strain_azim  = np.zeros((c_arr.shape[0]))
    for j in range(c_arr.shape[0]):
        if c_arr[j]<0:
            strain_azim[j] = c_arr[j] + 180.
        else:
            strain_azim[j] = c_arr[j]
    strain_azim = 180. - strain_azim
    station_azim = np.array(data["station-azim"])

    cc = corrcoef(station_azim*2, strain_azim*2)
    print ("correlation coefficient: ",cc)


