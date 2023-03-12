#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fitting-dat.py

@author: tomoya
"""
import numpy as np
import sys
import math as m
import matplotlib.pyplot as plt
import pyproj
from scipy.optimize import curve_fit


mu = 10.0 ## Rigidity [GPa]
lamb = mu ## Lame constant under poisson's ratio of 0.25
rho = 2700.0 ## Density [km/m3]
K = (3.0*lamb+2.0*mu)/3.0 ## bulk modulus [GPa]
S11 = (lamb + mu)/(mu*(2.0*mu + 3.0*lamb))
S12 = -1.0*lamb / (2.0*mu*(2.0*mu+3.0*lamb))
E   = (mu*(3.*lamb + 2.*mu))/(lamb + mu)

def sinx2(x, x0, a, b):
    y = a*np.sin(x-x0)*np.sin(x-x0)+b
    return y

def funcP(x, x0, l, m):
    y = -1.0*E*( np.cos(x-x0)*np.cos(x-x0) + 2.0*(lamb+2.0*mu)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0)) + (2.0*l+4.0*m)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0)) + 2.0*l*(S11*np.sin(x-x0)*np.sin(x-x0) + 2.0*S12*np.cos(x-x0)*np.cos(x-x0)) ) / (2.0*(lamb+2.0*mu))#-1.0*E* (np.cos(x-x0)*np.cos(x-x0) + 2.0*(lamb+2.0*mu)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0))  + (2*l+4*m)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0)) + 2*l*(S11*np.sin(x-x0)*np.sin(x-x0)+2*S12*np.cos(x-x0)*np.cos(x-x0)))/(lamb+2*mu)
    return y

def funcSV(x, x0, l, m, n):
    y = -1.0 * E * ( np.cos(x-x0)*np.cos(x-x0) + 2.0*mu*(S11*np.sin(x-x0)*np.sin(x-x0) + S12*np.cos(x-x0)*np.cos(x-x0)) + (l*0.5 + m)*(S11+S12) +0.5*l*(S12-S11)- (0.5 * n + l - m) * S12 ) / mu #- 1.0*E*( np.cos(x-x0)*np.cos(x-x0) + 2.0*(lamb+2.0*mu)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0)) + (2.0*l+4.0*m)*(S11*np.cos(x-x0)*np.cos(x-x0) + S12*np.sin(x-x0)*np.sin(x-x0)) + 2.0*l*(S11*np.sin(x-x0)*np.sin(x-x0) + 2.0*S12*np.cos(x-x0)*np.cos(x-x0)) ) / (2.0*(lamb+2.0*mu))
    
    return y


if __name__ == "__main__":
    params = sys.argv
    length = len(params)
    
    if(length != 4):
        quit()

    grs80 = pyproj.Geod(ellps='GRS80')

    inptfile = params[1]
    data = np.genfromtxt(inptfile,
                         dtype={'names':('theta','dat','std'),'formats':('f8','f8','f8')})
    
    inp_x  = data['theta']
    inp_y  = data['dat']/(2.*np.pi)
    inp_ey = data['std']/(2.*np.pi)

    inptfile_2 = params[2]
    data_2 = np.genfromtxt(inptfile_2,
                           dtype={'names':('name','lati_1','long_1','lati_2','long_2'),'formats':('S20','f8','f8','f8','f8')})
    
    stnname = data_2['name']
    lat1 = data_2['lati_1']
    lng1 = data_2['long_1']
    lat2 = data_2['lati_2']
    lng2 = data_2['long_2']

    OutFile = params[3]

    ##fitting initial parameter
    x0_init = np.radians(1.)
    l_init = -1E5#/(2.*np.pi)
    m_init = -1E5#/(2.*np.pi)
    n_init = -1E5#/(2.*np.pi)
    r_init = 0.5
    ##fitting
    inp_r = np.radians(inp_x)

    param_bounds = ((0,-1E6/(2.*np.pi),-1E6/(2.*np.pi)),(np.radians(360),1E6/(2.*np.pi),1E6/(2.*np.pi)))
    popt, pcov = curve_fit(funcP, inp_r, inp_y, p0=[x0_init, l_init, m_init])
    
    print (stnname, popt[0], '{:e}'.format(popt[1]), '{:e}'.format(popt[2]))#,'{:e}'.format(popt[3]))#, popt[4])    
