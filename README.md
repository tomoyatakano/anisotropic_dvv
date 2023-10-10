# anisotropic_dvv

This repository provides the data and scripts for the analyses in "Anisotropic seismic velocity variations in response to different orientations of tidal deformations 
". 

Citation: Takano, T., Nishimura, T., & Nakahara, H. (2023). Anisotropic seismic velocity variations in response to different orientations of tidal deformations. Geophysical Journal International, ggad386.

data/tide.pickle: 
pickle file of tidal strain at each azimuth for each volcano

data/result_table2.pickle:
pickle file of the parameters of fitted function

codes/fitting-dat.py:
script for estimating third-order elastic constants.

codes/azimuth-tide.c:
script for estimating tidal strain at each azimuth

codes/circular_corr.py:
script for computing correlation coeficient based on a circular statistics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8038742.svg)](https://doi.org/10.5281/zenodo.8038742)
