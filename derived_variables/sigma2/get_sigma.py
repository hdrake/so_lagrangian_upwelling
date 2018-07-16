#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc

def get_sigma2(file,lat,lon):

    ncFile = nc.Dataset(file)
    sigma2 = ncFile.variables['sigma2'][0,:,lat,lon]
    
    return sigma2