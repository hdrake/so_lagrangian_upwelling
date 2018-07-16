#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import numpy as np
import os
from jdcal import gcal2jd

dir = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/sigma0/'
dir = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/for_CMS_SO_Edited/'

dir = '/archive/hfd/CM2.5/CM2.5_A_Control-1860_Y03-MBLING-BLING/5day/sigma2/'

dir = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/for_CMS_SO_Edited/'

files = os.listdir(dir)
files = np.sort(files)

count = 0
for file in files:
    filenum = file.split('.')[1]
    #filenum = file.split('_')[2]
    if count > 0:
        if np.sum(gcal2jd(int(filenum[0:4]),int(filenum[4:6]),int(filenum[6:8])))-np.sum(gcal2jd(int(lastnum[0:4]),int(lastnum[4:6]),int(lastnum[6:8])))!= 5:
            print filenum
    else:
        print dir+file
    count+=1
    lastnum = filenum