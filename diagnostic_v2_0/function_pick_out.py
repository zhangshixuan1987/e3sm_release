#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:24:40 2017

@author: zhunguo, guozhun@uwm.edu, guozhun@lasg.iap.ac.cn 
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import os
from subprocess import call


def pick_out(ncases, cases,years, nsite,lats, lons,area, filepath,casedir):
# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data
 print(ncases)
# inptrs = [ncases]

 if not os.path.exists("data"):
        os.mkdir("data")


 for im in range(0,ncases):
    
     infile=filepath[im]+cases[im]+'/run/'+cases[im]+'.cam.h0.'+years[im]+'-01.nc'

     print(infile)
     print(im)
     inptrs = Dataset(infile,'r')       # pointer to file1
     lat=inptrs.variables['lat'][:]
     nlat=len(lat)
     lon=inptrs.variables['lon'][:] 
     nlon=len(lon)
#     idx_cols=[[0 for i in range(5)] for j in range(nsite)] 
     nh=[0 for k in range(nsite)]
#     print(idx_cols)
     cols=[0,1,2,3,4]
     sits=np.linspace(0,nsite-1,nsite)

     txtfile=filepath[im]+cases[im]+'/run/diff*.asc'
     os.system('mkdir '+ casedir+'/txt/')
     os.system('cp -f '+ txtfile+ ' '+ casedir+'/txt/')

     os.system('rm -f ./data/'+cases[im]+'_site_location.nc')
     outf =Dataset('./data/'+cases[im]+'_site_location.nc','w')
     outf.createDimension("sit",nsite)
     outf.createDimension("col",5)
#     outf.variables['sit'][:]=sits
#     outf.variables['col'][:]=cols
     outf.createVariable('idx_cols','i',('sit','col'))
     outf.createVariable('n','i',('sit'))
     outf.variables['n'][:]=0
     outf.variables['idx_cols'][:,:]=0

# ========================================================================== 
# find out the cols and their numbers
#    the axis of site is stored in idx_cols(site,n)

     for i in range(0,nlat):
         for j in range(0,nsite): 
             if (lon[i] >= lons[j]-area) & (lon[i] < lons[j]+area) & (lat[i]>=lats[j]-area) & (lat[i] < lats[j]+area): 
                 outf.variables['idx_cols'][j,nh[j]]=i 
                 outf.variables['n'][j]=nh[j]+1
     outf.close()
