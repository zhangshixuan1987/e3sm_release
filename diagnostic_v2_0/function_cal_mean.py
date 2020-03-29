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

def cal_mean(ncases, cases,years, nsite,lats, lons,area, filepath):
# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data

 for im in range (0, ncases ):
#  call('rm -f  ',filepath[im]+cases[im]+'/'+cases[im]+'_*_climo.nc')
  datalocal = filepath[im] +cases[im]+'/run/'
  print(datalocal)
#  ncea_str='/global/common/sw/cray/cnl7/haswell/nco/4.7.9/gcc/8.2.0/unbt25h/bin/ncea '
  ncea_str='/blues/gpfs/home/software/spack-0.10.1/opt/spack/linux-centos7-x86_64/intel-17.0.4/nco-4.7.4-x4y66ep2ydoyegnckicvv5ljwrheniun/bin/ncea '

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_DJF_climo.nc'
  infile=datalocal+cases[im]+'.cam.h0.0002-01.nc '+ datalocal+cases[im]+'.cam.h0.0002-02.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-12.nc '
  os.system(ncea_str+infile +' -O '+outfile)

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_DJF_clm.nc'
  infile=datalocal+cases[im]+'.clm2.h0.0002-01.nc '+ datalocal+cases[im]+'.clm2.h0.0002-02.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-12.nc '
  os.system(ncea_str+infile +' -O '+outfile)


  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc'
  infile=datalocal+cases[im]+'.cam.h0.'+years[im]+'-03.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-04.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-05.nc '
  os.system(ncea_str+infile +' -O '+outfile)

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_MAM_clm.nc'
  infile=datalocal+cases[im]+'.clm2.h0.'+years[im]+'-03.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-04.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-05.nc '
  os.system(ncea_str+infile +' -O '+outfile)


  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_JJA_climo.nc'
  infile=datalocal+cases[im]+'.cam.h0.'+years[im]+'-06.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-07.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-08.nc '
  os.system(ncea_str+infile +' -O '+outfile)

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_JJA_clm.nc'
  infile=datalocal+cases[im]+'.clm2.h0.'+years[im]+'-06.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-07.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-08.nc '
  os.system(ncea_str+infile +' -O '+outfile)


  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc'
  infile=datalocal+cases[im]+'.cam.h0.'+years[im]+'-09.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-10.nc '+ datalocal+cases[im]+'.cam.h0.'+years[im]+'-11.nc '
  os.system(ncea_str+infile +' -O '+outfile)

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_clm.nc'
  infile=datalocal+cases[im]+'.clm2.h0.'+years[im]+'-09.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-10.nc '+ datalocal+cases[im]+'.clm2.h0.'+years[im]+'-11.nc '
  os.system(ncea_str+infile +' -O '+outfile)


  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_ANN_climo.nc'
  infile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_JJA_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_DJF_climo.nc '
  os.system(ncea_str+infile +' -O '+outfile)

  outfile=filepath[im]+cases[im]+'/'+cases[im]+'_ANN_clm.nc'
  infile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_clm.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_JJA_clm.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_MAM_clm.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_DJF_clm.nc '
  os.system(ncea_str+infile +' -O '+outfile)



  ln_str='ln -s '
  os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_SON_budget_climo.nc' )
  os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_MAM_budget_climo.nc' )
  os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_JJA_budget_climo.nc' )
  os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_DJF_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_DJF_budget_climo.nc' )
  os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_ANN_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_ANN_budget_climo.nc' )



