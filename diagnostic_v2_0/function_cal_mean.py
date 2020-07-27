#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:24:40 2017

@author: zhunguo, guozhun@uwm.edu, guozhun@lasg.iap.ac.cn 
"""

from netCDF4 import Dataset
import numpy as np
import scipy as sp
import os

from subprocess import call

def cal_mean(ncases, cases,years,nyear, nsite,lats, lons,area, filepath):
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
     infile=' '
     for yr in range (0, nyear[im] ):
         infile=infile+datalocal+cases[im]+'.cam.h0.'+str(years[im]+1+yr).rjust(4,'0')+'-01.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+1+yr).rjust(4,'0')+'-02.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-12.nc '
     os.system(ncea_str+infile +' -O '+outfile)
     
     outfile=filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-03.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-04.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-05.nc '
     os.system(ncea_str+infile +' -O '+outfile)
   
     outfile=filepath[im]+cases[im]+'/'+cases[im]+'_JJA_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-06.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-07.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-08.nc '

     os.system(ncea_str+infile +' -O '+outfile)
   
     outfile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc'
     infile=' '
     for yr in range (0, nyear[im]):
         infile=infile+datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-09.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-10.nc '+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+yr).rjust(4,'0')+'-11.nc '

     os.system(ncea_str+infile +' -O '+outfile)
   
     outfile=filepath[im]+cases[im]+'/'+cases[im]+'_ANN_climo.nc'
     infile=filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_JJA_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc '+filepath[im]+cases[im]+'/'+cases[im]+'_DJF_climo.nc '
     os.system(ncea_str+infile +' -O '+outfile)
   
     ln_str='ln -s '
     os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_SON_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_SON_budget_climo.nc' )
     os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_MAM_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_MAM_budget_climo.nc' )
     os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_JJA_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_JJA_budget_climo.nc' )
     os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_DJF_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_DJF_budget_climo.nc' )
     os.system(ln_str+ filepath[im]+cases[im]+'/'+cases[im]+'_ANN_climo.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_ANN_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-10.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_10_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-10.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_10_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-11.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_11_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-11.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_11_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-12.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_12_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-12.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_12_budget_climo.nc' )
   
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+1).rjust(4,'0')+'-01.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_01_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+1).rjust(4,'0')+'-01.nc'  +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_01_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+1).rjust(4,'0')+'-02.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_02_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]+1).rjust(4,'0')+'-02.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_02_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-03.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_03_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-03.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_03_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-04.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_04_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-04.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_04_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-05.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_05_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-05.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_05_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-06.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_06_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-06.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_06_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-07.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_07_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-07.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_07_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-08.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_08_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-08.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_08_budget_climo.nc' )
   
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-09.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_09_climo.nc' )
     os.system(ln_str+ datalocal+cases[im]+'.cam.h0.'+str(years[im]).rjust(4,'0')+'-09.nc' +' '+ filepath[im]+cases[im]+'/'+cases[im]+'_09_budget_climo.nc' )



