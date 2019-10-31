# -*- coding: utf-8 -*-
"""
E3SM CLUBB Diagnostics package 

Main code to make 1) 2D plots,2) profiles, 3) budgets on selected stations, 
         and then build up  webpages  etc
    zhunguo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
"""

## ==========================================================
# Begin User Defined Settings
# User defined name used for this comparison, this will be the name 
#   given the directory for these diagnostics
casename="pic"
outdir="/global/u2/z/zhun/E3SM.wiki/python_scripts/"

filepath=["/global/cscratch1/sd/zhun/E3SM_simulations/",\
               "/global/cscratch1/sd/zhun/E3SM_simulations/" \
              ]
cases=[ \
               "zhun.A_WCYCL1850S_CMIP6_ifort18_taus_splat1.ne30_oEC.cori-knl",\
               "zhun.A_WCYCL1850S_CMIP6_ifort18_taus_C63.ne30_oEC.cori-knl"    \
              ]
# Give a short name for your experiment which will appears on plots

casenames=['splat=1','C63']

years=[\
                "0001","0001"]

# Observation Data
filepathobs="/global/project/projectdirs/m2689/zhun/amwg/obs_data_20140804/"

#------------------------------------------------------------------------
# Setting of plots.
ptype="png"   # eps, pdf, ps... are supported by this package
cseason="JJA" # Seasons, or others

#------------------------------------------------------------------------
# set to true when you first use it.
calmean=False     # make mean states
findout=False     # pick out the locations of your sites
draw2d=True       # This flag control 2D plots
drawlarge= True   # profiles for large-scale variable on your sites 
drawclubb=True    # profiles for standard clubb output
drawskw=True      # profiles for skewness functions
drawsilhs=False   # profiles for silhs variables
drawbgt=True      # budgets of 9 prognostic Eqs 

makeweb=True
maketar=True

#------------------------------------------------------------------------
# Please give the lat and lon of sites here.
# sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
#lats = [  20,  27, -20, -20,  -5,  -1, -60,  60,   2,   9,  56,  76,  45,   0,  10]
#lons = [ 190, 240, 275, 285, 355, 259, 340, 180, 140, 229, 311, 320, 180, 295,  90]
lats = [  20,  27, -20, -20,  -5,  -1,  60,   2,   9,  56,  45,   0,  10]
lons = [ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311,  180, 295,  90]

clevel = 500
area  = .5
# Note, this para helps to find out the 'ncol' within
# lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area

#========================================================================
#------------------------------------------------------------------------
# Do not need to change
ncases=len(cases)
nsite=len(lats)

import os
import function_cal_mean
import function_pick_out
import draw_plots_hoz_2D
import draw_plots_hoz_3D
import draw_large_scale
import draw_clubb_skew
import draw_silhs_standard 
import draw_clubb_standard
import draw_clubb_budget
import Common_functions
import Diagnostic_webpage


casedir=outdir+casename

if not os.path.exists(casedir):
    os.mkdir(casedir)

if calmean:
    print("Getting climatological mean")
    function_cal_mean.cal_mean(ncases, cases, years,nsite, lats, lons, area, filepath)

if findout:
    print("Find out the sites")
    function_pick_out.pick_out(ncases, cases, years, nsite, lats, lons, area, filepath)

if draw2d:
    print("Drawing 2d")
    plot2d=draw_plots_hoz_2D.draw_2D_plot(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)
    plot3d=draw_plots_hoz_3D.draw_3D_plot(ptype,clevel,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

if drawlarge:
    print("Large-scale variables on selected sites")
    plotlgs=draw_large_scale.large_scale_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

if drawclubb:
    print("CLUBB standard variables on selected sites")
    plotstd=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

if drawskw:
    print("CLUBB standard variables on selected sites")
    plotskw=draw_clubb_skew.clubb_skw_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

if drawsilhs:
    print("CLUBB standard variables on selected sites")
    plotsilhs=draw_silhs_standard.silhs_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

if drawbgt:
    print("CLUBB BUDGET on selected sites")
    plotbgt=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)

## ============================
# Make Web Page   

if makeweb:
    print("Making webpages")
    Diagnostic_webpage.main_web(casename,casedir)
    if (draw2d):
        plot2d.extend(plot3d[:])
        Diagnostic_webpage.sets_web(casename,casedir,plot2d,"2D",\
                                  "Horizontal Plots","354","268")

    if (drawlarge or drawclubb or drawskw or drawsilhs or drawbgt):
         for ire in range (0, nsite):
             plotclb=[]
             plotclb.append(plotlgs[ire])
             plotclb.append(plotstd[ire])
             plotclb.append(plotskw[ire])
             plotclb.append(plotbgt[ire*ncases])
             plotclb.append(plotbgt[ire*ncases+1])        
             Diagnostic_webpage.sets_web(casename,casedir,plotclb,str(lons[ire])+'E_'+str(lats[ire])+'N',\
                                  "Profiles on "+str(lons[ire])+'E_'+str(lats[ire])+'N',"454","318")

if maketar:
    print("Making tar file of case")
    Common_functions.make_tarfile(outdir+casename+'.tar',outdir+casename)
    
