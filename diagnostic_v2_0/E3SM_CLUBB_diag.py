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
casename="c7rvrs_n2p4_c81_clpf_fac1p5_a14_fnr" # A general case name
outdir="/lcrc/group/acme/zhun/plots/" # Location of plots

filepath=["/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
          ]
cases=[ \
        "anvil-centos7.clubb_silhs_v2_tau.c7rvrs_n2p4_c81_clpf_fac1p5_dsc2k_a14_fnr.ne16_ne16",\
        "anvil-centos7.clubb_silhs_v2_tau.c7rvrs_n2p4_c81_clpf_fac1p5_dsc2k.ne16_ne16", \
        "anvil-centos7.master_20191113.gust_polun_run2.ne16_ne16" \
              ]
# Give a short name for your experiment which will appears on plots

casenames=['fnr_min','def','zm_gust']

years=[\
         "0001",  "0001","0001"]
dpsc=[\
      "none","none","zm"]
# NOTE, dpsc has to be "none", if silhs was turned on. 

# Observation Data
#filepathobs="/global/project/projectdirs/m2689/zhun/amwg/obs_data_20140804/"
filepathobs="/blues/gpfs/home/zhun/amwg_diag_20140804/obs_data_20140804/"
#------------------------------------------------------------------------
# Setting of plots.
ptype="png"   # eps, pdf, ps... are supported by this package
cseason="ANN" # Seasons, or others

#------------------------------------------------------------------------
calmean=False      # make mean states
findout=True       # pick out the locations of your sites
draw2d=True        # This flag control 2D plots
drawlarge= True    # profiles for large-scale variable on your sites 
drawclubb= True    # profiles for standard clubb output
drawskw= True    # profiles for skewness functions
drawrain = True
drawbgt= True     # budgets of CLUBB prognostic Eqs 
drawe3smbgt= True #
# ONLY for SILHS
drawhf= False     # Tendency of holl filler 
drawsilhs=False    # profiles for silhs variables

makeweb=True      # Make a webpage?
maketar=True      # Tar them?

clevel = 500
area  = 1.
# Note, this para helps to find out the 'ncol' within
# lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area
#------------------------------------------------------------------------
# Please give the lat and lon of sites here.
# sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
lats = [  20,  27, -20, -20,  -5,  -1,  60,   2,   9,   56,  45,   0,  10,  20,   0,   5]
lons = [ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311, 180, 295,  90, 205, 325, 280]

#========================================================================

#------------------------------------------------------------------------
# Do not need to change
#------------------------------------------------------------------------

ncases=len(cases)
nsite=len(lats)

casedir=outdir+casename
print(casedir)

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
import draw_hollfiller
import draw_rain 
import draw_e3sm_budget
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
    function_pick_out.pick_out(ncases, cases, years, nsite, lats, lons, area, filepath,casedir)

if draw2d:
    print("Drawing 2d")
    plot2d=draw_plots_hoz_2D.draw_2D_plot(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)
    plot3d=draw_plots_hoz_3D.draw_3D_plot(ptype,clevel,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawlarge:
    print("Large-scale variables on selected sites")
    plotlgs=draw_large_scale.large_scale_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawclubb:
    print("CLUBB standard variables on selected sites")
    plotstd=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawskw:
    print("CLUBB standard variables on selected sites")
    plotskw=draw_clubb_skew.clubb_skw_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawsilhs:
    print("CLUBB standard variables on selected sites")
    plotsilhs=draw_silhs_standard.silhs_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawhf:
    print("Holl filler")
    plothf=draw_hollfiller.hollfiller_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawrain:
    print("Energy")
    plotrain=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawe3smbgt:
    print("e3sm_budget")
    plote3smbgt=draw_e3sm_budget.draw_e3sm_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,dpsc)

if drawbgt:
    print("CLUBB BUDGET on selected sites")
    varis = [ "wp2","wp3","up2","vp2"]
    cscale = [1, 1, 1,1]
    pname = "Budget1"
    plotbgt1=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)
    varis = [ "wprtp","wpthlp","rtp2","thlp2"]
    cscale = [ 1E3, 1, 1E6, 1]
    pname = "Budget2"
    plotbgt2=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)   
    varis = [ "um","vm","thlm","rtm"]
    cscale = [ 1, 1, 1, 1E3]
    pname = "Budget3"
    plotbgt3=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)  
    print(plotbgt3)

if makeweb:
    print("Making webpages")
    Diagnostic_webpage.main_web(casename,casedir)

    Diagnostic_webpage.sets_web(casename,casedir,"diff.*.asc","txt",\
                                "Gitdiff","1000","1000")

    if (draw2d):
        plot2d.extend(plot3d[:])
        Diagnostic_webpage.sets_web(casename,casedir,plot2d,"2D",\
				"Horizontal Plots","1000","1000")

    if (drawlarge or drawclubb or drawskw or drawsilhs or drawbgt or drawrain or drawhf or drawe3sm_budget):
         for ire in range (0, nsite):
             plotclb=[]
             plotclb.append(plotlgs[ire])
             plotclb.append(plotstd[ire])
             plotclb.append(plotskw[ire])
#             plotclb.append(plothf[ire])
             plotclb.append(plotrain[ire])
             plotclb.append(plote3smbgt[ire*ncases])
             plotclb.append(plote3smbgt[ire*ncases+1])
             plotclb.append(plote3smbgt[ire*ncases+2])
             plotclb.append(plotbgt1[ire*ncases])
             plotclb.append(plotbgt1[ire*ncases+1])
             plotclb.append(plotbgt1[ire*ncases+2])
             plotclb.append(plotbgt2[ire*ncases])        
             plotclb.append(plotbgt2[ire*ncases+1])
             plotclb.append(plotbgt2[ire*ncases+2])
             plotclb.append(plotbgt3[ire*ncases])         
             plotclb.append(plotbgt3[ire*ncases+1])
             plotclb.append(plotbgt3[ire*ncases+2])
             Diagnostic_webpage.sets_web(casename,casedir,plotclb,str(lons[ire])+'E_'+str(lats[ire])+'N',\
                                  "Profiles on "+str(lons[ire])+'E_'+str(lats[ire])+'N',"908","636")

if maketar:
    print("Making tar file of case")
    Common_functions.make_tarfile(outdir+casename+'.tar',outdir+casename)
    
