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
case="newzhun" # A general case name
outdir="/lcrc/group/acme/zhun/plots/" # Location of plots

filepath=["/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\
"/lcrc/group/acme/zhun/E3SM_simulations/",\

          ]
#cases=[ \
#       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try8.ne16_ne16",\   # n2p5 acc4 
##       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try6.ne16_ne16",\
#       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try14.ne16_ne16",\
#       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try15.ne16_ne16",\
#       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try3.ne16_ne16",\
##       "anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_test4.ne16_ne16",\
##       "anvil-centos7.clubb_silhs_v2_tau.ice_c50_rtclp10_bvm_facc1_kmax3o_berg2_4.ne16_ne16", \
#         "anvil-centos7.master_20191113.gust_polun_run3.ne16_ne16",\
#]

cases=[ \
      'anvil-centos7.new_zhun.LBAbest_level0_run3.ne16_ne16', \
      'anvil-centos7.new_zhun.LBAbest_newlatin_level0.ne16_ne16', \
      'anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_nogust.ne16_ne16', \
#      'anvil-centos7.clubb_silhs_v2_tau.facc1_kmax3o_berg2_dcs4h_acc_try12.ne16_ne16',\
      'anvil-centos7.master_20191113.gust_polun_run3.ne16_ne16',\
]

       
# Give a short name for your experiment which will appears on plots

#casenames=['arc4_prc3_1p2_wp3_wpxp1_n2p5','arc4_prc3_1p2_wp3_wpxp1_n2p45_c61p5','arc4_prc3_1p2_wp3_wpxp1_n2p5_c1c61p5','arc4_prc3','ZM']
casenames=['newlatin_repeat','bestLBA_newlatin','arc4_prc3_1p2_wp3_wpxp1_n2p4_nogust(bestLBA)','ZM']
years=[\
        "0001", "0001","0001", "0001","0001","0001"]
dpsc=[\
      "none","none","none","zm"]
# NOTE, dpsc,deep scheme, has to be "none", if silhs is turned on. 

# Observation Data
#filepathobs="/global/project/projectdirs/m2689/zhun/amwg/obs_data_20140804/"
filepathobs="/blues/gpfs/home/zhun/amwg_diag_20140804/obs_data_20140804/"
#------------------------------------------------------------------------
# Setting of plots.
ptype="png"   # eps, pdf, ps... are supported by this package
cseason="ANN" # Seasons, or others
casename=case+"_"+cseason

#------------------------------------------------------------------------
calmean=True      # make mean states
findout=True       # pick out the locations of your sites
draw2d=True        # This flag control 2D plots
drawlarge= True    # profiles for large-scale variable on your sites 
drawclubb= True    # profiles for standard clubb output
drawskw= True    # profiles for skewness functions
drawrain = True
drawbgt= True     # budgets of CLUBB prognostic Eqs 
drawe3smbgt= True #
drawmicrobgt= True
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
# sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
lats = [  20,  27, -20, -20,  -5,  -1,  60,   2,   9,   56,  45,   0,  10,  20,   0,   5,   9, -60]
lons = [ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311, 180, 295,  90, 205, 325, 280, 170, 340]

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
import draw_micro_budget
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

    pname = "std1"
    varis    = [ 'wp2','up2','vp2','rtp2','thlp2','wp3']
    cscale = [1, 1, 1, 1E6, 1, 1]
    plotstd1=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)

    pname = "std2"
    varis    = [ 'wprtp','wpthlp','wprcp','upwp','vpwp','wpthvp']
    cscale   = [     1E3,       1,    1E3,   1,     1,     1] 
    plotstd2=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)

    pname = "std3"
    varis    = [ 'wp2thlp','wp2thlp','wpthlp2','wprtp2','rcp2', 'wp2rcp']
    cscale   = [     1,       1,    1,   1E6,           1E6,    1E3] 
    plotstd3=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)

    varis    = [ 'wpthvp','wp2thvp','rtpthvp','thlpthvp','wp4','wprtpthlp']
    cscale   = [     1,           1,    1E3,       1,     1,     1E3] 
    pname = "std4"
    plotstd4=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,pname)


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
    varis=[ "AQRAIN","ANRAIN","ADRAIN","FREQR"]
    pname = "Rain"
    plotrain=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,pname)
    pname = "Snow"
    varis=[ "AQSNOW","ANSNOW","ADSNOW","FREQS"]
    plotsnow=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,pname)

if drawe3smbgt:
    print("e3sm_budget")
    plote3smbgt=draw_e3sm_budget.draw_e3sm_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,dpsc)

if drawmicrobgt:
    print("micro_budget")
    plotmicrobgt=draw_micro_budget.draw_micro_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,dpsc)

if drawbgt:
    print("CLUBB BUDGET on selected sites")
    varis   = [ 'wp2','wp3','up2','vp2']
    cscale  = [     1,    1,    1,    1]
    chscale = [    '1', '1',  '1',  '1']
    pname = "Budget1"
    plotbgt1=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)
    varis    = [ "wprtp","wpthlp",  "rtp2", "thlp2"]
    cscale   = [     1E7,     1E4,   1E11,     1E4]
    chscale  = [  '1E-7',  '1E-4', '1E-11',  '1E-4']

    pname = "Budget2"
    plotbgt2=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)   
    varis   = [  "um",   "vm",  "thlm",   "rtm"]
    cscale  = [   1E4,    1E4,     1E5,     1E8]
    chscale = ['1E-4', '1E-4',  '1E-5',  '1E-8']

    pname = "Budget3"
    plotbgt3=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)  

if makeweb:
    print("Making webpages")
    Diagnostic_webpage.main_web(casename,casedir)

    Diagnostic_webpage.sets_web(casename,casedir,"diff.*.asc","txt",\
                                "Gitdiff","1000","1000")

    if (draw2d):
        plot2d.extend(plot3d[:])
        Diagnostic_webpage.sets_web(casename,casedir,plot2d,"2D",\
				"Horizontal Plots","1000","1000")

    for ire in range (0, nsite):
        plotclb=[]
        if (drawlarge):
           plotclb.append(plotlgs[ire])
        if (drawclubb):  
           plotclb.append(plotstd1[ire])
           plotclb.append(plotstd2[ire])
           plotclb.append(plotstd3[ire])
           plotclb.append(plotstd4[ire])
        if (drawskw):
           plotclb.append(plotskw[ire])
        if (drawhf):
           plotclb.append(plothf[ire])
        if (drawrain):
           plotclb.append(plotrain[ire])
           plotclb.append(plotsnow[ire])

        if (drawmicrobgt):
           for im in range (0, ncases ):
               plotclb.append(plotmicrobgt[ire*ncases+im])

        if (drawe3smbgt):
           for im in range (0, ncases ):
               plotclb.append(plote3smbgt[ire*ncases+im])

        if (drawbgt):
           for im in range (0, ncases ):
               plotclb.append(plotbgt1[ire*ncases+im])
           for im in range (0, ncases ):
               plotclb.append(plotbgt2[ire*ncases+im])
           for im in range (0, ncases ):
               plotclb.append(plotbgt3[ire*ncases+im])

        Diagnostic_webpage.sets_web(casename,casedir,plotclb,str(lons[ire])+'E_'+str(lats[ire])+'N',\
                                  "Profiles on "+str(lons[ire])+'E_'+str(lats[ire])+'N',"908","636")

if maketar:
    print("Making tar file of case")
    Common_functions.make_tarfile(outdir+casename+'.tar',outdir+casename)
    
