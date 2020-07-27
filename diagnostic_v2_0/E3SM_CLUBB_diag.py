# -*- coding: utf-8 -*-
'''
E3SM CLUBB Diagnostics package 

Main code to make 1) 2D plots,2) profiles, 3) budgets on selected stations, 
         and then build up  webpages  etc
    zhunguo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
'''

## ==========================================================
# Begin User Defined Settings
# User defined name used for this comparison, this will be the name 
#   given the directory for these diagnostics
case='paper4' # A general case name
outdir='/lcrc/group/acme/zhun/plots/' # Location of plots

filepath=[ \
'/lcrc/group/acme/ac.zguo/E3SM_simulations/',\
'/lcrc/group/acme/ac.zguo/E3SM_simulations/',\
'/lcrc/group/acme/ac.zguo/E3SM_simulations/',\
'/lcrc/group/acme/ac.zguo/E3SM_simulations/',\
'/lcrc/group/acme/ac.zguo/E3SM_simulations/',\
'/lcrc/group/acme/zhun/E3SM_simulations/',\
          ]
cases=[ \
'anvil.EAMv1.F2010SC5-CMIP6_t1.ne30_ne30',\
'anvil-centos7.base2.wpxpri_3p3e4_1_3_0_12_C7ri.ne30_ne30',\
'anvil.EAMv1.FC5AV1C.ne30_ne30',\
'anvil-centos7.base2.wpxpri_3p3e4_1_3_0_12_C7ri_FC5AV1C.ne30_ne30',\
'anvil-centos7.base2.wpxpri_3p35e4_1_2p5_0_12_C7ri_F2010SC5-CMIP6.ne30_ne30',\
]

       
# Give a short name for your experiment which will appears on plots

casenames=[
'EAMv1.F2010SC5-CMIP6',\
'wpxpri_3p3e4_1_3_0_12_C7ri',\
'EAMv1.FC5AV1C',\
'wpxpri_3p3e4_1_3_0_12_C7ri_FC5AV1C',\
'wpxpri_3p35e4_1_2p5_0_12_C7ri_F2010SC5-CMIP6',\
]

years=[\
        1, 1, 1, 1,1,1]
nyear=[\
        10, 5, 1, 5,3,1]

dpsc=[\
      'none','none','none','none','none','none']
# NOTE, dpsc,deep scheme, has to be 'none', if silhs is turned on. 

# Observation Data
#filepathobs='/global/project/projectdirs/m2689/zhun/amwg/obs_data_20140804/'
filepathobs='/blues/gpfs/home/ac.zguo/amwg_diag_20140804/obs_data_20140804/'
#------------------------------------------------------------------------
# Setting of plots.
ptype         ='png'   # eps, pdf, ps, png, x11, ... are supported by this package
cseason       ='ANN' # Seasons, or others
casename      =case+'_'+cseason

#------------------------------------------------------------------------
calmean          = True       # make mean states
findout          = True       # pick out the locations of your sites
draw2d           = True       # 2D plots, SWCF etc.
drawlarge        = True       # profiles for large-scale variable on your sites 
drawclubb        = True       # profiles for standard clubb output
drawskw          = False       # profiles for skewness functions
drawrain         = True       # profiles for SNOW, Rain etc.
drawbgt          = True       # budgets of CLUBB prognostic Eqs 
drawe3smbgt      = True       # budgets of e3sm tendency
drawmicrobgt     = False       # budgets of MG2
drawaero         = False       # AERO for cloud brone
# ONLY for SILHS
drawhf           = False      # Tendency of holl filler 
drawsilhs        = False      # profiles for silhs variables

makeweb          =True        # Make a webpage?
maketar          =True        # Tar them?

clevel = 500
area  = 1.
# Note, this para helps to find out the 'ncol' within
# lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area
#------------------------------------------------------------------------
# Please give the lat and lon of sites here.
# sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   23   24   25   26   27   28   29   30  31  32 33
lats = [  20,  27, -20, -20,  -5,  -1,  60,   2,   9,   56,  45,   0,  10,  20,   0,   5,   9, -60,   0,   0, -45, -75,  30,  25 , 70 , 15,  17,  13,  36, 28, 29, 30, 27 , 30, 27]
lons = [ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311, 180, 295,  90, 205, 325, 280, 170, 340, 305,  25,  90,  90,  90, 105 , 90, 300, 300, 300, 263, 240,240, 240 , 242, 238, 238 ]


#========================================================================

#------------------------------------------------------------------------
# Do not need to change
#------------------------------------------------------------------------

ncases =len(cases)
nsite  =len(lats)

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
    print('Getting climatological mean')
    function_cal_mean.cal_mean(ncases, cases, years,nyear, nsite, lats, lons, area, filepath)

if findout:
    print('Find out the sites')
    function_pick_out.pick_out(ncases, cases, years, nsite, lats, lons, area, filepath,casedir)

if draw2d:
    print('Drawing 2d')
    plot2d=draw_plots_hoz_2D.draw_2D_plot(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)
    clevel=500
    plot3d=draw_plots_hoz_3D.draw_3D_plot(ptype,clevel,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawlarge:
    print('Drawing Large-scale variables on selected sites')
    plotlgs=draw_large_scale.large_scale_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawclubb:
    print('Drawing CLUBB standard variables on selected sites')

    pname = 'std1'
    varis    = [ 'wp2','up2','vp2','rtp2','thlp2','wp3']
    cscale   = [     1,    1,    1,   1E6,      1,    1]
    chscale  = [   '1',  '1',  '1','1E-6',    '1',  '1']
    plotstd1=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'std2'
    varis    = [ 'wprtp','wpthlp','wprcp','upwp','vpwp','rtpthlp']
    cscale   = [     1E3,       1,    1E3,     1,     1,     1E3] 
    chscale  = [  '1E-3',     '1', '1E-3',   '1',    '1',    '1E-3']

    plotstd2=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'std3'
    varis    = [ 'wp2thlp','wp2rtp','wpthlp2','wprtp2','rcp2', 'wp2rcp']
    cscale   = [         1,        1,       1,     1E6,   1E6,      1E3] 
    chscale  = [       '1',      '1',     '1',  '1E-6','1E-6',   '1E-3']

    plotstd3=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    varis    = [ 'wpthvp','wp2thvp','rtpthvp','thlpthvp','wp4','wprtpthlp']
    cscale   = [        1,        1,      1E3,         1,    1,        1E3] 
    chscale  = [      '1',      '1',   '1E-3',       '1',  '1',    '1-E-3']
    pname = 'std4'
    plotstd4=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

#    pname = 'Tau'
#    varis   = [ 'invrs_tau_bkgnd','invrs_tau_shear','invrs_tau_sfc','tau_no_N2_zm','tau_zm','tau_wp2_zm','tau_xp2_zm','tau_wp3_zm',  'bv_freq_sqd']
#    cscale  = [               1E3,              1E3,            1E3,           1E3,     1E3,         1E3,         1E3,         1E3,            1E3]
#    chscale = [            '1E-3',           '1E-3',         '1E-3',        '1E-3',  '1E-3',      '1E-3',      '1E-3',      '1E-3',         '1E-3']
#    plottau=draw_clubb_standard.clubb_std_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

if drawskw:
    print('Drawing CLUBB skewness functions on selected sites')
    plotskw=draw_clubb_skew.clubb_skw_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawsilhs:
    print('CLUBB standard variables on selected sites')
    plotsilhs=draw_silhs_standard.silhs_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawhf:
    print('Holl filler')
    plothf=draw_hollfiller.hollfiller_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir)

if drawrain:
    print('Drawing Rain and Snow properties')

    pname = 'Rain'
    varis   = [ 'AQRAIN','ANRAIN','ADRAIN','FREQR']
    cscale  = [      1E6,     1,       1E4,      1]
    chscale = [   '1E-6',  '1',    '1E-4',     '1']
    plotrain=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'Snow'
    varis   = [ 'AQSNOW','ANSNOW','ADSNOW','FREQS']
    cscale  = [      1E6,     1,       1E4,      1]
    chscale = [   '1E-6',  '1',    '1E-4',     '1']
    plotsnow=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'NUM'
    varis   = [ 'AWNC','AWNI','AREL','AREI']
    cscale  = [     1E-7,    1E-3,        1,         1]
    chscale = [    '1E7',   '1E3',      '1',       '1']
    plotsnum=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'RAINQM'
    varis   = [ 'RAINQM','NUMRAI']
    cscale  = [     1E-12,    1E-3]
    chscale = [    '1E12',   '1E3']
    plotsqm=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)


if drawaero:
    print('Drawing  aerosol related vars')

    pname = 'FREZ'
    varis   = [ 'DSTFREZIMM','BCFREZIMM','DSTFREZCNT','BCFREZCNT','DSTFREZDEP','BCFREZDEP']
    cscale  = [            1,          1,        1E12,       1E3 ,         1E8,        1E3]
    chscale = [          '1',        '1',     '1E-12',     '1E-3',      '1E-8',     '1E-3']
    plotsaero1=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'FREQAERO'
    varis   = ['FREQIMM','FREQCNT','FREQDEP','FREQMIX']
    cscale  = [      1,     1,       1E3,      1]
    chscale = [    '1',   '1',    '1E-3',     '1']
    plotsaero2=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    pname = 'ACN'
    varis   = [ 'bc_a1_num','bc_c1_num','dst_a1_num','dst_a3_num','dst_c1_num','dst_c3_num']
    cscale  = [            1,         1,           1,         1E3,          1,        1E3]
    chscale = [          '1',       '1',         '1',      '1E-3',         '1',     '1E-3']
    plotsaero3=draw_rain.rain_prf(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)


if drawe3smbgt:
    print('Drawing e3sm standard budgets')
    plote3smbgt=draw_e3sm_budget.draw_e3sm_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,dpsc)

if drawmicrobgt:
    print('Drawing MG budget')
    pname = 'mirco1'
    vname   = [ 'Liquid','Ice','Rain','Snow']
    varis   = [ 'MPDLIQ','MPDICE','QRSEDTEN','QSSEDTEN'] # We just need a unit
    cscale  = [      1E9,     1E9,       1E9,       1E9]
    chscale = [   '1E-9',  '1E-9',    '1E-9',    '1E-9']
    plotmicrobgt1=draw_micro_budget.draw_micro_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,vname,cscale,chscale,pname)

    pname = 'micro2'
    vname   = [ 'Vapor','NUMCLDLIQ','NUMCLDICE']
    varis   = [ 'QISEVAP','nnuccco',  'nnuccdo']  # We just need a unit
    cscale  = [      1E9,     1E-3,      1E-1]
    chscale = [   '1E-9',    '1E3',     '1E1']
    plotmicrobgt2=draw_micro_budget.draw_micro_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,vname,cscale,chscale,pname)

if drawbgt:
    print('Drawing CLUBB BUDGET')
    varis   = [ 'wp2','wp3','up2','vp2']
    cscale  = [     1,    1,    1,    1]
    chscale = [   '1',  '1',  '1',  '1']
    pname = 'Budget1'
    plotbgt1=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)

    varis    = [ 'wprtp','wpthlp',  'rtp2', 'thlp2']
    cscale   = [     1E7,     1E4,    1E11,     1E4]
    chscale  = [  '1E-7',  '1E-4', '1E-11',  '1E-4']
    pname = 'Budget2'
    plotbgt2=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)   

    varis   = [  'um',   'rtpthlp',  'thlm',   'rtm']
    cscale  = [   1E4,    1E4,     1E5,     1E8]
    chscale = ['1E-4', '1E-4',  '1E-5',  '1E-8']
    pname = 'Budget3'
    plotbgt3=draw_clubb_budget.draw_clubb_bgt(ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,cscale,chscale,pname)  

if makeweb:
    print('Making webpages')
    Diagnostic_webpage.main_web(casename,casedir)

    Diagnostic_webpage.sets_web(casename,casedir,'diff.*.asc','txt',\
                                'Gitdiff','1000','1000')

    if (draw2d):
        plot2d.extend(plot3d[:])

        Diagnostic_webpage.sets_web(casename,casedir,plot2d,'2D',\
				'Horizontal Plots','1000','1000')



    for ire in range (0, nsite):
        plotclb=[]
        if (drawlarge):
           plotclb.append(plotlgs[ire])
        if (drawclubb):  
           plotclb.append(plotstd1[ire])
           plotclb.append(plotstd2[ire])
           plotclb.append(plotstd3[ire])
           plotclb.append(plotstd4[ire])
#           plotclb.append(plottau[ire])

        if (drawskw):
           plotclb.append(plotskw[ire])
        if (drawhf):
           plotclb.append(plothf[ire])
        if (drawrain):
           plotclb.append(plotrain[ire])
           plotclb.append(plotsnow[ire])
           plotclb.append(plotsnum[ire])
           plotclb.append(plotsqm[ire])

        if (drawaero):
           plotclb.append(plotsaero1[ire])
           plotclb.append(plotsaero2[ire])
           plotclb.append(plotsaero3[ire])


        if (drawmicrobgt):
           for im in range (0, ncases ):
               plotclb.append(plotmicrobgt1[ire*ncases+im])
#               plotclb.append(plotmicrobgt2[ire*ncases+im])

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
                                  'Profiles on '+str(lons[ire])+'E_'+str(lats[ire])+'N','908','636')

if maketar:
    print('Making tar file of case')
    Common_functions.make_tarfile(outdir+casename+'.tar',outdir+casename)
    
