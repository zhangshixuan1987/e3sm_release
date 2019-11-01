This is a diagnostic package for E3SM_CLUBB_SILHS V2.1



Environments:  
Libraries needed, netcdf4, scipy, pyngl, pynio, ncl and matplotlib.

1 Install libs and Run diagnostic package:
 ./conda create -n ncl_to_python -c conda-forge xarray netcdf4 scipy pyngl pynio ncl matplotlib

 source activate ncl_to_python  # ncl_to_python could be replaced by other name.

 Then, 
 ~/.conda/envs/ncl_to_python/bin/python E3SM_CLUBB_diag.py   # RUN


2 Set some parameters:
 
 2.1 User defined name used for this comparison, this will be the name given the directory for these diagnostics
        casename="name"
        outdir="/home/zhun/E3SM_code/clubb_silhs_v2_tau/diagnostic_v2_0/"

2.2 Directory and casenames of output files live
 NOTE, make sure you have enough space in these Directories, because the climotoloical file will be generated in these directories
        filepath=["/lcrc/group/acme/zhun/E3SM_simulations/",\
              ]
 Real Case Names:
        cases=[ \
              "anvil.clubb_silhs_v2_tau.1year_test1.ne16_ne16",\
              ]

 Tell the package when you plan to draw:
        years=[\
              "0001","0001"]
 Give a short name for each cases, which will appear on the plots.
        casenames=['mergetau','mergesilhs']

2.3 Setting for Package.
        ptype="png"   # eps, pdf, ps... are supported by this package
        cseason="JJA" # Seasons
        clevel = 500
        area  = 1.
        Note! This para helps to find out the 'ncol' within lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area. If error message reads as "UnboundLocalError: local variable 'theunits' referenced before assignment", it usually means a large area is needed.

2.4 Flags 
    calmean=False     # make mean states, it should be 'True' when you draw plot first.
    findout=True      # pick out the locations of your sites, it should be 'True' when you draw plot first.
    draw2d=True       # Need 2D plots?
    drawlarge= True   # Need profiles for large-scale variable on your sites 
    drawclubb=True    # Need profiles for standard clubb output
    drawskw=True      # Need profiles for skewness functions
    drawsilhs=False   # Need profiles for silhs variables
    drawbgt=True      # Need budgets of 9 prognostic Eqs 
    makeweb=True      # Make a webpage?
    maketar=True      # Tar them? 

