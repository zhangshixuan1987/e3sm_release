'''
    RAIN parameters 
    zhunguo : guozhun@lasg.iap.ac.cn ; guozhun@uwm.edu
'''

import Ngl
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import os
from subprocess import call


def rain_prf (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir):

# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data
# inptrs = [ncases]
 if not os.path.exists(casedir):
        os.mkdir(casedir)

 _Font   = 25
 interp = 2
 extrap = False
 mkres = Ngl.Resources()
 mkres.gsMarkerIndex = 2
 mkres.gsMarkerColor = "Red"
 mkres.gsMarkerSizeF = 15.
 infiles  = ["" for x in range(ncases)]
 ncdfs    = ["" for x in range(ncases)]
 nregions = nsite

 varis    = [ "AQRAIN","ANRAIN","ADRAIN","FREQR"]
 varisobs = ["CC_ISBL", "OMEGA","SHUM","CLWC_ISBL", "THETA","RELHUM","U","CIWC_ISBL","T" ]
 nvaris = len(varis)
 cunits = ["%","mba/day","g/kg","g/kg","K", "%", "m/s", "g/kg", "m/s", "m/s","K","m" ]
 cscale = [1, 1, 1, 1, 1., 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
 cscaleobs = [100,        1,     1, 1000 , 1.,   1,     1,   1000,     1,1,1,1,1,1,1]
 obsdataset =["ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI","ERAI","ERAI"]
 plotrain=["" for x in range(nsite)]


 for ire in range (0, nsite):
     if not os.path.exists(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N'):
         os.mkdir(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N')

     plotname = casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N/Energy_'+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason
     plotrain[ire] = 'Energy_'+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason
     wks= Ngl.open_wks(ptype,plotname)

     Ngl.define_colormap(wks,"default")
     plot = []
     res     = Ngl.Resources()
     res.nglDraw              = False
     res.nglFrame             = False
     res.lgLabelFontHeightF     = .012                   # change font height
     res.lgPerimOn              = False                 # no box around
     res.vpWidthF         = 0.30                      # set width and height
     res.vpHeightF        = 0.30
     #res.vpXF             = 0.04
     # res.vpYF             = 0.30
     res.tmYLLabelFont  = 12
     res.tmXBLabelFont  = 12
     res.tmXBLabelFontHeightF = 0.005
     res.tmXBLabelFontThicknessF = 1.0
#     res.tmXBLabelAngleF = 45
     res.xyMarkLineMode      = "Lines"
     res.xyLineThicknesses = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.,2.,2.,2.,2,2,2,2,2,2,2]

     res.xyDashPatterns    = np.arange(0,24,1)
#     res.xyMarkers         = np.arange(16,40,1)
#     res.xyMarkerSizeF       = 0.005


     pres            = Ngl.Resources()
#     pres.nglMaximize = True
     pres.nglFrame = False
     pres.txFont = 12
     pres.nglPanelYWhiteSpacePercent = 5
     pres.nglPanelXWhiteSpacePercent = 5
     pres.nglPanelTop = 0.93

     for iv in range (0, nvaris):   
         if(iv == nvaris-1):
             res.pmLegendDisplayMode    = "ALWAYS"
             res.xyExplicitLegendLabels = casenames[:]
             res.pmLegendSide           = "top"             
             res.pmLegendParallelPosF   = 0.6               
             res.pmLegendOrthogonalPosF = -0.5                  
             res.pmLegendWidthF         = 0.10              
             res.pmLegendHeightF        = 0.10          
             res.lgLabelFontHeightF     = .02               
             res.lgLabelFontThicknessF  = 1.5
             res.lgPerimOn              = False
         else:
             res.pmLegendDisplayMode    = "NEVER"


         for im in range (0,ncases):
             ncdfs[im]  = './data/'+cases[im]+'_site_location.nc'
             infiles[im]= filepath[im]+cases[im]+'/'+cases[im]+'_'+cseason+'_climo.nc'
             inptrs = Dataset(infiles[im],'r')       # pointer to file1
             lat=inptrs.variables['lat'][:]
             nlat=len(lat)
             lon=inptrs.variables['lon'][:]
             nlon=len(lon)
             lev=inptrs.variables['lev'][:]
             nlev=len(lev)
             ncdf= Dataset(ncdfs[im],'r')
             n   =ncdf.variables['n'][:]
             idx_cols=ncdf.variables['idx_cols'][:,:]
             ncdf.close()
             if (im ==0):
                 A_field = np.zeros((ncases,nlev),np.float32)

             for subc in range( 0, n[ire]):
                 npoint=idx_cols[ire,n[subc]-1]-1
                 tmp=inptrs.variables[varis[iv]][0,:,npoint] 
                 theunits=str(cscale[iv])+inptrs.variables[varis[iv]].units
                 A_field[im,:] = (A_field[im,:]+tmp[:]/n[ire]).astype(np.float32 )
             A_field[im,:] = A_field[im,:] *cscale[iv]

             inptrs.close()
         res.tiMainString    =  varis[iv]+"  "+theunits
#         res.trXMinF = min(np.min(A_field[0, :]))
#         res.trXMaxF = max(np.max(A_field[0, :]))
         res.trYReverse        = True
         res.xyLineColors      = np.arange(2,20,1)
         res.xyMarkerColors    = np.arange(2,20,2)
         p = Ngl.xy(wks,A_field,lev,res)
         
#         res.trYReverse        = False
#         res.xyLineColors      = ["black"]
#         pt = Ngl.xy(wks,B,pre1,res)
#         Ngl.overlay(p,pt)

         plot.append(p)


     pres.txString   = "Rain at"+ str(lons[ire])+"E,"+str(lats[ire])+"N"
     txres = Ngl.Resources()
     txres.txFontHeightF = 0.020
     Ngl.text_ndc(wks,"Rain at"+ str(lons[ire])+"E,"+str(lats[ire])+"N",0.5,0.95,txres)

     Ngl.panel(wks,plot[:],[nvaris/2,2],pres)
     Ngl.frame(wks)
     Ngl.destroy(wks)

 return plotrain


     

