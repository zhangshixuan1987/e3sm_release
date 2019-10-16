'''
    Silhs variables 
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


def silhs_prf (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs):

# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data
# inptrs = [ncases]
 if not os.path.exists("pic"):
        os.mkdir("pic")

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

 varis    = [ "SILHS_CLUBB_PRECIP_FRAC","SILHS_CLUBB_ICE_SS_FRAC","HR","AWNC","AWNI","ADSNOW","ANSNOW","AQSNOW","RHW","QRS","QRL","AQRAIN" ]
 nvaris = len(varis)
 cunits = ["%","mba/day","g/kg","g/kg","K", "%", "m/s", "g/kg", "m/s", "m/s","K","m" ]
 cscale = [1, 1, 1, 1, 1., 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 ]

 plotsilhs=["" for x in range(nsite)]

 for ire in range (0, nsite):
     if not os.path.exists('pic/'+str(lons[ire])+'E_'+str(lats[ire])+'N'):
         os.mkdir('pic/'+str(lons[ire])+'E_'+str(lats[ire])+'N')

     plotname = './pic/'+str(lons[ire])+'E_'+str(lats[ire])+'N/silhs_'+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason
     plotsilhs[ire] = 'CLUBB_standard_'+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason
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
     res.txFontHeightF   = .01
     #res.vpXF             = 0.04
     # res.vpYF             = 0.30
     res.tmYLLabelFont  = 12
     res.tmXBLabelFont  = 12
     res.tmXBLabelFontHeightF = 0.005
     res.tmXBLabelFontThicknessF = 1.0
     res.tmXBLabelAngleF = 45
     res.xyMarkLineMode      = "Lines"
     res.xyLineThicknesses = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.,2.,2.,2.,2,2,2,2,2,2,2]

     res.xyDashPatterns    = np.arange(0,24,1)
#     res.xyMarkers         = np.arange(16,40,1)
#     res.xyMarkerSizeF       = 0.005


     pres            = Ngl.Resources()
#     pres.nglMaximize = True

     pres.txString   = str(lons[ire])+"E,"+str(lats[ire])+"N"
     pres.txFont = 12
     pres.nglPanelYWhiteSpacePercent = 5
     pres.nglPanelXWhiteSpacePercent = 5
     pres.nglPanelTop = 0.93

     txres               = Ngl.Resources()
     txres.txFontHeightF = 0.01
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
             lev=inptrs.variables['ilev'][:]
             nlev=len(lev)
             ncdf= Dataset(ncdfs[im],'r')
             n   =ncdf.variables['n'][:]
             idx_cols=ncdf.variables['idx_cols'][:,:]
             ncdf.close()
             if (im ==0):
                 A_field = np.zeros((ncases,nlev),np.float32)

             for subc in range( 0, n[ire]):
                 npoint=idx_cols[ire,n[subc]-1]-1
                 if (varis[iv] == 'THETA'):
                     tmp = inptrs.variables['T'][0,:,npoint]
                     hyam =inptrs.variables['hyam'][:]
                     hybm =inptrs.variables['hybm'][:]
                     ps=inptrs.variables['PS'][0,npoint] 
                     ps=ps
                     p0=inptrs.variables['P0']
                     pre = np.zeros((nlev),np.float32)
                     for il in range (0, nlev):
                         pre[il] = hyam[il]*p0 + hybm[il] * ps
                         tmp[il] = tmp[il] * (100000/pre[il])**0.286
                     theunits=str(cscale[iv])+inptrs.variables['T'].units

                 else:
                     tmp=inptrs.variables[varis[iv]][0,:,npoint] 
                     theunits=str(cscale[iv])+inptrs.variables[varis[iv]].units
                 A_field[im,:] = (A_field[im,:]+tmp[:]/n[ire]).astype(np.float32 )
             A_field[im,:] = A_field[im,:] *cscale[iv]

             inptrs.close()
         res.tiMainString    =  varis[iv]+"  x"+theunits
         res.trXMinF = min(np.min(A_field[0, :]))#,np.min(B))
         res.trXMaxF = max(np.max(A_field[0, :]))#,np.max(B))
         if(varis[iv] == "THETA"):
             res.trXMinF = 280.
             res.trXMaxF = 400.
         if(varis[iv] == "CLOUD" or varis[iv] =="RELHUM") :
             res.trXMinF = 0.
             res.trXMaxF = 100.
         if(varis[iv] == "T") :
             res.trXMinF = 180
             res.trXMaxF = 300
         if(varis[iv] == "U") :
             res.trXMinF = -40
             res.trXMaxF = 40
         res.trYReverse        = True
         res.xyLineColors      = np.arange(2,20,1)
         res.xyMarkerColors    = np.arange(2,20,2)
         p = Ngl.xy(wks,A_field,lev,res)
         
         plot.append(p)


     Ngl.panel(wks,plot[:],[nvaris/3,3],pres)
     txres = Ngl.Resources()
     txres.txFontHeightF = 0.020
     Ngl.text_ndc(wks,"SILHS VAR at"+ str(lons[ire])+"E,"+str(lats[ire])+"N",0.5,0.95,txres)


     Ngl.frame(wks)

     Ngl.destroy(wks)

 return plotsilhs


     

