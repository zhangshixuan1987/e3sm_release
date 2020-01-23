'''
    CLUBB budgets
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

 
def draw_e3sm_bgt (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,dpsc):

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

 varis = [ "DCQ","DCCLDLIQ","DCCLDICE","PTEQ","PTTEND","DTCOND"]
 varisobs = [ "CLOUD", "OMEGA","SHUM","CLWC_ISBL", "THATA","RELHUM"]
 nvaris = len(varis)
 cunits = ["%", "mba/day","g/kg","g/kg","K", "%", "mba/day", "K", "g/kg", "m/s", "m/s","K","m"]
 cscale = [1, 1, 1, 1, 1, 1, 1,1,1E3,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
 cscaleobs = [100., 100/86400., 1., 1000, 1., 1., 1, 1,1,1]
 obsdataset =["CCCM", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI","ERAI","ERAI"]

 plote3smbgt=["" for x in range(nsite*ncases)] 

 for ire in range (0, nsite):
     for im in range (0,ncases):
         if not os.path.exists(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N'):
             os.mkdir(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N')

         plotname = casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N/E3SM_Budgets_'+casenames[im]+"_"+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason
         plote3smbgt[im+ncases*ire] = 'E3SM_Budgets_'+casenames[im]+"_"+str(lons[ire])+"E_"+str(lats[ire])+"N_"+cseason

         wks= Ngl.open_wks(ptype,plotname)

         Ngl.define_colormap(wks,"radar_1")
         plot = []
         res     = Ngl.Resources()  
         res.nglDraw              = False
         res.nglFrame             = False
         res.lgLabelFontHeightF     = .012                   # change font height
         res.lgPerimOn              = False                 # no box around
         res.vpWidthF         = 0.30                      # set width and height
         res.vpHeightF        = 0.30

#         res.txFontHeightF   = .01
         #  res.vpXF             = 0.04
         # res.vpYF             = 0.30
         res.tmYLLabelFont  = 12
         res.tmXBLabelFont  = 12
         res.tmXBLabelFontHeightF = 0.005
         res.tmXBLabelFontThicknessF = 1.0
         res.xyMarkLineMode      = "MarkLines"
         res.xyLineThicknesses = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.,2.,2.,2.,2,2,2,2,2,2,2]
         res.xyLineColors      = np.arange(2,25,2)
         res.xyDashPatterns    = np.arange(0,24,1)
         res.xyMarkers         = np.arange(16,40,1)
         res.xyMarkerSizeF       = 0.005
         res.xyMarkerColors      = np.arange(2,25,2)
         res.pmLegendDisplayMode    = "ALWAYS"
         res.pmLegendSide           = "top"                 # Change location of
         res.pmLegendParallelPosF   = 0.75                  # move units right
         res.pmLegendOrthogonalPosF = -0.65                  # more neg = down
         res.pmLegendWidthF         = 0.1                  # Change width and
         res.pmLegendHeightF        = 0.1                  # height of legend.
         res.lgLabelFontHeightF     = 1.                   # change font height
         res.lgLabelFontThicknessF  = 1.
         res.lgPerimOn              = True
         res.tiYAxisString   = "PRESSURE"
     
#         res.nglLeftString     = varis[iv]
#         res.nglRightString    = cunits[iv]
         res.trYReverse        = True

         pres            = Ngl.Resources() 
#         pres.nglMaximize = True
         pres.wkWidth              = 2000
         pres.wkHeight             = 2000

         pres.nglFrame = False
         pres.txFont = 12
         pres.nglPanelYWhiteSpacePercent = 5
         pres.nglPanelXWhiteSpacePercent = 5
         pres.nglPanelTop = 0.93

         txres               = Ngl.Resources()
#         txres.txFontHeightF = 0.01

         for iv in range (0, nvaris):

             if (varis[iv] == "DCQ" ):
                if (dpsc[im] == "zm" ):
                   budget_ends = ["MPDQ", "RVMTEND_CLUBB","ZMDQ", "EVAPQZM"]
                else:
                   budget_ends = ["MPDQ", "RVMTEND_CLUBB"]
                nterms = len (budget_ends)
             if (varis[iv] == "DTCOND" ):
                if (dpsc[im] == "zm" ):
                   budget_ends = ["STEND_CLUBB", "MPDT", "DPDLFT","ZMDT", "EVAPTZM", "ZMMTT"]
                else:
                   budget_ends = ["STEND_CLUBB", "MPDT", "DPDLFT"]
                nterms = len (budget_ends)
             if (varis[iv] == "PTTEND") :
                 budget_ends = ["DTCOND", "QRS", "QRL",  "TTGW"]
                 nterms = len (budget_ends)
             if (varis[iv] == "PTEQ") :
                 if (dpsc[im] == "zm" ):
                    budget_ends = ["MPDQ", "RVMTEND_CLUBB","ZMDQ", "EVAPQZM"]
                 else:
                    budget_ends = ["MPDQ", "RVMTEND_CLUBB"]  
                 nterms = len (budget_ends)
             if (varis[iv] == "DCCLDLIQ") :
                 if (dpsc[im] == "zm" ):
                    budget_ends = ["MPDLIQ", "RCMTEND_CLUBB", "DPDLFLIQ","ZMDLIQ"]
                 else:
                    budget_ends = ["MPDLIQ", "RCMTEND_CLUBB", "DPDLFLIQ"]
                 nterms = len (budget_ends)
             if (varis[iv] == "DCCLDICE") :
                 if (dpsc[im] == "zm" ):
                    budget_ends = ["MPDICE", "RIMTEND_CLUBB", "DPDLFICE","ZMDICE"]
                 else:
                    budget_ends = ["MPDICE", "RIMTEND_CLUBB", "DPDLFICE"]
                 nterms = len (budget_ends)


             ncdfs[im]  = './data/'+cases[im]+'_site_location.nc'
             infiles[im]= filepath[im]+cases[im]+'/'+cases[im]+'_'+cseason+'_climo.nc'
             inptrs = Dataset(infiles[im],'r')       # pointer to file1
             lat=inptrs.variables['lat'][:]
             nlat=len(lat)
             lon=inptrs.variables['lon'][:]
             nlon=len(lon)
             ilev=inptrs.variables['lev'][:]
             nilev=len(ilev)
             ncdf= Dataset(ncdfs[im],'r')
             n   =ncdf.variables['n'][:]
             idx_cols=ncdf.variables['idx_cols'][:,:]
             ncdf.close()
             A_field = np.zeros((nterms,nilev),np.float32)
             theunits=str(cscale[iv])+"x"+inptrs.variables[varis[iv]].units
             res.tiMainString    =  varis[iv]+"  "+theunits 


             for it in range(0, nterms):
                 for subc in range( 0, n[ire]):
                     varis_bgt= budget_ends[it]
                     npoint=idx_cols[ire,n[subc]-1]-1
                     tmp=inptrs.variables[varis_bgt][0,:,npoint] #/n[ire]
                     tmp=tmp*cscale[iv]
                     if (varis_bgt == "MPDT" or varis_bgt == "STEND_CLUBB" ):
                        tmp=tmp/1004
                     A_field[it,:] = (A_field[it,:]+tmp[:]/n[ire]).astype(np.float32 )

             inptrs.close()
             res.xyExplicitLegendLabels =  budget_ends[:]
             p = Ngl.xy(wks,A_field,ilev,res)
             plot.append(p)

             xp=np.mod(iv,2)
             yp=int(iv/2)

         pres.txFontHeightF = 0.02
         pres.txString   = casenames[im]+"  BUDGET at" +str(lons[ire])+"E,"+str(lats[ire])+"N"

         Ngl.panel(wks,plot[:],[nvaris/2,2],pres)
         txres = Ngl.Resources()
         txres.txFontHeightF = 0.020
         Ngl.text_ndc(wks,casenames[im]+"  BUDGET at" +str(lons[ire])+"E,"+str(lats[ire])+"N",0.5,0.95,txres)

         Ngl.frame(wks)
         Ngl.destroy(wks) 

 return (plote3smbgt)

