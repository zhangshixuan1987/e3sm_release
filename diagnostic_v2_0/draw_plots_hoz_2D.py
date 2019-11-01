
'''
    Draw 2D plots
    Zhun Guo 
'''
import Ngl
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import os
import Common_functions
from subprocess import call


def draw_2D_plot (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir):

# ncases, the number of models
# cases, the name of models
# casename, the name of cases
# filepath, model output filepath
# filepathobs, filepath for observational data
# inptrs = [ncases]
 if not os.path.exists(casedir):
        os.mkdir(casedir)

 if not os.path.exists(casedir+"/2D"):
        os.mkdir(casedir+"/2D") 

 _Font   = 25
 interp = 2
 extrap = False
 mkres = Ngl.Resources()
 mkres.gsMarkerIndex = 2
 mkres.gsMarkerColor = "Red"
 mkres.gsMarkerSizeF = 15.   
 infiles  = ["" for x in range(ncases)] 
 ncdfs    = ["" for x in range(ncases)] 
 varis    = ["CLDTOT",     "SWCF","LWCF","PRECT"]
 varisobs = ["CLDTOT_CS2", "SWCF","LWCF","PRECT"]
 nvaris = len(varis)
 cunits = [""]
 cscale = [100,1,1,86400000, 86400,  1]
 cscaleobs =  [1,1,1,1,0.04]
# cntrs = [[0 for col in range(11)] for row in range(nvaris)]
 cntrs = np.zeros((nvaris,11),np.float32)

 obsdataset=  ["CLOUDSATCOSP","CERES","CERES","GPCP","ERAI", "ERAI", "ERAI", "ERAI"]
 
 plot2d=["" for x in range(nvaris)]
 for iv in range(0, nvaris):
# make plot for each field 
   if(varis[iv] == "CLDTOT"):
       cntrs[iv,:] = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 90]
   if(varis[iv] == "LWCF"):
       cntrs[iv,:] = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110] 
   if(varis[iv] =="SWCF"):
       cntrs[iv,:] = [-30, -40, -50, -60, -70, -80, -90, -100, -110, -120, -130]
   if(varis[iv]=="PRECT" or varis[iv]=="QFLX"):
       cntrs[iv,:] = [0.1, 0.2, 0.5, 1.0, 2, 3, 4, 5, 6, 7, 8]
   if(varis[iv] == "FLUT"):
       cntrs[iv,:] = [210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310]
   
    

#  Observational data
   if(obsdataset[iv] =="CCCM"):
       if(cseason == "ANN"):
           fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-"+cseason+".nc"
       else:
           fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-2010-"+cseason+".nc"
   else:
       if (varisobs[iv] =="PRECT"):
           fileobs = filepathobs+'/GPCP_'+cseason+'_climo.nc'
       else:
           fileobs = filepathobs + obsdataset[iv]+'_'+cseason+'_climo.nc'

#       infiles[im]=filepath[im]+cases[im]+'/run/'+cases[im]+'_'+cseason+'_climo.nc'


   inptrobs = Dataset(fileobs,'r') 
   latobs=inptrobs.variables['lat'][:]
   lonobs=inptrobs.variables['lon'][:]
   B=inptrobs.variables[varisobs[iv]][0,:,:] 
   B=B * cscaleobs[iv]

   #************************************************
   # create plot
   #************************************************
   plotname = casedir+'/2D/Horizontal_'+varis[iv]+'_'+cseason
   plot2d[iv] = 'Horizontal_'+varis[iv]+'_'+cseason
   wks= Ngl.open_wks(ptype,plotname)
 
#   ngl_define_colormap(wks,"prcp_1")
   plot = []

   txres               = Ngl.Resources()
   txres.txFontHeightF = 0.018

   pres            = Ngl.Resources()
#   pres.nglMaximize = True
   pres.nglPanelYWhiteSpacePercent = 5
   pres.nglPanelXWhiteSpacePercent = 5
   pres.nglPanelBottom = 0.02
   pres.nglPanelLabelBar        = True
   pres.pmLabelBarWidthF        = 0.8
   pres.nglFrame         = False


   res = Ngl.Resources()
   res.nglDraw         =  False            #-- don't draw plots
   res.nglFrame        =  False  
   res.cnFillOn     = True
   res.cnFillMode   = "RasterFill"
   res.cnLinesOn    = False
   res.nglMaximize = True
   res.mpFillOn     = True
   res.mpCenterLonF = 180
   res.tiMainFont                     = _Font
   res.tiMainFontHeightF              = 0.016
   res.tiXAxisString                  = ""
   res.tiXAxisFont                    = _Font
   res.tiXAxisFontHeightF             = 0.016
   res.tiYAxisString                  = ""
   res.tiYAxisFont                    = _Font
   res.tiYAxisOffsetXF                = 0.0
   res.tiYAxisFontHeightF             = 0.018        
   res.tmXBLabelFont = _Font
   res.tmYLLabelFont = _Font
   res.tiYAxisFont   = _Font
   res.vpWidthF         = 0.80                      # set width and height
   res.vpHeightF        = 0.40
   res.vpXF             = 0.04
   res.vpYF             = 0.30

 
#   res.nglStringFont                  = _Font
#   res.nglStringFontHeightF           = 0.04
#   res.nglRightString                 = ""#"Cloud Fraction"
#   res.nglScalarContour     = True         
   res.cnInfoLabelOn                  = False  
   res.cnFillOn                       = True
   res.cnLinesOn                      = False
   res.cnLineLabelsOn                 = False
   res.lbLabelBarOn                   = False
 
#   res.vcRefMagnitudeF = 5.
#   res.vcMinMagnitudeF = 1.
#   res.vcRefLengthF    = 0.04
#   res.vcRefAnnoOn     = True#False
#   res.vcRefAnnoZone   = 3
#   res.vcRefAnnoFontHeightF = 0.02
#   res.vcRefAnnoString2 =""
#   res.vcRefAnnoOrthogonalPosF   = -1.0   
#  res.vcRefAnnoArrowLineColor   = "blue"         # change ref vector color
#  res.vcRefAnnoArrowUseVecColor = False
#   res.vcMinDistanceF  = .05
#   res.vcMinFracLengthF         = .
#   res.vcRefAnnoParallelPosF    =  0.997
#   res.vcFillArrowsOn           = True
#   res.vcLineArrowThicknessF    =  3.0
#   res.vcLineArrowHeadMinSizeF   = 0.01
#   res.vcLineArrowHeadMaxSizeF   = 0.03
#   res.vcGlyphStyle              = "CurlyVector"     # turn on curley vectors 
#  res@vcGlyphStyle              ="Fillarrow"
#   res.vcMonoFillArrowFillColor = True
#   res.vcMonoLineArrowColor     = True
#   res.vcLineArrowColor          = "green"           # change vector color
#   res.vcFillArrowEdgeColor      ="white"
#   res.vcPositionMode            ="ArrowTail"
#   res.vcFillArrowHeadInteriorXF =0.1
#   res.vcFillArrowWidthF         =0.05           #default
#   res.vcFillArrowMinFracWidthF  =.5
#   res.vcFillArrowHeadMinFracXF  =.5
#   res.vcFillArrowHeadMinFracYF  =.5
#   res.vcFillArrowEdgeThicknessF = 2.0
   res.mpFillOn                   = False
   
   res.cnLevelSelectionMode = "ExplicitLevels"

   res.cnLevels      = cntrs[iv][:]

   for im in range(0, ncases):
       ncdfs[im]  = './data/'+cases[im]+'_site_location.nc' 
       infiles[im]= filepath[im]+cases[im]+'/'+cases[im]+'_'+cseason+'_climo.nc'
       inptrs = Dataset(infiles[im],'r')       # pointer to file1
       lat=inptrs.variables['lat'][:]
       nlat=len(lat)
       lon=inptrs.variables['lon'][:]
       nlon=len(lon)
       area=inptrs.variables['area'][:]

       area_wgt = np.zeros(nlat)
#       area_wgt[:] = gw[:]

       sits=np.linspace(0,nsite-1,nsite)
       ncdf= Dataset(ncdfs[im],'r')
       n   =ncdf.variables['n'][:]
       idx_cols=ncdf.variables['idx_cols'][:]
       if (varis[iv] == 'PRECT'):
           A = inptrs.variables['PRECC'][0,:]+inptrs.variables['PRECL'][0,:]
       else:
           A = inptrs.variables[varis[iv]][0,:]

       A_xy=A
       A_xy = A_xy * cscale[iv]
       ncdf.close()
       inptrs.close()

       if im == 0 :
           dsizes = len(A_xy)
           field_xy = [[0 for col in range(dsizes)] for row in range(ncases)] 
       
       field_xy[im][:] = A_xy
  
       res.lbLabelBarOn = False 
       if(np.mod(im,2)==0): 
           res.tiYAxisOn  = True
       else:
           res.tiYAxisOn  = False
  
       res.tiXAxisOn  = False 
       res.sfXArray     = lon
       res.sfYArray     = lat
       res.mpLimitMode  = "LatLon"
       res.mpMaxLonF    = max(lon) 
       res.mpMinLonF    = min(lon) 
       res.mpMinLatF    = min(lat) 
       res.mpMaxLatF    = max(lat) 
       res.tiMainString    =  casenames[im]+"  " +varis[iv] +"  GLB="+str(np.sum(A_xy[:]*area[:]/np.sum(area)))

       p = Ngl.contour_map(wks,A_xy,res)
       plot.append(p)

# observation 
#   res.nglLeftString = obsdataset[iv]
#  res@lbLabelBarOn = True 
#  res@lbOrientation        = "vertical"         # vertical label bars
   res.lbLabelFont          = 12
   res.tiYAxisOn  = True
   res.tiXAxisOn  = True
   res.tiXAxisFont = 12
   rad    = 4.0*np.arctan(1.0)/180.0
   re     = 6371220.0
   rr     = re*rad

   dlon   = abs(lonobs[2]-lonobs[1])*rr
   dx     = dlon* np.arccos(latobs*rad)
   jlat   = len(latobs )
   dy     = np.zeros(jlat,dtype=float)
                                                            # close enough
   dy[0]  = abs(lat[2]-lat[1])*rr
   dy[1:jlat-2]  = abs(lat[2:jlat-1]-lat[0:jlat-3])*rr*0.5
   dy[jlat-1]    = abs(lat[jlat-1]-lat[jlat-2])*rr
   area_wgt   = dx*dy  # 
   is_SE= False

   res.sfXArray     = lonobs
   res.sfYArray     = latobs
   res.mpLimitMode  = "LatLon"
   res.mpMaxLonF    = max(lonobs)
   res.mpMinLonF    = min(lonobs)
   res.mpMinLatF    = min(latobs)
   res.mpMaxLatF    = max(latobs)
   res.tiMainString    =  obsdataset[iv]+"  " +varis[iv] +"  GLB="+str(Common_functions.area_avg(B, area_wgt,is_SE))


   p =Ngl.contour_map(wks,B,res)
   plot.append(p)


   Ngl.panel(wks,plot[:],[ncases+1,1],pres) 

   Ngl.frame(wks)
   Ngl.destroy(wks)

#   Ngl.end()
 return plot2d


