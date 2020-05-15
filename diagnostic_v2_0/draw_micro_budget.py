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

 
def draw_micro_bgt (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs,casedir,varis,vname,cscale,chscale,pname):

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
 mkres.gsMarkerColor = 'Red'
 mkres.gsMarkerSizeF = 15.   
 infiles  = ['' for x in range(ncases)]
 ncdfs    = ['' for x in range(ncases)]
 nregions = nsite

 nvaris = len(varis)

 plotmicrobgt=['' for x in range(nsite*ncases)] 

 for ire in range (0, nsite):
     for im in range (0,ncases):
         if not os.path.exists(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N'):
             os.mkdir(casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N')

         plotname = casedir+'/'+str(lons[ire])+'E_'+str(lats[ire])+'N/'+pname+'_'+casenames[im]+'_'+str(lons[ire])+'E_'+str(lats[ire])+'N_'+cseason
         plotmicrobgt[im+ncases*ire] = pname+'_'+casenames[im]+'_'+str(lons[ire])+'E_'+str(lats[ire])+'N_'+cseason

         wks= Ngl.open_wks(ptype,plotname)

         Ngl.define_colormap(wks,'radar')
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
         res.tmXBLabelFontHeightF = 0.01
         res.tmXBLabelFontThicknessF = 1.0
         res.xyMarkLineMode      = 'MarkLines'
         res.xyLineThicknesses = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0,3.,3.,3.,3.,3,3,3,3,3,3,3]
         res.xyLineColors      = np.arange(2,16,1)
         res.xyDashPatterns    = np.arange(0,24,1)
         res.xyMarkers         = np.arange(16,40,1)
         res.xyMarkerSizeF       = 0.005
         res.xyMarkerColors      = np.arange(2,16,1)
         res.pmLegendDisplayMode    = 'ALWAYS'
         res.pmLegendSide           = 'top'                 # Change location of
         res.pmLegendParallelPosF   = 0.6                  # move units right
         res.pmLegendOrthogonalPosF = -0.55                  # more neg = down
         res.pmLegendWidthF         = 0.2       # Decrease width
         res.pmLegendHeightF        = 0.1       # Decrease height
         res.lgBoxMinorExtentF      = 0.1       # Shorten the legend lines
         res.lgLabelFontHeightF     = 0.015     # Change the font size
         res.lgPerimOn              = True
         res.tiYAxisString   = 'PRESSURE'
     
#         res.nglLeftString     = varis[iv]
#         res.nglRightString    = cunits[iv]
         res.trYReverse        = True

         pres            = Ngl.Resources() 
#         pres.nglMaximize = True

         pres.nglFrame = False
         pres.txFont = 12
         pres.nglPanelYWhiteSpacePercent = 5
         pres.nglPanelXWhiteSpacePercent = 5
         pres.nglPanelTop = 0.93

         txres               = Ngl.Resources()
#         txres.txFontHeightF = 0.01

         for iv in range (0, nvaris):

             if (varis[iv] == 'MPDLIQ' ):   # LIQ
                budget_ends = ['PRCO',  'PRAO', 'MNUCCCO', 'MNUCCTO', 'MSACWIO', 'PSACWSO', 'BERGSO','BERGO']
              # in fortran    prc*cld, pra*cld,mnuccc*cld,mnucct*cld,msacwi*cld,psacws*cld,bergs*cld, berg
              #               (-pra-prc-mnuccc-mnucct-msacwi- psacws-bergs)*lcldm-berg
                nterms = len (budget_ends)

             if (varis[iv] == 'MPDICE' ):    # ICE
                budget_ends = [ 'PRCIO', 'PRAIO', 'MSACWIO', 'MNUCCCO',  'MNUCCTO','mnudepo',   'BERGO',  'CMEIOUT','mnuccrio']
# in fortran                   prci*cld,prai*cld,msacwi*cld,mnuccc*cld, mnucct*cld,      berg, vap_dep + ice_sublim + mnuccd 
#                (mnuccc+mnucct+mnudep+msacwi)*lcldm+(-prci-prai)*icldm+(vap_dep+ice_sublim+mnuccd)+berg+mnuccri*precip_frac  
                nterms = len (budget_ends)

             if (varis[iv] == 'QRSEDTEN' ):  #  RAIN
                budget_ends = [ 'PRAO', 'PRCO',  'PRACSO',       'EVAPPREC', 'MNUCCRO','mnuccrio' ]
#                              pra*cld,prc*cld, psacs*prf, -pre*prf(nevapr),mnuccr*prf
#             (pra+prc)*lcldm+(pre-pracs- mnuccr-mnuccri)*precip_frac   
                nterms = len (budget_ends)

             if (varis[iv] == 'QSSEDTEN' ):  # SNOW
                budget_ends = [ 'PRAIO', 'PRCIO', 'PSACWSO', 'PRACSO','EVAPSNOW',  'MNUCCRO', 'BERGSO']
#                              prai*cld,prci*cld,psacws*cld,psacs*prf, -prds*prc, mnuccr*prf,bergs*cld 
#              (prai+prci)*icldm+(psacws+bergs)*lcldm+(prds+  pracs+mnuccr)*precip_frac
                nterms = len (budget_ends)

             if (varis[iv] == 'QISEVAP' ):  # Vapor
                budget_ends = [ 'EVAPPREC','EVAPSNOW','CMEIOUT','mnudepo' ]
               #          -pre*prf(nevapr),-prds*prc ,vap_dep + ice_sublim + mnuccd
#            -(pre+prds)*precip_frac-vap_dep-ice_sublim-mnuccd-mnudep*lcldm 
                nterms = len (budget_ends)

             if (varis[iv] == 'nnuccco' ):  # NUM of LIQ
                budget_ends = [ 'nnuccco', 'nnuccto', 'npsacwso', 'nsubco', 'nprao','nprc1o']
#                               nnuccc*cld,nnucct*cld,npsacws*cld,nsubc*cld,npra*cld,nprc1*cld
#                              (-nnuccc-nnucct-npsacws+nsubc-npra-nprc1)*lcldm
                nterms = len (budget_ends)
                 
             if (varis[iv] == 'nnuccdo' ):  # NUM of ICE
                budget_ends = [  'nnuccdo',  'nnuccto',  'tmpfrzo',  'nnudepo',  'nsacwio',  'nsubio',  'nprcio',  'npraio','nnuccrio','DETNICETND']
#                                nnuccd   ,nnucct*lcld,tmpfrz*lcld,nnudep*lcld,nsacwi*lcld,nsubi*icld,nprci*icld,nprai*icld,nnuccri*prf
#                                nnuccd+ (nnucct+tmpfrz+nnudep+nsacwi)*lcldm+(nsubi-nprci- nprai)*icldm+nnuccri*precip_frac
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
             theunits=str(chscale[iv])+'x'+inptrs.variables[varis[iv]].units
             res.tiMainString    =  vname[iv]+'  '+theunits 


             for it in range(0, nterms):
                 for subc in range( 0, n[ire]):
                     varis_bgt= budget_ends[it]
                     npoint=idx_cols[ire,n[subc]-1]-1
                     tmp=inptrs.variables[varis_bgt][0,:,npoint] #/n[ire]
                     tmp=tmp*cscale[iv]
                     lcldm=inptrs.variables['CLOUD'][0,:,npoint] 
                     icldm=lcldm
                     precip_frac=inptrs.variables['FREQR'][0,:,npoint]

                     if (varis_bgt == 'MPDT' or varis_bgt == 'STEND_CLUBB' ):
                        tmp=tmp/1004

                     if (varis[iv] == 'MPDLIQ'):  # LIQ
                       if (varis_bgt == 'PRCO' or varis_bgt ==  'PRCIO' or varis_bgt == 'PRAO' or varis_bgt == 'PRAICSO' or varis_bgt == 'PRAIO' \
                          or varis_bgt == 'MNUCCCO' or varis_bgt == 'MNUCCTO' or varis_bgt == 'MSACWIO' or varis_bgt == 'PSACWSO' or varis_bgt == 'BERGSO' 
                          or varis_bgt == 'BERGO'):
                          tmp=tmp *(-1) 

                     if (varis[iv] == 'MPDICE'):  # ICE
                        if ( varis_bgt == 'PRCIO' or varis_bgt == 'PRAIO' ):                 
                           tmp=tmp *(-1)

                     if (varis[iv] == 'QRSEDTEN'):  # RAIN
                        if ( varis_bgt == 'MNUCCRO'  or varis_bgt == 'PRACSO' or varis_bgt == 'EVAPPREC' or varis_bgt == 'mnuccrio'):
                           tmp=tmp*(-1)

                     if (varis[iv] == 'QSSEDTEN'):  # SNOW
                        if ( varis_bgt == 'EVAPSNOW' ):
                           tmp= tmp*(-1)

                     if (varis[iv] == 'QISEVAP'):  # Vapor
                        if ( varis_bgt == 'CMEIOUT' or varis_bgt == 'mnudepo' ):
                           tmp=-1* tmp

                     if (varis[iv] == 'nnuccco' ):  # NUM of LIQ
                        if ( varis_bgt == 'nnuccco' or varis_bgt == 'nnuccto' or varis_bgt == 'npsacwso' or varis_bgt == 'nprao' or varis_bgt == 'nprc1o'):
                           tmp=-1* tmp 

                     if (varis[iv] == 'nnuccdo' ):  # NUM of ICE
                        if ( varis_bgt == 'nprcio' or varis_bgt == 'npraio'):
                           tmp=-1* tmp

                     A_field[it,:] = (A_field[it,:]+tmp[:]/n[ire]).astype(np.float32 )

             inptrs.close()
             res.xyExplicitLegendLabels =  budget_ends[:]
             p = Ngl.xy(wks,A_field,ilev,res)
             plot.append(p)

             xp=np.mod(iv,2)
             yp=int(iv/2)

         pres.txFontHeightF = 0.02
         pres.txFont = _Font
         pres.txString   = casenames[im]+' Microphy BUDGET at' +str(lons[ire])+'E,'+str(lats[ire])+'N'

         if(np.mod(nvaris,2)==1):
            Ngl.panel(wks,plot[:],[(nvaris)/2+1,2],pres)
         else:
            Ngl.panel(wks,plot[:],[(nvaris)/2,2],pres)

         txres = Ngl.Resources()
         txres.txFontHeightF = 0.020
         txres.txFont = _Font
         Ngl.text_ndc(wks,casenames[im]+' Microphy BUDGET at' +str(lons[ire])+'E,'+str(lats[ire])+'N',0.5,0.95,txres)

         Ngl.frame(wks)
         Ngl.destroy(wks) 

 return (plotmicrobgt)
