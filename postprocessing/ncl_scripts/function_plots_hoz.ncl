
procedure plots_hoz_multiruns (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)
; ncases, the number of models
; casenames, the name of models
; casename, the name of cases
; filepath, model output filepath
; filepathobs, filepath for observational data

begin

 inptrs = new((/ncases/),  "file") ; pointers to files
 ncdfs   = inptrs

;=================================================================================================================
; Read sites info
;=================================================================================================================

; open files
 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+"_"+cseason+"_climo.nc"
 ;   print((/infile/))
    inptrs(im) = addfile(infile,"r")   
    ncdfs(im)  = addfile("./data/"+cases(im)+"_site_location.nc" ,"r")  
 end do

  mkres               = True
;  mkres@gsnDraw       = False
;  mkres@gsnFrame      = False  
  mkres@gsMarkerIndex = 2
  mkres@gsMarkerColor = "Red"
  mkres@gsMarkerSizeF = 15.        ; Increase marker sizes.

;=================================================================================================
 inptr = inptrs(0)


 delete(inptr)
ncdf = ncdfs(0)
n=ncdf->n
idx_cols=ncdf->ind_cols



    _Font              = 25
interp = 2
extrap = False
 clevel  = (/500,200,850/)
 nlevels =dimsizes (clevel) 
; fields to plot
 vars    = (/ "OMEGA","Z3"/)
 varsobs = (/  "OMEGA","Z3"/)
 nvars = dimsizes(vars)
 cunits = (/""/)
 cscale = (/864, 1,  1/)
 cscaleobs =  (/1. ,1,1/)
 cntrs = new ((/nvars,11/),float)
; observationl data set
 obsdataset =(/"ERAI", "ERAI", "ERAI", "ERAI"/)
 level =(/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,150.,100./)

 do i=0, nvars-1
   if(vars(i).eq."CLOUD") then
    cntrs(i, :) = (/0.01, 0.02, 0.05, 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9/)
   end if
   if(vars(i).eq."Q") then
    cntrs(i, :) = (/0.1, 0.5, 1, 2, 3, 5, 7, 9, 12, 15, 20/)
   end if
   if(vars(i).eq."RELHUM") then
    cntrs(i, :) = (/5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95/)
   end if
   if(vars(i).eq."OMEGA") then
    cntrs(i, :) = (/-0.15, -0.1, -0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, 0.1, 0.15/)*300
   end if
 end do

; open files
 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+"_"+cseason+"_climo.nc"
   ; print((/infile/))
    inptrs(im) = addfile(infile,"r")       ; pointer to file1
 end do 
 inptr=inptrs(0)
 lat   = dble2flt(inptr->lat)
 nlat  = dimsizes(lat)
 lon   = dble2flt(inptr->lon)
 nlon  = dimsizes(lon)
 lev   = inptr->lev
 nlev  = dimsizes(lev)
 clon = inptr->lon
 clat = inptr->lat
 delete(inptr)
 
 do iv=0, nvars-1
  do ire = 0, nlevels-1
; make plot for each field 
   do im=0, ncases-1
     inptr = inptrs(im)
     if (vars(iv) .eq. "PRECT" ) then
     A = dble2flt(inptr->$vars(iv)$(0,:))
     else
     A = dble2flt(inptr->$vars(iv)$(0,:, :))
     end if
  lev     = inptr->lev
  ilev    = inptr->ilev
  hyam    = inptr->hyam
  hybm    = inptr->hybm
  hyai    = inptr->hyai
  hybi    = inptr->hybi
  ps    = inptr->PS(0,:)                       ; get sfc pres
  p0    = (inptr->P0)/100  
     if (vars(iv) .eq. "PRECT" ) then
     A_xy=A
     else
     Z_xy=A 
     ;vinth2p(A,hyam,hybm,level,ps,interp,p0,1,extrap)
     A_xy=Z_xy({clevel(ire)},:);
     end if

     U = dble2flt(inptr->U(0,:, :))     
     U_xy=U
      ;vinth2p(U,hyam,hybm,level,ps,interp,p0,1,extrap)
     U_cl=U_xy({clevel(ire)},:)

     V = dble2flt(inptr->V(0,:,:))
     V_xy=V
     ;vinth2p(V,hyam,hybm,level,ps,interp,p0,1,extrap)
     V_cl=V_xy({clevel(ire)},:)

     A_xy = A_xy * cscale(iv)
    

     if(im.eq.0) then
      dsizes = dimsizes(A_xy)
      field_xy = new((/ncases, dsizes(0)/), float)
      field_xy!1 = "ncol"
      field_xy!0 = "model"
      delete(dsizes)
      u_fld=field_xy
      v_fld=field_xy
     end if
     
     field_xy(im, :) = ((/A_xy/))
     u_fld(im, :) = ((/U_cl/))
     v_fld(im, :) = ((/V_cl/))

     delete([/inptr,A,A_xy/])
   end do  ; end im

;  Observational data
   if(obsdataset(iv).eq."CCCM") then
     if(cseason.eq."ANN") then
       fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-"+cseason+".nc"
     else 
       fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-2010-"+cseason+".nc"
     end if
   else 
     if (obsdataset(iv).eq."PRECT")
     fileobs = filepathobs+"/GPCP_"+cseason+"_climo.nc"
     else
     fileobs = filepathobs + obsdataset(iv)+"_"+cseason+"_climo.nc"
     end if
   end if
   inptrobs = addfile(fileobs, "r")
   if(obsdataset(iv).eq."GPCP") then
     B = inptrobs->$varsobs(iv)$(0, :, :)
    fileobs1 = filepathobs + "ERAI_"+cseason+"_climo.nc"
    inptrobs1 = addfile(fileobs1, "r")
    Ua_cl =inptrobs1->U(0, {clevel(iv)}, :,:)
    Va_cl =inptrobs1->V(0, {clevel(iv)}, :,:)
    delete(inptrobs1 )
   else 
     B = inptrobs->$varsobs(iv)$(0, {clevel(iv)}, :, :)
     Ua_cl =inptrobs->U(0, {clevel(iv)}, :,:)
     Va_cl =inptrobs->V(0, {clevel(iv)}, :,:)
   end if
   B = B * cscaleobs(iv)
     
     if(obsdataset(iv).eq."CCCM") then
       B_xy = B
     else  
       B_xy = B
     end if

   dsizes = dimsizes(B_xy)
   fieldobs_xy = new((/dsizes(0), dsizes(1)/), float)
   fieldobs_xy!0 = "lat"
   fieldobs_xy!1 = "lon"
   fieldobs_xy&lon = B_xy&lon
   fieldobs_xy&lat = B_xy&lat
   fieldobs_xy&lon@units = "degrees_east"
   fieldobs_xy&lat@units = "degrees_north"
     
   fieldobs_xy(:, :) = ((/B_xy/))
 
   delete([/inptrobs,B,B_xy/])

;************************************************
; create plot
;************************************************
  plotname = "./pic/Horizontal_"+vars(iv)+"_"+cseason+clevel(ire)
   wks = gsn_open_wks(ptype,plotname)              ; open a ps file

  gsn_define_colormap(wks,"MPL_PuOr")
  plot = new(ncases+1, graphic)

  res                      = True               ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame    = False                             ; don't advance yet
  res@sfXArray     = clon                           ; cell center (lon)
  res@sfYArray     = clat                           ; cell center (lat)
  res@cnFillOn     = True                           ; do color fill
  res@cnFillMode   = "RasterFill"
  res@cnLinesOn    = False
 
  res@gsnMaximize = True
  res@mpFillOn     = True
  res@mpLimitMode  = "LatLon"
  res@mpMaxLonF    = max(clon);+1                     ; specify the plot domain
  res@mpMinLonF    = min(clon);-1                     ;                         
  res@mpMinLatF    = min(clat);-1                     ;
  res@mpMaxLatF    = max(clat);+1                     ;
  res@mpCenterLonF = 180
      res@tiMainString                   = ""
      res@tiMainFont                     = _Font
      res@tiMainFontHeightF              = 0.016

      res@tiXAxisString                  = ""
      res@tiXAxisFont                    = _Font
      res@tiXAxisFontHeightF             = 0.016

      res@tiYAxisString                  = ""
      res@tiYAxisFont                    = _Font
      res@tiYAxisOffsetXF                = 0.0
      res@tiYAxisFontHeightF             = 0.018

 ;     res@gsnCenterStringOrthogonalPosF  = 0.04              ; Position
      res@txFuncCode    = "~"

      res@tmXBLabelFont = _Font
      res@tmYLLabelFont = _Font
      res@gsnStringFont = _Font
      res@tiYAxisFont   = _Font

      res@gsnStringFont                  = _Font
      res@gsnStringFontHeightF           = 0.04
      res@gsnRightString                 = "";"Cloud Fraction"
      res@gsnScalarContour     = True               ; contours desired
      res@cnInfoLabelOn                  = False        ; Label needless
      res@cnFillOn                       = True
      res@cnLinesOn                      = False
      res@cnLineLabelsOn                 = False
      res@lbLabelBarOn                   = False


      res@vcRefMagnitudeF = 5.
      res@vcMinMagnitudeF = 1.
      res@vcRefLengthF    = 0.04
      res@vcRefAnnoOn     = True;False
      res@vcRefAnnoZone   = 3
      res@vcRefAnnoFontHeightF = 0.02
      res@vcRefAnnoString2 =""
      res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
;      res@vcRefAnnoArrowLineColor   = "blue"         ; change ref vector color
;      res@vcRefAnnoArrowUseVecColor = False

      res@vcMinDistanceF  = .05
;      res@vcMinFracLengthF         = .
      res@vcRefAnnoParallelPosF    =  0.997
      res@vcFillArrowsOn           = True
      res@vcLineArrowThicknessF    =  3.0
      res@vcLineArrowHeadMinSizeF   = 0.01
      res@vcLineArrowHeadMaxSizeF   = 0.03
      res@vcGlyphStyle              = "CurlyVector"     ; turn on curley vectors
;      res@vcGlyphStyle              ="Fillarrow"
      res@vcMonoFillArrowFillColor = True;False
      res@vcMonoLineArrowColor     = True;False
      res@vcLineArrowColor          = "green"           ; change vector color
      res@vcFillArrowEdgeColor      ="white"
      res@vcPositionMode            ="ArrowTail"
      res@vcFillArrowHeadInteriorXF =0.1
      res@vcFillArrowWidthF         =0.05           ;default
      res@vcFillArrowMinFracWidthF  =.5
      res@vcFillArrowHeadMinFracXF  =.5
      res@vcFillArrowHeadMinFracYF  =.5
      res@vcFillArrowEdgeThicknessF = 2.0
      res@mpFillOn                   = False

  do im=0, ncases-1 
     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     if(im.eq.ncases) then
       res@gsnLeftString = casenames(im)+" from "+obsdataset(iv)+" "+clevel(ire)
     else
       res@gsnLeftString = casenames(im)+" "+clevel(ire)
     end if
;     res@gsnRightString = "Mean: "+sprintf("%5.2f", global(im,iv,0))+" "+cunits(iv)
;     res@gsnSpreadColors = True
     delete(res@cnLevels)
     res@cnLevels             = cntrs(iv,:)
     res@lbLabelBarOn = False 
     if(mod(im,2).eq.0) then
        res@tiYAxisOn  = True
     else
        res@tiYAxisOn  = False
     end if
     res@tiXAxisOn  = False 
     plot(im) = gsn_csm_contour_map(wks,field_xy(im, :), res)       
     ;plot(im) =gsn_csm_vector_scalar_map_ce(wks,U_cl,V_cl,field_xy(im, :),res)
     ;  wmsetp("col", 2)    ; Draw in red.
     ;  wmsetp("wbs", 0.01)  ; Increase the size of the barb.
     ;  wmbarbmap(wks, lat, lon,  U_cl, V_cl)  ; Plot barbs.

;      gsn_polymarker(wks,plot(im),lons,lats,mkres)

;  frame(wks)
  end do  ; end im

  ; observation 
  delete ( res@sfXArray)
  delete ( res@sfYArray)
  res@gsnLeftString = obsdataset(iv)+" "+clevel(ire)
  delete(res@cnLevels)
  res@cnLevels             = cntrs(iv,:)
;  res@lbLabelBarOn = True 
;  res@lbOrientation        = "vertical"         ; vertical label bars
  res@lbLabelFont          = 12
  res@tiYAxisOn  = True
  res@tiXAxisOn  = True
  res@tiXAxisFont = 12
  plot(ncases) =   gsn_csm_vector_scalar_map_ce(wks,Ua_cl,Va_cl,fieldobs_xy(:, :),res)

  ;*********************************************
  ; panel plot
  ;*********************************************
  pres            = True
  pres@gsnMaximize = True
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.02
  pres@gsnPanelLabelBar        = True
  pres@pmLabelBarWidthF        = 0.8
  gsn_panel(wks, plot, (/(ncases+1),2/), pres)

  delete(field_xy)
  delete(fieldobs_xy)
  end do ; end ire
 end do  ; end iv

end 
