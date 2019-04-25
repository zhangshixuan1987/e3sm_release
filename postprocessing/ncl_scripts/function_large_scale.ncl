procedure large_scale_prof (ptype,cseason, ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)
; ncases, the number of models
; cases, the name of models
; casename, the name of cases
; filepath, model output filepath
; filepathobs, filepath for observational data

begin

 inptrs = new((/ncases/),  "file") ; pointers to files
 ncdfs   = inptrs

; region information 
 nregions = nsite   ;  6 regions
; fields to plot

 vars    = (/ "CLOUD"  , "OMEGA",   "Q",   "CLDLIQ", "THATA","RELHUM","U",   "CLDICE","T"/)
 varsobs = (/ "CC_ISBL", "OMEGA","SHUM","CLWC_ISBL", "THATA","RELHUM","U","CIWC_ISBL","T"/)
 nvars = dimsizes(vars)
 cunits    = (/"%","mba/day","g/kg","g/kg","K", "%", "m/s", "g/kg", "m/s", "m/s","K","m"/)
 cscale    = (/100,      864,  1000, 1000 , 1.,   1,     1,   1000,     1,1,1,1,1,1,1/)
 cscaleobs = (/100,        1,     1, 1000 , 1.,   1,     1,   1000,     1,1,1,1,1,1,1/)
 cntrs = new ((/nvars,11/),float)
        ; observationl data set
 obsdataset =(/"ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI","ERAI","ERAI"/)

;=================================================================================================================
; Read sites info
;=================================================================================================================



; open files
 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+"_"+cseason+"_climo.nc"
;    print((/infile/))
    inptrs(im) = addfile(infile,"r")       ; pointer to file1
    ncdfs(im)     = addfile("./data/"+cases(im)+"_site_location.nc" ,"r")  ; open output netCDF file
 end do
;=================================================================================================
 inptr = inptrs(0)
 lat   = inptr->lat
 nlat  = dimsizes(lat)
 lon   = inptr->lon
 nlon  = dimsizes(lon)
 lev   = inptr->lev
 nlev  = dimsizes(lev)
 delete(inptr)

do ire = 0, nsite-1

 do iv=0, nvars-1
   do im=0, ncases-1
     ncdf = ncdfs(im)
     n=ncdf->n
     idx_cols=ncdf->ind_cols
     inptr = inptrs(im)
    if (iv.eq.0) then
    print("CASE="+cases(im)+";num="+n(ire))
    end if
       A = dble2flt(inptr->T(0, :, idx_cols(ire,0)-1 ))
       A = 0.0
    do subc = 0, n(ire)-1
     if(vars(iv).ne."THATA") then
       A =(/ A+dble2flt((inptr->$vars(iv)$(0,:, idx_cols(ire,n(subc)-1)-1 )))/n(ire)/)
     else
       A =(/ A+dble2flt((inptr->T(0,:, idx_cols(ire,n(subc)-1)-1 )))/n(ire)/)
       dimsizes1 = dimsizes(A)
       hyam = inptr->hyam
       hybm = inptr->hybm
       ps = dble2flt(inptr->PS(0, idx_cols(ire,n(subc)-1)-1  ))
       p0 = dble2flt(inptr->P0)
       pre = new((/dimsizes1(0)/), "float")
       do ilev=0, dimsizes1(0)-1
         pre(ilev ) = doubletofloat(hyam(ilev))*doubletofloat(p0) + doubletofloat(hybm(ilev)) * ps( :)
         A(ilev) = A(ilev) * (100000/pre(ilev))^(0.286)
       end do
       delete([/dimsizes1,hyam,hybm,ps,p0,pre/])
     end if
    end do
     A = A * cscale(iv)


     if(im.eq.0) then
      dsizes = dimsizes(A)
      field_xz = new((/ncases, dsizes(0)/), float)
      field_xz!1 = "lev"
      field_xz!0 = "model"
      delete(dsizes)
     end if
     field_xz(im, :) = ((/A/))
     field_xz&lev=lev
     delete([/inptr,A,ncdf,idx_cols,n/])
   end do  ; end im
;  Observational data
   if(obsdataset(iv).eq."CCCM") then
     if(cseason.eq."ANN") then
       fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-"+cseason+".nc"
     else
       fileobs = "/Users/guoz/databank/CLD/CCCm/cccm_cloudfraction_2007-2010-"+cseason+".nc"
     end if
   else
     fileobs = filepathobs + obsdataset(iv)+"_"+cseason+"_climo.nc"
   end if
   inptrobs = addfile(fileobs, "r")
   if(obsdataset(iv).eq."CCCM") then
     B = inptrobs->$varsobs(iv)$(:, {lats(ire)}, {lons(ire)})
   else
     if(varsobs(iv).ne."THATA") then
       B = inptrobs->$varsobs(iv)$(0, :, {lats(ire)}, {lons(ire)})
     else
       B = inptrobs->T(0, :, {lats(ire)}, {lons(ire)})
       pre = inptrobs->lev
       dimsizes1 = dimsizes(pre)
       do ilev=0, dimsizes1-1
         B(ilev) = B(ilev) * (1000/pre(ilev))^0.286
       end do
       delete([/pre,dimsizes1/])
     end if
     Z3obs = inptrobs->Z3(0, :, {lats(ire)}, {lons(ire)})
   end if
   B = B * cscaleobs(iv)

   dsizes = dimsizes(B)
   fieldobs_xz = new((/dsizes(0)/), float)
   fieldobs_xz!0 = "lev"
     fieldobs_xz&lev = B&lev
   fieldobs_xz(:) = ((/B/))
   if(obsdataset(iv).ne."CCCM") then
     fieldobs_xz&lev = B&lev
     delete(Z3obs)
   end if

   delete([/inptrobs,B/])


;************************************************
; create plot
;************************************************
  if (iv.eq.0) then
   plotname = "./pic/Fields_"+lons(ire)+"E_"+lats(ire)+"N_"+cseason
   wks = gsn_open_wks(ptype,plotname)              ; open a ps file
   gsn_define_colormap(wks,"default")
   plot = new(nvars, graphic)

   res                      = True               ; plot mods desired
   res@gsnDraw              = False
   res@gsnFrame             = False

   res@lgLabelFontHeightF     = .012                   ; change font height
   res@lgPerimOn              = False                 ; no box around
   res@vpWidthF         = 0.30                      ; set width and height
   res@vpHeightF        = 0.30
   res@txFontHeightF   = .02
 ;  res@vpXF             = 0.04
  ; res@vpYF             = 0.30
   res@tmYLLabelFont  = 12
   res@tmXBLabelFont  = 12
   res@tmXBLabelFontHeightF = 0.01
   res@tmXBLabelFontThicknessF = 1.0
   res@tmXBLabelAngleF = 45

   res@xyLineThicknessF = 2.0;(/2.0, 2.0, 2.0, 2.0, 2.0, 2.0/)               ; make 2nd lines thicker
   res@xyLineColors      = ispan(2,20,1) 
   res@xyDashPatterns    = ispan(0,24,1)
  end if ; iv

  if(iv.eq.nvars-1) then
    res@pmLegendDisplayMode    = "ALWAYS"
    res@xyExplicitLegendLabels = casenames(0:ncases-1)
    res@pmLegendSide           = "top"                 ; Change location of
    res@pmLegendParallelPosF   = 0.6                  ; move units right
    res@pmLegendOrthogonalPosF = -0.5                  ; more neg = down
    res@pmLegendWidthF         = 0.10                 ; Change width and
    res@pmLegendHeightF        = 0.10                  ; height of legend.
    res@lgLabelFontHeightF     = .02                   ; change font height
    res@lgLabelFontThicknessF  = 1.5
    res@lgPerimOn              = False
  else
    res@pmLegendDisplayMode    = "NEVER"
  end if

  res@gsnLeftString     = vars(iv)
  res@gsnRightString  = cunits(iv)
  res@trYReverse        = True

  if(vars(iv).eq."THATA") then
   res@trXMinF = 280.
   res@trXMaxF = 400.
  else if(vars(iv).eq."CLOUD" .or. vars(iv).eq."RELHUM") then
   res@trXMinF = 0.
   res@trXMaxF = 100.
  else if(vars(iv).eq."T") then
   res@trXMinF = 180
   res@trXMaxF = 300
  else if(vars(iv).eq."U") then
   res@trXMinF = -40
   res@trXMaxF = 40

  else
   res@trXMinF = min((/min(field_xz(0:ncases-1, :)),min(fieldobs_xz(:))/))
   res@trXMaxF = max((/max(field_xz(0:ncases-1, :)),max(fieldobs_xz(:))/))
  end if
  end if
  end if
  end if

  plot(iv) = gsn_csm_xy (wks, field_xz(0:ncases-1, :), field_xz&lev, res)
   delete(res@trXMinF)
   delete(res@trXMaxF)


  ; observation 
  resobs = res
  delete(resobs@xyLineColors)
  delete(resobs@xyDashPatterns)
  resobs@xyLineThicknessF = 3.0               ; make 2nd lines thicker
  resobs@xyLineColors      = (/"black"/)          ; change line color
  resobs@xyDashPatterns = (/0/)
  delete(resobs@xyExplicitLegendLabels)
  resobs@xyExplicitLegendLabels= "ERAI"
  resobs@pmLegendParallelPosF   = 0.6                  ; move units right
  resobs@pmLegendOrthogonalPosF = -.75                ; more neg = down
  plotobs = gsn_csm_xy (wks, fieldobs_xz(:), fieldobs_xz&lev, resobs)

  overlay (plot(iv), plotobs)

  if(iv.eq.nvars-1) then
    ;*********************************************
    ; panel plot
    ;*********************************************
    pres            = True
;    pres@gsnMaximize = True/
    pres@txString   = lons(ire)+"E,"+lats(ire)+"N"
    pres@txFont = 12
    pres@gsnPanelYWhiteSpacePercent = 5
    pres@gsnPanelXWhiteSpacePercent = 5
    pres@gsnPanelBottom = 0.18
    gsn_panel(wks, plot, (/3,nvars/3/), pres)
  end if

  delete(field_xz)
  end do ; end iv

end do ; end ire

end
     

