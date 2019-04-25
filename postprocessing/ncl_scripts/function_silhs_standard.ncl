procedure SILHS_terms (ptype, cseason,ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)
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

 vars = (/ "SILHS_CLUBB_PRECIP_FRAC","SILHS_CLUBB_ICE_SS_FRAC","HR","AWNC","AWNI","ADSNOW","ANSNOW","AQSNOW","RHW","QRS","QRL","AQRAIN"/)

 nvars = dimsizes(vars)
 cunits = (/"%", "mba/day","g/kg","g/kg","K", "%", "mba/day", "K", "g/kg", "m/s", "m/s","K","m"/)
 cscale = (/1, 1, 1, 1, 1., 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
 cntrs = new ((/nvars,11/),float)
        ; observationl data set

;=================================================================================================================
; Read sites info
;=================================================================================================================



; open files
 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+"_"+cseason+"_climo.nc"
    print((/infile/))
    inptrs(im) = addfile(infile,"r")       ; pointer to file1
    ncdfs(im)     = addfile("./data/"+cases(im)+"_site_location.nc" ,"r")  ; open output netCDF file
 end do

 inptr = inptrs(0)
 lat = inptr->lat
 nlat = dimsizes(lat)
 lon = inptr->lon
 nlon = dimsizes(lon)
 lev = inptr->lev
 nlev = dimsizes(lev)
 delete(inptr)

do ire = 0, nsite-1

 do iv=0, nvars-1
   do im=0, ncases-1
     ncdf = ncdfs(im)
     n=ncdf->n
     idx_cols=ncdf->ind_cols
    print("num="+n(ire))
    delete (ncdf)
     inptr = inptrs(im)
       A = dble2flt(inptr->SILHS_CLUBB_PRECIP_FRAC(0, :, idx_cols(ire,0)-1 ))
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

cunits(iv)= inptr->$vars(iv)$@units
     if(im.eq.0) then
      dsizes = dimsizes(A)
      field_xz = new((/ncases, dsizes(0)/), float)
      field_xz!1 = "lev"
      field_xz!0 = "model"
      delete(dsizes)
     end if
     field_xz(im, :) = ((/A/))
     field_xz&lev=lev
     delete([/inptr,A,idx_cols,n/])
   end do  ; end im

;************************************************
; create plot
;************************************************
  if (iv.eq.0) then
   plotname = "./pic/CLUBB_SILHS_"+lons(ire)+"E_"+lats(ire)+"N_"+cseason
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
 ;  res@vpYF             = 0.30
   res@tmYLLabelFont  = 12
   res@tmXBLabelFont  = 12
   res@tmXBLabelFontHeightF = 0.01
   res@tmXBLabelFontThicknessF = 1.0

   res@tmXBLabelAngleF = 45
   delete(res@tmXBValues)
   delete(res@tmXBLabels)
   delete(res@xyLineThicknesses)
   delete(res@xyLineColors)
   delete(res@xyDashPatterns)
   delete(res@xyExplicitLegendLabels)

   res@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0, 20, 2.0/)               ; make 2nd lines thicker
   res@xyLineColors      = ispan (2,20,1) 
   res@xyDashPatterns = (/0, 1, 0, 1, 0, 1/)
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

  plot(iv) = gsn_csm_xy (wks, field_xz(0:ncases-1, :), field_xz&lev, res)

  ; observation 
  resobs = res
  delete(resobs@xyLineThicknesses)
  delete(resobs@xyLineColors)
  delete(resobs@xyDashPatterns)
  resobs@xyLineThicknesses = (/3.0/)               ; make 2nd lines thicker
  resobs@xyLineColors      = (/"black"/)          ; change line color
  resobs@xyDashPatterns = (/0/)
  delete(resobs@xyExplicitLegendLabels)
  resobs@xyExplicitLegendLabels= "Obs"
  resobs@pmLegendParallelPosF   = 0.6                  ; move units right
  resobs@pmLegendOrthogonalPosF = -.75                ; more neg = down
;  plotobs = gsn_csm_xy (wks, fieldobs_xz(:, 0), fieldobs_xz&lev, resobs)

;  overlay (plot(iv), plotobs)

  if(iv.eq.nvars-1) then
    ;*********************************************
    ; panel plot
    ;*********************************************
    pres            = True
;    pres@gsnMaximize = True
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
     

