procedure clubb_bgt (ptype, cseason,ncases, cases, casenames, nsite, lats, lons, filepath, filepathobs)
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

 vars = (/ "wp2","wp3","wprtp","wpthlp","rtp2","thlp2","up2","vp2","rtpthlp"/)
 varsobs = (/ "CLOUD", "OMEGA","SHUM","CLWC_ISBL", "THATA","RELHUM"/)
 nvars = dimsizes(vars)
 cunits = (/"%", "mba/day","g/kg","g/kg","K", "%", "mba/day", "K", "g/kg", "m/s", "m/s","K","m"/)
 cscale = (/1, 1, 1, 1, 1., 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
 cscaleobs =  (/100., 100/86400., 1., 1000, 1., 1., 1, 1,1,1/)
 cntrs = new ((/nvars,11/),float)
        ; observationl data set
 obsdataset =(/"CCCM", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI", "ERAI","ERAI","ERAI"/)

;=================================================================================================================
; Read sites info
;=================================================================================================================



; open files
 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+"_"+cseason+"_budget_climo.nc"
    print((/infile/))
    inptrs(im) = addfile(infile,"r")       ; pointer to file1
    ncdfs(im)     = addfile("./data/"+cases(im)+"_site_location.nc" ,"r")  ; open output netCDF file
 end do

 inptr = inptrs(0)
 lat = inptr->lat
 nlat = dimsizes(lat)
 lon = inptr->lon
 nlon = dimsizes(lon)
 lev = inptr->ilev
 nlev = dimsizes(lev)
; gw = inptr->gw
 delete(inptr)
do ire = 0, nsite-1
   do im=0, ncases-1
     plotname = "./pic/Budgets_"+cases(im)+lons(ire)+"E_"+lats(ire)+"N_"+cseason


       do iv=0, nvars-1
    if (vars(iv) .eq. "rtp2" .or. vars(iv) .eq. "thlp2") then 
    budget_ends = (/"_bt", "_ma", "_ta", "_tp", "_dp1", "_dp2", "_cl", "_pd", "_sf", "_forcing"/)
    nterms = dimsizes (budget_ends)
    end if
    if (vars(iv) .eq. "wprtp") then
    budget_ends = (/"_bt", "_ma", "_ta", "_tp", "_ac","_bp","_pr1","_pr2", "_pr3","_dp1","_mfl", "_cl", "_sicl","_pd", "_forcing"/)
    nterms = dimsizes (budget_ends)
    end if

    if (vars(iv) .eq. "wpthlp") then
    budget_ends = (/"_bt", "_ma", "_ta", "_tp", "_ac","_bp","_pr1","_pr2", "_pr3","_dp1","_mfl", "_cl", "_sicl", "_forcing"/)
    nterms = dimsizes (budget_ends)
    end if


    if (vars(iv) .eq. "rtpthlp") then
    budget_ends = (/"_bt", "_ma", "_ta", "_tp1","_tp2","_dp1","_dp2", "_cl", "_sf", "_forcing"/)
    nterms = dimsizes (budget_ends)
    end if

    if (vars(iv) .eq. "wp2") then
    budget_ends = (/"_bt", "_ma", "_ta", "_ac","_bp","_pr1","_pr2", "_pr3","_dp1","_dp2", "_cl", "_pd", "_sf"/)
    nterms = dimsizes (budget_ends)
    end if

    if (vars(iv) .eq. "wp3") then
    budget_ends = (/"_bt", "_ma", "_ta", "_tp", "_ac","_bp1","_pr1","_pr2","_dp1", "_cl"/)
    nterms = dimsizes (budget_ends)
    end if


    if (vars(iv) .eq. "up2" .or. vars(iv) .eq. "vp2") then
    budget_ends = (/"_bt", "_ma", "_ta", "_tp", "_dp1", "_dp2","_pr1","_pr2" ,"_cl", "_pd", "_sf"/)
    nterms = dimsizes (budget_ends)
    end if

     ncdf = ncdfs(im)
     n=ncdf->n
     idx_cols=ncdf->ind_cols
    print("num="+n(ire))
    delete (ncdf)
     inptr = inptrs(im)
       A = dble2flt(inptr->wp2(0, :, idx_cols(ire,0)-1 ))
       A = 0.0
       dsizes = dimsizes(A)
       A_field= new((/nterms, dsizes(0)/), float)
       A_field!1 = "lev"
       A_field!0 = "bgt"
       A_field&bgt= budget_ends 
       A_field = 0.0
       A_field&lev=lev

     do in = 0, nterms -1
      do subc = 0, n(ire)-1
        vars_bgt= vars(iv)+budget_ends(in)
        A_field(in,:) =(/ A_field(in,:)+dble2flt((inptr->$vars_bgt$(0,:, idx_cols(ire,n(subc)-1)-1 )))/n(ire)/)
      end do
     cunits(iv)= inptr->$vars_bgt$@units
     end do

;************************************************
; create plot
;************************************************
  if (iv.eq.0) then
   wks = gsn_open_wks(ptype,plotname)              ; open a ps file
   gsn_define_colormap(wks,"radar_1")
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
   res@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.,2.,2.,2.,2,2,2,2,2,2,2/)  
   res@xyLineColors      = ispan(2,25,2)
   res@xyDashPatterns = ispan(0,24,1)
  end if ; iv

    res@pmLegendDisplayMode    = "ALWAYS"
    res@xyExplicitLegendLabels =  budget_ends(:) 
    res@pmLegendSide           = "top"                 ; Change location of
    res@pmLegendParallelPosF   = 0.75                  ; move units right
    res@pmLegendOrthogonalPosF = -0.65                  ; more neg = down
    res@pmLegendWidthF         = 0.10                 ; Change width and
    res@pmLegendHeightF        = 0.15                  ; height of legend.
    res@lgLabelFontHeightF     = .01                   ; change font height
    res@lgLabelFontThicknessF  = 1.
    res@lgPerimOn              = False

  res@gsnLeftString     = vars(iv)
  res@gsnRightString  = cunits(iv)
  res@trYReverse        = True

  plot(iv) = gsn_csm_xy (wks, A_field(:, :), A_field&lev, res)
 delete(res@xyExplicitLegendLabels)
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
     delete([/inptr,A,idx_cols,n,A,A_field,budget_ends/])
  end do ; end iv


  end do ; im
end do ; end ire

end
     

