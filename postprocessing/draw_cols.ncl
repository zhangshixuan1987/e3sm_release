;========================================================================
; This is a diagnostic package for E3SM and CAM with SE dyn core.
; This package mainly draws vertical profiles of standard outputs, CLUBB terms,
; CLUBB budget temrs, skewness terms, etc. It also provides an opportunity
; to compare multiple experiments, even if they don't share a same horizontal resolution. 

; Horizontal_* spatial distributions at selected levels, e.g. 200, 500, 850 hPa. 
; Fields_locations_season for large scale features, such as U,V,T,Q
; CLUBB_SILHS_locations_season for SILHS variables
; Budgets_$cases_locations_season for detailed budgets
; CLUBB_skewfunctions_* for plots of skewness functions
;
; Please connect guozhun@uwm.edu for questions.
;------------------------------------------------------------------------

load  "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load  "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load  "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load  "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load  "./ncl_scripts/function_pick_out.ncl"
load  "./ncl_scripts/function_cal_mean.ncl"
load  "./ncl_scripts/function_large_scale.ncl"
load  "./ncl_scripts/function_clubb_standard.ncl"
load  "./ncl_scripts/function_clubb_budget.ncl"
load  "./ncl_scripts/function_clubb_skew.ncl" 
load  "./ncl_scripts/function_silhs_standard.ncl"
load  "./ncl_scripts/function_plots_hoz.ncl"
;------------------------------------------------------------------------
begin

;------------------------------------------------------------------------
; general setting
; cases' info, please change here.
filepath    = (/"/global/cscratch1/sd/zhun/csmruns/",\
                "/global/cscratch1/sd/zhun/csmruns/",\
                "/global/cscratch1/sd/zhun/csmruns/",\
                "/global/cscratch1/sd/zhun/csmruns/" \
              /)
cases       = (/"Reg_edison_FAMIPC5_ne30_ne30_CLUBB_SILHS"\
               ,"Reg2_edison_FAMIPC5_ne30_ne30_CLUBB_SILHS"\
               ,"Full09_edison_FAMIPC5_ne30_ne30_CLUBB_SILHS"\
               ,"TEST_edison_FAMIPC5_ne30_ne30_CLUBB_SILHS_L72_upgrade_cam_silhs_take2"\
              /)
; Give a short name for your experiment which will appears on plots
casenames   = (/"RR_NR_0.99_in","RR_NR_0.9_bel","RR_NR_0.9_both","BMG"/)

; Observation Data
filepathobs = (/"/global/project/projectdirs/m2689/zhun/amwg/obs_data_20140804/"/) 

;------------------------------------------------------------------------
; Setting of plots.
ptype="pdf"   ; eps, pdf, ps... are supported by this package
cseason="JJA" ; Seasons, or others 

;------------------------------------------------------------------------
; set to true when you first use it.
calmean="ftrue" ; make mean state
findout="ftrue" ; pick out the locations of your sites

;------------------------------------------------------------------------
; Please give the lat and lon of sites here. 
; If below list are changed, please switch findout="true"

; sites    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
;lats = (/  20,  27, -20, -20,  -5,  -1, -60,  60,   2,   9,  56,  76,  45,   0,  10/)
;lons = (/ 190, 240, 275, 285, 355, 259, 340, 180, 140, 229, 311, 320, 180, 295,  90/)
lats = (/  20,  27, -20, -20,  -5,  -1,  60,   2,   9,  56,  45,   0,  10/)
lons = (/ 190, 240, 275, 285, 355, 259,  180, 140, 229, 311,  180, 295,  90/)

area  = .5 
; Note, this para helps to find out the 'ncol' within 
; lats - area < lat(ncol) < lons + area .and. lons- area < lon(ncol) < lons + area

;------------------------------------------------------------------------
; no need to change
ncases = dimsizes(cases) 
nsite  = dimsizes(lats) 
print("Here we go!")
system("mkdir ./pic/")

;------------------------------------------------------------------------
; Step 1 
; For generate the mean state files, *_climo.nc (standard outputs) and 
; *_budget_climo.nc (budgets of CLUBB).  Just needed at the first run
if (calmean .eq. "true") then
cal_mean(ncases, cases, nsite, lats, lons, area, filepath)
end if

;------------------------------------------------------------------------
; Step 2 find out the locations of all sites from raw SE outputs
; Just need at the first run, info of sites will be stored in ./data/

if (findout .eq. "true" ) then
pick_out(ncases, cases, nsite, lats, lons, area, filepath)
end if

;------------------------------------------------------------------------
; Step 3 read $cases_site_location.nc, in which index of col ;
; and their number can be find.

; 3.1 For AMWG standard outputs, call function_plots_hoz.ncl 
plots_hoz_multiruns(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)

; 3.2 For AMWG standard outputs, profiles, call function_large_scale.ncl
large_scale_prof(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)

; 3.3 For CLUBB standard outputs, call function_clubb_standard.ncl
clubb_stat(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)

; 3.4 For SILHS and microphys, call function_silhs_standard.ncl
SILHS_terms(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)

; 3.5 CLUBB's budgets, call function_clubb_budget.ncl
clubb_bgt(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)

; 3.6 Skewness terms in h2 files. function_clubb_skew.ncl 
clubb_terms(ptype,cseason,ncases,cases, casenames, nsite, lats,lons,filepath, filepathobs)


end

