#!/bin/csh
date

set echo verbose


set fetch_code    = 0   # 0 = No, >0 = Yes
set compile_model = 1   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

####################################################################
# First, you should checkout the E3SM code following instructions
# here: 
#    https://kaizhangpnl.github.io/start.html#quick-start-guide
####################################################################
setenv CCSMTAG CLUBB_SILHS
setenv CCSMROOT $PWD

####################################################################
# Machine, compset, PE layout etc.
####################################################################

#setenv COMPSET FC5AV1C-04P2
setenv COMPSET FAMIPC5
setenv RESOLUTION ne30_ne30
setenv MACH      edison
setenv PTMP      /global/cscratch1/sd/$user/bld
#setenv ntasks 960
setenv ntasks 864
setenv nthrds 1

# Configuration parameters
setenv NUMSC 4
setenv MGVER 2 # Currently "1" and "2" are allowed

setenv CASESRC   nothing
setenv MYSRC     ${CCSMROOT}/mods_$CASESRC

# The run case
setenv CASE     TEST_${MACH}_${COMPSET}_${RESOLUTION}_${CCSMTAG}
# The executable case; if use single executable for multiple runs
setenv COMCASE  TEST_${MACH}_${COMPSET}_${RESOLUTION}_${CCSMTAG}

setenv CASEROOT  ${CCSMROOT}/cases/$CASE
setenv RUNDIR    /global/cscratch1/sd/$user/csmruns/$CASE

####################################################################
# Compile model
####################################################################
if ($compile_model > 0) then

   rm -rf $CASEROOT
   cd  $CCSMROOT/cime/scripts

   ./create_newcase --case $CASEROOT --mach $MACH --project m2689 \
                    --res $RESOLUTION --compset $COMPSET

   #====================================================================
   # set up case
   #====================================================================

   ###./create_newcase -list grids

   cd $CASEROOT

   ./xmlchange -file env_run.xml   -id RUNDIR  -val $RUNDIR
   ./xmlchange -file env_build.xml -id EXEROOT -val $PTMP/$COMCASE/bld/

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

   ./case.setup

   #====================================================================
   # my mods of source code
   #====================================================================
   #cd $CASEROOT
   #
   #ln -s ${MYSRC}/* SourceMods/src.cam    # put your mods in here
   #
   #./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append -val ' -cosp'

   ./xmlchange CAM_CONFIG_OPTS="-dyn se -phys cam5 -clubb_sgs -rad rrtmg -chem trop_mam3 -silent -nlev 72 -microphys mg$MGVER -psubcols $NUMSC -cppdefs '-DUWM_MISC -DSILHS'"
   #./xmlchange -file env_build.xml CAM_CONFIG_OPTS="-phys cam5 -clubb_sgs -microphys mg2 -chem linoz_mam4_resus_mom -rain_evap_to_coarse_aero -nlev 72 -bc_dep_to_snow_updates -psubcols 8 -cppdefs '-DUWM_MISC -DSILHS'"
   #./xmlchange -file env_build.xml CAM_CONFIG_OPTS="-phys cam5 -clubb_sgs -microphys mg2 -chem linoz_mam4_resus_mom -rain_evap_to_coarse_aero -bc_dep_to_snow_updates -psubcols 8 -cppdefs '-DUWM_MISC -DSILHS'"

   # Build the model
   cd $CASEROOT

   ./case.build

   ###./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then

#------------------
## set environment
#------------------

cd $CASEROOT

#./xmlchange  -file env_batch.xml  -id  JOB_WALLCLOCK_TIME   -val '00:30:00'
#./xmlchange  -file env_batch.xml  -id  JOB_QUEUE   -val 'debug'
./xmlchange  -file env_batch.xml  -id  JOB_WALLCLOCK_TIME   -val '18:00:00'
#./xmlchange  -file env_batch.xml  -id  JOB_WALLCLOCK_TIME   -val '01:30:00'
./xmlchange  -file env_batch.xml  -id  JOB_QUEUE   -val 'regular'

./xmlchange -file env_run.xml -id SAVE_TIMING_DIR -val '/global/homes/g/griffinb/'

./xmlchange  -file env_run.xml  -id  RUN_STARTDATE   -val '1979-01-01'
./xmlchange  -file env_run.xml  -id  RESUBMIT        -val '0'
##./xmlchange  -file env_run.xml  -id  CONTINUE_RUN    -val 'TRUE'
./xmlchange  -file env_run.xml  -id  STOP_N          -val '14'
./xmlchange  -file env_run.xml  -id  STOP_OPTION     -val 'nmonths'
./xmlchange  -file env_run.xml  -id  REST_N          -val '14'
./xmlchange  -file env_run.xml  -id  REST_OPTION     -val 'nmonths'
./xmlchange  -file env_run.xml  -id  DOUT_S          -val 'FALSE'

# A list of CLUBB variables
#clubb_vars_zt_list="'thlm', 'thvm', 'rtm', 'rcm', 'rvm', 'um', 'vm', 'um_ref','vm_ref','ug', 'vg', 'cloud_frac', 'cloud_cover', 'rcm_in_layer', 'rcm_in_cloud', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt', 'Lscale', 'Lscale_pert_1', 'Lscale_pert_2', 'T_in_K', 'rel_humidity', 'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down', 'tau_zt', 'Kh_zt', 'wp2thvp', 'wp2rcp', 'wprtpthlp', 'sigma_sqd_w_zt', 'rho', 'radht', 'radht_LW', 'radht_SW', 'Ncm', 'Nc_in_cloud', 'Nc_activated', 'snowslope', 'sed_rcm', 'rsat', 'rsati', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm', 'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_forcing', 'rtm_sdmp','rtm_mc', 'rtm_pd', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_bt', 'thlm_ma', 'thlm_ta', 'thlm_mfl', 'thlm_tacl', 'thlm_cl', 'thlm_forcing', 'thlm_sdmp','thlm_mc', 'thlm_old', 'thlm_without_ta', 'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl', 'rtm_old', 'rtm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl', 'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'vm_bt', 'vm_ma', 'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp', 'wp3_ac', 'wp3_bp1', 'wp3_bp2', 'wp3_pr1', 'wp3_pr2', 'wp3_dp1', 'wp3_cl', 'mixt_frac', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2', 'cloud_frac_1', 'cloud_frac_2', 'a3_coef_zt', 'wp3_on_wp2_zt', 'chi_1', 'chi_2', 'stdev_chi_1', 'stdev_chi_2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2', 'corr_chi_eta_1', 'corr_chi_eta_2', 'rrtthl', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'precip_frac', 'precip_frac_1', 'precip_frac_2', 'Ncnm', 'wp2_zt', 'thlp2_zt', 'wpthlp_zt', 'wprtp_zt', 'rtp2_zt', 'rtpthlp_zt', 'up2_zt', 'vp2_zt', 'upwp_zt', 'vpwp_zt', 'C11_Skw_fnc'"
#clubb_vars_zm_list="'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp4', 'up2', 'vp2', 'wpthvp', 'rtpthvp', 'thlpthvp', 'tau_zm', 'Kh_zm', 'wprcp', 'rc_coef', 'wm_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'rho_zm', 'sigma_sqd_w', 'Skw_velocity', 'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C7_Skw_fnc', 'C1_Skw_fnc', 'a3_coef', 'wp3_on_wp2', 'rcm_zm', 'rtm_zm', 'thlm_zm', 'cloud_frac_zm', 'rho_ds_zm', 'thv_ds_zm', 'em', 'mean_w_up', 'mean_w_down', 'shear', 'wp3_zm', 'Frad', 'Frad_LW', 'Frad_SW', 'Frad_LW_up', 'Frad_SW_up', 'Frad_LW_down', 'Frad_SW_down', 'Fprec', 'Fcsed', 'wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac', 'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_dp1', 'wp2_dp2', 'wp2_cl', 'wp2_pd', 'wp2_sf', 'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_cl', 'vp2_pd', 'vp2_sf', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1', 'up2_pr2', 'up2_cl', 'up2_pd', 'up2_sf', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp', 'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl', 'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma', 'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3', 'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'rtp2_bt', 'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf', 'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1', 'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt', 'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2', 'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'wpthlp_entermfl', 'wpthlp_exit_mfl', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max', 'wprtp_mfl_min', 'wprtp_mfl_max', 'Richardson_num', 'shear_sqd'"

cat >> user_nl_cam << EOF
&camexp

use_gw_convect = .false.

dtime = 1800
nhtfrq = 0,-24,0
mfilt = 1,5000,5000
ndens = 1,1,1,1,1,1
history_budget = .true.
microp_scheme = 'MG'
micro_mg_version = $MGVER
micro_mg_sub_version = 0
micro_mg_num_steps = 1
micro_mg_dcs = 390e-6
micro_mg_berg_eff_factor = 1.0
cldfrc2m_rhmini = 0.8
cldfrc2m_rhmaxi = 1.05

! following parameters overrided by hardwired values in initial implementation. 
! settings below to make it consistent 
clubb_C1 = 1.0D0
clubb_C2rt = 0.2D0
clubb_C2rtthl = 0.2D0
clubb_gamma_coef = 0.24D0
clubb_gamma_coefb = 0.37D0
clubb_c_K10 = 0.3D0
clubb_C14 = 0.5D0
clubb_C8 = 2.5D0

ice_supersat = .true.
macrop_scheme = 'CLUBB_SGS'
eddy_scheme = 'CLUBB_SGS'
shallow_scheme = 'CLUBB_SGS'
deep_scheme = 'off'
subcol_scheme = 'SILHS'
use_subcol_microp = .true.
microp_uniform = .true.
history_amwg = .true.
clubb_do_adv = .false.
clubb_expldiff = .false.
clubb_rainevap_turb = .false.
clubb_cloudtop_cooling = .false.
fincl1 = 'U:A','PS:A','T:A','V:A','OMEGA:A','Z3:A','PRECT:A',
'CLDLIQ:A', 'CLDICE:A', 'LWCF:A', 'SWCF:A', 'FLUT:A',
'TMQ:A', 'PRECC:A', 'PRECL:A', 'CME:A', 'PRODPREC:A',
'EVAPPREC:A','EVAPSNOW:A','ICWMRST:A','ICIMRST:A','PRAO:A',
'PRCO:A','QCSEVAP:A','QISEVAP:A','QVRES:A','CMEIOUT:A','VTRMI:A',
'VTRMC:A','QCSEDTEN:A','QISEDTEN:A','MNUCCCO:A','MNUCCTO:A',
'MNUCCDO:A','MNUCCDOhet:A','MSACWIO:A','PSACWSO:A','BERGSO:A',
'BERGO:A','MELTO:A','HOMOO:A','QCRESO:A','PRCIO:A','PRAIO:A',
'MELTSDT:A','FRZRDT:A','ADRAIN:A','ADSNOW:A','FREQR:A','FREQS:A',
'PE:A','APRL:A','PEFRAC:A','VPRCO:A','VPRAO:A','RACAU:A',
'QIRESO:A','QCRESO:A','PRACSO:A','MPDT:A','MPDQ:A','MPDLIQ:A',
'MPDICE:A','INEGCLPTEND', 'LNEGCLPTEND', 'VNEGCLPTEND',
'QCRAT:A',
'SL', 'Q', 'RHW', 'QRS', 'QRL', 'HR', 'FDL', 'SILHS_CLUBB_PRECIP_FRAC',
'SILHS_CLUBB_ICE_SS_FRAC'
fincl2 = 'CLDTOT', 'CLDST','CDNUMC','CLDLIQ','CLDICE','FLUT',
'LWCF','SWCF','PRECT'
end_restart = .true.
subcol_SILHS_weight = .true.
subcol_SILHS_numsubcol = $NUMSC
subcol_SILHS_corr_file_name = 'arm_97'
subcol_silhs_q_to_micro = .true. ! if .false. gridbox means are used instead of sample points
subcol_silhs_n_to_micro = .true. ! if .false. gridbox means are used instead of sample points
subcol_silhs_use_clear_col = .false.
subcol_SILHS_constrainmn = .false.
subcol_silhs_ncnp2_on_ncnm2 = 0.05,
hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip = 1.0,
hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip = 1.0,
hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip = 1.0,
hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip = 1.0,
hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip = 1.0,
hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip = 1.0
sol_facti_cloud_borne = 1.0D0
dust_emis_fact = 0.3D0
nucleate_ice_subgrid = 1.0
seasalt_emis_scale = 0.6

! CLUBB history!!!
clubb_history = .true.
clubb_rad_history = .false.
EOF

# Comment out if you're not sure it will be ready to commit
./case.submit

endif

