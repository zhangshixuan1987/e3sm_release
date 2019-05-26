module subcol_SILHS
   !---------------------------------------------------------------------------
   ! Purpose:
   !
   ! Implement a subcolumn scheme based on the Subgrid Importance Latin Hypercube 
   ! Sampling (SILHS) functionality of the CLUBB moist turbulence parameterization.
   !
   !---------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8, r4=>shr_kind_r4, i4=>shr_kind_i4
   use physics_types,    only: physics_state, physics_tend, physics_ptend
   use ppgrid,           only: pcols, psubcols, pver, pverp
   use constituents,     only: pcnst, cnst_get_ind
!   use abortutils,       only: endrun
   use shr_sys_mod,      only: endrun => shr_sys_abort
   use cam_logfile,      only: iulog
   use cam_history,      only: addfld, add_default, outfld, horiz_only
#ifdef CLUBB_SGS
#ifdef SILHS
   use clubb_api_module, only: &
        hmp2_ip_on_hmm2_ip_slope_type, &
        hmp2_ip_on_hmm2_ip_intrcpt_type
#endif
#endif
   use physconst,     only: cpair, gravit, latvap, latice, rair

   implicit none

      private
      save

   public :: subcol_register_SILHS  ! 
   public :: subcol_init_SILHS      ! Initialize 
   public :: subcol_gen_SILHS       ! Generate subcolumn fields by calling SILHS 
   public :: subcol_readnl_SILHS    ! SILHS namelist reader
   public :: subcol_ptend_avg_SILHS
#ifdef SILHS
   public :: subcol_SILHS_var_covar_driver, &
             subcol_SILHS_massless_droplet_destroyer, &
             subcol_SILHS_fill_holes_conserv
   private :: Abs_Temp_profile
   private :: StaticEng_profile
   ! Calc subcol mean ! Calc subcol variance
   private :: meansc
   private :: stdsc
   private :: fill_holes_sedimentation
#endif

   !-----
   ! Private module vars
   !-----
   ! Pbuf indicies
   integer :: thlm_idx, num_subcols_idx, pdf_params_idx, rcm_idx, rtm_idx, ice_supersat_idx, &
              alst_idx, cld_idx, qrain_idx, qsnow_idx, nrain_idx, nsnow_idx, ztodt_idx, &
              prec_pcw_idx, snow_pcw_idx, qcsedten_idx, qrsedten_idx, qisedten_idx, &
              qssedten_idx, vtrmc_idx, umr_idx, vtrmi_idx, ums_idx

   logical :: subcol_SILHS_weight  ! if set, sets up weights for averaging subcolumns for SILHS
   integer :: subcol_SILHS_numsubcol ! number of subcolumns for this run
   logical :: docldfracscaling = .false. ! Weight tendencies by cloud fraction

   character(len=256) :: subcol_SILHS_corr_file_path
   character(len=16)  :: subcol_SILHS_corr_file_name

   logical :: subcol_SILHS_q_to_micro, &
              subcol_SILHS_n_to_micro, &
              subcol_SILHS_use_clear_col, &
              subcol_SILHS_meanice, &
              subcol_SILHS_constrainmn

   logical :: subcol_SILHS_var_covar_src
   logical :: subcol_SILHS_destroy_massless_droplets

   real(r8) :: subcol_SILHS_ncnp2_on_ncnm2

   ! There may or may not be a better place to put this.
   real(r8), parameter :: p0_clubb = 100000._r8


!   real(r8) :: subcol_SILHS_c6rt, subcol_SILHS_c7, subcol_SILHS_c8, subcol_SILHS_c11, &
!               subcol_SILHS_c11b, subcol_SILHS_gamma_coef, &
!               subcol_SILHS_mult_coef, subcol_SILHS_mu

   real(r8) :: ztodt  ! model timestep
#ifdef CLUBB_SGS
#ifdef SILHS
    type(hmp2_ip_on_hmm2_ip_slope_type) :: hmp2_ip_on_hmm2_ip_slope    
    type(hmp2_ip_on_hmm2_ip_intrcpt_type) :: hmp2_ip_on_hmm2_ip_intrcpt
#endif
#endif

contains

   subroutine subcol_register_SILHS()

      !--------------------------------
      ! Register fields needed by SILHS in the physics buffer
      ! Currently, most fields needed by SILHS but calculated by CLUBB are registered
      ! by clubb in clubb_intr.F90.
      ! 
      !--------------------------------
#ifdef CLUBB_SGS
#ifdef SILHS
      use physics_buffer,  only: pbuf_add_field, dtype_r8
      use ppgrid,          only: pver, pverp, pcols, psubcols


     ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
      call pbuf_add_field('QRAIN',      'global',dtype_r8,(/pcols,pver/), qrain_idx)
      call pbuf_add_field('QSNOW',      'global',dtype_r8,(/pcols,pver/), qsnow_idx)
      call pbuf_add_field('NRAIN',      'global',dtype_r8,(/pcols,pver/), nrain_idx)
      call pbuf_add_field('NSNOW',      'global',dtype_r8,(/pcols,pver/), nsnow_idx)

#endif
#endif
   end subroutine subcol_register_SILHS


   subroutine subcol_readnl_SILHS(nlfile)
#ifdef CLUBB_SGS
#ifdef SILHS
      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: masterproc, masterprocid, mpicom
      use spmd_utils,      only: mpi_integer, mpi_logical, mpi_character, mpir8
      use clubb_api_module,only: core_rknd
#endif
#endif
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
#ifdef CLUBB_SGS
#ifdef SILHS
      namelist /subcol_SILHS_nl/ subcol_SILHS_weight, &
                                 subcol_SILHS_numsubcol, &
                                 subcol_SILHS_corr_file_path, &
                                 subcol_SILHS_corr_file_name, &
                                 subcol_SILHS_q_to_micro, &
                                 subcol_SILHS_n_to_micro, &
                                 subcol_SILHS_ncnp2_on_ncnm2, &
                                 hmp2_ip_on_hmm2_ip_slope, &
                                 hmp2_ip_on_hmm2_ip_intrcpt, &
                                 subcol_SILHS_meanice, &
                                 subcol_SILHS_use_clear_col, &
                                 subcol_SILHS_constrainmn
!                                 subcol_SILHS_c6rt, subcol_SILHS_c7, &
!                                 subcol_SILHS_c8, subcol_SILHS_c11, subcol_SILHS_c11b, &
!                                 subcol_SILHS_gamma_coef, subcol_SILHS_mult_coef, subcol_SILHS_mu

      !-----------------------------------------------------------------------------
      ! Set defaults
      subcol_SILHS_var_covar_src = .true.  ! TODO: put this in namelist
      subcol_SILHS_destroy_massless_droplets = .true.  ! TODO: put this in namelist

      ! Eric Raut changed a default.
      hmp2_ip_on_hmm2_ip_slope%Ni = 0.0_core_rknd
      hmp2_ip_on_hmm2_ip_intrcpt%Ni = 0.5_core_rknd

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'subcol_SILHS_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, subcol_SILHS_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun('subcol_readnl_SILHS: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpi_bcast(subcol_SILHS_weight,    1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_numsubcol, 1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_path, len(subcol_SILHS_corr_file_path), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_name, len(subcol_SILHS_corr_file_name), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_use_clear_col, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_constrainmn, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_meanice, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_q_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_n_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_var_covar_src,1,mpi_logical,masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_destroy_massless_droplets, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_ncnp2_on_ncnm2, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%rr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%Nr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%ri, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%Ni, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%rs, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_slope%Ns, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%rr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%Nr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%ri, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%Ni, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%rs, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(hmp2_ip_on_hmm2_ip_intrcpt%Ns, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c6rt, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c7, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c8, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11b, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_gamma_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mult_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mu, 1, mpir8, masterprocid, mpicom, ierr)

! SPMD
#endif
! SILHS
#endif
! CLUBB_SGS
#endif
   end subroutine subcol_readnl_SILHS


   subroutine subcol_init_SILHS(pbuf2d)

      !--------------------------------
      ! Read in parameters and initialize SILHS PDF fields.
      ! Set up indexes into Pbuf fields.
      ! Register history outputs.
      !--------------------------------

      use physics_buffer,          only: physics_buffer_desc, pbuf_get_field, &
                                         dtype_r8, pbuf_get_index
      use units,                   only: getunit, freeunit 
#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,        only: core_rknd, &
                                         l_prescribed_avg_deltaz, &
                                         l_diagnose_correlations, &
                                         l_calc_w_corr, l_use_cloud_cover, &
                                         l_fix_w_chi_eta_correlations, l_const_Nc_in_cloud, &
                                         pdf_dim, &
                                         setup_corr_varnce_array_api, &
                                         init_pdf_hydromet_arrays_api, &
                                         Ncnp2_on_Ncnm2, &
                                         set_clubb_debug_level_api

    use silhs_api_module, only :         l_lh_importance_sampling

#endif
#endif

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CLUBB_SGS
#ifdef SILHS

      integer :: iunit = 501 ! Default value, will get iunit from CAM 
      !character(len=*), parameter :: default_corr_case = "arm_97"
      character(len=*), parameter :: &
            cloud_file_ext  = "_corr_array_cloud.in", & ! File extensions for corr files
            below_file_ext  = "_corr_array_below.in"
      character(len=256) :: corr_file_path_cloud, corr_file_path_below

      ! To set up CLUBB hydromet indices
      integer :: &
          hydromet_dim, & ! Number of enabled hydrometeors
          iirr,         & ! Hydrometeor array index for rain water mixing ratio, rr
          iirs,         & ! Hydrometeor array index for snow mixing ratio, rs
          iiri,         & ! Hydrometeor array index for ice mixing ratio, ri
          iirg,         & ! Hydrometeor array index for graupel mixing ratio, rg
          iiNr,         & ! Hydrometeor array index for rain drop concentration, Nr
          iiNs,         & ! Hydrometeor array index for snow concentration, Ns
          iiNi,         & ! Hydrometeor array index for ice concentration, Ni
          iiNg            ! Hydrometeor array index for graupel concentration, Ng


      ! Set CLUBB's debug level
      ! This is called in module clubb_intr; no need to do it here.
!      call set_clubb_debug_level_api( 0 )

      !-------------------------------
      ! CLUBB-SILHS Parameters (global module variables)
      !-------------------------------

      l_fix_w_chi_eta_correlations = .true.
      l_lh_importance_sampling = .true.
      l_diagnose_correlations = .false.
      l_calc_w_corr = .false.
!      l_prescribed_avg_deltaz = .false.
      l_use_cloud_cover = .false.
      l_const_Nc_in_cloud = .true.

      ! Values from the namelist
      docldfracscaling = subcol_SILHS_use_clear_col

      ! Namelist "tuning" or set correlations
      ! KTC Todo: Move these to a tuning "in" file or into the namelist
      ! JHTODO: we might want these on CLUBB's API and ultimatively on a namelist for tuning
!      C6rt = subcol_SILHS_c6rt
!      C7 = subcol_SILHS_c7                                      ! to all ice clouds
!      C8 = subcol_SILHS_c8
!      C11 = subcol_SILHS_c11
!      C11b = subcol_SILHS_c11b
!      gamma_coef = subcol_SILHS_gamma_coef
!      mult_coef = subcol_SILHS_mult_coef
!      mu = subcol_SILHS_mu

      !call set_clubb_debug_level( 0 )  !#KTCtodo: Add a namelist variable to set debug level
     
      !-------------------------------
      ! Define physics buffer indexes
      !-------------------------------
      thlm_idx = pbuf_get_index('THLM')   
      num_subcols_idx = pbuf_get_index('num_subcols')
      pdf_params_idx = pbuf_get_index('PDF_PARAMS')
      rcm_idx = pbuf_get_index('RCM')
      rtm_idx = pbuf_get_index('RTM')
      cld_idx = pbuf_get_index('CLD')
      alst_idx = pbuf_get_index('ALST')  ! SILHS expects clubb's cloud_frac liq stratus fraction
      ztodt_idx = pbuf_get_index('ZTODT')
      ice_supersat_idx = pbuf_get_index('ISS_FRAC')
      prec_pcw_idx = pbuf_get_index('PREC_PCW')
      snow_pcw_idx = pbuf_get_index('SNOW_PCW')
      qcsedten_idx = pbuf_get_index('QCSEDTEN')
      qrsedten_idx = pbuf_get_index('QRSEDTEN')
      qisedten_idx = pbuf_get_index('QISEDTEN')
      qssedten_idx = pbuf_get_index('QSSEDTEN')
      vtrmc_idx = pbuf_get_index('VTRMC')
      umr_idx = pbuf_get_index('UMR')
      vtrmi_idx = pbuf_get_index('VTRMI')
      ums_idx = pbuf_get_index('UMS')
     
      !-------------------------------
      ! Set up SILHS hydrometeors #KTCtodo: move microphys specification to config time,
      !        Steve wants to set up a microphysics query so I can ask the microphysics
      !        scheme which hydrometeors to use. For the future.
      !-------------------------------
      iirr = 1
      iirs = 3
      iiri  = 5
      iirg = -1

      iiNr    = 2
      iiNs = 4
      iiNi    = 6
      iiNg = -1

      hydromet_dim = 6

 
      ! Set up pdf indices, hydromet indicies, hydromet arrays, and hydromet variance ratios
      call init_pdf_hydromet_arrays_api( 1.0_core_rknd, 1.0_core_rknd,  & ! intent(in)
                                         hydromet_dim,                  & ! intent(in)
                                         iirr, iiri, iirs, iirg,        & ! intent(in)
                                         iiNr, iiNi, iiNs, iiNg,        & ! intent(in)
                                         hmp2_ip_on_hmm2_ip_slope,      & ! optional(in)
                                         hmp2_ip_on_hmm2_ip_intrcpt )     ! optional(in)

      Ncnp2_on_Ncnm2 = subcol_SILHS_ncnp2_on_ncnm2

      !-------------------------------
      ! Set up hydrometeors and correlation arrays for SILHS
      !-------------------------------
      corr_file_path_cloud = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//cloud_file_ext
      corr_file_path_below = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//below_file_ext

      iunit = getunit()


      call setup_corr_varnce_array_api( corr_file_path_cloud, corr_file_path_below, &
                                    iunit )
      call freeunit(iunit) 

#endif
      !-------------------------------
      ! Register output fields from SILHS
      ! #KTCtodo: Remove these from the default output list
      !-------------------------------
      call addfld('SILHS_NCLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Cloud Number Concentration', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Number Concentration of Rain from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_OMEGA_SCOL', (/'psubcols', 'ilev    '/), 'I', 'Pa/s', &
           'Subcolumn vertical pressure velocity', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RCM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Number Conc from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Precipitating Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RT_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg ', &
           'Subcolumn Total Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_THLM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'K', &
           'Subcolumn liquid water pot temperature', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WEIGHT_SCOL', (/'psubcols'/), 'I', 'frac', &
           'Weights for each subcolumn', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm/s', &
           'Subcolumn vertical velocity from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)

      call addfld('NR_IN_LH', (/ 'lev' /), 'I', 'm^-3', &
                  'Num Rain Conc as input to SILHS')
     call addfld('RTM_CLUBB', (/ 'ilev' /), 'I', 'kg/kg', &
                  'Input total water mixing ratio')
     call addfld('THLM_CLUBB', (/ 'ilev' /), 'I', 'K', &
                  'Input liquid water potential temperature')
     call addfld('SILHS_QC_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud water mixing ratio')
     call addfld('SILHS_QI_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud ice mixing ratio')
     call addfld('SILHS_NC_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud water number concentration')
     call addfld('SILHS_NI_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud ice number concentration')
     call addfld('AKM_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Exact Kessler autoconversion')
     call addfld('AKM_LH_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Monte Carlo estimate of Kessler autoconversion')
     call addfld('INVS_EXNER', (/ 'lev' /), 'I', 'none', &
                  'inverse EXNER function from state in subcol_SILHS')
     call addfld('SILHS_ZTODT', horiz_only, 'I', 's', & 
                  'Length of Physics timestep (for debugging)')
     call addfld('SILHS_MSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean Cloud Ice across subcolumns')
     call addfld('SILHS_STDSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of Ice across subcolumns')
     call addfld('SILHS_MSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean Cloud Liquid across subcolumns')
     call addfld('SILHS_STDSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of Liquid across subcolumns')
     call addfld('SILHS_MSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                  'Mean water vapor across subcolumns')
     call addfld('SILHS_STDSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                  'Standard deviation of water vapor across subcolumns')
     call addfld('SILHS_EFF_CLDFRAC', (/ 'lev' /), 'A', 'frac', &
                  'Calculated cloud fraction from subcolumn liq or ice') 

     call addfld('SILHS_CLUBB_PRECIP_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Precipitation fraction from CLUBB (set_up_pdf_params_incl_hydromet)')
     call addfld('SILHS_CLUBB_ICE_SS_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Ice supersaturation fraction from CLUBB')

      !call add_default('SILHS_NCLD_SCOL', 1, ' ')
      !call add_default('SILHS_NRAIN_SCOL', 1, ' ')
      !call add_default('SILHS_OMEGA_SCOL', 1, ' ')
      !call add_default('SILHS_RCM_SCOL', 1, ' ')
      !call add_default('SILHS_RICLD_SCOL', 1, ' ')
      !call add_default('SILHS_NICLD_SCOL', 1, ' ')
      !call add_default('SILHS_RRAIN_SCOL', 1, ' ')
      !call add_default('SILHS_RT_SCOL', 1, ' ')
      !call add_default('SILHS_THLM_SCOL', 1, ' ')
      !call add_default('SILHS_WEIGHT_SCOL', 1, ' ')
      !call add_default('SILHS_WM_SCOL', 1, ' ')

      !call add_default('NR_IN_LH', 1, ' ')
      !call add_default('RTM_CLUBB', 1, ' ')
      !call add_default('THLM_CLUBB', 1, ' ')
      !call add_default('SILHS_QC_IN', 1, ' ')
      !call add_default('SILHS_QI_IN', 1, ' ')
      !call add_default('AKM_CLUBB', 1, ' ')
      !call add_default('AKM_LH_CLUBB', 1, ' ')
      !call add_default('INVS_EXNER', 1, ' ')

#endif
   end subroutine subcol_init_SILHS
   
   subroutine subcol_gen_SILHS(state, tend, state_sc, tend_sc, pbuf)
      !-------------------------------
      ! This is where the subcolumns are created, and the call to
      !      generate_silhs_sample_mod_api
      !    goes out. Variables needed to make this call are pulled from the 
      !    pbuf, from module data, and calculated based on the CAM state.
      !-------------------------------

      use physics_buffer,         only : physics_buffer_desc, pbuf_get_index, &
                                         pbuf_get_field
      use ppgrid,                 only : pver, pverp, pcols
      use ref_pres,               only : top_lev => trop_cloud_top_lev
      use time_manager,           only : get_nstep
      use subcol_utils,           only : subcol_set_subcols, subcol_set_weight
      use phys_control,           only : phys_getopts
      use spmd_utils,             only : masterproc
      use shr_const_mod,          only : SHR_CONST_PI, SHR_CONST_RHOFW

#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,       only : pdf_parameter, unpack_pdf_params_api, &
                                         num_pdf_params, &
                                         hydromet_dim, &

                                         Lscale, &

                                         setup_pdf_parameters_api, &

                                         l_stats_samp, &

                                         hydromet_pdf_parameter, &

                                         zm2zt_api, setup_grid_heights_api, gr, &

                                         iirr, iiNr, iirs, iiri, &
                                         iirg, iiNs, &
                                         iiNi, iiNg, &

                                         core_rknd, &

                                         w_tol_sqd, zero_threshold, cloud_frac_min, & ! rc_tol, &

                                         pdf_dim, &
                                         corr_array_n_cloud, &
                                         corr_array_n_below, &
                                         iiPDF_chi, iiPDF_rr, &
                                         iiPDF_w, iiPDF_Nr, &
                                         iiPDF_ri, iiPDF_Ni, &
                                         iiPDF_Ncn, iiPDF_rs, iiPDF_Ns
   
      use silhs_api_module, only :       generate_silhs_sample_api, & ! Ncn_to_Nc, &
                                         lh_clipped_variables_type, &
                                         clip_transform_silhs_output_api, &
                                         est_kessler_microphys_api
#endif
#endif
      
      ! CAM data structures
      type(physics_state), intent(inout) :: state
      type(physics_tend),  intent(inout) :: tend
      type(physics_state), intent(inout) :: state_sc        ! sub-column state
      type(physics_tend),  intent(inout) :: tend_sc         ! sub-column tend
      type(physics_buffer_desc), pointer :: pbuf(:)

#ifdef CLUBB_SGS
#ifdef SILHS
      !----------------
      ! Local variables
      !----------------
      logical, parameter :: &
                 l_implemented = .true.   ! Implemented in a host model
      logical, parameter :: rx_Nc = .false. ! Use NC calculated based on grid mean effective radius
      integer, parameter :: &
                 grid_type = 3            ! The 3 grid centered on momentum levels
      real(r8), parameter :: cldmin = 0.001_r8 ! To use when cld frac = 0.0 to be consistant with micro_mg
      real(r8), parameter :: min_num_conc = 1.0e-12_r8
      real(r8), parameter :: qsmall = 1.0e-18  ! Microphysics cut-off for cloud

      integer :: i, j, k, ngrdcol, ncol, lchnk, stncol
      integer :: ixcldice, ixnumice, ixq, ixcldliq, ixnumliq, ixrain, ixnumrain, ixsnow, ixnumsnow
      integer :: begin_height, end_height ! Output from setup_grid call
      real(r8) :: sfc_elevation  ! Surface elevation
      real(r8), dimension(pverp-top_lev+1) :: zt_g, zi_g ! Thermo & Momentum grids for clubb
      real(r8), dimension(pverp) :: scfrac     ! cloud fraction based on sc distributions
      real(r8) :: msc, std, maxcldfrac, maxsccldfrac
      real(r8) :: scale = 1.0_r8

      !----------------
      ! Required for set_up_pdf_params_incl_hydromet
      !----------------
      real(r8), dimension(pverp-top_lev+1) :: cld_frac_in  ! Cloud fraction
      type(hydromet_pdf_parameter), dimension(pverp-top_lev+1) :: &
                                    hydromet_pdf_params  ! Hydrometeor PDF parameters
      real(r8), dimension(:,:,:), allocatable :: &       ! Correlation matrix for pdf components
                                    corr_array_1, corr_array_2 
      real(r8), dimension(:,:), allocatable :: &
                                    mu_x_1, mu_x_2, &    ! Mean array for PDF components
                                    sigma_x_1, sigma_x_2 ! Std dev arr for PDF components
      real(r8), dimension(:,:,:), allocatable :: &       ! Transposed corr cholesky mtx
                                    corr_cholesky_mtx_1, corr_cholesky_mtx_2
      real(r8), dimension(pverp-top_lev+1) :: Nc_in_cloud
      real(r8), dimension(pverp-top_lev+1) :: ice_supersat_frac_in
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: hydrometp2


      !----------------
      ! Input to generate_silhs_sample
      !----------------
      integer :: iter                            ! CLUBB iteration 
      integer :: num_subcols                     ! Number of subcolumns
      integer, dimension(pcols) :: numsubcol_arr ! To set up the state struct
      integer, parameter :: sequence_length = 1  ! Number of timesteps btn subcol calls
      type(pdf_parameter), dimension(pverp-top_lev+1) :: pdf_params     ! PDF parameters
      real(r8), dimension(pverp, num_pdf_params) :: pdf_params_packed
      real(r8), dimension(pverp-top_lev+1) :: rho_ds_zt    ! Dry static density (kg/m^3) on thermo levs
      real(r8), dimension(pver)  :: dz_g         ! thickness of layer
      real(r8), dimension(pverp-top_lev+1) :: delta_zm     ! Difference in u wind altitudes
      real(r8), dimension(pverp-top_lev+1) :: invs_dzm     ! 1/delta_zm
      real(r8), dimension(pverp-top_lev+1) :: rcm_in       ! Cld water mixing ratio on CLUBB levs
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: hydromet  ! Hydrometeor species
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: wphydrometp  ! Hydrometeor flux
      real(r8), dimension(pverp-top_lev+1)              :: Ncm ! Mean cloud droplet concentration, <N_c>
      logical, parameter :: &  
         l_calc_weights_all_levs = .false. ! .false. if all time steps use the same
                                          !   weights at all vertical grid levels 
      logical :: & 
        l_calc_weights_all_levs_itime, & ! .true. if we calculate sample weights separately at all 
                                         !    grid levels at the current time step   
        l_rad_itime                      ! .true. if we calculate radiation at the current time step  
      
      !---------------
      !Output from generate_silhs_sample
      !--------------
      real(r8), allocatable, dimension(:,:,:) :: X_nl_all_levs ! Sample transformed to normal-lognormal
      real(r8), allocatable, dimension(:,:,:) :: X_nl_all_levs_raw ! Raw (unclipped) SILHS samples
      real(r8), allocatable, dimension(:,:)   :: lh_sample_point_weights ! Subcolumn weights
      integer,  allocatable, dimension(:,:)    :: X_mixt_comp_all_levs ! Which Mixture Component

      real(r8), allocatable, dimension(:,:) :: rc_all_points ! Calculate RCM from LH output
      real(r8), allocatable, dimension(:,:) :: rain_all_pts  ! Calculate Rain from LH output
      real(r8), allocatable, dimension(:,:) :: nrain_all_pts ! Calculate Rain Conc from LH
      real(r8), allocatable, dimension(:,:) :: snow_all_pts  ! Calculate Snow from LH output
      real(r8), allocatable, dimension(:,:) :: nsnow_all_pts ! Calculate Snow Conc from LH
      real(r8), allocatable, dimension(:,:) :: w_all_points  ! Calculate W from LH output
      ! real(r8), allocatable, dimension(:,:) :: RVM_lh_out    ! Vapor mixing ratio sent away
      real(r8), allocatable, dimension(:,:) :: ice_all_pts   ! Calculate Cld Ice from LH output
      real(r8), allocatable, dimension(:,:) :: nice_all_pts  ! Calculate Num cld ice from LH
      real(r8), allocatable, dimension(:,:) :: nclw_all_pts  ! Calculate Num cld wat from LH

      !----------------
      ! Output from clip_transform_silhs_output_api
      !----------------
      type(lh_clipped_variables_type), dimension(:,:), allocatable :: &
        lh_clipped_vars

      logical, parameter :: &
        l_use_Ncn_to_Nc = .true.  ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                                  ! Ncn_to_Nc might cause problems with the MG microphysics 
                                  ! since the changes made here (Nc-tendency) are not fed into 
                                  ! the microphysics
        

      !----------------
      ! Output to history
      !----------------
      ! V. Larson note: These variables are on the zt (full) levels: why do they
      ! have dimension pverp?  The pverp level corresponds to the CLUBB
      ! below-ground level.
      ! The variables in this paragraph are oriented like CAM variables (k=1 is
      ! the model top).
      ! They are flipped versions of CLUBB variables, for the entire chunk.
      real(r8), dimension(pcols*psubcols, pverp) :: RT_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: THL_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: OMEGA_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: WM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RVM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RCM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NCLW_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: ICE_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NICE_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RAIN_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NRAIN_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: SNOW_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NSNOW_lh_out

      real(r8), dimension(state_sc%psetcols) :: weights ! Subcol weights

      real(r8), dimension(pcols, pver) :: meansc_ice
      real(r8), dimension(pcols, pver) :: stdsc_ice

      real(r8), dimension(pcols, pver) :: meansc_liq
      real(r8), dimension(pcols, pver) :: stdsc_liq

      real(r8), dimension(pcols, pver) :: meansc_vap
      real(r8), dimension(pcols, pver) :: stdsc_vap
      real(r8), dimension(pcols, pver) :: grmn_eff_rad
      real(r8), dimension(pcols, pver) :: eff_cldfrac
      real(r8), dimension(pcols, pver) :: precip_frac_out

      real(r8) :: tmp_mean, diff_mean, rcubed

      !----------------
      ! Output from Est_Kessler_microphys
      !----------------
      real(r8), dimension(pverp-top_lev+1) :: lh_Akm     ! Monte Carlo estimate of Kessler Autoconversion
      real(r8), dimension(pverp-top_lev+1) :: AKm        ! Exact Kessler autoconversion
      real(r8), dimension(pverp-top_lev+1) :: AKstd      ! Exact Stdev of gba Kessler
      real(r8), dimension(pverp-top_lev+1) :: AKstd_cld  ! Exact w/in cloud stdev of gba Kessler
      real(r8), dimension(pverp-top_lev+1) :: AKm_rcm    ! Exact local gba Kessler auto based on rcm
      real(r8), dimension(pverp-top_lev+1) :: AKm_rcc    ! Exact local gba Kessler based on w/in cloud rc
      real(r8), dimension(pverp-top_lev+1) :: lh_rcm_avg ! LH estimate of grid box avg liquid water
      real(r8), dimension(pcols,pverp) :: lh_AKm_out, AKm_out

      !----------------
      ! Needed to update State
      !----------------
      real(r8), dimension(pver)  :: Temp_prof  ! Subcolumn LWPT converted to Abs Temp
      real(r8), dimension(pver)  :: SE_prof    ! Static Energy calculated from Abs Temp
      real(r8), dimension(pver)  :: No_cloud = 0.0_r8     ! Clear air condensate profile
      real(r8), dimension(pcols, pver)  :: invs_exner  ! inverse exner sent to conversion codw
                                                       ! pcols for output to history
      real(r8) :: eff_rad_coef = 1.0_r8/(4.0_r8/3.0_r8*SHR_CONST_RHOFW*SHR_CONST_PI)
      real(r8), dimension(pver) :: eff_rad_prof ! r^3 as calculated from grid mean MR & NC
     
      !----------------
      ! Pointers
      !----------------
      real(r8), pointer, dimension(:) :: num_subcol_ptr
      real(r8), pointer, dimension(:) :: ztodt_ptr
      real(r8), pointer, dimension(:,:) :: thlm      ! Mean temperature
      real(r8), pointer, dimension(:,:,:) :: pdf_params_ptr  ! Packed PDF_Params
      real(r8), pointer, dimension(:,:) :: ice_supersat_frac ! ice cloud fraction
      real(r8), pointer, dimension(:,:) :: rcm       ! CLUBB cld water mr
      real(r8), pointer, dimension(:,:) :: rtm       ! mean moisture mixing ratio
      real(r8), pointer, dimension(:,:) :: cld       ! CAM cloud fraction
      real(r8), pointer, dimension(:,:) :: alst      ! CLUBB liq cloud fraction
      real(r8), pointer, dimension(:,:) :: qrain     ! micro_mg rain from previous step
      real(r8), pointer, dimension(:,:) :: qsnow     
      real(r8), pointer, dimension(:,:) :: nrain     ! micro_mg rain num conc 
      real(r8), pointer, dimension(:,:) :: nsnow


      if (.not. allocated(state_sc%lat)) then
         call endrun('subcol_gen error: state_sc must be allocated before calling subcol_gen')
      end if

      ! Determine num of columns and which chunk we're working on and what timestep
      ngrdcol = state%ngrdcol
      ncol = state%ncol
      lchnk = state%lchnk
      iter = get_nstep() ! #KTCtodo: The model iteration is passed into SILHS without taking
                         !           substepping into account. I may need to change this in 
                         !           the future. Also, why does SILHS need an iter, but CLUBB
                         !           does not?
                         ! #ERDBG:   The model iteration number is not used in SILHS unless
                         !           sequence_length > 1, but nobody runs with that option.
      !----------------
      ! Establish associations between pointers and physics buffer fields
      ! (Do this now so that num_subcol_ptr is available for the state copy below)
      !----------------
      call pbuf_get_field(pbuf, thlm_idx, thlm)
      call pbuf_get_field(pbuf, num_subcols_idx, num_subcol_ptr)
      call pbuf_get_field(pbuf, ztodt_idx, ztodt_ptr)
      call pbuf_get_field(pbuf, pdf_params_idx, pdf_params_ptr)
      call pbuf_get_field(pbuf, ice_supersat_idx, ice_supersat_frac)
      call pbuf_get_field(pbuf, rcm_idx, rcm)
      call pbuf_get_field(pbuf, rtm_idx, rtm)
      call pbuf_get_field(pbuf, alst_idx, alst)
      call pbuf_get_field(pbuf, cld_idx, cld)
      call pbuf_get_field(pbuf, qrain_idx, qrain)
      call pbuf_get_field(pbuf, qsnow_idx, qsnow)
      call pbuf_get_field(pbuf, nrain_idx, nrain)
      call pbuf_get_field(pbuf, nsnow_idx, nsnow)

      !----------------
      ! Copy state and populate numbers and values of sub-columns
      !----------------
      ztodt = ztodt_ptr(1)
      numsubcol_arr(:) = 0  ! Start over each chunk
      numsubcol_arr(:ngrdcol) = subcol_SILHS_numsubcol ! Only set for valid grid columns
      call subcol_set_subcols(state, tend, numsubcol_arr, state_sc, tend_sc)

      !----------------
      ! Get indices for ice mass and number
      ! This is the same code from clubb_intr.F90
      !----------------
      call cnst_get_ind('Q', ixq)
      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('NUMICE', ixnumice)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('NUMLIQ', ixnumliq)
      call cnst_get_ind('RAINQM', ixrain, abort=.false.)
      call cnst_get_ind('NUMRAI', ixnumrain, abort=.false.)
      call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)
      call cnst_get_ind('NUMSNO', ixnumsnow, abort=.false.)

      ! The number of vertical grid levels used in CLUBB is pverp, which is originally
      ! set in the call to setup_clubb_core_api from subroutine clubb_ini_cam.  This
      ! is stored in CLUBB in the object gr%nz.  This isn't changed in CLUBB.
      ! However, when SILHS is used, SILHS only uses pverp - top_lev + 1 vertical grid
      ! levels and also uses the gr%nz object.  The value of gr%nz needs to be reset
      ! for SILHS here and then set again for CLUBB in subroutine clubb_tend_cam.
      gr%nz = pverp - top_lev + 1

      !----------------
      ! Loop over all the active grid columns in the chunk
      !----------------
      do i = 1, ngrdcol
      
         ! JHDBG: Big suspicion about that code
         ! V. Larson: I don't know what happens to arrays allocated with size
         ! num_subcols if num_subcols varies with the grid column.
         num_subcols = numsubcol_arr(i)
         stncol = 0         ! Each grid column needs to know how many subcolumns have gone by
         do k = 1, i-1
            ! stncol = stncol + numsubcol_arr(i-1)
            ! Eric Raut replaced i-1 with k in line immediately above.
            stncol = stncol + numsubcol_arr(k)
         enddo

         ! Setup the CLUBB vertical grid object. This must be done for each
         ! column as the z-distance between hybrid pressure levels can 
         ! change easily.
         sfc_elevation = state%zi(i,pverp)
         ! Define the CLUBB momentum grid (in height, units of m)
         do k = 1, pverp-top_lev+1
            zi_g(k) = state%zi(i,pverp-k+1)-sfc_elevation
         enddo
         ! Define the CLUBB thermodynamic grid (in units of m)
         do k = 1, pver-top_lev+1
            zt_g(k+1) = state%zm(i,pver-k+1)-state%zi(i,pverp)
         enddo
         ! Thermodynamic ghost point is below surface
         zt_g(1) = -1._r8*zt_g(2)
         ! Calculate the distance between grid levels on the host model grid,
         ! using host model grid indices.
         do k = top_lev, pver
            dz_g(k) = state%zi(i,k)-state%zi(i,k+1)
         enddo
         ! allocate grid object
         call setup_grid_heights_api( l_implemented, grid_type, &
                                      zi_g(2), zi_g(1), zi_g(1:pverp-top_lev+1), &
                                      zt_g(1:pverp-top_lev+1) )

         ! Pull pdf params out of 2-D real array
         ! The PDF parameters are passed out of CLUBB and into SILHS without
         ! being used at all in the rest of the host model code.  The arrays
         ! aren't flipped for the PDF parameters, and they don't need to be.
         call unpack_pdf_params_api(pdf_params_ptr(i,1:pverp-top_lev+1,:), pverp-top_lev+1, pdf_params)

         ! Inverse delta_zm is required for the 3-level L-scale averaging
         do k = 1, pver-top_lev+1
            delta_zm(k+1) = state%zi(i,pverp-k)-state%zi(i,pverp-k+1)
            invs_dzm(k+1) = 1.0_r8/delta_zm(k+1)
         enddo
         ! Handle CLUBB sub-sfc ghost point as done in clubb grid_class.F90
         delta_zm(1) = delta_zm(2) 
         invs_dzm(1) = invs_dzm(2)

         ! Compute dry static density on CLUBB vertical grid
         do k = 1, pver-top_lev+1
            rho_ds_zt(k+1) = (1._r8/gravit)*state%pdel(i,pver-k+1)/dz_g(pver-k+1)
         enddo
         ! CLUBB ghost point under the surface
         rho_ds_zt(1) = rho_ds_zt(2)

         ! Set up hydromet array, flipped from CAM vert grid to CLUBB
         do k = 1, pver-top_lev+1
            if ( iirr > 0 ) then
              ! If ixrain and family are greater than zero, then MG2 is
              ! being used, and rain and snow are part of state. Otherwise,
              ! diagnostic rain and snow from MG1 are used in hydromet.
               if (ixrain > 0) then
                  hydromet(k+1,iirr) = state%q(i,pver-k+1,ixrain)
               else
                  hydromet(k+1,iirr) = qrain(i,pver-k+1)
               endif
            endif
            if ( iiNr > 0 ) then
               if (ixnumrain > 0) then
                  hydromet(k+1,iiNr) = state%q(i,pver-k+1,ixnumrain)
               else
                  hydromet(k+1,iiNr) = nrain(i,pver-k+1)
               endif
            endif
            if ( iirs > 0 ) then
               if (ixsnow > 0) then
                  hydromet(k+1,iirs) = state%q(i,pver-k+1,ixsnow)
               else
                  hydromet(k+1,iirs) = qsnow(i,pver-k+1)
               endif
            endif
            if ( iiNs > 0 ) then
               if (ixnumsnow > 0) then
                  hydromet(k+1,iiNs) = state%q(i,pver-k+1,ixnumsnow)
               else
                  hydromet(k+1,iiNs) = nsnow(i,pver-k+1)
               endif
            endif
            if ( iiri > 0 ) then
               hydromet(k+1,iiri) = state%q(i,pver-k+1,ixcldice)
            endif
            if ( iiNi > 0 ) then
               hydromet(k+1,iiNi) = state%q(i,pver-k+1,ixnumice)
            endif
     
            Ncm(k+1) = state%q(i,pver-k+1,ixnumliq)

         enddo

         do k = 1, hydromet_dim ! ghost point below the surface
            hydromet(1,k) = hydromet(2,k)                  
         enddo

         Ncm(1) = Ncm(2)

         do k = top_lev, pver
            ! Calculate effective radius cubed, CAM-grid oriented for use in subcolumns
            eff_rad_prof(k) = eff_rad_coef*state%q(i,k,ixcldliq)/state%q(i,k,ixnumliq)
            ! Test a fixed effective radius
            ! eff_rad_prof(k) = 5.12e-16_r8 ! 8 microns
         enddo

         ! Allocate arrays for set_up_pdf_params_incl_hydromet
         allocate( corr_array_1(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( corr_array_2(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( mu_x_1(pdf_dim, pverp-top_lev+1) )
         allocate( mu_x_2(pdf_dim, pverp-top_lev+1) )
         allocate( sigma_x_1(pdf_dim, pverp-top_lev+1) )
         allocate( sigma_x_2(pdf_dim, pverp-top_lev+1) )
         allocate( corr_cholesky_mtx_1(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( corr_cholesky_mtx_2(pdf_dim, pdf_dim, pverp-top_lev+1) )
         ! Allocate arrays for SILHS output
         allocate( lh_sample_point_weights(pverp-top_lev+1,num_subcols) )
         allocate( X_mixt_comp_all_levs(pverp-top_lev+1,num_subcols) )
         allocate( X_nl_all_levs(pverp-top_lev+1,num_subcols,pdf_dim) )
         allocate( X_nl_all_levs_raw(pverp-top_lev+1,num_subcols,pdf_dim) )
         allocate( lh_clipped_vars(pverp-top_lev+1,num_subcols) )
         ! Allocate arrays for output to either history files or for updating state_sc
         allocate( rc_all_points(pverp-top_lev+1, num_subcols) )
         allocate( rain_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nrain_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( snow_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nsnow_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( w_all_points(pverp-top_lev+1, num_subcols) )
         ! allocate( RVM_lh_out(num_subcols, pverp) )  ! This one used only to update state
         allocate( ice_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nice_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nclw_all_pts(pverp-top_lev+1, num_subcols) )
         
         ! Convert from CAM vertical grid to CLUBB
         do k = 1, pverp-top_lev+1 
            rcm_in(k)  = rcm(i,pverp-k+1)
            ice_supersat_frac_in(k) = ice_supersat_frac(i,pverp-k+1)
         enddo
         do k = 1, pver-top_lev+1
            cld_frac_in(k+1) = alst(i,pver-k+1)
         enddo
         cld_frac_in(1) = cld_frac_in(2) ! Ghost pt below surface
         ! Calculate a clubb-specific exner function
         ! (This is grid mean, as pressure levels do not change in 
         !  the subcolumn state)
         invs_exner(i,:) = ((state%pmid(i,:)/p0_clubb)**(rair/cpair))

         ! Call setup_pdf_parameters to get the CLUBB PDF ready for SILHS
         ! Compute Num concentration of cloud nuclei
         Nc_in_cloud = Ncm / max( cld_frac_in, cloud_frac_min )

         ! The variable wphydrometp is only used when l_calc_w_corr is enabled.
         ! The l_calc_w_corr flag is turned off by default, so wphydrometp will
         ! simply be set to 0 to simplify matters.
         wphydrometp = 0.0_r8
     
         ! make the call
         call setup_pdf_parameters_api( pverp-top_lev+1, pdf_dim, ztodt, &    ! In
                                        Nc_in_cloud, rcm_in, cld_frac_in, &            ! In
                                        ice_supersat_frac_in, hydromet, wphydrometp, & ! In
                                        corr_array_n_cloud, corr_array_n_below, &      ! In
                                        pdf_params, l_stats_samp, &                    ! In
                                        hydrometp2, &                                  ! Out
                                        mu_x_1, mu_x_2, &                              ! Out
                                        sigma_x_1, sigma_x_2, &                        ! Out
                                        corr_array_1, corr_array_2, &                  ! Out
                                        corr_cholesky_mtx_1, corr_cholesky_mtx_2, &    ! Out
                                        hydromet_pdf_params )                          ! Out

         ! Calculate radiation only once in a while
         ! l_rad_itime = (mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1)  

         ! Calculate sample weights separately at all grid levels when
         ! radiation is not called  
         ! l_calc_weights_all_levs_itime = l_calc_weights_all_levs .and. .not.
         ! l_rad_itime  
         l_calc_weights_all_levs_itime = .false. ! subcol_utils cannot compute weighted avgs
                                                 !   when the weights vary with height.   
                                                 !   Don't set to true until this is fixed!!

         ! Lscale needs to be passed out of advance_clubb_core, saved to the
         ! pbuf, and then pulled out of the pbuf for use here.

         ! Let's generate some subcolumns!!!!!
         call generate_silhs_sample_api &
              ( iter, pdf_dim, num_subcols, sequence_length, pverp-top_lev+1, & ! In
                l_calc_weights_all_levs_itime, &                   ! In 
                pdf_params, delta_zm, rcm_in, Lscale(1:pverp-top_lev+1), &      ! In
                rho_ds_zt, mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, & ! In 
                corr_cholesky_mtx_1, corr_cholesky_mtx_2, &        ! In
                hydromet_pdf_params, &                             ! In
                X_nl_all_levs_raw, X_mixt_comp_all_levs, &         ! Out
                lh_sample_point_weights)                           ! Out

         ! Extract clipped variables from subcolumns
         call clip_transform_silhs_output_api( pverp-top_lev+1, num_subcols, &   ! In
                                               pdf_dim, hydromet_dim, & ! In
                                               X_mixt_comp_all_levs, & ! In
                                               X_nl_all_levs_raw, &    ! In
                                               pdf_params, l_use_Ncn_to_Nc, & ! In
                                               lh_clipped_vars, & ! Out
                                               X_nl_all_levs )    ! Out

         ! Test subcolumns by comparing to an estimate of kessler autoconversion
         call est_kessler_microphys_api &
              ( pverp-top_lev+1, num_subcols, pdf_dim, X_nl_all_levs, pdf_params, &
                rcm_in, cld_frac_in, X_mixt_comp_all_levs, lh_sample_point_weights, &
                lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, lh_rcm_avg)

         ! Calc column liquid water for output (rcm)
         rc_all_points = lh_clipped_vars(:,:)%rc

         if ( iiPDF_rr > 0 ) then
             ! Calc subcolumn precipitating liq water for output (rrm)
             rain_all_pts = real( X_nl_all_levs(:,:,iiPDF_rr), kind=r8 )
         end if

         if ( iiPDF_Nr > 0 ) then
             ! Calc subcolumn number rain conc for output (nrainm)
             nrain_all_pts = real( X_nl_all_levs(:,:,iiPDF_Nr), kind=r8 )
         end if

         if ( iiPDF_rs > 0 ) then
             ! Calc subcolumn precipitating snow      for output (rsm)
             snow_all_pts = real( X_nl_all_levs(:,:,iiPDF_rs), kind=r8 )
         end if

         if ( iiPDF_Ns > 0 ) then
             ! Calc subcolumn precipitating snow conc for output (Nsm)
             nsnow_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ns), kind=r8 )
         end if

         if ( iiPDF_ri > 0 ) then
             ! Calc subcolumn cloud ice mixing ratio
             ice_all_pts = real( X_nl_all_levs(:,:,iiPDF_ri), kind=r8)
         end if

         if ( iiPDF_Ni > 0 ) then
             ! Calc subcolumn cloud ice number
             nice_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ni), kind=r8)
         end if

         ! Calc subcolumn vert velocity for output (wm)
         w_all_points = real( X_nl_all_levs(:,:,iiPDF_w), kind=r8 )
         ! Calc cloud liq water number conc 
         nclw_all_pts = lh_clipped_vars(:,:)%Nc
         ! Calc mean liquid water potential temp for clear air
         !call THL_profile(pver, state%t(i,:), invs_exner(i,:), No_cloud, Temp_prof)

         ! Calc effective cloud fraction for testing
         eff_cldfrac(:,:) = 0.0_r8
         do k = top_lev, pver
            do j=1, num_subcols

               if ( ( rc_all_points(pverp-k+1,j) .gt. qsmall ) &
                      .or. ( ice_all_pts(pverp-k+1,j) .gt. qsmall ) ) then
                  eff_cldfrac(i,k) = eff_cldfrac(i,k)+lh_sample_point_weights(pverp-k+1,j)
               endif
            enddo 

            eff_cldfrac(i,k) = eff_cldfrac(i,k)/real(num_subcols, kind=r8)
         enddo

         ! Pack precip_frac for output
         do k = 2, pverp-top_lev+1
           precip_frac_out(i,pver-k+2) = hydromet_pdf_params(k)%precip_frac
         enddo

         ! Pack up weights for output
         do j = 1, num_subcols      
            if (subcol_SILHS_weight) then 
               weights(stncol+j) = lh_sample_point_weights(2,j) ! Using grid level 2 always won't work 
                                                                !   if weights vary with height.
            else
               weights(stncol+j) = 1._r8
            endif
         enddo

         ! Convert from CLUBB vertical grid to CAM grid for history output and
         ! Updating state variables
         do k = top_lev, pverp
            do j = 1, num_subcols
               RT_lh_out(    stncol+j,k ) = lh_clipped_vars(pverp-k+1,j)%rt
               RCM_lh_out(   stncol+j,k ) = rc_all_points(pverp-k+1,j)
               NCLW_lh_out(  stncol+j,k ) = nclw_all_pts(pverp-k+1,j)
               ICE_lh_out(   stncol+j,k ) = ice_all_pts(pverp-k+1,j)
               NICE_lh_out(  stncol+j,k ) = nice_all_pts(pverp-k+1,j)
!               RVM_lh_out(j,k) = RT_lh_out(stncol+j,k)-RCM_lh_out(stncol+j,k)-ICE_lh_out(stncol+j,k)
               RVM_lh_out(   stncol+j,k ) = lh_clipped_vars(pverp-k+1,j)%rv
               THL_lh_out(   stncol+j,k ) = lh_clipped_vars(pverp-k+1,j)%thl
               RAIN_lh_out(  stncol+j,k ) = rain_all_pts(pverp-k+1,j)
               NRAIN_lh_out( stncol+j,k ) = nrain_all_pts(pverp-k+1,j)
               SNOW_lh_out(  stncol+j,k ) = snow_all_pts(pverp-k+1,j)
               NSNOW_lh_out( stncol+j,k ) = nsnow_all_pts(pverp-k+1,j)
               WM_lh_out(    stncol+j,k ) = w_all_points(pverp-k+1,j)
               OMEGA_lh_out( stncol+j,k ) = -1._r8*WM_lh_out(stncol+j,k)*rho_ds_zt(pverp-k+1)*gravit
               AKm_out(i,k) = AKm(pverp-k+1)
               lh_AKm_out(i,k) = lh_AKm(pverp-k+1)
            enddo
         enddo

         ! Constrain the sample distribution of cloud water and ice to the same mean
         ! as the grid to prevent negative condensate errors
         if(subcol_SILHS_constrainmn) then
            call subcol_constrainmn( num_subcols, ICE_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldice), meansc_ice(i,:), stdsc_ice(i,:) )
            if ( ixrain > 0 ) &
            call subcol_constrainmn( num_subcols, RAIN_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixrain) )
            if ( ixsnow > 0 ) &
            call subcol_constrainmn( num_subcols, SNOW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixsnow) )
            call subcol_constrainmn( num_subcols, RCM_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldliq), meansc_liq(i,:), stdsc_liq(i,:) )
            call subcol_constrainmn( num_subcols, RVM_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixq), meansc_vap(i,:), stdsc_vap(i,:) )
            call subcol_constrainmn( num_subcols, NICE_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumice) )
            if ( ixnumrain > 0 ) &
            call subcol_constrainmn( num_subcols, NRAIN_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumrain) )
            if ( ixnumsnow > 0 ) &
            call subcol_constrainmn( num_subcols, NSNOW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumsnow) )
            call subcol_constrainmn( num_subcols, NCLW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumliq) )
            do k = top_lev, pver
               ! Look for exceptionally large values of condensate
               if(ANY(ICE_lh_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)) then
                  ! Clip the large values
                  where(ICE_lh_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)
                     ICE_lh_out(stncol+1:stncol+num_subcols,k) = 0.01_r8
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = 1.5e+7_r8
                  end where
                  ! Recalculate the weighted subcolumn mean
                  tmp_mean = meansc( ICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  ! Add the difference to each subcolumn
                  ICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                     ICE_lh_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  ! Recalculate the weight subcolumn mean for ice num conc
                  tmp_mean = meansc( NICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
                  ! Add the difference to each subcolumn
                  if(diff_mean.gt.0.0_r8) then
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                         NICE_lh_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  else ! just use the grid mean in each subcolumn
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                         state%q(i,k,ixnumice)
                  end if
                  ! Test adjusted means for debugging
                  tmp_mean = meansc( ICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  tmp_mean = meansc( NICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
               endif
            enddo ! k = top_lev, pver
         endif ! subcol_silhs_constrainm

         ! Code to update the state variables for interactive runs
         ! Set state variables
         do j = 1, numsubcol_arr(i)

            call Abs_Temp_profile( pver-top_lev+1, THL_lh_out(stncol+j,top_lev:pver), &
                                   invs_exner(i,top_lev:pver), RCM_lh_out(stncol+j,top_lev:pver), &
                                   Temp_prof(top_lev:pver) )
            state_sc%t(stncol+j,top_lev:pver) = Temp_prof(top_lev:pver)
            call StaticEng_profile( pver-top_lev+1, Temp_prof(top_lev:pver), &
                                    state%zm(i,top_lev:pver), state%phis(i), &
                                    SE_prof(top_lev:pver) )
            state_sc%s(stncol+j,top_lev:pver) = SE_prof(top_lev:pver)

            ! Vertical Velocity is not part of the energy conservation checks, but
            ! we need to be careful here, because the SILHS output VV is noisy.
            state_sc%omega(stncol+j,top_lev:pver) = OMEGA_lh_out(stncol+j,top_lev:pver)
            state_sc%q(stncol+j,top_lev:pver,ixq) = RVM_lh_out(stncol+j,top_lev:pver) 

            if( rx_Nc ) then
                stop "rx_Nc not enabled in subcol_gen_SILHS"
!               ! Test calculating num const based on grid mean eff radius
!                      where(eff_rad_prof.gt.0.0)
!                state_sc%q(stncol+j,:,ixnumliq) = (RCM_ADJ_out(stncol+j,1:pver) &
!                                                         /eff_rad_prof(:))*eff_rad_coef
!                      elsewhere
!                        state_sc%q(stncol+j,:,ixnumliq) = NCLW_ADJ_out(stncol+j,1:pver)
!                      end where
!                      NCLW_lh_out(stncol+j,1:pver) = state_sc%q(stncol+j,:,ixnumliq) 
!           else
!                       state_sc%q(stncol+j,:,ixnumliq) = NCLW_ADJ_out(stncol+j,1:pver)
!                       NCLW_lh_out(stncol+j,:) = NCLW_ADJ_out(stncol+j,:)
            endif


            if (subcol_SILHS_meanice) then
                stop "subcol_SILHS_meanice = T not currently available"
                state_sc%q(stncol+j,top_lev:pver,ixcldice) = state%q(i,top_lev:pver,ixcldice)
                state_sc%q(stncol+j,top_lev:pver,ixnumice) = state%q(i,top_lev:pver,ixnumice)
                state_sc%q(stncol+j,top_lev:pver,ixcldliq) = RCM_lh_out(stncol+j,top_lev:pver)
                state_sc%q(stncol+j,top_lev:pver,ixnumliq) = NCLW_lh_out(stncol+j,top_lev:pver)
            else
               if (subcol_SILHS_q_to_micro) then ! Send SILHS predicted constituents to microp
                   state_sc%q(stncol+j,top_lev:pver,ixcldliq) = RCM_lh_out(stncol+j,top_lev:pver)
                   state_sc%q(stncol+j,top_lev:pver,ixcldice) = ICE_lh_out(stncol+j,top_lev:pver)
                   if (ixrain > 0) &
                      state_sc%q(stncol+j,top_lev:pver,ixrain) = RAIN_lh_out(stncol+j,top_lev:pver)
                   if (ixsnow > 0) &
                      state_sc%q(stncol+j,top_lev:pver,ixsnow) = SNOW_lh_out(stncol+j,top_lev:pver)
               else            
                  state_sc%q(stncol+j,top_lev:pver,ixcldliq) = state%q(i,top_lev:pver,ixcldliq)
                  state_sc%q(stncol+j,top_lev:pver,ixcldice) = state%q(i,top_lev:pver,ixcldice)
                  if (ixrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixrain) = state%q(i,top_lev:pver,ixrain)
                  if (ixsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixsnow) = state%q(i,top_lev:pver,ixsnow)
               endif
               if (subcol_SILHS_n_to_micro) then ! Send SILHS predicted number conc to microp
                  state_sc%q(stncol+j,top_lev:pver,ixnumice) = NICE_lh_out(stncol+j,top_lev:pver)
                  state_sc%q(stncol+j,top_lev:pver,ixnumliq) = NCLW_lh_out(stncol+j,top_lev:pver)
                  if (ixnumrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumrain) = NRAIN_lh_out(stncol+j,top_lev:pver)
                  if (ixnumsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = NSNOW_lh_out(stncol+j,top_lev:pver)
               else            
                  state_sc%q(stncol+j,top_lev:pver,ixnumliq) = state%q(i,top_lev:pver,ixnumliq)
                  state_sc%q(stncol+j,top_lev:pver,ixnumice) = state%q(i,top_lev:pver,ixnumice)
                  if (ixnumrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumrain) = state%q(i,top_lev:pver,ixnumrain)
                  if (ixnumsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = state%q(i,top_lev:pver,ixnumsnow)
               endif
            endif ! meanice

            ! Change liq and ice (and rain and snow) num conc zeros to min values (1e-12)
            where (state_sc%q(stncol+j,top_lev:pver,ixnumliq) .lt. min_num_conc) 
               state_sc%q(stncol+j,top_lev:pver,ixnumliq) = min_num_conc
            end where
            where (state_sc%q(stncol+j,top_lev:pver,ixnumice) .lt. min_num_conc)
               state_sc%q(stncol+j,top_lev:pver,ixnumice) = min_num_conc
            end where
            if (ixnumrain > 0) then
               where(state_sc%q(stncol+j,top_lev:pver,ixnumrain) .lt. min_num_conc)
                  state_sc%q(stncol+j,top_lev:pver,ixnumrain) = min_num_conc
               end where
            endif
            if (ixnumsnow > 0) then
               where(state_sc%q(stncol+j,top_lev:pver,ixnumsnow) .lt. min_num_conc)
                  state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = min_num_conc
               end where
            endif
               
         enddo

         ! Only use weights if namelist variable turned on
         if (subcol_SILHS_weight) call subcol_set_weight(state_sc%lchnk, weights)

     
         ! Deallocate the dynamic arrays used
         deallocate( lh_sample_point_weights, X_mixt_comp_all_levs, &
                     X_nl_all_levs, X_nl_all_levs_raw, lh_clipped_vars, &
                     corr_array_1, corr_array_2, mu_x_1, mu_x_2, sigma_x_1, &
                     sigma_x_2, corr_cholesky_mtx_1, corr_cholesky_mtx_2 )
         ! deallocate( RVM_lh_out ) 
         deallocate( rc_all_points, rain_all_pts, nrain_all_pts, snow_all_pts, nsnow_all_pts, ice_all_pts, &
                     nice_all_pts, nclw_all_pts, w_all_points )
      enddo ! ngrdcol

      call outfld( 'SILHS_THLM_SCOL', THL_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RT_SCOL', RT_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_OMEGA_SCOL', OMEGA_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WM_SCOL', WM_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RCM_SCOL', RCM_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RICLD_SCOL', ICE_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NICLD_SCOL', NICE_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NCLD_SCOL', NCLW_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RRAIN_SCOL', RAIN_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NRAIN_SCOL', NRAIN_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WEIGHT_SCOL', weights, pcols*psubcols, lchnk )
      call outfld( 'NR_IN_LH', nrain, pcols, lchnk )
      call outfld( 'RTM_CLUBB', rtm, pcols, lchnk )
      call outfld( 'THLM_CLUBB', thlm, pcols, lchnk )
      call outfld( 'SILHS_QC_IN', state%q(:,:,ixcldliq), pcols, lchnk )
      call outfld( 'SILHS_QI_IN', state%q(:,:,ixcldice), pcols, lchnk )
      call outfld( 'SILHS_NC_IN', state%q(:,:,ixnumliq), pcols, lchnk )
      call outfld( 'SILHS_NI_IN', state%q(:,:,ixnumice), pcols, lchnk )
      call outfld( 'AKM_CLUBB', AKm_out, pcols, lchnk )
      call outfld( 'AKM_LH_CLUBB', lh_AKm_out, pcols, lchnk )
      call outfld( 'INVS_EXNER', invs_exner, pcols, lchnk )
      call outfld( 'SILHS_ZTODT', ztodt_ptr, pcols, lchnk )
      call outfld( 'SILHS_MSC_CLDICE', meansc_ice, pcols, lchnk )
      call outfld( 'SILHS_STDSC_CLDICE', stdsc_ice, pcols, lchnk )
      call outfld( 'SILHS_MSC_CLDLIQ', meansc_liq, pcols, lchnk )
      call outfld( 'SILHS_STDSC_CLDLIQ', stdsc_liq, pcols, lchnk )
      call outfld( 'SILHS_MSC_Q', meansc_vap, pcols, lchnk )
      call outfld( 'SILHS_STDSC_Q', stdsc_vap, pcols, lchnk )
      call outfld( 'SILHS_EFF_CLDFRAC', eff_cldfrac, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_PRECIP_FRAC', precip_frac_out, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_ICE_SS_FRAC', ice_supersat_frac, pcols, lchnk )

#endif
#endif
   end subroutine subcol_gen_SILHS

   subroutine subcol_ptend_avg_SILHS(ptend_sc, ngrdcol, lchnk, ptend)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_ptend_get_firstsubcol, subcol_ptend_avg_shr, &
                                  subcol_get_weight, subcol_get_filter, &
                                  is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      type(physics_ptend), intent(in)             :: ptend_sc        ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      type(physics_ptend), intent(inout)          :: ptend
      ! Because we can't get a state passed in here, we might have to use values from the 
      ! subcolumn generation. This would make any conservation checks invalid if this
      ! function is called after another parameterization... hmm.

       call subcol_ptend_avg_shr(ptend_sc, ngrdcol, lchnk, ptend, is_filter_set(), is_weight_set())

   end subroutine subcol_ptend_avg_SILHS

#ifdef SILHS
   real(r8) function meansc(arr_in, w_in, ns) result(val)
      real(r8), intent(in) :: ns                         ! Length of Array
      real(r8), dimension(int(ns)), intent(in) :: arr_in      ! Input array
      real(r8), dimension(int(ns)), intent(in) :: w_in        ! Weights
      real(r8) :: acc  ! accumulator
      integer :: i
      acc = 0
      val = 0
      do i=1,ns
         acc = acc + arr_in(i)*w_in(i)
      enddo
      val = acc/ns
   end function

   real(r8) function stdsc(arr_in, w_in, mn_in, ns) result(val)
      real(r8), intent(in) :: ns  ! Number of elements (subcolumns)
      real(r8), dimension(int(ns)), intent(in) :: arr_in, w_in  !Input array and weights
      real(r8), intent(in) :: mn_in   ! The mean of arr_in
      real(r8) :: accvar, var
      integer :: i
      accvar = 0
      do i=1,ns
         accvar = accvar + ((arr_in(i)-mn_in)**2)*w_in(i)
      enddo
      var = accvar/ns
      val = sqrt(var)
   end function

   subroutine Abs_Temp_profile(nz, LWPT_prof, ex_prof, rcm_prof, ABST_prof)

      use clubb_api_module,              only : thlm2T_in_K_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: LWPT_prof  ! Temp prof in LWPT
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: ABST_prof  ! Abs Temp prof
      integer :: i
 
      do i=1,nz
         ABST_prof(i) = thlm2T_in_K_api(LWPT_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine THL_profile(nz, ABST_prof, ex_prof, rcm_prof, THL_prof)

      use clubb_api_module,              only : T_in_K2thlm_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: ABST_prof  ! Abs Temp prof
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: THL_prof  ! LWPT prof
      integer :: i
 
      do i=1,nz
         THL_prof(i) = T_in_K2thlm_api(ABST_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine StaticEng_profile(nz, ABST_prof, zm_prof, zsfc, s_prof)
      integer,                 intent(in) :: nz
      real(r8), dimension(nz), intent(in) :: ABST_prof
      real(r8), dimension(nz), intent(in) :: zm_prof
      real(r8),                intent(in) :: zsfc
      real(r8), dimension(nz), intent(out) :: s_prof
      integer :: i

      do i=1,nz
         s_prof(i) = cpair*(ABST_prof(i)) + gravit*zm_prof(i)+zsfc
      enddo

   end subroutine

   subroutine subcol_constrainmn( num_subcols, samples, weights, grid_mean, mean_sc, std_sc )

      use ref_pres, only : top_lev => trop_cloud_top_lev

      ! Input/Output Variables
      integer, intent(in) :: num_subcols
      real(r8), dimension(num_subcols, pverp), intent(inout) :: samples
      real(r8), dimension(num_subcols), intent(in) :: weights
      real(r8), dimension(pverp), intent(in) :: grid_mean
      real(r8), dimension(pver), intent(out), optional :: mean_sc, std_sc

      ! Local Variables
      real(r8) :: meansc_loc, adj_rat
      integer :: k
   !------------------------------------------------------------------
      !----- Begin Code -----
      do k = top_lev, pver
         meansc_loc = meansc( samples(:,k), weights(:), real(num_subcols, r8) )

         if (present(mean_sc)) &
            mean_sc(k) = meansc_loc
         if (present(std_sc)) &
            std_sc(k) = stdsc( samples(:,k), weights(:), meansc_loc, &
                               real(num_subcols, r8) )

         if ( meansc_loc > 0.0_r8 ) then
            adj_rat = grid_mean(k)/meansc_loc
         else 
            ! If the mean is zero, then zero out all subcolumns to avoid
            ! negative samples
            adj_rat = 0.0_r8
         end if
         samples(:,k) = samples(:,k) * adj_rat
      end do
   end subroutine subcol_constrainmn

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !
   function clubb_flip_grid ( profile ) result( profile_flipped )

     ! Description:
     !   Swaps the elements in profile so they are in reverse order. CAM and
     !   CLUBB's grids are flipped with respect to each other.
     !
     !   Usage:
     !     clubb_var = clubb_flip_grid( cam_var )
     !     cam_var   = clubb_flip_grid( clubb_var )

     implicit none

     ! Input Variable
     real(r8), dimension(pverp), intent(in) :: profile

     ! Output Variable
     real(r8), dimension(pverp) :: profile_flipped

     ! Local Variable
     integer :: k

     do k = 1, pverp
       profile_flipped(k) = profile(pverp-k+1)
     end do ! k=1, pverp

     return
   end function clubb_flip_grid

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !

   subroutine subcol_SILHS_var_covar_driver &
              ( ztodt, state_sc, ptend_sc, &
                pbuf )

     ! This subroutine calculates microphysical effects on five variances and
     ! covariances: rtp2, thlp2, wprtp, wpthlp, and rtpthlp.
     !
     ! This code is experimental!!

     use physics_buffer,          only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
     use ref_pres,                only: top_lev => trop_cloud_top_lev
     use subcol_utils,            only: subcol_unpack, subcol_get_nsubcol, subcol_get_weight
     use clubb_api_module,        only: T_in_K2thlm_api, pdf_parameter, unpack_pdf_params_api
     use silhs_api_module,        only: lh_microphys_var_covar_driver_api

     implicit none

     ! Parameters
     !  This fill value is set to catch errors; it should not be read.
     real(r8), parameter                   :: fillvalue = -999._r8

     ! Input Variables
     real(r8), intent(in)                  :: ztodt        ! model time increment
     type(physics_state), intent(in)       :: state_sc     ! state for sub-columns
     type(physics_ptend), intent(in)       :: ptend_sc     ! ptend for sub-columns

     ! Pointers
     type(physics_buffer_desc), pointer    :: pbuf(:)

     ! Local Variables
     integer :: lchnk, ngrdcol, igrdcol, isubcol, ns, k
     integer, dimension(pcols) :: nsubcol
     real(r8), dimension(pcols*psubcols)       :: weights_packed
     real(r8), dimension(pcols,psubcols)       :: weights
     real(r8), dimension(pcols,psubcols,pverp) :: rc_all, rv_all, rt_all, w_all, thl_all
     real(r8), dimension(pcols,psubcols,pver ) :: s_all, t_all, zm_all, omega_all, pmid_all
     real(r8), dimension(pcols,psubcols)       :: phis_all
     real(r8), dimension(pcols,psubcols,pver ) :: stend, ttend
     real(r8), dimension(pcols,psubcols,pverp) :: thltend, qctend, qvtend

     real(r8), dimension(pcols,psubcols,pver)  :: dz_g, pdel_all, rho
     real(r8), dimension(pcols,psubcols,pverp) :: zi_all
 
     real(r8), pointer, dimension(:,:,:) :: pdf_params_ptr  ! Packed PDF_Params
     type(pdf_parameter), dimension(pverp-top_lev+1) :: pdf_params  ! PDF parameters [units vary]

     integer :: ixcldliq, ixq

     real(r8), dimension(pcols,psubcols,pver ) :: exner

     ! Inputs to lh_microphys_var_covar_driver
     real(r8), dimension(pcols,pverp,psubcols) :: rt_all_clubb, thl_all_clubb, w_all_clubb, &
                                                  qctend_clubb, qvtend_clubb, thltend_clubb
     real(r8), dimension(pcols,pverp-top_lev+1,psubcols) :: height_depndt_weights

     ! Outputs from lh_microphys_var_covar_driver
     real(r8), dimension(:,:), pointer :: rtp2_mc_zt, thlp2_mc_zt, wprtp_mc_zt, &
                                          wpthlp_mc_zt, rtpthlp_mc_zt

     ! pbuf indices
     integer :: &
       rtp2_mc_zt_idx,    &
       thlp2_mc_zt_idx,   &
       wprtp_mc_zt_idx,   &
       wpthlp_mc_zt_idx,  &
       rtpthlp_mc_zt_idx, &
       pdf_params_idx

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. subcol_SILHS_var_covar_src ) return

     lchnk = state_sc%lchnk
     ngrdcol  = state_sc%ngrdcol

     ! Obtain indices
     call cnst_get_ind('Q', ixq)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     rtp2_mc_zt_idx = pbuf_get_index('rtp2_mc_zt')
     thlp2_mc_zt_idx = pbuf_get_index('thlp2_mc_zt')
     wprtp_mc_zt_idx = pbuf_get_index('wprtp_mc_zt')
     wpthlp_mc_zt_idx = pbuf_get_index('wpthlp_mc_zt')
     rtpthlp_mc_zt_idx = pbuf_get_index('rtpthlp_mc_zt')
     pdf_params_idx = pbuf_get_index('PDF_PARAMS')

     ! Obtain pbuf fields for computation
     call pbuf_get_field(pbuf, pdf_params_idx, pdf_params_ptr)

     ! Obtain pbuf fields for output
     call pbuf_get_field(pbuf, rtp2_mc_zt_idx, rtp2_mc_zt)
     call pbuf_get_field(pbuf, thlp2_mc_zt_idx, thlp2_mc_zt)
     call pbuf_get_field(pbuf, wprtp_mc_zt_idx, wprtp_mc_zt)
     call pbuf_get_field(pbuf, wpthlp_mc_zt_idx, wpthlp_mc_zt)
     call pbuf_get_field(pbuf, rtpthlp_mc_zt_idx, rtpthlp_mc_zt)

     ! Unpack needed tendencies from subcolumn ptends
     call subcol_unpack(lchnk, ptend_sc%s(:,:), stend, fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixcldliq), qctend(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixq), qvtend(:,:,1:pver), fillvalue)

     ! Unpack sample point values from subcolumn states
     call subcol_unpack(lchnk, state_sc%q(:,:,ixcldliq), rc_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%q(:,:,ixq),      rv_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%omega(:,:),      omega_all (:,:,:),  fillvalue)
     call subcol_unpack(lchnk, state_sc%s(:,:),          s_all,              fillvalue)
     call subcol_unpack(lchnk, state_sc%zm,              zm_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%phis,            phis_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%zi,              zi_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%pdel,            pdel_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%pmid,            pmid_all,           fillvalue)

     ! Initialize fields to fillvalue.
     rt_all  = fillvalue
     thl_all = fillvalue
     w_all   = fillvalue
     qctend  = fillvalue
     qvtend  = fillvalue
     thltend = fillvalue

     ! How many subcolumns in each column?
     call subcol_get_nsubcol(lchnk, nsubcol)

     do igrdcol = 1, ngrdcol
        do isubcol = 1, nsubcol(igrdcol)

           rt_all(igrdcol,isubcol,top_lev:pver) = rc_all(igrdcol,isubcol,top_lev:pver) &
                                                  + rv_all(igrdcol,isubcol,top_lev:pver)

           ! Compute dry static density on CLUBB vertical grid
           do k = top_lev, pver
              dz_g(igrdcol,isubcol,k) = zi_all(igrdcol,isubcol,k) - zi_all(igrdcol,isubcol,k+1) ! thickness
              rho(igrdcol,isubcol,k) = (1._r8/gravit)*pdel_all(igrdcol,isubcol,k)/dz_g(igrdcol,isubcol,k)
           enddo

           ! Compute w from omega
           w_all(igrdcol,isubcol,top_lev:pver) = -omega_all(igrdcol,isubcol,top_lev:pver) &
                                                  / ( rho(igrdcol,isubcol,top_lev:pver) * gravit )

           ! Convert stend and s_all to ttend and t_all
           !  Note 1: With subcolumns, cpair is truly a constant (I think).
           !  Note 2: For tendencies, the extra terns zm and phis should
           !          not be included in the calculation.
           ttend(igrdcol,isubcol,top_lev:pver) = stend(igrdcol,isubcol,top_lev:pver) / cpair

           do k = top_lev, pver
              t_all(igrdcol,isubcol,k) = ( s_all(igrdcol,isubcol,k) &
                                           - gravit * zm_all(igrdcol,isubcol,k) &
                                           - phis_all(igrdcol,isubcol) ) / cpair
           enddo ! k = 1, pver

           ! This formula is taken from earlier in this file.
           exner(igrdcol,isubcol,top_lev:pver) &
           = ( pmid_all(igrdcol,isubcol,top_lev:pver) / p0_clubb )**(rair/cpair)

           ! Note: all tendencies or all means should be used in the call to
           !       T_in_K2thlm_api (with the exception of exner)
           do k = top_lev, pver
              thltend(igrdcol,isubcol,k) &
              = T_in_K2thlm_api( ttend(igrdcol,isubcol,k), exner(igrdcol,isubcol,k), &
                                 qctend(igrdcol,isubcol,k) )
              thl_all(igrdcol,isubcol,k) &
              = T_in_K2thlm_api( t_all(igrdcol,isubcol,k), exner(igrdcol,isubcol,k), &
                                 rc_all(igrdcol,isubcol,k) )
           enddo ! k = 1, pver

           ! Add ghost points
           rt_all (igrdcol,isubcol,pverp) = rt_all (igrdcol,isubcol,pver)
           thl_all(igrdcol,isubcol,pverp) = thl_all(igrdcol,isubcol,pver)
           w_all  (igrdcol,isubcol,pverp) = w_all  (igrdcol,isubcol,pver)
           qctend (igrdcol,isubcol,pverp) = qctend (igrdcol,isubcol,pver)
           qvtend (igrdcol,isubcol,pverp) = qvtend (igrdcol,isubcol,pver)
           thltend(igrdcol,isubcol,pverp) = thltend(igrdcol,isubcol,pver)

           ! Flip inputs to CLUBB's grid. Note the dimension ordering change.
           rt_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( rt_all(igrdcol,isubcol,1:pverp) )
           thl_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( thl_all(igrdcol,isubcol,1:pverp) )
           w_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( w_all(igrdcol,isubcol,1:pverp) )
           qctend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( qctend(igrdcol,isubcol,1:pverp) )
           qvtend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( qvtend(igrdcol,isubcol,1:pverp) )
           thltend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( thltend(igrdcol,isubcol,1:pverp) )

        enddo ! isubcol = 1, nsubcol(igrdcol)
     enddo ! igrdcol = 1, ngrdcol

     ! Obtain weights
     call subcol_get_weight(lchnk, weights_packed)
     call subcol_unpack(lchnk, weights_packed, weights, fillvalue)

     ! Call lh_microphys_var_covar_driver for each column
     do igrdcol=1, ngrdcol
       ns = nsubcol(igrdcol)

       ! Obtain PDF parameters
       ! The PDF parameters are passed out of CLUBB and into SILHS without
       ! being used at all in the rest of the host model code.  The arrays
       ! aren't flipped for the PDF parameters, and they don't need to be.
       call unpack_pdf_params_api(pdf_params_ptr(igrdcol,1:pverp-top_lev+1,:), pverp-top_lev+1, pdf_params)

       ! This code assumes that the weights are height independent.
       ! It will have to change once the weights vary with altitude!
       ! I'm not sure whether the grid will need to be flipped.
       do k = 1, pverp-top_lev+1
          height_depndt_weights(igrdcol,k,1:ns) = weights(igrdcol,1:ns)
       end do

       ! Make the call!!!!!
       call lh_microphys_var_covar_driver_api &
            ( pverp-top_lev+1, ns, ztodt, height_depndt_weights(igrdcol,1:pverp-top_lev+1,1:ns), pdf_params, &
              rt_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), thl_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              w_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), qctend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              qvtend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), thltend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              rtp2_mc_zt(igrdcol,1:pverp-top_lev+1), thlp2_mc_zt(igrdcol,1:pverp-top_lev+1), &
              wprtp_mc_zt(igrdcol,1:pverp-top_lev+1), wpthlp_mc_zt(igrdcol,1:pverp-top_lev+1), &
              rtpthlp_mc_zt(igrdcol,1:pverp-top_lev+1) )

       ! The *_mc_zt microphysics tendencies are passed out of SILHS and back
       ! to CLUBB without being used at all in the rest of the host model code.
       ! The arrays aren't flipped for the *_mc_zt microphysics tendencies, and
       ! they don't need to be.

       ! CLUBB used pverp vertical levels, but SILHS only uses
       ! pverp - top_lev + 1 vertical levels.
       ! Fill the upper levels with 0s when necessary.
       if ( pverp > pverp-top_lev+1 ) then
          rtp2_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          thlp2_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          wprtp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          wpthlp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          rtpthlp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
       endif ! pverp > pverp-top_lev+1

     enddo ! igrdcol = 1, ngrdcol

     return
   end subroutine subcol_SILHS_var_covar_driver

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !

   subroutine subcol_SILHS_massless_droplet_destroyer &
              ( ztodt, state, &
                ptend )

     ! This subroutine eradicates cloud droplets in grid boxes with no cloud
     ! mass. This code is actually not at all related to SILHS. It should be
     ! moved out of this cage they call "subcol_SILHS.F90" someday.
     !
     ! This code is experimental!!

     use micro_mg_utils, only: qsmall
     use ref_pres,       only: top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     real(r8), intent(in)                  :: ztodt     ! model time increment
     type(physics_state), intent(in)       :: state     ! state for columns

     ! Input/Output Variables
     type(physics_ptend), intent(inout)    :: ptend     ! ptend for columns

     ! Local Variables
     integer :: icol, k

     integer :: ixcldliq, ixnumliq

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. subcol_SILHS_destroy_massless_droplets ) return

     ! Indices!
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call cnst_get_ind('NUMLIQ', ixnumliq)

     ! These "labels" in loops are really cool. We should start doing this in
     ! CLUBB.
     col_loop: do icol=1, state%ncol
       ! If updated qc (after microphysics) is zero, then ensure updated nc is also zero!!
       vert_loop: do k = top_lev, pver
         if ( state%q(icol,k,ixcldliq) + (ztodt*ptend%q(icol,k,ixcldliq)) < qsmall ) then
           ptend%lq(ixnumliq) = .true. ! This is probably already true, but it doesn't
                                       ! hurt to set it.
           ptend%q(icol,k,ixnumliq) = -(state%q(icol,k,ixnumliq) / ztodt)
         end if
       end do vert_loop
     end do col_loop

     return
   end subroutine subcol_SILHS_massless_droplet_destroyer

   !============================================================================
   subroutine subcol_SILHS_fill_holes_conserv( state, dt, ptend, pbuf )

     ! The William F. Buckley Jr. Conservative Hole Filler.

     ! Description:
     ! Stops holes from forming in a hydrometeor mixing ratio by reducing the
     ! microphysics tendency of that hydrometeor mixing ratio which would
     ! otherwise cause that hydrometeor mixing ratio to have a negative value
     ! once the microphysics tendency is applied.  This code is used to prevent
     ! holes in water mass, not number concentration.
     !
     ! This subroutine is called after microphysics has completed and after
     ! microphysics fields from subcolumns have been averaged back to grid
     ! columns, but before the grid-column microphysics tendencies have been
     ! applied in physics_update.  This code is meant for use with the SILHS
     ! subcolumn approach.  This code needs to be applied to grid columns, not
     ! subcolumns.
     !
     ! This code adjusts the tendencies (ptend) before they are used to update
     ! the grid mean fields (state variables).
     !
     ! The column-integrated total water needs to be conserved during
     ! microphysics.  The conserved amount includes the amount of water that
     ! precipitated to the ground from sedimentation during microphysics.
     ! The conservation equation for each grid column is:
     !
     ! SUM(k=top_lev:pver) ( rv_start(k) + rc_start(k) + rr_start(k)
     !                       + ri_start(k) + rs_start(k) ) * pdel(k) / g
     ! = SUM(k=top_lev:pver) ( rv(k) + rc(k) + rr(k) + ri(k) + rs(k) )
     !                       * pdel(k) / g
     !   + prect * dt * 1000;
     !
     ! where rv_start, rc_start, rr_start, ri_start, and rs_start are water
     ! vapor, cloud water, rain water, cloud ice, and snow mixing ratios before
     ! microphysics is called; rv, rc, rr, ri, and rs are water vapor, cloud
     ! water, rain water, cloud ice, and snow mixing ratios after being updated
     ! by microphysics; pdel is the pressure difference between vertical levels,
     ! g is gravity, and prect * dt * 1000 is the total amount of water (from
     ! all precipitating hydrometeors) that sedimented to the ground during
     ! microphysics (dt is the timestep used for microphysics).  The units of
     ! column-integrated total water are kg (water) / m^2.
     !
     ! All the updated hydrometeor fields are related to the hydrometeor fields
     ! at the start by:
     !
     ! rv(k) = rv_start(k) + rv_tend(k) * dt;
     ! rc(k) = rc_start(k) + rc_tend(k) * dt;
     ! rr(k) = rr_start(k) + rr_tend(k) * dt;
     ! ri(k) = ri_start(k) + ri_tend(k) * dt; and
     ! rs(k) = rs_start(k) + rs_tend(k) * dt;
     !
     ! where rv_tend, rc_tend, rr_tend, ri_tend, and rs_tend are water vapor,
     ! cloud water, rain water, cloud ice, and snow mixing ratio tendencies
     ! from microphysics, which includes the sum of microphysics process rates
     ! and sedimentation.  When these equations are applied to the equation
     ! for column-integrated total water, that equation becomes:
     !
     ! SUM(k=top_lev:pver) ( rv_tend(k) + rc_tend(k) + rr_tend(k)
     !                       + ri_tend(k) + rs_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! As stated above, the hydrometeor tendencies are the sum of tendencies
     ! from microphysics process rates and tendencies from sedimentation:
     !
     ! rv_tend(k) = rv_mc_tend(k);
     ! rc_tend(k) = rc_mc_tend(k) + rc_sed_tend(k);
     ! rr_tend(k) = rr_mc_tend(k) + rr_sed_tend(k);
     ! ri_tend(k) = ri_mc_tend(k) + ri_sed_tend(k); and
     ! rs_tend(k) = rs_mc_tend(k) + rs_sed_tend(k);
     !
     ! where rv_mc_tend, rc_mc_tend, rr_mc_tend, ri_mc_tend, and rs_mc_tend are
     ! the tendencies of water vapor, cloud water, rain water, cloud ice, and
     ! snow from microphysics process rates, and rc_sed_tend, rr_sed_tend,
     ! ri_sed_tend, and rs_sed_tend are the tendencies of cloud water,
     ! rain water, cloud ice, and snow from sedimentation.  When these equations
     ! are applied to the equation for column-integrated total water, that
     ! equation becomes:
     !
     ! SUM(k=top_lev:pver) ( rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     !                       + ri_mc_tend(k) + rs_mc_tend(k) )
     !                     * dt * pdel(k) / g
     ! + SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) + ri_sed_tend(k)
     !                         + rs_sed_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! At any vertical level, the tendencies from microphysics process rates
     ! (mc_tend variables) must balance:
     !
     ! rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     ! + ri_mc_tend(k) + rs_mc_tend(k) = 0; for all k from top_lev to pver.
     !
     ! The column-integrated total water equation can be applied to
     ! sedimentation:
     !
     ! SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) + ri_sed_tend(k)
     !                       + rs_sed_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! The total precipitation rate, prect, can be split into liquid
     ! precipitation rate, precl, and frozen precipitation rate, preci:
     !
     ! prect = precl + preci.
     !
     ! The microphysics code outputs prect and preci, so precl can be calculated
     ! by precl = prect - preci.  The column-integrated total water equation can
     ! be split into:
     !
     ! SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
     !                     * dt * pdel(k) / g
     ! + precl * dt * 1000 = 0; and
     !
     ! SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
     !                     * dt * pdel(k) / g
     ! + preci * dt * 1000 = 0.
     !
     ! Overall, the conservation methods used in this subroutine are:
     !
     ! 1) When adjusting the tendencies from microphysics process rates,
     !    conserve:
     !
     !    rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     !    + ri_mc_tend(k) + rs_mc_tend(k) = 0; for all k from top_lev to pver.
     !
     ! 2) When adjusting the tendencies from microphysics process rates, adjust
     !    dry static energy appropriately.  The change in dry static energy
     !    is necessary because of phase changes.  This "puts back" the extra dry
     !    static energy that was "taken out" when an excessive phase-changing
     !    process rate was produced by microphysics.
     !
     ! 3) When adjusting the hydrometeor tendency from sedimentation of a
     !    liquid hydrometeor (cloud water or rain water), conserve:
     !
     !    SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + precl * dt * 1000 = 0.
     !
     ! 4) When adjusting the hydrometeor tendency from sedimentation of a
     !    frozen hydrometeor (cloud ice or snow), conserve:
     !
     !    SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + preci * dt * 1000 = 0.

     !----------------------------------------------------------------------

     use physics_buffer, only: &
         physics_buffer_desc, &
         pbuf_get_field

     use ppgrid, only: &
         pcols

     use constituents, only: &
         qmin

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     type(physics_state), intent(in) :: state     ! Physics state variables
     real(r8), intent(in) :: dt                   ! Time step duration

     ! Input/Output Variables
     type(physics_ptend),  intent(inout) :: ptend  ! Parameterization tendencies
     type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer

     ! Local Variables
     real(r8), dimension(pcols,pver) :: &
       rv_start, & ! Water vapor mixing ratio at start of microphysics  [kg/kg]
       rc_start, & ! Cloud water mixing ratio at start of microphysics  [kg/kg]
       rr_start, & ! Rain water mixing ratio at start of microphysics   [kg/kg]
       ri_start, & ! Cloud ice mixing ratio at start of microphysics    [kg/kg]
       rs_start    ! Snow mixing ratio at start of microphysics         [kg/kg]

     real(r8), dimension(pcols,pver) :: &
       rv_tend, & ! Water vapor mixing ratio tendency  [kg/kg/s]
       rc_tend, & ! Cloud water mixing ratio tendency  [kg/kg/s]
       rr_tend, & ! Rain water mixing ratio tendency   [kg/kg/s]
       ri_tend, & ! Cloud ice mixing ratio tendency    [kg/kg/s]
       rs_tend, & ! Snow mixing ratio tendency         [kg/kg/s]
       stend      ! Dry static energy tendency         [J/kg/s]

     real(r8), dimension(:), pointer :: &
       prect, & ! Total precipitation rate (surface)        [m/s]
       preci    ! Ice-phase precipitation rate (surface)    [m/s]

     real(r8), dimension(:,:), pointer :: &
       rc_sed_tend, & ! Mean cloud water sedimentation tendency    [kg/kg/s]
       rr_sed_tend, & ! Mean rain water sedimentation tendency     [kg/kg/s]
       ri_sed_tend, & ! Mean cloud ice sedimentation tendency      [kg/kg/s]
       rs_sed_tend, & ! Mean snow sedimentation tendency           [kg/kg/s]
       vtrmc,       & ! Mean cloud water sedimentation velocity    [m/s]
       umr,         & ! Mean rain water sedimentation velocity     [m/s]
       vtrmi,       & ! Mean cloud ice sedimentation velocity      [m/s]
       ums            ! Mean snow sedimentation velocity           [m/s]

     real(r8), dimension(pcols,pver) :: &
       rv_mc_tend, & ! Water vapor mixing ratio microphysics tendency  [kg/kg/s]
       rc_mc_tend, & ! Cloud water mixing ratio microphysics tendency  [kg/kg/s]
       rr_mc_tend, & ! Rain water mixing ratio microphysics tendency   [kg/kg/s]
       ri_mc_tend, & ! Cloud ice mixing ratio microphysics tendency    [kg/kg/s]
       rs_mc_tend    ! Snow mixing ratio microphysics tendency         [kg/kg/s]

     real(r8) :: &
       rv_curr, & ! Current water vapor mixing ratio    [kg/kg]
       rc_curr, & ! Current cloud water mixing ratio    [kg/kg]
       rr_curr, & ! Current rain water mixing ratio     [kg/kg]
       ri_curr, & ! Current cloud ice mixing ratio      [kg/kg]
       rs_curr    ! Current snow mixing ratio           [kg/kg]

     logical :: &
       l_pos_rv_mc_tend, & ! Flag for positive water vapor mixing ratio mc tend.
       l_pos_rc_mc_tend, & ! Flag for positive cloud water mixing ratio mc tend.
       l_pos_rr_mc_tend, & ! Flag for positive rain water mixing ratio mc tend.
       l_pos_ri_mc_tend, & ! Flag for positive cloud ice mixing ratio mc tend.
       l_pos_rs_mc_tend    ! Flag for positive snow mixing ratio mc tend.

     real(r8) :: &
       mc_tend_max_mag,     & ! Max. allowable mag. of (neg.) mc tend [kg/kg/s]
       mc_tend_correction,  & ! Amnt. correction necessary to mc tend [kg/kg/s]
       total_mc_positive,   & ! Total of all positive mc tendencies   [kg/kg/s]
       mc_correction_ratio    ! Ratio: mc_tend_correction/total_mc_positive [-]

     real(r8), dimension(pcols) :: &
       precl    ! Liquid-phase precipitation rate (surface)        [m/s]

     ! Budgeting terms for hole filling.
     ! These variables are for use in stats output.
     real(r8), dimension(pcols,pver) :: &
       rv_hf_tend, & ! Water vapor mixing ratio hole-filling tendency  [kg/kg/s]
       rc_hf_tend, & ! Cloud water mixing ratio hole-filling tendency  [kg/kg/s]
       rr_hf_tend, & ! Rain water mixing ratio hole-filling tendency   [kg/kg/s]
       ri_hf_tend, & ! Cloud ice mixing ratio hole-filling tendency    [kg/kg/s]
       rs_hf_tend, & ! Snow mixing ratio hole-filling tendency         [kg/kg/s]
       s_hf_tend     ! Dry static energy hole-filling tendency         [J/kg/s]

     integer :: ixcldice, ixcldliq, ixrain, ixsnow ! Hydrometeor indices

     integer :: ncol  ! Number of grid columns

     integer :: icol, k  ! Loop indices

     ! Flag to perform hole filling after the original sedimentation tendency
     ! is added back on to the new microphysics process tendency.  This calls
     ! the sedimentation hole filler.
     logical, parameter :: &
       l_sed_hole_fill = .true.

     logical, parameter :: &
       l_check_conservation = .true. ! Flag to perform water conservation check

     ! Vertically-integrated grand total water (rv+rc+rr+ri+rs)    [kg/m^2]
     real(r8), dimension(pcols) :: &
       grand_total_water_column_start,  & ! Column integral at start
       grand_total_water_column_finish    ! Column integral at finish

     ! Vertically-integrated total water energy    [J/m^2]
     real(r8), dimension(pcols) :: &
       total_energy_column_start,  & ! Column integral at start
       total_energy_column_finish    ! Column integral at finish

     real(r8), dimension(pcols) :: &
       tot_water_rel_err,  & ! Relative error: vert-integrated grand total water
       tot_energy_rel_err    ! Relative error: vert-integrated total energy

     real(r8), parameter :: &
       err_thresh = 1.0e-14_r8  ! Threshold of relative error


     ! Get indices for hydrometeor fields.
     call cnst_get_ind('CLDICE', ixcldice)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call cnst_get_ind('RAINQM', ixrain, abort=.false.)
     call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)

     ! Get the number of grid columns.
     ncol = state%ncol

     ! Get fields from the pbuf.
     call pbuf_get_field(pbuf, prec_pcw_idx, prect)
     call pbuf_get_field(pbuf, snow_pcw_idx, preci)
     call pbuf_get_field(pbuf, qcsedten_idx, rc_sed_tend)
     call pbuf_get_field(pbuf, qrsedten_idx, rr_sed_tend)
     call pbuf_get_field(pbuf, qisedten_idx, ri_sed_tend)
     call pbuf_get_field(pbuf, qssedten_idx, rs_sed_tend)
     call pbuf_get_field(pbuf, vtrmc_idx, vtrmc)
     call pbuf_get_field(pbuf, umr_idx, umr)
     call pbuf_get_field(pbuf, vtrmi_idx, vtrmi)
     call pbuf_get_field(pbuf, ums_idx, ums)

     ! Calculate liquid precipitation rate (precl) from the total precipitation
     ! rate (prect) and the frozen preciptation rate (preci).  This should never
     ! be negative, but just to be safe, threshold at 0.
     precl = max( prect - preci, 0.0_r8 )

     ! Perform total water and total energy conservation checks.
     if ( l_check_conservation ) then

        ! Calculate total water in each column.
        ! This calculation is the vertically-integrated grand total water
        ! in each grid column before microphysics began.
        do icol = 1, ncol
           grand_total_water_column_start(icol) = 0.0_r8
           do k = top_lev, pver
              grand_total_water_column_start(icol) &
              = grand_total_water_column_start(icol) &
                + ( state%q(icol,k,1) + state%q(icol,k,ixcldliq) &
                    + state%q(icol,k,ixcldice) ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 grand_total_water_column_start(icol) &
                 = grand_total_water_column_start(icol) &
                   + state%q(icol,k,ixrain) * state%pdel(icol,k) / gravit
              endif
              if ( ixsnow > 0 ) then
                 grand_total_water_column_start(icol) &
                 = grand_total_water_column_start(icol) &
                   + state%q(icol,k,ixsnow) * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
        enddo ! icol = 1, ncol

        ! Calculate total energy in each column.
        ! This calculation is the vertically-integrated total energy in each
        ! grid column before microphysics began.  Since, the microphysics code
        ! does not directly change kinetic energy, 0.5 * ( u^2 + v^2 ), it can
        ! be skipped as part of the energy conservation check.
        !do icol = 1, ncol
        !   total_energy_column_start(icol) = 0.0_r8
        !   do k = top_lev, pver
        !      total_energy_column_start(icol) &
        !      = total_energy_column_start(icol) &
        !        + ( state%s(icol,k) &
        !            + ( latvap + latice ) * state%q(icol,k,1) &
        !            + latice * state%q(icol,k,ixcldliq) ) &
        !          * state%pdel(icol,k) / gravit
        !      if ( ixrain > 0 ) then
        !         total_energy_column_start(icol) &
        !         = total_energy_column_start(icol) &
        !           + latice * state%q(icol,k,ixrain) &
        !             * state%pdel(icol,k) / gravit
        !      endif
        !   enddo ! k = top_lev, pver
        !enddo ! icol = 1, ncol
        ! This calculation is the vertically-integrated total energy in each
        ! grid column after microphysics, but at the start of hole-filling.
        do icol = 1, ncol
           total_energy_column_start(icol) = 0.0_r8
           do k = top_lev, pver
              total_energy_column_start(icol) &
              = total_energy_column_start(icol) &
                + ( state%s(icol,k) + ptend%s(icol,k) * dt &
                    + ( latvap + latice ) &
                      * ( state%q(icol,k,1) + ptend%q(icol,k,1) * dt ) &
                    + latice * ( state%q(icol,k,ixcldliq) &
                                 + ptend%q(icol,k,ixcldliq) * dt ) ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 total_energy_column_start(icol) &
                 = total_energy_column_start(icol) &
                   + latice * ( state%q(icol,k,ixrain) &
                                + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
              total_energy_column_start(icol) &
              = total_energy_column_start(icol) &
                + latice * precl(icol) * dt * 1000.0_r8
           enddo ! k = top_lev, pver
        enddo ! icol = 1, ncol

     endif ! l_check_conservation

     ! The fields within state haven't been updated yet, since this is before
     ! the call to physics_update.
     rv_start = state%q(:,:,1)
     rc_start = state%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_start = state%q(:,:,ixrain)
     endif
     ri_start = state%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_start = state%q(:,:,ixsnow)
     endif

     ! Unpack the current total tendencies for hydrometeor mixing ratio fields.
     rv_tend = ptend%q(:,:,1)
     rc_tend = ptend%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_tend = ptend%q(:,:,ixrain)
     endif
     ri_tend = ptend%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_tend = ptend%q(:,:,ixsnow)
     endif

     ! Unpack the current tendency for dry static energy.
     stend = ptend%s

     ! The total hydrometeor tendencies are the sum of microphysics process
     ! rates and sedimentation rates.  Calculate the microphysics process
     ! tendencies by subtracting the sedimentation tendencies from the overall
     ! tendencies.
     rv_mc_tend = rv_tend
     rc_mc_tend = rc_tend - rc_sed_tend
     if ( ixrain > 0 ) then
        rr_mc_tend = rr_tend - rr_sed_tend
     endif
     ri_mc_tend = ri_tend - ri_sed_tend
     if ( ixsnow > 0 ) then
        rs_mc_tend = rs_tend - rs_sed_tend
     endif

     ! Loop over all columns, performing any tendency adjustments one column
     ! at a time.
     do icol = 1, ncol

        ! Loop over all vertical levels, performing any microphysics process
        ! tendency adjustments one level at a time.
        do k = top_lev, pver

           ! Find which hydrometeors have positive microphysics process
           ! tendencies at this level.
           if ( rv_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_rv_mc_tend = .true.
           else
              l_pos_rv_mc_tend = .false.
           endif
           if ( rc_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_rc_mc_tend = .true.
           else
              l_pos_rc_mc_tend = .false.
           endif
           if ( ixrain > 0 ) then
              if ( rr_mc_tend(icol,k) >= 0.0_r8 ) then
                 l_pos_rr_mc_tend = .true.
              else
                 l_pos_rr_mc_tend = .false.
              endif
           endif
           if ( ri_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_ri_mc_tend = .true.
           else
              l_pos_ri_mc_tend = .false.
           endif
           if ( ixsnow > 0 ) then
              if ( rs_mc_tend(icol,k) >= 0.0_r8 ) then
                 l_pos_rs_mc_tend = .true.
              else
                 l_pos_rs_mc_tend = .false.
              endif
           endif

           !!! Check for holes in water vapor mixing ratio
           if ( .not. l_pos_rv_mc_tend ) then

              ! Calculate the water vapor mixing ratio as it would be with the
              ! current microphysics process tendency.
              rv_curr = rv_start(icol,k) + rv_mc_tend(icol,k) * dt

              if ( rv_curr < qmin(1) ) then

                 ! Microphysics processes are causing a hole in water vapor
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) water
                 ! vapor microphysics process tendency.
                 mc_tend_max_mag = ( qmin(1) - rv_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the water vapor mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rv_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principal, this should never be greater than 1, but this
                 ! is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative water vapor mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of water vapor).
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latvap * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latvap * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - ( latvap + latice ) &
                        * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to water vapor cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - ( latvap + latice ) &
                        * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new water vapor mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rv_mc_tend(icol,k) &
                 = rv_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rv_curr < qmin(1)

           endif ! .not. l_pos_rv_mc_tend

           !!! Check for holes in cloud water mixing ratio
           if ( .not. l_pos_rc_mc_tend ) then

              ! Calculate the cloud water mixing ratio as it would be with the
              ! current microphysics process tendency.
              rc_curr = rc_start(icol,k) + rc_mc_tend(icol,k) * dt

              if ( rc_curr < qmin(ixcldliq) ) then

                 ! Microphysics processes are causing a hole in cloud water
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) cloud
                 ! water microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixcldliq) - rc_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the cloud water mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rc_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principal, this should never be greater than 1, but this
                 ! is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative cloud water mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of cloud water).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to cloud water heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latvap * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to cloud water does not change
                    ! dry static energy.
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to cloud water cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to cloud water cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new cloud water mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rc_mc_tend(icol,k) &
                 = rc_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rc_curr < qmin(ixcldliq)

           endif ! .not. l_pos_rc_mc_tend

           !!! Check for holes in rain water mixing ratio
           if ( ixrain > 0 .and. ( .not. l_pos_rr_mc_tend ) ) then

              ! Calculate the rain water mixing ratio as it would be with the
              ! current microphysics process tendency.
              rr_curr = rr_start(icol,k) + rr_mc_tend(icol,k) * dt

              if ( rr_curr < qmin(ixrain) ) then

                 ! Microphysics processes are causing a hole in rain water
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) rain
                 ! water microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixrain) - rr_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the rain water mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rr_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principal, this should never be greater than 1, but this
                 ! is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative rain water mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of rain water).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to rain water heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latvap * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to rain water does not change
                    ! dry static energy.
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to rain water cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to rain water cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new rain water mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rr_mc_tend(icol,k) &
                 = rr_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rr_curr < qmin(ixrain)

           endif ! ixrain > 0 .and. ( .not. l_pos_rr_mc_tend )

           !!! Check for holes in cloud ice mixing ratio
           if ( .not. l_pos_ri_mc_tend ) then

              ! Calculate the cloud ice mixing ratio as it would be with the
              ! current microphysics process tendency.
              ri_curr = ri_start(icol,k) + ri_mc_tend(icol,k) * dt

              if ( ri_curr < qmin(ixcldice) ) then

                 ! Microphysics processes are causing a hole in cloud ice
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) cloud
                 ! ice microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixcldice) - ri_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the cloud ice mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - ri_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principal, this should never be greater than 1, but this
                 ! is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative cloud ice mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of cloud ice).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + ( latvap + latice ) &
                        * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to cloud ice does not change dry static
                    ! energy.
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new cloud ice mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 ri_mc_tend(icol,k) &
                 = ri_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! ri_curr < qmin(ixcldice)

           endif ! .not. l_pos_ri_mc_tend

           !!! Check for holes in snow mixing ratio
           if ( ixsnow > 0 .and. ( .not. l_pos_rs_mc_tend ) ) then

              ! Calculate the snow mixing ratio as it would be with the
              ! current microphysics process tendency.
              rs_curr = rs_start(icol,k) + rs_mc_tend(icol,k) * dt

              if ( rs_curr < qmin(ixsnow) ) then

                 ! Microphysics processes are causing a hole in snow mixing
                 ! ratio.

                 ! Calculate the maximum allowable magnitude of (negative) snow
                 ! microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixsnow) - rs_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the snow mixing ratio microphysics process tendency.
                 ! This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rs_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principal, this should never be greater than 1, but this
                 ! is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative snow mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of snow).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + ( latvap + latice ) &
                        * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to snow does not change dry static
                    ! energy.
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new snow mixing ratio microphysics process
                 ! tendency.  This should be equal to the maximum magnitude
                 ! (negative) amount allowed, mc_tend_max_mag.
                 rs_mc_tend(icol,k) &
                 = rs_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rs_curr < qmin(ixsnow)

           endif ! ixsnow > 0 .and. ( .not. l_pos_rs_mc_tend )

        enddo ! k = 1, pver

     enddo ! icol = 1, ncol

     ! Calculate the new overall tendencies by adding the sedimentation
     ! tendencies back onto the new microphysics process tendencies.
     rv_tend = rv_mc_tend
     rc_tend = rc_mc_tend + rc_sed_tend
     if ( ixrain > 0 ) then
        rr_tend = rr_mc_tend + rr_sed_tend
     endif
     ri_tend = ri_mc_tend + ri_sed_tend
     if ( ixsnow > 0 ) then
        rs_tend = rs_mc_tend + rs_sed_tend
     endif

     ! Now that the original sedimentation tendency has been added to the
     ! new microphysics process tendency, the new total microphysics tendency
     ! can still cause holes to form.  These holes can be filled conservatively
     ! using the sedimentation hole filler.
     if ( l_sed_hole_fill ) then

        ! Call the sedimentation hole filler for cloud water mixing ratio.
        ! This can update rc_tend and precl.
        call fill_holes_sedimentation( dt, ncol, rc_start, state%pdel, &
                                       vtrmc, state%zi, qmin(ixcldliq), &
                                       rc_tend, precl )

        ! Call the sedimentation hole filler for rain water mixing ratio.
        ! This can update rr_tend and precl.
        if ( ixrain > 0 ) then
           call fill_holes_sedimentation( dt, ncol, rr_start, state%pdel, &
                                          umr, state%zi, qmin(ixrain), &
                                          rr_tend, precl )
        endif ! ixrain > 0

        ! Call the sedimentation hole filler for cloud ice mixing ratio.
        ! This can update ri_tend and preci.
        call fill_holes_sedimentation( dt, ncol, ri_start, state%pdel, &
                                       vtrmi, state%zi, qmin(ixcldice), &
                                       ri_tend, preci )

        ! Call the sedimentation hole filler for snow mixing ratio.
        ! This can update rs_tend and preci.
        if ( ixsnow > 0 ) then
           call fill_holes_sedimentation( dt, ncol, rs_start, state%pdel, &
                                          ums, state%zi, qmin(ixsnow), &
                                          rs_tend, preci )
        endif ! ixsnow > 0

        ! Update the total precipitation rate (prect) from the updated liquid
        ! precipitation rate (precl) and the updated frozen preciptation rate
        ! (preci).
        prect = precl + preci

     endif ! l_sed_hole_fill

     ! The updated total microphysics tendencies after hole filling have not
     ! been used to update ptend yet, so record the budget terms for hole
     ! filling first.
     rv_hf_tend = rv_tend - ptend%q(:,:,1)
     rc_hf_tend = rc_tend - ptend%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_hf_tend = rr_tend - ptend%q(:,:,ixrain)
     endif ! ixrain > 0
     ri_hf_tend = ri_tend - ptend%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_hf_tend = rs_tend - ptend%q(:,:,ixsnow)
     endif ! ixsnow > 0

     ! The updated dry static energy tendency after hole filling has not been
     ! used to update ptend yet, so record the budget term for hole filling
     ! first.
     s_hf_tend = stend - ptend%s

     ! Pack the current total tendencies for hydrometeor mixing ratio fields.
     ptend%q(:,:,1) = rv_tend
     ptend%q(:,:,ixcldliq) = rc_tend
     if ( ixrain > 0 ) then
        ptend%q(:,:,ixrain) = rr_tend
     endif
     ptend%q(:,:,ixcldice) = ri_tend
     if ( ixsnow > 0 ) then
        ptend%q(:,:,ixsnow) = rs_tend
     endif

     ! Pack the current tendency for dry static energy.
     ptend%s = stend

     ! Perform total water and total energy conservation checks.
     if ( l_check_conservation ) then

        ! Calculate total water in each grid column.
        ! This calculation is the vertically-integrated grand total water
        ! in each grid column updated for all microphysics and hole filling.
        ! This includes the amount that precipitated to the surface.
        do icol = 1, ncol
           grand_total_water_column_finish(icol) = 0.0_r8
           do k = top_lev, pver
              grand_total_water_column_finish(icol) &
              = grand_total_water_column_finish(icol) &
                + ( state%q(icol,k,1) &
                    + ptend%q(icol,k,1) * dt &
                    + state%q(icol,k,ixcldliq) &
                    + ptend%q(icol,k,ixcldliq) * dt &
                    + state%q(icol,k,ixcldice) &
                    + ptend%q(icol,k,ixcldice) * dt ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 grand_total_water_column_finish(icol) &
                 = grand_total_water_column_finish(icol) &
                   + ( state%q(icol,k,ixrain) + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
              if ( ixsnow > 0 ) then
                 grand_total_water_column_finish(icol) &
                 = grand_total_water_column_finish(icol) &
                   + ( state%q(icol,k,ixsnow) + ptend%q(icol,k,ixsnow) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
           grand_total_water_column_finish(icol) &
           = grand_total_water_column_finish(icol) &
             + prect(icol) * dt * 1000.0_r8
        enddo ! icol = 1, ncol

        ! Calculate total energy in each column.
        ! This calculation is the vertically-integrated total energy in each
        ! grid column updated for all microphysics and hole filling.  This
        ! includes the amount that precipitated to the surface.  Since, the
        ! microphysics code does not directly change kinetic energy,
        ! 0.5 * ( u^2 + v^2 ), it can be skipped as part of the energy
        ! conservation check.
        do icol = 1, ncol
           total_energy_column_finish(icol) = 0.0_r8
           do k = top_lev, pver
              total_energy_column_finish(icol) &
              = total_energy_column_finish(icol) &
                + ( state%s(icol,k) + ptend%s(icol,k) * dt &
                    + ( latvap + latice ) &
                      * ( state%q(icol,k,1) + ptend%q(icol,k,1) * dt ) &
                    + latice * ( state%q(icol,k,ixcldliq) &
                                 + ptend%q(icol,k,ixcldliq) * dt ) ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 total_energy_column_finish(icol) &
                 = total_energy_column_finish(icol) &
                   + latice * ( state%q(icol,k,ixrain) &
                                + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
              total_energy_column_finish(icol) &
              = total_energy_column_finish(icol) &
                + latice * precl(icol) * dt * 1000.0_r8
           enddo ! k = top_lev, pver
        enddo ! icol = 1, ncol

        ! Calculate the total relative error in each grid column.
        do icol = 1, ncol

           tot_water_rel_err(icol) &
           = abs( ( grand_total_water_column_finish(icol) &
                    - grand_total_water_column_start(icol) ) ) &
             / min( grand_total_water_column_finish(icol), &
                    grand_total_water_column_start(icol) )

           tot_energy_rel_err(icol) &
           = abs( ( total_energy_column_finish(icol) &
                    - total_energy_column_start(icol) ) ) &
             / min( total_energy_column_finish(icol), &
                    total_energy_column_start(icol) )

        enddo ! icol = 1, ncol

        ! Print an error message if any total water relative error is found to
        ! be greater than the threshold.
        if ( any( tot_water_rel_err >= err_thresh ) ) then
           print *, "Water conservation error reported in hole filling"
           do icol = 1, ncol
              if ( tot_water_rel_err(icol) >= err_thresh ) then
                 print *, "Column = ", icol, &
                          "Relative error = ", tot_water_rel_err(icol), &
                          "Column-integrated grand total water at start = ", &
                          grand_total_water_column_start(icol), &
                          "Column-integrated grand total water at finish = ", &
                          grand_total_water_column_finish(icol)
              endif ! tot_water_rel_err(icol) >= err_thresh
           enddo ! icol = 1, ncol
        endif ! any( tot_water_rel_err >= err_thresh )

        ! Print an error message if any total energy relative error is found to
        ! be greater than the threshold.
        if ( any( tot_energy_rel_err >= err_thresh ) ) then
           print *, "Energy conservation error reported in hole filling"
           do icol = 1, ncol
              if ( tot_energy_rel_err(icol) >= err_thresh ) then
                 print *, "Column = ", icol, &
                          "Relative error = ", tot_energy_rel_err(icol), &
                          "Column-integrated total energy at start = ", &
                          total_energy_column_start(icol), &
                          "Column-integrated total energy at finish = ", &
                          total_energy_column_finish(icol)
              endif ! tot_energy_rel_err(icol) >= err_thresh
           enddo ! icol = 1, ncol
        endif ! any( tot_energy_rel_err >= err_thresh )

     endif ! l_check_conservation


     return

   end subroutine subcol_SILHS_fill_holes_conserv

   !============================================================================
   subroutine fill_holes_sedimentation( dt, ncol, hm_start, pdel, &
                                        fallspeed_m_per_s, zi, qmin_hm, &
                                        hm_tend, prec )

     ! Description:
     ! After hydrometeor tendencies from microphysics processes were adjusted
     ! so that holes don't form in a hydrometeor field from microphysics
     ! processes, the sedimentation tendency was added back on to produce an
     ! updated total microphysics tendency.  The first-order "up-gravity"
     ! sedimentation method that was originally used is positive definite.
     ! However, since the microphysics process tendencies were altered so that
     ! holes don't form, it is possible that adding the old sedimentation
     ! tendencies back onto the new microphysics process tendencies could
     ! produce new total microphysics tendencies that cause holes to form.
     !
     ! In this subroutine, holes in a hydrometeor field are checked for after
     ! the updated microphysics tendency is applied.  If any are found, they are
     ! filled from positive hydrometeor mass found at grid levels below where
     ! the hole is found.  The levels that are used to fill are within range
     ! based on fallspeed of the hydrometeor.  If the level that contains the
     ! hole is within proximity to the surface, then the water that sedimented
     ! to the surface can be used to fill the hole, as well.

     !----------------------------------------------------------------------

     use ppgrid, only: &
         pcols

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     real(r8), intent(in) :: dt                   ! Time step duration

     integer, intent(in) :: ncol                  ! Number of grid columns

     real(r8), dimension(pcols,pver), intent(in) :: &
       hm_start, & ! Hydrometeor mixing ratio at start of microphysics  [kg/kg]
       pdel        ! Pressure difference between grid levels            [Pa]

     real(r8), dimension(pcols,pver), intent(in) :: &
       fallspeed_m_per_s    ! Hydrometeor mixing ratio fall speed     [m/s]

     real(r8), dimension(pcols,pverp), intent(in) :: &
       zi    ! Height of momentum (interface) grid levels    [m]

     real(r8), intent(in) :: &
       qmin_hm    ! Minimum threshold value of hydrometeor mixing ratio  [kg/kg]

     ! Input/Output Variables
     real(r8), dimension(pcols,pver), intent(inout) :: &
       hm_tend    ! Hydrometeor mixing ratio tendency  [kg/kg/s]

     real(r8), dimension(pcols), intent(inout) :: &
       prec    ! Precipitation rate (surface)        [m/s]

     ! Local Variables
     real(r8), dimension(pver) :: &
       hm_update, & ! Hydrometeor mixing ratio; start of sed. hole fill  [kg/kg]
       hm_curr      ! Current value of hydrometeor mixing ratio          [kg/kg]

     real(r8) :: &
       total_hole,          & ! Total mass of hole in hydrometeor       [kg/m^2]
       total_fill_mass,     & ! Total mass available to fill hole       [kg/m^2]
       hole_fillmass_ratio, & ! Ratio: total_hole / total_fill_mass     [-]
       fallspeed_Pa_per_s,  & ! Hydrometeor mixing ratio fall speed     [Pa/s]
       total_fall_Pa,       & ! Pressure "distance" hydrometeor fell    [Pa]
       sum_pdel               ! Sum of pdel over levels                 [Pa]

     logical, dimension(pver) :: &
       l_pos_hm  ! Flag for a hydrometeor having a positive (>= qmin_hm) value

     ! Flag for whether surface precipitation mass needs to be included in
     ! the total_fill_mass for hole filling.
     logical :: l_reached_surface

     integer :: icol  ! Grid column index

     integer :: k, idx  ! Vertical grid level indices

     ! Index of the lowest vertical grid level that needs to be included in the
     ! total_fill_mass for hole filling.
     integer :: lowest_level_idx


     ! Loop over all columns, performing any adjustments one column at a time.
     do icol = 1, ncol

        ! Calculate the updated value of the hydrometeor field based on the
        ! updated microphysics tendency.  Since the original sedimentation
        ! tendency has been added to the updated microphysics process tendency
        ! to produce the updated total microphysics tendency (hm_tend), the
        ! updated value of the hydrometeor field (hm_update) could be negative.
        hm_update = hm_start(icol,:) + hm_tend(icol,:) * dt
        hm_curr = hm_update

        ! Check for any holes in the vertical profile
        if ( any( hm_curr(top_lev:pver) < qmin_hm ) ) then

           ! At least one hole is found in this hydrometeor species in this
           ! grid column.  The holes must be filled conservatively.

           ! Check which levels have values of the hydrometeor that are at or
           ! above the minimum threshold value.
           do k = top_lev, pver
              if ( hm_curr(k) >= qmin_hm ) then
                l_pos_hm(k) = .true.
              else ! hm_curr < qmin_hm
                l_pos_hm(k) = .false.
              endif ! hm_curr >= qmin_hm
           enddo ! k = top_lev, pver

           do k = pver, top_lev, -1

              if ( .not. l_pos_hm(k) ) then

                 ! A hole is found in the hydrometeor at this grid level.
                 if ( k == pver ) then

                    ! A hole is found at the lowermost level.
                    ! The only place the hydrometeor could have sedimented
                    ! to is the surface, so fill from only the surface.

                    ! Calculate the total hydrometeor mass of the hole that
                    ! needs to be filled.
                    ! The value of the hydrometeor mixing ratio is negative,
                    ! but total_hole is positive.
                    total_hole &
                    = ( qmin_hm - hm_curr(k) ) * pdel(icol,k) / gravit

                    ! Calculate the available amount of hydrometeor mass to
                    ! fill the hole.
                    total_fill_mass = prec(icol) * dt * 1000.0_r8

                    ! Calculate the ratio of total hole to total fill mass.
                    ! This should not exceed 1, but use thresholding to be safe.
                    hole_fillmass_ratio &
                    = min( total_hole / max( total_fill_mass, 1.0e-30 ), &
                           1.0_r8 )

                    ! Modify (reduce) the amount of surface precipitation in
                    ! order to fill the hole.  Since dt and 1000 are constants,
                    ! the only variable that needs to be modified
                    ! proportionately is prec.
                    prec(icol) = prec(icol) * ( 1.0_r8 - hole_fillmass_ratio )

                    ! Update the value of the hydrometeor at the level where the
                    ! hole was found.  Mathematically, as long as the available
                    ! mass was able to fill the entire hole, the new value of
                    ! the hydrometeor mixing ratio (hm_curr) should be 0.
                    hm_curr(k) &
                    = hm_curr(k) &
                      + hole_fillmass_ratio * total_fill_mass &
                        * gravit / pdel(icol,k)

                 else ! top_lev <= k < pver

                    ! Calculate the hydrometeor fallspeed in Pa/s.
                    ! In MG2, the equation for this is given by:
                    !
                    ! fallspeed([Pa/s]) = g * rho * fallspeed([m/s]).
                    !
                    ! The value of rho is typically calculated from the
                    ! hydrostatic approximation:
                    !
                    ! rho = - ( 1 / g ) * dp/dz.
                    !
                    ! The equation for fallspeed in Pa/s becomes:
                    !
                    ! fallspeed([Pa/s]) = - dp/dz * fallspeed([m/s]).
                    fallspeed_Pa_per_s &
                    = fallspeed_m_per_s(icol,k) &
                      * pdel(icol,k) / ( zi(icol,k) - zi(icol,k+1) )

                    ! Calculate the fall "distance" in Pa.
                    total_fall_Pa = fallspeed_Pa_per_s * dt

                    ! Find the index of the vertical level that the hydrometeor
                    ! sedimented to in one timestep.  It must sediment at least
                    ! one level.
                    sum_pdel = 0.0_r8
                    idx = k + 1
                    do
                       ! Update the total pressure difference between the
                       ! level of origin and the current level.
                       sum_pdel = sum_pdel + pdel(icol,idx)
                       if ( sum_pdel >= total_fall_Pa ) then
                          ! The total pressure difference between the level of
                          ! origin and the current level exceeds the total
                          ! hydrometeor fall "distance" (in Pa).
                          lowest_level_idx = idx
                          l_reached_surface = .false.
                          exit
                       else ! sum_pdel < total_fall_Pa
                          ! The total hydrometeor fall "distance" (in Pa)
                          ! exceeds the total pressure difference between the
                          ! level of origin and the current level.
                          if ( idx == pver ) then
                             ! The lowest level of the model has been reached.
                             ! The hydrometeor sedimented to the surface.
                             lowest_level_idx = pver
                             l_reached_surface = .true.
                             exit
                          else ! idx < pver
                             ! Increment idx and keep going.
                             idx = idx + 1
                          endif ! idx == pver
                       endif ! sum_pdel >= total_fall_Pa
                    enddo

                    ! Calculate the total hydrometeor mass of the hole that
                    ! needs to be filled.
                    ! The value of the hydrometeor mixing ratio is negative,
                    ! but total_hole is positive.
                    total_hole &
                    = ( qmin_hm - hm_curr(k) ) * pdel(icol,k) / gravit

                    ! Calculate the available amount of hydrometeor mass to
                    ! fill the hole.
                    total_fill_mass = 0.0_r8
                    if ( l_reached_surface ) then
                       ! The hydrometeor sedimented to the surface, so
                       ! automatically loop down to pver and include the
                       ! surface mass.
                       do idx = k+1, pver, 1
                          if ( l_pos_hm(idx) ) then
                             total_fill_mass &
                             = total_fill_mass &
                               + ( hm_curr(idx) - qmin_hm ) &
                                 * pdel(icol,idx) / gravit
                          endif ! l_pos_hm(idx)
                       enddo ! idx = k+1, pver, 1
                       ! Contribution to total fill mass from the surface.
                       total_fill_mass &
                       = total_fill_mass + prec(icol) * dt * 1000.0_r8
                    else ! .not. l_reached_surface
                       ! The hydrometeor sedimented to lowest_level_idx.
                       idx = k + 1
                       do
                          if ( l_pos_hm(idx) ) then
                             total_fill_mass &
                             = total_fill_mass &
                               + ( hm_curr(idx) - qmin_hm ) &
                                 * pdel(icol,idx) / gravit
                          endif ! l_pos_hm(idx)
                          if ( idx >= lowest_level_idx ) then
                             ! Check if enough mass has been gathered in
                             ! total_fill_mass to fill the hole.
                             if ( total_fill_mass >= total_hole ) then
                                ! There has been enough total_fill_mass
                                ! gathered to completely fill the hole.
                                lowest_level_idx = idx
                                exit
                             else ! total_fill_mass < total_hole
                                ! Even though lowest_level_idx has been reached,
                                ! more total_fill_mass needs to be added in
                                ! order to completely fill the hole, so keep
                                ! going.
                                if ( idx == pver ) then
                                   ! The lowest vertical level has already been
                                   ! reached, so go to the surface.
                                   lowest_level_idx = pver
                                   l_reached_surface = .true.
                                   ! Contribution to total fill mass from the
                                   ! surface.
                                   total_fill_mass &
                                   = total_fill_mass &
                                     + prec(icol) * dt * 1000.0_r8
                                   exit
                                else ! idx < pver
                                   ! Haven't reached pver yet, so increment
                                   ! and keep going.
                                   idx = idx + 1
                                endif ! idx == pver
                             endif ! total_fill_mass >= total_hole
                          else ! idx < lowest_level_idx
                             ! Haven't reached lowest_level_idx yet, so increment
                             ! and keep going.
                             idx = idx + 1
                          endif ! idx >= lowest_level_idx
                       enddo
                    endif ! l_reached_surface

                    ! Calculate the ratio of total hole to total fill mass.
                    ! This should not exceed 1, but use thresholding to be safe.
                    hole_fillmass_ratio &
                    = min( total_hole / max( total_fill_mass, 1.0e-30 ), &
                           1.0_r8 )

                    ! Modify (reduce) the amount of the hydrometeor at levels
                    ! that were used to fill the hole.
                    do idx = k+1, lowest_level_idx
                       if ( l_pos_hm(idx) ) then
                          ! Since pver at a grid level does not change and
                          ! gravit is constant, the only variable that needs to
                          ! be modified proportionately is hm_curr.
                          hm_curr(idx) &
                          = qmin_hm &
                            + ( hm_curr(idx) - qmin_hm ) &
                              * ( 1.0_r8 - hole_fillmass_ratio )
                       endif ! l_pos_hm(idx)
                    enddo ! idx = k+1, lowest_level_idx

                    if ( l_reached_surface ) then
                       ! Modify (reduce) the amount of surface precipitation in
                       ! order to fill the hole.
                       ! Since dt and 1000 are constants, the only variable that
                       ! needs to be modified proportionately is prec.
                       prec(icol) &
                       = prec(icol) * ( 1.0_r8 - hole_fillmass_ratio )
                    endif ! l_reached_surface

                    ! Update the value of the hydrometeor at the level where the
                    ! hole was found.  Mathematically, as long as the available
                    ! mass was able to fill the entire hole, the new value of
                    ! the hydrometeor mixing ratio (hm_curr) should be 0.
                    hm_curr(k) &
                    = hm_curr(k) &
                      + hole_fillmass_ratio * total_fill_mass &
                        * gravit / pdel(icol,k)

                 endif ! k == pver

              endif ! .not. l_pos_hm(k)

           enddo ! k = pver, top_lev, -1

        endif ! any( hm_curr(top_lev:pver) < qmin_hm )

        ! Update the value of total microphysics tendency after hole filling.
        hm_tend(icol,:) = hm_tend(icol,:) + ( hm_curr - hm_update ) / dt

     enddo ! icol = 1, ncol


     return

   end subroutine fill_holes_sedimentation

   !============================================================================
#endif
   
end module subcol_SILHS 
