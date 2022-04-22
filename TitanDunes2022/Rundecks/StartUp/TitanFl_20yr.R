P4SM40_likeTitan_jb_05302018.R GISS Model E  2004 modelE                     rar     07/15/2009
!! P4qsM20.R GISS Model E  2004 modelE 65m q-flux ocn      rar     07/15/2009

!! P4qsM20: P4SM40 with 65m q-flux ocean ("!!"-lines are for P4qsM20.R)
P4SM40_likeTitan_jb_05302018: modelE as frozen in July 2009 without gravity wave drag
modelE 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology  (see OSST/SICE)
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)
radiation: planet (SOCRATES)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define USE_ENT
#define NEW_IO
#define PLANET_PARAMS likeTitan    ! define to avoid tropwmo warnings about tropopause
#define USE_PLANET_RAD
#define GISS_RAD_OFF


End Preprocessor Options

Object modules: (in order of decreasing priority)
     ! resolution-specific source codes

Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL40p STRAT_DUM          ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M FFT72


IO_DRV                              ! new i/o

    ! GISS dynamics
! ATMDYN MOMEN2ND                     ! atmospheric dynamics
ATMDYN_TITAN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV TQUS_DRV                    ! advection of Q/tracers

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                            ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
Zenith	   			    ! shared zenith calculations
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OCEAN OCNML                         ! ocean modules

planet_rad planet_alb lw_control sw_control ! planet radiation source files

Components:
shared MPI_Support solvers giss_LSM
dd2d
socrates
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
    ! resolution dependent files
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy           ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
! AIC=NCARIC.72x46.D7712010_ext.nc     ! initial conditions (atm.)      needs GIC, ISTART=2
! GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground)
! AIC=AIC.RES_M20A.D771201.nc          ! initial conditions (atm.)      needs GIC, ISTART=2
! AIC=AIC.RES_F40.D771201.nc           ! initial conditions (atm)
AIC=/home/maxqc/titan_input/1JAN1970.rsfTitan_20yr_1949.nc        ! New atmo initial cond for Titan 8/18/17
! Initial atmospheric temperature, specific humidity, wind, surface pressure.
! If rundeck parameter initial_psurf_from_topo=1, initial surface pressure is
! reset to be consistent with orography.
! AIC=planet/Mars/AIC.coldmars.nc
! GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground)
!
SOILIC=planet/desert_world/soilic_drysoil.nc
!

! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
! For q-flux ocean, replace lines above by the next 2 lines & set KOCEAN=1, ISTART=8
!! AIC=1JAN1961.rsfE4M20.MXL65m   ! = end of preliminary run with KOCEAN=0,Kvflxo=1
!! OHT=OTSPEC.E4M20.MXL65m.1956-1960 ! ocean horizontal heat transport
! OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (needed for post processing)
! TOPO=Z72X46N.cor4_nocasp.nc       ! topography
TOPO=/home/maxqc/titan_input/titan_topo_flat.nc   ! flat topography
! SOIL=S4X50093.ext.nc              ! soil bdy.conds
SOIL=planet/desert_world/soil_allsand.nc
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
! VEG=V72X46.1.cor2_no_crops.ext.nc ! veg. fractions
VEG=/home/maxqc/titan_input/veg_allbare_alb22.nc  ! all bare soil with albedo = 0.20
! CROPS=CROPS2007_72X46N.cor4_nocasp.nc       ! crops history
soil_textures=soil_textures_top30cm
SOILCARB_global=soilcarb_top30cm_4x5.nc
CDN=CD4X500S.ext.nc               ! surf.drag coefficient
REG=REG4X5                        ! special regions-diag
! RVR=RD_modelE_M.nc                ! river direction file
! NAMERVR=RD_modelE_M.names.txt     ! named river outlets
! TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
TOP_INDEX=/home/maxqc/titan_input/titan_srfcsgma.nc   ! Titan topography std. deviations
! GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution
    ! resolution independent files
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table

!ISCCP=ISCCP.tautables
!GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
GHG=/home/maxqc/titan_input/GHG.Titan.txt     ! fixed 50000 ppm CH4

! *** Aerosols and O3 are deactivated by default in planet rundecks ***
!RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr ! no volcanic aerosols
!RADN8=cloud.epsilon4.72x46 ! map of CLDEPS only applicable to modern Earth
!dH2O=dH2O_by_CH4_monthly ! no stratospheric generation of H2O by CH4
!DUSTaer=dust_mass_CakmurMillerJGR06_72x46x20x7x12.nc
!BC_dep=BC.Dry+Wet.depositions.ann.nc
! updated aerosols need MADAER=3
!TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h.nc
!TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h.nc
!TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h.nc
!TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h.nc
!TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h.nc
!TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h.nc
!O3file=o3_2005_shindelltrop_72x46x49_1850-1997.nc
MSU_wts=MSU.RSS.weights.data

Label and Namelist:
P4SM40_likeTitan_jb_05302018 (ModelE1 4x5, 20 lyrs, 1850 atm/ocn)

&&PARAMETERS

! calculate initial surface pressure for hydrostatic consistency with topography
initial_psurf_from_topo=1

! minimum and maximum allowed column mass (kg/m2) for error checking
mincolmass=0.0
maxcolmass=200000.0

planetName = 'likeTitan'
eccentricity = 0.056       ! Saturn J2000 values from Allison & Ferrier 2004 unpublished
obliquity  = 26.4d0        ! degrees, Titan  J2000 values from Allison & Ferrier 2004 unpublished
! Earths sidereal rotation period is 86400 * (365/366) = 86163.934426229375 s
siderealRotationPeriod=1374315.0d0       ! seconds (x 15.95 days)
siderealOrbitalPeriod=929523600.d0                   ! seconds (x 29.475 years)
quantizeYearLength='False'
planet_s0=14.82            ! Solar irradiance W/m^2.


! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)
ocn_cycl=1      ! =0 if ocean varies from year to year

variable_lk=1   ! variable lakes

wsn_max=2.   ! restrict snow depth to 2 m-h2o (if 0. snow depth is NOT restricted)
! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! Input files for planet radiation (SOCRATES)
solar_spec_dir='/home/maxqc/ModelE_Support/socrates/stellar_spectra'
spectral_dir='/home/maxqc/ModelE_Support/socrates/spectral_files'
solar_spec='sun'
!spectral_file_lw='sp_lw_ga7/sp_lw_ga7_dsa'
!spectral_file_sw='sp_sw_ga7/sp_sw_ga7_dsa'
spectral_file_lw='sp_lw_dsa_titan/sp_lw_10_dsa_titan'           ! For Titan
spectral_file_sw='sp_sw_dsa_titan/sp_sw_14_dsa_titan_sun'
l_cloud=0                       ! These Titan spectral files do not contain water optical properties
!aer_opt_prop_lw='sp_lw_ga7/aer_lw_ga7.nc'
!aer_opt_prop_sw='sp_sw_ga7/aer_sw_ga7.nc'

! tuning param.: this setting works for 1850; use U00wtrX=1.28 for 1979

! Tuning parameters as of 2016/10/24 (D.S. Amundsen)
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=0.50   ! below 850mb and MC regions; then tune this to get rad.balance
radiusl_multiplier=0.97

! Cloud inhomogeneity correction
KCLDEP=1    ! use a constant value for CLDEPS
EPSCON=0.12 ! use CLDEPS=0.12

! Stellar constant at 1 AU (default = 1360.67 W/m2)
!planet_s0=1360.67
! S0X has been removed and will have no effect on results

minGroundTemperature=-300.0d0

CO2X=0.
CH4X=1.0
O2X=0.
NO2X=0.
N2OX=0.
CFC11X=0.
CFC12X=0.
N2CX=0.
XGHGX=0.
YGHGX=0.
SO2X=0.
O3X=0. ! turn of O3, to turn on delete this line and add O3file

H2OstratX=1.

H2ObyCH4=0.     ! deactivates strat.H2O generated by CH4
l_uniform_ghg=1 ! deactivates prescribed GHG distributions
KSIALB=0        ! 6-band albedo (Hansen) (=-1 no land ice fixup, 1 Lacis scheme)

!madaer=3        ! updated aerosols
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
crops_yr=-1
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! dont include 1st indirect effect
cc_cdncx=0.        ! dont include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
!dalbsnX=.024      ! no BC deposition
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD

! parameters that control the Shapiro filter
! DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
! DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.      ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
! DT=450.       ! 8/13/14: added 4 lines after NIsurf=2
NIsurf=2        ! increase as layer 1 gets thinner
nrad=1
! DT=225.         ! Changed per advice of Max (7/17/14) to fix aadvtx crashes
! DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
! DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
!
! With smaller radius (Titan 2575 km), need a smaller timestep
! Now, DTsrc / DT  has to be an EVEN integer, so 1800 / 150 is allowed
DT=225.0          ! Changed per advice of Max (10/21/14)to fix ADVECM crashes
DT_XUfilter=225.0 ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225.0 ! Shapiro filter on V in E-W direction; usually same as DT


! parameters that affect at most diagn. output:
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
iCOPY=2          ! saving acc + rsf
isccp_diags=0    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
&&END_PARAMETERS

 &INPUTZ
 YEARI=1950,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1970,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=8,IRANDI=0, YEARE=1950,MONTHE=1,DATEE=1,HOURE=1,
!! suggested settings for P4qsF40:
!! YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!! YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
!! ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
/
