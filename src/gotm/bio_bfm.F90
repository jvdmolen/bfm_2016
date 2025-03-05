!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_bfm --- BFM bio model \label{sec:bio_bfm}
!
! !INTERFACE:
   module bio_bfm
!
! !DESCRIPTION:
!
!
! !USES:
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_bfm, pointers_gotm_bfm,            &
          var_info_bfm, set_env_bio_bfm, do_bio_bfm,  &
          allocate_memory_bfm,reset_diagonal,         &
          test_on_all_negative_states,                &
          test_mass_conservation,                     &
          test_on_negative_states, end_bio_bfm,       &
          do_bfm_river_loads,assign_adv_rates,        &
          CalcVertFluxAtLev,calc_sigma_depth,         &
          DeriveFromGotm,return_bfm,test_structure, &
          test_model_states,test_model_rates,print_structure
!
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the template bio module
!
! !INTERFACE:
   subroutine init_bio_bfm(nlev,out_unit)
   use bio_var,only:numc,numbc,numc_diag,numc_flux,numbc_flux,numcc, &
                    numbc_diag
#ifdef INCLUDE_DIAGNOS_PRF
   use bio_var,only:numbc_prf,nprf
#endif
!
! !DESCRIPTION:
!  Here, the main communication of array dimensions between GOTM
!  and BFM is done.
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_D2_BOX_FLUX, NO_D3_BOX_FLUX,&
                  NO_STATES
#ifdef INCLUDE_DIAGNOS_PRF
   use mem,only:  NO_BOXES_PRF,NO_D3_BOX_DIAGNOSS_PRF
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: nlev
   integer,          intent(in)   :: out_unit
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_bfm'


   ! BFM  --> GOTM
   numc  = NO_D3_BOX_STATES
   numbc = NO_D2_BOX_STATES
   numc_diag  = NO_D3_BOX_DIAGNOSS
   numbc_diag = NO_D2_BOX_DIAGNOSS
   numc_flux  = NO_D3_BOX_FLUX
   numbc_flux = NO_D2_BOX_FLUX
   ! numcc is the number of transported variables
   numcc = numc

   ! GOTM --> BFM
   NO_BOXES_X  = 1
   NO_BOXES_Y  = 1
   NO_BOXES_Z  = nlev
   NO_BOXES    = NO_BOXES_X * NO_BOXES_Y * NO_BOXES_Z
   NO_BOXES_XY = NO_BOXES_X * NO_BOXES_Y
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY
#ifdef INCLUDE_DIAGNOS_PRF
   numbc_prf= NO_D3_BOX_DIAGNOSS_PRF
   nprf= NO_BOXES_PRF
#endif
   !LOGUNIT = out_unit

   LEVEL3 'pelagic variables =',numc
   LEVEL3 'pelagic transported variables =',numcc
   LEVEL3 'benthic variables =',numbc
   LEVEL3 'pelagic variables prepared for output',numc_diag
   LEVEL3 'benthic variables prepared for output',numbc_diag
#ifdef INCLUDE_DIAGNOS_PRF
   LEVEL3 'benthic profile variables prepared for output',numbc_prf
#endif
   LEVEL3 'NO_BOXES_X=',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y=',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z=',NO_BOXES_Z
   LEVEL3 'NO_BOXES=',NO_BOXES
   LEVEL3 'NO_BOXES_XY=',NO_BOXES_XY
#ifdef INCLUDE_DIAGNOS_PRF
   LEVEL3 'NO_BOXES_PRF=',NO_BOXES_PRF
#endif
   LEVEL3 'NO_STATES=',NO_STATES
   LEVEL3 'Step 1 of GOTM <-> BFM initialisation done ...'

!  sfl=_ZERO_
!  sfl_read=_ZERO_
   return

   end subroutine init_bio_bfm
!EOC
!-----------------------------------------------------------------------

! init_var_bfm is found ../share
!

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize BFM and GETM shared memory
!
! !INTERFACE:
   subroutine pointers_gotm_bfm()
!
! !DESCRIPTION:
! Allocate pointers to GOTM memory
! Used by init_var_bfm
!
! !USES:
   use bio_var,only:bio_setup,cc,pp,dd,pelvar_type,diag, &
                ccb,ppb,ddb,benvar_type,diagb,numc_diag,numbc_diag
   use mem, only: D3STATE,D3SOURCE,D3SINK,D3STATETYPE, &
                  D3DIAGNOS,D2STATE,D2SOURCE,D2SINK,   &
                  D2STATETYPE,NO_BOXES,NO_BOXES_XY,    &
                  D2DIAGNOS,NO_D2_BOX_STATES,          &
                  NO_D2_BOX_DIAGNOSS
#ifdef INCLUDE_DIAGNOS_PRF
      use mem,only:D3DIAGNOS_PRF,NO_BOXES_PRF,NO_D3_BOX_DIAGNOSS_PRF
      use bio_var,only:numbc_prf,diagb_prf,nprf
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!

   !---------------------------------------------
   ! Pelagic pointers
   !---------------------------------------------
!JM   D3STATE  => cc(:,1:NO_BOXES)
!JM   D3SOURCE => pp(:,:,1:NO_BOXES)
!JM   D3SINK   => dd(:,:,1:NO_BOXES)
   D3STATE  => cc(1:NO_BOXES,:)
   D3SOURCE => pp(1:NO_BOXES,:,:)
   D3SINK   => dd(1:NO_BOXES,:,:)
   D3STATETYPE => pelvar_type
   if (numc_diag > 0) D3DIAGNOS => diag(1:NO_BOXES,:)

   !---------------------------------------------
   ! Benthic pointers
   !---------------------------------------------
   if (bio_setup >=2 ) then
!JM      D2STATE  => ccb(:,1:NO_BOXES_XY)
!JM      D2SOURCE => ppb(:,:,1:NO_BOXES_XY)
!JM      D2SINK   => ddb(:,:,1:NO_BOXES_XY)
      D2STATE  => ccb(1:NO_BOXES_XY,:)
      D2SOURCE => ppb(1:NO_BOXES_XY,:,:)
      D2SINK   => ddb(1:NO_BOXES_XY,:,:)
      D2STATETYPE => benvar_type
      if (numbc_diag>0) D2DIAGNOS => diagb(1:NO_BOXES_XY,:)
#ifdef INCLUDE_DIAGNOS_PRF
      if (numbc_prf>0) D3DIAGNOS_PRF => diagb_prf(1:NO_BOXES_PRF,:)
#endif
   else
      ! allocate memory anyhow to avoid problems with BFM allocation
!JM      allocate(D2STATE(1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
!JM      allocate(D2SOURCE(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
!JM      allocate(D2SINK(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
!JM      allocate(D2STATETYPE(1:NO_D2_BOX_STATES ))
      allocate(D2STATE(1:NO_BOXES_XY,1:NO_D2_BOX_STATES))
      allocate(D2SOURCE(1:NO_BOXES_XY,1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES))
      allocate(D2SINK(1:NO_BOXES_XY,1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES))
      allocate(D2STATETYPE(1:NO_D2_BOX_STATES ))
      if (numbc_diag>0)  &
         allocate(D2DIAGNOS(1:NO_BOXES_XY,1:NO_D2_BOX_DIAGNOSS))
      if (numbc_prf>0)  &
         allocate(D2DIAGNOS(1:NO_BOXES_PRF,1:NO_D3_BOX_DIAGNOSS_PRF))
   end if

   end subroutine pointers_gotm_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_bfm()
!
! !DESCRIPTION:
!  This subroutine provides information on the variables. To be used
!  when storing data in NetCDF files.
!
! !USES:
   use mem
   use bfm_output,only:setup_bio_output
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call setup_bio_output
   call set_var_info_bfm
   return
   end subroutine var_info_bfm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine set_env_bio_bfm(nlev,dt,bathy_dep,u_taub,h,t,s,rho,nuh,dl, &
          uwind_gotm,vwind_gotm,ulevel,vlevel, I_0,julianday,dry_z_getm,abioshade)
!JMC ifdef SPM=true, then abioshade is filled from kd_spm in spm module
!
! !DESCRIPTION
!
! !USES
! BFM modules
 use bio_var, only: nlev_local=>nlev,I_0_local => I_0,dt_local=>dt,h_l,nuh_l, &
                    bio_julianday=>julianday,bathy_dep_local=>bathy_dep
 use mem,       ONLY: NO_BOXES, ESS, ERHO,SUNQ, &
                     EUWIND,EVWIND,ETW, ESW, Wind,R9x,EUCURR_LEVEL,EVCURR_LEVEL,    &
                     dry_z,ETAUB,OCDepth,Depth, EIR,EIRr, ABIO_eps
 use mem_Param,  ONLY: use_function_for_eps0,p_eps0,p_poro


IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(IN)                  :: nlev
   REALTYPE, intent(in)                :: dt
   REALTYPE, intent(in)                :: bathy_dep
   REALTYPE, intent(in)                :: u_taub
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
   REALTYPE, intent(in)                :: rho(0:nlev)
   REALTYPE, intent(in)                :: s(0:nlev)
   REALTYPE, intent(in)                :: nuh(0:nlev)
   REALTYPE, intent(in)                :: I_0
   REALTYPE, intent(in)                :: dl
   integer, intent(in)                 :: julianday
   REALTYPE, intent(in)                :: uwind_gotm,vwind_gotm,ulevel,vlevel
   REALTYPE,intent(in),optional        :: dry_z_getm
   REALTYPE,intent(in),optional        :: abioshade(0:nlev)
!
! !OUTPUT PARAMETERS:
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer             :: n
   REALTYPE            :: psilt
!EOP
!-----------------------------------------------------------------------
!BOC

!   LEVEL2 'calculating environmental forcings for the BFM'
!LEVEL2 'set_env_bio_bfm, rho',rho
   !---------------------------------------------
   ! Assign depths of layers
   ! temperature and salinity
   !---------------------------------------------

    nlev_local=nlev
    dt_local=dt
    bathy_dep_local = bathy_dep
    I_0_local = I_0
    bio_julianday=julianday
    dry_z=_ONE_
    if (present(dry_z_getm)) dry_z(1)=dry_z_getm

    Depth(:) = h(1:nlev)
    OCDepth(NO_BOXES) =depth(NO_BOXES)
    EUWIND(1)=uwind_gotm
    EVWIND(1)=vwind_gotm
    Wind= sqrt( EUWIND(1) * EUWIND(1) + EVWIND(1) *EVWIND(1))
    ETAUB(1)=u_taub
!JM pass currents
    EUCURR_LEVEL(1)=ulevel
    EVCURR_LEVEL(1)=vlevel
    do n=NO_BOXES-1,1,-1
       OCDepth(n)=depth(n)+ OCDepth(n+1)
    enddo
   ETW(1:nlev) = t(1:nlev)
   ESW(1:nlev) = s(1:nlev)
   ERHO(1:nlev) = rho(1:nlev)
   psilt=(p_poro(1) - 0.38662 )/ 0.00415
   ESS(:)=R9x(:)
   nuh_l=nuh
   h_l=h
!LEVEL2 'set_env_bio_bfm, I_0',I_0
!stop
   !---------------------------------------------
   ! Compute extinction coefficient
   !---------------------------------------------

#ifdef SPM
     p_eps0=0.0
     ABIO_eps(1:nlev) = abioshade(1:nlev)
#else
     select case (use_function_for_eps0)
        case (1) ; p_eps0=1.17692308-0.0307692308*min(35.0,ESW(1)) !northsea
        case (2) ; p_eps0=1.17692308-0.0307692308*max(30.0,min(35.0,ESW(1))) !wz
     end select
#endif


!  if (bioshade_feedback) then
     ! 0= Special calculation  of vertical extinction only controlled by
     !   biological constituents for use in gotm
     ! 1= full calulation used to caluclate bilolgical vertical extinction
     ! 2= Special calculation  of vertical extinction only controlled by
     !   biological constituents and silt when dilt inis include in the BFM part.
!    call  CalcVerticalExtinction(0)
!    bioshade(1)=1.0
!    do i=nlev,2,-1
!      bioshade(i-1) = bioshade(i)*exp(-xEPS(i)*Depth(i))
!    end do
!    bioshade(1:nlev) =  bioshade(1:nlev)*exp(-xEPS(:)*Depth(:)*0.5)
!  endif

   !---------------------------------------------
   ! Notice that irradiance in the BFM is in
   ! uE/m2/s and is defined at the top of each
   ! layer (the derivation of the middle-layer
   ! EIR for production is done in the
   ! Phytoplankton routines)
   !---------------------------------------------
   EIRr=I_0
   SUNQ=dl
!  call  CalcVerticalExtinction(1)
!  EIR(nlev) = max(p_small,p_PAR*I_0/E2W)
!  do i=nlev,2,-1
!    EIR(i-1) = EIR(i)*exp(-xEPS(i)*Depth(i)*dry_z(1))
!  end do

   !---------------------------------------------
   ! bioshade is instead derived in the
   ! middle of the layer and it's non-dimensional
   !---------------------------------------------

   end subroutine set_env_bio_bfm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the BFM model
!
! !INTERFACE
   subroutine do_bio_bfm(first)
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!
! !USES
   use bio_var,only:bio_setup,surface_flux_method,sfl,bfl,sfl_N3n,sfl_N4n, &
                    surface_flux_method,test_mode
   use mem_param,only:AssignAirPelFluxesInBFMFlag,AssignPelBenFluxesInBFMFlag, &
                      CalcBenthicFlag
   use mem, only:     iiC,iiN,iiP,iiS,iiL, &
                  ppR2c, ppRZc, ppR6c, ppR6n, ppR6p, ppR6s, NO_BOXES_Z,   &
                  ppR1c, ppR1n, ppR1p,ppO3c,NO_BOXES,   &
                  ppO2o,ppN1p,ppN3n,ppN4n,ppN5s,ppN6r,  &
                  NO_D3_BOX_STATES, Depth, &              
                  ppPhytoPlankton,iiPhytoPlankton, &
                  PELBOTTOM, PELSURFACE, &
                  jK3G4n,jK23K13n,jK34K24n,flN3O4n,jK23G4n
use mem, only: ERHO
use bfm_output,only:var_ave
!JM #IFDEF INCLUDE_MACROPHYT
!JM #ENDIF
!JM #IFDEF INCLUDE_DAAN
!JM     use  Daan
!JM #ENDIF

   use constants,  only: SEC_PER_DAY
   use gotm_error_msg, only:get_d3_model_flag_from_getm

   IMPLICIT NONE
!
   logical,intent(in)          :: first


   logical                     :: d3_model_flag
   logical,save                :: start=.true.
   integer                     :: i,k
   REALTYPE                    :: topm3psec,Nloss,Hloss
!See "USE association" above
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from template by Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

!LEVEL1 'do_bio_bfm'

   !---------------------------------------------
   ! Compute BFM terms
   !---------------------------------------------
!JM #IFDEF INCLUDE_DAAN
!JM    !DJG routines from Daan.f90 for varying constants:
!JM    call change_fish_predation()
!JM #ENDIF
   call SiltDynamics
!LEVEL1 'ecologydynamics'
!LEVEL1 'ERHO',ERHO
!stop
   call EcologyDynamics
!LEVEL1 'benthicsilt'
!stop
   call BenthicSiltDist(0,0.0D+00)
!LEVEL1 'get_d3_model_flag'
!stop
   call get_d3_model_flag_from_getm(d3_model_flag)
   Nloss=_ZERO_;Hloss =_ZERO_;
   if ( .not.(d3_model_flag.or.start.or.test_mode)) then
      if (CalcBenthicFlag==3)  then
        Nloss=0.0;
        Nloss=jK3G4n(1)-jK23K13n(1)-jK34K24n(1)+sum(flN3O4n(:)*Depth(:));
!       Hloss=0.0;
        Hloss= -0.5 * jK3G4n(1)- 0.5*sum(flN3O4n(:)*Depth(:)) +(jK23K13n(1)+jK23G4n(1)) -0.5* jK23G4n(1)
!       Hloss= -0.5 * jK3G4n(1)- 0.5*sum(flN3O4n(:)*Depth(:)) +jK13K3n(1)
      elseif (CalcBenthicFlag==0)  then
        Nloss=sum(flN3O4n(:)*Depth(:))
        Hloss=sum(flN3O4n(:)*Depth(:))
      endif
   endif
   start=.false.
!LEVEL1 'surface fluxes',bio_setup
!stop
   !---------------------------------------------
   ! Surface fluxes
   !---------------------------------------------
   if ( bio_setup ==2 ) return
!  topm3psec=_ONE_/Depth(NO_BOXES_Z)/ SEC_PER_DAY
   topm3psec=_ONE_/SEC_PER_DAY
   sfl=_ZERO_;
   if ( .NOT. AssignAirPelFluxesInBFMFlag ) then
!JM     sfl(ppO2o) =   PELSURFACE(ppO2o,1) *topm3psec
!JM     if ( ppO3C > 0 ) sfl(ppO3c) =   PELSURFACE(ppO3c,1) *topm3psec
     sfl(ppO2o) =   PELSURFACE(1,ppO2o) *topm3psec
     if ( ppO3C > 0 ) sfl(ppO3c) =   PELSURFACE(1,ppO3c) *topm3psec
   endif
   select case (surface_flux_method)
        case (-1)! absolutely nothing
        case (0) ! constant
           sfl(ppN3n) =   0.12  *topm3psec
           sfl(ppN4n) =   0.09  *topm3psec
           !this is called here to test track when the 1d model is used in 1D-mode.
           !In this case get sfl(pptrN3n) the same value as slf(ppN3n)
           ! It work only if d3_model_flag ==flase and if trakcing is active.
!          call fill_sfl(d3_model_flag,ppN3n,numc,sfl)
!          call fill_sfl(d3_model_flag,ppN4n,numc,sfl)
           sfl(ppN1p) =   _ZERO_  !0.0
        case (1) ! from file via sfl_read
           ! fluxes are in mmol m-2 d-1
           sfl(ppN3n) =    sfl_N3n  *topm3psec
           sfl(ppN4n) =    sfl_N4n  *topm3psec
!          call fill_sfl(d3_model_flag,ppN3n,numc,sfl)
!          call fill_sfl(d3_model_flag,ppN4n,numc,sfl)
        case (3) ! sfl array filled externally - for 3D models
         ! option 3 works only in 1D-mode!!!!!!!
           sfl(ppN3n)= Nloss * topm3psec
!          call flux(iiPel,NO_BOXES_Z,ppN3n,ppN3n,Nloss/Depth(NO_BOXES_Z))
           PELSURFACE(ppN3n,1)=Nloss
!          if ( ppO3h> 0) sfl(ppO3h)= Hloss * topm3psec
!          call flux(iiPel,NO_BOXES_Z,ppO3h,ppO3h,Hloss/Depth(NO_BOXES_Z))
        case default
   end select

   !---------------------------------------------
   ! Bottom fluxes
   !---------------------------------------------
   topm3psec=1.0/Depth(1)/ SEC_PER_DAY
   bfl=_ZERO_;
   if ((bio_setup == 3 ) .and. ( .NOT.AssignPelBenFluxesInBFMFlag)) then

      bfl(ppRZc) = PELBOTTOM(1,ppRZc)*topm3psec
      bfl(ppR2c) = PELBOTTOM(1,ppR2c)*topm3psec
      bfl(ppR6c) = PELBOTTOM(1,ppR6c)*topm3psec
      bfl(ppR6n) = PELBOTTOM(1,ppR6n)*topm3psec
      bfl(ppR6p) = PELBOTTOM(1,ppR6p)*topm3psec
      bfl(ppR6s) = PELBOTTOM(1,ppR6s)*topm3psec

      bfl(ppR1c) =  PELBOTTOM(1,ppR1c)*topm3psec
      bfl(ppR1n) =  PELBOTTOM(1,ppR1n)*topm3psec
      bfl(ppR1p) =  PELBOTTOM(1,ppR1p)*topm3psec

      bfl(ppO2o) = PELBOTTOM(1,ppO2o)*topm3psec
      bfl(ppN1p) = PELBOTTOM(1,ppN1p)*topm3psec
      bfl(ppN3n) = PELBOTTOM(1,ppN3n)*topm3psec
      bfl(ppN4n) = PELBOTTOM(1,ppN4n)*topm3psec
      bfl(ppN5s) = PELBOTTOM(1,ppN5s)*topm3psec
      bfl(ppN6r) = PELBOTTOM(1,ppN6r)*topm3psec

      do i=1,iiPhytoPlankton
        k=ppPhytoPlankton(i,iiC)
        bfl(k) = PELBOTTOM(1,k)*topm3psec
        k=ppPhytoPlankton(i,iiN)
        bfl(k) = PELBOTTOM(1,k)*topm3psec
        k=ppPhytoPlankton(i,iiP)
        bfl(k) = PELBOTTOM(1,k)*topm3psec
        k=ppPhytoPlankton(i,iiL)
        bfl(k) = PELBOTTOM(1,k)*topm3psec
        k=ppPhytoPlankton(i,iiS)
        if ( k > 0 ) bfl(k) = PELBOTTOM(1,k)*topm3psec
      enddo
   endif
!LEVEL1 'end do_bio_bfm'
!LEVEL1 'var_ave',var_ave
!stop
   end subroutine do_bio_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Assign  adv_rates
!
! !INTERFACE
   subroutine CalcVertFluxAtLev(statenr,lev,nlev,Depth,dt,out)
!
! !DESCRIPTION
! !USES
   use bio_var,only: cc_before_transport,cc
   use mem,only: PELBOTTOM
   use constants,  only: SEC_PER_DAY
   IMPLICIT NONE
!
   integer,intent(IN)          :: statenr
   integer,intent(IN)          :: lev
   integer,intent(IN)          :: nlev
   REALTYPE,intent(IN)         :: Depth(1:nlev)
   REALTYPE,intent(IN)         :: dt
   REALTYPE,intent(OUT)        :: out

   out=0.0D+00
!  if ( lev > 0.and.allocated(cc_before_transport) ) then
!      ! rate calculate in time unit of physical model (secs)
!      out= sum((cc_before_transport(statenr,lev+1:nlev) &
!              -cc(statenr,lev+1:nlev))*Depth(lev+1:nlev))/dt
!  elseif ( lev==0 ) then
!      ! rate calculate transferred from time unit of eco model
!      ! to the one  of physical model (secs)
!      out=PELBOTTOM(statenr,1)/SEC_PER_DAY;
!  else
!     out=_ZERO_
!  endif
!  return
   end subroutine CalcVertFluxAtLev
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Assign  adv_rates
!
! !INTERFACE
   subroutine assign_adv_rates(dt)
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!
! !USES
   use bio_var,only:bio_setup,ws,llws,c1dimz,rel_max_sedi_rate,pelvar_type,SILTTRANSPORT
   use mem, only: NO_BOXES_Z,   &
                  NO_D3_BOX_STATES, OCDepth,Depth,              &
                  D3DIAGNOS,iiPELSINKREF
   use bfm_output,only:var_names,var_ave
   use constants,  only: SEC_PER_DAY
   use mem_globalfun,only: insw_vector
use global_mem,only:LOGUNIT

   IMPLICIT NONE
!
   REALTYPE,intent(in)          :: dt

   logical                     :: ll_larger
   logical,save                :: first=.true.
   integer                     :: i,j,ldep,k
   REALTYPE                    :: corr(1:NO_BOXES_Z)
   REALTYPE                    ::r
   character(len=22)           ::onem,pnem

!write(LOGUNIT,*)'start assign_adv_rates, var_ave',var_ave
   !---------------------------------------------
   ! Transfer sinking velocities (m/d -> m/s)
   !---------------------------------------------
   if ( bio_setup ==2 ) return
   ldep=1
   do j=1,NO_D3_BOX_STATES
     i= iiPELSINKREF(j)
     if ( i > 0 ) then
       ll_larger=(maxval(abs(D3DIAGNOS(1:NO_BOXES_Z,i)))>0.001)
       c1dimz(NO_BOXES_Z)=0.0;
       if ( ll_larger) then
        c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*D3DIAGNOS(1:NO_BOXES_Z-1,i)&
                            +Depth(1:NO_BOXES_Z-1)*D3DIAGNOS(2:NO_BOXES_Z,i))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))
        if (pelvar_type(j)>=SILTTRANSPORT) then
           corr=1.0D+00
        else
           corr=min(OCDepth(ldep)*rel_max_sedi_rate, &
             abs(c1dimz(1:NO_BOXES_Z)))/ (1.0D-80+abs(c1dimz(1:NO_BOXES_Z)))
           corr=corr*OCDepth(ldep)/(5.0+OCdepth(ldep))
        endif
!JM           ws(j,1:NO_BOXES_Z) = -c1dimz(1:NO_BOXES_Z)/SEC_PER_DAY *corr
           ws(1:NO_BOXES_Z,j) = -c1dimz(1:NO_BOXES_Z)/SEC_PER_DAY *corr
       else
!JM          ws(j,1:NO_BOXES_Z) = 0.0;
          ws(1:NO_BOXES_Z,j) = 0.0;
       endif
       llws(j)=ll_larger
!JM       ws(j,0)= ws(j,1)
       ws(0,j)= ws(1,j)
     elseif (i< -NO_D3_BOX_STATES) then
       llws(j)=.false.
     elseif (-i==j) then
       ! only when iiPelSinkRef is equal to the negative value of its
       ! is possible to define straight a ws .
       llws(j)=.false.
!JM       ws(j,0)= ws(j,1)
       ws(0,j)= ws(1,j)
!      if ( OCDEPTH(1)>1.5)llws(j)=.true.
       llws(j)=.true.
     elseif (i<0) then
       if ( -i>=j) then
         STDERR "i=",i
         stop 'error:assign_adv_rates'
       endif
       i=-i
!JM       ws(j,0:NO_BOXES_Z) =ws(i,0:NO_BOXES_Z)
       ws(0:NO_BOXES_Z,j) =ws(0:NO_BOXES_Z,i)
       llws(j)=llws(i)
     else
       llws(j)=.false.
     endif
   enddo
!write(LOGUNIT,*)'mid assign_adv_rates, var_ave',var_ave
   if (first) then
     ldep=0;k=0
     do j=1,NO_D3_BOX_STATES
       if (k==0) then
           k=1 ; LEVEL3 "Vertical advective transport is present for"
       endif
        i= iiPELSINKREF(j)
        onem=var_names(j); ldep=len_trim(onem)
        if ( i > 0 ) then
          LEVEL4 onem(1:ldep)
        elseif (i <-NO_D3_BOX_STATES) then
          LEVEL4 onem(1:ldep)," excluded outside PelGlobal.F90"
        elseif (-i==j) then
          LEVEL4 onem(1:ldep)," defined outside PelGlobal.F90"
        elseif (i<0) then
          pnem=var_names(-i); i=len_trim(pnem)
          LEVEL4 onem(1:ldep)," coupled to ",pnem(1:i)
        endif
        if (pelvar_type(j)>=SILTTRANSPORT) &
          LEVEL4 onem(1:ldep), " Diffusive transport coupled to advective transport"  
     enddo
     first=.false.
     if ( i.ne.0) then  
        LEVEL3 "End of Vertical transport";LEVEL3 " "
     endif
   endif

!write(LOGUNIT,*)'end assign_adv_rates, var_ave',var_ave
!stop
   return

   end subroutine assign_adv_rates
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Reset diagonal id a 3d array
!
! !INTERFACE:
   subroutine reset_diagonal(n,pp)
!
! !DESCRIPTION:
!    Reset of the diagonal
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                    :: n
   REALTYPE,dimension(:,:,:),intent(inout) :: pp
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES:
   integer                   :: i
!EOP
!-----------------------------------------------------------------------
!BOC
     do i=1,n
!JM       pp(i,i,:) = _ZERO_
       pp(:,i,i) = _ZERO_
     end do

   return
   end subroutine reset_diagonal
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocate_bfm
!
! !INTERFACE:
        subroutine allocate_memory_bfm(nlev)
! !USES:
       use bio_var,only:bio_setup,numc_diag,diag,ccb,ppb,ddb,benvar_type,diagb,&
                    numc,numbc,numbc_diag, &
                    cc_before_transport
       use bfm_output,only:stBenFluxS,stBenFluxE,var_ids
       use mem,only:flx_option
#ifdef INCLUDE_DIAGNOS_PRF
      use bio_var,only:numbc_prf,diagb_prf,nprf
#endif
! !INPUT PARAMETERS:
        implicit none
        integer,intent(IN)            ::nlev
!
! !LOCAL VARIABLES:
   integer                   :: rc=0,n,i
   logical                   :: with_sedimentation_flux=.false.
!
! !DESCRIPTION:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  28-04-2006  Piet Ruardij Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

   if ( numc_diag > 0 ) then
     allocate(diag(0:nlev,1:numc_diag),stat=rc)
    if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (cc)'
     diag=_ZERO_                                                     !BFM
   endif


   !fi sedimentation flux is calculated.. a copy of the cc-array is needed
!  if ( with_sedimentation_flux ) then
     allocate(cc_before_transport(1:numc,0:nlev),stat=rc)
     if (rc /= 0) STOP 'init_bio: Error allocating (cc_before_transport)'
     cc_before_transport=_ZERO_
!  endif

   if (bio_setup >= 2) then                                         !BFM
     ! allocate benthic state variables                             !BFM
!JM     allocate(ccb(1:numbc,0:1),stat=rc)                             !BFM
     allocate(ccb(0:1,1:numbc),stat=rc)                             !BFM
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ccb)'!BFM
!JM     allocate(ppb(1:numbc,1:numbc,0:1),stat=rc)                     !BFM
     allocate(ppb(0:1,1:numbc,1:numbc),stat=rc)                     !BFM
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ppb)'!BFM
!JHM     allocate(ddb(1:numbc,1:numbc,0:1),stat=rc)                     !BFM
     allocate(ddb(0:1,1:numbc,1:numbc),stat=rc)                     !BFM
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ppb)'!BFM

     ccb=_ZERO_                                                     !BFM
     ppb=_ZERO_                                                     !BFM
     ddb=_ZERO_                                                     !BFM

     ! allocate variable holding type and save attributes           !BFM
     allocate(benvar_type(1:numbc),stat=rc)                         !BFM
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (benvar_type)'!BFM
     benvar_type = 0

     if ( numbc_diag > 0 ) then
       allocate(diagb(0:1,1:numbc_diag),stat=rc)
       if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (diagb)'
       diagb=_ZERO_                                                 !BFM

#ifdef INCLUDE_DIAGNOS_PRF
      if (numbc_prf>0)  then
         allocate(diagb_prf(0:nprf,1:numbc_prf),stat=rc)
         if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (diagb_prf)'
         diagb_prf=_ZERO_                                           !BFM
      endif
#endif

     endif

   !---------------------------------------------
   ! Create pointers
   !---------------------------------------------
     call pointers_gotm_bfm()
     call AllocateMem
     call InitBoxParams

   end if

   i=0
   STDERR "stBenFluxS,stBenFluxE",stBenFluxS,stBenFluxE
   do n=stBenFluxS,stBenFluxE
     i=i+1
     if (var_ids(n) > 0.and.flx_option(i) ==20) with_sedimentation_flux=.true. 
   end do



 end subroutine allocate_memory_bfm
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test negative concentrations
!
! !INTERFACE:
       subroutine test_on_negative_states( statenr,lldeep, h, nlev, after, c1, &
                                       error, msg,counter )
!
! !USES:
       use bfm_output,only:var_names
       use gotm_error_msg, only:set_warning_for_getm
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer,intent(in)                      :: statenr
       integer,intent(in)                      :: nlev
       logical,intent(in)                      :: lldeep
       REALTYPE,intent(in)                     :: h(0:nlev)
       character(len=*),intent(IN)             :: after
!
! !OUTPUT PARAMETERS:
       REALTYPE,intent(inout)                  :: c1(0:nlev)
       integer,intent(OUT)                     :: error
       character(len=*),intent(OUT)            :: msg
       integer,intent(INOUT),optional          :: counter
!          Array c1 is modified if ncecessry
!
! !LOCAL VARIABLES:
        integer              ::k
        integer              ::i,n
        REALTYPE             ::r,s
        REALTYPE             ::sumbefore,sumafter
        character(len=22)    ::onem

! !DESCRIPTION:
!   Routine to check for negative values.
!   Negative values are corrected with the aveage of neighbouring
!   grid points. A warning is given.
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!      created by P. Ruardij 21-06-2006
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
       error=0
       ! in this way NaN are directly found!
       if (minval(c1(1:nlev)) .ge. 0.00D+00) then           !BFM
          continue
       else
          sumbefore=sum(c1(1:nlev)*h(1:nlev))
          if ( .not.lldeep ) then
            s=_ZERO_;do i=1,nlev
              r=c1(i);if (r<_ZERO_) s=s-r*h(i)
            enddo
            r=max(_ZERO_,sumbefore/sum(h(1:nlev)))
            c1(1:nlev)=r
            sumafter=1.0D-80+sum(h(1:nlev))*r
            r=max(0.0,s/sumafter)
            ! recalculate percentage shift
            r =100.0D+00*s
            if (present(counter).and.r< 0.001D+00) then
              counter=counter+1
            elseif ( r > 1.0D-10) then
              onem=var_names(statenr); i=len_trim(onem)
              call set_warning_for_getm()
              write(msg,'(A,'': Negative values after '' ,A,''averaging n>=1 shift='',F10.3,''%'')') &
                    onem(1:i),after,r
                i=len_trim(msg)
                STDERR msg(1:i)
                call set_warning_for_getm()
             endif
          else
            k=0                                                !BFM
            n=0
            do i = 1,nlev                                      !BFM
              if ( c1(i).lt._ZERO_) then                   !BFM
                  write(onem,'(I2,'':'',G12.3,''/'')') i,c1(i)
                  if (index(onem,'-')>0 )n=n+1
                  k=-i                                          !BFM
                  if ( i == 1 ) then
                     if ( c1(i+1) > _ZERO_ )  then
                       c1(i)=0.1* c1(i+1)
                       k=i
                     else
                       c1(i)=_ZERO_
                       k=i
                     endif
                  elseif ( i == nlev ) then
                     if ( c1(i-1) > _ZERO_ )  then
                       c1(i)=0.1D+00* c1(i-1)
                       k=i
                     else
                       c1(i)=_ZERO_
                       k=i
                     endif
                  else if ( (c1(i-1) > _ZERO_) .and. ( c1(i+1)>_ZERO_ ) ) then
                     k=i
                     c1(i)=(c1(i-1)+c1(i+1)) * 0.1D+00
                  else if ( c1(i-1) >= _ZERO_ ) then
                       c1(i)=0.1* c1(i-1)
                       k=i
                  else if ( c1(i+1) >= _ZERO_ ) then
                       c1(i)=0.1* c1(i+1)
                       k=i
                  endif
               endif
               if ( error.ge.0) error=k
            end do                                    !BFM
            sumafter=1.0D-80+sum(c1(1:nlev)*h(1:nlev))
            r=max(_ZERO_,sumbefore/sumafter)
!           if ( r.eq.0.0) STDERR 'SUMFEFORE=', sumbefore ,sumafter
!           STDERR 'nega:',sumbefore,sumafter
            ! recalculate percentage shift
            c1=c1*min(1.0D+00,r);r =100.0D+00*(1.0-r)
            if (present(counter).and.r< 0.001D+00) then
              counter=counter+n
            elseif ( n > 0 ) then
              onem=var_names(statenr); i=len_trim(onem)
              write(msg,'(A,'': Negative values after '' ,A,'' n='',I2,'' shift='',F10.3,''%'')') &
                    onem(1:i),after,n,r
            endif
         endif
       endif

     end subroutine test_on_negative_states
!EOC
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test negative concentrations for all_states
!
! !INTERFACE:
       subroutine test_on_all_negative_states (lldo, lldeep, &
         b_setup, h, nstates,nlev, warning_level, messafter, ccx, enderror )
!
! !USES:
       use bio_var,only:pelvar_type,ALLTRANSPORT,c1dimz
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       logical,intent(in)                      :: lldo
       logical,intent(in)                      :: lldeep
       integer,intent(in)                      :: b_setup
       integer,intent(in)                      :: nstates
       integer,intent(in)                      :: nlev
       integer,intent(in)                      :: warning_level
       REALTYPE,intent(in)                     :: h(0:nlev)
       character(len=*),intent(IN)             :: messafter
!
! !OUTPUT PARAMETERS:
!          Array ccx is modified if ncecessry
!JM       REALTYPE,intent(inout)                  :: ccx(1:nstates,0:nlev)
       REALTYPE,intent(inout)                  :: ccx(0:nlev,1:nstates)
       integer,intent(OUT)                     :: enderror
! !LOCAL VARIABLES:
        integer              ::j,i
        character(len=180)   ::msg
        integer              :: counter
        integer              :: error
           error=0;enderror=0;
           if (lldo) then
              do j=1,nstates
               if (b_setup /= 2 ) then
                 if (pelvar_type(j)>=ALLTRANSPORT) then
!JM                   c1dimz=ccx(j,:)
                   c1dimz=ccx(:,j)
                   counter=0
                   call test_on_negative_states ( j,lldeep,h,nlev, messafter, &
                   c1dimz, error,msg,counter )
                   ! counter.gt.0: neglectable small errors
                   i=error*max(0,1-counter)
                   if ((pelvar_type(j)>=ALLTRANSPORT.and.warning_level>1) &
                     .and.i.ne.0) then
                      STDERR msg(1:len_trim(msg)) ;enderror=1;
                   endif
                   if (error.lt.0) then
                       enderror=error;return
                   endif
!JM                   if (error.gt.0) ccx(j,:)=c1dimz
                   if (error.gt.0) ccx(:,j)=c1dimz
                 endif
               endif
             enddo
           endif
     end subroutine test_on_all_negative_states


!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test mass conservation concentrations
!
! !INTERFACE:
!
       subroutine test_mass_conservation(lldo,mode,after,totsysn_old,totsysp_old)
! !USES:
       use gotm_error_msg, only:set_warning_for_getm
       use mem, only: totsysn,totsysp
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       logical,intent(in)                      :: lldo
       integer,intent(in)                      :: mode
       character(len=*),intent(IN)             :: after
!
!
! !INPUT/OUTPUT PARAMETERS:
       REALTYPE,intent(INOUT)                   :: totsysn_old
       REALTYPE,intent(INOUT)                   :: totsysp_old
!
! !LOCAL VARIABLES:
        REALTYPE             ::r
        integer              ::i
        character(len=160)   ::msg=''

        if (lldo ) then
            if ( mode.le.2) then
              !initialize
              if ( mode.eq.1) call CheckMassConservationNPSDynamics
              totsysn_old=totsysn(1);totsysp_old=totsysp(1)
            else
              ! CheckMassConservation need only be called if this routine is called if state vars
              ! are changed by physical processes.
              if ( mode.eq.3) call CheckMassConservationNPSDynamics
              r= (totsysn(1)-totsysn_old)/abs(totsysn(1)) * 100.0
              if ( r> 0.0001) then
                  call set_warning_for_getm
                  write(msg,'(''No Mass conservation for N after call to '',A,'': '',F7.4,''%'')') &
                        after,r
                  i=len_trim(msg); STDERR msg(1:i)

              endif
              r= (totsysp(1)-totsysp_old)/abs(totsysp(1)) * 100.0
              if ( r> 0.0001) then
                 call set_warning_for_getm
                 write(msg,'(''No Mass conservation for P after call to '',A,'': '',F7.4,''%'')') &
                        after,r
                  i=len_trim(msg); STDERR msg(1:i)
              endif
              totsysn_old=totsysn(1);totsysp_old=totsysp(1)
            endif
        endif
    end subroutine test_mass_conservation


!EOC
!-------------------------------------------------------------------------




!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine calc_sigma_depth( nlev,ddu,maxdepth,arr)
!
! !DESCRIPTION:
!  Calculate Sigma depth.
!  This routine is used to calculate the profiles of the benthic
!  nutrient model.
!  This routine is a simplification of
!  the calculation used in gotm/getm
!
! !USES:
   IMPLICIT NONE
! !INPUT PARAMETERS:
     integer,intent(IN)           :: nlev
     REALTYPE,intent(IN)          :: ddu
     REALTYPE,intent(IN)          :: maxdepth
! !OUTPUT PARAMETERS:
     REALTYPE,intent(OUT)        :: arr(1:nlev)
! !LOCAL PARAMETERS:
     REALTYPE                    :: r,s
     integer                     :: i
!
!
! !REVISION HISTORY:
!  Original by Piet Ruardij
!
!EOP
!-----------------------------------------------------------------------
!BOC

   r =0.0
   do i=1,nlev
      s= maxdepth*(1.0 - tanh(ddu*float(nlev-i)/float(nlev))/tanh(ddu))
      arr(i)=(s+r) * 0.5
      r=s
   enddo
   return
   end subroutine calc_sigma_depth
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_bfm
!
! !DESCRIPTION:
!  Nothing done here --- supplied for completeness
!  with GOTM bio structure.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_bfm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DeriveFromGotm
!
! !INTERFACE:
   subroutine DeriveFromGotm(mode,nlev,arr)
!
! !DESCRIPTION:
!  Nothing done here --- supplied for completeness
!  with GOTM bio structure.
!
! !USES:
   use turbulence,   only: num,eps

   IMPLICIT NONE
! !INPUT PARAMETERS:
     integer,intent(IN)           :: mode
     integer,intent(IN)           :: nlev
! !OUTPUT PARAMETERS:
     REALTYPE,intent(OUT)        :: arr(1:nlev)
! !REVISION HISTORY:
     select case (mode)
      case(1)   ! Calculate Shear Rate (1/s)
          arr(1:nlev)= sqrt(eps(1:nlev)/(1.0D-80+num(1:nlev)))
     end select


!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine DeriveFromGotm
!EOC

   subroutine print_structure(mode,nr,layer)
   use mem,only:D3SOURCE,D3SINK,D2SOURCE,D2SINK, PELSURFACE,PELBOTTOM, &
                D3DIAGNOS,D2DIAGNOS,NO_D3_BOX_DIAGNOSS,NO_D2_BOX_DIAGNOSS, &
                NO_D2_BOX_STATES,NO_D3_BOX_STATES
    
   implicit none
   integer,intent(IN)          ::mode,nr,layer

   integer                     ::keep(20)
   integer,save                ::keep_o(20)
   integer                     ::i,n,ikeep,nra
   integer,save                ::ikeep_o=-1,icount=0,istep=0
   logical                     ::l
   ikeep=0
   istep=istep+1
   if (ikeep_o==-1) then
       ikeep_o=0
       keep_o=0
   endif
   select case (mode)
     case(2); n=NO_D2_BOX_STATES
     case(3); n=NO_D3_BOX_STATES
   end select

   nra=abs(nr)
   do i=1,n
     if ( nr>0) then
     select case (mode)
!JM       case(2); l=(abs(D2SOURCE(nra,i,layer))>1.0D-80)
!JM       case(3); l=(abs(D3SOURCE(nra,i,layer))>1.0D-80)
       case(2); l=(abs(D2SOURCE(layer,nra,i))>1.0D-80)
       case(3); l=(abs(D3SOURCE(layer,nra,i))>1.0D-80)
     end select
     if (l) then
       ikeep=ikeep+1;keep(ikeep)=i
     endif
     endif
     if ( nr<0) then
     select case (mode)
!JM       case(2); l=(abs(D2SINK(nra,i,layer))>1.0D-80)
!JM       case(3); l=(abs(D3SINK(nra,i,layer))>1.0D-80)
       case(2); l=(abs(D2SINK(layer,nra,i))>1.0D-80)
       case(3); l=(abs(D3SINK(layer,nra,i))>1.0D-80)
     end select
     if (l) then
       ikeep=ikeep+1;keep(ikeep)=-i
     endif
     endif
   enddo
     
   if (ikeep >0.and.ikeep.ne.ikeep_o) then
     write(stderr,'(''flows '',I3,'':'',(I7,20I4))') nra,istep, keep(1:ikeep)
   else
      l=.false.
      do i=1,ikeep
         l= (keep(i).ne.keep_o(i))
      enddo
     if (l) then
         write(stderr,'(''flows '',I3,'':'',(I7,20I4))') nra,istep,keep(1:ikeep)
     endif
   endif
   keep_o=keep;ikeep_o=ikeep
   end subroutine print_structure

   subroutine test_structure(test_3d,test_2d,test_part)
   use bio_var,only:test_mode
   use mem,only:D3SOURCE,D3SINK,D2SOURCE,D2SINK, PELSURFACE,PELBOTTOM, &
                D3DIAGNOS,D2DIAGNOS,NO_D3_BOX_DIAGNOSS,NO_D2_BOX_DIAGNOSS, &
                NO_D2_BOX_STATES,NO_D3_BOX_STATES

   implicit none
   integer,intent(IN)          ::test_3d,test_2d,test_part

   integer                     ::i3d,j3d,i2d,j2d
   integer,save                ::l=0,k,m,n
   real                        ::rk
   REALTYPE                    ::value

   test_mode=test_3d.ne.0.or.test_2d.ne.0
   if (l==0.and.test_mode)  &
    STDERR "test_mode: loss processes are NOT compensated with surface input"
   j3d=mod(test_3d,10);i3d=(test_3d-j3d)/10
   j2d=mod(test_2d,10);i2d=(test_2d-j2d)/10
   if (l==0) STDERR "test_structure",test_3d,test_2d,i3d,j3d,i2d,j2d
   if (test_part>=200) then
      m=mod(test_part,100);k=(test_part-m)/100;
   elseif (test_part<0) then
       k=test_part
   else
      m=0;k=0
   endif
   if ( j3d.ge.1.and.j3d.le.3) then
     if (k.ne.0) then
        if (k>0) then 
           n=int(m*rk);m=int((m-1)*rk)+1
           rk=real(NO_D3_BOX_STATES )/real(k)
        else ;k=-k;m=k;n=k; endif
!JM        if( j3d.eq.1.or.j3d.eq.3) D3SOURCE(m:n,:,:)=_ZERO_
!JM        if (j3d.ge.2) D3SINK(m:n,:,:)=_ZERO_
        if( j3d.eq.1.or.j3d.eq.3) D3SOURCE(:,m:n,:)=_ZERO_
        if (j3d.ge.2) D3SINK(:,m:n,:)=_ZERO_
        if (l.eq.0) STDERR "test_structure Reset D3SOURCE D3SINK",rk,m,n
      else
        if( j3d.eq.1.or.j3d.eq.3) D3SOURCE(:,:,:)=_ZERO_
        if (j3d.ge.2) D3SINK(:,:,:)=_ZERO_
        if (l.eq.0) STDERR "test_structure Reset D3SOURCE D3SINK"
      endif
   endif
   if ( j2d.ge.1.and.j3d.le.3) then
     if (k.ne.0) then
        if (k>0) then 
            n=int(m*rk);m=int((m-1)*rk)+1
            rk=real(NO_D2_BOX_STATES )/real(k)
        else ;k=-k;m=k;n=k; endif
!JM        if( j2d.eq.1.or.j2d.eq.3) D2SOURCE(m:n,:,:)=_ZERO_
!JM        if (j2d.ge.2) D2SINK(m:n,:,:)=_ZERO_
        if( j2d.eq.1.or.j2d.eq.3) D2SOURCE(:,m:n,:)=_ZERO_
        if (j2d.ge.2) D2SINK(:,m:n,:)=_ZERO_
        if (l.eq.0) STDERR "test_structure Reset D2SOURCE D2SINK",rk,m,n
      else
        if( j2d.eq.1.or.j2d.eq.3) D2SOURCE(:,:,:)=_ZERO_
        if (j2d.ge.2) D2SINK(:,:,:)=_ZERO_
        if (l.eq.0) STDERR "test_structure Reset D2SOURCE D2SINK"
      endif
   endif
   if ( i3d.eq.2.or.i3d.eq.3) then
     select case  (i3d)
       case (2);value=-1.0D+11
       case (3);value=1.0D+11
       case (4);value=_ZERO_
     end select
     if (k>0) then
        rk=real(NO_D3_BOX_DIAGNOSS )/real(k)
        n=int(m*rk);m=int((m-1)*rk)+1
        D3DIAGNOS=_ZERO_;
        D3DIAGNOS(:,m:n)= value
        if (l.eq.0) STDERR "test_structure Reset D3DIAGNOS on ",value,rk,m,n
     else
        D3DIAGNOS=value
        if (l.eq.0) STDERR "test_structure Reset D3DIAGNOS on ",value
     endif
   endif
   if ( i2d.eq.2.or.i2d.eq.3) then
     select case  (i2d)
       case (2);value=-1.0D+11
       case (3);value=1.0D+11
       case (4);value=_ZERO_
     end select
     if (k>0) then
        rk=real(NO_D2_BOX_DIAGNOSS )/real(k)
        n=int(m*rk);m=int((m-1)*rk)+1
        if (NO_D2_BOX_DIAGNOSS-k<n)n=NO_D2_BOX_DIAGNOSS
        D2DIAGNOS=_ZERO_;
        D2DIAGNOS(:,m:n)= value
        if (l.eq.0) STDERR "test_structure Reset D2DIAGNOS on ",value,m,n
     else
        D2DIAGNOS=value
        if (l.eq.0) STDERR "test_structure Reset D2DIAGNOS on ",value
     endif
   endif
   l=1;
   end subroutine test_structure

   subroutine test_model_states(mode,nsel,numc,kmax,cc,h)

   implicit none
   integer,intent(IN)                                 ::mode
   integer,intent(IN)                                 ::nsel
   integer,intent(IN)                                 ::numc
   integer,intent(IN)                                 ::kmax
   REALTYPE,intent(IN),dimension(1:numc,0:kmax)       ::cc
   REALTYPE,intent(IN),dimension(0:kmax),optional      ::h

   integer           ::i,from,ito
   REALTYPE          ::r
   REALTYPE,save     ::old_fr=-1.0
   integer,save      ::old_exp=-1
   REALTYPE          ::n_fr
   integer           ::n_exp
   n_fr=_ZERO_; n_exp=_ZERO_
   if (nsel.eq.0) then
     from=1 ; ito=numc;
   else
     from=nsel;ito=nsel;
   endif
   do i=from,ito
     if ( kmax.gt.1) then
 !JM       r=sum(cc(i,1:kmax)*h(1:kmax))
        r=sum(cc(1:kmax,i)*h(1:kmax))
    else
!JM        r=cc(i,1)
        r=cc(1,i)
     endif
     n_fr=n_fr+ fraction(r)
     n_exp=n_exp +exponent(r)
   enddo
   if (old_exp>0.and.mode.gt.0) then
     if ( n_exp.ne.old_exp) STDERR 'test_model_states change in exp mode=',mode
     if ( abs(n_fr -old_fr)/old_fr > 1.0D-06 ) STDERR 'test_model_states change in frac'
   endif
   old_fr=n_fr; old_exp=n_exp
   end subroutine test_model_states

   subroutine test_model_rates(mode,nsel,numc,kmax,h)
   use bio_var, only: bio_setup,pp,dd,numbc,ccb,ppb,ddb,c1dimz
   implicit none
   integer,intent(IN)                                 ::mode
   integer,intent(IN)                                 ::nsel
   integer,intent(IN)                                 ::numc
   integer,intent(IN)                                 ::kmax
   REALTYPE,intent(IN),dimension(0:kmax),optional      ::h


   integer           ::i,from,ito
   REALTYPE          ::r
   REALTYPE,save     ::old_fr=-1.0
   integer,save      ::old_exp=-1
   REALTYPE          ::n_fr
   integer           ::n_exp
   n_fr=_ZERO_; n_exp=_ZERO_
   if (nsel.eq.0) then
     from=1 ; ito=numc;
   else
     from=nsel;ito=nsel;
   endif
   do i=from,ito
     if ( kmax.gt.1) then
!JM        r=sum(sum(pp(i,:,1:kmax)-dd(i,:,1:kmax),2)*h(1:kmax))
        r=sum(sum(pp(1:kmax,i,:)-dd(1:kmax,i,:),2)*h(1:kmax))
     else
!JM        r=sum(ppb(i,:,1)-ddb(i,:,1))
        r=sum(ppb(1,i,:)-ddb(1,i,:))
     endif
     n_fr=n_fr+ fraction(r)
     n_exp=n_exp +exponent(r)
   enddo
   if (old_exp>0.and.mode.gt.0) then
     if ( n_exp.ne.old_exp) STDERR 'test_model_rates change in exp mode=',mode
     if ( abs(n_fr -old_fr)/old_fr > 1.0D-06 ) STDERR 'test_model_rates change in frac'
   endif
   old_fr=n_fr; old_exp=n_exp
   end subroutine test_model_rates

!-----------------------------------------------------------------------


   end module bio_bfm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License
! www.gnu.org
!-----------------------------------------------------------------------
