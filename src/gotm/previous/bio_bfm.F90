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
          test_on_negative_states, end_bio_bfm,       &
          do_bfm_river_loads,assign_adv_rates,        &
          CalcVertFluxAtLev,calc_sigma_depth,         &
          DeriveFromGotm,return_bfm,test_structure,   &
          test_mass_conservation, &
          print_structure,output_warning_shiftings, &
          do_nan_test_after_adv,set_pointers_from_gotm_to_bfm,&
          warning_test_mass_conservation

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
   subroutine init_bio_bfm(nlev)
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
   nprf= NO_D3_BOX_DIAGNOSS_PRF
   numbc_prf= NO_BOXES_PRF
#endif

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

! call to init_var_bfm is found  in bio.F90
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
   use bio_var,only:cc,pp,dd,pelvar_type,diag, &
                ccb,ppb,ddb,benvar_type,diagb,numc_diag,numbc_diag
   use mem, only: D3STATE,D3SOURCE,D3SINK,D3STATETYPE, &
                  D3DIAGNOS,D2STATE,D2SOURCE,D2SINK,   &
                  D2STATETYPE,NO_BOXES,NO_BOXES_XY,    &
                  D2DIAGNOS
#ifdef INCLUDE_DIAGNOS_PRF
      use mem,only:D3DIAGNOS_PRF,NO_BOXES_PRF
      use bio_var,only:numbc_prf,diagb_prf
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!

   !---------------------------------------------
   ! Fixed Pelagic pointers
   !---------------------------------------------
   D3STATETYPE => pelvar_type
   if (numc_diag > 0) D3DIAGNOS=> diag(:,1:NO_BOXES)
#ifdef INCLUDE_DIAGNOS_PRF
   if (numbc_prf>0) D3DIAGNOS_PRF => diagb_prf(1:NO_BOXES_PRF,:)
#endif
   !---------------------------------------------
   ! Fixed Benthic pointers
   !---------------------------------------------
   D2STATETYPE => benvar_type
   if (numbc_diag>0) D2DIAGNOS => diagb(:,1:NO_BOXES_XY)

   !---------------------------------------------
   ! Pelagic pointers
   !---------------------------------------------
   D3STATE  => cc(1:NO_BOXES,:)
   D3SOURCE => pp(1:NO_BOXES,:,:)
   D3SINK   => dd(1:NO_BOXES,:,:)
   !---------------------------------------------
   ! Benthic pointers
   !---------------------------------------------
   D2STATE  => ccb(1:NO_BOXES_XY,:)
   D2SOURCE => ppb(1:NO_BOXES_XY,:,:)
   D2SINK   => ddb(1:NO_BOXES_XY,:,:)
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
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call setup_bio_output
   call set_var_info_bfm
   return
   end subroutine var_info_bfm
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine set_pointers_from_gotm_to_bfm( &
        check,nlev,dt,h,t,s,rho,nuh)
!
! !DESCRIPTION
!
! !USES
! BFM modules
 use bio_var, only: nlev_local=>nlev,dt_local=>dt,h_l,nuh_l
 use mem,ONLY: ESS, ERHO, ETW,ESW, R9x,Depth
 use bfm_output,only:unselect_var_from_output

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(IN)                  :: check
   integer,intent(IN)                  :: nlev
   REALTYPE, intent(in)                :: dt
   REALTYPE, intent(in),target         :: h(0:nlev)
   REALTYPE, intent(in),target         :: t(0:nlev)
   REALTYPE, intent(in),target         :: rho(0:nlev)
   REALTYPE, intent(in),target         :: s(0:nlev)
   REALTYPE, intent(in),target         :: nuh(0:nlev)
!
! !OUTPUT PARAMETERS:
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !---------------------------------------------
   ! Assign depths of layers
   ! temperature and salinity
   !---------------------------------------------
    nlev_local=nlev
    dt_local=dt
    Depth=> h(1:nlev)
    ETW => t(1:nlev)
    ESW => s(1:nlev)
    ERHO => rho(1:nlev)
    ESS=>R9x(:)
    nuh_l=>nuh(0:nlev)
    h_l=>h(0:nlev)
    if ( check==1) then
      call unselect_var_from_output('Depth')
      call unselect_var_from_output('EUCURR_LEVEL')
      call unselect_var_from_output('EVCURR_LEVEL')
      call unselect_var_from_output('ETW')
      call unselect_var_from_output('ESW')
      call unselect_var_from_output('ERHO')
      call unselect_var_from_output('ESS')
    endif
   !---------------------------------------------
   ! Notice that irradiance in the BFM is in
   ! uE/m2/s and is defined at the top of each
   ! layer (the derivation of the middle-layer
   ! EIR for production is done in the
   ! Phytoplankton routines)
   !---------------------------------------------
   end subroutine set_pointers_from_gotm_to_bfm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine set_env_bio_bfm(nlev,bathy_dep,u_taub,dl, &
          uwind_gotm,vwind_gotm,ulevel,vlevel,I_0,julianday,dry_z_getm)
! !DESCRIPTION
!
! !USES
! BFM modules
   use bio_var,only: bathy_dep_local=>bathy_dep,bio_julianday=>julianday, &
                                                         I_0_local=>I_0
   use mem,only:Wind,EUWIND,EVWIND,EUCURR_LEVEL,EVCURR_LEVEL, &
     ETAUB,dry_z,EIRr,SUNQ,OCDepth,Depth,NO_BOXES

   implicit none
   integer,intent(in)             :: nlev
   REALTYPE, intent(in)           :: bathy_dep
   REALTYPE, intent(in)           :: u_taub
   REALTYPE, intent(in)           :: dl
   REALTYPE, intent(in)           :: I_0
   REALTYPE, intent(in)           :: uwind_gotm,vwind_gotm
   REALTYPE, intent(in)           :: ulevel,vlevel
   integer(8),intent(in)          :: julianday
   REALTYPE,intent(in),optional   :: dry_z_getm

   integer             :: n
   bio_julianday=julianday

   EUWIND(1)=uwind_gotm
   EvWIND(1)=vwind_gotm
   Wind= sqrt( EUWIND(1) * EUWIND(1) + EVWIND(1) *EVWIND(1))
   ETAUB(1)=u_taub
   dry_z(1)=_ONE_
   if (present(dry_z_getm)) dry_z(1)=dry_z_getm
   bathy_dep_local = bathy_dep
   I_0_local= I_0
   EIRr=I_0
   SUNQ=dl
!JM pass currents
    EUCURR_LEVEL=ulevel
    EVCURR_LEVEL=vlevel

   OCDepth(NO_BOXES) =depth(NO_BOXES)
   do n=NO_BOXES-1,1,-1
     OCDepth(n)=depth(n)+ OCDepth(n+1)
   enddo
   end subroutine set_env_bio_bfm
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the BFM model
!
! !INTERFACE
   subroutine do_bio_bfm(first,check_rates)
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!
! !USES
   use bio_var,only:bio_setup,surface_flux_method,sfl,bfl,sfl_N3n,sfl_N4n, &
                    surface_flux_method,numc ! ,diagb
!  use bio,only: ModelCallNr
   use mem, only:     iiC,iiN,iiP,iiS,iiL, &
       NO_BOXES_XY,ppO3c,ppO2o,ppN1p,ppN3n,ppN4n,ppN5s,Depth, &
       PELBOTTOM, PELSURFACE, ESW ! ,dry_z

   use mem_param,only:AssignAirPelFluxesInBFMFlag,AssignPelBenFluxesInBFMFlag, &
     use_function_for_eps0,p_eps0

   use constants,  only: SEC_PER_DAY
   use gotm_error_msg, only:get_d3_model_flag_from_getm
   use global_interface,only:DefinitiveLossGain
   IMPLICIT NONE
!
   logical,intent(in)          :: first
   logical,intent(in)          :: check_rates


   logical                     :: d3_model_flag
   logical,save                :: start=.true.
   integer                     :: k
   REALTYPE                    :: topm3psec
   REALTYPE,dimension(NO_BOXES_XY) :: ChangePel,ChangeBen
 
   interface
     function CompensateLoss(ppnutrient,valuePel,valueBen,multiply)
     use global_mem,only:RLEN
     use mem,only: NO_BOXES_XY
     integer,intent(IN)                           :: ppnutrient
     real(RLEN),dimension(NO_BOXES_XY),intent(IN) :: valuePel
     real(RLEN),dimension(NO_BOXES_XY),intent(IN) :: valueBen
     integer,intent(IN)                           :: multiply
     real(RLEN)                                   :: CompensateLoss
     end function CompensateLoss
   end interface

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
   !---------------------------------------------
   ! Compute extinction coefficient
   !---------------------------------------------

#ifdef SPM
     p_eps0=0.0D+00
!    ABIO_eps(1:nlev) = abioshade(1:nlev)
#else
     select case (use_function_for_eps0)
       !northsea
       case (1) ; p_eps0=1.17692308D+00-0.0307692308D+00*min(35.0D+00,ESW(1))
       !wz
       case (2) ; p_eps0=1.17692308D+00-0.0307692308D+00 &
                                          &*max(30.0D+00,min(35.0D+00,ESW(1)))
     end select
#endif

   !---------------------------------------------
   ! Compute BFM terms
   !---------------------------------------------
   if (.not.check_rates)call SiltDynamics
   call EcologyDynamics
   call get_d3_model_flag_from_getm(d3_model_flag)
   ChangePel=0.0D+00;ChangeBen=0.0D+00
   if ( .not.(d3_model_flag.or.start)) then
      call DefinitiveLossGain(iiN,ChangePel,ChangeBen)
   endif
   start=.false.
   !---------------------------------------------
   ! Surface fluxes
   !---------------------------------------------
   if ( bio_setup ==2 ) return
   topm3psec=_ONE_/SEC_PER_DAY
   sfl=_ZERO_
   if ( .NOT. AssignAirPelFluxesInBFMFlag ) then
     sfl(ppO2o) =   PELSURFACE(ppO2o,1) *topm3psec
     if ( ppO3C > 0 ) sfl(ppO3c) =   PELSURFACE(ppO3c,1) *topm3psec
   endif
   select case (mod(surface_flux_method,10))
     case (-1)! absolutely nothing
     case (0) ! constant
       sfl(ppN3n) =   0.12D+00  *topm3psec
       sfl(ppN4n) =   0.09D+00  *topm3psec
       !this is called here to test track when the 1d model is used in 1D-mode.
       !In this case get sfl(pptrN3n) the same value as slf(ppN3n)
       ! It work only if d3_model_flag ==flase and if trakcing is active.
!       call fill_sfl(d3_model_flag,ppN3n,numc,sfl)
!       call fill_sfl(d3_model_flag,ppN4n,numc,sfl)
       sfl(ppN1p) =   _ZERO_  !0.0
     case (1) ! from file via sfl_read
       ! fluxes are in mmol m-2 d-1
       sfl(ppN3n) =    sfl_N3n  *topm3psec
       sfl(ppN4n) =    sfl_N4n  *topm3psec
!      call fill_sfl(d3_model_flag,ppN3n,numc,sfl)
!       call fill_sfl(d3_model_flag,ppN4n,numc,sfl)
     case (3) ! sfl array filled externally - for 3D models
       ! option 3 works only in 1D-mode!!!!!!!
        sfl(ppN3n)= CompensateLoss(ppN3n,ChangePel,ChangeBen, &
                                            (surface_flux_method-3)/10)
        sfl(ppN3n)=sfl(ppN3n)  * topm3psec
        PELSURFACE(ppN3n,1)=(ChangePel(1)+ChangeBen(1))
        call DefinitiveLossGain(iiP,ChangePel,ChangeBen)
        PELSURFACE(ppN1p,1)=(ChangePel(1)+ChangeBen(1))
        sfl(ppN1p)= (ChangePel(1)+ChangeBen(1))* topm3psec 
        call DefinitiveLossGain(iiS,ChangePel,ChangeBen)
        PELSURFACE(ppN5s,1)=(ChangePel(1)+ChangeBen(1))
        sfl(ppN5s)= (ChangePel(1)+ChangeBen(1))* topm3psec 
!       if ( ppO3h> 0) sfl(ppO3h)= Hloss * topm3psec
!       call flux(iiPel,NO_BOXES_Z,ppO3h,ppO3h,Hloss/Depth(NO_BOXES_Z))
     case default
   end select

   !---------------------------------------------
   ! Bottom fluxes
   !---------------------------------------------
   topm3psec=1.0/Depth(1)/ SEC_PER_DAY
   bfl=_ZERO_
   if ((bio_setup == 3 ) .and. ( .NOT.AssignPelBenFluxesInBFMFlag)) then

      do k=1,numc
        bfl(k)=PELBOTTOM(k,1)*topm3psec
      enddo
   endif
   end subroutine do_bio_bfm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Assign  adv_rates
!
! !INTERFACE
!  subroutine CalcVertFluxAtLev(statenr,lev,nlev,Depth,dt,out)
!
! !DESCRIPTION
! !USES
!  use bio_var,only: cc_before_transport,cc
!  use mem,only: PELBOTTOM
!  use constants,  only: SEC_PER_DAY
!  IMPLICIT NONE
!
!  integer,intent(IN)          :: statenr
!  integer,intent(IN)          :: lev
!  integer,intent(IN)          :: nlev
!  REALTYPE,intent(IN)         :: Depth(1:nlev)
!  REALTYPE,intent(IN)         :: dt
!  REALTYPE,intent(OUT)        :: out

!  out=0.0D+00
!  if ( lev > 0.and.allocated(cc_before_transport) ) then
!      ! rate calculate in time unit of physical model (secs)
!      out= sum((cc_before_transport(statenr,lev+1:nlev) &
!              -cc(statenr,lev+1:nlev))*Depth(lev+1:nlev))/dt
!  elseif ( lev==0 ) then
!      ! rate calculate transferred from time unit of eco model
!      ! to the one  of physical model (secs)
!      out=PELBOTTOM(statenr,1)/SEC_PER_DAY
!  else
!     out=_ZERO_
!  endif
!  return
!  end subroutine CalcVertFluxAtLev
!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Assign  adv_rates
!
! !INTERFACE
   subroutine assign_adv_rates
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!
! !USES
   use bio_var,only:bio_setup,ws,llws,c1dimz, &
               pelvar_type,SILTTRANSPORT
   use mem, only: NO_BOXES_Z,NO_D3_BOX_STATES,Depth,ppR6c
   use mem_BFM_Sinking,only:CopyInfoOnSinkingtoGOTM
   use bfm_output,only:var_names
   use constants,  only: SEC_PER_DAY

   IMPLICIT NONE
!

   logical                     :: ll_larger
   logical,save                :: first=.true.
   integer                     :: i,j,ldep,k
   REALTYPE                    :: corr(1:NO_BOXES_Z)
   character(len=22)           ::onem,pnem

   !---------------------------------------------
   ! Transfer sinking velocities (m/d -> m/s)
   !---------------------------------------------
   if ( bio_setup ==2 ) return
   ldep=1
   do j=1,NO_D3_BOX_STATES
     ! sedi vars in BFM copied in c1dimz and  reunited from m/d to sec/d
     call CopyInfoOnSinkingtoGOTM(j,i,c1dimz)
     if ( i > 0 ) then
       ll_larger=(maxval(abs(c1dimz(1:NO_BOXES_Z)))>1.0D-7)
       ! Reclculate from the average in a layer to the
       ! sinking rates at boundarys between the layers
       if ( ll_larger) then
         c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*c1dimz(1:NO_BOXES_Z-1)&
           +Depth(1:NO_BOXES_Z-1)*c1dimz(2:NO_BOXES_Z)) &
           /(Depth(1:NO_BOXES_Z-1)+Depth(2:NO_BOXES_Z-1))
         c1dimz(NO_BOXES_Z)=0.0D+00
         !sedimentation rate biological constituents  will
         !always be lower than of silt
         corr=min(0.002D+00, &
              abs(c1dimz(1:NO_BOXES_Z)))/ (1.0D-80+abs(c1dimz(1:NO_BOXES_Z)))
         ws(1:NO_BOXES_Z,j) = -c1dimz(1:NO_BOXES_Z)*corr
       else
        ws(1:NO_BOXES_Z,j) = 0.0D+00
       endif
       llws(j)=ll_larger
       ws(0,j)= ws(1,j)
     elseif (i< -NO_D3_BOX_STATES) then
       llws(j)=.false.
     elseif (-i==j) then
       ! only when iiPelSinkRef is equal to the negative value of its
       ! is possible to define straight a ws .
       ws(0,j)= ws(1,j)
       llws(j)=.true.
     elseif (i<0) then
       if ( -i>=j) then
         STDERR "i=",i
         stop 'error:assign_adv_rates'
       endif
       i=-i
       ws(0:NO_BOXES_Z,j) =ws(0:NO_BOXES_Z,i)
       llws(j)=llws(i)
     else
       llws(j)=.false.
     endif
   enddo

   if (first) then
     ldep=0;k=0
     do j=1,NO_D3_BOX_STATES
       if (k==0) then
           k=1 ; LEVEL3 "Vertical advective transport is present for"
       endif
       call CopyInfoOnSinkingtoGOTM(j,i,c1dimz)
       onem=var_names(j); ldep=len_trim(onem)
       if ( i > 0 ) then
         LEVEL4 onem(1:ldep)
       elseif (i <-NO_D3_BOX_STATES) then
          LEVEL4 onem(1:ldep)," No vertical advective transport"
       elseif (-i==j) then
          LEVEL4 onem(1:ldep)," defined outside CoupleInfoOnNSinkingToGotm.F90"
       elseif (i<0) then
          pnem=var_names(-i)
          LEVEL4 onem(1:ldep)," coupled to ",pnem(1:len_trim(pnem))
       endif
       if (pelvar_type(j)>=SILTTRANSPORT) LEVEL4 onem(1:ldep), &
           " Diffusive transport coupled to advective transport"
     enddo
     if ( i.ne.0) then
        LEVEL3 "End of Vertical transport";LEVEL3 " "
     endif
     first=.false.
   endif

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
       pp(:,i,i) = _ZERO_
     end do

   return
   end subroutine reset_diagonal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocate_bfm
!
! !INTERFACE:
   subroutine allocate_memory_bfm(nlev)
! !USES:
     use bio_var,only:numc_diag,diag,ccb,ppb,ddb,benvar_type,diagb,&
        numc,numbc,numbc_diag, cc_before_transport
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
     allocate(diag(1:numc_diag,0:nlev),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (cc)'
     diag=_ZERO_                                                     !BFM
   endif

   !fi sedimentation flux is calculated.. a copy of the cc-array is needed
!  if ( with_sedimentation_flux ) then
     allocate(cc_before_transport(0:nlev,1:numc),stat=rc)
     if (rc /= 0) STOP 'init_bio: Error allocating (cc_before_transport)'
     cc_before_transport=_ZERO_
!  endif

   ! allocate benthic state variables                             !BFM
   allocate(ccb(0:1,1:numbc),stat=rc)                         !BFM
   if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ccb)'!BFM
   allocate(ppb(0:1,1:numbc,1:numbc),stat=rc)                     !BFM
   if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ppb)'!BFM
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
     allocate(diagb(1:numbc_diag,0:1),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (diagb)'
     diagb=_ZERO_                                                 !BFM
   endif

#ifdef INCLUDE_DIAGNOS_PRF
    if (numbc_prf>0)  then
      allocate(diagb_prf(0:numbc_prf,1:nprf),stat=rc)
      if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (diagb_prf)'
      diagb_prf=_ZERO_                                           !BFM
    endif
#endif


   !---------------------------------------------
   ! Create pointers
   !---------------------------------------------
   call pointers_gotm_bfm()
   call Allocate3dMem
   call Allocate2dMem
   call AllocateOtMem
   call InitBoxParams

   i=0
   do n=stBenFluxS,stBenFluxE
     i=i+1
     if (var_ids(n) > 0.and.flx_option(i) ==20) with_sedimentation_flux=.true.
   end do



 end subroutine allocate_memory_bfm
!EOC
!-------------------------------------------------------------------------
     function do_nan_test_after_adv(statenr,nlev,c0,c1,ws)
! !USES:
       use bfm_output,only:var_names
!      use gotm_error_msg, only:set_warning_for_getm
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer,intent(in)                      :: statenr
       integer,intent(in)                      :: nlev
       REALTYPE,intent(in)                     :: c0(0:nlev)
       REALTYPE,intent(in)                     :: c1(0:nlev)
       REALTYPE,intent(in)                     :: ws(0:nlev)
       integer                                 :: do_nan_test_after_adv
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
        integer              ::k,error
        integer              ::i,n
        character(len=22)    ::onem

        k=0;error=0
        do n=nlev,0,-1
          if (isnan(c1(n)).or.isnan(ws(n))) then
            k=n
            onem=var_names(statenr);i=len_trim(onem)
            STDERR "do_nan_test_after_adv:"
            STDERR "statevar,layer,ws(k):",onem(1:i),k,ws(k)
            STDERR "old - new-value:",c0(k),c1(k)
!           write(msg,"(""Nan in element "",I2,""of statevar "",A10)") &
!                                                               k,onem(1:i)
            error=-1
            STDERR 'ws=',ws(1:nlev)
            STDERR onem(1:i),c0(1:nlev)
            return
          endif
        enddo
        do_nan_test_after_adv=error
        end function do_nan_test_after_adv
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
          if (.not.lldeep ) then
            s=_ZERO_;do i=1,nlev
              r=c1(i);if (r<_ZERO_) s=s-r*h(i)
            enddo
            r=max(_ZERO_,sumbefore/sum(h(1:nlev)))
            c1(1:nlev)=r
            sumafter=1.0D-80+sum(h(1:nlev))*r
            r=max(0.0D+00,s/sumafter)
            ! recalculate percentage shift
            r =100.0D+00*s
            if (present(counter).and.r< 0.001D+00) then
              counter=counter+1
            elseif ( r > 1.0D-10) then
              onem=var_names(statenr)
              write(msg,'(A,'': Negative values after '' &
                &,A,''averaging n>=1 shift='',F10.3,''%'')') &
                onem(1:len_trim(onem)),after,r
              STDERR msg(1:len_trim(msg))
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
            ! recalculate percentage shift
            c1=c1*min(1.0D+00,r);r =100.0D+00*(1.0-r)
            if (present(counter).and.r< 0.001D+00) then
              counter=counter+n
            elseif ( n > 0 ) then
              onem=var_names(statenr)
              write(msg,'(A,'': Negative values after '' &
                    &,A,'' n='',I2,'' shift='',F10.3,''%'')') &
                    onem(1:len_trim(onem)),after,n,r
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
       REALTYPE,intent(inout)                  :: ccx(0:nlev,1:nstates)
       integer,intent(OUT)                     :: enderror
! !LOCAL VARIABLES:
        integer              ::j,i
        character(len=180)   ::msg
        integer              :: counter
        integer              :: error
           error=0;enderror=0
           if (lldo) then
              do j=1,nstates
               if (b_setup /= 2 ) then
                 if (pelvar_type(j)>=ALLTRANSPORT) then
                   c1dimz=ccx(:,j)
                   counter=0
                   call test_on_negative_states ( j,lldeep,h,nlev, messafter, &
                   c1dimz, error,msg,counter )
                   ! counter.gt.0: neglectable small errors
                   i=error*max(0,1-counter)
                   if ((pelvar_type(j)>=ALLTRANSPORT.and.warning_level>1) &
                     .and.i.ne.0) then
                      STDERR msg(1:len_trim(msg)) ;enderror=1
                   endif
                   if (error.lt.0) then
                       enderror=error;return
                   endif
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
   subroutine test_mass_conservation(lldo,mode,after,dt,totsys_old)
! !USES:
   use mem, only: totsysc,totsysn,totsysp,totsyss,totsysh,NO_BOXES_XY, &
   jtBENPELc,jtBENPELn,jtBENPELp,jtBENPELs,iiN,iiC,iiP,iiS, &
   iiBen,iiPel,totbenc,totbenn,totbenp,totbens,totbenh,dry_point
   use global_interface,only:DefinitiveLossGain
   use mem_CheckMassConservation,only: &
          CheckMassConservationNPS,CheckMassConservationC
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in)                      :: lldo
   integer,intent(in)                      :: mode
   character(len=*),intent(IN)             :: after
   REALTYPE,intent(IN)                     :: dt
!
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(:,:,:),intent(INOUT)    :: totsys_old
!
! !LOCAL VARIABLES:
   REALTYPE                        ::r,t,oldmass
   REALTYPE,dimension(NO_BOXES_XY) ::ChangePel,ChangeBen
   integer                         ::iXY
   integer,parameter               ::kkC=1,kkN=2,kkP=3,kkSi=4,kkH=5
!  character(len=160)   ::msg=''
   integer              :: yes,yesT

   if (lldo ) then
     if ( mode.le.2) then
       !initialize
       if (mode.eq.1) then
         call CheckMassConservationNPS;call CheckMassConservationC
       endif
       do iXY=1,NO_BOXES_XY
        totsys_old(kkC,1,iXY)=totsysc(iXY);totsys_old(kkC,2,iXY)=totbenc(iXY)
        totsys_old(kkN,1,iXY)=totsysn(iXY);totsys_old(kkN,2,iXY)=totbenn(iXY)
        totsys_old(kkP,1,iXY)=totsysp(iXY);totsys_old(kkP,2,iXY)=totbenp(iXY)
        totsys_old(kkSi,1,iXY)=totsyss(iXY);totsys_old(kkSi,2,iXY)=totbens(iXY)
        totsys_old(kkH,1,iXY)=totsysh(iXY);totsys_old(kkH,2,iXY)=totbenh(iXY)
       enddo
     else
       ! CheckMassConservation need only be called if this routine
       ! is called if state vars are not changed by physical processes.

       if (mode.eq.3) then
         call CheckMassConservationNPS;call CheckMassConservationC
       endif
       call DefinitiveLossGain(iiC,changePel,changeBen)
       yesT=0
       do iXY=1,NO_BOXES_XY
         yes=0
         if (dry_point(iXY).gt.1.0D-80) then
           r= calc_diff_mass_conservation(totsys_old(kkC,1,iXY) &
           -totsys_old(kkC,2,iXY),totsysc(iXY)-totbenc(iXY),changeBen(iXY),dt)
           call warning_test_mass_conservation(r,after,"(dry) Pelagic C",yes)
           yesT=yes
           if (yes==2) then
             r= calc_diff_mass_conservation(totsys_old(kkC,2,iXY), &
             totbenc(iXY), jtBENPELc(iXY)+changeBen(iXY),dt)
             call warning_test_mass_conservation(r,after,"Ben C",yes)
             if (yes==2) call OutputAfterConservationTest( &
                                  yes,totsys_old(kkC,2,iXY),r,iiBen,iiC)
           endif
         endif
       enddo

       if (yesT==0) then
       do iXY=1,NO_BOXES_XY
         yes=0
         r= calc_diff_mass_conservation(totsys_old(kkC,1,iXY),totsysc(iXY), &
         changePel(iXY)+changeBen(iXY),dt)
         call warning_test_mass_conservation(r,after,"C",yes)
         yesT=yes
         if (yes==2) then
           r= calc_diff_mass_conservation(totsys_old(kkC,2,iXY), &
           totbenc(iXY), jtBENPELc(iXY)+changeBen(iXY),dt)
           call warning_test_mass_conservation(r,after,"Ben C",yes)
           if (yes==2) call OutputAfterConservationTest( &
                                  yes,totsys_old(kkC,2,iXY),r,iiBen,iiC)
         endif
         if (yesT==2.and.yes==0) then
           oldmass=totsys_old(kkC,1,iXY)-totsys_old(kkC,2,iXY)
           r= calc_diff_mass_conservation(oldmass,totsysc(iXY)-totbenc(iXY),&
              changePel(iXY)-jtBENPELc(iXY),dt)
           call warning_test_mass_conservation(r,after,"Pel C",yes)
           if (yes==2)call OutputAfterConservationTest(yes,oldmass,r,iiPel,iiC)
         endif
       enddo
       endif

       call DefinitiveLossGain(iiN,changePel,changeBen)
       do iXY=1,NO_BOXES_XY
         yes=0
         r= calc_diff_mass_conservation(totsys_old(kkN,1,iXY),totsysn(iXY), &
         changePel(iXY)+changeBen(iXY),dt)
         call warning_test_mass_conservation(r,after,"N",yes)
         yesT=yes
         if (yes==2) then
           r= calc_diff_mass_conservation(totsys_old(kkN,2,iXY), &
           totbenn(iXY), jtBENPELn(iXY)+changeBen(iXY),dt)
           call warning_test_mass_conservation(r,after,"Ben N",yes)
           if ( yes==2) call OutputAfterConservationTest( &
                                  yes,totsys_old(kkN,2,iXY),r,iiBen,iiN)
         endif
         if (yesT==2.and.yes==0) then
           oldmass=totsys_old(kkN,1,iXY)-totsys_old(kkN,2,iXY)
           t=changePel(iXY)-jtBENPELn(iXY)
           r= calc_diff_mass_conservation(oldmass,totsysn(iXY)-totbenn(iXY),&
              t,dt)
           call warning_test_mass_conservation(r,after,"Pel N",yes)
           if (yes==2)call OutputAfterConservationTest(yes,oldmass,t,iiPel,iiN)
         endif
       enddo

       do iXY=1,NO_BOXES_XY
         yes=0
         r= calc_diff_mass_conservation(totsys_old(kkP,1,iXY),totsysp(iXY), &
                                                                 _ZERO_, dt)
         call warning_test_mass_conservation(r,after,"P",yes)
         yesT=yes
         if (yes==2) then
           r= calc_diff_mass_conservation(totsys_old(kkP,2,iXY), &
                                 totbenp(iXY), jtBENPELp(iXY),dt)
           call warning_test_mass_conservation(r,after,"Ben P",yes)
           if (yes==2) call OutputAfterConservationTest( &
                    yes,totsys_old(kkP,2,iXY),jtBENPELp(iXY),iiBen,iiP)
         endif
         if (yesT==2.and.yes==0) then
           oldmass=totsys_old(kkP,1,iXY)-totsys_old(kkP,2,iXY)
           r= calc_diff_mass_conservation( &
              oldmass,totsysp(iXY)-totbenp(iXY), -jtBENPELp(iXY),dt)
           call warning_test_mass_conservation(r,after,"Pel P",yes)
           if (yes==2)call OutputAfterConservationTest(yes,oldmass,r,iiPel,iiP)
         endif
       enddo

       do iXY=1,NO_BOXES_XY
         yes=0
         r= calc_diff_mass_conservation(totsys_old(kkSi,1,iXY),totsyss(iXY), &
                                                                  _ZERO_, dt)
         call warning_test_mass_conservation(r,after,"Si",yes)
         yesT=yes
         if (yes==2) then
           r= calc_diff_mass_conservation(totsys_old(kkSi,2,iXY), &
                                  totbens(iXY), jtBENPELs(iXY),dt)
           call warning_test_mass_conservation(r,after,"Ben Si",yes)
           if (yes==2) call OutputAfterConservationTest( &
                    yes,totsys_old(kkSi,2,iXY),jtBENPELs(iXY),iiBen,iiS)
         endif
         if (yesT==2.and.yes==0) then
           oldmass=totsys_old(kkSi,1,iXY)-totsys_old(kkSi,2,iXY)
           r= calc_diff_mass_conservation( &
              oldmass,totsyss(iXY)-totbens(iXY), -jtBENPELs(iXY),dt)
           call warning_test_mass_conservation(r,after,"Pel Si",yes)
           if (yes==2)call OutputAfterConservationTest(yes,oldmass,r,iiPel,iiS)
         endif
       enddo

!      do iXY=1,NO_BOXES_XY
!        r= calc_diff_mass_conservation( &
!            totsys_old(kkH,1,iXY),totsysh(iXY),_ZERO_, dt)
!        call warning_test_mass_conservation(r,after,"Alkalinity",yes)
!      enddo

       do iXY=1,NO_BOXES_XY
         totsys_old(kkN,1,iXY)=totsysn(iXY)
         totsys_old(kkP,1,iXY)=totsysp(iXY)
         totsys_old(kkSi,1,iXY)=totsyss(iXY)
         totsys_old(kkC,1,iXY)=totsysc(1)
         totsys_old(kkH,1,iXY)=totsysh(1)
       enddo
     endif
   endif
   end subroutine test_mass_conservation
   function calc_diff_mass_conservation(oldmass,newmass,rate,dt)
   use bio_var,only:secs_pr_day
   implicit none
   REALTYPE,intent(IN)             ::oldmass,newmass,rate,dt
   REALTYPE                        ::correct,t,rate_diff
   REALTYPE                        ::calc_diff_mass_conservation

   correct=newmass+rate*dt/secs_pr_day
   t=min(correct,oldmass)
   rate_diff= (correct-oldmass)/abs(t) * 100.0
   calc_diff_mass_conservation=rate_diff*secs_pr_day/dt
   end function calc_diff_mass_conservation

    subroutine warning_test_mass_conservation(value,after,nutrient,yes)
    use gotm_error_msg, only:set_warning_for_getm
    implicit none
    REALTYPE,intent(IN)            :: value
    character(len=*),intent(IN)    :: after
    character(len=*),intent(IN)    :: nutrient
    integer,intent(INOUT)          :: yes
    REALTYPE                       :: border,absvalue
    character(len=200)             :: msg
    character(len=1)               :: direction

!   if ((index('P',nutrient)==1)) STDERR "warning_test",value,yes
    yes=0
    absvalue=abs(value)
    if (absvalue< 0.01) return
    border=99.99
    if (absvalue>border) then
      direction='>';yes=2
    else
      border=0.05;yes=1
      if (absvalue> border) then
        yes=2
        direction='='
        write(msg,'(''No Mass conservation for '',A,'' after call to ''&
          &,A,'': '',A1,F8.4,''%'')') nutrient,after,direction,value
        STDERR trim(msg)
         call set_warning_for_getm
        return
      elseif (absvalue<0.02) then
        yes=0;border=0.02;direction='<'
      else
        yes=1;direction='<'
      endif
    endif
    write(msg,'(''No Mass conservation for '',A,'' after call to ''&
                 &,A,'': '',A1,F7.4,''%'')') nutrient,after,direction,border
    STDERR trim(msg)
    call set_warning_for_getm
    end subroutine warning_test_mass_conservation
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

   r =0.0D+00
   do i=1,nlev
      s= maxdepth*(1.0D+00 - tanh(ddu*float(nlev-i)/float(nlev))/tanh(ddu))
      arr(i)=(s+r) * 0.5D+00
      r=s
   enddo
   return
   end subroutine calc_sigma_depth
!EOC
!-----------------------------------------------------------------------
!!BOP
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
!-----------------------------------------------------------------------
!BOC

  subroutine output_warning_shiftings(warning_level,shiftcounter)
  use controlled_messages,only:controlled_output
  implicit none
  integer,intent(IN)              ::warning_level,shiftcounter

  character(len=80)               ::mess

  if (shiftcounter > 0) then
    if (warning_level >2) then
      write(mess,'(A,I3)')  &
        'number of masses shiftings between layers < 0.001%=',shiftcounter
    else
      write(mess,'(A,I3)')  &
        'number of shiftings of masses between vertical layers=',shiftcounter
     endif
   endif
   call controlled_output(shiftcounter,mess)
   end subroutine output_warning_shiftings

   end module bio_bfm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License
! www.gnu.org
!-----------------------------------------------------------------------
