#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenOxygen(mode)
!
! DESCRIPTION
!   Description of first order oxic processes in the sediment
!       and computation of oxygen penetration depth
!
!
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenOxygenDynamics(mode)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D1m, G2o, D2m
  ! The following global scalar vars are used: InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : shiftD1m
  ! The following Benthic 1-d global boxvars are used: ETW_Ben, irrenh, &
  ! rrBTo, jG2K3o, jG2K7o, O2o_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, &
  ! CalcBenthicFlag

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,dummy,ZERO_KELVIN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: D1m, G2o, D2m, irri_bio
#endif
  use mem, ONLY: ppD1m, ppO2o,ppG2o, ppD2m, iiPel,iiBen,KNO3,M3n, &
   shiftD1m,Wind,Depth_Ben,ETW_Ben,O2o_Ben,  &
    jcrrBTo,rrBTo,reBTo, jG2K3o, jG2K7o,jO2Y2o,irrenh,dry_z, G2_xavail_o,  & 
    LocalDelta,max_change_per_step, InitializeModel,  &
    NO_BOXES_XY, BoxNumberXY,flux_vector
  use constants,  ONLY: SEC_PER_DAY, ONE_PER_DAY, BENTHIC_BIO,STANDARD,EQUATION
  use constants,  ONLY: ANY
  use mem_Param,  ONLY: p_poro, p_d_tot, CalcBenthicFlag
  use botflux,onlY:addbotflux_vector
  use mem_BenOxygen

  use bennut_interface, ONLY: CalculateFromSet
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m,p_dry_ben,p_clDxm
  use gotm_error_msg,only:set_warning_for_getm
  use LimitRates, ONLY:LimitChange_vector
  use Constants, ONLY:NEGATIVE



! use mem,ONLY:Output2d_1,Output2d_2,Output2d_3, Output2d_4


  IMPLICIT NONE
  integer,intent(IN)          :: mode
!  
!
! !AUTHORS
!   P. Ruardij
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M.Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: rh
  real(RLEN),dimension(NO_BOXES_XY)  :: diff
  real(RLEN),dimension(NO_BOXES_XY)  :: zmG2o,reO2o,remO2o
  real(RLEN),dimension(NO_BOXES_XY)  :: D1mNew,D1m_max
  real(RLEN),dimension(NO_BOXES_XY)  :: G2oNew
  real(RLEN),dimension(NO_BOXES_XY)  :: jG2O2o

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Emperical equation derived from Broecker and Peng (1973)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call SurfaceO2Diffusion(0,ETW_Ben(1),Wind,Depth_Ben(1), &
                      rh(1),diff(1))

    where ( dry_z <=p_clDxm .and. p_dry_ben) 
      diff=diff/p_poro*irrenh
    elsewhere
      diff = SEC_PER_DAY* 1.0D-9* (10.0D+00)**((- 984.26D+00/( -ZERO_KELVIN+ &
      ETW_Ben(:))+ 3.672D+00))*irrenh* p_exsaf
    endwhere

  if ( mode.eq.1) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Determine maximum value of D1m after 1 day when limiting the enhancement 
    ! to p_sumD1m per day
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     D1m_max=D1m+p_sumD1m*LocalDelta
     remO2o =   2.0D+00* diff* p_poro* O2o_Ben(:)/  &
                               (D1m_max)**2 *(D1m_max*p_poro)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Recalculate total consumption from /m2 to /m3 pw:
    ! In case of of "netative respirations' (due to production of the benthic 
    ! diatoms) the consequences are kep in limits by limiting enhanvement of D1m 
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    reO2o  =  ( rrBTo(:)-reBTo(:)+ jG2K3o(:)+ jG2K7o(:))
    zmG2o  =  max(reO2o,remO2o) /( D1m(:)* p_poro)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Determine new thickness:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    D1mNew  =   sqrt(  2.0D+00* diff* p_poro* O2o_Ben(:)/( NZERO+ zmG2o))
    D1mNew  = min(D1mNew ,p_d_tot-2.0 * p_clD1D2m);

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate rate of change of thickness of the aerobic layer:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    shiftD1m(:)  =  ( max(  p_mD1m,  D1mNew)- D1m(:))/ ONE_PER_DAY

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Damping the change of D1m in case of large changes
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if ( InitializeModel== 0) then
!   shiftD1m(:) = shiftD1m(:)* (D1m(:)/( D1m(:)+ &
!     abs(shiftD1m(:))))**(p_xdampingD1m)*( p_chD1m/( p_chD1m+ D1m(:)))
      shiftD1m(:) = shiftD1m(:)*(p_chD1m/(p_chD1m+ D1m(:)))
      if (D1m(1).lt.p_mD1m) write(LOGUNIT,*)'BOI shiftD1m',D1m,shiftD1m
 
      do BoxNumberXY=1,NO_BOXES_XY
       rh(BoxNumberXY) =max(ZERO, CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
                       STANDARD, D1m(BoxNumberXY), dummy)/(NZERO+M3n(BoxNumberXY)))
 !     if ( rh(BoxNumberXY) < ZERO.and.M3n(BoxNumberXY)>NZERO) then
 !       write(LOGUNIT,*) "BenOxygen proportion M3n(D1m)/M3n(0..D2m)=", &
 !                 rh(BoxNumberXY),' ShiftD1m=',shiftD1m,' M3n=',M3n(BoxNumberXY)
 !       call set_warning_for_getm
 !     endif
      enddo

      ! At D1m where above production of NO3 (nitrification) take place and 
      ! below NO3 consumption ( denitrification) represent the point 
      ! with (most of time) the highest concentration. 
      ! if M3n at D1m is smalller than the average or is even negative there 
      ! is no production of NO3 (=nitrification) and we have to deal with 
      ! nearly anoxic sutuations. In these cases we limit the changes of D1m.
      where (shiftD1m(:)> ZERO) &
      shiftD1m(:)= shiftD1m(:) * max(-DONE,min(DONE,rh*2.0D+00));
    end if
      if (D1m(1).lt.p_mD1m) write(LOGUNIT,*)'BOII shiftD1m,rh',D1m,shiftD1m,rh


    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Damping the change of D1m in case of too thick D1m
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    rh= min( shiftD1m(:),max(ZERO,p_d_tot-p_chD1m-D1m(:)))
    if (D1m(1).lt.p_mD1m) write(LOGUNIT,*)'BOIII shiftD1m',D1m,rh
  
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! recalculate the new D1mNew at the actual time step:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    D1mNew  =   D1m(:)+ rh* LocalDelta

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! calculate the consumption which belongs to the corrected D1mNew
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    zmG2o  =  ( 2.0D+00* diff* p_poro* O2o_Ben(:))/( D1mNew* D1mNew)
  
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! New oxygen conc. in the sediment:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    G2oNew = D1mNew*( O2o_Ben(:)- 0.66667D+00* zmG2o* D1mNew* D1mNew/( 2.0D+00* &
      diff* p_poro))* p_poro

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! flux to pelagic: correct flux for rate of change of G2o
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    jG2O2o = -reO2o-( G2oNew- G2o(:))/  LocalDelta

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !  Assign fluxes
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    if ( InitializeModel== 0) then
      shiftD1m(:)=rh
      call flux_vector( iiBen, ppD1m,ppD1m, rh )
      call addbotflux_vector(ANY,iiBen,ppG2o,iiPel,ppO2o,jG2o2o)

      if ( CalcBenthicFlag== BENTHIC_BIO) then
        ! Compute shifting of the denitrification layer here in case of running 
        ! only the benthic submodel and NOT the benthic nutrient model.
        rh  =   p_d_tot- D2m(:)
        call flux_vector( iiBen, ppD2m,ppD2m, shiftD1m(:)* rh/( rh+ 0.01D+00) )
      end if
    else
        G2o(:)=G2oNew
    endif
  elseif ( mode.eq.0 ) then
    ! this call is before the biology is caluclated
    ! because of irrenh depnds on activity of the biology we use
    ! here the value of the previous step  which is keep in 2d-state irri_bio
    where ( dry_z <=p_clDxm .and. p_dry_ben) 
      diff=diff/p_poro*irri_bio
    elsewhere
      diff = SEC_PER_DAY* 1.0D-9* (10.0D+00)**((- 984.26D+00/( -ZERO_KELVIN+ &
      ETW_Ben(:))+ 3.672D+00))*irri_bio* p_exsaf
    endwhere
    !Determine max.allowed change in D1m
    rh=-10.0D+00
    call LimitChange_vector(NEGATIVE,rh,D1m,max_change_per_step)
    !Let change D1m such that the local consmption becomes high as possible..
    rh=max(p_mD1m,D1m +rh*LocalDelta)
    zmG2o  =  ( 2.0D+00* diff* p_poro* O2o_Ben(:))/( rh* rh)
    !Calculate G2_xavail_o for use in in the call to LimitChange'
    G2_xavail_o=G2o+zmG2o*rh*p_poro*LocalDelta
  elseif(mode.eq.2) then
    !Calculate for actual D1m and O2_Ben mineralization
    zmG2o  =  ( 2.0D+00* diff* p_poro* O2o_Ben(:))/( D1m* D1m)
    !Calculate what mineralization should be if all oxygen should be taken of
    ! oxiclayer instead of form the oxygen in the adjacent layer.
    rh=zmG2o+jO2Y2o/(p_poro *D1m)
!   write(LOGUNIT,*) 'rh=',rh
    rh =rh*(D1m*D1m)/(2.00D+00*p_poro*O2o_Ben) 
    !Calculate irrenh on basis 
    irrenh=irrenh* max(DONE,rh/diff)
  endif
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
