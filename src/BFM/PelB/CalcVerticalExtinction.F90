#include "DEBUG.h"
#include "INCLUDE.h"

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcVerticalExtinction
!
! DESCRIPTION
!   Calculates the vertical extinction.
!     
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcVerticalExtinction(mode)
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6c
  ! The following box states are used (NOT in fluxes): PhytoPlankton
  ! The following Pelagic 1-d global boxvars are modified : xEPS
  ! The following Pelagic 1-d global boxvars  are used: ABIO_eps, ESS
  ! The following groupmember vars  are used: iiPhytoPlankton
  ! The following constituent constants  are used: iiC, iiL
  ! The following 0-d global parameters are used: p_eps0, &
  ! p_epsR6, p_epsESS, ChlLightFlag, p_epsChla
  ! The following 1-d global parameter vars are used: p_qchlc
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: R6c
! use mem, ONLY:P1l,P2l,P3l,P4l,P5l,P6l,P1c,P2c,P3c,P4c,P5c,P6c
#endif
  use mem, ONLY: xEPS, xEPS_0, xEPS_ESS, xEPS_Chl, ABIO_eps, ESS,  &
      iiPhytoPlankton, iiC, iiL, NO_BOXES,NO_BOXES_XY, PhytoPlankton
  use mem_Param, ONLY: p_eps0, p_epsR6, p_epsESS, ChlLightFlag, p_epsChla
  use mem_Phyto,ONLY:p_qchlc 
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
#ifdef INCLUDE_MACROPHYT
 use mem,ONLY: iiMacroContent,MacroContent
 use mem_Param, ONLY: CalcMacroPhyto 
 use constants,only:GET
 use mem_MacroPhyto,only:save_MacroPhyt_status=>save_status, &
        farm_surface,surface,recalc_ben_to_pel,layerMsm,Depth_layer_Msc
#endif

! !INPUT PARAMETERS:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer,intent(IN)   :: mode
!
! !AUTHORS
!   ERSEM-team
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i,iout
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
#ifdef INCLUDE_MACROPHYT
  real(RLEN), dimension(NO_BOXES)    ::r3
  real(RLEN), dimension(NO_BOXES_XY) ::mtm,psur
#endif

  select case (mode)
    case(0)
        xEPS_0(:)= ZERO
        xEPS_ESS(:)= ZERO
        xEPS_Chl(:)= p_epsR6* R6c(:)
    case (1)
      select case ( p_eps0  == ZERO)
        case( .TRUE. )
          xEPS_ESS(:)  =   ABIO_eps(:)
          xEPS_Chl(:)  =   p_epsR6* R6c(:)
        case( .FALSE. )
          xEPS_0(:)    =   p_eps0 
          xEPS_ESS(:)  =   p_epsESS* ESS(:)
          xEPS_Chl(:)  =   p_epsR6* R6c(:)
      end select
    case(2)
        xEPS_0(:)= ZERO
        xEPS_ESS(:)  = p_epsESS* ESS(:)
        xEPS_Chl(:)  = p_epsR6* R6c(:)
  end select
  call findnan(xEPS_Chl,NO_BOXES,iout)
  if ( iout.gt.0) then
     write(LOGUNIT,*)'xEPS_Chl=NaN p_eps0,mode,iout,R6c,ESS:',p_eps0,mode,iout,R6c(iout),ESS(iout)
  endif

  select case ( ChlLightFlag)
    case ( 1 )
      do i = 1 , ( iiPhytoPlankton)
        lcl_PhytoPlankton =>    PhytoPlankton(i,iiC)
        xEPS_Chl(:)  =   xEPS_Chl(:)+ p_epsChla * p_qchlc(i)* lcl_PhytoPlankton
      end do

    case ( 2 )
      do i = 1 , ( iiPhytoPlankton)
        lcl_PhytoPlankton =>    PhytoPlankton(i,iiL)
        xEPS_Chl(:)  =   xEPS_Chl(:)+ p_epsChla * lcl_PhytoPlankton
      end do
  end select
#ifdef INCLUDE_MACROPHYT
!write(LOGUNIT,*)'extinction: macrophyto'
      if ( CalcMacroPhyto) then
        do i=1,iiMacroContent
!write(LOGUNIT,*)"i,iiMacroContent",i,iiMacroContent
          where (save_MacroPhyt_status(i,:)==1) 
            mtm=Depth_layer_Msc(i,:) 
            psur=farm_surface(i,:)/surface
          endwhere
!write(LOGUNIT,*)'set pointer, iiL',iiL
          lcl_PhytoPlankton =>    MacroContent(i,iiL)
!write(LOGUNIT,*)'recalc_ben, GET',GET
!write(LOGUNIT,*)'NO_BOXES, NO_BOXES_XY',NO_BOXES,NO_BOXES_XY
!write(LOGUNIT,*)'shape',size(layerMsm(i,:))
!write(LOGUNIT,*)'shape2',size(lcl_PhytoPlankton)
!write(LOGUNIT,*)'shape3',size(mtm),size(psur),size(r3)
!write(LOGUNIT,*)'layerMsm',layerMsm(i,:)
!write(LOGUNIT,*)'lcl_Ph',lcl_PhytoPlankton
!write(LOGUNIT,*)'mtm',mtm
!write(LOGUNIT,*)'psur',psur
!write(LOGUNIT,*)'r3',r3
!JM          call recalc_ben_to_pel(GET,layerMsm(i,:),lcl_PhytoPlankton, &
!JM                                                     mtm,psur,r3)
          if (sum(save_MacroPhyt_status(i,:)).gt.0) &
            call recalc_ben_to_pel(GET,layerMsm(i,:),lcl_PhytoPlankton, &
                                                     mtm,psur,r3)
!         if (save_MacroPhyt_status(i,1)) &
!         write(LOGUNIT,*) "macroP:",r3(layerMsm(i,1)),lcl_phytoPlankton(1)
!write(LOGUNIT,*)'save_status'
          where (save_MacroPhyt_status(i,:)==1) &
                     xEPS_Chl(:)=xEPS_Chl(:)+r3* p_epsChla
         enddo
      endif
!write(LOGUNIT,*)'extinction: after macrophyto'
#endif
 
  xEPS= xEPS_0 + xEPS_ESS +xEPS_Chl

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
