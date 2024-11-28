#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Alkalinity
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine AlkalinityDynamics(control)
!
! !USES:

#ifdef INCLUDE_PELCO2
use global_mem, ONLY:RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#endif
  use mem, ONLY: ppO3h, ppN6r,ppN3n,ppB1n,NO_BOXES, iiPel,iiN,  &
      ppPhytoPlankton,iiPhytoPlankton, &
      flN3O4n,flN4N3n, flux_vector, Depth,iiTotal,iiConsumption
  use mem_CO2, ONLY: p_qhK4K3n,p_qhATo
  use mem_param,  ONLY: p_qro, p_qon_dentri 
  use SourceFunctions,ONLY:Source_D3_withstate,Source_D3_withgroup 
  use botflux,only:getbotflux_3D

!
! !AUTHORS
!   16 March 1999 Original version by H. Thomas
!
!
!
! !REVISION_HISTORY
!   
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
  IMPLICIT NONE

  integer,intent(IN)                      :: control

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)          :: rN3n,rN6r
  integer                                 :: iout
  !
  ! correction of the alkalinity
  !
  !Use only rates defined for nitrate:
  rN3n= flN4N3n - Source_D3_withgroup( &
      ppN3n,ppPhytoPlankton,iiPhytoPlankton,iiN,iiConsumption) &
      -Source_D3_withstate(ppN3n,ppB1n,iiConsumption)   

  !cycle for N6r is not defined : correct total flux due production and consumption
  ! for input from the sediment
  rN6r=Source_D3_withstate(ppN6r,ppN6r,iiTotal)  
  call flux_vector( iiPel, ppO3h,ppO3h,  - (rN3n-flN3O4n)*p_qhK4K3n & 
                                        + (rN6r/p_qro-flN3O4n/p_qon_dentri)*p_qhATo)
#endif
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
