
#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelForcingForBen
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine FilterLimPart(BoxNumber,efilPART)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
  use mem, ONLY: ESS, iiPel,iiPhytoPlankton, ppPhytoPlankton, PhytoPlankton 
  use mem, ONLY: NO_BOXES, NO_BOXES_X, NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY,iiP1,iiP5,iiC
  use mem_Param,  ONLY: CalcPhytoPlankton
!  
!
! !AUTHORS
!   Piet Ruardij
!
!
! !REVISION_HISTORY
!   Created at Wed Jun 16 02:04:44 PM CEST 2004
!
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
!------------------------------------------------------------------ 
!BOC
!
  IMPLICIT NONE
  integer,intent(in)                                   :: BoxNumber
  real(RLEN),intent(OUT)                               :: efilPART


  real(RLEN), dimension(:), pointer                ::lcl_Plankton
  real(RLEN)                                       ::part,phyto
  real(RLEN)                                       ::parphyt=8.0e6
  real(RLEN)                                       ::parsil=2.38e+5
  real(RLEN)                                       ::OCP=2.0E10

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


    phyto=ZERO
    lcl_Plankton=>PhytoPlankton(iiP1,iiC);   
    phyto=phyto+lcl_Plankton(BoxNumber)
    lcl_Plankton=>PhytoPlankton(iiP5,iiC);   
    phyto=phyto+lcl_Plankton(BoxNumber)
    part=parphyt*phyto +ESS(BoxNumber)*parsil
    efilPART=2.0*OCP/(2.0*OCP+part)
end
