!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitBoxParams
!
! DESCRIPTION
!   INitialization of for Box1DParameters

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine InitBoxParams
!
! USES:
  use global_mem
  use mem, ONLY: NO_BOXES, NO_BOXES_XY, ppK1p,ppK5s
  use mem_Param, ONLY:  p_pK1_ae,p_pK5_ae,p_poro,p_1d_poro

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   ---
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer:: i
  integer:: status
! real(RLEN) :: kd
! real(RLEN) :: rho=2.65

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ! Set Values for 1D-parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (.NOT. allocated(p_pK1_ae) ) allocate(p_pK1_ae(1:NO_BOXES_XY),stat=status)
  if (.NOT. allocated(p_pK5_ae) ) allocate(p_pK5_ae(1:NO_BOXES_XY),stat=status)

  if (.NOT. allocated(p_poro) ) allocate(p_poro(1:NO_BOXES_XY),stat=status)
  do i=1,NO_BOXES_XY
    p_poro(i)= p_1d_poro
    call  CalcAdsorptionFromPoro(ppK1p,p_poro(i),p_pK1_ae(i))
    call  CalcAdsorptionFromPoro(ppK5s,p_poro(i),p_pK5_ae(i))
!   kd = 4.03087  * (p_poro(i) - 0.38662 )/ 0.00415 
!   p_pK1_ae(i) = kd*(1.0-p_poro(i))/p_poro(i)*rho
  enddo

  end subroutine

  subroutine CalcAdsorptionFromPoro(mode,poro,output )
! USES:
  use global_mem,ONLY:RLEN
  use mem,Only: ppK1p,ppK5s

  implicit none
  integer               :: mode
  real(RLEN),INTENT(IN) :: poro
  real(RLEN),INTENT(OUT):: output

  real(RLEN),parameter  :: rho=2.65
  real(RLEN)            :: kd
  
  select case (mode)
    case (ppK1p)
      kd = 4.03087  * (poro - 0.38662 )/ 0.00415 
      output = kd*(1.0-poro)/poro*rho
    case (ppK5s)
      kd = 4.03087  * (poro - 0.38662 )/ 0.2856 
      output = kd*(1.0-poro)/poro*rho
  end select
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
