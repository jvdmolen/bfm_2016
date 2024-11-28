#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalculateThermoVars
!
! DESCRIPTION
!   !	This submodel calls all other submodels
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine CalculateThermoVarsBFM
!
! !USES:
  ! The following 0-d global parameters are used: CalcPelagicFlag, &
  ! CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO
  use mem,  ONLY: ETW, Depth, SdTdzTH, TMLd, BMLd
  use mem, ONLY: BoxNumberXY, NO_BOXES_Z,NO_BOXES_XY,PelBoxAbove
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:ResetTotMassVar
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: ResetTotMassVar
  use mem_Param, ONLY: p_mdTdz

!  
!
! !AUTHORS
!   Piet Ruardij
!
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer                            :: f,t
  integer                            :: iiTMLd, iiBMLd,iimdTdz

  DO BoxNumberXY=1,NO_BOXES_XY
    f=PelBoxAbove(BoxNumberXY)
    t=f+NO_BOXES_Z-1
    call CalculateThermovar(NO_BOXES_Z,p_mdTdz,Depth(f:t),ETW(f:t), &
            iiBMLd,BMLd(BoxNumberXY),iiTMLd,TMLd(BoxNumberXY),SdTdzTh(BoxNumberXY))
  enddo
  end

  subroutine CalculateThermovar(n,mdTdz,z,T,iiBMLd,BMLd,iiTMLd,TMLd,SdTdzTh)
    use global_mem, ONLY:RLEN,LOGUNIT,ZERO
  IMPLICIT NONE
  integer,intent(IN)                 :: n
  real(RLEN),intent(IN),dimension(1:n)  :: z
  real(RLEN),intent(IN),dimension(1:n)  :: T
  real(RLEN),intent(IN)              :: mdTdz
  integer,intent(OUT)                ::iiBMLD,iiTMLd 
  real(RLEN),intent(OUT)             ::BMLD,TMLd,SdTdzTh 
  integer                            :: i,jf,jt,iimdtdz
  real(RLEN),dimension(1:n)          :: dTdz, dT2dz2,dz
  integer,dimension(1)               :: jjmaxdT2Dz2,jjmaxdTdz, jjmindT2dz2

    jf=1
    jt=n
    dz(jf:jt-1)=0.5*(z(jf:jt-1)+z(1+jf:jt))
    dTdz(jf:jt-1)=(T(1+jf:jt)-T(jf:jt-1))/ dz(jf:jt-1)
    iiTMLd=0
    if ( maxval((dTdz(jf:jt-1))) >mdTdz)  then
      ! Find interface beween two ;layer with the highest dTdZ
      jjmaxdTdz=maxloc((dTdz(jf:jt-1))); iimdTdz=jjmaxdTdz(1)+jf-1;
      ! Calculated dt2/dz2
      dz(jf:jt-2)=0.5*(dz(jf:jt-2)+dz(1+jf:jt-1))
      dT2dz2(jf:jt-2)=(dTdz(1+jf:jt-1)-dTdz(jf:jt-2))/dz 
      ! Determine in lower part of the gradient where the maximum is:
      ! depth of the the top of the bottom mixed layer
      jjmaxdT2Dz2=maxloc((dT2dz2(jf:iimdTdz)));iiBMLd=jjmaxdT2Dz2(1)+jf-1;
      ! Determine in upper part of the gradient where the minimum is:
      ! depth of the top/surface mixed layer
      jjmindT2dz2=minloc((dT2dz2(iimdTdz:jt-1))); iiTMLd=jjmindT2dz2(1)+jjmaxdTdz(1)-1
      if (iiTMLd>1) then
        TMLd=-sum(z(iiTMLd+1:n))
        BMLd=TMLd-sum(z(iiBMLd+1:iiTMLd))
        SdTdzTh=sum(dTdz(iiBMLd:iiTMLd))
      endif
    endif
    if (iiTMLd<=1) then
      SdTdzTh=ZERO;iiBMLD=0;iiTMLd=0
      BMLd=ZERO;TMLd=ZERO
    endif
  end

