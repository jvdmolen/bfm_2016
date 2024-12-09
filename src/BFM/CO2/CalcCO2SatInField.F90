#include "DEBUG.h"


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcCO2SatInField.f90
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcCO2SatInField(k,n,ERHO,ETW,ESW,cc)

!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem, ONLY: ppN1p,ppN5s,ppO3h,ppO3c
  use CO2System,ONLY: CalcCO2System,HPlusBasis
  use mem_CO2
  use constants,ONLY:MW_C

!
!
! !AUTHORS
!   H. Thomas   and  P. Ruardij
!
! !REVISION_HISTORY
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij and M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation
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
  integer                    :: k
  integer                    :: n
  real(RLEN),intent(IN)      :: ERHO(1:k)
  real(RLEN),intent(IN)      :: ETW(1:k)
  real(RLEN),intent(IN)      :: ESW(1:k)
!JM  real(RLEN),intent(INOUT)   :: cc(1:n,1:k)
  real(RLEN),intent(INOUT)   :: cc(1:k,1:n)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,external    ::D3toD1
  integer             ::error
  integer             ::i
  integer             ::j
  integer             ::m
  real(RLEN)          ::Ac,Ac_max
  real(RLEN)          ::DIC
  real(RLEN)          ::dumCO2
  real(RLEN)          ::dumCO3
  real(RLEN)          ::dumHCO3
  real(RLEN)          ::dumCac
  real(RLEN)          ::dumpH

  m=ppO3c
  j=ppO3h
  if ( k > 0 ) then
    do i=1,k
!JM          Ac   =   cc(j,i)
          Ac   =   cc(i,j)
          Ac_max=   (Ac+HplusBASIS)*0.5
          Ac=min(Ac_max,Ac)
!JM          error= CalcCO2System(MethodCalcCO2,ESW(i),ETW(i),ERHO(i),&
!JM                   cc(ppN1p,i),cc(ppN5s,i),Ac,&
!JM                   dumCO2,dumHCO3,dumCO3,dumCac,&
!JM                   pCO2_in=pCO2_Air,DIC_out=DIC,pH_out=dumpH)
          error= CalcCO2System(MethodCalcCO2,ESW(i),ETW(i),ERHO(i),&
                   cc(i,ppN1p),cc(i,ppN5s),Ac,&
                   dumCO2,dumHCO3,dumCO3,dumCac,&
                   pCO2_in=pCO2_Air,DIC_out=DIC,pH_out=dumpH)
          if ( error.eq.0) then
!JM            cc(m,i)=DIC* MW_C
!JM            cc(j,i)=Ac
            cc(i,m)=DIC* MW_C
            cc(i,j)=Ac
          else
            write(LOGUNIT,*) 'error,DIC',error,DIC
            write(LOGUNIT,*) 'i=',i
            write(LOGUNIT,*) 'AC=',Ac
!JM            write(LOGUNIT,*) 'N1p=',cc(ppN1p,i)
!JM            write(LOGUNIT,*) 'N5s=',cc(ppN5s,i)
            write(LOGUNIT,*) 'N1p=',cc(i,ppN1p)
            write(LOGUNIT,*) 'N5s=',cc(i,ppN5s)
            write(LOGUNIT,*) 'ESW=',ESW(i)
            write(LOGUNIT,*) 'ETW=',ETW(i)
            write(LOGUNIT,*) 'ERHO=',ERHO(i)
            write(LOGUNIT,*) 'pCO2=',pCO2_Air
          endif
    enddo
  endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
