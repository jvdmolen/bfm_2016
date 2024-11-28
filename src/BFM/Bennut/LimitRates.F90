#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: LimitRates
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
      module LimitRates
! !USES:
        USE global_mem,      ONLY:RLEN,ZERO,DONE,NZERO
        use mem_globalfun,   ONLY: insw


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  functions 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     contains


!EOP
!-------------------------------------------------------------------------!
!BOP

!
! FUNCTION
!   LimitShift
!
! DESCRIPTION
!   
!	function LimitShift
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine LimitShift(jK10K0x,K0x,K10x,max_shift,p)
!
! !AUTHORS
!   Piet Ruardij   
!
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        IMPLICIT  NONE
        REAL(RLEN),intent(INOUT)         ::jK10K0x ! Specification
        REAL(RLEN),intent(IN)            ::K0x ! Specification
        REAL(RLEN),intent(IN)            ::K10x ! Specification
        REAL(RLEN),intent(IN)            ::max_shift ! Specification
        REAL(RLEN),intent(OUT),optional  :: p ! Specification

        REAL(RLEN)       :: r,s 
 
        r= max(ZERO,insw(jK10K0x)*max(ZERO,K10x ) &
                                   +insw(-jK10K0x)* max(ZERO,K0x))
        if ( r.eq.ZERO ) then
          if (present(p)) then
             p=ZERO
          else
             jK10K0x=ZERO
          endif
        elseif (present(p)) then
          p=max_shift/(abs(jK10K0x/r)+max_shift);
        else
          jK10K0x=jK10K0x*max_shift/(abs(jK10K0x/r)+max_shift);
        endif
        return
        end subroutine  LimitShift
!EOP
!-------------------------------------------------------------------------!
!BOP

!
! FUNCTION
!   LimitChange.f90
!
! FILE
!   LimitChange.f90
!
! DESCRIPTION
!   
!	function LimitChange
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine LimitChange(mode,jK0x,K0x,max_change,p)
!
! !AUTHORS
!   Piet Ruardij   
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
        IMPLICIT  NONE
        INTEGER,intent(IN)        :: mode
        REAL(RLEN),intent(INOUT)  :: jK0x       ! Specification
        REAL(RLEN),intent(IN)     :: K0x        ! Specification
        REAL(RLEN),intent(IN)     :: max_change ! Specification
        REAL(RLEN),intent(OUT),optional  :: p ! Specification
 
        REAL(RLEN)                ::r
        LOGICAL                   ::ll
        ll=( mode ==1 .or. ( mode==2.and.jK0x <-NZERO) &
                                        .or.(mode==3.and.jK0x> NZERO)) 
        r=DONE
        if ( K0x <=ZERO) then
            if (ll) r=ZERO;
        elseif (ll) then
            r=max_change/(abs(jK0x)/K0x+max_change)
        endif
        if (present(p)) then
           p=r
        else
           jK0x=jK0x*r
        endif
        return
      end subroutine LimitChange
!EOP
!-------------------------------------------------------------------------!
!BOP
! FUNCTION
!   LimitChange_vector
!
!
! DESCRIPTION
!   
!	function LimitChange
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine LimitChange_vector(mode,jK0x,K0x,max_change,p)
!
! !AUTHORS
!   Piet Ruardij   
!
! CHANGE_LOG
!   
!
! COPYING
!   
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
        IMPLICIT  NONE
        INTEGER,intent(IN)             :: mode
        REAL(RLEN),intent(INOUT)       ::jK0x(:)    ! Specification
        REAL(RLEN),intent(IN)          ::K0x(:)     ! Specification
        REAL(RLEN),intent(IN)          ::max_change ! Specification
        REAL(RLEN),optional,intent(OUT)::p(:)       ! Specification

        REAl(RLEN)                :: r(size(jK0x(:)))
        logical                   :: ll(size(jK0x(:)))
        logical                   :: mm(size(jK0x(:)))
        mm=abs(jK0x)<NZERO
        select case  (mode)
          case(1); ll=.TRUE.
          case(2); ll=jK0x<-NZERO
          case(3); ll=jK0x> NZERO
        end select
        r=DONE
        where ( K0x <NZERO .or. mm ) 
          r=ZERO
        elsewhere (ll) 
          r=max_change/(abs(jK0x)/K0x+max_change)
        end where
        if (present(p)) then
           p=r
        else
           jK0x=jK0x*r;
        endif

        return
      end subroutine LimitChange_vector
!EOP
!-------------------------------------------------------------------------!
!BOP
! FUNCTION
!   LimitChange_vector
!
!
! DESCRIPTION
!   
!	function LimitChange
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine DoubleLimitChange_vector(mode,jK0x,K0x, &
                                 double_change,max_change,p)
!
! !AUTHORS
!   Piet Ruardij   
!
! CHANGE_LOG
!   
!
! COPYING
!   
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
        IMPLICIT  NONE
        INTEGER,intent(IN)             :: mode
        REAL(RLEN),intent(INOUT)       ::jK0x(:)    ! Specification
        REAL(RLEN),intent(IN)          ::K0x(:)     ! Specification
        REAL(RLEN),intent(IN)          ::double_change(:) ! Specification
        REAL(RLEN),intent(IN)          ::max_change ! Specification
        REAL(RLEN),optional,intent(OUT)::p(:)       ! Specification

        REAl(RLEN)                :: r(size(jK0x(:)))
        REAl(RLEN)                :: t(size(jK0x(:)))
        logical                   :: ll(size(jK0x(:)))
        logical                   :: mm(size(jK0x(:)))
        mm=abs(jK0x)<NZERO
        select case  (mode)
          case(1); ll=.TRUE.
          case(2); ll=jK0x<-NZERO
          case(3); ll=jK0x> NZERO
        end select
        r=DONE
        where ( K0x <NZERO .or. mm ) 
          r=ZERO
        elsewhere (ll) 
          t=double_change*max_change
          r=t/(abs(jK0x)/K0x+t)
        end where
        if (present(p)) then
           p=r
        else
           jK0x=jK0x*r;
        endif

        return
      end subroutine DoubleLimitChange_vector
!EOP
!-------------------------------------------------------------------------!
!BOP

!
! FUNCTION
!   LimitShift_m3
!
! DESCRIPTION
!   
!	function LimitShift
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine LimitShift_m3(jK10K0x,K0x,K10x,eK0m,eK10m,eshiftm,max_shiftt,p)
!
! !AUTHORS
!   Piet Ruardij   
!
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        IMPLICIT  NONE
        REAL(RLEN),intent(INOUT)         ::jK10K0x ! flux from K10->K)
        REAL(RLEN),intent(IN)            ::K0x ! mass in layer "0"  mmol/m2
        REAL(RLEN),intent(IN)            ::K10x !mass in layer "10"  mmol/m2
        REAL(RLEN),intent(IN)            ::eK0m ! transfer /m3-porewater ->/m2  in layer "0"
        REAL(RLEN),intent(IN)            ::eK10m ! transfer /m3-porewater ->/m2  in layer "0"
        REAL(RLEN),intent(IN)            ::eshiftm !  change of layers per timestep
        REAL(RLEN),intent(IN)            ::max_shiftt !  maxshift per time step
        REAL(RLEN),intent(OUT),optional  :: p ! Specification

        REAL(RLEN)       :: N 
        REAL(RLEN)       :: N0x
        REAL(RLEN)       :: N10x
        REAL(RLEN)       :: jN0x
        REAL(RLEN)       :: jN10x
        REAL(RLEN)       :: j
 
        N0x  =max(NZERO,K0x/eK0m)
        N10x =max(NZERO,K10x/eK10m)
        jN0x =(K0x+jK10K0x)/(eK0m+eshiftm) -K0x/eK0m
        jN10x=(K10x-jK10K0x)/(eK10m-eshiftm) -K10x/eK10m
        j    =insw(jK10K0x)
        N    = j*N10x +(1-j)* N0x
        j    = jN10x  +(1-j)*jN0x
        if (present(p)) then
          p=max_shiftt/(abs(j/N)+max_shiftt);
        else
          jK10K0x=jK10K0x*max_shiftt/(abs(j/N)+max_shiftt);
        endif
        return
        end subroutine  LimitShift_m3
!EOP
!------------------------------------------------------------------------

      end module LimitRates

