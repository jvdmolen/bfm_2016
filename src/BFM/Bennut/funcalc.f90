!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   funcalc.f90
!
! FILE
!   funcalc.f90
!
! DESCRIPTION
!   
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
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
      REAL(RLEN) FUNCTION funcalc(mode,chterm,coeff,basis,x)
        USE global_mem, ONLY:RLEN,LOGUNIT,DONE
        USE bennut_type
        USE bennut_variables,only:lowlevel_error
        USE constants
        USE bennut_interface,ONLY:BESSK1, BESSK0, BESSI1, BESSI0, &
          QGAUS_EXP
        use mem_globalfun,   ONLY: insw
        IMPLICIT  NONE
        integer,intent(IN) ::mode          ! Specification
        integer,intent(IN) ::chterm        ! Specification
        type (ty_coeff),intent(IN) ::coeff ! Specification
        REAL(RLEN),intent(IN) ::basis      ! Specification
        REAL(RLEN),intent(IN) ::x          ! Specification

        REAL(RLEN),parameter :: pi=3.141592D+00
        REAL(RLEN) ::r
        REAL(RLEN) ::rx
        REAL(RLEN) ::s
        integer ::i
        integer ::j
        integer ::term

        rx=x-basis
        term=coeff%ia
        if ( chterm .ne. 0 )term=max(term-2,min(term,0))
        select case (term)
        case (CONSTANT_TERM)
          select case (mode)
            case (INTEGRAL)
               funcalc=rx
            case (EQUATION,DERIVATIVE,SDERIVATIVE)
               funcalc=DBLE(max(1+mode,0))
            case default
               write(LOGUNIT,*) 'mode,term=',mode,term
               stop 'funcalc CONSTANT_TERM '
          end select
        case (LINEAR_TERM:)
          !integration of terms......
          select case (mode)
            case (INTEGRAL)
              funcalc=(rx**(term+1))/(term+1)
            case (EQUATION)
               funcalc=rx**(term)
            case (DERIVATIVE,SDERIVATIVE)
              !calculation for differention
              j=term+mode
              if (j.lt.0) then
                funcalc=ZERO
              else
                r=term
                if (mode == SDERIVATIVE .and. term > 2) r=(term-1)*r
                if (j > 0) r=r*(rx**(j))
                funcalc=r
              endif
            case default
               stop 'error funcalc FIRST_ORDER_TERM'
          end select
        case (EXPONENTIAL_TERM,ZERO_EXPONENTIAL_TERM)
          r=exp(coeff%labda(1)*rx)
          select case (mode)
            case (INTEGRAL)
              funcalc=r/coeff%labda(1)
            case (EQUATION)
              funcalc=r
            case (DERIVATIVE,SDERIVATIVE)
               funcalc=(coeff%labda(1)**(-mode))*r
            case (EXPONENTIAL_INTEGRAL)
              funcalc=ZERO
            case default
              stop 'funcalc EXPONENTIAL_TERM'
          end select 
        case (DIST_EXPONENTIAL_TERM)
          r=exp(coeff%labda(1)*rx)
          s=insw(coeff%labda(1))
          select case (mode)
            case (INTEGRAL)
              funcalc=r/coeff%labda(1)-rx*s
            case (EQUATION)
              funcalc=r-s
            case (DERIVATIVE,SDERIVATIVE)
               funcalc=(coeff%labda(1)**(-mode))*r
            case (EXPONENTIAL_INTEGRAL)
              funcalc=ZERO
            case default
              stop 'funcalc EXPONENTIAL_TERM'
          end select 
        case (BESSELI_EXP_TERM)
          r=exp(coeff%labda(1)*rx)
          s=r*coeff%labda(2)
          select case (mode)
            case (EXPONENTIAL_INTEGRAL)
              r=ZERO
              if (rx.ne.ZERO) then
                if (coeff%labda(1).le.ZERO) then
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                     2.0D+00*coeff%labda(1),bessi0)
                else
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                     2.0D+00*coeff%labda(1),bessi1)
                endif
              endif
              funcalc=r;
            case (INTEGRAL)
              r=ZERO
              if (rx.ne.ZERO) then
                if ( coeff%labda(1).lt.ZERO) then
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                                      ZERO,bessi0)
                else
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                                      ZERO,bessi1)
                endif
              endif
              funcalc=r
            case (EQUATION)
              if ( coeff%labda(1).lt.ZERO) then
                funcalc=bessi0(s) 
              else
                funcalc=bessi1(s) 
              endif
            case (DERIVATIVE)
              if ( coeff%labda(1).lt.ZERO) then
                funcalc=bessi1(s)*coeff%labda(1)*s
              else
                funcalc=bessi0(s)*coeff%labda(1)**2*s*rx
              endif
            case default
              stop 'funcalc BESSELI_EXP_TERM'
            end select 
        case (BESSELK_EXP_TERM)
          r=exp(coeff%labda(1)*rx)
          s=r*coeff%labda(2)
          select case (mode)
            case (EXPONENTIAL_INTEGRAL)
              r=ZERO
              if (rx.ne.ZERO) then
                if ( coeff%labda(1).lt.ZERO) then
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                     2.0D+00*coeff%labda(1),bessk0)
                else
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                     2.0D+00*coeff%labda(1),bessk1)
                endif
              endif
              funcalc=r;
            case (INTEGRAL)
              r=ZERO
              if (rx.ne.ZERO) then
               if ( coeff%labda(1).lt.ZERO) then
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                                      ZERO,bessk0)
                else
                  r=qgaus_exp(ZERO,rx,coeff%labda(1),coeff%labda(2), &
                                                      ZERO,bessk1)
                endif
              endif
              funcalc=r
            case (EQUATION)
              if ( coeff%labda(1).lt.ZERO) then
                funcalc=bessk0(s) 
              else
                funcalc=bessk1(s) 
              endif
            case (DERIVATIVE)
              if ( coeff%labda(1).lt.ZERO) then
                funcalc=-bessk1(s)*coeff%labda(1)*s
              else
                funcalc=bessk0(s)*coeff%labda(1)**2*s*rx
              endif
            case default
              stop 'funcalc BESSELK_EXP_TERM'
          end select 
        case default
           write(LOGUNIT,*) 'term=',term,'funcalc no such term'
           lowlevel_error=1
           funcalc=DONE
        end select 
        return
      end
