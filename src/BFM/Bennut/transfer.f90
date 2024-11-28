!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   transfer.f90
!
! FILE
!   transfer.f90
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
      REAL(RLEN) FUNCTION transfer(mode,coeff,input,diff)
        USE global_mem, ONLY:LOGUNIT,ZERO
        USE bennut_constants
        USE bennut_type
        USE constants
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        type (ty_coeff),intent(IN) ::coeff ! Specification
        REAL(RLEN),intent(IN) ::diff ! Specification
        REAL(RLEN),intent(IN) ::input ! Specification

        REAL(RLEN) ::r
        integer ::term 

        term=coeff%ia
        if (input.eq.ZERO) then
          transfer=ZERO
          return
        endif

        select case (term)
          case (ZERO_EXPONENTIAL_TERM,DIST_EXPONENTIAL_TERM)
            r=-coeff%labda(1)**2*diff
            if (mode.eq.PARA2COEFF) then
              transfer=input/r
            else
              transfer=input*r
            endif
          case (QUADRATIC_TERM:)
            r=-real(term)*diff
!           r=-2.D+00*diff
            if (mode.eq.PARA2COEFF) then
              transfer=input/r
            else
              transfer=input*r
            endif
          case ( CONSTANT_TERM,LINEAR_TERM) 
            transfer=input
          case (BESSELI_EXP_TERM, BESSELK_EXP_TERM)
            if (mode.eq.PARA2COEFF) then
              transfer=sqrt(input/diff)/abs(coeff%labda(1))/2.0D+00
            elseif (mode.eq.LABDA_2) then
              r=coeff%labda(1) 
              transfer = 2.0 / abs( r ) * sqrt(coeff%labda(2) / diff);
            elseif (mode.eq.LABDA_1) then
              transfer = coeff%labda(1) / 2.0D+00 
            else
              write (LOGUNIT,*) 'mode=',mode
              write (LOGUNIT,*) 'term=',term
              stop 'stop transfer'
            endif
          case default
            write (LOGUNIT,*) 'mode=',mode
            write (LOGUNIT,*) 'term=',term
            stop 'stop ransfer '
        end select

        return
      end

