#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   FixProportionCoeff0
!
! FILE
!   FixProportionCoeff.f90
!
! DESCRIPTION
!   
!	function FixProportionCoeff
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine FixProportionCoeff(NUTR,coeff1,coeff2,value1,value2)
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
        USE global_mem,      ONLY:RLEN,ZERO,NZERO,dummy
        use constants, ONLY: PARAMETER, INPUT_TERM, START_ADD_TERM, &
                             INPUT_ADD_TERM 
        use bennut_interface, ONLY: CompleteSet
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
!
!
        IMPLICIT  NONE
        INTEGER,intent(IN)        :: NUTR
        INTEGER,intent(IN)        :: coeff1
        INTEGER,intent(IN)        :: coeff2
        REAL(RLEN),intent(IN) ::value1
        REAL(RLEN),intent(IN) ::value2

        REAL(RLEN)            :: tmp
        select case ( abs(value2)+abs(value1)> 1.0D-20)
           case( .FALSE. )
              call CompleteSet( NUTR, INPUT_TERM, coeff2, PARAMETER, &
                                                     dummy, value=ZERO)
           case( .TRUE. )
             tmp=sign(abs(value2)+NZERO,value2)
             call CompleteSet( NUTR, START_ADD_TERM, coeff2, PARAMETER,&
                                               dummy, mfac=+1.0D+00/tmp)
             
             tmp=sign(abs(value1)+NZERO,value1)
             call CompleteSet( NUTR, INPUT_ADD_TERM, coeff1, PARAMETER,&
                                               dummy, mfac=-1.0D+00/tmp)
       end select
 
       return
      end

