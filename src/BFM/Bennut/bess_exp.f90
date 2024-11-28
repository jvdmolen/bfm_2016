!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   bess_exp.f90
!
! FILE
!   bess_exp.f90
!
! DESCRIPTION
!   FILE: Solve_coupled.F  nutrient dynamics in the benthos 
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
      REAL(RLEN) FUNCTION BESS_EXP(X,LAMBDA_1,LAMBDA_2,LAMBDA_3,FN)
        USE global_mem, ONLY:RLEN
        USE constants
        IMPLICIT  NONE
        REAL(RLEN),intent(IN)    ::x ! Specification
        REAL(RLEN),intent(IN)    ::LAMBDA_1 ! Specification
        REAL(RLEN),intent(IN)    ::LAMBDA_2 ! Specification
        REAL(RLEN),intent(IN)    ::LAMBDA_3 ! Specification

        INTERFACE                             ! Specification
          REAL(RLEN) FUNCTION FN(X)          ! Specification
          USE global_mem, ONLY:RLEN
          REAL(RLEN),INTENT(IN)   ::X              ! Specification
        END FUNCTION                         ! Specification
        END  INTERFACE                       ! Specification     

        real(RLEN)  :: r

        r=LAMBDA_2*exp(X*LAMBDA_1)
        if (LAMBDA_3.eq.0.0) then
           BESS_EXP=FN(r)
        else
           BESS_EXP=exp(X*LAMBDA_3)* FN(r);
        endif
        RETURN
      end

