#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   RecalcPenetrationDepth.f90
!
! FILE
!   RecalcPenetrationDepth.f90
!
! DESCRIPTION
!
! In subroutine RecalcPenetrationDepth the penetration depth of detritus is
! recalculated when detritius is sedimentated on top of the sediment in
! such a way that the amount of detritus in the lower layers keep the same mass.
!       The new thickness is calculated Newton's rule of approximation:
!           x(i+1)=x(i)-f(x)/f'(x for f(x)==0
!
! !INTERFACE
    subroutine RecalcPenetrationDepth(D1m,Dtm, Dxm,Dhm, input, mass,newDxm )
!
! !AUTHORS
!   Piet Ruardij
!
! !USES:
     use global_mem, ONLY:RLEN,NZERO,DONE,LOGUNIT,ZERO
     use mem,ONLY:NO_BOXES_XY
     use mem_Param,  ONLY: p_clDxm
     use mem_globalfun, ONLY:exp_limit_scalar,IntegralExpDist, &
           SetExpDist,insw
     use constants,only:MAX_VAL_EXPFUN
     use BFM_ERROR_MSG,only:set_warning_for_getm


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

        REAL(RLEN),intent(IN)     ::D1m
        REAL(RLEN),intent(IN)     ::Dtm
        REAL(RLEN),intent(IN)     ::Dxm
        REAL(RLEN),intent(IN)     ::Dhm
        REAL(RLEN),intent(IN)     ::input
        REAL(RLEN),intent(IN)     :: mass
        REAL(RLEN),intent(OUT)    :: newDxm

        REAL(RLEN),parameter      ::p_cm_Dxm=1000.0,p_ml_Dxm=100.0,p_layer=0.3
!       REAL(RLEN),parameter      ::p_clDxm=0.015
        REAL(RLEN)                ::alpha
        REAL(RLEN)                ::old
        REAL(RLEN)                ::fx
        REAL(RLEN)                ::dfx
        REAL(RLEN)                ::c,mass0,massinf
        REAL(RLEN)                ::newalpha,lone,lsig

        alpha=SetExpDist(Dxm,Dlm=p_clDxm,Dmm=p_cm_Dxm)
!       alpha=abs(alpha)
        if (alpha*input< ZERO .and.abs(Dxm).gt.p_ml_dxm) alpha=-alpha;
        ! calculated mass present below D1m
!       old = mass * exp(-alpha *D1m)/alpha
        mass0= mass /IntegralExpDist(-alpha,D1m)
        massinf=mass0*IntegralExpDist(-alpha,p_layer)
        newalpha= alpha

        !alpha>0 : lone=0 ,lsig=-1
        lone=DONE-insw(alpha)
        lsig=-DONE+2.0D+00*lone
        old = massinf * (exp(-alpha *D1m) -lone)

        !start iteration with old_alpha
        newalpha= alpha
        if ( abs(input) > NZERO) then
          c=DONE
          do while ( c > 1.D-5)
            !calc mass below D1m with new input
            fx = (massinf+input) * ( exp(-newalpha *D1m) -lone)
            !determine first derivative
            dfx  =  ( lsig * D1m  ) *(massinf+input) * exp(-newalpha *D1m) 
            !keep old value
            c =newalpha
            ! calc new alpha for (fx-old)==0
            newalpha = newalpha - (fx-old) /dfx
            ! use c to calculate iteration precision
            c = abs(c -newalpha)/c
          enddo
       else
         newDxm=alpha;return
       endif

       !limit maximal step
       dfx=c
       c=DONE/max(0.05D+00*(-lsig)*alpha,-lsig*newalpha)
       c=sign(min(max(c,p_clDxm),p_cm_Dxm),alpha)
        newDxm=c
        return
        end
