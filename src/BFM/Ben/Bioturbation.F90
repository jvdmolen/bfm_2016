#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Bioturbation
!
! DESCRIPTION
!   Describes sediment reworking by benthic organisms and their effect on
!   vertical transport of dissolved (bioirrigation) and particulate 
!   (bioturbation) matter
! 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BioturbationDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D6m, D7m, D8m, D9m
  ! The following Benthic-states are used (NOT in fluxes): Y2c, Y5c, Y1c, Y4c
  ! The following Benthic 1-d global boxvars are modified : turenh
  ! The following Benthic 1-d global boxvars got a value: irrenh, &
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,DONE
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: D1m,D6m, D7m, D8m, D9m, Y2c, Y5c, Y1c, Y3c, Y4c
#endif
  use mem, ONLY: ppD6m, ppD7m, ppD8m, ppD9m,iiY1,iiY2,iiY4,iiY5, &
    max_change_per_step, LocalDelta, irrenh, turenh,pxturinD1, ETW_Ben, &
    ppirri_bio,irri_bio,NO_BOXES_XY, iiBen, flux_vector,jnetYIc,jnetY3c
  use mem_Bioturbation
  use constants,only:NEGATIVE
  use LimitRates, ONLY:LimitChange_vector
  use mem_Param,  ONLY: p_clDxm
! use mem_BenOrganism, ONLY:p_cm
! use mem,ONLY: LocalDelta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: MM_vector,insw_vector,SetExpDist_vector
!  
!
! !AUTHORS
!   W. Ebenhoeh and C. Kohlmeier.
!
!
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
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: Ytur
  real(RLEN),dimension(NO_BOXES_XY)  :: Yirr
  real(RLEN),dimension(NO_BOXES_XY)  :: d,h,r
  real(RLEN),dimension(NO_BOXES_XY)  :: alpha
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! JM more Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: jnetY1,jnetY2,jnetY4,jnetY5
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !JM assign local variables to avoid pointer errors in insw_vector calls
  jnetY1=jnetYIc(iiY1,:)
  jnetY2=jnetYIc(iiY2,:)
  jnetY4=jnetYIc(iiY4,:)
  jnetY5=jnetYIc(iiY5,:)

  d=p_cturm

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature Response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! et  =   (p_q10)**(( ETW_Ben(:)- 10.0D+00)* 0.1D+00)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioturbation. ''turenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!JM  r=   Y1c(:)*insw_vector(jnetYIc(iiY1,:))+ &
!JM       Y4c(:)*insw_vector(jnetYIc(iiY4,:))
  r=   Y1c(:)*insw_vector(jnetY1)+ &
       Y4c(:)*insw_vector(jnetY4)
!JM  Ytur=r+   Y2c(:)*insw_vector(jnetYIc(iiY2,:)) &
!JM     + Y5c(:)*insw_vector(jnetYIc(iiY5,:)) 
  Ytur=r+   Y2c(:)*insw_vector(jnetY2) &
     + Y5c(:)*insw_vector(jnetY5) 

  turenh(:)  =DONE+   p_cmtur* MM_vector(  Ytur-r,  p_chtur)
  h  =DONE+   p_cmtur* MM_vector(  r,  p_chtur)
  turenh(:)=turenh(:)+h
  pxturinD1=h/turenh
  alpha=SetExpDist_vector(D6m,Dlm=p_clDxm)
  !less effective bioturbation when gradient (-alpha*exp(-p_cturm*alpha) of detritus becomes small
  r= p_Etur* turenh/(NZERO+D6m) *((DONE-pxturinD1)*(DONE-exp(-alpha*d)) + &
                        pxturinD1 *(DONE-exp(-alpha*D1m(:)))) 
  call flux_vector(iiBen, ppD6m,ppD6m,r)

  alpha=SetExpDist_vector(D7m,Dlm=p_clDxm)
  r= p_Etur* turenh/(NZERO+D7m) *((DONE-pxturinD1)*(DONE-exp(-alpha*d)) + &
                        pxturinD1 *(DONE-exp(-alpha*D1m(:))))
  call flux_vector(iiBen, ppD7m,ppD7m,r)

  alpha=SetExpDist_vector(D8m,Dlm=p_clDxm)
  r= p_Etur* turenh/(NZERO+D8m) *((DONE-pxturinD1)*(DONE-exp(-alpha*d)) + &
                        pxturinD1 *(DONE-exp(-alpha*D1m(:)))) 
  call flux_vector(iiBen, ppD8m,ppD8m,r)

  alpha=SetExpDist_vector(D9m,Dlm=p_clDxm)
  r= p_Etur* turenh/(NZERO+D9m) *((DONE-pxturinD1)*(DONE-exp(-alpha*d)) + &
                        pxturinD1 *(DONE-exp(-alpha*D1m(:)))) 
  call flux_vector(iiBen, ppD9m,ppD9m,r)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioirrigation. ''irrenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (sw_irr==1.or.sw_irr==3 ) then
    ! bacckground bioirrigation is caused by the wholes which are
    ! by the different types of benthic organisms.
!JM    Yirr=    Y2c(:)*insw_vector(jnetYIc(iiY2,:)) &
!JM    +        Y5c(:)*insw_vector(jnetYIc(iiY5,:)) &
!JM    +p_irrY4*Y4c(:)*insw_vector(jnetYIc(iiY4,:)) &
    Yirr=    Y2c(:)*insw_vector(jnetY2) &
    +        Y5c(:)*insw_vector(jnetY5) &
    +p_irrY4*Y4c(:)*insw_vector(jnetY4) &
    +p_irrY3*Y3c(:)*insw_vector(jnetY3c)
  
    irrenh(:)  =   DONE+ p_cmirr *MM_vector(  Yirr,  p_chirr)
    !irrigation is completely partly based on activtiy of deposit feeders which
    ! need oxyge from the pelagic to handle aboxic bacteria
    if (sw_irr==3) call BenOxygenDynamics(2)
  elseif (sw_irr==2) then
    !irrigation is completely based on activtiy of deposit feeders which
    ! need oxyge from the pelagic to handle aboxic bacteria
    irrenh(:)=DONE
    call BenOxygenDynamics(2)

  endif
  if (p_limDIm > ZERO) then
    r=-D6m*D6m/(p_limDim+D6m)*DONE ! adaptation rate 1 perday
    call LimitChange_vector(NEGATIVE,r,D6m,max_change_per_step)
    call flux_vector(iiBen, ppD6m,ppD6m, -r)
    r=-D7m*D7m/(p_limDim+D7m)*DONE ! adaptation rate 1 perday
    call LimitChange_vector(NEGATIVE,r,D7m,max_change_per_step)
    call flux_vector(iiBen, ppD7m,ppD7m, -r)
    r=-D8m*D8m/(p_limDim+D8m)*DONE ! adaptation rate 1 perday
    call LimitChange_vector(NEGATIVE,r,D8m,max_change_per_step)
    call flux_vector(iiBen, ppD8m,ppD8m, -r)
    r=-D9m*D9m/(p_limDim+D9m)*DONE ! adaptation rate 1 perday
    call LimitChange_vector(NEGATIVE,r,D9m,max_change_per_step)
    call flux_vector(iiBen, ppD9m,ppD9m, -r)
  endif

  call flux_vector(iiBen,ppirri_bio,ppirri_bio,(irrenh-irri_bio)/LocalDelta)

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
