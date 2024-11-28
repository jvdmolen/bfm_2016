#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenDenitriDepth
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine CalcSulphideHorizon(BoxNumberXY,D1m,D2m,Dend,error,newD2m)
!
! !USES:

  ! ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE,dummy
  use mem_globalfun,   ONLY: rtsafe
  use mem_BenDenitriDepth, ONLY:boxnr

  implicit none
  integer ,intent(IN)              :: BoxNumberXY
  real(RLEN),intent(IN)            :: D1m
  real(RLEN),intent(IN)            :: D2m
  real(RLEN),intent(IN)            :: Dend
  integer,intent(OUT)              :: error
  real(RLEN),intent(OUT)           :: newD2m

  real(RLEN)                       ::Dxm
  external                         ::CalcShm


  Dxm=(D1m+D2m)*0.5

  boxnr=BoxNumberXY
  error = rtsafe( CalcSHm, D1m, Dend, 1.0D-6, newD2m)

  return
  end

  subroutine CalcSHm(Dxm,f,df)
  use mem,only: KNO3,KRED
  use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE,dummy
  use constants,ONLY: EQUATION,STANDARD,GET,DERIVATIVE,RFLUX
  use bennut_interface, ONLY: CalculateFromSet
  use mem_BenDenitriDepth, ONLY:boxnr,p_xmulti

  implicit none
  real(RLEN),intent(IN)            :: Dxm
  real(RLEN),intent(out)           :: f
  real(RLEN),intent(out)           :: df

  f=p_xmulti *CalculateFromSet(KNO3(boxnr),EQUATION,STANDARD,Dxm,dummy) &
             -CalculateFromSet(KRED(boxnr),EQUATION,STANDARD,Dxm,dummy) 
  df=p_xmulti *CalculateFromSet(KNO3(boxnr),DERIVATIVE,STANDARD,Dxm,dummy) &
             -CalculateFromSet(KRED(boxnr),DERIVATIVE,STANDARD,Dxm,dummy) 
  end
