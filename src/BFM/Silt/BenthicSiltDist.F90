#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicSiltDist
!
! DESCRIPTION
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicSiltDist(mode,rate)
!
! !USES:


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NZERO,DONE,ZERO
  use mem,only:DSM,QSx,jrESS,turenh,ppDSM,ppQSx,D1m, &
               flux_vector,NO_BOXES_XY,pxturinD1,iiBen,max_change_per_step,LocalDelta
  use mem_Bioturbation,only:p_ETur,p_cturm
  use mem_globalfun,   ONLY: SetExpDist_vector
  use Constants, ONLY:NEGATIVE,POSITIVE
  use mem_Param,  ONLY: p_clDxm,p_d_tot
  use bio_var,ONLY:sw_BenSilt
  use LimitRates, ONLY:LimitChange_vector


  implicit none
  integer,intent(IN)                                       :: mode
  real(RLEN),dimension(NO_BOXES_XY),intent(IN)             ::rate
  
  real(RLEN),dimension(NO_BOXES_XY)                         ::alpha
  real(RLEN),dimension(NO_BOXES_XY)                         ::r,m
  reaL(RLEN),dimension(NO_BOXES_XY)                         ::thlayer

  if (sw_BenSilt==0 ) return
  if ( mode.eq.0) then
     alpha=SetExpDist_vector(DSm,Dlm=p_clDxm)
     !less effective bioturbation when gradient (-alpha*exp(-p_cturm*alpha) 
     ! of detritus becomes small
     r=p_Etur*turenh*((DONE-pxturinD1)/(NZERO+DSm)*(DONE-exp(-alpha*p_cturm))+ &
                        pxturinD1 /(NZERO+DSm)*(DONE-exp(-alpha*D1m(:)))) 
     call flux_vector(iiBen, ppDSm,ppDSm,r)

     call CouplingBioSilt(2,thlayer)
     thlayer=abs(thlayer*LocalDelta)
     r=log(DONE -0.5*(DONE-exp(-thlayer/DSm)))*(-DSm)
     where (jrESS.gt.ZERO) 
      r=jrESS/QSx*(-DSm)
     elsewhere
       m=QSx-jrESS*LocalDelta
       r=jrESS/m*(-Dsm)
     endwhere
     r=r-DONE*max(ZERO,DSM-p_d_tot)
     call LimitChange_vector(NEGATIVE,r,DSm,max_change_per_step)
     call flux_vector(iiBen, ppDSm,ppDSm,r)
     call flux_vector(iiBen, ppQSx,ppQSx,-jrESS)

! else
  endif

  end
