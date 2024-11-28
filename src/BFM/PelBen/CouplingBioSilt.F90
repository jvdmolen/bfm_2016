#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   This process calculate the variables which determine
!   the coupling between siltdynmaics and the ecology  
!
!
! !INTERFACE
  subroutine CouplingBioSilt(mode0,out)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE
      use mem,ONLY:NO_BOXES_XY,jrESS,BP1c,DSm,QSx
      use mem_Param,ONLY:p_poro,p_d_tot
      use mem_benphyto,ONLY:p_mxSiltBPc
      use bio_var,only:sw_BenSilt
      use mem_Silt,ONLY: p_SampleDepth
      use mem_globalfun,   ONLY: IntegralExpDist_vector,SetExpDist_vector

      implicit none
      integer,intent(IN)                             :: mode0
      REAL(RLEN),dimension(NO_BOXES_XY),intent(OUT)  ::out
      
      REAL(RLEN),dimension(NO_BOXES_XY)              :: psilt
      REAL(RLEN),dimension(NO_BOXES_XY)              :: silt
      REAL(RLEN),dimension(NO_BOXES_XY)              :: factor
      REAL(RLEN),dimension(NO_BOXES_XY)              :: alpha
      REAL(RLEN),dimension(NO_BOXES_XY)              :: r_xy1,r_xy2
      integer                                        :: mode

      mode=abs(mode0)
      factor=p_mxSiltBPc/(BP1c+p_mxSiltBPc)
      if (sw_BenSilt<2.or.mode0<0) then
        psilt=max(0.1,(p_poro - 0.38662 )/ 0.415*factor)
      else
        alpha=SetExpDist_vector(Dsm,Dmm=p_d_tot)
        r_xy1=p_d_tot; r_xy2=p_sampleDepth
        psilt= QSx/ IntegralExpDist_vector(-alpha, r_xy1)* &
              IntegralExpDist_vector( -alpha,r_xy2)/(2.65E+09*r_xy2)
      endif
      silt= 2.65e+09 *psilt ! 2.65e+9 mg silt per m3
        ! rrESS is reuspension/d-->
        !               in one timestep rrESS*LocalDelta is resuspended.
        ! mg Silt/m2/d*d/(mg silt/m2)*1m
      select case (mode)
        case (0) ; out=psilt
        case (1) ; out=silt
        case (2) ; out=jrESS/silt
      end select
      end
