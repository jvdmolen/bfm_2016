#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CorrectConcNearBed
!
! DESCRIPTION
!   Description of the anoxic diagenitic processes in the sediment
!
!   Calculate uptake by assuming a Rouse profile within the pelagic cell above the bed.
!
!   Rouse profile requires the assumption that:
!   -settling and resuspension of diatoms from the sea bed are in equilibrium
!   -the above process is not significantly affected by the take-up by filterfeeders
!
!   The first assumption would seem reasonable if, in the model, the flux of diatoms to the bed 
!   (other than by takeup by filter feeders) is set to zero (implying that pickup and deposition 
!   are equal), and changes in vertical diffusion in time are slow enough to allow the settling 
!   velocity to maintian a near-equilibrium concentration profile.
!
!   The second assumption is reasonable if the consumption by filter feeders is much smaller 
!   than ws*Ca, where ws is the settling velocity and Ca the near-bed diatom concentration.
!
!   Close to the sea bed, the Rouse profile can be approximated by a power law.
!
!   Then, the concentration C_f at height z_f above the bed where filter feeders are feeding 
!   can be approximated by (assuming linear increase of eddy diffusivity in the cell above the bed):
!
!   C_f = C(1) * ( z_f/(0.5*h(1)) )**-b
!
!   b = ws/(kappa*u_taub)
!
!   C(1) is concentration in centre of first cell above bed
!   h(1) is layer width of first cell above bed
!   ws is settling velocity of diatoms
!   kappa=0.4 is von Karman's constant
!   u_taub is shear velocity (from GETM/GOTM) 
!
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
      subroutine CorrectConcNearBed_vector(DepthLayer, Sedi, fto, &
                      p_p_max,n, correction)
!
! !USES:
     use global_mem,ONLY:RLEN,ZERO,NZERO,DONE
     use constants, only: SEC_PER_DAY
     use turbulence,  ONLY: kappa
     use mem,         ONLY: ETAUB,NO_BOXES_XY

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     !    Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      IMPLICIT NONE
! !INPUT PARAMETERS:
     integer,intent(IN)                    :: n
     REAL(RLEN), intent(IN),dimension(n)   ::DepthLayer
     REAL(RLEN), intent(IN),dimension(n)   ::Sedi
     REAL(RLEN), intent(IN)                ::fto
     REAL(RLEN), intent(IN)                ::p_p_max
!OUTPUT PARAMETERS:
     REAL(RLEN),dimension(n),intent(OUT)   ::correction
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Local Variables
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      REAL(RLEN),dimension(NO_BOXES_XY)                ::b
      REAL(RLEN),dimension(NO_BOXES_XY)                ::f
      REAL(RLEN),dimension(NO_BOXES_XY)                ::r
      REAL(RLEN),dimension(NO_BOXES_XY)                ::s
      REAL(RLEN),dimension(NO_BOXES_XY)                ::d
      REAL(RLEN),parameter                             ::p_small=1.0D-10
    
      s=ZERO;d=ZERO;
      f=  min(fto,0.5*DepthLayer)
      where ( abs(Sedi).gt.NZERO .and. DepthLayer.gt.NZERO) 
!JM        b = min(100.0D+00,Sedi/(p_small+SEC_PER_DAY*kappa*ETAUB(1))) 
!JM        s=max(DONE,(f/max(p_small,DepthLayer-f))**(-b))        
        b = (Sedi/SEC_PER_DAY)/(p_small+kappa*ETAUB(1))
        s=(f/(0.5D0*Depthlayer))**(-b)      
        correction=s;
     elsewhere (Sedi.gt.0)
        correction=DONE;
     elsewhere
        correction=ZERO;
     endwhere
!JM     where (correction *f/DepthLayer .gt.p_p_max) 
!JM        correction=p_p_max *DepthLayer/f
!JM     endwhere
     where (correction.gt.p_p_max) 
        correction=p_p_max
     endwhere
     return
     end

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: CorrectConcNearBed
!
! DESCRIPTION
!   Description of the anoxic diagenitic processes in the sediment
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
      subroutine CorrectConcNearBed(depthlayer, sedi, fto,p_p_max, &
                                                        correction)
!
! !USES:
     use global_mem,only:RLEN,LOGUNIT,ZERO,NZERO,DONE
     use constants, only: SEC_PER_DAY
     use turbulence,  ONLY: kappa
     use mem,         ONLY: ETAUB,NO_BOXES_XY

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     !    Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      IMPLICIT NONE
! !INPUT PARAMETERS:
     REAL(RLEN), intent(IN)   ::DepthLayer
     REAL(RLEN), intent(IN)   ::Sedi
     REAL(RLEN), intent(IN)   ::fto
     REAL(RLEN), intent(IN)  ::p_p_max
!OUTPUT PARAMETERS:
     REAL(RLEN),intent(OUT)   ::correction

!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Local Variables
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      REAL(RLEN):: b
      REAL(RLEN):: f
      REAL(RLEN):: r
      REAL(RLEN):: s
      REAL(RLEN):: d
      REAL(RLEN),parameter:: p_small=1.0D-10
    
      s=ZERO;d=ZERO;
      if ( Depthlayer.gt.NZERO ) then
        f=  min(fto,0.5*DepthLayer)
        if ( abs(Sedi).gt.NZERO) then 
!JM        b = min(100.0D+00,Sedi/(p_small+SEC_PER_DAY*kappa*ETAUB(1))) 
!JM        s=max(DONE,(f/max(1.D-10,DepthLayer-f))**(-b))        
          b = (Sedi/SEC_PER_DAY)/(p_small+kappa*ETAUB(1))
          s=(f/(0.5D0*Depthlayer))**(-b)      
          correction=s;
        else
           correction=DONE;
        endif
      else
           correction=ZERO;
      endif
!JM     if (correction *f/DepthLayer .gt.p_p_max) then
!JM        correction=p_p_max *DepthLayer/f
!JM     endif
     if (correction.gt.p_p_max) then
        correction=p_p_max
     endif
     return
     end

!EOC


