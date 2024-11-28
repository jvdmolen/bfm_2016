#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!   !    This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!        - dissolution of biogenic silica
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelChemDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: N4n, N3n, O2o, &
  ! N6r, R6s, N5s, P1s
  ! The following Pelagic 1-d global boxvars are modified : flN3O4n
  ! The following Pelagic 1-d global boxvars  are used: ETW, flPIR6s
  ! The following 0-d global parameters are used: p_qon_nitri, p_qro, &
  ! p_qon_dentri
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NZERO,ZERO,ZERO_KELVIN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: N4n, N3n, O2o, N6r, R6s,R2c ,P6c,R3c
#endif
  use mem, ONLY: ppN4n, ppN3n,ppN1p,ppN5s,  ppO2o, ppN6r, &
     ppR6n,ppR6p,ppR6s, ppO3c,ppR1c,ppR1n,ppR1p,ppR3c,ppR2c,&
     iiPhytoPlankton,ppPhytoPlankton,iiN,iiP,iiS,iiPel,   &
    flN3O4n, flN3N4n,flN4N3n,ETW, flPIR6s,flPTN6r,NO_BOXES,&
    flux_vector, sN4N3n,pH, &
    flPIR6n,flPIR1n,flPIR6p,flPIR1p,flR3R2c,flnDIp,flnDIn,flR1N4n,flR1O3c, &
    Depth,  max_change_per_step,iiTotal
  use mem_Param,  ONLY: p_qon_nitri, p_qro, p_qon_dentri,CalcPhytoPlankton
  use constants,  ONLY: MW_C,POSITIVE
  use mem_PelChem
  use mem_PelBac,only:p_version_PelBac=>p_version
  use LimitRates, ONLY:LimitChange_vector
! use mem_phyto,only:p_qnRc,p_xdiv,p_lN3N4n,p_qun
!use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4,iiP2,Source

#ifdef INCLUDE_PELCO2
  use mem,ONLY:CO2
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector, eTq_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: MM_vector, eTq_vector, insw_vector
  use BFM_ERROR_MSG,only:set_warning_for_getm

!  
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
!
!
! !REVISION_HISTORY
!   !
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
  integer                         :: i,iout
  integer                         :: j
  real                            :: rscalar
  real(RLEN),dimension(NO_BOXES)  :: fN6O2r
  real(RLEN),dimension(NO_BOXES)  :: eo
  real(RLEN),dimension(NO_BOXES)  :: er
  real(RLEN),dimension(NO_BOXES)  :: fR6N5s
  real(RLEN),dimension(NO_BOXES)  :: sOs
  real(RLEN),dimension(NO_BOXES)  :: o,r,s,limit_1,limit_2


  !
  ! Regulating factors
  !

  eo  =   MM_vector(  max(NZERO,O2o(:)),  p_clO2o)
  er  =   MM_vector(  max(NZERO,N6r(:)),  p_clN6r)

  !fluxes due to breakdown of urea
   call flux_vector(iiPel,ppR1c,ppO3c,flR1O3c)
   call flux_vector(iiPel,ppR1n,ppN4n,flR1N4n)
  !
  ! Nitrification in the water
  !
  
  select case (p_version_PelBac)
     case (5)      ;flN4N3n(:)=sN4N3n * N4n
     case default  ; flN4N3n(:)=max(ZERO,p_sN4N3* N4n(:)*  &
                      eTq_vector( ETW, p_q10N4N3)* eo) * N4n/(0.5+N4n)
  end select

  if (p_version_PelBac.ne.5) then
#ifdef INCLUDE_PELCO2
    !limitation at hight pH == low [CO2] in oxic layer : no growth possible
    ! prim prod due to nitrification recalculated to mmolC/m3 porewater
    r  = p_cyn * flN4N3n *MW_C   
    call LimitChange_vector(POSITIVE,r,CO2,max_change_per_step,o)
    flN4N3n=flN4N3n*o
#endif
  endif
  call flux_vector( iiPel, ppO2o,ppO2o,-( flN4N3n(:)* p_qon_nitri) )

  call flux_vector( iiPel, ppN4n,ppN3n,   flN4N3n(:) )
  !
  ! Denitrification in the water
  !  this process takes only place in cases of "negative oyxgen concentration"

  r= p_pN3O4*flPTN6r/p_qro*p_qon_dentri
  s= p_pN3O4*flPTN6r
  !reclalculate to o-units
  !Check only on negative fluxes. and limit if necessary.......
  call LimitChange_vector(POSITIVE,r,N3n,max_change_per_step,limit_1)
  call LimitChange_vector(POSITIVE,s,N6r,max_change_per_step,limit_2)
  limit_1=min(limit_1,limit_2)
  flN3O4n=limit_1*r
  call flux_vector( iiPel, ppN3n,ppN3n,-flN3O4n(:) )
  call flux_vector( iiPel, ppN6r,ppN6r, -limit_1*s)
  
  rscalar=sum(flN3O4n*Depth)
  if (rscalar.gt.ZERO) then
    j=sum(insw_vector(-(O2o-p_clO2o)));
    if ( j.gt.0) then
      i=sum(insw_vector((flN3O4n)));
      write(LOGUNIT,'(A,I2,A,F6.2,A,I2,A)') '  pelagic denitrification in',i, &
       ' , low oxygen (<',p_clO2o,') in ',j,' layer(s)'
      call set_warning_for_getm
    endif
  endif

  !
  ! Reoxidation of reduction equivalents
  !
  if ( p_sOS.eq.ZERO) then
     r=7.0
#ifdef INCLUDE_PELCO2
     where ( pH.gt.ZERO)  r=pH
#endif
       !equation of Millero Environ Sci Technol. 1987,21,439-443
       ! chosen for constant ion strength of 0.7 (seawater) : 0.3681=0.44*0.7**0.5
       sOS=10.0**(11.78-3.0e3/(-ZERO_KELVIN+ETW)+0.16*r+0.3681) * 24.E-6
       sOs=sOs*O2o
  else
      sOs  =   p_sOS* eo * er
  endif
  fN6O2r  =   sOS* N6r(:)

  r=fN6O2r/p_qro
  call LimitChange_vector(POSITIVE,r,O2o,max_change_per_step,o)
  call LimitChange_vector(POSITIVE,fN6O2r,N6r,max_change_per_step,r)
  fN6O2r=fN6O2r*min(o,r,N6r/(N6r+0.001))
  call flux_vector( iiPel, ppN6r,ppN6r,-( fN6O2r) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( fN6O2r/ p_qro) )

  !
  ! Regeneration of dissolved silica
  !

  fR6N5s  =   p_sR6N5* eTq_vector( ETW(:), p_q10R6N5)* R6s(:)  &
                                                   *100.0/(R2c(:)+100.0) 
  call flux_vector( iiPel, ppR6s,ppN5s, fR6N5s )

  do i=1,iiPhytoPlankton
    if ( CalcPhytoPlankton(i) ) then
      call flux_vector( iiPel, ppPhytoPlankton(i,iiN),ppR6n, flPIR6n(i,:) )
      call flux_vector( iiPel, ppPhytoPlankton(i,iiN),ppR1n, flPIR1n(i,:) )
      call flux_vector( iiPel, ppPhytoPlankton(i,iiP),ppR6p, flPIR6p(i,:) )
      call flux_vector( iiPel, ppPhytoPlankton(i,iiP),ppR1p, flPIR1p(i,:) )
      j=ppPhytoPlankton(i,iiS)
      if ( j.gt.0)  call flux_vector( iiPel, j,       ppR6s, flPIR6s(i,:) )
    endif
  enddo
  call findnan(flR3R2c,NO_BOXES,iout)
  if ( iout>0) write(LOGUNIT,*) 'PelChem lR3R2=NAN',iout,P6c(iout),flR3R2c(iout) 
  call LimitChange_vector(POSITIVE,flR3R2c,R3c,max_change_per_step)
  call flux_vector(iiPel,ppR3c,ppR2c,flR3R2c(:))
  call flux_vector(iiPel,ppN3n,ppN4n,flN3N4n(:))

! flnDIn=Source_D3_vector(ppN4n,iiTotal) + Source_D3_vector(ppN3n,iiTotal)
! flnDIp=Source_D3_vector(ppN1p,iiTotal)


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
