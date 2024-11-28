!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleInterface
!
! DESCRIPTION
!   Definition of Explicit Interfaces

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE global_interface

  !

  !
!
! !AUTHORS
!   mfstep/ERSEM team, especially  Momme Butenschoen:
!
! !REVISION_HISTORY
!   ---
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
  INTERFACE

  subroutine MesoZooDynamics(zoo, ppzooc, ppzoon, ppzoop)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: zoo
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop
  end subroutine MesoZooDynamics
  end INTERFACE

  INTERFACE

  subroutine MicroZooDynamics(zoo, ppzooc, ppzoon, ppzoop)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: zoo
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop
  end subroutine MicroZooDynamics
  end INTERFACE

  INTERFACE

  subroutine PhytoDynamics(phyto, ppphytoc, ppphyton, ppphytop, ppphytos, &
    ppphytol)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol
  end subroutine PhytoDynamics
  end INTERFACE

  INTERFACE

  subroutine PhotoAvailableRadiationDynamics(phyto, ppphytoc, ppphyton, &
    ppphytop, ppphytos, ppphytol)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol
  end subroutine PhotoAvailableRadiationDynamics
  end INTERFACE

  INTERFACE

  subroutine LightAdaptationDynamics(phyto, ppphytoc, ppphyton, ppphytop, &
    ppphytos, ppphytol)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol
  end subroutine LightAdaptationDynamics
  end INTERFACE

  INTERFACE

  subroutine BenOrganismDynamics(y, ppyc, ppyn, ppyp)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: y
  integer,intent(IN) :: ppyc
  integer,intent(IN) :: ppyn
  integer,intent(IN) :: ppyp
  end subroutine BenOrganismDynamics
  end INTERFACE

  INTERFACE

  subroutine BenBacDynamics(hx, pphxc, pphxn, pphxp)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: hx
  integer,intent(IN) :: pphxc
  integer,intent(IN) :: pphxn
  integer,intent(IN) :: pphxp
  end subroutine BenBacDynamics
  end INTERFACE

  INTERFACE

  FUNCTION eTq(temp,p_q10)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: temp
  real(RLEN),intent(IN) :: p_q10
  real(RLEN) :: eTq
  end FUNCTION eTq
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcChlorophylla()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcChlorophylla
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcVerticalExtinction(mode)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: mode
  end SUBROUTINE CalcVerticalExtinction
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcLight
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcLight
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcOxygenSaturation()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcOxygenSaturation
  end INTERFACE

  INTERFACE
  subroutine SurfaceO2Diffusion(mode,t,Wind,depth,rdiff,flux)
  use global_mem, ONLY:RLEN
  IMPLICIT NONE
  integer,INTENT(IN)        :: mode
  real(RLEN), INTENT(IN)    :: t
  real(RLEN), INTENT(IN)    :: Wind
  real(RLEN), INTENT(IN)    :: depth
  real(RLEN), INTENT(IN)    :: rdiff
  real(RLEN), INTENT(OUT)   :: flux
  end subroutine SurfaceO2Diffusion
  end INTERFACE

  INTERFACE
  subroutine SurfaceCO2Diffusion(mode,t,s,rho,Wind,depth,rdiff,flux)
  use global_mem, ONLY:RLEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer,INTENT(IN)        :: mode
  real(RLEN), INTENT(IN)    :: t
  real(RLEN), INTENT(IN)    :: s
  real(RLEN), INTENT(IN)    :: rho
  real(RLEN), INTENT(IN)    :: Wind
  real(RLEN), INTENT(IN)    :: depth
  real(RLEN), INTENT(IN)    :: rdiff
  real(RLEN), INTENT(OUT)   :: flux
  end subroutine SurfaceCO2Diffusion
  end INTERFACE


   interface
   subroutine SurfaceCO2Processes(BoxNumber,BoxNumberXY)
    IMPLICIT NONE
    integer, INTENT(IN)    :: BoxNumber
    integer, INTENT(IN)    :: BoxNumberXY
  end subroutine SurfaceCO2Processes
  end interface


  INTERFACE

  SUBROUTINE ResetTotMassVar()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE ResetTotMassVar
  end INTERFACE

  INTERFACE

  FUNCTION BoxAbove(box_no)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: box_no
  integer :: BoxAbove
  end FUNCTION BoxAbove
  end INTERFACE

  INTERFACE

  FUNCTION BoxBeneath(box_no)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: box_no
  integer :: BoxBeneath
  end FUNCTION BoxBeneath
  end INTERFACE

  INTERFACE

  FUNCTION CalcSchmidtNumberOx(Temp)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: Temp
  real(RLEN) :: CalcSchmidtNumberOx
  end FUNCTION CalcSchmidtNumberOx
  end INTERFACE

  INTERFACE

  FUNCTION CalcSchmidtNumberCO2(Temp)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: Temp
  real(RLEN) :: CalcSchmidtNumberCO2
  end FUNCTION CalcSchmidtNumberCO2
  end INTERFACE

 interface
    subroutine RecalcPenetrationDepth(D1m,Dtm, Dxm,Dhm, input, mass,newDxm )
        use global_mem, ONLY:RLEN
        REAL(RLEN),intent(IN)     ::D1m
        REAL(RLEN),intent(IN)     ::Dtm
        REAL(RLEN),intent(IN)     ::Dxm
        REAL(RLEN),intent(IN)     ::Dhm
        REAL(RLEN),intent(IN)     ::input
        REAL(RLEN),intent(IN)     :: mass
        REAL(RLEN),intent(OUT)    :: newDxm
    end subroutine RecalcPenetrationDepth
  end interface

  interface
     subroutine CorrectConcNearBed_vector(depthlayer, sedi, fto, &
                      p_p_max,n, correction)
     use global_mem,ONLY:RLEN
     integer,intent(IN)                    :: n
     REAL(RLEN), intent(IN),dimension(n)   ::DepthLayer
     REAL(RLEN), intent(IN),dimension(n)   ::Sedi
     REAL(RLEN), intent(IN)                ::fto
     REAL(RLEN), intent(IN)                ::p_p_max
     REAL(RLEN),dimension(n),intent(OUT)   ::correction
     end subroutine CorrectConcNearBed_vector
  end interface

  interface
     subroutine CorrectConcNearBed(depthlayer, sedi, fto,p_p_max, &
                                                        correction)
     use global_mem,only:RLEN
     IMPLICIT NONE
     REAL(RLEN), intent(IN)   ::DepthLayer
     REAL(RLEN), intent(IN)   ::Sedi
     REAL(RLEN), intent(IN)   ::fto
     REAL(RLEN), intent(IN)  ::p_p_max
     REAL(RLEN),intent(OUT)   ::correction
     end subroutine CorrectConcNearBed
  end interface

  interface
  subroutine BenPhytoDynamics(bp, ppbpc, ppbpn, ppbpp, ppbps, ppbpl)
    use global_mem, ONLY:RLEN
    IMPLICIT NONE
    integer,intent(IN)  :: bp
    integer,intent(IN) :: ppbpc
    integer,intent(IN) :: ppbpn
    integer,intent(IN) :: ppbpp
    integer,intent(IN) :: ppbps
    integer,intent(IN) :: ppbpl
  endsubroutine BenPhytoDynamics
  end interface

  interface
  subroutine FindNaNInRates(iiSys,ppState,message)
  use global_mem, ONLY:RLEN
  IMPLICIT NONE
  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message
  end subroutine FindNaNInRates
  end interface

  INTERFACE
  subroutine BenQ1TransportDynamics(KQ1,mode,ppR1x,R1_Benx,p_qnUc)
    use global_mem, ONLY:RLEN
    use mem,ONLY:NO_BOXES_XY
    IMPLICIT NONE
    integer,dimension(NO_BOXES_XY),intent(IN)           ::KQ1
    integer,intent(IN)                                   ::mode
    integer,intent(IN),optional                          ::ppR1x
    real(RLEN),optional                                  ::p_qnUc
    real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional::R1_Benx
  end subroutine BenQ1TransportDynamics
  end INTERFACE


  INTERFACE
  subroutine BenQ1TransportCalculate(KQ1,mode,R1_Benx,p_qnUc)
    use global_mem, ONLY:RLEN
    use mem,ONLY:NO_BOXES_XY
    IMPLICIT NONE
    integer,dimension(NO_BOXES_XY),intent(INOUT)         ::KQ1
    integer,intent(IN)                                   ::mode
    real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional::R1_Benx
    real(RLEN),optional                                  ::p_qnUc
  end subroutine BenQ1TransportCalculate
  end INTERFACE

  interface
  subroutine PhaeocystisCalc(mode,phyto,output,input1,param)
    use global_mem, ONLY:RLEN
    use mem, ONLY: NO_BOXES
    implicit none
    integer,intent(IN)                         :: mode
    integer,intent(IN)                         :: phyto
    real(RLEN),dimension(NO_BOXES),intent(INOUT) :: output
    real(RLEN),dimension(NO_BOXES),intent(IN)  :: input1
    real(RLEN),intent(IN)  :: param
  end subroutine PhaeocystisCalc
  end interface

  INTERFACE
  subroutine PhaeocystisCalc_1l(mode,phyto,output,input1,param)
  use global_mem, ONLY:RLEN
  use mem,ONLY:NO_BOXES_XY
  implicit none
  integer,intent(IN)                            :: mode
  integer,intent(IN)                            :: phyto
  real(RLEN),intent(OUT),dimension(NO_BOXES_XY) :: output
  real(RLEN),intent(IN),dimension(NO_BOXES_XY)  :: input1
  real(RLEN),intent(IN)                         :: param
  end subroutine PhaeocystisCalc_1l
  end INTERFACE

  INTERFACE
  subroutine FilterFeederDynamics(y,ppyc,ppyn,ppyp)
    IMPLICIT NONE
    integer,intent(IN)  :: y
    integer,intent(IN) :: ppyc
    integer,intent(IN) :: ppyn
    integer,intent(IN) :: ppyp
  end subroutine FilterFeederDynamics
  end interface

  interface
  SUBROUTINE CalcBndyConcentration(area,side,numc,kmax,mode,depend, &
                           multi,rho,salt,temp, cc_old,cc_inout,error)
  use global_mem, ONLY:RLEN
  IMPLICIT NONE
  character(len=*),intent(IN)            :: area
  character(len=*),intent(IN)            :: side
  integer,intent(IN)                     :: numc,kmax,mode,depend
  real(RLEN),intent(IN)                  :: multi
  real(RLEN),intent(IN)                  :: temp(0:kmax)
  real(RLEN),intent(IN)                  :: salt(0:kmax)
  real(RLEN),intent(IN)                  :: rho(0:kmax)
  real(RLEN),intent(IN)                  :: cc_old(1:numc,0:kmax)
  real(RLEN),intent(INOUT)               :: cc_inout(1:numc,0:kmax)
  integer,intent(OUT)                    :: error
  end SUBROUTINE CalcBndyConcentration
  end interface

  INTERFACE
  SUBROUTINE CalcRiverConcentration(mode,kmax,numc,ETW,cc_old_river,cc_river)
  use global_mem, ONLY:RLEN
  IMPLICIT NONE
  integer,intent(IN)         :: mode
  integer,intent(IN)         :: kmax
  integer,intent(IN)         :: numc
  real(RLEN),intent(IN)      :: ETW(1:kmax)
  real(RLEN),intent(IN)      :: cc_old_river(1:numc,1:kmax)
  real(RLEN),intent(INOUT)   :: cc_river(1:numc,1:kmax)
  end SUBROUTINE CalcRiverConcentration
  end INTERFACE

  INTERFACE
  subroutine check_if_in_output(text,put_in_list)
  IMPLICIT NONE
  character(len=*),intent(IN)  :: text
  logical,intent(OUT)          :: put_in_list
  end subroutine check_if_in_output
  end INTERFACE

  INTERFACE
  subroutine findnan( vector,n,iout)
      use global_mem, only:RLEN,LOGUNIT
      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      integer,intent(OUT)              :: iout
  end subroutine findnan
  end INTERFACE

  INTERFACE
       function CorrectForDiffBoundLayer(mol_diff,r_single,qu,T,Cinf,r,cells)
       use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
       use mem,only: NO_BOXES

       implicit none
       real(RLEN),intent(IN)          ::mol_diff
       real(RLEN),intent(IN)          ::r_single
       real(RLEN),intent(IN)          ::qu
       real(RLEN),intent(IN)          ::T(NO_BOXES)
       real(RLEN),intent(IN)          ::Cinf(NO_BOXES)
       real(RLEN),intent(IN),optional ::r(NO_BOXES)
       real(RLEN),intent(IN),optional ::cells(NO_BOXES)
       real(RLEN)                     ::CorrectForDiffBoundLayer(NO_BOXES)
       end function CorrectForDiffBoundLayer
  end INTERFACE
  INTERFACE
      function CalcSh(E, diameter, diffusion)
! !USES:
      use global_mem, ONLY:RLEN,DONE
      use mem,ONLY:NO_BOXES

! !INPUT:
      implicit none
      REAL(RLEN),intent(IN)         ::E(NO_BOXES)
      REAL(RLEN),intent(IN)         ::diameter(NO_BOXES)
      REAL(RLEN),intent(IN)         ::diffusion(NO_BOXES)
      real(RLEN)                    ::CalcSh(NO_BOXES)
      end function CalcSh
  end INTERFACE
  INTERFACE
       subroutine CalcWaterProp(n,Salt,Temp,Rho,Mu)
       use global_mem,ONLY:RLEN
       IMPLICIT NONE
       integer,intent(IN)                              ::n
       real(RLEN),dimension(n),intent(IN)              ::Salt
       real(RLEN),dimension(n),intent(IN)              ::Temp
       real(RLEN),dimension(n),intent(OUT),optional    ::Rho
       real(RLEN),dimension(n),intent(OUT),optional    ::Mu
       end subroutine CalcWaterProp
  end INTERFACE

  INTERFACE
function FindMidPointExp(Dm,p_from,p_to)
      use global_mem,ONLY:RLEN
      use mem,only: NO_BOXES_XY
      implicit none
      real(RLEN),intent(IN)          :: Dm
      real(RLEN),intent(IN)          :: p_from
      real(RLEN),intent(IN)          :: p_to
      real(RLEN)                     :: FindMidPointExp
      end function FindMidPointExp
  end INTERFACE
  INTERFACE
function FindMidPointExp_vector(Dm,p_from,d1_from,p_to,d1_to)
      use global_mem,ONLY:RLEN
      use mem,only: NO_BOXES_XY
      implicit none
      real(RLEN),dimension(NO_BOXES_XY),intent(IN)          :: Dm
      real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional :: d1_from
      real(RLEN),intent(IN),optional                        :: p_from
      real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional :: d1_to
      real(RLEN),intent(IN),optional                        :: p_to
      real(RLEN),dimension(NO_BOXES_XY)                     :: FindMidPointExp_vector
      end function FindMidPointExp_vector
  end INTERFACE

  INTERFACE
         subroutine TopPredLosses(iisys,n,spec_loss,p_ex, &
                   state_c,state_n,state_p, &
                   flO3c,flUrean,flUreac,flN1p,flR6c,flR6n,flR6p)
         USE global_mem, ONLY:RLEN 
         implicit none
         integer, intent(IN)           ::iisys
         integer, intent(IN)           ::n
         real(RLEN),intent(IN)         ::spec_loss(n)
         real(RLEN),intent(IN)         ::p_ex
         real(RLEN), intent(IN)        ::state_c(n)
         real(RLEN), intent(IN)        ::state_n(n)
         real(RLEN), intent(IN)        ::state_p(n)
         real(RLEN), intent(INOUT)     ::flO3c(n)
         real(RLEN), intent(INOUT)     ::flUrean(n)
         real(RLEN), intent(INOUT)     ::flUreac(n)
         real(RLEN), intent(INOUT)     ::flN1p(n)
         real(RLEN), intent(INOUT)     ::flR6c(n)
         real(RLEN), intent(INOUT)     ::flR6n(n)
         real(RLEN), intent(INOUT)     ::flR6p(n)
         end subroutine TopPredLosses
  end INTERFACE

  INTERFACE
    function CalcPelMassInM2(ppX)
     USE global_mem, ONLY:RLEN
     use mem,only:NO_BOXES_XY
     integer,intent(IN)      ::ppX
     real(RLEN),dimension(NO_BOXES_XY)  ::CalcPelMassINM2
    end function CalcPelMassInM2
  end INTERFACE




#ifdef INCLUDE_MACROPHYT
  INTERFACE

  subroutine MacroPhytoDynamics(macro, ppmcsc, ppmcac, ppmcan, ppmcap, ppmcal)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: macro
  integer,intent(IN) :: ppmcsc
  integer,intent(IN) :: ppmcac
  integer,intent(IN) :: ppmcan
  integer,intent(IN) :: ppmcap
  integer,intent(IN) :: ppmcal
  end subroutine MacroPhytoDynamics
  end INTERFACE
#endif
  end MODULE global_interface
