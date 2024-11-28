#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: AttachRatetoTEP
!
! DESCRIPTION
!   This process describes the dynamics of the special Phaeocyctis
!    features in the ERSEM model. The differences in behaviour
!    are expressed by differences in parameter-values only.
!
!
!

!
! !INTERFACE
  module AttachRatetoTEP
!
! !USES:

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,ZERO_KELVIN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   contains

  function RelAttachRateToTEP(mode,n,size_cell,E, &
               radius_macro,number_macros,ref)
  use global_interface,ONLY: CalcSh

  use mem,only:ESW,ETW
  use constants,only: SEC_PER_DAY
  implicit none
  integer,intent(IN)                          ::mode
  integer,intent(IN)                          ::n
  real(RLEN),dimension(n),intent(IN)          ::size_cell                        
  real(RLEN),dimension(n),intent(IN)          ::E
  real(RLEN),intent(IN),dimension(n),optional ::radius_macro
  real(RLEN),intent(IN),dimension(n),optional ::number_macros
  real(RLEN),intent(IN),dimension(n),optional ::ref
  real(RLEN),dimension(n)                     ::RelAttachRateToTEP                        

  real(RLEN),dimension(n)                     ::Sh
  real(RLEN),dimension(n)                     ::diff
  real(RLEN),dimension(n)                     ::strange_density
  real(RLEN),dimension(n)                     ::rx_any,sphere
  real(RLEN),parameter                        ::rPI= 3.1415926535897D+00
  real(RLEN),parameter                        ::r4d3=1.3333333333333D+00
  real(RLEN),parameter                        ::four=4.0E+00

  
  diff=CalcDiffSphere(n,ESW,ETW,size_cell)
  sh=CalcSH(E,size_cell,diff)
  select case (mode)
    case (1)
      rx_any=four*rPI*diff*radius_macro*number_macros*SEC_PER_DAY
    case (2)
      strange_density=(r4d3*rPI* size_cell**3)
      ! as a porxy for the "relative" weight of a cell is the volume taken...
      rx_any=Sh*four*rPI*diff/strange_density ! unit 1/M3/d
      ! dimension less number: difference in attach rate relative to reference
      if (present(ref)) rx_any=rx_any/ref 
  end select
  RelAttachRateToTEP=rx_any                        
  return
  end function RelAttachRateToTEP

  function CalcDiffSPhere(n,ESW,ETW,size_cell)
    use global_interface,ONLY:CalcWaterProp

    integer,intent(IN)                 ::n
    real(RLEN),dimension(n),intent(IN) ::ESW
    real(RLEN),dimension(n),intent(IN) ::ETW
    real(RLEN),dimension(n),intent(IN) ::size_cell
    real(RLEN),dimension(n)            :: CalcDiffSphere !output:m2/s

    real(RLEN),parameter          ::K_B=3.806e-23 !Boltmann:m2⋅kg/(s2⋅K)
    real(RLEN),parameter          ::rPI= 3.1415926535897D+00
    real(RLEN),parameter          ::six= 6.0D+00
    real(RLEN),dimension(n)::Mu !viscosity

    call CalcWaterProp(n,ESW,ETW,Mu=Mu)
    CalcDiffSphere= K_B* (ETW- ZERO_KELVIN)/(six*rPI*Mu*size_cell)
  end function CalcDiffSPhere
  end module AttachRatetoTEP

