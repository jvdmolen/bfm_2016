!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenQ1Transport
!
! DESCRIPTION
!   This file holds all alternative parameter values for parameters
!	defined in the process BenAnoxic.p.

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenQ1Transport
  use global_mem,ONLY: RLEN,LOGUNIT,NML_OPEN,NML_READ,NMLUNIT,ALLOC, &
                      error_msg_prn
  use mem,  ONLY: NO_BOXES_XY,iiBenLabileDetritus
!
  IMPLICIT NONE
  public
  real(RLEN),dimension(:,:),allocatable    :: jQIBTx
  real(RLEN),dimension(:,:),allocatable    :: jBTQIx
  real(RLEN),dimension(:),allocatable      :: Q1x
  real(RLEN),dimension(:),allocatable      :: Q11x
  real(RLEN),dimension(:),allocatable      :: Q21x

  integer             ::NUTR,layer,constituent
  integer,parameter   ::iiUrea=1,iiLOC=2
  real(RLEN)          :: rto,r50
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenQ1Transport PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  All parameter names and values were copied from the .p-file:
  real(RLEN)  :: p_p =0.0D+00       ! adsorption distribution coefficient
  real(RLEN)  :: p_slQ1=0.001D+00  ! minimum specific bacteria uptake per biomass of Q1 (arth. reasons)
  real(RLEN)  :: p_shQ1=20.00D+00  ! maximum specific bacteria uptake per biomass of Q1 (arth. reasons)
  real(RLEN)  :: p_pxtr=0.001D+00 ! part of Q2 which is tranported. Q2 is EPS produced by benPhyto
                                   ! which keep the benphyto sticking to the sediment and cannot
                                   ! be considered completely as a dissolved nutrient.

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenQ1Transport
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenQ1Transport()

  implicit none
  integer           ::status
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenQ1Transport_parameters/  p_p, p_slQ1,p_shQ1,p_pxtr
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  allocate(jQIBTx(1:iiBenLabileDetritus,1:NO_BOXES_XY),stat=status)
  if (status /= 0) call error_msg_prn(ALLOC,"BenQ1Transport","jQiBTx")
  allocate(jBTQIx(1:iiBenLabileDetritus,1:NO_BOXES_XY),stat=status)
  if (status /= 0) call error_msg_prn(ALLOC,"BenQ1Transport","jBTQix")
  allocate(Q1x(1:NO_BOXES_XY),stat=status)
  if (status /= 0) call error_msg_prn(ALLOC,"BenQ1Transport","Q1x")
  allocate(Q11x(1:NO_BOXES_XY),stat=status)
  if (status /= 0) call error_msg_prn(ALLOC,"BenQ1Transport","Q11x")
  allocate(Q21x(1:NO_BOXES_XY),stat=status)
  if (status /= 0) call error_msg_prn(ALLOC,"BenQ1Transport","Q21x")


  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenQ1Transport parameters.."
   open(NMLUNIT,file='BenQ1Transport.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=BenQ1Transport_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenQ1Transport_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenQ1Transport.f90","BenQ1Transport.nml")
101 call error_msg_prn(NML_READ,"InitBenQ1Transport.f90","BenQ1Transport_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine InitBenQ1Transport
  end module mem_BenQ1Transport
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
