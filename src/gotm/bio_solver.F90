!$Id: ode_solvers.F90,v 1.9 2005-12-13 14:12:45 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_solver --- integration of biological models
!
! !INTERFACE:
   module bio_solver
!
! !DESCRIPTION:
! This module is a wrapper for all the ode solvers.
! It is meant to simplify the exchange of information between the bio 
! models and the memory system.
! IMPORTANT NOTE: process_model is always called with cc for consistency
! with GOTM but that array is not passed in the call to BFM
! 4th order schemes not implemented yet!
!
! !USES:
#ifdef BFM_GOTM
   use constants, only: SEC_PER_DAY
   use bio_var, only: bio_setup,pp,dd,numbc,ccb,ppb,ddb
   use bio_bfm, only: reset_diagonal
!  use mem,only:ppY2p
#endif
!
!  default: all is private.
   private
   integer,parameter :: NMAXVAR=100 ! max number of variables
   integer,parameter :: NLEVB=1     ! number of benthic levels
   logical,save                        :: first=.TRUE.
   integer                             :: ci,i,j
!  REALTYPE                            :: rhs0
!  REALTYPE                            :: ppsum0 !,ddsum0
   REALTYPE                            :: a(1:NMAXVAR,1:NMAXVAR)
   REALTYPE                            :: r(1:NMAXVAR)
   REALTYPE,allocatable,dimension(:,:) :: cc0,cc1,ppsum,ddsum, &
                                          rhs,rhs1,rhs2,rhs3
#ifdef BFM_GOTM
   REALTYPE,allocatable,dimension(:,:) :: ccb0,ccb1,ppbsum,ddbsum
#else
   REALTYPE,allocatable,dimension(:,:,:) :: pp,dd
#endif
!  integer                               ::iout
!
! !PUBLIC MEMBER FUNCTIONS:
  public ode_solver
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!                    : Marcello Vichi & Piet Ruardij
!
! !PRIVATE DATA MEMBERS:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: General ODE solver \label{sec:ode-solver}
!
! !INTERFACE:
   subroutine ode_solver(solver,numc,nlev,dt,cc)
!
! !DESCRIPTION:
! Here, 10 different numerical solvers for the right hand sides of the
! biogeochemical models are implemented for computing the ordinary
! differential equations (ODEs) which are calculated as the second
! step of the operational split method for the complete biogeochemical
! models. The remaining ODE is
! \begin{equation}\label{ODESystem}
!   \partial_t c_i 
! =   P_i(\vec{c}) -D_i(\vec{c}), \;\; i = 1,\ldots,I, 
! \end{equation}
! with $c_i$ denoting the concentrations of state variables.
! The right hand side denotes the reaction terms,
! which are composed of contributions
! $d_{i,j}(\vec{c})$, which represent reactive fluxes from
! $c_i$ to $c_j$, and in turn, $p_{i,j}(\vec{c})$ are reactive fluxes from
! $c_j$ received by $c_i$, see equation (\ref{eq:am:a}).
! 
! These methods are:
!
! \begin{enumerate}
! \item First-order explicit (not unconditionally positive)
! \item Second order explicit Runge-Kutta (not unconditionally positive)
! \item Fourth-order explicit Runge-Kutta (not unconditionally positive)
! \item First-order Patankar (not conservative)
! \item Second-order Patankar-Runge-Kutta (not conservative)
! \item Fourth-order Patankar-Runge-Kutta (does not work, not conservative)
! \item First-order Modified Patankar (conservative and positive)
! \item Second-order Modified Patankar-Runge-Kutta (conservative and positive)
! \item Fourth-order Modified Patankar-Runge-Kutta 
!       (does not work, conservative and positive)
! \item First-order Extended Modified Patankar 
!       (stoichiometrically conservative and positive)
! \item Second-order Extended Modified Patankar-Runge-Kutta 
!       (stoichiometrically conservative and positive)
! \end{enumerate}
!
! The schemes 1 - 5 and 7 - 8 have been described in detail by
! \cite{Burchardetal2003b}. Later, \cite{Bruggemanetal2005} could
! show that the Modified Patankar schemes 7 - 8 are only conservative
! for one limiting nutrient and therefore they developed the
! Extended Modified Patankar (EMP) schemes 10 and 11 which are also
! stoichiometrically conservative. Patankar and Modified Patankar
! schemes of fourth order have not yet been developed, such that 
! choices 6 and 9 do not work yet.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: solver,nlev,numc
   REALTYPE, intent(in)                :: dt
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                            :: dt_local
!EOP
!-----------------------------------------------------------------------
!BOC
LEVEL1 'ode_solver'

! Allocate global variables
   if (.NOT. allocated(cc0)) then
      allocate(cc0(0:nlev,1:numc))
      allocate(cc1(0:nlev,1:numc))
      allocate(rhs(0:nlev,1:numc))
      allocate(rhs1(0:nlev,1:numc))
      allocate(rhs2(0:nlev,1:numc))
      allocate(rhs3(0:nlev,1:numc))
      allocate(ppsum(0:nlev,1:numc))
      allocate(ddsum(0:nlev,1:numc))
#ifdef BFM_GOTM
!JM      allocate(ccb0(1:numbc,0:nlev))
!JM      allocate(ccb1(1:numbc,0:nlev))
!JM      allocate(ppbsum(1:numbc,0:nlev))
!JM      allocate(ddbsum(1:numbc,0:nlev))
      allocate(ccb0(0:nlev,1:numbc))
      allocate(ccb1(0:nlev,1:numbc))
      allocate(ppbsum(0:nlev,1:numbc))
      allocate(ddbsum(0:nlev,1:numbc))
      STDERR "bio_solver numc=",numc
      STDERR "bio_solver nlev=",nlev
      STDERR "bio_solver dt=",dt
      STDERR "bio_solver solver=",solver
#else
      allocate(pp(0:nlev,1:numc,1:numc))
      allocate(dd(0:nlev,1:numc,1:numc))
#endif
   end if
!#ifdef BFM_GOTM
!   dt_local=dt
!#else
   dt_local=dt/SEC_PER_DAY
!#endif
LEVEL1 'calling solver',solver
LEVEL1 'dt,SEC_PER_DAY',dt,SEC_PER_DAY
LEVEL1 'dt_local,numc,nlev',dt_local,numc,nlev

   select case (solver)
      case (1)
         call euler_forward(dt_local,numc,nlev,cc)
      case (2)
         call runge_kutta_2(dt_local,numc,nlev,cc)
      case (3)
         call runge_kutta_4(dt_local,numc,nlev,cc)
      case (4)
         call patankar(dt_local,numc,nlev,cc)
      case (5)
         call patankar_runge_kutta_2(dt_local,numc,nlev,cc)
      case (6)
         call patankar_runge_kutta_4(dt_local,numc,nlev,cc)
      case (7)
         call modified_patankar(dt_local,numc,nlev,cc)
      case (8)
         call modified_patankar_2(dt_local,numc,nlev,cc)
      case (9)
         call modified_patankar_4(dt_local,numc,nlev,cc)
      case (10)
         call emp_1(dt_local,numc,nlev,cc)
      case (11)
         call emp_2(dt_local,numc,nlev,cc)
      case default
         stop "bio: no valid solver method specified in bio.inp !"
   end select
!  call findnega(ccb(ppY2p,1),1,iout)
!  if ( iout.gt.0) then
!    STDERR  "solver 1:Negative value at end  "; STOP "in bio_solver"
!  endif

   return
   end subroutine ode_solver
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Euler-forward scheme
!
! !INTERFACE:
   subroutine euler_forward(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the first-order Euler-forward (E1) scheme is coded, with one 
! evaluation of the right-hand sides per time step:
! \begin{equation}\label{eq:am:euler}
! \begin{array}{rcl}
! c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
! - D_i\left(\underline{c}^n\right) \right\}.
! \end{array}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
  REALTYPE, intent(inout)              :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
! REALTYPE :: rhs
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
      ! reset diagonal terms only
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)
!  call findnega(ccb(ppY2p,1),1,iout)
!  if ( iout.gt.0) then
!    STDERR  "solver 1:Negative value after ode_solver "; STOP "in bio_solver"
!  endif


#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
     cc(:,:) = cc(:,:) + dt*sum(pp(:,:,:)-dd(:,:,:),1)
   end if
   ! integrate benthic variables 
!  rhs=ccb(ppY2p,1)
   if (bio_setup>1) then
      ccb(:,:) = ccb(:,:) + dt*sum(ppb(:,:,:)-ddb(:,:,:),1)
   end if
!  call findnega(ccb(ppY2p,1),1,iout)
!  if ( iout.gt.0) then
!    STDERR  "solver 2:Negative value after ode_solver "
!    STDERR  "prev value ccb(ppY2p,:)=",rhs
!    STDERR  "ccb(ppY2p,:)=",ccb(ppY2p,:)
!    STDERR  "sum(ppb(:,ppY2p,1)):",sum(ppb(:,ppY2p,1))*dt
!    STDERR  "sum(ddb(:,ppY2p,1)):",sum(ddb(:,ppY2p,1))*dt
!    STDERR  "ppb(ppY2p,:,1):",ppb(ppY2p,:,1)*dt
!    STDERR  "ppb(:,ppY2p,1)",ppb(:,ppY2p,1)*dt
!    STDERR  "ddb(ppY2p,:,1):",ddb(ppY2p,:,1)*dt
!    STDERR  "ddb(:,ppY2p,1)",ddb(:,ppY2p,1)*dt; STOP "in bio_solver"
!  endif
#else
   do ci=1,nlev
      do i=1,numc
         rhs=0.
         do j=1,numc
            rhs=rhs+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc(ci,i)=cc(ci,i)+dt*rhs
      end do
   end do
#endif

   return
   end subroutine euler_forward
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Runge-Kutta scheme
!
! !INTERFACE:
   subroutine runge_kutta_2(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the second-order Runge-Kutta (RK2) scheme is coded, with two
! evaluations of the right hand side per time step:
! \begin{equation}\label{eq:am:RK}
! \left.
! \begin{array}{rcl}
! c_i^{(1)} &=&
! \displaystyle
! c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
! - D_i\left(\underline{c}^n\right) \right), \\ \\
!  c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \dfrac{\Delta t}{2}
! \left\{P_i\left(\underline{c}^n\right) + P_i\left(\underline{c}^{(1)}\right)
!  - D_i\left(\underline{c}^n\right) - D_i\left(\underline{c}^{(1)}\right)
! \right\}.
! \end{array}
! \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE :: rhs(0:nlev,1:numc),rhs1(1:numc)
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     ! store the initial pelagic state
     cc0 = cc
     do ci=1,nlev
        do i=1,numc
           rhs(ci,i)=0.
           do j=1,numc
              rhs(ci,i)=rhs(ci,i)+pp(ci,i,j)-dd(ci,i,j)
           end do
           cc(ci,i)=cc0(ci,i)+dt*rhs(ci,i)
        end do
     end do
   end if ! bio_setup/=2

   if (bio_setup>1) then
     ! store the initial benthic state
     ccb0 = ccb
     do ci=1,NLEVB
        do i=1,numbc
           rhs(ci,i)=0.
           do j=1,numbc
              rhs(ci,i)=rhs(ci,i)+pp(ci,i,j)-dd(ci,i,j)
           end do
           ccb(ci,i)=ccb0(ci,i)+dt*rhs(ci,i)
        end do
     end do
   end if ! bio_setup>1
   
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif
   call process_model(first,numc,nlev,cc,pp,dd)

   if (bio_setup/=2) then
     do ci=1,nlev
        do i=1,numc
           rhs1(i)=0.
           do j=1,numc
              rhs1(i)=rhs1(i)+pp(ci,i,j)-dd(ci,i,j)
           end do
           cc(ci,i)=cc0(ci,i)+dt*0.5*(rhs(ci,i)+rhs1(i))
        end do
     end do
   end if ! bio_setup/=2

   if (bio_setup>1) then
     do ci=1,NLEVB
        do i=1,numbc
           rhs1(i)=0.
           do j=1,numbc
              rhs1(i)=rhs1(i)+ppb(ci,i,j)-ddb(ci,i,j)
           end do
           ccb(ci,i)=ccb0(ci,i)+dt*0.5*(rhs(ci,i)+rhs1(i))
        end do
     end do
   end if ! bio_setup>1

#else
   do ci=1,nlev
      do i=1,numc
         rhs(ci,i)=0.
         do j=1,numc
            rhs(ci,i)=rhs(ci,i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc1(ci,i)=cc(ci,i)+dt*rhs(ci,i)
      end do
   end do

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif
   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         rhs1(i)=0.
         do j=1,numc
            rhs1(i)=rhs1(i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc(ci,i)=cc(ci,i)+dt*0.5*(rhs(ci,i)+rhs1(i))
      end do
   end do
#endif

   return
   end subroutine runge_kutta_2
!EOC

!MAV NOT DONE
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Runge-Kutta scheme
!
! !INTERFACE:
   subroutine runge_kutta_4(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the fourth-order Runge-Kutta (RK4) scheme is coded, 
! with four evaluations
! of the right hand sides per time step:
! \begin{equation}\label{eq2}
! \left.
! \begin{array}{rcl}
! c_i^{(1)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\underline c^n\right)
! -D_i\left(\underline c^n\right)\right\} \\ \\
! c_i^{(2)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\frac12\left(\underline c^n+\underline 
! c^{(1)}\right)\right)-D_i\left(\frac12\left(\underline c^n+\underline 
! c^{(1)}\right)\right)\right\} \\ \\
! c_i^{(3)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\frac12\left(\underline c^n+\underline 
! c^{(2)}\right)\right)-D_i\left(\frac12\left(\underline c^n+\underline 
! c^{(2)}\right)\right)\right\} \\ \\
! c_i^{(4)} &=&
! \displaystyle
! c_i^n+\Delta t \left\{P_i\left(\underline c^{(3)}\right)-D_i\left(\underline 
! c^{(3)}\right)\right\} \\ \\
! c_i^{n+1} &=&
! \displaystyle
!  \frac16 \left\{c_i^{(1)}+2c_i^{(2)}+2c_i^{(3)}+c_i^{(4)}   \right\}.
! \end{array}
! \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE :: rhs(0:nlev,1:numc),rhs1(0:nlev,1:numc)
  REALTYPE :: rhs2(0:nlev,1:numc),rhs3(0:nlev,1:numc)
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd)

   do ci=1,nlev
      do i=1,numc
         rhs(ci,i)=0.
         do j=1,numc
            rhs(ci,i)=rhs(ci,i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc1(ci,i)=cc(ci,i)+dt*rhs(ci,i)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         rhs1(ci,i)=0.
         do j=1,numc
            rhs1(ci,i)=rhs1(ci,i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc1(ci,i)=cc(ci,i)+dt*rhs1(ci,i)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         rhs2(ci,i)=0.
         do j=1,numc
            rhs2(ci,i)=rhs2(ci,i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc1(ci,i)=cc(ci,i)+dt*rhs2(ci,i)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         rhs3(ci,i)=0.
         do j=1,numc
            rhs3(ci,i)=rhs3(ci,i)+pp(ci,i,j)-dd(ci,i,j)
         end do
         cc(ci,i)=cc(ci,i)+dt*1./3.*(0.5*rhs(ci,i)+rhs1(ci,i)+rhs2(ci,i)+0.5*rhs3(ci,i))
      end do
   end do

   return
   end subroutine runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Patankar scheme
!
! !INTERFACE:
   subroutine patankar(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the first-order Patankar-Euler scheme (PE1) scheme is coded,
! with one evaluation of the right hand sides per time step:
! \begin{equation}\label{eq:am:patankar}
! \begin{array}{rcl}
!   c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \Delta t \left\{P_i\left(\underline{c}^n\right)
!                                       - D_i\left(\underline{c}^n\right)
!                                         \frac{c_i^{n+1}}{c_i^n} \right\}.
! \end{array}
! \end{equation}
!
! !USES:
   use bio_bfm, only : test_on_negative_states

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
! REALTYPE :: ppsum ,ddsum
! integer  :: j
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     ! reset diagonal terms only
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     ! reset diagonal terms only
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif
   call process_model(first,numc,nlev,cc,pp,dd)


#ifdef BFM_GOTM
   if (bio_setup/=2) then
     cc = (cc+dt*sum(pp(:,1:numc,1:numc),1))/(1.+dt*sum(dd(:,1:numc,1:numc),1)/(1.0D-80+cc))
   end if
   ! compute benthic variables 
    if (bio_setup>1) then
      ccb = (ccb+dt*sum(ppb(:,:,:),1))/(1.+dt*sum(ddb(:,:,:),1)/(1.0D-80+ccb))
    end if
#else
   do ci=1,nlev
      do i=1,numc
         ppsum=0.
         ddsum=0.
         do j=1,numc
            ppsum=ppsum+pp(ci,i,j)
            ddsum=ddsum+dd(ci,i,j)
         end do
         cc(ci,i)=(cc(ci,i)+dt*ppsum)/(1.+dt*ddsum/cc(ci,i))
      end do
   end do
#endif

   return
   end subroutine patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine patankar_runge_kutta_2(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the second-order Patankar-Runge-Kutta (PRK2) scheme is coded,
! with two evaluations of the right hand sides per time step:
! 
! \begin{equation}\label{eq:am:PRK}
!   \left.
!   \begin{array}{rcl}
!     c_i^{(1)} &=&
! \displaystyle
!  c_i^n  + \Delta t
!                   \left\{P_i\left(\underline{c}^n\right)
!                       - D_i\left(\underline{c}^n\right)
!                         \dfrac{c_i^{(1)}}{c_i^n}\right\},
!                   \\ \\
!     c_i^{n+1} &=&
! \displaystyle
!  c_i^n  + \dfrac{\Delta t}{2}
!                   \left\{P_i\left(\underline{c}^n\right)
!                         + P_i\left(\underline{c}^{(1)}\right)
!                         - \left( D_i\left(\underline{c}^n\right)
!                         + D_i\left(\underline{c}^{(1)}\right)\right)
!                         \dfrac{c_i^{n+1}}{c_i^{(1)}}
!                   \right\}.
!   \end{array}
!   \right\}
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE :: ppsum(0:nlev,1:numc),ddsum(0:nlev,1:numc)
#ifdef BFM_GOTM
  REALTYPE :: ppbsum(0:nlev,1:numc),ddbsum(0:nlev,1:numc)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
        do i=1,numc
           ppsum(ci,i)=0.
           ddsum(ci,i)=0.
           do j=1,numc
              ppsum(ci,i)=ppsum(ci,i)+pp(ci,i,j)
              ddsum(ci,i)=ddsum(ci,i)+dd(ci,i,j)
           end do
           cc(ci,i)=(cc0(ci,i)+dt*ppsum(ci,i))/(1.+dt*ddsum(ci,i)/cc0(ci,i))
        end do
     end do
   end if
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,NLEVB
        do i=1,numbc
           ppbsum(ci,i)=0.
           ddbsum(ci,i)=0.
           do j=1,numbc
              ppbsum(ci,i)=ppbsum(ci,i)+ppb(ci,i,j)
              ddbsum(ci,i)=ddbsum(ci,i)+ddb(ci,i,j)
           end do
           ccb(ci,i)=(ccb0(ci,i)+dt*ppbsum(ci,i))/(1.+dt*ddbsum(ci,i)/ccb0(ci,i))
        end do
     end do
   end if

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

  ! call process_model(first,numc,nlev) ! TODO check if this line should exist
                                        !     and the correct arguments 

   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
        do i=1,numc
           do j=1,numc
              ppsum(ci,i)=ppsum(ci,i)+pp(ci,i,j)
              ddsum(ci,i)=ddsum(ci,i)+dd(ci,i,j)
           end do
           cc(ci,i)=(cc0(ci,i)+0.5*dt*ppsum(ci,i))/(1.+0.5*dt*ddsum(ci,i)/cc(ci,i))
        end do
     end do
   end if
   ! integrate benthic variables 
   if (bio_setup/=2) then
     do ci=1,NLEVB
        do i=1,numbc
           do j=1,numbc
              ppbsum(ci,i)=ppbsum(ci,i)+ppb(ci,i,j)
              ddbsum(ci,i)=ddbsum(ci,i)+ddb(ci,i,j)
           end do
           ccb(ci,i)=(ccb0(ci,i)+0.5*dt*ppbsum(ci,i))/(1.+0.5*dt*ddbsum(ci,i)/ccb(ci,i))
        end do
     end do
   end if

#else
   do ci=1,nlev
      do i=1,numc
         ppsum(ci,i)=0.
         ddsum(ci,i)=0.
         do j=1,numc
            ppsum(ci,i)=ppsum(ci,i)+pp(ci,i,j)
            ddsum(ci,i)=ddsum(ci,i)+dd(ci,i,j)
         end do
         cc1(ci,i)=(cc(ci,i)+dt*ppsum(ci,i))/(1.+dt*ddsum(ci,i)/cc(ci,i))
      end do
   end do

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         do j=1,numc
            ppsum(ci,i)=ppsum(ci,i)+pp(ci,i,j)
            ddsum(ci,i)=ddsum(ci,i)+dd(ci,i,j)
         end do
         cc(ci,i)=(cc(ci,i)+0.5*dt*ppsum(ci,i))/(1.+0.5*dt*ddsum(ci,i)/cc1(ci,i))
      end do
   end do
#endif

   return
   end subroutine patankar_runge_kutta_2
!EOC

!MAV NOT DONE
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine patankar_runge_kutta_4(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! This subroutine should become the fourth-order Patankar Runge-Kutta
! scheme, but it does not yet work.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE :: ppsum1(0:nlev,1:numc),ddsum1(0:nlev,1:numc)
   REALTYPE :: ppsum2(0:nlev,1:numc),ddsum2(0:nlev,1:numc)
   REALTYPE :: ppsum3(0:nlev,1:numc),ddsum3(0:nlev,1:numc)
#ifdef BFM_GOTM
!  REALTYPE :: ppbsum1(1:numbc,0:nlev) !,ddbsum1(1:numbc,0:nlev)
!  REALTYPE :: ppbsum2(1:numbc,0:nlev) !,ddbsum2(1:numbc,0:nlev)
!  REALTYPE :: ppbsum3(1:numbc,0:nlev) !,ddbsum3(1:numbc,0:nlev)
#endif

!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd)

   do ci=1,nlev
      do i=1,numc
         ppsum(ci,i)=0.
         ddsum(ci,i)=0.
         do j=1,numc
            ppsum(ci,i)=ppsum(ci,i)+pp(ci,i,j)
            ddsum(ci,i)=ddsum(ci,i)+dd(ci,i,j)
         end do
         cc1(ci,i)=(cc(ci,i)+dt*ppsum(ci,i))/(1.+dt*ddsum(ci,i)/cc(ci,i))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         ppsum1(ci,i)=0.
         ddsum1(ci,i)=0.
         do j=1,numc
            ppsum1(ci,i)=ppsum1(ci,i)+pp(ci,i,j)
            ddsum1(ci,i)=ddsum1(ci,i)+dd(ci,i,j)
         end do
         cc1(ci,i)=(cc(ci,i)+dt*ppsum1(ci,i))/(1.+dt*ddsum1(ci,i)/cc1(ci,i))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         ppsum2(ci,i)=0.
         ddsum2(ci,i)=0.
         do j=1,numc
            ppsum2(ci,i)=ppsum2(ci,i)+pp(ci,i,j)
            ddsum2(ci,i)=ddsum2(ci,i)+dd(ci,i,j)
         end do
         cc1(ci,i)=(cc(ci,i)+dt*ppsum2(ci,i))/(1.+dt*ddsum2(ci,i)/cc1(ci,i))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd)

   do ci=1,nlev
      do i=1,numc
         ppsum3(ci,i)=0.
         ddsum3(ci,i)=0.
         do j=1,numc
            ppsum3(ci,i)=ppsum3(ci,i)+pp(ci,i,j)
            ddsum3(ci,i)=ddsum3(ci,i)+dd(ci,i,j)
         end do
         ppsum(ci,i)=1./3.*(0.5*ppsum(ci,i)+ppsum1(ci,i)+ppsum2(ci,i)+0.5*ppsum3(ci,i))
         ddsum(ci,i)=1./3.*(0.5*ddsum(ci,i)+ddsum1(ci,i)+ddsum2(ci,i)+0.5*ddsum3(ci,i))
         cc(ci,i)=(cc(ci,i)+dt*ppsum(ci,i))/(1.+dt*ddsum(ci,i)/cc1(ci,i))
      end do
   end do

   return
   end subroutine patankar_runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Modified Patankar scheme
!
! !INTERFACE:
   subroutine modified_patankar(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the first-order Modified Patankar-Euler scheme (MPE1) scheme is coded,
! with one evaluation of the right hand side per time step:
! \begin{equation}\label{eq:am:MP}
! \begin{array}{rcl}
!   c_i^{n+1} &=&
! \displaystyle
!  c_i^n
!               + \Delta t \left\{ \sum\limits_{\stackrel{j=1}{j \not= i}}^I
!                p_{i,j}\left(\underline{c}^n\right) \dfrac{c_j^{n+1}}{c_j^n}
!                                         + p_{i,i}\left(\underline{c}^n\right)
!                           - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
!                                         \dfrac{c_i^{n+1}}{c_i^n} \right\}.
! \end{array}
! \end{equation}
! 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
   ! reset diagonal terms only
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     ! reset diagonal terms only
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
#endif
     do ci=1,nlev
        do i=1,numc
           a(i,i)=0.
           do j=1,numc
              a(i,i)=a(i,i)+dd(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/(1.0D-80+cc(j,ci))
           end do
           a(i,i)=dt*a(i,i)/(1.0D-80+cc(ci,i))
           a(i,i)=_ONE_+a(i,i)
           r(i)=cc(ci,i)+dt*pp(i,ci,i)
        end do
        call matrix(numc,a,r,cc(ci,:))
     end do
   end if

#ifdef BFM_GOTM
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,NLEVB
        do i=1,numbc
           a(i,i)=0.
           do j=1,numbc
              a(i,i)=a(i,i)+ddb(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*ppb(ci,i,j)/(1.0D-80+ccb(j,ci))
           end do
           a(i,i)=dt*a(i,i)/(1.0D-80+ccb(ci,i))
           a(i,i)=_ONE_+a(i,i)
           r(i)=ccb(ci,i)+dt*ppb(i,ci,i)
        end do
        call matrix(numbc,a,r,ccb(ci,:))
     end do
    end if
#endif

   return
   end subroutine modified_patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Modified Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine modified_patankar_2(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the second-order Modified Patankar-Runge-Kutta (MPRK2) scheme is coded,
! with two evaluations of the right hand sides per time step:
! 
! \begin{equation}\label{eq:am:MPRK}
!   \left. \begin{array}{rcl}
!     c_i^{(1)} &=&
! \displaystyle
! c_i^n  + \Delta t
! \left\{
! \sum\limits_{\stackrel{j=1}{j \not= i}}^I p_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_j^{(1)}}{c_j^n}
! + p_{i,i}\left(\underline{c}^n\right)
! - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_i^{(1)}}{c_i^n}
! \right\},
! \\ \\
! c_i^{n+1} &=&
! \displaystyle
! c_i^n  + \dfrac{\Delta t}{2}
!                   \left\{
!                     \sum\limits_{\stackrel{j=1}{j \not= i}}^I
!                       \left(p_{i,j}\left(\underline{c}^n\right)
!                           + p_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_j^{n+1}}{c_j^{(1)}}
!                       + p_{i,i}\left(\underline{c}^n\right)
!                       + p_{i,i}\left(\underline{c}^{(1)}\right)
! \right.\\ \\
!               & &
! \displaystyle
! \left.\phantom{c_i^n  + \dfrac{\Delta t}{2} }
!                   - \sum_{j=1}^I
!                       \left(d_{i,j}\left(\underline{c}^n\right)
!                           + d_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_i^{n+1}}{c_i^{(1)}}
!                   \right\}.
!   \end{array}
!   \right\}
! \end{equation}
! 
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE :: pp1(1:numc,0:nlev,1:numc),dd1(1:numc,0:nlev,1:numc)
#ifdef BFM_GOTM
   REALTYPE :: ppb1(1:numbc,1:numbc,0:nlev),ddb1(1:numbc,1:numbc,0:nlev)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif
   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
     cc0 = cc
     do ci=1,nlev
        do i=1,numc
           a(i,i)=_ZERO_
           do j=1,numc
              a(i,i)=a(i,i)+dd(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc0(j,ci)
           end do
           a(i,i)=dt*a(i,i)/cc0(ci,i)
           a(i,i)=1.+a(i,i)
           r(i)=cc0(ci,i)+dt*pp(i,ci,i)
        end do
        call matrix(numc,a,r,cc(ci,:))
     end do
   end if

   ! integrate benthic variables 
   if (bio_setup>1) then
     ccb0 = ccb
     do ci=1,NLEVB
        do i=1,numbc
           a(i,i)=_ZERO_
           do j=1,numbc
              a(i,i)=a(i,i)+ddb(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*ppb(ci,i,j)/ccb0(j,ci)
           end do
           a(i,i)=dt*a(i,i)/ccb0(ci,i)
           a(i,i)=1.+a(i,i)
           r(i)=ccb0(ci,i)+dt*ppb(i,ci,i)
        end do
        call matrix(numbc,a,r,ccb(ci,:))
     end do
   end if

#else
   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp(i,ci,i)
      end do
      call matrix(numc,a,r,cc1(ci,:))
   end do
#endif

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc1,pp1,dd1)

   pp=0.5*(pp+pp1)
   dd=0.5*(dd+dd1)
#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
        do i=1,numc
           a(i,i)=0.
           do j=1,numc
              a(i,i)=a(i,i)+dd(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc(j,ci)
           end do
           a(i,i)=dt*a(i,i)/cc(ci,i)
           a(i,i)=1.+a(i,i)
           r(i)=cc0(ci,i)+dt*pp(i,ci,i)
        end do
        call matrix(numc,a,r,cc0(ci,:))
     end do
    end if
   ppb=0.5*(ppb+ppb1)
   ddb=0.5*(ddb+ddb1)
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,NLEVB
        do i=1,numbc
           a(i,i)=0.
           do j=1,numbc
              a(i,i)=a(i,i)+ddb(ci,i,j)
              if (i.ne.j) a(i,j)=-dt*ppb(ci,i,j)/ccb(j,ci)
           end do
           a(i,i)=dt*a(i,i)/ccb(ci,i)
           a(i,i)=1.+a(i,i)
           r(i)=ccb0(ci,i)+dt*ppb(i,ci,i)
        end do
        call matrix(numbc,a,r,ccb0(ci,:))
     end do
   end if

#else
   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp(i,ci,i)
      end do
      call matrix(numc,a,r,cc(ci,:))
   end do
#endif

   return
   end subroutine modified_patankar_2
!EOC

!MAV NOT DONE
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Modified Patankar-Runge-Kutta scheme
!
! !INTERFACE:
   subroutine modified_patankar_4(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! This subroutine should become the fourth-order Modified Patankar Runge-Kutta
! scheme, but it does not yet work.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE :: pp1(1:numc,0:nlev,1:numc) ,dd1(1:numc,0:nlev,1:numc)
   REALTYPE :: pp2(1:numc,0:nlev,1:numc) ,dd2(1:numc,0:nlev,1:numc)
   REALTYPE :: pp3(1:numc,0:nlev,1:numc) ,dd3(1:numc,0:nlev,1:numc)
#ifdef BFM_GOTM
!  REALTYPE :: ppb1(1:numbc,1:numbc,0:nlev) !,ddb1(1:numbc,1:numbc,0:nlev)
!  REALTYPE :: ppb2(1:numbc,1:numbc,0:nlev) !,ddb2(1:numbc,1:numbc,0:nlev)
!  REALTYPE :: ppb3(1:numbc,1:numbc,0:nlev) !,ddb3(1:numbc,1:numbc,0:nlev)
#endif

!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp(i,ci,i)
      end do
      call matrix(numc,a,r,cc1(ci,:))
   end do

   call process_model(first,numc,nlev,cc1,pp1,dd1)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd1(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp1(ci,i,j)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp1(i,ci,i)
      end do
      call matrix(numc,a,r,cc1(ci,:))
   end do

   call process_model(first,numc,nlev,cc1,pp2,dd2)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd2(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp2(ci,i,j)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp2(i,ci,i)
      end do
      call matrix(numc,a,r,cc1(ci,:))
   end do

   call process_model(first,numc,nlev,cc1,pp3,dd3)

   pp=1./3.*(0.5*pp+pp1+pp2+0.5*pp3)
   dd=1./3.*(0.5*dd+dd1+dd2+0.5*dd3)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(ci,i,j)
            if (i.ne.j) a(i,j)=-dt*pp(ci,i,j)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(ci,i)
         a(i,i)=1.+a(i,i)
         r(i)=cc(ci,i)+dt*pp(i,ci,i)
      end do
      call matrix(numc,a,r,cc(ci,:))
   end do

   return
   end subroutine modified_patankar_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order Extended Modified Patankar scheme
!
! !INTERFACE:
   subroutine emp_1(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the first-order Extended Modified Patankar scheme for
! biogeochemical models is coded, with one evaluation of the right-hand
! side per time step:
!
! \begin{eqnarray}
!    \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{f}(t^n,\vec{c}^n)\prod\limits_{j \in J^n} {\frac{c_j^{n + 1} }{c_j^n}}} \nonumber \\
!    & & \mbox{with } J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\}
! \end{eqnarray}
!
! This system of non-linear implicit equations is solved in auxiliary subroutine
! findp\_bisection, using the fact this system can be reduced to a polynomial in one
! unknown, and additionally using the restrictions imposed by the requirement of positivity.
! For more details, see \cite{Bruggemanetal2005}.
!
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
  REALTYPE, intent(inout)              :: cc(0:nlev,1:numc)
!
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  REALTYPE :: pi, derivative(1:numc)
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
#endif
      do ci=1,nlev
         derivative(:) = sum(pp(ci,:,:),2)-sum(dd(ci,:,:),2)
         call findp_bisection(numc, cc(ci,:), derivative(:), dt, 1.d-9, pi)
         cc(ci,:) = cc(ci,:) + dt*derivative(:)*pi
      end do
#ifdef BFM_GOTM
   end if

   if (bio_setup>1) then
     do ci=1,NLEVB
       derivative(:) = sum(ppb(ci,:,:),2)-sum(ddb(ci,:,:),2)
       call findp_bisection(numbc, ccb(ci,:), derivative(:), dt, 1.d-9, pi)
       ccb(ci,:) = ccb(ci,:) + dt*derivative(:)*pi
     end do
   end if
#endif
 
   return
   end subroutine emp_1
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Extended Modified Patankar scheme
!
! !INTERFACE:
   subroutine emp_2(dt,numc,nlev,cc)
!
! !DESCRIPTION:
! Here, the second-order Extended Modified Patankar scheme for
! biogeochemical models is coded, with two evaluations of the right-hand
! side per time step:
!
! \begin{eqnarray}
!   \vec{c}^{(1)} & = & \vec{c}^n  + \Delta t \: \vec{f}(t^n ,\vec{c}^n )\prod\limits_{j \in J^n } {\frac{{c_j^{(1)} }}{{c_j^n }}} \nonumber \\
!   \vec{c}^{n + 1} & = & \vec{c}^n  + \frac{{\Delta t}}{2}\left( {\vec{f}(t^n ,\vec{c}^n ) + \vec{f}(t^{n + 1} ,\vec{c}^{(1)} )} \right)\prod\limits_{k \in K^n } {\frac{{c_k^{n + 1} }}{{c_k^{(1)} }}} 
! \end{eqnarray}
!
! where
!
! \begin{eqnarray}
!   J^n & = & \left\{ {i:f_i (t^n ,\vec{c}^n ) < 0, i \in \{ 1,...,I\} } \right\} \nonumber \\
!   K^n & = & \left\{ {i:f_i (t^n ,\vec{c}^n ) + f_i (t^{n+1} ,\vec{c}^{(1)} ) < 0, i \in \{ 1,...,I\} } \right\}.
! \end{eqnarray}
!
! The first step is identical to a step with the first-order EMP scheme. The second step mathmatically identical to
! a step with the first-order scheme if we rewrite it as
!
! \begin{eqnarray}
!   \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)})\prod\limits_{k \in 
!     K^n} {\frac{c_k^{n + 1} }{c_k^n }}} \nonumber \\
!   & & \mbox{with }\vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)}) = \frac{1}{2}\left( {\vec{f}(t^n,\vec{c}^n) + \vec{f}(t^{n + 1},\vec{c}^{(1)})} \right)\prod\limits_{k \in K^n} {\frac{c_k^n }{c_k^{(1)} }}.
! \end{eqnarray}
!
! Therefore, this scheme can be implemented as two consecutive steps with the first-order scheme, the second using
! $\vec{h}(t^n,t^{n + 1},\vec{c}^n,\vec{c}^{(1)})$. The non-linear problem of each consecutive step is solved
! in auxiliary subroutine findp\_bisection.
!
! For more details, see \cite{Bruggemanetal2005}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
!
! !INPUT/OUTPUT PARAMETER:
  REALTYPE, intent(inout)              :: cc(0:nlev,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  REALTYPE :: pi, rhs(0:nlev,1:numc), cc_med(0:nlev,1:numc)
#ifdef BFM_GOTM
  REALTYPE :: rhsb(0:nlev,1:numc), ccb_med(0:nlev,1:numc)
#endif

!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif
   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate pelagic variables 
   if (bio_setup/=2) then
#endif
   do ci=1,nlev
      rhs(ci,:) = sum(pp(ci,:,:),2) - sum(dd(ci,:,:),2)
      call findp_bisection(numc, cc(ci,:), rhs(ci,:), dt, 1.d-9, pi)
      cc_med(ci,:) = cc(ci,:) + dt*rhs(ci,:)*pi
   end do
#ifdef BFM_GOTM
   end if

   if (bio_setup>1) then
     do ci=1,NLEVB
      rhsb(ci,:) = sum(ppb(ci,:,:),2) - sum(ddb(ci,:,:),2)
      call findp_bisection(numbc, ccb(ci,:), rhsb(ci,:), dt, 1.d-9, pi)
      ccb_med(ci,:) = ccb(ci,:) + dt*rhsb(ci,:)*pi
     end do
   end if
#endif

#ifdef BFM_GOTM
   if (bio_setup/=2) then
     call reset_diagonal(numc,pp)
     call reset_diagonal(numc,dd)
   end if
   if (bio_setup>1) then
     call reset_diagonal(numbc,ppb)
     call reset_diagonal(numbc,ddb)
   endif
#endif

   call process_model(first,numc,nlev,cc,pp,dd)

#ifdef BFM_GOTM
   ! integrate 2nd round pelagic variables 
   if (bio_setup/=2) then
#endif
      do ci=1,nlev
         rhs(ci,:) = 0.5 * (rhs(ci,:) + sum(pp(ci,:,:),2) - sum(dd(ci,:,:),2))

         ! Correct for the state variables that will be included in 'p'.
         do i=1,numc
            if (rhs(ci,i) .lt. 0.) rhs(ci,:) = rhs(ci,:) * cc(ci,i)/cc_med(ci,i)
         end do

         call findp_bisection(numc, cc(ci,:), rhs(ci,:), dt, 1.d-9, pi)

         cc(ci,:) = cc(ci,:) + dt*rhs(ci,:)*pi
      end do ! ci (z-levels)
#ifdef BFM_GOTM
      ! reset diagonal terms only
   end if

   if (bio_setup>1) then
      do ci=1,NLEVB
         rhsb(ci,:) = 0.5 * (rhsb(ci,:) + sum(ppb(ci,:,:),2) - sum(ddb(ci,:,:),2))

         ! Correct for the state variables that will be included in 'p'.
         do i=1,numbc
            if (rhsb(ci,i) .lt. 0.) rhsb(ci,:) = rhsb(ci,:) * ccb(ci,i)/ccb_med(ci,i)
         end do

         call findp_bisection(numbc, ccb(ci,:), rhsb(ci,:), dt, 1.d-9, pi)

         ccb(ci,:) = ccb(ci,:) + dt*rhsb(ci,:)*pi
      end do ! ci (z-levels)
      ! reset diagonal terms only
   end if
#endif

   return
   end subroutine emp_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculation of the EMP product term 'p'
!
! !INTERFACE:
   subroutine findp_bisection(numc, cc, derivative, dt, accuracy, pi)
!
! !DESCRIPTION:
! Auxiliary subroutine for finding the Extended Modified Patankar 
! product term $p$ with the bisection technique.
!
! This subroutine solves the non-linear problem
!
! \begin{eqnarray}
!    \lefteqn{\vec{c}^{n + 1} = \vec{c}^n + \Delta t \: \vec{f}(t^n,\vec{c}^n)\prod\limits_{j \in J^n} {\frac{c_j^{n + 1} }{c_j^n}}} \nonumber \\
!    & & \mbox{with } J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\}
! \end{eqnarray}
!
! using the fact that it can be reduced to the problem of finding the root of a polynomial
! in one unknown $p: = \prod\limits_{j \in J^n}{{c_j^{n + 1}}/{c_j^n}}$:
!
! \begin{equation}
!   g(p) = \prod\limits_{j \in J^n} {\left( {1 + \frac{\Delta t \: f_j(t^n,\vec{c}^n)}{c_j^n }p} \right)} - p = 0,
! \end{equation}
!
! with
!
! \begin{equation}
!   J^n = \left\{ {i:f_i (t^n,\vec{c}^n) < 0,i \in \{1,...,I\}} \right\},
! \end{equation}
!
! Additionally, it makes use of the the positivity requirement $\vec{c}^{n+1}_i>0\ \forall\ i$, which
! imposes restriction
!
! \begin{equation}
!   p \in \left( 0, \min \left( {1,\mathop {\min }\limits_{j \in J^n} \left( { - 
!                   \frac{c_j^n }{\Delta t \: f_j (t^n,\vec{c}^n)}} \right)} \right) \right).
! \end{equation}
!
! It has been proved that there exists exactly one $p$ for which the above is 
! true, see \cite{Bruggemanetal2005}.
! The resulting problem is solved using the bisection scheme, which is guaranteed to converge.
!
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)   :: numc
   REALTYPE, intent(in)  :: cc(1:numc), derivative(1:numc)
   REALTYPE, intent(in)  :: dt, accuracy
!
! !OUTPUT PARAMETER:
   REALTYPE, intent(out) :: pi
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE :: pileft, piright, fnow
   REALTYPE :: relderivative(1:numc)
   integer  :: iter, i, potnegcount
!EOP
!-----------------------------------------------------------------------
!BOC
! Sort the supplied derivatives (find out which are negative).
   potnegcount = 0
   piright = 1.
   do i=1,numc

      if (derivative(i).lt.0.) then
!        State variable could become zero or less; include it in the
!        J set of the EMP scheme.
         if (cc(i).eq.0.) write (*,*) "Error: state variable ",i," is zero and has negative derivative!"
         potnegcount = potnegcount+1
         relderivative(potnegcount) = dt*derivative(i)/cc(i)

!        Derivative is negative, and therefore places an upper bound on pi.
         if (-1./relderivative(potnegcount).lt.piright) piright = -1./relderivative(potnegcount)
     end if

   end do

   if (potnegcount.eq.0) then
!     All derivatives are positive, just do Euler.
      pi = 1.0
      return
   end if

   pileft = 0.      ! polynomial(0) = 1

!  Determine maximum number of bisection iterations from
!  requested accuracy.
!  maxiter = -int(ceiling(dlog10(accuracy)/dlog10(2.D0)))

   do iter=1,20
!     New pi to test is middle of current pi-domain.
      pi = 0.5*(piright+pileft)

!     Calculate polynomial value.
      fnow = 1.
      do i=1,potnegcount
         fnow = fnow*(1.+relderivative(i)*pi)
      end do

      if (fnow>pi) then
!        Polynomial(pi)>0; we have a new left bound for pi.
         pileft = pi
      elseif (fnow<pi) then
!       Polynomial(pi)<0; we have a new right bound for pi.
        piright = pi
      else
!       Freak occurrence: polynomial(pi)=0, we happened to pinpoint
!       the exact pi.
        exit
      end if
!     Check if we now pi accurately enough (accuracy refers to the
!     number of decimals we know).
      if ((piright-pileft)/pi<accuracy) exit
   end do

!  Low pi values imply very large negative relative derivative. This happens
!  for stiff systems (or very high delta_t), and for non-positive systems.
!  Then EMP is not suitable (it will stall state variable values), so warn user.
   if (pi.lt.1.d-4) then
     write (*,*) "Warning: small pi=",pi," in Extended Modified Patankar slows down system!"
!    write (*,*) "relative derivatives: ",derivative(:)*dt/cc(:)
!    write (*,*) "You system may be stiff or non-positive, or you time step is too large."
!    stop "ode_solvers::findp_bisection"
   end if

   return

   end subroutine findp_bisection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE:
   subroutine matrix(n,a,r,c)
!
! !DESCRIPTION:
! This is a Gaussian solver for multi-dimensional linear equations.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! INPUT/OUTPUT PARAMETERS:
  REALTYPE,dimension(:,:),intent(INOUT)              :: a
  REALTYPE,dimension(:),intent(INOUT)                :: r
!
! OUTPUT PARAMETERS:
  REALTYPE,dimension(:),intent(OUT)                :: c
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
! Alternative using vectors:
!  do i=1,n
!     a(i,i+1:n)=a(i,i+1:n)/a(i,i)
!     a(i,i)=1.0
!     do k=i+1,n
!        r(k)=r(k)-a(k,i)*r(i)
!        a(k,i+1:n)=a(k,i+1:n)-a(k,i)*a(i,i+1:n)
!     end do
!  end do
!  do i=n,1,-1
!     c(i)=r(i)-sum(a(i,i+1:n)*c(i+1:n)) 
!  end do

   do i=1,n
      r(i)=r(i)/a(i,i)
      do j=n,i,-1
         a(i,j)=a(i,j)/(1.0D-80+a(i,i))
      end do
      do k=i+1,n
         r(k)=r(k)-a(k,i)*r(i)
         do j=i+1,n
            a(k,j)=a(k,j)-a(k,i)*a(i,j)
         end do
      end do
   end do

   do i=n,1,-1
      c(i)=r(i) 
      do j=i+1,n
         c(i)=c(i)-a(i,j)*c(j)
      end do
   end do

   return
   end subroutine matrix
!EOC

  end module bio_solver

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team 
! under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
