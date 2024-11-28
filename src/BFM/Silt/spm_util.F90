!$Id$
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spm_util --- suspended matter model \label{sec:spm}
!
! !INTERFACE:
   module spm_util
!
! !DESCRIPTION:
!
! Contains more or less general purpose routines for SPM calculations
!-----------------------------------------------------------------------
!
! !USES:
!   use spm_var
!
!  default: all is public.
   public
!
! !PUBLIC MEMBER FUNCTIONS:
!
! !PRIVATE MEMBER FUNCTIONS:
!   private 
!
! !REVISION HISTORY:
!  Original author(s): Johan Van Der Molen & Karsten Bolding
!
! Nov 2013: taken from  gotm-3.3.2_withspm_may2011, made independent of spm_var 
!           and simplified to apply to a single size fraction
!
!  $Log$
!
! !LOCAL VARIABLES:
   REALTYPE, parameter       :: g=9.81, pi=3.141592654 

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculates reference concentrations
!
! !INTERFACE:
   subroutine reference_concentration(cref,z_a,refconc_mode,tau_mode,tauc_mode, &
   &                                  tauw_mode,taucrit_mode,tau_w,tau_c,tau_bed,nu,phi_wc,u_orb, &
   &                                  u_star,rho,rho_s,grainsize,u,v,h1,Tz,eta)
!
! !DESCRIPTION:
!  Calls routine to calculate combined bed-shear stress.
!
!  Calculates critical shields number following Soulsby, 1997 p. 106
!
!   theta_crit=0.3/(1.0+1.2*dstar)+0.055*(1.0-exp(-0.02*dstar))
!
!  Calculates critical bed-shear stress
!
!     tau_crit=g*(rho_s-rho)*grainsize*theta_crit
!
!  Calculates Non-dimensional bed-shear stress:
!
!    T_a=max(0, (tau_bed-taucrit_mode*tau_crit)/tau_crit )
!
!  taucrit_mode is a user-set switch to include (1) or exclude (0) a 
!  threshold for motion.
!
!  Calculates near-bed reference concentrations for all size fractions
!  1 option:
!  refconc_mode=1: Smith & McLean (1977) method, see Soulsby (1997), p. 140.
!                    -calculates dimensionless bed-shear stress T_a
!                     (subroutine) (also critical shear stress tau_crit)
!                    -calculates reference level for coarsest fraction:
!                      z_a=26.3*tau_crit*T_a/(rho*g*(rho_s/rho-1)) 
!                          +grainsize/12.
!                    -calculates reference concentration
!                      cref=rho_s*c2*T_a/(1+c3*T_a) 
!                    -limit total concentration to 65% (solid sand)
!
!  Refs:
!  Smith, J.D., McLean, S.R. (1977). Spatially averaged flow over a wavy
!     surface. Journal of Geophysical Research 82, 1735-1743.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)     :: refconc_mode,tau_mode,tauc_mode,tauw_mode,taucrit_mode

   REALTYPE, intent(out)   :: tau_w,tau_c,tau_bed
   REALTYPE, intent(in)    :: rho,rho_s,grainsize,nu,u,v,h1,Tz,u_orb
   REALTYPE, intent(out)   :: cref,z_a,u_star,eta
!
! !REVISION HISTORY:
!  Original author(s): Johan van der Molen
!  Spring 2008: JM: made changes to ensure that all variables (e.g. reference height) 
!                   are calculated per size fraction instead of referring to 
!                   the coarsest fraction
!
! !LOCAL VARIABLES
   integer              :: n
   REALTYPE, parameter :: c1=0.015, c2=0.00156, c3=0.0024
   REALTYPE            :: dstar,theta_crit,tau_crit,T_a,phi_wc,z0
!EOP
!-----------------------------------------------------------------------
!BOC

   select case (refconc_mode)
   case (1)   ! Smith & McLean 1977, Soulsby 1997, p. 140
    ! calculate dstar
       dstar=( (g*(rho_s/rho-1)/(nu**2.0))**(1.0/3.0) )*grainsize
    ! calculate critical shear stress for initiation of motion
       theta_crit=0.3/(1.0+1.2*dstar)+0.055*(1.0-exp(-0.02*dstar))
       tau_crit=g*(rho_s-rho)*grainsize*theta_crit
    ! calculate combined bed-shear stress
       call combined_bed_shear_stress(tau_w,tau_c,tau_bed,tau_mode,tauc_mode,tauw_mode, &
       &                              rho,rho_s,grainsize,theta_crit,dstar,z0,phi_wc, &
       &                              u_orb,u_star,u,v,h1,Tz,eta &
       &                              )
    ! dimensionless bedshear stress, no pickup if tau_bed<tau_crit
       T_a=max(_ZERO_, (tau_bed-taucrit_mode*tau_crit)/tau_crit ) 
    ! calculate reference level
       z_a=26.3*tau_crit*T_a/(rho*g*(rho_s/rho-1)) &
       &     +z0
    ! calculate reference concentration
       cref=rho_s*c2*T_a/(1+c3*T_a)

   case default
     stop
!JM error message here
   end select

   ! limit to 65% volume concentration, which is the limit of 'solid' sediment
   cref=min(cref,0.65*rho_s)

   return
   end subroutine reference_concentration
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculates wave parameters 
!
! !INTERFACE:
   subroutine wave_parameters(Hs,Tz,u,v,phi_w,phi_wc,u_orb,depth)
!
! !DESCRIPTION:
! Calculates wave parameters: wave-current angle, wave number, orbital 
!  velocity, orbital excursion.
! Currently set up to ignore current effect on wave number; code present
! to include that effect.
!
! Subroutine calls to calculate angle between waves and currents.
! Function call to solve dispersion relation and calculate wave number k.
! Orbital velocity:
!
! u_orb=abs(pi*Hs/(T*sinh(2*pi*depth/lambda)))
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
  REALTYPE, INTENT(IN)  :: Hs, Tz, u, v, phi_w, depth
  REALTYPE, INTENT(OUT) :: u_orb,phi_wc

! !REVISION HISTORY:
!  Original author(s): Johan van der Molen
!
! !LOCAL VARIABLES
   REALTYPE :: omega,k,lwm,tp,a_delta
   REALTYPE, parameter :: pi=3.141592654, d360=360.0
   REALTYPE :: phi_c
!EOP
!-----------------------------------------------------------------------
!BOC
    omega=2*pi/Tz

    !calculate wave-current angle
    phi_c=atan2(v,u)              ! wave angle wrt x-axis -pi<phi_c<pi
    if (phi_c .lt. _ZERO_) then       ! degrees wrt North
       phi_c=180.*(pi/2-phi_c)/pi
    else
       phi_c=180.*(2*pi-(phi_c-pi/2))/pi
    end if
    phi_c=mod(phi_c,d360)             ! reduce to <360
    phi_wc=abs(phi_c-phi_w)          ! angle between waves and currents
    if (phi_wc>180.) phi_wc=360.-phi_wc           ! reduce to <180
    phi_wc=phi_wc*pi/180.

   ! calculate wave number, wave length and period
    k=dispersion_relation(omega,depth,_ZERO_) 
    lwm=2*pi/k
    tp=Tz
! the alternative would be to include effects current on wave number as below
! here vr is depth-averaged current
!    k=dispersion_relation(omega,depth,vr*cos(phi_wc))    
!    lwm=2*pi/k
!    tp=Tz/(1-vr*Tz*cos(phi_wc)/lwm)

   ! calculate near-bed orbital excursion
    a_delta=abs(Hs/(2*sinh(2*pi*depth/lwm)))

   ! calculate near-bed peak orbital velocity
    u_orb=abs(pi*Hs/(tp*sinh(2*pi*depth/lwm)))

   return
   end subroutine wave_parameters
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculates combined bed-shear stress by waves and currents
!
! !INTERFACE:
   subroutine combined_bed_shear_stress(tau_wt,tau_ct,tau_bedt,tau_mode,tauc_mode,tauw_mode, &
   &                                    rho,rho_s,grainsize,theta_crit,dstar,z0,phi_wc, &
   &                                    u_orb,u_star,u,v,h1,Tz,eta &
   &                                    )
!
! !DESCRIPTION:
! Calculates skin-friction bed-shear stress by currents
!
! 1 option: 
! tauc_mode=1:   law of wall
!     Skin-friction roughness length for coarsest size class:
!
!        z0=grainsize/12.
!  
!     Skin-friction shear velocity from first grid cell above bottom:
!
!       u_star=sqrt(u1**2+v1**2)*kappa/log((h1/2)/z0)
!
!     Current-friction factor:
!
!        tau_c=rho*(u_star**2)
!
!  Calculates skin-friction bed-shear stress  component due to waves.
!  1 option:
!  tauw_mode=1:  Soulsby, 1997, p. 78.
!     Skin-friction roughness length for coarsest size class:
!
!        z0=grainsize/12.
!  
!     Wave-friction factor:
!
!        fw=c1*(u_orb*wper/(z0*2*pi))**c2
!
!     Wave bed-shear stress:
!
!        tauw=0.5*rho_w*fw*(u_orb**2)
!
! Calculates combined bed-shear stress by waves and currents
! 2 options:
! tau_mode=1: vector addition of current and wave components
!
!          tau_bed=sqrt(tau_c**2+tau_w**2+2*tau_c*tau_w*cos(phi_wc))
!
! tau_mode=2: Soulsby's method (Soulsby, 1997, p. 92)
!
!          tau_wc=tau_c*(1.0+1.2*(tau_w/(tau_c+tau_w))**3.2)
!          tau_bed=sqrt((tau_wc+tau_w*cos(phi_wc))**2+(tau_w*sin(phi_wc))**2)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS
   INTEGER, INTENT(IN)     :: tau_mode,tauc_mode,tauw_mode
   REALTYPE, INTENT(IN)    :: grainsize,phi_wc,rho,rho_s,theta_crit,dstar,u_orb
   REALTYPE, INTENT(IN)    :: u,v,h1,Tz
   REALTYPE, intent(out)   :: tau_wt,tau_ct,tau_bedt,eta
   REALTYPE, intent(out)   :: z0,u_star

! !OUTPUT PARAMETERS

! ! LOCAL PARAMETERS

!
! !REVISION HISTORY:
!  Original author(s): Johan van der Molen
!  Spring 2008: JM: changed z0 to refer to average bed composition instead of coarsest
!                   fraction
!
! !LOCAL PARAMETERS
   REALTYPE, parameter  :: c1=1.39, c2=-0.52, eps=0.0001, a_r=0.045 !0.923
   REALTYPE, parameter  :: c3=5.0, c4=30.0, g=9.81
   REALTYPE             :: omega        ! radial frequency
   REALTYPE             :: z0f, z0t, z0tot, &
   &                       lambda, u_start, fwt, eta_w, lambda_w, eta_c, &
   &                       fw,lambda_c,tau_wc,tau_w,tau_c,tau_bed,tau_wct
   REALTYPE             :: eta_washout=0.15
   INTEGER, parameter   :: niter=100
   LOGICAL              :: washout

!EOP
!-----------------------------------------------------------------------
!BOC
   ! skin friction roughness length
   z0=grainsize/12.

   ! current skin-friction shear stress
   select case(tauc_mode)
     case(1)   ! law of wall
       u_star=sqrt(u**2+v**2)*0.4/log((h1/2)/z0)
       tau_c=rho*(u_star**2)            
     case default
!JM error message here
       stop
   end select

   ! wave skin-friction shear stress
   select case(tauw_mode)
     case(1)         ! Soulsby
       if (u_orb.gt.eps) then
          fw=c1*(u_orb*Tz/(z0*2*pi))**c2     ! Soulsby, p.78
          tau_w=0.5*rho*fw*(u_orb**2)
       else
          tau_w=_ZERO_
       endif
     case default
       stop
!JM error message here
   end select

   ! combined skin-friction shear stress
   select case (tau_mode)
      case (1)   ! vector addition
         tau_bed=sqrt(tau_c**2+tau_w**2+2*tau_c*tau_w*cos(phi_wc))
      case (2)   ! Soulsby (1995,1997)
         tau_wc=tau_c*(1.0+1.2*(tau_w/(max(tau_c+tau_w,1E-10)))**3.2)
         tau_bed=sqrt((tau_wc+tau_w*cos(phi_wc))**2+(tau_w*sin(phi_wc))**2)
  end select

   ! calculate bed forms based on skin friction shear stresses
   call do_bedforms(eta_w,lambda_w,eta_c,lambda_c,rho,rho_s,grainsize, &
   &                      theta_crit,dstar,tau_bed,tau_w,tau_c,u_orb,Tz,washout)

   ! form-related bed-shear stress
   if (eta_w.gt.eta_c) then
     eta=eta_w
     lambda=lambda_w
   else
     eta=eta_c
     lambda=lambda_c
   endif
   if (eta.gt.eps .and. lambda.gt.eps) then
     z0f=a_r*(eta**2)/lambda
   else
     if (washout) eta=eta_washout
     z0f=_ZERO_
   endif

   ! transport-related bed-shear stress (Wilson, 1989)
   z0t=c3*tau_wc/(c4*g*(rho_s-rho))
   ! total bed shear stress
   z0tot=z0+z0f+z0t

   ! current total shear stress
   select case(tauc_mode)
     case(1)   ! law of wall
       u_start=sqrt(u**2+v**2)*0.4/log((h1/2)/z0tot)
       tau_ct=rho*(u_start**2)            
     case default
!JM error message here
       stop
   end select

   ! wave total shear stress
   select case(tauw_mode)
     case(1)         ! Soulsby
       if (u_orb.gt.eps) then
          fwt=c1*(u_orb*Tz/(z0tot*2*pi))**c2     ! Soulsby, p.78
          tau_wt=0.5*rho*fwt*(u_orb**2)
       else
          tau_wt=_ZERO_
       endif
     case default
       stop
!JM error message here
   end select

   ! combined total bed-shear stress
   select case (tau_mode)
      case (1)   ! vector addition
         tau_bedt=sqrt(tau_ct**2+tau_wt**2+2*tau_ct*tau_wt*cos(phi_wc))
      case (2)   ! Soulsby (1995,1997)
         tau_wct=tau_ct*(1.0+1.2*(tau_wt/(max(tau_ct+tau_wt,1E-10)))**3.2)
         tau_bedt=sqrt((tau_wct+tau_wt*cos(phi_wc))**2+(tau_wt*sin(phi_wc))**2)
   end select

   return
   end subroutine combined_bed_shear_stress
!EOC

!-----------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Calculate vertical sediment diffusion coefficient
!!
!! !INTERFACE:
!   subroutine sediment_diffusivity(sd_mode)
!!
!! !DESCRIPTION:
!! Calculates vertical sediment diffusion coefficient
!! 2 options:
!! sd_mode=1: linear
!!       K_s=kappa*u_star*z_a
!! sd_mode=2: parabolic
!!       K_s=kappa*u_star*z_a*(1-z_a/depth)
!!
!! Calculates u_star from combined bed-shear stress.
!!
!! !USES:
!   IMPLICIT NONE
!!
!! !INPUT PARAMETERS
!   INTEGER, INTENT(IN) :: sd_mode
!!
!! !REVISION HISTORY:
!!  Original author(s): Johan van der Molen
!!
!! !LOCAL PARAMETERS
!   REALTYPE, PARAMETER :: kappa=0.4
!!EOP
!!-----------------------------------------------------------------------
!!BOC
!   u_startot=sqrt(tau_bed/rho(1))
!
!   select case(sd_mode)
!     case(1)   ! linear
!       K_s=kappa*u_startot*z_a
!     case(2)   ! parabolic
!       K_s=kappa*u_startot*z_a*(1-z_a/depth)
!   end select
!
!   return
!   end subroutine sediment_diffusivity
!!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Solves dispersion relation using Newton Raphson method
!
! !INTERFACE:
   function dispersion_relation(omega,h,vcosphi)
!
! !DESCRIPTION:
! Solves dispersion relation iteratively using Newton Raphson method
! Capable of handling modification of wave number by presence of currents
! by setting vcosphi=v*cos(phi) with v current magnitude and phi 
! wave-current angle. vcosphi=0 ignores the wave-current effect.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS
   REALTYPE, INTENT(IN) :: omega,h,vcosphi

! !OUTPUT PARAMETERS
   REALTYPE :: dispersion_relation
!
! !REVISION HISTORY:
!  Original author(s): Johan van der Molen
!  Oct. 2006: added modification by currents
!
! !LOCAL PARAMETERS
   REALTYPE, parameter  :: niter=20, tol=0.001
   REALTYPE             :: k=0.01,eps,ff,dff,kplus,sigma
   integer              :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   n=0
   eps=10*tol

   do
      if (n.ge.niter .or. eps.lt.tol) EXIT
      n=n+1
      sigma=omega-k*vcosphi                       ! modified frequency
      ff=g*k*tanh(k*h)-sigma**2                   ! dispersion relation
      dff=g*tanh(k*h)+g*k*h*(1-(tanh(k*h))**2)+ & ! derivative
      &    2.0*vcosphi*(omega-k*vcosphi)
      kplus=k-ff/dff                              ! new estimate
      eps=abs(kplus-k)                            ! improvement
      k=kplus
   end do

   dispersion_relation=k

   return
   end function dispersion_relation
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate bed form dimensions
!
! !INTERFACE:
   subroutine do_bedforms(eta_w,lambda_w,eta_c,lambda_c,rho,rho_s,grainsize, &
   &                      theta_crit,dstar,tau_bed,tau_w,tau_c,u_orb,Tz,washout)
!
! !DESCRIPTION:
! Calculates bed form dimensions according to Aldridge et al. (in press)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS
!   INTEGER, INTENT(IN) :: 
  REALTYPE, INTENT(IN)    :: rho,rho_s,grainsize,theta_crit,dstar,tau_bed,tau_w, &
  &                          tau_c,u_orb,Tz
  REALTYPE, INTENT(INOUT) :: eta_w,lambda_w,eta_c,lambda_c
  LOGICAL, INTENT(OUT)    :: washout
!
! !REVISION HISTORY:
!  Original author(s): Johan van der Molen
!
! !LOCAL PARAMETERS
   INTEGER, PARAMETER  :: wave_bedform_method=1
   REALTYPE, PARAMETER :: theta_washout=0.83, thousand=1000.0, seven=7.0
   REALTYPE, PARAMETER :: max_silt_size=59.0E-6, gz=9.81
   REALTYPE, PARAMETER :: c1=1.8, c2=1.5, c3=4.0, c4=0.6
   REALTYPE, PARAMETER :: c5=0.22, c6=-0.16, c7=0.16, c8=-0.04
   REALTYPE, PARAMETER :: c9=0.48, c10=1.5, c11=4.0, c12=0.8, c13=-1.5
   REALTYPE, PARAMETER :: c14=0.28, c15=1.5, c16=4.0, c17=0.6, c18=-1.0
   REALTYPE, PARAMETER :: a3d=0.22E-3, a2d=0.3E-3, fh3d=0.55, fl3d=0.73
   REALTYPE, PARAMETER :: c_19=0.275, c_20=0.022, c_21=0.42, c_22=1.97, &
   &                      c_23=0.44, c_24=0.21, c_25=0.001, c_26=2.5
   REALTYPE            :: theta_wc, theta_w, theta_c, theta_b, A_w
   REALTYPE            :: term1, term2
   REALTYPE            :: a2d3d,coreta3d,corlambda3d,etastar,lambdastar,rmob, &
   &                      gs, u_orb10, u_rms, eta_wo, lambda_wo
!EOP
!-----------------------------------------------------------------------
!BOC
   eta_c=_ZERO_
   eta_w=_ZERO_
   lambda_c=_ZERO_
   lambda_w=_ZERO_

   ! if the dominant class is not a silt fraction, calculate ripple height
   if (grainsize.gt.max_silt_size) then
     theta_wc=tau_bed/(g*(rho_s-rho)*grainsize)
     theta_w=tau_w/(g*(rho_s-rho)*grainsize)
     theta_c=tau_c/(g*(rho_s-rho)*grainsize)

     A_w=u_orb*Tz/(2*pi)  ! orbital excursion

     washout=theta_wc.gt.theta_washout   ! sheet flow conditions?

     if ((.not.washout) .and. theta_wc.gt.theta_crit) then

       ! current ripples
       if (theta_c.gt.theta_crit) then
         lambda_c=thousand*grainsize
         eta_c=thousand*grainsize/seven/2.0  !JM division by two gives better results West Gabbard
       endif

       ! wave ripples
       if (theta_w.gt.theta_crit) then
         select case(wave_bedform_method)
           case (1)             ! Grant & Madsen, 1982

             theta_b=c1*theta_crit*((dstar**c2)/c3)**c4
             if (theta_w.lt.theta_b) then
               eta_w=c5*A_w*(theta_w/(theta_crit))**c6
               lambda_w=eta_w/(c7*(theta_w/theta_crit)**c8)
             else
               term1=(dstar**c10)/c11
               term2=theta_w/theta_crit
               eta_w=c9*a_w*(term1**c12)*(term2**c13)
               term1=((dstar**c15)/c16)**c17
               term2=(theta_w/theta_crit)**c18
               lambda_w=eta_w/(c14*term1*term2)
             endif
             ! post-storm Large Wave Ripples not implemented, may be added here

           case default
! error message here
             stop
         end select

       endif ! wave
!       write(*,*)"No sheet flow!"
     else
!       write(*,*)"Sheet flow!"
!       stop
     endif   ! washout
   endif     ! silt

   return
   end subroutine do_bedforms
!EOC

!-----------------------------------------------------------------------

   end module spm_util

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
