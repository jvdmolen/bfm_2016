#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcWaterProp
!
! DESCRIPTION
!

!
! !INTERFACE
        subroutine CalcWaterProp(n,Salt,T,Rho,Mu)
        use global_mem,ONLY:RLEN,DONE

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
!
  integer,intent(IN)                              ::n
  real(RLEN),dimension(n),intent(IN)              ::Salt
  real(RLEN),dimension(n),intent(IN)              ::T
  real(RLEN),dimension(n),intent(out),optional    ::Rho
  real(RLEN),dimension(n),intent(out),optional    ::Mu

  real(RLEN),dimension(5)                         ::d_a
  real(RLEN),dimension(5)                         ::d_b
  real(RLEN),dimension(10)                        ::v_a
  real(RLEN),dimension(n)                         ::Rho_w
  real(RLEN),dimension(n)                         ::D_Rho
  real(RLEN),dimension(n)                         ::S
  real(RLEN),dimension(n)                         ::A
  real(RLEN),dimension(n)                         ::B
  integer                                         ::check


  check=0
  S=Salt/1.0D+03
  if (present(Rho)) then
    check=1
    d_a(1) = 9.9992293295D+02
    d_a(2) = 2.0341179217D-02
    d_a(3) =-6.1624591598D-03
    d_a(4) = 2.2614664708D-05
    d_a(5) =-4.6570659168D-08
    d_b(1) = 8.0200240891D+02
    d_b(2) =-2.0005183488D+00
    d_b(3) = 1.6771024982D-02
    d_b(4) =-3.0600536746D-05
    d_b(5) =-1.6132224742D-05

!orginal matlab:
!  rho_w = a(1) + a(2)*T + a(3)*T.^2 + a(4)*T.^3 + a(5)*T.^4;
!  D_rho = b(1)*s + b(2)*s.*T + b(3)*s.*T.^2 + b(4)*s.*T.^3 + b(5)*s^2.*T.^2;

    rho_w= d_a(1) + T*(d_a(2) + T*(d_a(3) +T* d_a(4) + T* d_a(5)));
    D_rho= S*(d_b(1) +T*(d_b(2)*T +T*(d_b(3) +T*d_b(4))) +S*T*T*d_b(5))
    Rho   = rho_w + D_rho;
  endif
  if (present(Mu)) then
     check=1
     v_a(1)= 1.5700386464D-01
     v_a(2)= 6.4992620050D+01
     v_a(3)=-9.1296496657D+01
     v_a(4)= 4.2844324477D-05
     v_a(5)= 1.5409136040D+00
     v_a(6)= 1.9981117208D-02
     v_a(7)=-9.5203865864D-05
     v_a(8)= 7.9739318223D+00
     v_a(9)=-7.5614568881D-02
     v_a(10)=4.7237011074D-04

! Original matlab equations:
! mu_w = a(4) + 1./(a(1)*(T+a(2)).^2+a(3));
! A  = a(5) + a(6) * T + a(7) * T.^2;
! B  = a(8) + a(9) * T + a(10)* T.^2;
! mu = mu_w.*(1 + A.*S + B.*S.^2);

   A  = v_a(5) + T*(v_a(6) + T*v_a(7))
   B  = v_a(8) + T*(v_a(9) + T*v_a(10))
   Mu=(v_a(4) + DONE/(v_a(1)*(T+v_a(2))**2+v_a(3)))*(DONE+S*(A+S*B));
  endif
  if ( check==0) then
     stop "error no known optional output parameter given"
  endif
  end
