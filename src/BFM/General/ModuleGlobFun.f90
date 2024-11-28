!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleGlobFun
!
! DESCRIPTION
!   List of general model functions

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE mem_globalfun
!
! !USES:
  USE global_mem, ONLY:RLEN, ZERO, DONE,BASETEMP

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   --------
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  ! SHARED GLOBAL FUNCTIONS (must be below contains)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  contains

FUNCTION EXP_LIMIT_SCALAR(input)
  use constants, ONLY: DBL_MAX, DBL_MIN, MAX_VAL_EXPFUN, MIN_VAL_EXPFUN

        real(RLEN),intent(IN) ::input
        real(RLEN) ::EXP_LIMIT_scalar


        integer  :: check

        check=0
        if (abs(input) > MAX_VAL_EXPFUN ) check=1 
        if (abs(input) < MIN_VAL_EXPFUN ) check=-1 ;
        ! EXP_LIMIT function is defined between 1.0e+10*DBL_MIN to 1.0e-10 DBL_MAX
        if (check == 0) then
            EXP_LIMIT_SCALAR=exp(input)
        elseif(input>0.0D+00) then
            EXP_LIMIT_SCALAR=0.5D+00*(-DBL_MIN*1.0D+10*real(check-1) &
                                  +DBL_MAX*1.0D-10*real(check+1))
        else
            EXP_LIMIT_SCALAR=0.5D+00*(-DBL_MAX*1.0D-10* real(check-1) &
                                  +DBL_MIN*1.0D+10*real(check+1))
        endif 
        end function EXP_LIMIT_SCALAR
FUNCTION EXP_LIMIT(input)
  use constants, ONLY: DBL_MAX, DBL_MIN, MAX_VAL_EXPFUN, MIN_VAL_EXPFUN

        real(RLEN),intent(IN) ::input(:)
        real(RLEN) ::EXP_LIMIT(size(input))

!       real(RLEN) :: check(size(input))

        integer  :: check(size(input))

        check=0
        where (abs(input) > MAX_VAL_EXPFUN ) check=1 
        where (abs(input) < MIN_VAL_EXPFUN ) check=-1 ;
        where (check == 0)
            EXP_LIMIT=exp(input)
        ! EXP_LIMIT function is defined between 1.0e+10*DBL_MIN to 1.0e-10 DBL_MAX
        elsewhere (input>0.0D+00)
            EXP_LIMIT=0.5D+00*(-DBL_MIN*1.0D+10*real(check-1) &
                           +DBL_MAX*1.0D-10*real(check+1))
        elsewhere 
            EXP_LIMIT=0.5D+00*(-DBL_MAX*1.0D-10*real(check-1) &
                           +DBL_MIN*1.0D+10 *real(check+1))
        endwhere

        end function EXP_LIMIT
FUNCTION INSW_VECTOR(input)
        real(RLEN),intent(IN) ::input(:)
        real(RLEN) ::INSW_VECTOR(size(input))

        INSW_VECTOR =ZERO
        where (input > ZERO ) INSW_VECTOR=DONE 

        end function INSW_VECTOR
FUNCTION MM_VECTOR(vector,param)
        real(RLEN),intent(IN) ::param
        real(RLEN)            ::vector(:)
        real(RLEN)            ::MM_VECTOR(size(vector))

        MM_VECTOR= VECTOR / ( VECTOR+  PARAM)

        end function MM_VECTOR
FUNCTION MM_POWER_VECTOR(vector,param,pow)
        real(RLEN),intent(IN) ::param
        integer               ::pow
        real(RLEN)            ::vector(:)
        real(RLEN)            ::MM_POWER_VECTOR(size(vector))

        MM_POWER_VECTOR= VECTOR**pow / ( VECTOR**pow+  PARAM**pow)

        end function MM_POWER_VECTOR
FUNCTION eramp_VECTOR(x, m)
        real(RLEN),intent(IN) ::x(:)
        real(RLEN),intent(IN) ::m
        real(RLEN)            ::eramp_VECTOR(size(x))

        eramp_VECTOR =ZERO
        where (x > 0 ) 
           eramp_VECTOR=DONE 
           where (X< M) eramp_VECTOR=X/M;
        endwhere

        end function

FUNCTION PartQ(p, d_a, d_b, d_m)
        implicit none
        real(RLEN),intent(IN) ::p
        real(RLEN),intent(IN) ::d_a
        real(RLEN),intent(IN) ::d_b
        real(RLEN),intent(IN) ::d_m
        real(RLEN)            ::PartQ

        REAl(RLEN)  ::alpha
        REAl(RLEN)  ::c1
        REAL(RLEN)  ::b1 
        REAL(RLEN)  ::a1 
        REAL(RLEN)  ::norm 
        REAL(RLEN)  ::r 

        alpha=SetExpDist(p)
        c1 = min(abs(p) * (-DONE * log(1.0D-20)), d_m);
        b1 = min(d_b, c1);
        a1 = min(d_a, b1);
        r=ZERO
        if ( d_a ==ZERO ) r=DONE

        if (c1 > ZERO .and. p /= ZERO) then
           norm=IntegralExpDist(-alpha,c1)
           PartQ=(IntegralExpDist(-alpha,b1)-IntegralExpDist(-alpha,a1))/norm
        else 
           PartQ =r
        endif

        end function

FUNCTION PartQ_vector(p, d_a, d_b, d_m)
        real(RLEN),intent(IN) ::p(:)
        real(RLEN),intent(IN) ::d_a(:)
        real(RLEN),intent(IN) ::d_b(:)
        real(RLEN),intent(IN) ::d_m
        real(RLEN)            ::PartQ_VECTOR(size(p))

        real(RLEN),dimension(:),allocatable ::alpha
        real(RLEN),dimension(:),allocatable ::c1
        REAL(RLEN),dimension(:),allocatable ::b1 
        REAL(RLEN),dimension(:),allocatable ::a1 
        REAL(RLEN),dimension(:),allocatable ::norm 
        REAL(RLEN),dimension(:),allocatable ::r 

        ALLOCATE(alpha(size(p)))
        ALLOCATE(c1(size(p)))
        ALLOCATE(b1(size(p)))
        ALLOCATE(a1(size(p)))
        ALLOCATE(norm(size(p)))
        ALLOCATE(r(size(p)))

        alpha=SetExpDist_vector(p);
        c1 = min(abs(p) * (-DONE * log(1.0D-20)), d_m);
        b1 = min(d_b, c1);
        a1 = min(d_a, b1);
        r=ZERO
        where ( d_a ==ZERO ) r=DONE

        where (c1 > ZERO .and. p /= ZERO)
           norm=IntegralExpDist_vector(-alpha,c1)
           PartQ_VECTOR=(IntegralExpDist_vector(-alpha,b1)- &
                      IntegralExpDist_vector(-alpha,a1))/norm
        elsewhere 
           PartQ_VECTOR =r
        endwhere

        DEALLOCATE(c1,b1,a1,norm,r)
        end function
function eTq_VECTOR(temp,p_q10)

        IMPLICIT NONE
        real(RLEN)            :: temp(:)
        real(RLEN),intent(IN) :: p_q10
        real(RLEN)            ::eTq_VECTOR(size(temp))

        if ( abs(p_q10-DONE) < 1.0D-6) then
          eTq_VECTOR=  DONE
        else
          eTq_VECTOR=  exp(  log(  p_q10)*( temp- BASETEMP)/ 10.0D+00)
        endif
        end function eTq_VECTOR

function IntegralExp_vector(alfa,x)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: alfa(:)
        real(RLEN),intent(IN)       :: x(:)
        real(RLEN)                 :: IntegralExp_vector(size(alfa))
        
          IntegralExp_vector=(exp(alfa * x) -DONE)/alfa
        end function IntegralExp_vector

function IntegralExpDist_vector(alfa,x)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: alfa(:)
        real(RLEN),intent(IN)       :: x(:)
        real(RLEN)                 :: IntegralExpDist_vector(size(alfa))
        
        where (alfa.lt.ZERO)
          IntegralExpDist_vector=(exp_limit(alfa * x) -DONE)/alfa
        elsewhere
          IntegralExpDist_vector=(exp_limit(alfa*x)-alfa*x-DONE)/alfa
        endwhere
        end function IntegralExpDist_vector
function IntegralExp(alfa,x)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: alfa
        real(RLEN),intent(IN)       :: x
        real(RLEN)                 :: IntegralExp
        
          IntegralExp=(exp_limit_scalar(alfa * x) -DONE)/alfa
        end function IntegralExp
function IntegralExpDist(alfa,x)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: alfa
        real(RLEN),intent(IN)       :: x
        real(RLEN)                 :: IntegralExpDist
        
        if ( alfa .lt.ZERO) then
          IntegralExpDist=(exp_limit_scalar(alfa * x) -DONE)/alfa
        else
          IntegralExpDist=(exp_limit_scalar(alfa*x)-alfa*x-DONE)/alfa
        endif
        end function IntegralExpDist

function SetExpDist(Dm,Dlm,Dmm)
      IMPLICIT none
      real(RLEN),intent(IN)          :: Dm
      real(RLEN),intent(IN),optional :: Dlm
      real(RLEN),intent(IN),optional :: Dmm
      real(RLEN)                     :: SetExpDist

      real(RLEN)                     :: r

      r=Dm
      if (present(Dlm)) r=sign(max(abs(r),Dlm),Dm)
      if (present(Dmm)) r=sign(min(abs(r),Dmm),Dm)
      SetExpDist=DONE/r

      return
      end function SetExpDist

function SetExpDist_vector(Dm,Dlm,Dmm)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: Dm(:)
        real(RLEN),intent(IN),optional :: Dlm
        real(RLEN),intent(IN),optional :: Dmm
        real(RLEN)                  :: SetExpDist_vector(size(Dm))

        real(RLEN)                     :: r(size(Dm))

        r=Dm
        if (present(Dlm)) r=sign(max(abs(r),Dlm),Dm)
        if (present(Dmm)) r=sign(min(abs(r),Dmm),Dm)
        SetExpDist_vector=DONE/r
      return
      end function SetExpDist_vector

function ExpDist(alpha,Dm)
      IMPLICIT none
      real(RLEN),intent(IN)          :: alpha
      real(RLEN),intent(IN)          :: Dm
      real(RLEN)                     :: ExpDist

      if (alpha.lt.ZERO) then
          ExpDist=exp_limit_scalar(alpha*Dm)
      else
          ExpDist=exp_limit_scalar(alpha*Dm)-DONE
      endif
      end function ExpDist

FUNCTION INSW(input)
        real(RLEN),intent(IN) ::input
        real(RLEN) ::INSW

        INSW =ZERO
        if (input > ZERO ) INSW=DONE 

        end function INSW
!#RTSAFE.FOR 
      function rtsafe(funcd,x1,x2,xacc,xout)
      IMPLICIT NONE
      real(RLEN),intent(IN) ::X1
      real(RLEN),intent(IN) ::X2
      real(RLEN),intent(IN) ::XACC
      real(RLEN),intent(OUT)::XOUT
      real(RLEN)            ::rtsafe 


      INTERFACE                           ! Specification
        SUBROUTINE FUNCD(X,F,DF )         ! Specification
          USE global_mem, ONLY:RLEN
          REAL(RLEN),INTENT(IN)   ::X     ! Specification
          REAL(RLEN),INTENT(OUT)  ::F     ! Specification
          REAL(RLEN),INTENT(OUT)  ::DF    ! Specification
        END SUBROUTINE FUNCD                ! Specification
      END INTERFACE                       ! Specification    
      

      real(RLEN),PARAMETER        ::MAXIT=100
      real(RLEN)                  :: F,FL,DF,FH
      real(RLEN)                  :: XL,XH,TEMP
      real(RLEN)                  :: DXOLD,DX
      integer                     :: j
      logical                     :: ready

      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL*FH.GE.0.) then
        rtsafe=1;return
      ELSEIF(FL.LT.0.)THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        TEMP=FL
        FL=FH
        FH=TEMP
      ENDIF
      XOUT=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(XOUT,F,DF)
!     DO J=1,MAXIT
      j=0;
      ready=.FALSE.
      do while (.not.ready .and. J<MAXIT)
        J=J+1
        IF(((XOUT-XH)*DF-F)*((XOUT-XL)*DF-F).GE.0. &
            .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          XOUT=XL+DX
          ready=(XL.EQ.XOUT)
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=XOUT
          XOUT=XOUT-DX
          ready= (TEMP.EQ.XOUT)
        ENDIF
        ready=(ABS(DX).LT.XACC) 
        if ( .not.ready) then
          CALL FUNCD(XOUT,F,DF)
          IF(F.LT.0.) THEN
            XL=XOUT
            FL=F
          ELSE
            XH=XOUT
            FH=F
          ENDIF
        endif
     enddo
     if ( j.ge.MAXIT) then
       rtsafe=2;RETURN
     else
       rtsafe=0;return
     ENDIF
     end function rtsafe

  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
