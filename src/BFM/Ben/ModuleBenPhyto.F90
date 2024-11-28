#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   Parameter values for the phytoplankton groups
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenPhyto
!
! !USES:

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO, &
          LOGUNIT,NML_OPEN,NML_READ,NMLUNIT,error_msg_prn
  use mem,ONLY: iiBenPhyto,iiBP1,NO_BOXES_XY,CoupledtoBDc
  use mem_Param,only:ChlLightFlag

!  
!
! !AUTHORS
!   the ERSEM group, Marcello Vichi, JWB, HBB
!
!
!
! !REVISION_HISTORY
!   !
!
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenPhyto PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !
  !  ---------------- Physiological parameters -----------------
  !
  !
  !  ---------------- Nutrient parameters in phytoplankton -----------------
  !
  integer     :: p_useparams(iiBenPhyto)  ! excretion as sugars or as TEP
  real(RLEN)  :: p_xresus(iiBenPhyto)  ! part of stress output as ecreted c  

  real(RLEN)  :: p_xladm(iiBenPhyto)  ! distance wich diatom can move a day
                                      ! and extend the light layer.
  real(RLEN)  :: p_hBPc(iiBenPhyto) !mmax biomass per m3 porewater	
                 !diatoms form mats .at low densites of Y2c Y2 is not cpable to to brak this mats
                 ! only at high densities  when holes are for surface to deep the surface for diatoms
                 ! is limited and cannot for large mats.
  real(RLEN)  :: p_px_silt_in_pw(iiBenPhyto) !fraction of silt in sediment present inporewater.
  real  :: p_mxSiltBpc ! at high benphyto  conc %silt in upper layer is lowered
                       !with a monod function  which this parmeter
  integer     :: sw_precision ! limitation of diffsuion of nutients in dense
                                          !diatomlayer
  !
  !  ------------- Chlorophyll parameters -----------
  !  skel: Skeletonema costatum pav: Pavlova lutheri
  !  syn: Synechoccus sp. (significant alpha decrease with irradiance)
  !  gyr: Gyrodinium sp. iso: Isochrysis galbana
  !              skel     iso      syn      gyr


  public InitBenPhyto,CalculateBenPhyto,CalculateBenPhyto_vector
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenPhyto()

  implicit none
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  namelist /BenPhyto_parameters/p_useparams,p_xresus, &
                  p_xladm,p_hBPc,p_px_silt_in_pw,p_mxSiltBpc,sw_precision
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   if ( ChlLightFlag.ne.2) then
     write(LOGUNIT,*)'BenPhyto cannot be used if ChlLightFlag<>2 ( See Param.nml)!'
     goto 101
   endif
   CoupledtoBDc=0
   p_px_silt_in_pw=0.1
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenPhyto parameters.."
   open(NMLUNIT,file='BenPhyto.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=BenPhyto_parameters,err=101)
   close(NMLUNIT)
   write(LOGUNIT,*) "#  Namelist is:"
   write(LOGUNIT,nml=BenPhyto_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenPhyto.f90","BenPhyto.nml")
101 call error_msg_prn(NML_READ,"InitBenPhyto.f90","BenPhyto_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBenPhyto

  function CalculateBenPhyto_vector(constituent,mode0,from,xto,inverse,param)
  use mem,ONLY: Dfm,Dcm,iiC,iiL
  use constants, ONLY: INTEGRAL,AVERAGE,EQUATION

  implicit none
  integer,intent(IN)                           :: constituent
  integer,intent(IN)                           :: mode0
  real(RLEN),dimension(NO_BOXES_XY)            :: CalculateBenPhyto_vector
  real(RLEN),dimension(NO_BOXES_XY),intent(IN) :: from
  real(RLEN),dimension(NO_BOXES_XY),intent(IN) :: xto
  real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional :: inverse
  real(RLEN),intent(IN),optional :: param

  integer                                      :: mode
  real(RLEN),dimension(NO_BOXES_XY)            ::  r,d
  real(RLEN),dimension(NO_BOXES_XY)            ::  s
  real(RLEN),parameter                         ::  pi=3.14159265D+00
  real(RLEN),parameter                         ::  two=2.0D+00
  ! mode=EQUATION:
  ! mode=INTEGRAL integral  from->xto
  ! mode= AVERAGE: modus depth between frsom->xto
  ! mode= -INTEGRAL DEPTH WHERE INTEGRAL =xto
  ! mode=-PARAMETER calc Dfm where INTEGRAL between 0 and xto ==from
! write(LOGUNIT,*) "CalculateBenPhyto_vector"
  mode=mode0
  if (present(inverse).and.mode0.eq.INTEGRAL)mode=-mode0
  if (.not.present(param)) then
     select case (constituent)
       case(iiC); d=Dfm+NZERO
       case(iiL); d=Dcm+NZERO
     end select
  else
     d=param
  endif
   select case (mode)
      case (EQUATION)
         r=d**2
         CalculateBenPhyto_vector=r/(r+from**2)
      case (INTEGRAL,AVERAGE)
         s=ZERO
         r=max(ZERO,min(DONE,atan(xto/d)*two/pi))
         where (from.gt.ZERO) s=max(ZERO,min(DONE,atan(from/d)*two/pi))
         select case (mode)
            case (INTEGRAL) ; CalculateBenPhyto_vector=r-s 
            case (AVERAGE)  
               CalculateBenPhyto_vector= &
                min((from+xto)/two,max(from,abs(d*tan((s+r)*pi))))
         end select
     case (-INTEGRAL)
           CalculateBenPhyto_vector =min(xto,max(from,d*tan(inverse*pi/two)))
      case default ; write(LOGUNIT,*) 'calculateBenPhyto_vector mode=',mode
                     stop 'stopped'
    end select
  end function CalculateBenPhyto_vector

  function CalculateBenPhyto(constituent,mode0,nr,from,xto,inverse,param)
  use mem,ONLY: Dfm,Dcm,iiC,iiL
  use constants, ONLY: INTEGRAL,EQUATION,AVERAGE

  implicit none
  integer,intent(IN)                           :: constituent
  integer,intent(IN)                           :: mode0
  integer,intent(IN)                           :: nr
  real(RLEN)                     :: CalculateBenPhyto
  real(RLEN),intent(IN)          :: from
  real(RLEN),intent(IN)          :: xto
  real(RLEN),intent(IN),optional :: inverse
  real(RLEN),intent(IN),optional :: param

  integer               :: mode
  real(RLEN)            ::  r
  real(RLEN)            ::  s,d
  real(RLEN),parameter  ::  pi=3.14159265D+00
  real(RLEN),parameter  ::  two=2.0D+00

  mode=mode0
  if (present(inverse).and.mode0.eq.INTEGRAL)mode=-mode0
  if (.not.present(param)) then
     select case (constituent)
       case(iiC); d=Dfm(nr)+NZERO
       case(iiL); d=Dcm(nr)+NZERO
     end select
  else
     d=param
  endif
   select case (mode)
     case (EQUATION)
       s=d**2
       CalculateBenPhyto=s/(s+from**2)
     case (INTEGRAL,AVERAGE)
       s=ZERO
       r=max(ZERO,min(DONE,atan(xto/d)*two/pi))
       if (from.gt.ZERO) s=max(ZERO,min(DONE,atan(from/d)*two/pi))
       select case (mode)
         case (INTEGRAL) ; CalculateBenPhyto=r-s 
         case (AVERAGE)  ; CalculateBenPhyto= &
                min((from+xto)/two,max(from,abs(d*tan((s+r)*pi))))
       end select
     case (-INTEGRAL)
           CalculateBenPhyto=min(xto,max(from,d*tan(inverse*pi/two)))
     case default ; write(LOGUNIT,*) 'CalculateBenPhyto mode=',mode0,mode
                    stop 'stopped'

   end select
   end function CalculateBenPhyto

  subroutine CalculateAverageNutrientConcInLightLayer
  !This routine assume that there is only one type of BenthicPhytoplankton
  !present. (type seq.nummor Benphyt== iiBP1)

  use mem,ONLY: BoxNumberXY,D1m,Dlm,ShiftDlm,iiBen,&
           ppKp4n,ppKn4n,ppKp3n,ppKp1p,ppKp5s,ppQpun,ETW_Ben, &
           Kp4n,Kn4n,Kp3n,Kp1p,Kp5s,Qpun,  K3n,K4n,K1p,K5s, &
           KNH4,KNO3,KPO4,KSiO3,KQun, flux,iiL
  use mem_BenSilica, ONLY: p_chM5s, p_cvM5s,p_q10
  use mem_Param,  ONLY: p_poro,p_pK5_ae,p_d_tot
  use Constants, only:INTEGRAL,MASS
  use bennut_interface,   ONLY: CalculateFromSet
  use global_interface,   ONLY: eTq
! use mem_globalfun,   ONLY: insw
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  implicit none
! integer               ::  bp=1
  integer               ::  i
  real(RLEN)            ::  r1,r3,r4,r5,rkn,t1,t3,t4,t5,tkn,tkp
  real(RLEN)            ::  s,st
  real(RLEN)            ::  p,d0,dc,step
  real(RLEN)            ::  chM5s,Dlm_new
  real(RLEN),external   ::  GetDelta

! nrphyto=p_useparams(bp)
  do BoxNumberXY=1,NO_BOXES_XY
     chM5s = (p_chM5s+ &
         p_cvM5s*( eTq( ETW_Ben(BoxNumberXY), p_q10)- DONE)) &
                   *p_poro(BoxNumberXY)*(DONE+p_pK5_ae(BoxNumberXY))
     Dlm_new=Dlm(BoxNumberXY)+ShiftDlm(BoxNumberXY)*GetDelta()
     if ( sw_precision ==0 ) then
       t3=max(ZERO,CalculateFromSet(KNO3(BoxNumberXY),INTEGRAL, &
                                                  MASS,ZERO,Dlm_new))
       t4=max(ZERO,CalculateFromSet(KNH4(BoxNumberXY),INTEGRAL, &
                                                  MASS,ZERO,Dlm_new))
       t1=max(ZERO,CalculateFromSet(KPO4(BoxNumberXY),INTEGRAL, &
                                                  MASS,ZERO,Dlm_new))
       t5=CalculateFromSet(KSiO3(BoxNumberXY),INTEGRAL,MASS,ZERO,Dlm_new)
       t5=chM5s*Dlm_new-t5
       tkn=max(ZERO,CalculateFromSet(KQun(BoxNumberXY),INTEGRAL, &
                                                  MASS,ZERO,Dlm_new))
     else 
       t3=ZERO;t4=ZERO;t1=ZERO;t5=ZERO;p=ZERO;tkn=ZERO;tkp=ZERO
       d0=ZERO
       st=CalculateBenPhyto(iiL,INTEGRAL,BoxNumberXY,ZERO,Dlm_new)
!      write(logunit,*) 'CANCILL st',st,Dlm,Dlm_new,ShiftDlm(BoxNumberXY)
       step=st/sw_precision
       d0=ZERO
       do i=1,sw_precision
         p=p+step
         dc=CalculateBenPhyto(iiL,INTEGRAL,BoxNumberXY,ZERO,Dlm_new,inverse=p)
!        write(logunit,*) 'CANCILL d1',d1
         r3=CalculateFromSet(KNO3(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         r4=CalculateFromSet(KNH4(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         r1=CalculateFromSet(KPO4(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         r5=CalculateFromSet(KSiO3(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         s=CalculateBenPhyto(iiL,INTEGRAL,BoxNumberXY, d0,dc)
!        rkc=CalculateFromSet(KQ1c(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         rkn=CalculateFromSet(KQun(BoxNumberXY),INTEGRAL,MASS,d0,dc)
         t3=t3+r3*s
         t4=t4+r4*s
         t1=t1+r1*s
         t5=t5+(chM5s*(dc-d0)-r5)*s
!        tkc=tkc+rkc*s
         tkn=tkn+rkn*s
         d0=dc
       enddo
       t3=max(ZERO,t3/st)
       t4=max(ZERO,t4/st)
       t1=max(ZERO,t1/st)
       t5=max(ZERO,t5/st)
!      tkc=max(ZERO,tkc/st)
       tkn=max(ZERO,tkn/st)
     endif
!    write(LOGUNIT,*) ' CalculateAverageNutrientConcInLightLayer'
     if (isnan(t4)) then
       write(LOGUNIT,*) 'MBP t1 is NaN', t3,t4,t1,t5,Dlm_new,s,st
       write(LOGUNIT,*) 'MBP t3 is NaN', K3n,K4n,K1p,K5s,Dlm_new,s,st
     endif
     !Here state variable are only used to store information for the next time 
     !step when these concentrations are used as an estimation of the average 
     !external nutrient concentration used to estimate the nutrient uptake in 
     !benthic diatoms in benthic diatoms.See BFM/Ben/limitNutrientUptake.F90 
     if ( isnan(Kp3n(boxNumberXY)))Kp3n(boxNumberXY)=ZERO
     if ( isnan(Kp4n(boxNumberXY)))Kp4n(boxNumberXY)=ZERO
     if ( isnan(Kp1p(boxNumberXY)))Kp1p(boxNumberXY)=ZERO
     if ( isnan(Kp5s(boxNumberXY)))Kp5s(boxNumberXY)=ZERO
     if ( isnan(Qpun(boxNumberXY)))Qpun(boxNumberXY)=ZERO
     s=max(ZERO,CalculateFromSet(KNH4(BoxNumberXY),INTEGRAL, &
                                                MASS,Dlm_new,D1m(BoxNumberXY)))
     call flux(BoxNumberXY,iiBen,ppKp3n,ppKp3n,( &
                                            t3-Kp3n(BoxNumberXY))/GetDelta())
     call flux(BoxNumberXY,iiBen,ppKp4n,ppKp4n, &
                                            (t4-Kp4n(BoxNumberXY))/GetDelta())
     call flux(BoxNumberXY,iiBen,ppKn4n,ppKn4n, &
                                            ( s-Kn4n(BoxNumberXY))/GetDelta())
     call flux(BoxNumberXY,iiBen,ppKp1p,ppKp1p,( &
                                            t1-Kp1p(BoxNumberXY))/GetDelta())
     call flux(BoxNumberXY,iiBen,ppQpun,ppQpun,( &
                                            tkn-Qpun(BoxNumberXY))/GetDelta())

     if ( D1m(BoxNumberXY) > Dlm_new) then
       r5=CalculateFromSet(KSiO3(BoxNumberXY),INTEGRAL,MASS,ZERO, &
                                   D1m(BoxNumberXY))
!      Output2d_1(BoxNumberXY)=r5
       r5=(chM5s*D1m(BoxNumberXY)-r5)
!      Output2d_2(BoxNumberXY)=r5
!      Output2d_3(BoxNumberXY)=chm5s*D1m(BoxNumberXY)
       t5=max(ZERO,min(t5,r5))
     endif
     call flux(BoxNumberXY,iiBen,ppKp5s,ppKp5s, &
                                   (t5-Kp5s(BoxNumberXY))/GetDelta())
   enddo
  end subroutine CalculateAverageNutrientConcInLightLayer
  
   subroutine test_BenPhyto_status(mode,out)

   use mem,ONLY: CoupledtoBDc,iiPhytoPlankton,BenPhyto,iiC
   use mem_Param,ONLY:CalcPelagicFlag,CalcPhytoPlankton

      implicit none
      integer,intent(IN)           :: mode
      integer,intent(OUT),optional :: out
      real(RLEN)                   :: scalar2
      integer                      :: nrphyto,i

      if (.not.present(out)) then
        do i=1,iiBenPhyto
          nrphyto=p_useparams(i)
          scalar2=maxval(BenPhyto(i,iiC))
          CoupledtoBDc(nrphyto)=0
          if (p_xresus(i)>ZERO .and. CalcPhytoPlankton(nrphyto) &
                           .and.CalcPelagicFlag) then
             CoupledtoBDc(nrphyto)=i
          ! if there is no coupled if Benphyto to a pelagic var
          ! CoupledtoBDc(nrphyto)=0
          ! if Benphyto-biomass is very small the number phytoplankton
          ! types are added  to CoupledtoBdc . This means:
          ! a. NO benphytodynamics is calculated.
          ! b. However sedimentation of resuspended diatoms is still possible.
          !    (See PelBen/Settling.F90)
             if (scalar2 < 1.0D-5) &
                 CoupledtoBDc(nrphyto)=CoupledtoBDc(nrphyto)+iiPhytoPlankton
          endif
        enddo
      elseif (mode.gt.0) then
        out=-1
        nrphyto=p_useparams(mode)
        if ( nrphyto.gt.0) then
          out=1
          if( CoupledtoBDc(nrphyto)>iiPhytoPlankton) out=0
        endif
      else
        out=0
        do i=1,iiBenPhyto
          nrphyto=p_useparams(i)
          if ( nrphyto.gt.0) then
            if (out.eq.0.and.CoupledtoBDc(nrphyto)<=iiPhytoPlankton) out=1
          endif
        enddo
      endif
   end subroutine test_BenPhyto_status

  end module mem_BenPhyto
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
