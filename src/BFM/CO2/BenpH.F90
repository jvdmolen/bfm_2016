#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenpH
!
! DESCRIPTION
!   Calculate pH in the oxic and in the denitrification+anoxic layer
!

! !INTERFACE
  subroutine BenpHDynamics
!

#ifdef INCLUDE_BENCO2

! !USES:


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
  use mem,  ONLY: G23h,G13h, G3h, G23c,G13c, G3c,O3c, D1m, D2m,D2STATE
  use mem, ONLY: PelBoxAbove,BoxNumber,BoxNumberXY, &
     NO_BOXES_XY,ppG23h,ppG13h,ppG3h, &
     DICae,  pHae, pCO2ae, DICan,  pHan,pHdn, pCO2an,  ETW, LocalDelta, &
     HCO3ae, CO3ae, CO2ae,HCO3an, CO3an, CO2an,CAcae,CAcan, iiBen,&
    ESW, ERHO, K1p, K5s,K15s,Acae, Acan,K11p,K21p,D1m,D2m,flux,O3h_Ben
  USE BFM_ERROR_MSG, ONLY: BFM_ERROR,set_warning_for_getm

  use CO2System,ONLY: CalcCO2System,HplusBASIS
  use mem_Param,  ONLY: p_d_tot,p_poro,p_d_tot_2,p_clD1D2m,p_pK1_ae

  use mem_BenthicNutrient3, ONLY: p_pAn2Ni
  use mem_BenSilica, ONLY: p_clK5D1m => p_clD1m,p_clK5D2m => p_clD2m, &
        p_chK5D2m => p_chD2m,p_pK5_ae=>p_p_ae,p_pK5_an=>p_p_an, &
        p_chM5s,p_cvM5s,p_qK5_10=>p_q10
  use mem_BenPhosphate, ONLY: p_pK1_an=>p_p_an
  use mem_BenAlkalinity,ONLY: flag_ae_reset,flag_an_reset,p_pG3h=>p_p
  use mem_BenCO2Transport,ONLY:p_pG3c=>p_p
  use constants,only:MW_C

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
!   the Free Software Foundation
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
  real(RLEN)  :: r,s
  real(RLEN)  :: cD1m,cD2m,dn
  real(RLEN)  :: cM5s,cM15s
  real(RLEN)  :: cM1p,cM11p
  real(RLEN)  :: DICdn,AcDn,cMd1p
  integer     :: error

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)

    ! In order to calculate the PH in the bottom before the nutrients gradients in the sediment
    ! are determined we need to !estimate porewater concentration from the state variables.

    !Silica
    cD2m  =   min(  max(  D2m(BoxNumberXY),   p_clK5D2m),  p_chK5D2m)
    cD1m  =   min(  max(  D1m(BoxNumberXY),   p_clK5D1m),  cD2m- p_clD1D2m)
    dn=    (DONE-p_pAn2Ni) * cD2m +p_pAn2Ni*p_d_tot
    cM5s=max(ZERO,K5s(BoxNumberXY)/dn/( DONE+p_pK5_ae) / p_poro(BoxNumberXY))
    cM15s=max(ZERO,K15s(BoxNumberXY)/(p_d_tot_2-cD2m)/  &
                                  (DONE +p_pK5_an)/p_poro(BoxNumberXY))
    cM15s=(cM5s*(dn-cD1m)+cM15s*(p_d_tot_2-dn))/ (p_d_tot_2-cD1m)
    !Phosphate
    cM1p = K1p(BoxNumberXY)/ p_poro(BoxNumberXY)/( &
       p_pK1_ae(BoxNumberXY)+ DONE)/( D1m(BoxNumberXY))
    cMd1p= K11p(BoxNumberXY)/(DONE+p_pK1_ae(BoxNumberXY)) &
          /(D2m(BoxNumberXY)-D1m(BoxNumberXY))/ p_poro(BoxNumberXY)
    cM11p=( K11p(BoxNumberXY)/(DONE+p_pK1_ae(BoxNumberXY))+ &
      K21p(BoxNumberXY)/(DONE+p_pK1_an)) /(p_d_tot_2-D1m(BoxNumberXY))/ &
                                                   p_poro(BoxNumberXY)
    DICae(BoxNumberXY) = G3c(BoxNumberXY)/ MW_C/ p_poro(BoxNumberXY)/ &
        ( p_pG3c+ DONE)/( D1m(BoxNumberXY))
    DICdn = G13c(BoxNumberXY) / MW_C/ p_poro(BoxNumberXY)/ &
       ( p_pG3c+ DONE)/ ( D2m(BoxNumberXY)- D1m(BoxNumberXY))
    DICan(BoxNumberXY) = (G13c(BoxNumberXY)+G23c(BoxNumberXY)) &
        / MW_C/ p_poro(BoxNumberXY)/( p_pG3c+ DONE)/ &
                                       ( p_d_tot_2- D1m(BoxNumberXY))
    Acae(BoxNumberXY) = G3h(BoxNumberXY)/ p_poro(BoxNumberXY)/( &
      p_pG3h+ DONE)/D1m(BoxNumberXY)
    Acdn = G13h(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_pG3h+ DONE) &
                     /( D2m(BoxNumberXY)- D1m(BoxNumberXY))
    Acan(BoxNumberXY) = (G23h(BoxNumberXY)+G13h(BoxNumberXY)) / &
      p_poro(BoxNumberXY)/( p_pG3h+ DONE)/( p_d_tot_2- D1m(BoxNumberXY))

    if (G3h(BoxNumberXY).lt.ZERO ) then
        G3h(BoxNumberXY)=(O3h_Ben(BoxNumberXY)+max(ZERO,Acan(BoxNumberXY)) ) &
              *0.5* p_poro(BoxNumberXY)*( p_pG3h+ DONE)*D1m(BoxNumberXY)
        write(LOGUNIT,*) "G3h<0", 'Reset to' ,G3h(BoxNumberXY)
        call set_warning_for_getm()
    endif


    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the belonging concentrations for the state variables
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     error= CalcCO2System(2,ESW(BoxNumber),&
       ETW(BoxNumber),ERHO(BoxNumber),cM1p,cM5s,Acae(BoxNumberXY),&
       CO2ae(BoxNumberXY),HCO3ae(BoxNumberXY),CO3ae(BoxNumberXY), &
       CAcae(BoxNumberXY), DIC_in=DICae(BoxNumberXY),&
       pCO2_out=pCO2ae(BoxNumberXY), pH_out=pHae(BoxNumberXY))
     if ( error > 0 ) then
          write(LOGUNIT,*)"oBpH: pH outside range"
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:ESW',ESW(BoxNumber)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:ETW',ETW(BoxNumber)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:ERHO',ERHO(BoxNumber)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:DICae',DICae(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:M1p',cM1p
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:M5s',cM5s
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:Acae',Acae(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:D1m',D1m(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',3G13.6)') 'BpH:G3h,G13h,G23h',&
                           G3h(BoxNumberXY),G13h(BoxNumberXY),G23h(BoxNumberXY)
          write(LOGUNIT,'('' pBpH:Hae='',G13.6)') pHae(BoxNumberXY)
          write(LOGUNIT,*) "BenpHDynamics pHae outside range 2-11"
          if ( DICae(BoxNumberXY).lt.Acae(BoxNumberXY).and.flag_ae_reset==1) then
            G3h(BoxNumberXY)=G3c(BoxNumberXY)/MW_C
            write(LOGUNIT,*) 'G3h reset '
          endif
          call set_warning_for_getm()
          pHae(BoxNumberXY)=15.0D+00
      endif

     error= CalcCO2System(2,ESW(BoxNumber),&
                 ETW(BoxNumber),ERHO(BoxNumber),cMd1p,cM15s,Acdn, &
                 CO2an(BoxNumberXY),HCO3an(BoxNumberXY),CO3an(BoxNumberXY), &
                 CAcan(BoxNumberXY),DIC_in=DICdn,&
                 pCO2_out=pCO2an(BoxNumberXY),pH_out=pHdn(BoxNumberXY))
      if ( error > 0 ) then
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:DIC',O3c(BoxNumber)/MW_C
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:DICdn',DICdn
          write(LOGUNIT,'(A,'' ='',G13.6,G13.6)') 'BpH:Acan,G13h', &
                                           Acdn,G13h(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:Md1p',cMd1p
          write(LOGUNIT,'(A,'' ='',2G13.6)') 'BpH:M15s',cM15s,K15s
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:D1m',D1m(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:D2m',D2m(BoxNumberXY)
          write(LOGUNIT,'('' pBpH:pHan='',G13.6)') pHdn(BoxNumberXY)
          write(LOGUNIT,*) "BenpHDynamics:pHdn outside range 2-11"
          pHdn(BoxNumberXY)=8.0
          error= CalcCO2System(2,ESW(BoxNumber),&
            ETW(BoxNumber),ERHO(BoxNumber),cMd1p,cM15s,Acdn, &
            CO2an(BoxNumberXY),HCO3an(BoxNumberXY),CO3an(BoxNumberXY), &
            CAcan(BoxNumberXY),pH_in=pHdn(BoxNumberXY),DIC_in=DICdn,&
            pCO2_out=pCO2an(BoxNumberXY))
            G13h(BoxNumberXY)=Acdn  * p_poro(BoxNumberXY)*( p_pG3h+ DONE) &
                     *( D2m(BoxNumberXY)- D1m(BoxNumberXY))
            write(LOGUNIT,*) 'G13h reset to ',G13h
            call set_warning_for_getm()
      endif
     error= CalcCO2System(2,ESW(BoxNumber),&
                 ETW(BoxNumber),ERHO(BoxNumber),cM11p,cM15s,Acan(BoxNumberXY), &
                 CO2an(BoxNumberXY),HCO3an(BoxNumberXY),CO3an(BoxNumberXY), &
                 CAcan(BoxNumberXY),DIC_in=DICan(BoxNumberXY),&
                 pCO2_out=pCO2an(BoxNumberXY),pH_out=pHan(BoxNumberXY))
      if ( error > 0 ) then
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:DIC',O3c(BoxNumber)/MW_C
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:DICan',DICan(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6,G13.6)') 'BpH:Acan,G3h', &
                                           Acan(BoxNumberXY),G3h(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:M11p',cM11p
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:M15s',cM15s
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:D1m',D1m(BoxNumberXY)
          write(LOGUNIT,'(A,'' ='',G13.6)') 'BpH:D2m',D2m(BoxNumberXY)
          write(LOGUNIT,'('' pBpH:pHan='',G13.6)') pHan(BoxNumberXY)
          write(LOGUNIT,*) "BenpHDynamics:pHan outside range 2-11"
          if ( DICan(BoxNumberXY).lt.Acan(BoxNumberXY).and.flag_an_reset==1)then
            s=max((p_d_tot_2-D1m(BoxNumberXY))*HplusBasis, &
                                 0.5*(G23c(BoxNumberXY)+G13c(BoxNumberXY))/MW_C)
            r=(G13h(BoxNumberXY)+G23h(BoxNumberXY)-s)/MW_C
            G23h(BoxNumberXY)=G23c(BoxNumberXY)/MW_C
            r=min(0.5,r/(NZERO+G23h(BoxNumberXY)))/LocalDelta*G23h(BoxNumberXY)
            call flux(BoxNumberXY,iiBen,ppG23h,ppG23h,-r)
            r=(G13h(BoxNumberXY)+G23h(BoxNumberXY)-s)/MW_C
            G13h(BoxNumberXY)=G13c(BoxNumberXY)/MW_C
            r=min(0.5,r/(NZERO+G13h(BoxNumberXY)))/LocalDelta*G13h(BoxNumberXY)
            call flux(BoxNumberXY,iiBen,ppG13h,ppG13h,-r)
            write(LOGUNIT,*) 'G13h,G23h reset'
          endif
          pHan(BoxNumberXY)=-1
          call set_warning_for_getm()
      else
       pHan(BoxNumberXY)=pHdn(BoxNumberXY)*(D2m(BoxNumberXY)-D1m(BoxNumberXY)) &
       +pHan(BoxNumberXY)*(p_d_tot-D2m(BoxNumberXY))/(p_d_tot-D1m(BoxNumberXY))
      endif
  enddo
#endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
