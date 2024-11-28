#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   BenProfiles.f90
!
! FILE
!   LimitChange.f90
!
! DESCRIPTION
!   
!	function BenProfiles
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine BenProfiles
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
#ifdef INCLUDE_DIAGNOS_PRF
        use global_mem,      ONLY:RLEN,LOGUNIT,ZERO,DONE,NZERO
        use mem, ONLY:BoxNumberZ, NO_BOXES_PRF, BoxNumberX, NO_BOXES_X, &
          BoxNumberY,NO_BOXES_Y,NO_BOXES_XY,BoxNumber,BoxNumberXY, &
          PrQun,PrM1p,PrM3n,PrM4n,PrM5s,PrM6r,PrsM1p,PreM5s,PrBDc,PrBChlC, &
          KQ1c,KQun,KPO4,KNO3,KNH4,KRED,KPO4sh,KSiO3,KSiO3eq,BP1c,iiBP1, &
          BP1l,seddepth,prmdepth,ETW_Ben,iiC,iiL

        use mem_Param, ONLY: p_poro,CalcBenPhyto
        use constants, ONLY: INTEGRAL,STANDARD
        use bennut_interface, ONLY:CalculateFromSet
        use mem_BenSilica, ONLY: p_chM5s, p_cvM5s,p_q10
        use mem_BenPhyto, ONLY: CalculateBenPhyto
        use bio_var,only:diagb_prf
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!       The following global functions are used:eTq
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        use global_interface,   ONLY: eTq,check_if_in_output

!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
       IMPLICIT  NONE
       REAL(RLEN)    ::h,h2,r,s,d,chM5s,rp,sp
       logical::start,llM1p,llsM1p,llM3n,llM4n,lleM5s,llM5s,llM6r, &
                       llQun,llBDc,llChlC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! user defined external functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       integer, external  :: D3toD1
       integer, external  :: D2toD1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
         call check_if_in_output('PrM1p', llM1p)
         call check_if_in_output('PrsM1p',llsM1p)
         call check_if_in_output('PrM3n', llM3n)
         call check_if_in_output('PrM4n', llM4n)
         call check_if_in_output('PrM5s', llM5s)
         call check_if_in_output('PreM5s',lleM5s)
         call check_if_in_output('PrM6r', llM6r)
         call check_if_in_output('PrQun', llQun)
         call check_if_in_output('PrBDc', llBDc)
         call check_if_in_output('PrBChlC', llChlC)
         ! Check if PrpH, profile of pH in the sediment,
         ! will be outputted in the run.
         ! if yes, calculation of profiles of phosphate and silicate
         ! are needed for calculation of the pH profile.
         ! See further in CO2/BenCO2Profiles.F90
         call check_if_in_output('PrpH',  start)
         lleM5s=lleM5s.or.start;llM1p=llM1p.or.start

      sp=ZERO;
      s=ZERO;
      DO BoxNumberZ= 1,NO_BOXES_PRF
      DO BoxNumberY=1,NO_BOXES_Y
        DO BoxNumberX=1,NO_BOXES_X
          BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
          BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)
 
          r=abs(seddepth(BoxNumber))
          d=r-s
          if (llM1p) PrM1p(BoxNumber)= &
           CalculateFromSet(KPO4(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          if (llsM1p.and.KPO4sh(BoxNumberXY).gt.0) PrsM1p(BoxNumber)= &
            CalculateFromSet(KPO4sh(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          if (llM3n) PrM3n(BoxNumber)= &
            CalculateFromSet(KNO3(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          if (llM4n) PrM4n(BoxNumber)= &
            CalculateFromSet(KNH4(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          if (llM5s.or.lleM5s) then
            chM5s = p_chM5s+ &
             p_cvM5s*( eTq( ETW_Ben(BoxNumberXY), p_q10)- DONE)
            if (lleM5s.and.KSiO3eq(BoxNumberXY).gt.0)  &
                PreM5s(BoxNumber)= chM5s- CalculateFromSet( &
                          KSiO3eq(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
            if (llM5s.and.KSiO3(BoxNumberXY).gt.0)  &
                PrM5s(BoxNumber)= chM5s- CalculateFromSet( &
                          KSiO3(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          endif
          if (llM6r) PrM6r(BoxNumber)= &
            CalculateFromSet(KRED(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d
          if (llQun) PrQun(BoxNumber)= &
            CalculateFromSet(KQun(BoxNumberXY),INTEGRAL,STANDARD,s,r)/d

          s=r;
          if ((llBDc.or.llChlc).and.CalcBenPhyto(1)) then 
             rp=abs(prmdepth(BoxNumber))
             d=rp-sp
             if (llBDc.or.llChlC) then 
               h= CalculateBenPhyto(iiC,INTEGRAL,BoxNumberXY, sp,rp)
               h=h*BP1c(BoxNumberXY)/d;
               if (llBDc) PrBDc(BoxNumber)=h
             endif
             if (llChlC) then 
               h2= CalculateBenPhyto(iiL,INTEGRAL,BoxNumberXY, sp,rp)
               h2=h2*BP1l(BoxNumberXY)/d;
               PrBCHlC(BoxNumber)=h2/(NZERO+h)
             endif
          sp=rp;
          endif
           

         ENDDO
      ENDDO
      ENDDO
      
      return
#endif
      end

