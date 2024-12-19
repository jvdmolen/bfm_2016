#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PhaecoystisCalc
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
  subroutine PhaeocystisCalc(mode,phyto,output,input1,param)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: R3c,N3n,N4n,N1p,Pcc,R1c,R1p
#endif
  use mem, ONLY: ETW, ESW,EIR,Nun, &
    NO_BOXES, iiPel, flux_vector,iiC,iiN,iiP,PhytoPlankton,&
    ppP2c,ppP2n,ppP2p,ppP2l,ppP6c,ppP6n,ppP6p,ppP6l,ppPcc,ppR1c,ppR3c,&
    iiP6,ppR6c,ppR1c,dry_z,iiP6,Depth,PI_dw,&
    flPIR6n,flPIR1n, flPIR6p,flPIR1p,flR3R2c,OCDepth,qnPc,qpPc,qlPc,max_change_per_step
  use constants,  ONLY:POSITIVE,MW_H,MW_O,MW_C
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n,p_pe_R1p,p_qpPhc
  use mem_Diffusion,ONLY:p_diff_N1,p_diff_N3,p_diff_N4,p_diff_urea
  use mem_Phaeo
  use mem_Phyto,  ONLY: p_xqn,p_xqp,p_qnRc,p_qpRc,p_lqnlc,p_qplc,p_Ke, &
    p_qun,p_qup,p_qchlc,p_lN3N4n,p_lureaN4n,p_lN1,p_xsize_m,p_thdo


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use LimitRates, ONLY:LimitChange_vector
  use mem_globalfun,   ONLY: exp_limit,exp_limit_scalar,insw_vector
  use global_interface,only:CorrectForDiffBoundLayer
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  IMPLICIT NONE
  real(RLEN),external  ::GetDelta

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)                         :: mode
  integer,intent(IN)                         :: phyto
  real(RLEN),dimension(NO_BOXES),intent(INOUT) :: output
  real(RLEN),dimension(NO_BOXES),intent(IN)  :: input1
  real(RLEN),intent(IN)  :: param

!
!
! !AUTHORS
!   ERSEM group !     P. Ruardij (NIOZ)
!
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES):: phytoc
  real(RLEN),dimension(NO_BOXES):: phyton
  real(RLEN),dimension(NO_BOXES):: phytop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer    :: iout
  real(RLEN)                      :: mol_diff,l_qu
  real(RLEN),dimension(NO_BOXES)  :: nrCells
  real(RLEN),dimension(NO_BOXES)  :: nrCols
  real(RLEN),dimension(NO_BOXES)  :: ColSize
  real(RLEN),dimension(NO_BOXES)  :: et,tN
  real(RLEN),dimension(NO_BOXES)  :: C0
  real(RLEN),dimension(NO_BOXES)  :: sdo
  real(RLEN),dimension(NO_BOXES)  :: rump,rumn
  real(RLEN),dimension(NO_BOXES)  :: runc,sunc
  real(RLEN),dimension(NO_BOXES)  :: Radius,Volume,Content
  real(RLEN),dimension(NO_BOXES)  :: pxSurvivalStruct,pxSurvivalDepth
  real(RLEN),dimension(NO_BOXES)  :: pxSurvivalSalt
  real(RLEN),dimension(NO_BOXES)  :: optimal_growth
  real(RLEN),dimension(NO_BOXES)  :: eun3,euR1n,eup,c3n,cun,cup
  real(RLEN),dimension(NO_BOXES)  :: pe_R1n,pe_R1p,pe_R1c
  real(RLEN),dimension(NO_BOXES)  :: cx_any !any concentration
  real(RLEN),dimension(NO_BOXES)  :: rx_any
  real(RLEN),dimension(NO_BOXES)  :: sx_any !any relative rate
  real(RLEN),dimension(NO_BOXES)  :: ex_any !any  dimensionless number (0->inf)
  real(RLEN),dimension(NO_BOXES)  :: px_any !any  dimensionless number (0->1)
  real(RLEN),dimension(NO_BOXES)  :: xpure_R3c,x_buoancy
  real(RLEN),dimension(NO_BOXES)  :: exp_arg

  real(RLEN)                      :: rn_1,rn_2,actual_depth
  real(RLEN),parameter            :: low_phaeo=1.0D-10
  real(RLEN),parameter            :: low_colsize=1.99D+00


  phytoc = PhytoPlankton(phyto,iiC)
  actual_depth=OCDepth(1)*dry_z(1)

  select case (mode)
    !CALC_SEDIMENTATION:relative sedimentation rate of colonies
    !          Colonies sinks due to size or include more C than the maximum
    !CALC_MORTALITY: Mortality of degradation of colonies and mortality in cells
    case (CALC_SEDIMENTATION,CALC_MORTALITY)
!write(LOGUNIT,*)'calc_sedimentation',mode,CALC_SEDIMENTATION,CALC_MORTALITY
!if (.false.) then
      output=ZERO
      ColSize=DONE
      if ((phyto.ne.iiP6).or.(param.le.NZERO))return
      where (phytoc>low_phaeo)  &
         ColSize=max(DONE,phytoc/(NZERO+Pcc))
      call CalcPhaeoColonyParameters(NO_BOXES,p_xsize_m(phyto),ColSize, &
                                                         Radius,Volume)
      Content=p_wP6c
      nrCols=Pcc/p_wP6c
!write(LOGUNIT,*)'phaeocalc 1,MW_C,nrCols,NZERO',MW_C,nrCols,NZERO

      where (phytoc>low_phaeo.and.Colsize> low_colsize) Content= &
          ColSize*PI_dw(:,phyto)+(R3c+R3c/MW_C*(2.0*MW_H+MW_O))/(NZERO+nrCols)
!write(LOGUNIT,*)'phaeocalc 2'
      call CalcSinking(NO_BOXES,Radius,Content,ETW,ESW,rx_any)
!write(LOGUNIT,*)'phaeocalc 3'

      ex_any=min((qpPc(:,phyto)- p_qplc(phyto))/(p_qpRc(phyto)-p_qplc(phyto)), &
         (qnPc(:,phyto)- p_lqnlc(phyto))/(p_qnRc(phyto)- p_lqnlc(phyto)))
!write(LOGUNIT,*)'phaeocalc 4'

      ex_any=2.0* p_thdo(phyto)* (DONE/( ex_any+ p_thdo(phyto)))
!JM      px_any= exp_limit(- qlPc(:,phyto)/ p_qchlc(phyto)/ p_Ke(phyto)* EIR)
      exp_arg=- qlPc(:,phyto)/ p_qchlc(phyto)/ p_Ke(phyto)* EIR
      px_any=exp_limit(exp_arg)
      ! More light  more buancy less sedimentation
      ! More nutrient limitation less buoancy more sedimentation
      x_buoancy=max(ZERO,ex_any,px_any)*insw_vector(ColSize-DONE)
!     write(LOGUNIT,*),ex_any(NO_BOXES),qpPc(phyto,NO_BOXES),px_any(NO_BOXES),x_buoancy(NO_BOXES)
      rx_any=rx_any*min(DONE,x_buoancy)
      if (mode==CALC_MORTALITY) then
        !Mortality of degradation of colonies and mortality in cells
        pxSurvivalStruct=DONE -rx_any/p_rsP6m
        pxSurvivalDepth=(actual_depth-p_mxDepthm)/p_mxDepthm
        ! survival impossible at low salinities.
        pxSurvivalSalt= (ESW-p_mSal)/p_mSal
        ex_any=min(pxSurvivalStruct,pxSurvivalDepth,pxSurvivalSalt)
!-------------------------------------------------------------------------
!       call findlarge(rx_any,NO_BOXES,1000.0D+00,iout)
!       if ( iout.gt.0) then
!         px_any=exp_limit(-ex_any*p_steep_mort)
!         write(LOGUNIT,*) 'PHaeo',rx_any(iout),pxSurvivalStruct(iout), &
!                                     ex_any(iout),px_any(iout)
!       endif
!-------------------------------------------------------------------------
        ex_any=max(ZERO,exp_limit(-p_steep_mort*ex_any) &
                                           -exp_limit_scalar(-p_steep_mort))
        !total mortality is limited at low phyto concentration to avoid NaN's
        cx_any=DONE;sx_any=ZERO
        where (phytoc>low_phaeo)
           sx_any=ex_any* phytoc/(phytoc+low_phaeo)*param
        end where
        call LimitChange_vector(POSITIVE,sx_any,cx_any,max_change_per_step)
        output=sx_any
      else
        output=rx_any !m/d (sinking rate
      endif
!endif
!write(LOGUNIT,*)'phaeocalc 5',output
!stop
    case (CALC_MORTALITY_CELLS_IN_COLONY)  !20  mortality of cells in colony
      ! after a direct return the value in output keep its value!
      output=ZERO
      if (param.lt.NZERO) return
      if (phyto.ne.iiP6.or.p_smxInC<NZERO )return
      if (param.lt.1.5) then
        et=input1
        !Calculate number of colonies per m3
        nrCols=(low_phaeo+Pcc)/p_wP6c
        !Calculate number of cells per m3
        nrCells=phytoc/p_wP6c
        ! Calculate of Cel per colony
        ColSize=max(DONE,nrCells/nrCols)
        where (ColSize.gt.low_colsize)
         output=et*p_smxInC*ColSize/(ColSize+p_chnxInC)
        endwhere
      elseif (param.lt.2.5) then
         tN=input1
         output=2.0* p_chnxInc* (DONE/( tN+ p_chnxInc))
      else
        stop 'wrong option'
      end if
    case (COLONY_DEGRADATION) !
      if (phyto.ne.iiP6 )return
      sx_any=input1
      where(phytoc<low_phaeo) &
         sx_any=DONE*max(ZERO,Pcc-phytoc)/(NZERO+Pcc)/GetDelta()
      call flux_vector(iiPel,ppPcc,ppPcc,-sx_any*Pcc)
      output=sx_any
    case (CALC_LOC_DET_FLUX)
      !calculate fluxes from Carbon & nutrients collected in the interstial room between the colonies
      if (phyto.ne.iiP6 )return
      phyton = PhytoPlankton(phyto,iiN)
      phytop = PhytoPlankton(phyto,iiP)
      sdo=input1
      cx_any=ZERO;xpure_R3c=ZERO
      pe_R1c=p_pe_R1c;pe_R1n=p_pe_R1n;pe_R1p=p_pe_R1p

       if (sw_detloc<=1) then
         select case (sw_detloc)
            case (-1) ;pe_R1c=DONE;    pe_R1n=DONE;    pe_R1p=DONE
            case (0)  ;pe_R1c=p_pe_R1c;pe_R1n=p_pe_R1n;pe_R1p=p_pe_R1p
            case (1)  ;pe_R1c=ZERO;    pe_R1n=ZERO;    pe_R1p=ZERO
         end select
       else
!        ! under phosphate limitation phae-detritus stick to gether and
!        ! bacteria cannot penetrate well in the material--.
!        ! Hence at complete lim all detr->R6c
!        ex_any=max(NZERO,min(DONE,(qnPc(:,phyto)-p_qnlc(phyto))/ &
!                                       (p_qnRc(phyto)-p_qnlc(phyto))))
!        pe_R1c=p_pe_R1c*ex_any;pe_R1n=p_pe_R1n*ex_any;pe_R1p=p_pe_R1p*ex_any
!        flR3R2c(:)=flR3R2c(:)+l_flR3R2c*pe_R1c/p_pe_R1c
!        l_flR3R2c=l_flR3R2c*(DONE-pe_R1c)/p_pe_R1c
       endif
       cx_any=max(ZERO,phyton-qnPc(:,phyto)*phytoc)
       flPIR1n(:,phyto)=flPIR1n(:,phyto) + sdo* pe_R1n* cx_any
       flPIR6n(:,phyto)=flPIR6n(:,phyto) + sdo* (DONE-pe_R1n)* cx_any
       cx_any=max(ZERO,phytop-qpPc(:,phyto)*phytoc)
       flPIR1p(:,phyto)=flPIR1p(:,phyto) + sdo* pe_R1p* cx_any
       flPIR6p(:,phyto)=flPIR6p(:,phyto) + sdo* (DONE-pe_R1p)* cx_any
       call flux_vector(iiPel,ppR3c,ppR1c,sdo* pe_R1c* R3c)
       call flux_vector(iiPel,ppR3c,ppR6c,sdo*(DONE-pe_R1c)*R3c)
       rx_any=ZERO;ColSize=DONE
      where (phytoc>low_phaeo)  &
         ColSize=max(DONE,phytoc/(NZERO+Pcc))
       where (R3c>low_phaeo.and.ColSize<low_colsize)  &
           rx_any=1000.00* R3c/(NZERO+phytoc)*R3c
       call LimitChange_vector(POSITIVE,rx_any,R3c,max_change_per_step)
       flR3R2c(:)=flR3R2c(:)+rx_any
!      call findnan(rx_any,NO_BOXES,iout)
!      if ( iout>0) write(LOGUNIT,*) 'Phaeo flR3R3=NaN',iout,rx_any(iout)
       output=(sdo*(DONE-pe_R1c)*xpure_R3c)*(p_pe_R1c-pe_R1c)/p_pe_R1c
    case (LIMIT_BOUND)
        output=ZERO
        where (input1.gt.NZERO.and.output.gt.low_phaeo) &
           output= output -min(output/(NZERO+input1),max_change_per_step) &
                                               *output*GetDelta()
    case (NEW_COLONIES)  !2 transfer of new cells which form colonies
      !in very shallow areas (<5meter) no colonies can be formed:
      !there is too much physical stress...
      sunc=input1;runc=sunc*phytoc ;ex_any=ZERO
      if (actual_depth>p_mxDepthm )then

        eun3 =p_lN3N4n(phyto) /(NZERO+p_lN3N4n(phyto)+N4n)
        c3n=eun3*N3n
        euR1n =p_lureaN4n(phyto) /(NZERO+p_lureaN4n(phyto) +N4n)
        cun=euR1n * Nun
        eup  =p_lN1(phyto) /(NZERO+p_lN1(phyto) + N1p)
        cup=eup*max(ZERO,R1p(:)- R1c(:) * p_qpPhc)  ! rich R1 :

!     !only activity when input( =sum) is larger zero

        rumn(:)  = p_qun(iiP6)* (N4n+c3n(:)+cun)  !* phytoc*fr_lim_PI_n(:,phyto)
        rump(:)  = p_qup(iiP6)* (N1p(:)+cup)  !*fr_lim_PI_p(:,phyto)

        where (runc >NZERO)
          optimal_growth=min(DONE,min(rumn/(p_xqn(phyto)*p_qnRc(phyto)), &
                         rump/(p_xqp(phyto)*p_qpRc(phyto)))/max(NZERO,sunc))
        end where
!       Output2d_1(1)=runc(NO_BOXES)
!       Output2d_2(1)=optimal_growth(NO_BOXES)

        where (runc >NZERO)
          optimal_growth=optimal_growth*exp(-p_steep_init*abs(N1p-p_cmP2P6p) &
                                                        /min(p_cmP2P6p,N1p))
          cx_any=(N4n+c3n(:)+cun)
          ex_any=exp(-p_steep_init*abs(cx_any-p_cmP2P6n )/min(p_cmP2P6n,cx_any))
          ex_any= max(ZERO,optimal_growth,optimal_growth*ex_any)* p_sP2P6  &
            *max(ZERO,DONE-exp(ETW-16.5))
        end where
      endif
!     Output2d_2(1)=optimal_growth(NO_BOXES)
!     Output2d_4(1)=ex_any(NO_BOXES)

      call flux_vector(iiPel,ppP2c,ppP6c,ex_any*runc)
      call flux_vector(iiPel,ppP2n,ppP6n,ex_any*runc*p_qnRc(phyto))
      call flux_vector(iiPel,ppP2p,ppP6p,ex_any*runc*p_qpRc(phyto))
      call flux_vector(iiPel,ppP2l,ppP6l,ex_any*runc*p_qchlc(phyto))
      call flux_vector(iiPel,ppPcc,ppPcc,ex_any*runc)
      output=ex_any
    !11,13,14,15,nutrient  limitation  of cells in colonies
    case (CALC_REL_PHOSPHATE_UPTAKE,CALC_REL_NITRATE_UPTAKE, &
                         CALC_REL_AMMONIUM_UPTAKE,CALC_REL_UREA_UPTAKE)
      output=DONE;mol_diff=ZERO
      if (phyto.ne.iiP6 )return
      select case (mode)
        case(CALC_REL_PHOSPHATE_UPTAKE); mol_diff=p_diff_N1;l_qu=p_qup(phyto)
        case(CALC_REL_NITRATE_UPTAKE)  ; mol_diff=p_diff_N3;l_qu=p_qun(phyto)
        case(CALC_REL_AMMONIUM_UPTAKE) ; mol_diff=p_diff_N4;l_qu=p_qun(phyto)
        case(CALC_REL_UREA_UPTAKE)    ;mol_diff=p_diff_urea;l_qu=p_qun(phyto)
      end select
      !Calculate number of colonies per m3
      nrCols=(low_phaeo+Pcc)/p_wP6c
      !Calculate number of cells per m3
      nrCells=phytoc/p_wP6c
      ! Calculate of Cel per colony
      ColSize=max(DONE,nrCells/nrCols)
      ! output: Radius and Volume of colony
      call CalcPhaeoColonyParameters(NO_BOXES,  &
                                       p_xsize_m(phyto),ColSize,Radius,Volume)
      !-------------------------------------------------------------------------
      !	Start of Calculation of limitation of nutrient uptake
      ! in "input1" externe nutrient concentration
      C0=input1
      cx_any=C0
      where (ColSize.gt.low_colsize)
        !Calc concentration at outside of  a single cell.
        cx_any=  CorrectForDiffBoundLayer(mol_diff, &
                             p_xsize_m(phyto),l_qu,ETW,input1)
        !Calc concentration at outside of colony
        C0= CorrectForDiffBoundLayer(mol_diff, &
                             p_xsize_m(phyto),l_qu,ETW, input1,Radius,ColSize)
        ! use of quotient of concentration at outside of single and colony
        !as a measure for the limitations of the nutrient uptake
        output=(NZERO+C0)/(NZERO+cx_any)
      end where
      output=min(DONE,output)
    case (CALC_NET_NITROGEN_UPTAKE)
      output=ZERO
      if (phyto.eq.iiP6) output=input1*p_qnR3c
    case (CALC_NET_PHOSPHATE_UPTAKE)
      output=ZERO
      if (phyto.eq.iiP6) output=input1*p_qpR3c
    case (CALC_REL_RADIUS) !relative sedimentation rate of colonies
      output=DONE
      if (phyto.ne.iiP6.or.param.ne.NZERO)return
      ColSize=DONE
      where (phytoc>low_phaeo)  &
         ColSize=max(DONE,phytoc/(NZERO+Pcc))
      call CalcPhaeoColonyParameters(NO_BOXES,p_xsize_m(phyto),ColSize, &
                                                           Radius,Volume)
      output=max(DONE,Radius/p_xsize_m(phyto))
    case (CALC_FOOD_MESOZOO,CALC_FOOD_MICROZOO, &
                             CALC_GRAZING_MESOZOO,CALC_GRAZING_MICROZOO)
      output=input1;ex_any=ZERO
      if (phyto.ne.iiP6.or.param.eq.ZERO )return
      if (sum(phytoc*Depth)/OCDepth(1)<low_phaeo ) then
        output=ZERO;return
      endif
      rn_1=ZERO
      select case (mode)
        case (CALC_FOOD_MESOZOO,CALC_GRAZING_MESOZOO)
           rn_1=p_cuZ4a(1); rn_2=p_cuZ4a(2)
        case (CALC_FOOD_MICROZOO,CALC_GRAZING_MICROZOO)
           rn_1=p_cuZ5a(1); rn_2=p_cuZ5a(2)
      end select
      ! Calculate number of Cells per colony
      ColSize=max(DONE,phytoc/(NZERO+Pcc))
      ! Only for the given size range the results of this calculation ==0
      where (phytoc > low_phaeo.and.Pcc.gt.low_phaeo)
        ex_any=min(DONE,(abs(ColSize-rn_1) +abs(rn_2-ColSize)-(rn_2-rn_1))/ColSize)
        ex_any=max(ZERO,DONE-ex_any)
      endwhere
      select case (mode)
         case (CALC_GRAZING_MESOZOO,CALC_GRAZING_MICROZOO)
           rx_any=ZERO
           where(ex_any>NZERO)rx_any=input1/(NZERO+phytoc)*Pcc/ColSize
           call flux_vector(iiPel,ppPcc,ppPcc,-rx_any)
           output=rx_any
         case DEFAULT
            output=input1*ex_any
      end select
  end select

!endif !skip option
!write(LOGUNIT,*)'phaeocalc end'
!stop

  end subroutine PhaeocystisCalc
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! !INTERFACE
  subroutine PhaeocystisCalc_1l(mode,phyto,output,input1,param)
!
! !USES:
  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#endif
  use mem, ONLY: NO_BOXES_XY, iiPel, &
    ppPcc,Pcc,iiP6,P6c,PelBoxAbove,BoxNumberXY,BoxNumber
  use mem_Phaeo
  use mem_phyto,only:p_xsize_m
  use botflux,only:openbotflux_vector
  use constants,only:NEGATIVE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  implicit none
  integer,intent(IN)                            :: mode
  integer,intent(IN)                            :: phyto
  real(RLEN),intent(OUT),dimension(NO_BOXES_XY) :: output
  real(RLEN),intent(IN),dimension(NO_BOXES_XY)  :: input1
  real(RLEN),intent(IN)                         :: param

  real(RLEN),external  ::GetDelta
  real(RLEN),dimension(NO_BOXES_XY)          ::lphytoc
  real(RLEN),dimension(NO_BOXES_XY)          ::lPcc
  real(RLEN),dimension(NO_BOXES_XY)          ::s
  real(RLEN),dimension(NO_BOXES_XY)          ::q
  real(RLEN),dimension(NO_BOXES_XY)          ::ex_any
  real(RLEN),dimension(NO_BOXES_XY)          ::rx_any
  real(RLEN),dimension(NO_BOXES_XY)          ::Radius,Volume
  real(RLEN),dimension(NO_BOXES_XY)          ::ColSize
  real(RLEN),dimension(NO_BOXES_XY)          ::nrCols
  real(RLEN),dimension(NO_BOXES_XY)          ::nrCells
  real(RLEN)                                 ::r2(2)
! real(RLEN)                      :: delta
  real(RLEN),parameter            :: rPI= 3.1415926535897D+00
  real(RLEN),parameter            :: low_phaeo=1.0D-10

  DO BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)
    lphytoc(BoxNumberXY)=P6c(BoxNumber)
    lPcc(BoxNumberXY)=Pcc(BoxNumber)
  enddo

  select case (mode)
    case (CALC_FOOD_YOUNG_FILTERFEEDER,CALC_FOOD_FILTERFEEDER, &
       CALC_GRAZING_FILTERFEEDER,CALC_GRAZING_YOUNG_FILTERFEEDER)
      output=input1 ;ex_any=ZERO
      if (phyto.ne.iiP6.or.param.eq.ZERO )return
      select case (mode)
        case (CALC_FOOD_FILTERFEEDER,CALC_GRAZING_FILTERFEEDER)
           r2=p_cuY3a  !read range min->max
        case (CALC_FOOD_YOUNG_FILTERFEEDER,CALC_GRAZING_YOUNG_FILTERFEEDER)
           r2=p_cuYy3a  !read range min->max
      end select
      ! Calculate number of Cells per colony
      ColSize=DONE
      where (lphytoc>low_phaeo.and.lPcc>low_phaeo) &
        ColSize=max(DONE,lphytoc/(NZERO+lPcc))
      ! Only for the given size range the results of this calculation ==0
      ex_any=min(DONE,(abs(ColSize-r2(1)) +abs(r2(2)-ColSize)-(r2(2)-r2(1)))/ColSize)
      ex_any=max(ZERO,DONE-ex_any)
      select case (mode)
        case (CALC_GRAZING_FILTERFEEDER,CALC_GRAZING_YOUNG_FILTERFEEDER)
          rx_any=ZERO
          where (ex_any>NZERO) rx_any=input1*lPcc/(NZERO+lphytoc)/ColSize
          call openbotflux_vector(NEGATIVE,iiPel,ppPcc,-rx_any)
          output=rx_any
      case DEFAULT
         output=input1*ex_any
      end select
    case (CALC_LIMIT_FILTERCAP)
      ! Hampering of filtering of filterfeeders by Phaeocystis
      output=DONE
      if (phyto.ne.iiP6.or.param.eq.ZERO)return
     !Calculate number of colonies per m3
      nrCols=1
      where (lphytoc>low_phaeo.and.lPcc>low_phaeo) &
                           nrCols=(low_phaeo+lPcc)/p_wP6c
      !Calculate number of cells per m3
      nrCells=lphytoc/p_wP6c
      ! Calculate of Cel per colony
      ColSize=max(DONE,nrCells/nrCols)
      call CalcPhaeoColonyParameters(NO_BOXES_XY, &
                            p_xsize_m(phyto),ColSize,Radius,Volume)
      s=4.0*rPI*Radius*Radius
      !Calculate total surface
      q=s*nrCols
      !Calculate filterlimitations according
      !Smaal,AC && Twisk,F(1997), JournalOF Exp. Mar. Biologyand Ecology,209:33-46
      !(unit of Smaal) mm2/ml ====  m2/m3/ (unit of BFM)
      ! use only to calculate the fraction of the maximum filtration
      output=max(ZERO,0.774D+00-0.12D+00*q)/0.774D+00
  end select
  end subroutine PhaeocystisCalc_1l
  subroutine CalcPhaeoColonyParameters(n,RadiusCell,ColSize,Radius,Volume)
  use global_mem, ONLY:RLEN,ZERO,DONE
  implicit none
  integer,intent(IN)           :: n
  real(RLEN),intent(IN)        ::RadiusCell      ! raidus single   cell
  real(RLEN),intent(IN)        ::ColSize(1:n)
  real(RLEN),intent(OUT)       ::Radius(1:n)
  real(RLEN),intent(OUT)       ::Volume(1:n)

  real(RLEN),parameter            :: rPI= 3.1415926535897D+00
  real(RLEN),parameter            :: mult=4.0D+00/3.00D+00
  real(RLEN),parameter            :: ppower=3.00D+00
  real(RLEN),parameter            :: low_colsize=1.99D+00

  Radius=ZERO
  where (ColSize.gt.low_colsize)
    !Calculate the volume of one colony according Rousseau(1990) and
    !transfer from mm3-->m3
    Volume=1.0D-9*10.0D+00**((log10(ColSize)-3.67D+00)/0.51D+00)
    !Calculate from volume -->diameter
    Radius=(Volume/(mult*rPI))**(DONE/ppower)
  endwhere
  where(Radius.lt.RadiusCell)
     Radius=RadiusCell;Volume=mult*rPI*Radius**ppower
  endwhere
  end
  subroutine CalcSinking(n,Radius,Content,ETW,ESW,SinkingRate)
  use global_mem, ONLY:RLEN,LOGUNIT
  use constants, ONLY:SEC_PER_DAY
  use global_interface,ONLY:CalcWaterProp
  implicit none
  integer,intent(IN)           ::n
  real(RLEN),intent(IN)        ::Radius(1:n)  !m
  real(RLEN),intent(IN)        ::Content(1:n) !mg dryweight /m3
  real(RLEN),intent(IN)        ::ETW(1:n)
  real(RLEN),intent(IN)        ::ESW(1:n)
  real(RLEN),intent(OUT)       ::SinkingRate(1:n) !m/day

  real(RLEN)                   :: viscosity(1:n)
  real(RLEN),parameter         :: rPI= 3.1415926535897D+00
! real(RLEN),parameter         :: q_xdry_c=2.5
  real(RLEN),parameter         :: g=9.81
  real(RLEN),parameter         :: multi1=2.0D+00/9.00
  real(RLEN),parameter         :: multi2=4.0D+00/3.00

     !viscosity in kg /M3/s
     call CalcWaterProp(n,ESW,ETW,Mu=viscosity)
     viscosity=viscosity*1.0D+6
     !viscosity in mg /M3/s
     SinkingRate=multi1*Content*g &
                            /(multi2*rPI*Radius)/viscosity*SEC_PER_DAY

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

