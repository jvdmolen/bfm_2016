#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BeLimitNutrietnUptake
!
! DESCRIPTION
!
!
!
!
! !INTERFACE
     subroutine BenLimitNutrientDynamics
!
! !USES:

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Modules (use of ONLY is strongly encouraged!)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     use global_mem, ONLY:RLEN,ZERO,DONE,NZERO,LOGUNIT
#ifdef NOPOINTERS
     use mem,  ONLY: D2STATE
#else
     use mem, ONLY: H1c,HNc,Hac,K1p,K4n,Qun
#endif
     use mem, ONLY:NO_BOXES_XY, &
          iiPhytoPlankton,iiBenPhyto,BenPhyto,iiH1,iiHN,iiC, &
          CoupledtoBDc,fr_lim_HI_o,fr_lim_HI_n,fr_lim_BPI_n,fr_lim_HI_p, &
                  fr_lim_BPI_p,cNIBTc,fr_lim_Ha_n,fr_lim_Ha_o

     use mem_Param,only:p_qon_nitri
     use mem_BenPhyto,only:p_useparams
     use mem_Phyto,  only: p_quPn=>p_qun,p_quPp=>p_qup
     use mem_BenBac, only: p_quHn=>p_qun,p_quHp=>p_qup,p_purH=>p_pur, &
            p_qnHc=>p_qnc,p_qpHc=>p_qpc
     use constants,only:MW_C

     use mem_BenNBac,only: p_qupA,p_qupB,p_qunA,p_qunB,p_qnNc, &
                           p_qnHNc=>p_qnc,p_qpHNc=>p_qpc


     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Implicit typing is never allowe
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     implicit none

     logical                            :: active_Phyto
     integer                            :: i,j,nrphyto
     integer,dimension(NO_BOXES_XY)     :: iOne,Nt
     real(RLEN),dimension(NO_BOXES_XY)  :: bphytoc,Hbc,prodH1_n,prodHN_n, &
                    prodHa_n, prodH1_p,prodHN_p,prodH1_o,prodHN_o,prodHa_o
     real(RLEN),dimension(NO_BOXES_XY)  :: prod_o,prod_n,prod_p,masstc
     real(RLEN),dimension(iiBenPhyto,NO_BOXES_XY)  :: prodBP_n,prodBP_p

     masstc=ZERO;Nt=0;iOne=1;bphytoc=ZERO
     active_Phyto=.FALSE.
     do i=1,iiBenPhyto
        nrphyto=p_useparams(i)
        j= CoupledtoBDc(nrphyto)
        ! if j> iiPhytoplankton there Benphyto is coupled to a Pelagic phyto
        ! however the mass of the benthic component is very low. 
        if( j>0.and.j<= iiPhytoPlankton) then
             active_Phyto=.TRUE. 
             bphytoc=max(ZERO,BenPhyto(i,iiC))
             masstc=masstc+bphytoc
             Nt= Nt+iOne
        endif
     enddo
     masstc=masstc+H1c
     Nt= Nt+iOne
     masstc=masstc+HNc
     Nt= Nt+2*iOne
     cNIBTc=masstc/Nt*0.1

     prod_n=NZERO
     prod_p=NZERO
     prod_o=NZERO
! Phytoplankton
      if (active_Phyto) then
        do i=1,iiBenPhyto
          nrphyto=p_useparams(i)
          if(CoupledtoBDc(nrphyto)>=0) then
            bphytoc=max(ZERO,BenPhyto(i,iiC))
            prodBP_n(i,:)=p_quPn(nrphyto) *(cNIBTc+bphytoc)
            prodBP_p(i,:)=p_quPp(nrphyto) *(cNIBTc+bphytoc)
            prod_n=prod_n+prodBP_n(i,:)
            prod_p=prod_p+prodBP_p(i,:)
          else
            prodBP_n(i,:)=ZERO;
            prodBP_p(i,:)=ZERO;
          endif
        enddo
      else
        prodBP_n(:,:)=ZERO;
        prodBP_p(:,:)=ZERO;
      endif

      prodH1_n=p_quHn*(cNIBTc+H1c)
      prodH1_p=p_quHp*(cNIBTc+H1c)
      prodH1_o=min(prodH1_n*K4n/p_qnHc,prodH1_p*K1p/p_qpHc) &
                                   /(1-p_purH)*p_purH/MW_C
      prod_n=prod_n+prodH1_n
      prod_p=prod_p+prodH1_p
      prod_o=prod_o+prodH1_o

      Hbc=max(ZERO,HNc-Hac)
      prodHN_n=p_qunB*(cNIBTc+Hbc) + p_qunA*(cNIBTc+Hac) 
      prodHa_n=p_qunA*(cNIBTc+Hac) 
      prodHN_p=p_qupB*(cNIBTc+Hbc) + p_qupA*(cNIBTc+Hac) 
      prodHN_o=min(prodHN_n*(Qun+K4n)/p_qnHNc, &
                          prodHN_p*K1p/p_qpHNc)/p_qnNc*p_qon_nitri/MW_C
      prodHa_o=min(prodHa_n*(Qun+K4n)/p_qnHNc, &
               p_qupA*(cNIBTc+Hac) *K1p/p_qpHNc)/p_qnNc*p_qon_nitri/MW_C
      prod_n=prod_n+prodHN_n
      prod_p=prod_p+prodHN_p
      prod_o=prod_o+ prodHN_o
!     write(LOGUNIT,*) prod_o,prodHN_o,prodH1_o

      fr_lim_Ha_n=DONE
      fr_lim_HI_n=DONE
      fr_lim_HI_p=DONE
      fr_lim_HI_o=DONE
      fr_lim_BPI_n=DONE
      fr_lim_BPI_p=DONE
      fr_lim_HI_n(iiH1,:)=prodH1_n/prod_n
      fr_lim_Ha_n=prodHa_n/prod_n
      if ( fr_lim_Ha_n(1)<0.0) then
         write(LOGUNIT,*) 'fr_lim_Ha_N<0.0',prodHa_n,p_qunA,cNIBTc,Hac
         write(LOGUNIT,*) bphytoc,H1c,HNc
      endif
      fr_lim_HI_p(iiH1,:)=prodH1_p/prod_p
      fr_lim_HI_n(iiHN,:)=prodHN_n/prod_n
      fr_lim_HI_p(iiHN,:)=prodHN_p/prod_p
      fr_lim_HI_o(iiH1,:)=prodH1_o/prod_o
      fr_lim_HI_o(iiHN,:)=prodHN_o/prod_o
      fr_lim_Ha_o=prodHa_o/prod_o
      if (active_Phyto) then
        do i=1,iiBenPhyto
          fr_lim_BPI_n(i,:)=prodBP_n(i,:)/prod_n
          fr_lim_BPI_p(i,:)=prodBP_p(i,:)/prod_p
        enddo
      endif
     end
