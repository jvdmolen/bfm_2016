#include "DEBUG.h"
#include "INCLUDE.h"

         subroutine TopPredLosses(iisys,n,spec_loss,p_ex, &
                   state_c,state_n,state_p, &
                   flO3c,flUrean,flUreac,flN1p,flR6c,flR6n,flR6p)
         USE global_mem, ONLY:RLEN, ZERO, DONE,NZERO
         use mem_Param, only:p_peZ_R1c, p_peZ_R1n, p_peZ_R1p
         use constants, only:p_qnUc
         implicit none
         integer, intent(IN)           ::iisys
         integer, intent(IN)           ::n
         real(RLEN),intent(IN)         ::spec_loss(n)
         real(RLEN),intent(IN)         ::p_ex
         real(RLEN), intent(IN)        ::state_c(n)
         real(RLEN), intent(IN)        ::state_n(n)
         real(RLEN), intent(IN)        ::state_p(n)
         real(RLEN), intent(INOUT)     ::flO3c(n)
         real(RLEN), intent(INOUT)     ::flUrean(n)
         real(RLEN), intent(INOUT)     ::flUreac(n)
         real(RLEN), intent(INOUT)     ::flN1p(n)
         real(RLEN), intent(INOUT)     ::flR6c(n)
         real(RLEN), intent(INOUT)     ::flR6n(n)
         real(RLEN), intent(INOUT)     ::flR6p(n)

         real(RLEN)    ::slow_loss_c(n),fast_loss_c(n)
         real(RLEN)    ::slow_loss_n(n),fast_loss_n(n)
         real(RLEN)    ::slow_loss_p(n),fast_loss_p(n)

         slow_loss_c=p_ex*(DONE-p_peZ_R1c)*spec_loss*state_c
         slow_loss_n=p_ex*(DONE-p_peZ_R1n)*spec_loss*state_n
         slow_loss_p=p_ex*(DONE-p_peZ_R1p)*spec_loss*state_p

         fast_loss_n=spec_loss*state_n-slow_loss_n
         fast_loss_p=spec_loss*state_p-slow_loss_p
         fast_loss_c=spec_loss*state_c-slow_loss_c-fast_loss_n/p_qnUc
         flO3c=flO3c+fast_loss_c
         flUrean=flUrean+fast_loss_n
         flUreac=flUreac+fast_loss_n/p_qnUc
         flN1p=flN1p+fast_loss_p

         flR6c=flR6c+ slow_loss_c
         flR6n=flR6n+ slow_loss_n
         flR6p=flR6p+ slow_loss_p
         return
         end
