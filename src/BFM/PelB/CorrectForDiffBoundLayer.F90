#include "DEBUG.h"
#include "INCLUDE.h"
       function CorrectForDiffBoundLayer(mol_diff,radix_single, &
                                               qu_mg,T,Cinf,radix,cells)
       use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
       use mem,only: NO_BOXES
       use constants,only:SEC_PER_DAY
       use mem_globalfun,ONLY: eTq_vector
       use mem_Param,  ONLY:p_q10diff
       use bio_bfm, ONLY: DeriveFromGotm
       use global_interface,ONLY: CalcSh

       implicit none
       real(RLEN),intent(IN)         ::mol_diff
       real(RLEN),intent(IN)         ::radix_single
       real(RLEN),intent(IN)         ::qu_mg
       real(RLEN),intent(IN)         ::T(NO_BOXES)
       real(RLEN),intent(IN)         ::Cinf(NO_BOXES)
       real(RLEN),intent(IN),optional::radix(NO_BOXES)
       real(RLEN),intent(IN),optional::cells(NO_BOXES)
       real(RLEN)                    ::CorrectForDiffBoundLayer(NO_BOXES)


       real(RLEN),parameter          ::p_w=14.2E-9
       real(RLEN),parameter          ::p_radix=3.5e-6
       real(RLEN),parameter          ::rPI= 3.1415926535897D+00
       real(RLEN)                    ::qu_cell
       real(RLEN)                    ::diff(NO_BOXES)
       real(RLEN)                    ::Q(NO_BOXES)
       real(RLEN)                    ::cell(NO_BOXES)
       real(RLEN)                    ::lradix(NO_BOXES)
       real(RLEN)                    ::Sh(NO_BOXES)
       real(RLEN)                    ::E(NO_BOXES)

       qu_cell=p_w*(radix_single/p_radix)**3 * qu_mg
       diff=mol_diff* eTq_vector( T, p_q10diff)
       cell=DONE; if (present(cells)) cell=cells
       lradix=radix_single; if (present(radix)) lradix=radix
       call DeriveFromGotm(1,NO_BOXES,E)
       ! CalculateSherwoodNumber
       Sh=CalcSh(E,lradix,diff)
       ! the area integrated flux Q
       Q=Sh*rPI*4.0D+00*lradix*diff
       !assume  Q= qu Cinf
       CorrectForDiffBoundLayer=Q*Cinf/(Q+cell*qu_cell)
       return
       end function CorrectForDiffBoundLayer
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! DESCRIPTION
! Calculate Sherwood Number
! !INTERFACE
      function CalcSh(E, radix, diffusion)
!
! !USES:
      use global_mem, ONLY:RLEN,DONE
      use mem,ONLY:NO_BOXES
      use constants,ONLY:SEC_PER_DAY

! !INPUT:
      implicit none
      REAL(RLEN),intent(IN)         ::E(NO_BOXES)
      REAL(RLEN),intent(IN)         ::radix(NO_BOXES)
      REAL(RLEN),intent(IN)         ::diffusion(NO_BOXES)
      real(RLEN)                    ::CalcSh(NO_BOXES)

      real(RLEN)                   :: Pe(NO_BOXES)
      real(RLEN)                   :: sh(NO_BOXES)

      !dimensions: m2    (1/s)  (m2/d)  /(d->s)
      Pe= radix**2.0D+00 * E /(diffusion/SEC_PER_DAY)

      where ( Pe < 0.01D+00 )
        sh= ( DONE + 0.29D+00 *  sqrt(Pe))
      elsewhere ( Pe < 100.0D+00 )
        sh= ( 1.014D+00 + 0.51D+00 *  sqrt(Pe))
      elsewhere
        sh= ( 0.55D+00 * Pe**0.3333D+00)
      endwhere
      CalcSh=sh

     end function CalcSh



