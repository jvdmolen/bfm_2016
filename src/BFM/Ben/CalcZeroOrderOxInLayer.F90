#include "DEBUG.h"
#include "INCLUDE.h"

     subroutine CalcZeroOrderOxInLayer(G2o,D1m,from,xto,Oo)
     use global_mem, ONLY:RLEN,DONE,NZERO,ZERO
     use mem,ONLY:NO_BOXES_XY

     implicit none
     real(RLEN),dimension(NO_BOXES_XY),INTENT(IN) ::G2o
     real(RLEN),dimension(NO_BOXES_XY),INTENT(IN) ::D1m
     real(RLEN),dimension(NO_BOXES_XY),INTENT(IN) ::from,xto
     real(RLEN),dimension(NO_BOXES_XY),INTENT(out)::Oo

     real(RLEN),dimension(NO_BOXES_XY) ::a,b,c,z0,h

     z0=min(D1m,xto)
     a= 3.0D+00 *G2o/D1m**3
     b=-6.0D+00 *G2o/D1m**2
     c= 3.0D+00 *G2o/D1m

     h=ZERO
     where (from.gt.ZERO) h=a/3.0D+00*from**3+b*0.5D+00*from**2+c*from
     Oo=max(ZERO,a/3.0D+00*z0**3+b*0.5D+00*z0**2+c*z0 -h)

     return
     end
