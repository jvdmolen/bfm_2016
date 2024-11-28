#include "DEBUG.h"
#include "INCLUDE.h"
         subroutine CalcSiltResuspension()
         use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
         use mem,ONLY: R9x,Q9x,Qp9x,NO_BOXES,NO_BOXES_Z,NO_BOXES_XY, &
                 PelBoxABove,Depth

         implicit none
         integer    :: i,j,k
         real(RLEN) ::new_Q9x
         real(RLEN),external   :: GetDelta

         do i=1,NO_BOXES_XY
           j=PelBoxAbove(i)
           k=j+NO_BOXES_Z-1
           Qp9x(i)=Q9x(i)
           Q9x(i)=sum(R9x(j:k)*Depth(j:k))
         enddo
         end
