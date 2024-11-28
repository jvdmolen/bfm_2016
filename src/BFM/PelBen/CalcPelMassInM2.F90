     function CalcPelMassInM2(ppX)
     USE global_mem, ONLY:RLEN,LOGUNIT
     use mem,only:NO_BOXES,NO_BOXES_XY,D3STATE,NO_BOXES_Z,DEPTH,PelBoxAbove
     integer,intent(IN)      ::ppX
     real(RLEN),dimension(NO_BOXES_XY)  ::CalcPelMassINM2

     integer                ::j,k,n
     real(RLEN),dimension(NO_BOXES_XY)  ::r

     do i=j,NO_BOXES_XY
       k=PelBoxAbove(j)
       n=k+NO_BOXES_Z-1
       r(i)=sum(D3STATE(ppX,k:n)*Depth(k:n))
     enddo
     CalcPelMassInM2=r
     return
     end
