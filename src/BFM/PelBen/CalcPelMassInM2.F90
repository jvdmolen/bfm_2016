     function CalcPelMassInM2(ppX)
     USE global_mem, ONLY:RLEN,LOGUNIT
     use mem,only:NO_BOXES,NO_BOXES_XY,D3STATE,NO_BOXES_Z,DEPTH,PelBoxAbove
     integer,intent(IN)      ::ppX
     real(RLEN),dimension(NO_BOXES_XY)  ::CalcPelMassINM2

     integer                ::j,k,n
     real(RLEN),dimension(NO_BOXES_XY)  ::r
     real(RLEN),allocatable  ::sum_arg(:)     !JM added
     integer   :: ppX_local

     ppX_local=ppX
!write(LOGUNIT,*) 'calcpelmassinm2',ppX_local
!stop

!JM bug, j not defined     do i=j,NO_BOXES_XY
!JM bug       k=PelBoxAbove(j)
     do i=1,NO_BOXES_XY
       k=PelBoxAbove(i)
      n=k+NO_BOXES_Z-1
!write(LOGUNIT,*) 'i,k,n',i,k,n
!stop
!JM       r(i)=sum(D3STATE(k:n,ppX)*Depth(k:n))
       allocate(sum_arg(k:n))
!JM       sum_arg=D3STATE(k:n,ppX)
       sum_arg=D3STATE(k:n,ppX_local)
       sum_arg=sum_arg*Depth(k:n)
       r(i)=sum(sum_arg)
       deallocate(sum_arg)
     enddo
     CalcPelMassInM2=r
     return
     end
