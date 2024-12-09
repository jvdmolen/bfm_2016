
      subroutine checkD2Rates(mode,nr,out)
      use global_mem, only:RLEN, ZERO,LOGUNIT,ALLOC,error_msg_prn
      use mem,only:D2SOURCE,D2SINK,NO_D2_BOX_STATES
      implicit none
      integer,intent(IN)          ::mode
      integer,intent(IN)          ::nr
      integer,intent(OUT)         ::out

      integer,save                ::init=0
      integer,pointer,dimension(:):: llsource
      integer,pointer,dimension(:):: llsink
      integer                     :: status,llout,i,j

      if ( mode.eq.0) then
        if ( init.eq.0) then
          allocate(llsource(1:NO_D2_BOX_STATES),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllCMem","llsource")
          allocate(llsink(1:NO_D2_BOX_STATES),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllCMem","llsource")
          init=1
        endif
        llsource = 0
        llsink = 0
        do i=1,NO_D2_BOX_STATES
          if (D2SINK(1,nr,i) > ZERO.and.i.ne.nr) llsink(i)=1
          if (D2SOURCE(1,nr,i) > ZERO.and.i.ne.nr) llsource(i)=1
        enddo
        out=sum(llsink)+sum(llsource)
      elseif (init.eq.0) then
         return
      elseif ( mode.eq.1) then
        llout=0
        do i=1,NO_D2_BOX_STATES
          if (D2SINK(1,nr,i) > ZERO.and.llsink(i)==0.and.i.ne.nr ) llout=llout+1
          if (D2SINK(1,nr,i) == ZERO.and.llsink(i)>0.and.i.ne.nr ) llout=llout+1
          if (D2SOURCE(1,nr,i) > ZERO.and.llsource(i)==0.and.i.ne.nr) llout=llout+1
          if (D2SOURCE(1,nr,i) == ZERO.and.llsource(i)>0.and.i.ne.nr) llout=llout+1
        enddo
        out=llout
      elseif ( mode.eq.2) then
        j=0
        do i=1,NO_D2_BOX_STATES
          if (D2SINK(1,nr,i) > ZERO.and.llsink(i)==0 ) then
             if (j==0) write(LOGUNIT,&
               '(''CheckD2Rate:flux-> only defined in actual step'',$)')
             write(LOGUNIT,'(I4,$)'),i
             j=1
          endif
        enddo
        if (j==1)write(LOGUNIT,*) ''
       j=0
        do i=1,NO_D2_BOX_STATES
          if (D2SINK(1,nr,i) == ZERO.and.llsink(i)>0.and.i.ne.nr ) then 
             if (j==0) write(LOGUNIT,&
             '(''CheckD2Rate:flux-> only defined in previous step'',$)')
             write(LOGUNIT,'(I4,$)'),i
             j=1
          endif
        enddo
        if (j==1) write(LOGUNIT,*) ''
        j=0
        do i=1,NO_D2_BOX_STATES
          if (D2SOURCE(1,nr,i) > ZERO.and.llsource(i)==0)  then
            if (j==0)write(LOGUNIT,&
            '(''CheckD2Rate:flux<- only defined in actual step'',$)') 
            write(LOGUNIT,'(I4,$)'),i
            j=1
          endif
        enddo
        if (j==1)write(LOGUNIT,*) ''
        j=0
        do i=1,NO_D2_BOX_STATES
          if (D2SOURCE(1,nr,i) == ZERO.and.llsource(i)>0.and.i.ne.nr) then
             if (j==0)write(LOGUNIT,&
             '(''CheckD2Rate:flux<- only defined in previous step'',$)')
             write(LOGUNIT,'(I4,$)'),i
            j=1
          endif
        enddo
        if (j==1)write(LOGUNIT,*) ''

         out=0
      endif
      return
      end
