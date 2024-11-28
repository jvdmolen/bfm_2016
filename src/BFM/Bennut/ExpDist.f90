
      function ExpDist(Dm,Dlm,Dmm)
      use global_mem,ONLY: RLEN,DONE
      IMPLICIT none
!JM      real(RLEN),intent(IN)          :: Dm(:)
!JM      real(RLEN),intent(IN),optional :: Dlm(:)
!JM      real(RLEN),intent(IN),optional :: Dmm(:)
      real(RLEN),intent(IN)          :: Dm
      real(RLEN),intent(IN),optional :: Dlm
      real(RLEN),intent(IN),optional :: Dmm
      real(RLEN)                     :: ExpDist

      real(RLEN)                     :: r

      r=Dm
!JM      if (present(Dlm)) r=sign(max(abs(r),Dlm,r))
!JM      if (present(Dmm)) r=sign(min(abs(r),Dmm,r))
      if (present(Dlm)) then
         if ( max(abs(r),Dlm,r) .ge. 0 ) then
           r=1
         else
           r=-1
         endif
      endif
      if (present(Dmm)) then
          if ( min(abs(r),Dmm,r) .ge. 0 ) then
           r=1
         else
           r=-1
         endif
      endif
      ExpDist=-DONE/r

      return
      end function ExpDist

