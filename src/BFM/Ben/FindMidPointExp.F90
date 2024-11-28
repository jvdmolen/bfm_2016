      function FindMidPointExp(Dm,p_from,p_to)
      use global_mem,ONLY:RLEN,DONE
      use mem_globalfun,   ONLY: SetExpDist
      use mem_Param,only:p_clDxm,p_d_tot 
      implicit none
      real(RLEN),intent(IN)          :: Dm
      real(RLEN),intent(IN)          :: p_from
      real(RLEN),intent(IN)          :: p_to
      real(RLEN)                     :: FindMidPointExp
      !local
      real(RLEN)                     ::alpha

      alpha  =   SetExpDist(Dm,Dlm=p_clDxm,Dmm=p_d_tot)
      FindMidPointExp= &
           log(DONE -0.5*(DONE-exp(-(p_to-p_from)*alpha)))/(-alpha)+p_from
      return
      end

      function FindMidPointExp_vector(Dm,p_from,d1_from,p_to,d1_to)
      use global_mem,ONLY:RLEN,DONE
      use mem_globalfun,   ONLY: SetExpDist_vector
      use mem_Param,only:p_clDxm,p_d_tot 
      use mem,only: NO_BOXES_XY
      implicit none
      real(RLEN),dimension(NO_BOXES_XY),intent(IN)          :: Dm
      real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional :: d1_from
      real(RLEN),intent(IN),optional                        :: p_from
      real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional :: d1_to
      real(RLEN),intent(IN),optional                        :: p_to
      real(RLEN),dimension(NO_BOXES_XY)   :: FindMidPointExp_vector
      !local
      real(RLEN),dimension(NO_BOXES_XY)   ::rfrom,rto,alpha

      if (present(d1_from)) then
        rfrom=d1_from
      elseif (present(p_from)) then
        rfrom=p_from
      else
       stop 'error one of the 2 from-parameters is not present in ncall'
      endif
      if (present(d1_to)) then
        rto=d1_to
      elseif (present(p_to)) then
        rto=p_to
      else
       stop 'error one of the 2 to-parameters is not present in ncall'
      endif

      alpha  =   SetExpDist_vector(Dm,Dlm=p_clDxm,Dmm=p_d_tot)
      FindMidPointExp_vector= &
           log(DONE -0.5*(DONE-exp(-(rto-rfrom)*alpha)))/(-alpha)+rfrom
      return
      end

