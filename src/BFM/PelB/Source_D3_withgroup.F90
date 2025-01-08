          function Source_D3_withgroup(iistate,ppgroup,lgroup,type,mode)
          use global_mem, only: RLEN, ZERO,DONE
          use mem,only:D3SINK,D3SOURCE,iiProduction,iiConsumption,iiTotal
          use constants, only:SEC_PER_DAY

          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::lgroup
          integer, intent(IN)             ::type
          integer, intent(IN)             ::mode
          interface                         ! Specification
            integer function ppgroup(n,iiC) ! Specification
            integer,intent(IN)   ::n,iiC    ! Specification
            end function                    ! Specification
          end  interface   
          real(RLEN) :: Source_D3_withgroup(size(D3SOURCE,DIM=1))
          real(RLEN) :: fill(size(D3SOURCE,DIM=1))
          integer    :: i,j,l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          do i=1,lgroup
            j=ppgroup(i,type)
            if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l* D3SINK(:,iistate,j)
            if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D3SOURCE(:,iistate,j)
          enddo
          Source_D3_withgroup=fill*SEC_PER_DAY
        end function Source_D3_withgroup

