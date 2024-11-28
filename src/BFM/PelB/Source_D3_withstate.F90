          function Source_D3_withstate(iistate,jjstate,mode)
          use global_mem, only: RLEN, ZERO,DONE
          use constants, only: SEC_PER_DAY
          use mem,only:D3SINK,D3SOURCE,iiProduction,iiConsumption,iiTotal

          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::jjstate
          integer, intent(IN)             ::mode
          real(RLEN) :: Source_D3_withstate(size(D3SOURCE,DIM=3))
          real(RLEN) :: fill(size(D3SOURCE,DIM=3))
          real(RLEN) :: l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l*D3SINK(iistate,jjstate,:)
          if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D3SOURCE(iistate,jjstate,:)
          Source_D3_withstate=fill*SEC_PER_DAY
        end function Source_D3_withstate

