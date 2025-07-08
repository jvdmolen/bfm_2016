!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: coupling_getm_bfm
!
! !INTERFACE:
   module coupling_getm_bfm
!
! !BFM:
!    This is special routine for BFM  
!    Ths routine include code for:
!       1. to fill values into the diagnostic values.
!       2. to caluclate averages per time step for the diagnostic variables.
!       3. calculate the sums which are flagged for output
!                   (see more in Genreal/GlobalDefsBFM.model)
!       4. river input calculations:
!            -set limits on input when no data are available
!            -set boundary condition for tracking-state variables.
!
! !USES:
    use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
    use domain, only: az,au,av
    use domain, ONLY:H,imin,imax,jmin,jmax
    use exceptions, only:getm_error

    character(len=80),private                 :: msg="",sub
! !PUBLIC DATA MEMBERS:
   character(len=64),allocatable,dimension(:) :: tracked_river_name
   character(len=64)                          :: tracer_type
   REALTYPE                                   :: lon_min 
   REALTYPE                                   :: lat_min=-9999.0 
   REALTYPE                                   :: lon_max 
   REALTYPE                                   :: lat_max 
#ifdef INCLUDE_DIAGNOS_PRF
   integer                                    :: diag_end_sections(3,7)
#else
   integer                                    :: diag_end_sections(3,6)
#endif
   integer                                    :: start_tracking_in_jul
   logical,public                             :: read_poro=.false.
   integer,parameter                          :: DIAG_RESET=0,DIAG_ADD=1,&
                                                   DIAG_AVERAGE=2,DIAG_INFO=3

!JM   public fill_diagn_bfm_vars, init_pel_co2, &
!          getm_bfm_bennut_calc_initial, set_2d_grid_parameters, &
!          init_2d_grid,unlabeled_var_index, &
!          assign_river_inputs_to_bio_states,check_reset_tracking, &
!          check_3d_track_rates,make_uv_flux_output,make_river_flux_output, &
!          output_2d_grid_parameters
   public fill_diagn_bfm_vars, init_pel_co2, &
          getm_bfm_bennut_calc_initial, set_2d_grid_parameters, &
          init_2d_grid,unlabeled_var_index, &
          check_reset_tracking, &
          check_3d_track_rates,make_uv_flux_output, &
          output_2d_grid_parameters
!EOP
!-----------------------------------------------------------------------

   contains
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  fill_diagn_bfm_vars
! 
! !INTERFACE:
   subroutine fill_diagn_bfm_vars(mode,llwrite, ig,jg, h,nlev,dt )
!
! !USES:
#ifdef BFM_GOTM
!JM  use bio_var, only: stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS
!JM  use bio_var, only: stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE
  use bfm_output, only: stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS
  use bfm_output, only: stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE
  use variables_bio_3d, only:counter_reset,counter_ave,n_cc3d_out,n_ccb3d_out,flag_out,bio_missing
#ifdef INCLUDE_DIAGNOS_PRF
   use variables_bio_3d, only: ccb3d_prf,n_ccb3d_prf
!JM   use bio_var, only: numbc_prf,nprf,stPRFDiagS,stPRFDiagE
   use bio_var, only: numbc_prf,nprf
   use bfm_output, only: stPRFDiagS,stPRFDiagE
#endif
#endif
!
! !INPUT PARAMETERS: 
   implicit none 
   logical,intent(IN)                        :: llwrite
   integer,intent(IN)                        :: mode
   integer,intent(IN)                        :: ig
   integer,intent(IN)                        :: jg
   integer,intent(IN)                        :: nlev
   REALTYPE,intent(IN),dimension(0:nlev)     :: h
   REALTYPE,intent(IN)                       :: dt
!
!
! !DESCRIPTION: 
!      With this routine all the output prepared for all non-state variables 
!      and for state variables of which the average value hato be collected
!      between 2 subsequent output steps.
!
!      All values are stored in the array cc3d_out :
!      The squenbces in the array is:
!       1, average values:
!                      a. state variables
!                      b. diagnostic variables which are calculated in the model.
!                      c. flux variables.
!      2. normal values:
!                      b. diagnostic variables
!                      c. normal variables. 
! !BUGS:  
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES: 
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef BFM_GOTM
  integer                         ::k
  integer                         ::l

   select case (mode)
    case(DIAG_RESET)
        counter_ave(imin:imax,jmin:jmax)=0; 
    case(DIAG_ADD)
        k=0
        ! add value for average output:
        call fill_diag_array(1,.true., stPelStateS,stPelStateE,ig,jg, nlev,dt,k )
        call fill_diag_array(2,.true., stPelDiagS,stPelDiagE,ig,jg, nlev,dt, k )
        call fill_diag_array(3,.true., stPelFluxS,stPelFluxE,ig,jg, nlev,dt, k )
        if ( llwrite) then
          ! set value for normal output after the average output:
          call fill_diag_array(2,.false., stPelDiagS,stPelDiagE,ig,jg, nlev,dt,k )
          call fill_diag_array(3,.false., stPelFluxS,stPelFluxE,ig,jg, nlev,dt,k )
        endif
        if ( k > n_cc3d_out) then
           STDERR "Error: k > n_cc3d_out"
           STDERR "k=",k,"n_cc3d_out=",n_cc3d_out
        endif
        k=0
        ! add value for average output:
        call fill_diag_array(4,.true., stBenStateS,stBenStateE,ig,jg, nlev,dt,k )
        call fill_diag_array(5,.true., stBenDiagS,stBenDiagE,ig,jg, nlev,dt,k )
        call fill_diag_array(6,.true., stBenFluxS,stBenFluxE,ig,jg, nlev,dt,k )
        if ( llwrite) then
          ! set value for normal output after the average output:
          call fill_diag_array(5,.false., stBenDiagS,stBenDiagE,ig,jg, nlev,dt,k )
          call fill_diag_array(6,.false., stBenFluxS,stBenFluxE,ig,jg, nlev,dt,k )
        endif
        if ( k > n_ccb3d_out) then
           STDERR "Error: k > n_ccb3d_out"
           STDERR "k=",k,"n_ccb3d_out=",n_ccb3d_out
        endif
#ifdef INCLUDE_DIAGNOS_PRF
        k=0
        call fill_diag_array(7,.true., stPRFDiagS,stPRFDiagE,ig,jg, nprf,dt,k )
        if ( llwrite) & 
           call fill_diag_array(7,.false., stPRFDiagS,stPRFDiagE,ig,jg, nprf,dt,k )
        if ( k > n_ccb3d_prf) then
           STDERR "Error: k > n_ccb3d_prf"
           STDERR "k=",k,"n_ccb3d_prf=",n_ccb3d_prf
        endif
#endif
        counter_ave(ig,jg)=counter_ave(ig,jg)+1
        flag_out=1
    case(DIAG_AVERAGE)
        ! calculate the average value: divide values for average output by countere_ave
        k=0
        call fill_diag_array(11,.true., stPelStateS,stPelFluxE,0,0, nlev,dt, k )
        call fill_diag_array(11,.false.,stPelDiagS, stPelFluxE,0,0, nlev,dt, k )
        k=0
        call fill_diag_array(14,.true., stBenStateS,stBenFluxE,0,0, nlev,dt, k )
        call fill_diag_array(14,.false.,stBenDiagS, stBenFluxE,0,0, nlev,dt, k )
#ifdef INCLUDE_DIAGNOS_PRF
        k=0
        call fill_diag_array(15,.true., stPRFDiagS,stPRFDiagE,0,0, nprf,dt, k )
        call fill_diag_array(15,.false.,stPRFDiagS,stPRFDiagE,0,0, nprf,dt, k )
#endif
    case(DIAG_INFO)
        ! calculate the average value: divide values for average output by countere_ave
        LEVEL2 "Detailed info about diagnostic Variables selected for output"
        k=0
        diag_end_sections(1,1) =k;l=k;
        call fill_diag_array(21,.true., stPelStateS,stPelStateE,0,0, nlev,dt,k )
        if ( k-l> 0 ) LEVEL2 "section 1:Averaged Pelagic States start at",l+1
        diag_end_sections(1,2) =k;l=k;
        call fill_diag_array(21,.true., stPelDiagS,stPelDiagE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 2:Averaged Pelagic Diagnostic Variabeles start at",l+1
        diag_end_sections(1,3) =k;l=k;
        call fill_diag_array(22,.true., stPelFluxS,stPelFluxE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 3:Averaged Pelagic Calculations  start at",l+1
        diag_end_sections(1,4) =k;l=k;
        call fill_diag_array(21,.false.,stPelDiagS,stPelDiagE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 4:Sampled Pelagic Diagnostic Variables start at",l+1
        diag_end_sections(1,5) =k;l=k;
        call fill_diag_array(22,.false.,stPelFluxS,stPelFluxE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 5:Sampled Pelagic Calculations  start at",l+1
        diag_end_sections(1,6) =k;l=k;
        k=0
        diag_end_sections(2,1) =k;l=k;
        call fill_diag_array(24,.true., stBenStateS,stBenStateE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 1:Averaged Benthic States start at",l+1
        diag_end_sections(2,2) =k;l=k;
        call fill_diag_array(24,.true., stBenDiagS,stBenDiagE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 2:Averaged Benthic Diagnostic Variabeles start at",l+1
        diag_end_sections(2,3) =k;l=k;
        call fill_diag_array(25,.true., stBenFluxS,stBenFluxE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 3:Averaged Benthic Calculations  start at",l+1
        diag_end_sections(2,4) =k;l=k;
        call fill_diag_array(24,.false.,stBenDiagS,stBenDiagE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 4:Sampled Benthic Diagnostic Variables start at",l+1
        diag_end_sections(2,5) =k;l=k;
        call fill_diag_array(25,.false.,stBenFluxS,stBenFluxE,0,0, nlev,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 5:Sampled Benthic Calculations  start at",l+1
        diag_end_sections(2,6) =k;l=k;
#ifdef INCLUDE_DIAGNOS_PRF
        k=0
        call fill_diag_array(27,.true., stPRFDiagS,stPRFDiagE,0,0, nprf,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 2:Averaged BenProfile Diagnostic Variabeles start at",l+1
        diag_end_sections(3,3) =k;l=k;
        call fill_diag_array(28,.false.,stPRFDiagS,stPRFDiagE,0,0, nprf,dt, k )
        if ( k-l> 0 ) LEVEL2 "section 4:Sampled BenProfile Diagnostic Variables start at",l+1
        diag_end_sections(3,6) =k;l=k;
#endif
    end select

#endif
   end subroutine fill_diagn_bfm_vars

!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: fill_diag_array 
! 
! !INTERFACE:
!
   subroutine fill_diag_array(mode,ave, from,to,ig,jg,nlev,dt,k)
! !USES:
#ifdef BFM_GOTM
   use variables_bio_3d,only:cc3d_out,ccb3d_out,n_cc3d_out,n_ccb3d_out,flag_out
   use variables_bio_3d,only:counter_ave,bio_missing
!JM   use bio_var, only: diag, diagb,var_ids,var_ave,var_names
   use bio_var, only: diag, diagb
   use bfm_output, only: var_ids,var_ave,var_names
   use bio_var,only: c1dimz,cc,ccb,numbc
#ifdef INCLUDE_DIAGNOS_PRF
   use variables_bio_3d, only: ccb3d_prf,n_ccb3d_prf
   use bio_var, only: diagb_prf,numbc_prf
#endif
#endif
!
! !INPUT PARAMETERS: 
!
   implicit none 
   integer,intent(IN)                        :: mode
   integer,intent(IN)                        :: from
   logical,intent(IN)                        :: ave
   integer,intent(IN)                        :: to
   integer,intent(IN)                        :: ig
   integer,intent(IN)                        :: jg
   integer,intent(IN)                        :: nlev
   REALTYPE,intent(IN)                       :: dt
! !INPUT/OUTPUT PARAMETERS:
!
   integer,intent(INOUT)                     :: k

! !DESCRIPTION: 
!   fill 4d array with columns calculation in an efficient way
!
! !LOCAL VARIABLES:
!
#ifdef BFM_GOTM
   integer                        :: i
   integer                        :: j
   integer                        :: l
   integer                        :: m
   integer                        :: n
   integer                        :: k2
   logical                        :: llcalc

! !BUGS:  
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES: 
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
!  2006-05-22   Piet Ruardij  Initial code.
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (mode)
     case(21,22,24,25,27,28)
        j=0
        do i=from,to
          j=j+1
          if (ave) then
            if ( (var_ids(i) /= 0 )  .and. var_ave(i) ) then 
              k=k+1
              select case(mode)
                case(21,22)
                  LEVEL3 '3d averaged: ',trim(var_names(i)),' seq_nrs:',k
                case(24,25)
                  LEVEL3 '2d averaged: ',trim(var_names(i)),' seq_nrs:',k
#ifdef INCLUDE_DIAGNOS_PRF
                case(27,28)
                  LEVEL3 '2d averaged: ',trim(var_names(i)),' seq_nrs:',k
#endif
              end select
            endif
          elseif ( (var_ids(i) /= 0) .and. ( .not.var_ave(i)) ) then 
            k=k+1
            select case(mode)
              case(21,22)
                  LEVEL3 '3d sampled: ',trim(var_names(i)),' seq_nrs:',k
                  if ( mode.eq.22) then
                      call test_flux_output(1,j,var_ave(i),l)
                      if ( l.eq.1) var_ids(i)=0
                  endif
              case(24,25)
                  LEVEL3 '2d sampled: ',trim(var_names(i)),' seq_nrs:',k
                  if ( mode.eq.25) then
                      call test_flux_output(2,j,var_ave(i),l)
                      if ( l.eq.1) var_ids(i)=0
                  endif
#ifdef INCLUDE_DIAGNOS_PRF
              case(27,28)
                  LEVEL3 '3d-prf sampled: ',trim(var_names(i)),' seq_nrs:',k
                  if ( mode.eq.25) then
                      call test_flux_output(2,j,var_ave(i),l)
                      if ( l.eq.1) var_ids(i)=0
                  endif
#endif
            end select
          endif
        enddo
      case(11,14,15)
        if ( flag_out ==0) return
        do i=from,to
          if ( ave ) then
            if ( (var_ids(i) /= 0 )  .and. var_ave(i) ) then 
              k=k+1
              select case(mode)
                case(11)
                  do n=jmin,jmax
                    do m=imin,imax
                      if (az(m,n) .ge. 1 .and. counter_ave(m,n).gt.0 ) then
                         cc3d_out(m,n,0:nlev,k) = cc3d_out(m,n,0:nlev,k)/counter_ave(m,n)
                      else
                         cc3d_out(m,n,0:nlev,k)=bio_missing
                      endif
                    enddo
                  enddo 
                case(14)
                  do n=jmin,jmax
                    do m=imin,imax
                      if (az(m,n) .ge. 1.and. counter_ave(m,n).gt.0 ) then
                           ccb3d_out(m,n,0:1,k) = ccb3d_out(m,n,0:1,k)/counter_ave(m,n)
                      else
                           ccb3d_out(m,n,0:1,k)=bio_missing
                      endif
                    enddo
                  enddo 
#ifdef INCLUDE_DIAGNOS_PRF
                case(15)
                  do n=jmin,jmax
                    do m=imin,imax
                      if (az(m,n) .ge. 1.and. counter_ave(m,n).gt.0 ) then
                        ccb3d_prf(m,n,0:nlev,k)=ccb3d_prf(m,n,0:nlev,k)/counter_ave(m,n)
                      else
                        ccb3d_prf(m,n,0:nlev,k)=bio_missing
                      endif
                    enddo
                  enddo 
#endif
              end select
            endif
          elseif ( (var_ids(i) /= 0) .and. ( .not.var_ave(i)) ) then 
            k=k+1
            select case(mode)
              case(11)
                do n=jmin,jmax
                  do m=imin,imax
                    if (az(m,n).lt.1 .or. counter_ave(m,n)==0 ) &
                          cc3d_out(m,n,0:nlev,k)=bio_missing
                  enddo
                enddo 
              case(14)
                do n=jmin,jmax
                  do m=imin,imax
                    if (az(m,n).lt.1 .or. counter_ave(m,n) ==0 ) &
                          ccb3d_out(m,n,0:1,k)=bio_missing
                  enddo
                enddo 
#ifdef INCLUDE_DIAGNOS_PRF
              case(15)
                do n=jmin,jmax
                  do m=imin,imax
                    if (az(m,n).lt.1 .or. counter_ave(m,n) ==0 ) &
                          ccb3d_prf(m,n,0:nlev,k)=bio_missing
                  enddo
                enddo 
#endif
            end select
          endif
        enddo
      case default 
        j=0
        do i=from,to
          j=j+1
          if ( (var_ids(i) /= 0 ).and. (var_ave(i).eqv.ave) ) then
            k=k+1
            if ( ( counter_ave(ig,jg).gt.0 )  .and.  ave) then
            ! add value for average output 
               select case(mode)
               case(1)
!JM                    cc3d_out(ig,jg,0:nlev,k) = cc3d_out(ig,jg,0:nlev,k)+cc(j,0:nlev)
                    cc3d_out(ig,jg,0:nlev,k) = cc3d_out(ig,jg,0:nlev,k)+cc(0:nlev,j)
               case(2)
!JM                    cc3d_out(ig,jg,0:nlev,k) = cc3d_out(ig,jg,:,k)+diag(j,0:nlev)
                    cc3d_out(ig,jg,0:nlev,k) = cc3d_out(ig,jg,:,k)+diag(0:nlev,j)
               case(3)
                   call make_flux_output(1,j,0,nlev,dt,c1dimz,llcalc)
                   if (llcalc) cc3d_out(ig,jg,0:nlev,k) = cc3d_out(ig,jg,0:nlev,k)+ c1dimz(0:nlev)
               case(4)
!JM                    ccb3d_out(ig,jg,0:1,k) = ccb3d_out(ig,jg,0:1,k)+ccb(j,0:1)
                    ccb3d_out(ig,jg,0:1,k) = ccb3d_out(ig,jg,0:1,k)+ccb(0:1,j)
               case(5)
!JM                    ccb3d_out(ig,jg,0:1,k) = ccb3d_out(ig,jg,0:1,k)+diagb(j,0:1)
                    ccb3d_out(ig,jg,0:1,k) = ccb3d_out(ig,jg,0:1,k)+diagb(0:1,j)
               case(6)
                   call make_flux_output(2,j,0,nlev,dt,c1dimz,llcalc)
                   if (llcalc) ccb3d_out(ig,jg,0:1,k) = ccb3d_out(ig,jg,0:1,k)+ c1dimz(0:1)
#ifdef INCLUDE_DIAGNOS_PRF
               case(7)
!JM                    ccb3d_prf(ig,jg,0:nlev,k) = ccb3d_prf(ig,jg,0:nlev,k)+diagb_prf(j,0:nlev)
                    ccb3d_prf(ig,jg,0:nlev,k) = ccb3d_prf(ig,jg,0:nlev,k)+diagb_prf(0:nlev,j)
#endif
               end select
            else
            ! set value for normal output or reset value for average output on zero 
            ! when counter_ave ==0 . 
               select case(mode)
               case(1)
!JM                    cc3d_out(ig,jg,0:nlev,k) = cc(j,0:nlev)
                    cc3d_out(ig,jg,0:nlev,k) = cc(0:nlev,j)
               case(2)
!JM                    cc3d_out(ig,jg,0:nlev,k) = diag(j,0:nlev)
                    cc3d_out(ig,jg,0:nlev,k) = diag(0:nlev,j)
               case(3)
                   call make_flux_output(1,j,0,nlev,dt,c1dimz,llcalc)
                   if ( llcalc) cc3d_out(ig,jg,0:nlev,k) = c1dimz(0:nlev)
               case(4)
!JM                    ccb3d_out(ig,jg,0:1,k) = ccb(j,0:1)
                    ccb3d_out(ig,jg,0:1,k) = ccb(0:1,j)
!                  ccb3d_out(k,ig,jg,0:1) = reshape(ccb(0:1,j), (/numbc,2/), order=(/2,1/))
               case(5)
                 !
                    print *, "ig:", ig, " jg:", jg, " k:", k, " j:", j
                    print *, "ccb3d_out dimensions:", size(ccb3d_out, 1), size(ccb3d_out, 2), size(ccb3d_out, 3), size(ccb3d_out, 4)
                    print *, "diagb dimensions:", size(diagb, 1), size(diagb, 2)

!JM                    ccb3d_out(ig,jg,0:1,k) = diagb(j,0:1)
                    ccb3d_out(ig,jg,0:1,k) = diagb(0:1,j)
               case(6)
                   call make_flux_output(2,j,0,nlev,dt,c1dimz,llcalc)
                   if ( llcalc) ccb3d_out(ig,jg,0:1,k) = c1dimz(0:1)
#ifdef INCLUDE_DIAGNOS_PRF
               case(7)
!JM                    ccb3d_prf(ig,jg,0:nlev,k) = diagb_prf(j,0:nlev)
                    ccb3d_prf(ig,jg,0:nlev,k) = diagb_prf(0:nlev,j)
#endif
               end select
            endif
          endif
        enddo
      end select
#endif
  end subroutine fill_diag_array
!EOC
!-----------------------------------------------------------------------

!BOP
!
! !ROUTINE: init_pel_co2
!
! !INTERFACE:
     subroutine init_pel_co2(mode,n,k)
!
! !USES:
#ifdef BFM_GOTM
    use mem, only:ETW,ESW,ERHO,InitializeModel,ppO3c
    use variables_3d, only: T,S,hn,rho
    use variables_bio_3d, only: ccb3d,cc3d
#endif

! !INPUT PARAMETERS:
    implicit none
    integer,intent(IN)          ::mode
    integer,intent(IN)          ::n
    integer,intent(IN)          ::k
!EOP
!-------------------------------------------------------------------------
!BOC

#ifdef BFM_GOTM
     integer                    :: i
     integer                    :: j
     REALTYPE                   :: cc(1:k,1:n)

     if ( mode >= 1 ) then
      Initializemodel=1
      do j=jmin,jmax
         do i=imin,imax
           if (az(i,j) .ge. 1 ) then
             ETW=T(i,j,1:kmax)
             ESW=S(i,j,1:kmax)
             ERHO=rho(i,j,1)
             cc(1:k,1:n)=cc3d(i,j,1:k,1:n)
             call CalcCO2SatInField(k,n,ERHO,ETW,ESW,cc);
             cc3d(i,j,1:k,ppO3c)=cc(1:k,ppO3c)
           endif
         enddo
       end do
       Initializemodel=0
     endif
     return
#endif
     end subroutine init_pel_co2
!-----------------------------------------------------------------------

!BOP
!
! !ROUTINE: getm_bfm_bennut_calc_initial
!
! !INTERFACE:
     subroutine getm_bfm_bennut_calc_initial(calc_init_bennut_states,read_poro0,hotstart_bio)
!
! !USES:
#ifdef BFM_GOTM
    use bio_var, only: numc,ccb,ppb,ddb,numbc,bio_setup
    use mem, only:ETW,ESW,ETAUB,Depth,ppESS,InitializeModel
    use mem_Param, ONLY:  p_poro
    use bio_bfm, only:reset_diagonal
    use variables_3d, only: T,S,hn
    use variables_bio_3d, only: ccb3d
    use gotm_error_msg, only: set_d3_model_flag_for_gotm, output_gotm_error,get_warning_for_getm
#endif
!
!
!
! !INPUT PARAMETERS:
    implicit none
      integer,intent(IN)          ::calc_init_bennut_states
      logical,intent(IN)          ::read_poro0
      logical,intent(IN)          ::hotstart_bio
#ifdef BFM_GOTM
!
! !DESCRIPTION:  read 2d-field for initialization of benthic states vars.
!
! !REVISION HISTORY:
!
!  15072006 Piet Raurdij
!
!EOP
!-------------------------------------------------------------------------
!BOC
     integer                    :: i
     integer                    :: j
     logical                    :: warning_flag
     LOGICAL                    :: error_flag
     LOGICAL                    :: r_p
     REALTYPE                   :: r


     if ( (.not.hotstart_bio).and.bio_setup /=2 ) then
           r_p=read_poro0
           call init_pel_co2(calc_init_bennut_states,numc,kmax)
     endif

     call output_2d_grid_parameters(0)

     if ( (.not.hotstart_bio) .and. calc_init_bennut_states >= 1 ) then
        if (calc_init_bennut_states ==1 ) then
           LEVEL4 "Benthic nutrient State vars are calculated with an assumption of an"
           LEVEL4 "equilibrium between input (nutrient regeneration) and output"
           LEVEL4 "(flux to water column and definitive loss processes)"
        endif
        if (calc_init_bennut_states >=1 ) then
           LEVEL4 "Benthic inorganic carbon State vars are initalized with an assumption of an"
           LEVEL4 "equilibrium between input (respiration) and output (flux to the air)"
        endif
        LEVEL2 'checking benhtic system...'
        do j=jmin,jmax
           do i=imin,imax
              if (az(i,j) .ge. 1 ) then
                 ! Only the the value of the lowest layer is needed
                 call set_2d_grid_parameters(r_p,igrid=i,jgrid=j)
                 Depth(1)=hn(i,j,1)
                 ETW(1)=T(i,j,1)
                 ETAUB(1)=0.003
                 if ( Depth(1) <=_ZERO_) then
                   write(msg,'(''ETW='',F10.2,'' Depth='',F10.2,'' '',I4,I4,I4)') &
                                    ETW(1),Depth(1),i+ioff,j+joff,az(i,j)
                   LEVEL4 msg
                 endif
                 ccb(:,1:numbc)=ccb3d(i,j,:,1:numbc)
!JM above compiles but crashes                 ccb(:,1:numbc)=reshape(ccb3d(1:numbc,i,j,:),(/numbc,1,1,2/), order=(/4,2,3,1/))
!                call InitBenthicNutrientDynamics(calc_init_bennut_states)
                 ccb3d(i,j,:,1:numbc)=ccb(:,1:numbc)
!JM above compiles but crashes                 ccb3d(1:numbc,i,j,:)=reshape(ccb(:,1:numbc), (/numbc,2/), order=(/2,1/))
                 call reset_diagonal(numbc,ppb)
                 call reset_diagonal(numbc,ddb)
                 call get_warning_for_getm(warning_flag )
                 if ( warning_flag ) then
                    r=sum(hn(i,j,:))
                    write(msg,'(''i,j,='',I4,''('',I4,'') '',I4,''('',I4,'') Depth='',F10.3)')&
                          i+ioff,i,j+joff,j,r
                   LEVEL3 msg
                 endif
                 call output_gotm_error( error_flag, sub, msg)
                 if ( error_flag ) then
                     LEVEL1 'i,j=',i+ioff,j+joff
                     call getm_error( sub, msg)
                 endif
              endif
           enddo
         end do
        LEVEL4 "end of benthic initialization------------------------------------"
      endif

     return
#endif
     end subroutine getm_bfm_bennut_calc_initial
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_2d_grid_parameters
!
! !INTERFACE:
 subroutine set_2d_grid_parameters(fields_present,file,bio_missing, &
                                                         igrid,jgrid,test_presence)
!
! !IUSES:
#ifdef BFM_GOTM
 use variables_bio_3d,ONLY: p_pK1_ae_2d,p_pK5_ae_2d,p_poro_2d,p_poro_default,&
                        p_pK1_ae_default,p_pK5_ae_default
 use mem_Param, ONLY:  p_pK1_ae,p_pK5_ae,p_poro
 use mem,ONLY:         ppK1p,ppK5s
#endif
!

! !INPUT PARAMETERS:
      implicit none
      logical,intent(IN)                        ::fields_present
      character(len=*),intent(IN),optional      ::file
      REALTYPE,intent(IN),optional              ::bio_missing
      integer,intent(IN),optional               ::igrid
      integer,intent(IN),optional               ::jgrid
      logical,intent(OUT),optional              ::test_presence

#ifdef BFM_GOTM
!
! !LOCAL VARIABLES
      integer                ::rc
      integer                ::nc=0
      integer                ::i
      integer                ::j
      integer                ::status
      logical                ::mode2
      REALTYPE               :: kd
      REALTYPE               :: rho=2.65
      REALTYPE               :: vr(2)

!
! !DESCRIPTION:
!         read porosities from netcdf-file
!         calculatest adsorption ceoof for Phosphate in aerobicx layer from
!         porosties
!
! !REVISION HISTORY:
!
!  01072006 Created by Piet Ruardij
!
!EOP
!-------------------------------------------------------------------------
!BOC

  if ( present(file)) then
    if (fields_present) then
         call getm_error('set_2d_grid_parameters', &
                                        'set_2d_grid_parameters routine already initialized')
         return
    endif
    ! initialization
    p_poro_default=p_poro(1)
    STDERR 'p_poro_default=',p_poro_default
    p_pK1_ae_default=p_pK1_ae(1)
    p_pK5_ae_default=p_pK5_ae(1)
    test_presence=.false.
    if (len_trim(file) .ne. 0 ) then
      inquire(file=file,exist=test_presence)
      if ( .not.test_presence) then
        call getm_error('set_2d_grid_parameters','error:file does not exist name='//file)
      endif
    endif
    allocate(p_poro_2d(E2DFIELD),stat=rc)
    if (rc /= 0) call getm_error('set_2d_grid_parameters',' Error allocating memory (p_poro_2d)')
    allocate(p_pK1_ae_2d(E2DFIELD),stat=rc)
    if (rc /= 0) call getm_error('set_2d_grid_parameters',' Error allocating memory (p_p_K1_ae)')
    allocate(p_pK5_ae_2d(E2DFIELD),stat=rc)
    if (rc /= 0) call getm_error('set_2d_grid_parameters',' Error allocating memory (p_p_K5_ae)')
    p_poro_2d=bio_missing;
    if (.not.test_presence ) then
       p_poro_2d=p_poro_default
    STDERR 'p_poro_default=',p_poro_default
       p_pK1_ae_2d=p_pK1_ae_default
       p_pK5_ae_2d=p_pK5_ae_default
       return
    else
      p_pK1_ae_2d=bio_missing
      p_pK5_ae_2d=bio_missing

       call init_2dbio_ncdf( 11,file,'Porosity',nc,status,p_poro_2d)
       if ( status.ne.0) call getm_error("set_2d_grid_parameters()",  &
                             'Could not find name in '//file )

       do j=jmin,jmax
          do i=imin,imax
             if (az(i,j) .ge. 1 ) then
                if ( p_poro_2d(i,j) > 0 ) then
                   p_poro_2d(i,j)=max(0.39D+00,p_poro_2d(i,j))
                else
                   p_poro_2d(i,j)=0.39D+00
                endif
                 call  CalcAdsorptionFromPoro(ppK1p,p_poro_2d(i,j),p_pK1_ae_2d(i,j))
                 call  CalcAdsorptionFromPoro(ppK5s,p_poro_2d(i,j),p_pK5_ae_2d(i,j))
             else
                p_poro_2d(i,j)=bio_missing
             endif
          enddo
       enddo
!      STDERR  'p_poro',p_poro_2d(imin:imax,jmin:jmax)
!      STDERR  'az',az(imin:imax,jmin:jmax)
    endif
  elseif ( present(igrid) .and. present(jgrid) ) then
     if ( fields_present ) then
        if ( p_poro_2d(igrid,jgrid) > 0 ) then
           p_pK1_ae(1)=p_pK1_ae_2d(igrid,jgrid)
           p_pK5_ae(1)=p_pK5_ae_2d(igrid,jgrid)
           p_poro(1)=p_poro_2d(igrid,jgrid)
           return
        endif
     else
        p_pK1_ae(1)=p_pK1_ae_default
        p_pK5_ae(1)=p_pK5_ae_default
        p_poro(1)=p_poro_default
        return
     endif
  else
     LEVEL3 'wrong use of set_2d_grid_parameters if call is NOT according:'
     LEVEL3 'call set_2d_grid_parameters(file=''filename'')'
     LEVEL3 'call set_2d_grid_parameters(igrid=i,jgrid=j)'
     call getm_error('set_2d_grid_parameters','')
  endif
#endif

 end subroutine set_2d_grid_parameters
!-------------------------------------------------------------------------
 subroutine output_2d_grid_parameters(mode)
!
! !IUSES:
#ifdef BFM_GOTM
     use variables_bio_3d,ONLY: p_pK1_ae_2d,p_pK5_ae_2d,p_poro_2d

     implicit none
     integer,intent(IN)   :: mode
     REALTYPE             :: vr(2)

#ifdef DEBUG
     STDERR 'do_output_2d_grid_paramters'
#endif

     vr(1)=_ZERO_;vr(2)=1.0
     call save_bfm_2d('p_poro_2d',p_poro_2d,vr,'-','porosity' )
     vr(1)=1.0;vr(2)=800.0
     call save_bfm_2d('p_pK1_ae',p_pK1_ae_2d,vr,'-','adsorption coeff,phosphate')
     call save_bfm_2d('p_pK5_ae',p_pK5_ae_2d,vr,'-','adsorption coeff,silicate')
#endif 
 end subroutine output_2d_grid_parameters

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_2d_grid
!
! !INTERFACE:
     subroutine init_2d_grid(file)
!
! !USES:
#ifdef BFM_GOTM
    use string_functions, only: empty,getseq_number
!JM    use bio_var, only: var_names,stBenStateS,stBenStateE
    use bfm_output, only: var_names,stBenStateS,stBenStateE
    use mem, only:ETW, Depth
    use variables_bio_3d, only: ccb3d
    use string_functions, only: empty,getseq_number
#endif

!
! !INPUT PARAMETERS:
     implicit none
     character(len=*)              ::file
#ifdef BFM_GOTM
!
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:  read 2d-field for initialization of benthic states vars.
!
! !REVISION HISTORY:
!
!  15072006 Piet Ruardij
!
!EOP
!-------------------------------------------------------------------------
!BOC
     integer                    ::i,k,m,nn
     integer                    ::l
     integer                    ::j
     integer                    ::nc=0
     integer                    ::status=0
     integer, dimension(:), allocatable :: flag_read
     character(len=10)          ::text
     character(len=80)          ::out
     REALTYPE                   :: ffb(I2DFIELD)


     if ( empty(file) ) return
     allocate(flag_read(1:stBenStateE-stBenStateS+1),stat=k)   ! pel. biological fields of diagnos.
     if (k /= 0) stop 'init_2d-grid: Error allocating memory (flag_read)'
     flag_read=0
     l=1
     j=0
     LEVEL1 '---Benthic Initialization from: ',file,':'
     do i=stBenStateS,stBenStateE
        j=j+1
        call init_2dbio_ncdf(l,file,var_names(i),nc,status,ffb)
        if ( status ==0 ) then
          flag_read(j)=1
          ccb3d(:,:,1,j)=ffb
        endif
        l=0
     enddo
     call init_2dbio_ncdf(-1,file,'',nc,status,ffb)
     j=0
     do i=stBenStateS,stBenStateE
       j=j+1
       if ( flag_read(j) ==1 ) then
          k=len_trim(var_names(i))-1
          text=var_names(i)(1:k)
          if ( getseq_number(text,var_names,stBenStateE,.FALSE.,i+1).gt.0) then
            write(text,'(A,A1)') var_names(i)(1:k),'n'
            m= getseq_number(text,var_names,stBenStateE,.TRUE.,stBenStateS)
            if ( m> 0 ) then
               m=m-stBenStateS+1
               if ( flag_read(m)==0 ) then
                  flag_read(m)=2
                  ccb3d(:,:,1,m)= ccb3d(:,:,1,j) * 0.01275
                  write(out,'(A,'' derived from '',A)') text ,var_names(i)(1:k+1)
                  nn=len_trim(out);LEVEL3 out(1:nn)
               endif
            endif
            write(text,'(A,A1)') var_names(i)(1:k),'p'
            m= getseq_number(text,var_names,stBenStateE,.TRUE.,stBenStateS)
            if ( m> 0 ) then
               m=m-stBenStateS+1
               if ( flag_read(m)==0 ) then
                  flag_read(m)=2
                  ccb3d(:,:,1,m)= ccb3d(:,:,1,j) * 0.07862e-3
                  write(out,'(A,'' derived from '',A)') text ,var_names(i)(1:k+1)
                  nn=len_trim(out);LEVEL3 out(1:nn)
               endif
            endif
            write(text,'(A,A1)') var_names(i)(1:k),'s'
            m= getseq_number(text,var_names,stBenStateE,.TRUE.,stBenStateS)
            if ( m> 0 ) then
               m=m-stBenStateS+1
               if ( flag_read(m)==0 ) then
                  flag_read(m)=2
                  ccb3d(:,:,1,m)= ccb3d(:,:,1,j) * 0.0145
                  write(out,'(A,'' derived from '',A)') text ,var_names(i)(1:k+1)
                  nn=len_trim(out);LEVEL3 out(1:nn)
               endif
            endif
          endif
       endif
     enddo
     j=0
     m=0
     do i=stBenStateS,stBenStateE
        j=j+1
        if ( flag_read(j) ==0 ) then
          if ( m==0 ) then
           LEVEL2 'Warning:'
           LEVEL2 'the following vars. were not directly/indirectly initialized from file:'
           m=1
          endif
          k=len_trim(var_names(i))
          LEVEL3 var_names(i)(1:k)
        endif
     enddo
     deallocate(flag_read)

#endif
     end subroutine init_2d_grid

  logical function  check_reset_tracking(mode,nr,ju,jutr) !AN
#ifdef BFM_GOTM
     use mem, ONLY: ii3dTrack,ii2dTrack,iiTrack
     use bio_var, only: numc,numbc
#endif
     implicit none 
     integer,intent(IN)      ::mode
     integer,intent(IN)      ::nr
     integer,intent(IN)      ::ju
     integer,intent(IN)      ::jutr

#ifdef BFM_GOTM
     if (jutr<ju.or.iiTrack==0) then
       check_reset_tracking=.false.
     elseif ( mode ==2 ) then
       check_reset_tracking=nr.gt.numbc-ii2dTrack    
     elseif ( mode ==3 ) then
       check_reset_tracking=nr.gt.numc-ii3dTrack    
     else
       check_reset_tracking=.false.
     endif
#else
       check_reset_tracking=.false.
#endif
  end function check_reset_tracking

  integer function  unlabeled_var_index(mode,nr)
#ifdef BFM_GOTM
     use mem, ONLY: ii3dTrack,ii2dTrack,iiTrack,nr_2d_track,nr_3d_track
     use bio_var, only: numc,numbc
#endif
     implicit none 
     integer,intent(IN)      ::mode
     integer,intent(IN)      ::nr
     
#ifdef BFM_GOTM
     integer                 ::out

     out=nr
     if ( iiTrack==0 ) then
     elseif ( mode ==2 ) then
       if (nr.gt.numbc-ii2dTrack) out=nr_2d_track(nr)  
     elseif ( mode ==3 ) then
       if (nr.gt.numc -ii3dTrack) out=nr_3d_track(nr)    
     endif
     if ( out.lt.0) out=nr
     unlabeled_var_index=out
     return
#else
     unlabeled_var_index=0
#endif
  end function unlabeled_var_index

!  integer function  assign_river_inputs_to_bio_states(ncid,river_name,river_nr,state_nr,iout)
!#ifdef BFM_GOTM
!     use netcdf
!     use string_functions, ONLY: getseq_number
!     use bio_var, only: numc,var_names
!     use mem, ONLY: ii3dTrack,nr_3d_track,iiTrack
!     use trace_bdy,only: trace
!     use time, only:JulDay
!#endif
!
!
!     implicit none
!     integer,intent(IN)           ::ncid
!     integer,intent(IN)           ::river_nr
!     integer,intent(IN)           ::state_nr
!     integer,intent(OUT)          ::iout
!     character(len=*),intent(IN)  ::river_name
!
!#ifdef BFM_GOTM
!     integer                                    :: nr_real_states
!     integer                                    :: unit=311
!     integer                                    :: rc
!     integer                                    :: n
!     integer                                    :: m
!     integer                                    :: err=0
!     integer                                    :: nriver=0
!     character(len=128)                         :: name=""
!     character(len=20)                          :: river_part
!     REAL_4B                                    :: att_wrk
!     REALTYPE                                   :: lon_river
!     REALTYPE                                   :: lat_river
!     logical                                    :: use_tracer_type=.false.
!     logical                                    :: exist=.false.
!     integer                                    :: yy,mm,dd
!     character(len=1)                           :: c1,c2
!  
!     namelist /trace_nml/ tracer_type
!
!     nr_real_states=numc
!     river_part="";if (len_trim(river_name)> 0 .and.river_nr.gt.0 ) river_part =trim(river_name)//'_'
!     if ( iiTrack.ge.1) nr_real_states=numc-ii3dTrack;
!     if ( ncid==0 ) then
!       if ( index(var_names(1),'@riv_')>0 ) trace=.true.
!       if (trace) then
!         inquire(file='tracer_info.dat',exist=exist)
!         tracer_type=''
!         if (exist) then
!              !open(unit,file='tracer_info.dat',action='read',status='old',iostat=err)
!              if ( err.ne.0) &
!                 call getm_error("assign_river_inputs_to_bio_states","error when opening tracer_info.dat")
!              read(unit,nml=trace_nml,iostat=err)
!              if ( err.ne.0) &
!                 call getm_error("assign_river_inputs_to_bio_states", &
!                                        "error when reading namelist from tracer_info.dat")
!              close(unit)
!         endif
!         use_tracer_type=len_trim(tracer_type).ne.0
!         LEVEL3 'TRACER_RUN'
!         LEVEL3 'tracer_type=',tracer_type
!       elseif ( iiTrack.ge.1.and.ii3dTrack.gt.1 ) then
!         LEVEL3 'TRACK_RUN'
!         open(unit,file='track_info.dat',action='read',status='old',iostat=err)
!         if ( err.ne.0) &
!           call getm_error("assign_river_inputs_to_bio_states","error when opening track_info.dat")
!         read(unit,*) name
!         ! get a date from the first line in track_info.dat
!         if( scan(name,'/-') .eq.0) &
!          call getm_error("assign_river_inputs_to_bio_states","No valid date on first line of track_info.dat") 
!          ! reclaculate the date in julian days  and put in integer_var n
!          read(name,'(i4,a1,i2,a1,i2)') yy,c1,mm,c2,dd; call JulDay(yy,mm,dd,n)
!         start_tracking_in_jul=n 
!         read(unit,*) nriver
!         if ( nriver >0 ) then
!           allocate(tracked_river_name(nriver),stat=rc) ! NetCDF name of river which will be tracked
!           if (rc /= 0) stop 'assign_river_inputs_to_bio_states:Error allocating memory (tracked_river_name)'
!           do n=1,nriver
!              read(unit,*) tracked_river_name(n)
!              LEVEL4 "A nutrient or the tracer in river ",trim(tracked_river_name(n))," will be tracked"
!           enddo
!         else
!           read(unit,*) lon_min,lat_min
!           read(unit,*) lon_max,lat_max
!           write(name,1000) lon_min,lat_min,lon_max,lat_max
! 1000      format('All rivers in the area ',F5.1,',',F6.1,' to ',F5.1,',',F6.1,' will be tracked')
!           LEVEL4  trim(name)
!         endif
!         close(unit)
!       endif
!     elseif ( iiTrack.eq.0.or.state_nr<=nr_real_states) then
!        iout=-9999
!        if (.not. trace ) then
!          name=trim(river_part)//trim(var_names(state_nr))
!        else
!          n=state_nr-river_nr;if ( n.eq.0) iout=0;
!          if (n .eq.0.or.(.not.use_tracer_type)) goto 90 ; !OK: no tracer  variable
!          name=trim(river_part)//trim(tracer_type)
!        endif
!        err =  nf90_inq_varid(ncid,trim(name),iout)
!        if (err .ne. NF90_NOERR) iout = -9999
!        if ( iout.ne. -9999 ) LEVEL4 trim(river_name),': ',trim(var_names(state_nr))
!     elseif ( ii3dTrack>1) then
!        ! tracking only-----------------------------------------
!        m=0;n=nr_3d_track(state_nr) ;iout=-2; 
!        ! n> 0 : tracked state variable, n =-1: tracer.
!        if ( n.gt. 0 .or.n.eq.-1 ) then
!          if (allocated(tracked_river_name)) then
!            ! check if river is track if this is the case  m will 
!            !     get the sequence number of the river.
!            m=size(tracked_river_name)
!            m=getseq_number(river_name,tracked_river_name(1:m),m,.TRUE.)
!            err=NF90_NOERR
!          elseif (lat_min> -9998.0) then
!            ! if river point is in the area  m will get the value 1
!            err =  nf90_inq_varid(ncid,trim(river_name),m)
!            if (err .eq. NF90_NOERR) then
!              err=nf90_get_att(ncid,m,'longitude',att_wrk);
!              if (err .eq. NF90_NOERR) then
!                lon_river=att_wrk
!                err=nf90_get_att(ncid,m,'latitude',att_wrk);
!                if (err .eq. NF90_NOERR) then
!                  lat_river=att_wrk;m=0
!                  if (lat_river-lat_max< -1.0D-6 .and. lat_river-lat_min >1.0D-6 &
!                  .and. lon_river-lon_max< -1.0D-6 .and. lon_river-lon_min >1.0D-6 ) m=1
!                endif
!              endif
!            endif
!          endif
!        endif
!        if ( err.eq.NF90_NOERR) then
!          if (m>0) then
!            if (n>0) then
!              ! tracked state variable
!              name=trim(river_part)//trim(var_names(n))
!              err =  nf90_inq_varid(ncid,trim(name),iout)
!              if ( err.eq.NF90_NOERR) then
!                ! there is an timeseries found which will be used also for the tracked variable
!                write(name,'(A,''-'',A,'' in '',A,'' will be tracked by adding this input to '',A)') &
!                 trim(river_name), trim(var_names(n)),trim(river_name),trim(var_names(state_nr))
!                 LEVEL4 trim(name) 
!              else
!               ! there is a tracked river found however no timeseries
!               ! instead of a timeseries of the limits will be used which are used for the
!               ! the full non_traced state variable.
!               iout=-1
!              endif
!            elseif ( n.eq.-1) then
!              ! sepecial river tracer state variable (tr_t)
!              ! iout=0 controls river input such that the conv. of the river is set on zero.
!              write(name,'(''Water from '',A,'' will be traced by setting the input for '',A,'' on 1'')') &
!                   trim(river_name), trim(river_name)
!              LEVEL4 trim(name) 
!              iout=0
!            endif
!          endif
!        endif
!     endif
!
!  90  assign_river_inputs_to_bio_states=err;
!#else
!     iout=0
!     assign_river_inputs_to_bio_states=0;
!#endif
!  end function assign_river_inputs_to_bio_states

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
      subroutine  make_nettransport_flux_output(mode,dt,iistate,from,to,ave,k,iicheck)
    
!
! !USES:
#ifdef BFM_GOTM
     use variables_bio_3d,only: cc3d,ffp,ccb3d_out,counter_reset
     use variables_3d,only: hn
     use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff   !getm
     use domain, only: az                                       !getm
!JM     use bio_var,only: var_ids,var_ave                          !bfm
     use bfm_output,only: var_ids,var_ave                          !bfm
     use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &  !bfm
              flx_ostates,flx_SS,flx_cal_ben_start,flx_option   !bfm
#endif
!
! !INPUT PARAMETERS:
     implicit none
     integer,intent(IN)          :: mode
     REALTYPE,intent(IN)         :: dt
     integer,intent(IN)          :: iistate
     integer,intent(IN)          :: from
     integer,intent(IN)          :: to
     logical,intent(IN)          :: ave
     integer,intent(INOUT)       :: k
     integer,intent(OUT)         :: iicheck
   
!
!
!
#ifdef BFM_GOTM
! LOCAL PARAMETERS:
      integer                                 :: i
      integer                                 :: j
      integer                                 :: l
      integer                                 :: n
      integer                                 :: m
      integer                                 :: z
      logical                                 :: ll_state_found
      REALTYPE                                :: h
      REALTYPE                                :: r
!
!
!EOP
!-------------------------------------------------------------------------
!BOC

      j=flx_cal_ben_start;
      l=10000;iicheck=0
      do i=from,to  !loop over the ben_glux vars to be stored.
       j=j+1
       if ( var_ids(i).ne.0) then
         k=k+1
         if (flx_option(j).ne.60.and.flx_option(j).ne.61 ) then
           !shortcut: no action
         elseif (mode.eq.0) then
           if (.not.counter_reset) return
           !initialization of output vars: happens at start/ after writing output.
           if (var_ave(i).eqv.ave ) then
              do n=jmin,jmax
              do m=imin,imax
                if (az(m,n).ge.1) ccb3d_out(m,n,:,k)=_ZERO_
              enddo
              enddo 
            endif
            l=min(l,flx_states(flx_calc_nr(j)))
            iicheck=1
         else
           ! filling of output vars
           ! filling hapens at every delta if ave==.true
           ! filling happens once if there is sampled only at outdelt ( ave=.false.)
           if (var_ave(i).eqv.ave ) then
             ll_state_found=.false.
             do l=flx_calc_nr(j-1)+1,flx_calc_nr(j)
               if (ll_state_found) then
               ! In this way conmparing indices is stopped as soon as one is found.
               elseif ( flx_states(l) == iistate) then
                 ll_state_found=.true.
                 do n=jmin,jmax
                 do m=imin,imax
                   if (az(m,n).ge.1) then
                      r=sum((ffp(m,n,1:kmax)-cc3d(m,n,1:kmax,iistate))*hn(m,n,1:kmax))/dt*86400.0 ! FIX: check hn and ffp index
                      if ( flx_option(j)==61 ) r=max(_ZERO_,r)
                      ccb3d_out(m,n,1,k)= ccb3d_out(m,n,1,k)+ r
                   endif
                 enddo
                 enddo 
               endif
             enddo
           endif
         endif
       endif
      enddo
      ! iicheck=1 : sum(s) on states* uv_fluxes present
      ! iicheck=2 : sum(s) on uv_fluxes ( water transport) present
      if ( mode ==0.and.iicheck==1.and.l==0 ) iicheck=2
#else
      iicheck=0
#endif
     end subroutine make_nettransport_flux_output
!-----------------------------------------------------------------------

!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
      subroutine  make_uv_flux_output(mode,iistate,from,to,ave,k,iicheck)
    
!
! !USES:
#ifdef BFM_GOTM
     use variables_bio_3d, only: cut,cvt,cc3d_out,counter_reset     !getm
     use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff   !getm
     use domain, only: dxu,dxv,dyu,dyv                          !getm
     use domain, only: az                                       !getm
!JM     use bio_var,only: var_ids,var_ave                          !bfm
     use bfm_output,only: var_ids,var_ave                          !bfm
     use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &  !bfm
              flx_ostates,flx_SS,flx_cal_ben_start,flx_option   !bfm
#endif
!
! !INPUT PARAMETERS:
     implicit none
     integer,intent(IN)          :: mode
     integer,intent(IN)          :: iistate
     integer,intent(IN)          :: from
     integer,intent(IN)          :: to
     logical,intent(IN)          :: ave
     integer,intent(INOUT)       :: k
     integer,intent(OUT)       :: iicheck
   
!
!
!
#ifdef BFM_GOTM
! LOCAL PARAMETERS:
      integer                                 :: i
      integer                                 :: j
      integer                                 :: l
      integer                                 :: n
      integer                                 :: m
      integer                                 :: z
      logical                                 :: ll_state_found
      REALTYPE                                :: h
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
!LEVEL1 "iistate",iistate
      j=0;l=10000;iicheck=0
      do i=from,to
       j=j+1
       if ( var_ids(i).ne.0) then
         k=k+1
         if (flx_option(j).ne.30.and.flx_option(j).ne.31.and.  &
           flx_option(j).ne.40.and.flx_option(j).ne.41 ) then
          !shortcut: no action
         elseif (mode.eq.0) then
           if (.not.counter_reset) return
           !initialization of output vars: happens at start/ after writing output.
           if (var_ave(i).eqv.ave ) then
             do n=jmin,jmax
             do m=imin,imax
               if (az(m,n).ge.1) cc3d_out(m,n,:,k)=_ZERO_
             enddo
             enddo 
           endif
           l=min(l,flx_states(flx_calc_nr(j)))
           iicheck=1
         else
           ! filling of output vars
           ! filling hapens at every delta if ave==.true
           ! filling happens once if there is sampled only at outdelt ( ave=.false.)
           if (var_ave(i).eqv.ave ) then
             ll_state_found=.false.
             do l=flx_calc_nr(j-1)+1,flx_calc_nr(j)
               if (ll_state_found) then
               ! In this way conmparing indices is stopped as soon as one is found.
               elseif ( flx_states(l) == iistate) then
                 ll_state_found=.true.
                 select case (flx_option(j))
                   case(30)
                     do n=jmin,jmax
                     do m=imin,imax
                       if (az(m,n).ge.1) then
                         do z=1,kmax
                           cc3d_out(m,n,z,k)=&
                             cc3d_out(m,n,z,k)+cut(m,n,z)*dyu(m,n)  
!if (m.eq.2 .and. n.eq.2 .and. z.eq.25) then
!LEVEL1 "m,n,z,k",m,n,z,k
!LEVEL1 "cut, dyu",cut(m,n,z),dyu(m,n)
!LEVEL1 "cc3d_out",cc3d_out(m,n,z,k)
!endif
                         enddo
                       endif
                     enddo
                     enddo 
                   case(31)
                     do n=jmin,jmax
                     do m=imin,imax
                       if (az(m,n).ge.1 )  then
                         do z=1,kmax
                           cc3d_out(m,n,z,k)=cc3d_out(m,n,z,k) &
                            +max(_ZERO_,cut(m,n,z) )*dyu(m,n))  
                         enddo
                       endif
                     enddo
                     enddo 
                   case(40)
                     do n=jmin,jmax
                     do m=imin,imax
                       if (az(m,n).ge.1)  then
                         do z=1,kmax
                            cc3d_out(m,n,z,k)=cc3d_out(m,n,z,k)+cvt(m,n,z)*dxv(m,n) 
!if (m.eq.2 .and. n.eq.2 .and. z.eq.25) then
!LEVEL1 "m,n,z,k",m,n,z,k
!LEVEL1 "cvt, dxu",cvt(m,n,z),dxv(m,n)
!LEVEL1 "cc3d_out",cc3d_out(m,n,z,k)
!endif
                         enddo
                       endif
                     enddo
                     enddo 
                   case(41)
                     do n=jmin,jmax
                     do m=imin,imax
                       if (az(m,n).ge.1 )  then
                         do z=1,kmax
                             cc3d_out(m,n,z,k)= &
                               cc3d_out(m,n,z,k)+max(_ZERO_,cvt(m,n,z) )*dxv(m,n))
                         enddo
                       endif
                     enddo
                     enddo 
                 end select
               endif
             enddo
           endif
         endif
       endif
      enddo
      ! iicheck=1 : sum(s) on states* uv_fluxes present
      ! iicheck=2 : sum(s) on uv_fluxes ( water transport) present
      if ( mode ==0.and.iicheck==1.and.l==0 ) iicheck=2
!if (iistate .ne. 0) then
!LEVEL1 "end make_uv_flux_output"
!stop 
!endif
#else
      iicheck=0
#endif
     end subroutine make_uv_flux_output
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
      subroutine  test_flux_output(mode,nr0,ave,iout)
    
!
! !USES:
#ifdef BFM_GOTM
     use mem, only: flx_cal_ben_start,flx_option
#endif
! !INPUT PARAMETERS:
     implicit none
     integer,intent(IN)          :: mode
     integer,intent(IN)          :: nr0
     logical,intent(IN)          :: ave
     integer,intent(OUT)         :: iout
#ifdef BFM_GOTM
! LOCAL PARAMETERS:
      integer                                 :: nr

      iout=0;nr=nr0
      if ( mode == 2 ) nr=nr+flx_cal_ben_start
      if ( flx_option(nr).ge.30.and.(.not.ave)) then
        STDERR 'Warning: a 3d-calculation (u,v,riv) is set in the var_save output section'
        STDERR 'Warning: instead of ave_save section.'
        STDERR 'Warning: Therefore output of these variables will suppressed'
        iout=1
      endif

#else
       iout=0
#endif
     end subroutine  test_flux_output
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!      subroutine  make_river_flux_output(height,dt,rSurface,rivcon,n,ig,jg)
!    
!
!! !USES:
!#ifdef BFM_GOTM
!     use variables_bio_3d, only: cut,cvt,ccb3d_out,counter_ave      !getm
!     use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff   !getm
!     use domain, only: dxu,dxv,dyu,dyv                          !getm
!     use domain, only: az                                       !getm
!     use bio_var,only: var_ids,var_ave,stBenFluxS,stBenFluxE    !bfm
!     use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &  !bfm
!              flx_ostates,flx_SS,flx_cal_ben_start,flx_option   !bfm
!     use time, only:secsprday
!#endif
!!
!! !INPUT PARAMETERS:
!     implicit none
!     REALTYPE,intent(IN)         :: height
!     REALTYPE,intent(IN)         :: dt
!     REALTYPE,intent(IN)         :: rSurface !reciproke of surface
!     REALTYPE,intent(IN)         :: rivcon(n)
!     
!     integer,intent(IN)          :: n
!     integer,intent(IN)          :: ig
!     integer,intent(IN)          :: jg
!!
!!
!!
!#ifdef BFM_GOTM
!! LOCAL PARAMETERS:
!      integer                                 :: i
!      integer                                 :: j
!      integer                                 :: l
!      integer                                 :: k
!      REALTYPE                                :: r
!      REALTYPE                                :: s
!!
!!
!!EOP
!!-------------------------------------------------------------------------
!!BOC
!      ! diag_end_sections is initiliazed in the call to fill_diagn_bfm_vars(0,......
!      ! this routine is called in init_getm_bio
!      ! OUTPUT of is written in a "benthic" var. (D2FIELD) .. (one number per river)
!      ! 
!      k=diag_end_sections(2,3)
!      j=flx_cal_ben_start;
!      do i=stBenFluxS,stBenFluxE
!       j=j+1
!       if ( var_ids(i)==0) then
!         !shortcut: no action
!       else 
!         k=k+1;
!         r=_ZERO_;
!         if (flx_option(j).eq.50.or.flx_option(j).eq.52) then
!           if ( dt.eq._ZERO_ ) then
!             s=_ZERO_
!           elseif ( flx_option(j).eq.50) then
!             ! Calculate vol per second 
!             s=height/(dt*rSurface)
!           else
!             ! Calculate  input per m3 per day
!             s=height/dt*secsprday
!           endif
!           ! filling of output vars
!           ! filling hapens at every delta if ave==.true
!           do l=flx_calc_nr(j-1)+1,flx_calc_nr(j)
!             select case (flx_states(l))
!               case (0)     ; r=s
!               case default ; r=r+s* rivcon(flx_states(l))
!             end select
!           enddo
!           select case (counter_ave(ig,jg))
!             case (0)     ; ccb3d_out(k,ig,jg,1)=r
!             case default ; ccb3d_out(k,ig,jg,1)=ccb3d_out(k,ig,jg,1)+r
!           end select
!         endif
!       endif
!      enddo
!#endif
!     end subroutine make_river_flux_output
!!EOC
!!-----------------------------------------------------------------------
!!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
   subroutine check_3d_track_rates(where)
! !USES:
#ifdef BFM_GOTM
   use bio_var, only: ccb,ppb,ddb,numbc,cc,numc
   use track, only: check_track_states
   use variables_bio_3d, only: cc3d
#endif

! !INPUT PARAMETERS:
   implicit none
   character*(*),intent(IN) :: where

! LOCAL PARAMETERS:
#ifdef BFM_GOTM
   integer          :: i
   integer          :: j
   REALTYPE         :: cc_copy(1:kmax,1:numc)
 
!EOP
!-------------------------------------------------------------------------
!BOC
   cc_copy=cc
   do i=imin,imax
       do j=jmin,jmax
          if (az(i,j) .ge. 1 )  then
            cc =cc3d(i,j,:,:)
            call check_track_states(-1,where)
          endif
       enddo
   enddo 
   cc=cc_copy
   return
#endif
   end subroutine check_3d_track_rates
!EOP
!-------------------------------------------------------------------------
!BOC


end module coupling_getm_bfm

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2008 - BFM                                              !
!-----------------------------------------------------------------------

