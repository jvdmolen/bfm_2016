!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: output_restart_bio
!
! !INTERFACE:
   module output_restart_bio
!
! !BFM:
!      routine to make a hotstart initialization of a model which include BFM
!      - assign all state variables with initial values derived from an netcdf-file
!      - do some checking.
!      - if state variable does not exist, initial values of 1d-gotm is taken.
!        (See bio_bfm.nml)
!        In this way old-restart files can be used if the model is extended
!        with new state variables and no complete reintialization is needed
!
! !USES:
    use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff,H
    use domain, only: az,au,av
    use exceptions, only:getm_error
    use coupling_getm_bfm,only:check_reset_tracking,read_poro,&
                  set_2d_grid_parameters, start_tracking_in_jul,init_pel_co2

    character(len=80),private                 :: msg="",sub
    REALTYPE,parameter                        :: pretty_small=1.0D-20

    public restart_file_bio

!EOP
!-----------------------------------------------------------------------

   contains
!-----------------------------------------------------------------------

!!-------------------------------------------------------------------------
!EOC
!
! !ROUTINE: restart_file_bio
!
! !INTERFACE:
     subroutine restart_file_bio(mode,fmt_file,iounit,julianday)
!
! !USES:
#ifdef BFM_GOTM
     use bio, only: bio_calc
!JM     use bio_var, only: numc,numbc,bio_setup,var_names,cc,ccb
!JM     use bio_var, only:stPelStateS,stBenStateS,stPelStateE,stBenStateE
     use bio_var, only: numc,numbc,bio_setup,cc,ccb
     use bfm_output, only:var_names,stPelStateS,stBenStateS,stPelStateE,stBenStateE
     use variables_bio_3d,  only: cc3d,ccb3d,ffp,ffb,bio_missing
     use variables_3d,  only: T,S
     use string_functions, ONLY: getseq_number
!AN     use track, only: check_ben_states_on_neg_values_for_track, &
!AN                      check_track_states
     use mem,ONLY:iiMicroZooPlankton,ppMicroZooPlankton,iiMesoZooPlankton, &
     ppMesoZooPlankton,ppBenOrganisms,iiBenOrganisms,ppBenBacteria,iiBenBacteria
     use mem,only: iiC,iiN,iiP,iiPel,iiBen,ppO3c,ppO3h,ppG3c,ppG13c,ppG23c, &
       ppG3h,ppG13h,ppG23h
!JM     use mem,only: ppQ6c,ppQ6n,ppQ6p,ppQ6s,ppQ16c,ppQ16n,ppQ16p,ppQ16s, &
!JM                   ppD6m,ppD7m,ppD8m,ppD9m,ppD1m,ppDtm
     use mem,only: ppQ6c,ppQ6n,ppQ6p,ppQ6s, & !JM not in bfm2016 ppQ16c,ppQ16n,ppQ16p,ppQ16s, &
                   ppD1m!AN ,ppDtm
     use mem,only: ppHNn,ppHnp
     use mem,only: ppD1m,ppD2m,ppP6c,ppPCc
     use mem,only: ppDcm,ppDfm
     use mem_Param, ONLY: p_poro,p_d_tot,p_clDxm
     use mem_Bioturbation,ONLY:p_cturm
#endif
     use halo_zones, only: update_2d_halo,update_3d_halo,wait_halo
     use halo_zones, only: D_TAG
!
! !INPUT PARAMETERS:

     implicit none
     integer,intent(IN)      :: mode
     integer,intent(IN)      :: fmt_file
     integer,intent(IN)      :: iounit
     integer,intent(IN)      :: julianday
! !LOCAL PARAMETERS:

#ifdef BFM_GOTM
     integer                                     :: i
     integer                                     :: j
     integer                                     :: n,k,l
     integer                                     :: rc
     integer                                     :: numc_hotstart
     integer                                     :: numbc_hotstart
     character(len=64)                           :: string
     character(len=64),dimension(:),allocatable  :: var_names_hotstart
     integer                                     :: flag(numc)
     integer                                     :: flagb(numbc)
     REALTYPE                                    :: ffb_b(I2DFIELD,0:1)
     REALTYPE                                    :: cc_copy(0:kmax,1:numc)
     REALTYPE                                    :: ccb_copy(0:1,1:numbc)
     REALTYPE                                    :: ccb1d(0:1)
     REALTYPE                                    :: cc1d(0:kmax)
     REALTYPE                                    :: c_in_layer_1
!    REALTYPE                                    :: nxx,vxx
!    integer                                     :: ix,jx
     logical                                     :: reset,lldone,llcheck
!
!
!EOP
!-------------------------------------------------------------------------
     if (mode.eq.WRITING.and.fmt_file.eq.BINARY) then

       if(bio_calc) then
         LEVEL3 'saving bio variables'
         if ( bio_setup /= 2 ) then
           LEVEL4 'saving 3d/pelagic variabls'
           write(iounit) numc
           write(iounit) var_names(stPelStateS:stPelStateE)
           do i=1,numc
             write(iounit) cc3d(:,:,:,i)
           enddo
          endif
          if ( bio_setup  > 2 ) then
           LEVEL4 'saving 2d/benthic variabls'
           write(iounit) numbc
           write(iounit) var_names(stBenStateS:stBenStateE)
           do i=1,numbc
             write(iounit) ccb3d(:,:,:,i)
           enddo
          endif
       end if
     elseif (mode.eq.READING) then
       if (bio_calc) then
          LEVEL3 'reading bio variables'
          if ( bio_setup /= 2 ) then
              LEVEL3 'read 3d/pelagic variables'
              if ( fmt_file.eq.BINARY) then
                read(iounit) numc_hotstart
                allocate(var_names_hotstart(1:numc_hotstart),stat=rc)
                if (rc /= 0) STOP 'reading_file_bio: Error allocating (var_names_hotstart)'
                read(iounit) var_names_hotstart
                flag=0
                do i=1,numc_hotstart
                  read(iounit) ffp
                  j=getseq_number(var_names_hotstart(i), &
                              var_names(stPelStateS:stPelStateE),numc,.TRUE.)
                  reset=check_reset_tracking(3,j,julianday, &
                                                       start_tracking_in_jul)
                  if ( j>0.and.(.not.reset)) then
                     cc3d(:,:,:,j)=ffp;flag(j)=1
                  endif
                enddo
              elseif (fmt_file.eq.NETCDF) then
                call read_restart_bio_ncdf(1,0,numc_hotstart,string)
                allocate(var_names_hotstart(1:numc_hotstart),stat=rc)
                if (rc /= 0)  &
                  STOP 'reading_file_bio: Error allocating (var_names_hotstart)'
                do i=1,numc_hotstart
                   call read_restart_bio_ncdf(3,i,n,var_names_hotstart(i))
                enddo
                flag=0
                do i=1,numc_hotstart
                   call read_restart_bio_ncdf(5,i,n,string)
                   j=getseq_number(var_names_hotstart(i), &
                               var_names(stPelStateS:stPelStateE),numc,.TRUE.)
                   reset=check_reset_tracking &
                           (3,j,julianday,start_tracking_in_jul)
                   if ( j>0.and.(.not.reset)) then
                     cc3d(:,:,:,j)=ffp(:,:,:);flag(j)=1
                   endif
                enddo
              endif
              l=0;lldone=.false.
              do i=imin,imax
                do j=jmin,jmax
                  n=count(cc3d(i,j,1,1:numc)<  bio_missing+_ONE_)
                  if (az(i,j) .ge. 1.and. n>0 ) then
                    write(msg, &
                      '(A,'' i='',I3,'' ('',I3'') j='',I3,'' (''I3,'')'')') &
                    'missing intial pelagic data  gridpoint',i+ioff,i,j+joff,j
                    LEVEL3 trim(msg )
                    cc3d(i,j,:,1:numc)=cc(:,1:numc)
                    T(i,j,:)=bio_missing
                    S(i,j,:)=bio_missing ;l=l+1;lldone=.true.
                    LEVEL3 'Initial values defined in bio_bfm.nml are used'
                  endif
                enddo
              enddo
              if (.not.lldone) then
                do i=imin,imax
                do j=jmin,jmax
                  if (az(i,j) .ge. 1.and. S(i,j,1)< pretty_small ) then
                    write(msg, &
                      '(A,'' i='',I3,'' ('',I3'') j='',I3,'' (''I3,'')'')') &
                    'missing intial pelagic S and T gridpoint',i+ioff,i,j+joff,j
                    LEVEL3 trim(msg )
                    T(i,j,:)=bio_missing
                    S(i,j,:)=bio_missing ;l=l+1
                    LEVEL3 'Initial values of near grid points  are used'
                  endif
                enddo
                enddo
              endif
!JM it's not clear what this while loop adds as passing through this bit of code
!JM multiple times does not do anything more than passing once.
!JM also the loop can hang. So disabled
!JM              do while (l>0)
                do i=imin,imax
                do j=jmin,jmax
                  if (S(i,j,1)< bio_missing+_ONE_) then
                    lldone=.false.
                    do k=-1,1
                      do n=-1,1
                        if ((i+k.ge.1.and.j+n.ge.1).and. &
                        ((k.eq.0.and.n.ne.0).or.(n.eq.0.and.k.ne.0)) ) then
                        llcheck=(cc3d(i+k,j+n,1,1)>0.0) &
                                    .and.(S(i+k,j+n,1)>bio_missing+_ONE_)
                        if (llcheck.and.(.not.lldone)) then
                           T(i,j,:)=T(i+k,j+n,:)
                           S(i,j,:)=S(i+k,j+n,:)
                           lldone=.true. ;l=l-1
                        endif
                        endif
                      enddo
                    enddo
                  endif
                enddo
                enddo
!JM              enddo

              call calc_GrpNP(iiPel,cc,flag,iiC,iiN,iiMicroZooPlankton,&
                                            numc,kmax,ppMicroZooPlankton)
              call calc_GrpNP(iiPel,cc,flag,iiC,iiP,iiMicroZooPlankton,&
                                            numc,kmax,ppMicroZooPlankton)
              call calc_GrpNP(iiPel,cc,flag,iiC,iiN,iiMesoZooPlankton,&
                                            numc,kmax,ppMesoZooPlankton)
              call calc_GrpNP(iiPel,cc,flag,iiC,iiP,iiMesoZooPlankton,&
                                            numc,kmax,ppMesoZooPlankton)
              if ( ppPCc> 0 .and. ppP6c >0 ) then
                if ( flag(ppPcc)==0 .and. flag(ppP6c) == 1) then
                  LEVEL3 'No initial values present for ppPCc'
                  LEVEL3 'Make initial values for ppPCc equal to ppP6c'
                  do i=imin,imax
                     do j=jmin,jmax
                       if (az(i,j).ge.1)&
                         cc3d(i,j,0:kmax,ppPCc)= cc3d(i,j,0:kmax,ppP6c)
                     enddo
                  enddo
                  flag(ppPcc)=1
                endif
              endif
              if (ppo3c >0) then
                if (flag(ppO3c)==0) then
                  LEVEL3 'No initial values present for ppO3c'
                  LEVEL3 'Pelagic CO2-DIC intialized'
                  cc1d=cc(0:kmax,ppO3h)
                  do i=imin,imax
                     do j=jmin,jmax
                       if (az(i,j) .ge. 1 ) cc3d(i,j,0:kmax,ppO3h)=cc1d(0:kmax)
                     enddo
                  enddo
                  call init_pel_co2(2,numc,kmax) ; flag(ppO3c)=1; flag(ppO3h)=1
                endif
              endif
              do n=1,numc
                if (flag(n)==0 ) then
                  LEVEL3 'pelstates: No initial data present:', &
                       trim(var_names(n+stPelStateS-1)), &
                       ' initialized with default (1D)values'
                  cc1d=cc(0:kmax,n)
                  do i=imin,imax
                     do j=jmin,jmax
                       if (az(i,j) .ge. 1 ) cc3d(i,j,0:kmax,n)=cc1d(0:kmax)
                     enddo
                   enddo
                endif
              enddo
              do i=imin,imax
                do j=jmin,jmax
                  if (az(i,j)>0  .and. cc3d(i,j,1,1)<_ZERO_) then
                    LEVEL3 'cc3d(i,j,1,1)=',cc3d(i,j,1,1)
                    LEVEL3 'Apparently bathymetry changed for point:', i,j
                    LEVEL3 'Pelagic values for this point initialized with default (1D) values'
                    cc3d(i,j,:,:)=cc
                  endif
                enddo
              enddo
              do n=1,numc
                ffp(:,:,:)=cc3d(:,:,:,n)
                call update_3d_halo(ffp,ffp,az,imin,jmin,imax,jmax,kmax,D_TAG)
                call wait_halo(D_TAG)
                cc3d(:,:,:,n)=ffp(:,:,0:)
              enddo
              deallocate(var_names_hotstart)
          endif
          if ( bio_setup  > 2 ) then
            LEVEL3 'read 2d/benthic variables'
            if (fmt_file.eq.BINARY) then
              read(iounit) numbc_hotstart
              allocate(var_names_hotstart(1:numbc_hotstart),stat=rc)
              if (rc /= 0)  &
                STOP 'reading_file_bio: Error allocating (var_names_hotstart)'
              read(iounit) var_names_hotstart
              flagb=0
              do i=1,numbc_hotstart
                read(iounit) ffb_b
                j=getseq_number(var_names_hotstart(i),&
                           var_names(stBenStateS:stBenStateE),numbc,.TRUE.)
                reset=check_reset_tracking(2,j,julianday,start_tracking_in_jul)
                if ( j>0.and.(.not.reset)) then
                    ccb3d(:,:,:,j)=ffb_b ;flagb(j)=1
                endif
              enddo
            elseif (fmt_file.eq.NETCDF) then
              call read_restart_bio_ncdf(2,0,numbc_hotstart,string)
              allocate(var_names_hotstart(1:numbc_hotstart),stat=rc)
              if (rc /= 0) STOP  &
                   'reading_file_bio: Error allocating (var_names_hotstart)'
              do i=1,numbc_hotstart
                call read_restart_bio_ncdf(4,i,n,var_names_hotstart(i))
              enddo
              flagb=0
              do i=1,numbc_hotstart
                call read_restart_bio_ncdf(6,i,n,string)
                j=getseq_number(var_names_hotstart(i), &
                        var_names(stBenStateS:stBenStateE),numbc,.TRUE.)
                reset=check_reset_tracking(2,j,julianday,start_tracking_in_jul)
                if ( j>0.and.(.not.reset)) then
                  ccb3d(:,:,1,j)=ffb;flagb(j)=1
                endif
              enddo
            endif
            l=numbc/2
            do i=imin,imax
              do j=jmin,jmax
                n=count(ccb3d(i,j,1,1:numbc)<  bio_missing+_ONE_)
!JM                if (az(i,j) .ge. 1.and. n >l ) then      ! goes wrong for Dcm=(0,-9999)
                if (az(i,j) .ge. 1.and. n >0 ) then
                  write(msg, &
                    '(A,'' i='',I3,'' ('',I3'') j='',I3,'' (''I3,'')'')') &
                  'missing intial benthic data  gridpoint',i+ioff,i,j+joff,j
                  LEVEL3 trim(msg )
                  do n=1,numbc
                    ccb1d=ccb(0:1,n)
                    ccb3d(i,j,:,n)=ccb1d
                  enddo
                  LEVEL3 'Initial values defined in bio_bfm.nml are used'
                endif
              enddo
            enddo
            ccb_copy(0:1,1:numbc)=ccb(0:1,1:numbc)
            if ( ppG3c >0 ) then
              if ( flagb(ppG3c)==0) then
                LEVEL3 'No initial calues present for ppG3c'
                LEVEL3 'Benthic CO2-DIC initialized'
                do i=imin,imax
                  do j=jmin,jmax
                    if (az(i,j) .ge. 1 ) then
                      call set_2d_grid_parameters(read_poro,igrid=i,jgrid=j)
                      ccb(0:1,1:numbc)=ccb3d(i,j,0:1,1:numbc)
                      c_in_layer_1=cc3d(i,j,1,ppO3c)
                      ccb(1,ppG3c)=c_in_layer_1*p_poro(1)*ccb(1,ppD1m)
                      ccb(1,ppG13c)= &
                          c_in_layer_1*p_poro(1)*(ccb(1,ppD2m)-ccb(1,ppD1m))
                      ccb(1,ppG23c)= &
                          c_in_layer_1*p_poro(1)*(p_d_tot-ccb(1,ppD2m))
                      c_in_layer_1=cc3d(i,j,1,ppO3h)
                      ccb(1,ppG3h)=c_in_layer_1*p_poro(1)*ccb(1,ppD1m)
                      ccb(1,ppG13h)= &
                          c_in_layer_1*p_poro(1)*(ccb(1,ppD2m)-ccb(1,ppD1m))
                      ccb(1,ppG23h)= &
                          c_in_layer_1*p_poro(1)*(p_d_tot-ccb(1,ppD2m))
                      ccb3d(i,j,0:1,1:numbc)=ccb(0:1,1:numbc)
                    endif
                  enddo
                enddo
                flagb(ppG3c)=1; flagb(ppG3h)=1;flagb(ppG13c)=1;flagb(ppG13h)=1
                flagb(ppG23c)=1; flagb(ppG23h)=1
              endif
            endif
            rc=0
            do i=imin,imax
              do j=jmin,jmax
                if (az(i,j) .eq. 1 ) then
                  k=0
                  do n=1,numbc
                    if ( flagb(n).eq.1.and.ccb3d(i,j,1,n) < pretty_small) k=-1
                  enddo
                  if ( k<0) then
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiC,&
                               iiBenBacteria,numbc,1,ppBenBacteria,ig=i,jg=j)
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiN,&
                               iiBenBacteria,numbc,1,ppBenBacteria,ig=i,jg=j)
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiP,&
                                iiBenBacteria,numbc,1,ppBenBacteria,ig=i,jg=j)
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiC,&
                              iiBenOrganisms,numbc,1,ppBenOrganisms,ig=i,jg=j)
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiN,&
                              iiBenOrganisms,numbc,1,ppBenOrganisms,ig=i,jg=j)
                    call calc_GrpNP(iiBen,ccb_copy,flagb,iiC,iiP,&
                             iiBenOrganisms,numbc,1,ppBenOrganisms,ig=i,jg=j)
!                   STDERR 'output_restart: Hnn,Hnp', &
!                                &        ccb_copy(1,ppHNn),ccb_copy(1,ppHNp)
                    ccb(0:1,1:numbc)=ccb3d(i,j,0:1,1:numbc)
                  endif
                endif
              enddo
            enddo
!%JM not in bfm2016            if (flagb(ppQ16c)==0) then
!%JM              do i=imin,imax
!%JM                do j=jmin,jmax
!%JM                  if (az(i,j).ge.1) then
!%JM!JM                    ccb1d(:)=min(p_d_tot,max(p_clDxm,ccb3d(i,j,:,ppD6m)))
!%JM                    ccb1d(:)=min(p_d_tot,p_clDxm)
!%JM                    ccb1d(:)=ccb3d(i,j,:,ppQ6c)/(_ONE_-exp(-p_d_tot/ccb1d(:))) &
!%JM                                   *(_ONE_-exp(-ccb3d(i,j,:,ppD1m)/ccb1d(:)))
!%JM                    ccb3d(i,j,:,ppQ16c)=ccb3d(i,j,:,ppQ6c)-ccb1d(:)  
!%JM                    ccb3d(i,j,:,ppQ6c)=ccb1d(:)  
!%JM!                   ccb3d(ppD6m,i,j,:)=ccb3d(ppD6m,i,j,:)-ccb3d(ppD1m,i,j,:)
!%JM!JM                    ccb1d(:)=min(p_d_tot,max(p_clDxm,ccb3d(ppD7m,i,j,:)))
!%JM                    ccb1d(:)=min(p_d_tot,p_clDxm)
!%JM                    ccb1d(:)=ccb3d(i,j,:,ppQ6n)/(_ONE_-exp(-p_d_tot/ccb1d(:))) &
!%JM                                     *(_ONE_-exp(-ccb3d(i,j,:,ppD1m)/ccb1d(:)))
!%JM                    ccb3d(i,j,:,ppQ16n)=ccb3d(i,j,:,ppQ6n)-ccb1d(:)  
!%JM                    ccb3d(i,j,:,ppQ6n)=ccb1d(:)  
!%JM!                   ccb3d(ppD7m,i,j,:)=ccb3d(ppD7m,i,j,:)-ccb3d(ppD1m,i,j,:)
!%JM!JM                    ccb1d(:)=min(p_d_tot,max(p_clDxm,ccb3d(ppD8m,i,j,:)))
!%JM                    ccb1d(:)=min(p_d_tot,p_clDxm)
!%JM                    ccb1d(:)=ccb3d(i,j,:,ppQ6p)/(_ONE_-exp(-p_d_tot/ccb1d(:))) &
!%JM                                    *(_ONE_-exp(-ccb3d(i,j,:,ppD1m)/ccb1d(:)))
!%JM                    ccb3d(i,j,:,ppQ16p)=ccb3d(i,j,:,ppQ6p)-ccb1d(:)  
!%JM                    ccb3d(i,j,:,ppQ6p)=ccb1d(:)  
!%JM!                   ccb3d(ppD8m,i,j,:)=ccb3d(ppD8m,i,j,:)-ccb3d(ppD1m,i,j,:)
!%JM!JM                    ccb1d(:)=min(p_d_tot,max(p_clDxm,ccb3d(ppD9m,i,j,:)))
!%JM                    ccb1d(:)=min(p_d_tot,p_clDxm)
!%JM                    ccb1d(:)=ccb3d(i,j,:,ppQ6s)/(_ONE_-exp(-p_d_tot/ccb1d(:))) &
!%JM                                     *(_ONE_-exp(-ccb3d(i,j,:,ppD1m)/ccb1d(:)))
!%JM                    ccb3d(i,j,:,ppQ16s)=ccb3d(i,j,:,ppQ6s)-ccb1d(:)  
!%JM                    ccb3d(i,j,:,ppQ6s)=ccb1d(:)  
!%JM!                   ccb3d(ppD9m,i,j,:)=ccb3d(ppD9m,i,j,:)-ccb3d(ppD1m,i,j,:)
!%JM                  endif
!%JM                enddo
!%JM              enddo
!%JM              flagb(ppQ16c)=1;flagb(ppQ16n)=1;flagb(ppQ16p)=1;flagb(ppQ16s)=1;
!%JM            endif
!AN            if (flagb(ppDtm)==0) then
!AN              do i=imin,imax
!AN                do j=jmin,jmax
!AN                  if (az(i,j).ge.1) then
!AN                    ccb3d(ppDtm,i,j,:)=abs(p_cturm)  
!AN                  endif
!AN                enddo
!AN                flagb(ppDtm)=1
!AN              enddo
!AN            endif
            if ( ppDcm> 0 .and. ppDfm >0 ) then
              if ( flagb(ppDcm)==0 .and. flagb(ppDfm) == 1) then
                LEVEL3 'No initial values present for ppDcm'
                LEVEL3 'Make initial values for ppDcm equal to ppDfm'
                do i=imin,imax
                   do j=jmin,jmax
                     if (az(i,j).ge.1)&
                       ccb3d(i,j,:,ppDcm)= ccb3d(i,j,:,ppDfm)
                   enddo
                enddo
                flagb(ppDcm)=1
              endif
            endif
            do n=1,numbc
!               STDERR trim(var_names(stBenStateS-1+n)),flagb(n)
              if (flagb(n)==0 ) then
                LEVEL3 'ben_states: No initial data present: ', &
                 trim(var_names(n+stBenStateS-1)), &
                                        ' initialized with default (1D)values'
                do i=imin,imax
                  do j=jmin,jmax
                    if (az(i,j) .eq. 1 ) then
!                     call check_ben_states_on_neg_values_for_track &
!                                           (n,numbc,1,i,j,ccb_copy)
                      ccb1d(0:1)=ccb_copy(0:1,n)
                      ccb3d(i,j,0:1,n)=ccb1d(0:1)
                    endif
                  enddo
                enddo
              endif
            enddo
            do i=imin,imax
              do j=jmin,jmax
                if(az(i,j)==1  .and. ccb3d(i,j,1,1)<_ZERO_) then
                    LEVEL4 'Apparently bathymetry changed for point:', i,j
                    LEVEL4 'Benthic values for this initialized with default (1D) values'
                    ccb3d(i,j,:,:)=ccb
                endif
              enddo
            enddo
!            do n=1,numbc
            do i=imin,imax
              do j=jmin,jmax
                if(az(i,j)==2 ) then
                  ccb3d(i,j,0:1,:)=-9999.0
!                 if (ccb3d(n,i,j,1) < -9998.0) then
!                      ccb1d(0:1)=ccb_copy(n,0:1)
!                      ccb3d(n,i,j,0:1)=ccb1d(0:1)
!                 endif
                endif
              enddo
            enddo
            deallocate(var_names_hotstart)
          endif
          if (fmt_file.eq.NETCDF) call read_restart_bio_ncdf(7,0,n,string)
       end if
     endif
#endif

  end subroutine restart_file_bio

  subroutine calc_GrpNP(iiSys,cx,flag,iiC,iiN,iiM,n,nlev,ppM,ig,jg)
#ifdef BFM_GOTM
!JM  use bio_var, only: var_names,stBenStateS
  use bfm_output, only: var_names,stBenStateS
  use variables_bio_3d, only: cc3d,ccb3d
  use mem,only: iiPel,iiBen
#endif

  implicit none
  integer,intent(IN)            ::iiSys
  REALTYPE,intent(IN)           ::cx(0:nlev,1:n)
  integer,intent(INOUT)         ::flag(1:n)
  integer,intent(IN)            ::iiC
  integer,intent(IN)            ::iiN
  integer,intent(IN)            ::iiM
  integer,intent(IN)            ::n
  integer,intent(IN)            ::nlev
  integer,intent(IN),optional   ::ig
  integer,intent(IN),optional   ::jg


  INTERFACE                             ! Specification
   integer FUNCTION ppM(i,ll)          ! Specification
     integer,INTENT(IN)   ::i           ! Specification
     integer,INTENT(IN)   ::ll         ! Specification
   END FUNCTION                         ! Specification
  END  INTERFACE                        ! Specification

#ifdef BFM_GOTM
  integer                               :: i,j,k,l,plus
  REALTYPE                              :: r

  do i=1,iiM
    k=ppM(i,iiC)
    j=ppM(i,iiN)
    if ( j>0 ) then
      plus=0
      select case (iiSYS)
        case (iiPel)
          if ( flag(j)==0 ) then
            write(msg,'(A,'' was not initialized'')') trim(var_names(j+plus)) ; LEVEL4 trim(msg)
            write(msg,'(A,'' initialized assuming a fixed ../C quotum of '',G12.6)')  &
                 trim(var_names(j+plus)), cx(nlev,j)/(1.0D-80+cx(nlev,k))
            flag(j)=1
            forall (l=1:nlev) &
               cc3d(:,:,l,j)=cc3d(:,:,l,k)*cx(l,j)/(1.0D-80+cx(l,k))
             LEVEL4 trim(msg )
          endif
       case(iiBen)
         if ( flag(j)==1 ) then
          plus=stBenStateS-1
          r=1.0D+80
          do l=1,nlev
             r=min(r,ccb3d(ig,jg,l,j))
          enddo
!         STDERR "init_ben_orgranism",nlev,i,k,j,ig,jg
          if ( r<pretty_small) then
            write(msg,'(A,'' was smaller than '',G12.6,'' i='',I3,'' ('',I3'') j='',I3,'' (''I3,'')'')') &
                            trim(var_names(j+plus)),pretty_small,ig+ioff,ig,jg+joff,jg
            LEVEL4 trim(msg)
            select case (i==iiC)
              case (.true.)
               write(msg,'(A,'' set on '',G12.6)')  trim(var_names(j+plus)),pretty_small
               forall (l=1:nlev) ccb3d(ig,jg,l,j)=pretty_small
              case (.false.)
               write(msg,'(A,'' reinitialized assuming a fixed ../C quotum of '',G12.6)')  &
                 trim(var_names(j+plus)), cx(nlev,j)/(1.0D-80+cx(nlev,k))
               forall (l=1:nlev) &
                  ccb3d(ig,jg,l,j)=max(_ZERO_,ccb3d(ig,jg,l,k))*cx(l,j)/(1.0D-80+cx(l,k))
            end select
            LEVEL4 trim(msg)
          endif
        endif
       end select
    endif
   enddo
#endif
  end subroutine calc_GrpNP

end module output_restart_bio

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2008 - BFM                                              !
!-----------------------------------------------------------------------

