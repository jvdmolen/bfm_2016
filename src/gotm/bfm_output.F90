!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_bfm --- BFM bio model \label{sec:bio_bfm}
!
! !INTERFACE:
     module bfm_output
!
! !DESCRIPTION:
!
!
! !USES:
!  default: all is private.
!     use bio_var, only: stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS, &
!                        stBenFluxS,stPelStateE,stPelDiagE,stPelFluxE,stBenStateE, &
!                        stBenDiagE,stBenFluxE,stSYSE
!#ifdef INCLUDE_DIAGNOS_PRF
!     use bio_var, only: stPRFDiagS,stPRFFluxS,stPRFDiagE,stPRFFluxE
!#endif

!
! !PUBLIC MEMBER FUNCTIONS
! store arrays for average computations (pelagic and benthic)
     logical , dimension(:) ,   allocatable           :: var_ave
     integer, dimension(:), allocatable    :: var_ids
     character(len=64), dimension(:), allocatable :: var_names
     character(len=64), dimension(:), allocatable :: var_units
     character(len=92), dimension(:), allocatable :: var_long

     REALTYPE, dimension(:,:),   allocatable, target   :: cc_ave
     REALTYPE, dimension(:,:),   allocatable, target   :: ccb_ave
     REALTYPE                                          :: ave_count
     logical                                           :: reset_count

#ifdef INCLUDE_DIAGNOS_PRF
     REALTYPE, dimension(:,:),   allocatable, target :: ccb_ave_prf
#endif

     logical      :: write_results
     ! Start and End markers for variable and diagnostics storage
     integer      :: stPelStateS
     integer      :: stPelDiagS
     integer      :: stPelFluxS
     integer      :: stBenStateS
     integer      :: stBenDiagS
     integer      :: stBenFluxS
     integer      :: stPelStateE
     integer      :: stPelDiagE
     integer      :: stPelFluxE
     integer      :: stBenStateE
     integer      :: stBenDiagE
     integer      :: stBenFluxE
     integer      :: stSYSE
#ifdef INCLUDE_DIAGNOS_PRF
     integer      :: stPRFDiagS
     integer      :: stPRFFluxS
     integer      :: stPRFDiagE
     integer      :: stPRFFluxE
#endif



!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

     contains

     subroutine setup_bio_output
     use global_mem,only:LOGUNIT
     use bio_var,only:numc,numc_diag,numc_flux,numbc,numbc_diag,numbc_flux
#ifdef INCLUDE_DIAGNOS_PRF
     use bio_var, only: numbc_prf
#endif
     implicit none
     integer                     ::k
     integer                     ::rc
     k=numc+numc_diag+numc_flux+numbc+numbc_diag+numbc_flux
#ifdef INCLUDE_DIAGNOS_PRF
     k=k+numbc_prf
#endif
     write(LOGUNIT,*) 'prepare_bio_output k=',k
     allocate(var_ave(1:k),stat=rc)
     if (rc /= 0) stop 'setup_bio_output(): Error allocating var_ave)'
     var_ave=.false.

     allocate(var_ids(1:k),stat=rc)
     if (rc /= 0) stop 'setup_bio_output(): Error allocating var_ids)'
     var_ids=0

     allocate(var_names(1:k),stat=rc)
     if (rc /= 0) stop 'setup_bio_output(): Error allocating var_names)'

     allocate(var_units(1:k),stat=rc)
     if (rc /= 0) stop 'setup_bio_output(): Error allocating var_units)'

     allocate(var_long(1:k),stat=rc)
     if (rc /= 0) stop 'setup_bio_output(): Error allocating var_long)'

     end subroutine setup_bio_output


     subroutine prepare_bio_output(mode)
     use global_mem,only:LOGUNIT
     use bio_var,only: nlev,dt, &
              cc,ccb,diag,diagb,c1dimz,bio_setup, &
              adv1d_courant,adv1d_number,numc
#ifdef INCLUDE_DIAGNOS_PRF
     use bio_var, only: numbc_prf,diagb_prf 
#endif

      implicit none
      integer,intent(IN)                    ::mode

      integer                     ::i
      integer                     ::j
      integer                     ::k
      integer                     ::rc
      logical                     ::llcalc
      character(len=90)           ::msg=""

      select case (mode)
        case(0)   ! initialization
!LEVEL1 'bfm output case 0',bio_setup,stPelStateS,stPelFluxE
!LEVEL1 'var_ave',var_ave
!stop
          i=count(var_ave(stPelStateS:stPelFluxE))
!LEVEL1 'i',i
          if ( (i> 0) .and. bio_setup/=2) then
!LEVEL1 'allocating'
            allocate(cc_ave(0:nlev,1:i),stat=rc)
!LEVEL1 'rc',rc
            if (rc /= 0) stop 'init_bio(): Error allocating cc_ave)'
             cc_ave=0
          endif
!LEVEL1 allocated(cc_ave)
          i=count(var_ave(stBenStateS:stBenFluxE))
          if ( ( i> 0) .and. bio_setup>1) then
            allocate(ccb_ave(0:1,1:i),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating cc_ave)'
            ccb_ave=0
          endif
#ifdef INCLUDE_DIAGNOS_PRF
          i=count(var_ave(stPRFDiagS:stPRFDiagE))
          if ( (i> 0) .and. bio_setup>1) then
            allocate(ccb_ave_prf(0:numbc_prf,1:i),stat=rc)
            if (rc/= 0) stop 'init_bio(): Error allocating ccb_ave_prf)'
            ccb_ave_prf=0
          endif
#endif

          ave_count=0.0
!LEVEL1 'bfm output case 0 end'
!JM: gaat fout, want var_ave is overal .false.
!stop
        case(1)  ! prepare for printing
          do j=1,numc
            if ( adv1d_courant(j)> 0.0 ) then
              write(msg,'(A,'':adv_center:MaxCourantN.='',G12.4, &
                &'' iter.='',F5.1)') trim(var_names(j)), &
                adv1d_courant(j)/ave_count,adv1d_number(j)/ave_count
              i=len_trim(msg);write(stderr,'(''        '',A)') msg(1:i)
              adv1d_courant(j)=0.0;adv1d_number(j)=0.0
             endif
           enddo
           if (bio_setup/=2.and.allocated(cc_ave)) cc_ave(:,:)=cc_ave(:,:)/ave_count
           if (bio_setup>1.and.allocated(ccb_ave)) ccb_ave=ccb_ave/ave_count
#ifdef INCLUDE_DIAGNOS_PRF
           if (bio_setup>1.and.allocated(ccb_ave_prf)) &
                                       ccb_ave_prf=ccb_ave_prf/ave_count
#endif
           ave_count=0.0
        case(10) ! Start of new time-step
           ave_count=ave_count+1.0
!LEVEL1 'case(10),stPelStateS,stPelStateE',stPelStateS,stPelStateE
!LEVEL1 'var_ave',var_ave
        case(11) ! add pel value
!LEVEL1 'case(11),stPelStateS,stPelStateE',stPelStateS,stPelStateE
!LEVEL1 'var_ave',var_ave
!stop
           k=0
           j=0
           if (stPelStateE==0 .or. bio_setup==2) return
           do i=stPelStateS,stPelStateE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                if ( ave_count< 1.5 ) then
                   cc_ave(:,k)=cc(:,j)
                else
                   cc_ave(:,k)=cc_ave(:,k)+cc(:,j)
                endif
             endif
           enddo
           j=0
           do i=stPelDiagS,stPelDiagE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                if ( ave_count< 1.5 ) then
!JM                   cc_ave(:,k)=diag(j,:)
                   cc_ave(:,k)=diag(:,j)
                else
!JM                   cc_ave(:,k)=cc_ave(:,k)+diag(j,:)
                   cc_ave(:,k)=cc_ave(:,k)+diag(:,j)
                endif
              endif
           enddo
           j=0
           do i=stPelFluxS,stPelFluxE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                call make_flux_output(1,j,0,nlev,dt,c1dimz,llcalc)
                if ( ave_count< 1.5 ) then
                   cc_ave(:,k)=c1dimz
                else
                   cc_ave(:,k)=cc_ave(:,k)+c1dimz
                endif
              endif
           enddo
        case(12) ! add ben value
           k=0
           j=0
           if (stBenStateE==0 .or. bio_setup==1) return
           do i=stBenStateS,stBenStateE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                if ( ave_count< 1.5 ) then
                   ccb_ave(0:1,k)=ccb(0:1,j)
                else
                   ccb_ave(0:1,k)=ccb_ave(0:1,k)+ccb(0:1,j)
                endif
              endif
           enddo
           j=0
           do i=stBenDiagS,stBenDiagE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                if ( ave_count< 1.5 ) then
!JM                   ccb_ave(0:1,k)=diagb(j,0:1)
                   ccb_ave(0:1,k)=diagb(0:1,j)
                else
!JM                   ccb_ave(0:1,k)=ccb_ave(0:1,k)+diagb(j,0:1)
                   ccb_ave(0:1,k)=ccb_ave(0:1,k)+diagb(0:1,j)
                endif
              endif
           enddo
           j=0
           do i=stBenFluxS,stBenFluxE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
!JM bug?                call make_flux_output(2,j,0,nlev,dt,c1dimz,llcalc)
                call make_flux_output(2,j,0,1,dt,c1dimz,llcalc)
                if ( ave_count< 1.5 ) then
                   ccb_ave(0:1,k)=c1dimz(0:1)
                else
                   ccb_ave(:,k)=ccb_ave(:,k)+c1dimz(0:1)
                endif
              endif
           enddo
#ifdef INCLUDE_DIAGNOS_PRF
           do i=stPRFDiagS,stPRFDiagE
             j=j+1
             if ( var_ave(i) ) then
                k=k+1
                if ( ave_count< 1.5 ) then
!JM                   ccb_ave_prf(k,:)=diagb_prf(:,j)
                   ccb_ave_prf(:,k)=diagb_prf(:,j)
                else
!JM                   ccb_ave_prf(k,:)=ccb_ave_prf(k,:)+diagb_prf(:,j)
                   ccb_ave_prf(:,k)=ccb_ave_prf(:,k)+diagb_prf(:,j)
                endif
              endif
           enddo
#endif
     end select
     return
     end subroutine prepare_bio_output

     subroutine select_vars_for_output(namlst)
     use global_mem,only:LOGUNIT
     use string_functions, ONLY: getseq_number,empty,set_nco_style

     implicit none
     integer,intent(IN)                ::namlst
     logical                           ::nco_style_vars
     !Maximum no variables which can be saved
     integer,parameter    :: NSAVE=201  
     character(len=64),dimension(NSAVE):: var_save
     character(len=64),dimension(NSAVE):: ave_save
     integer                           ::icontrol,i,ierr,k

     namelist /bfm_save_nml/ nco_style_vars,var_save, ave_save

     nco_style_vars=.false.
     var_save=""
     ave_save=""
     var_ave=.false.
     LEVEL2 "Read variables names which will be saved in netcdf-file"
     read(namlst,nml=bfm_save_nml,err=100,iostat=ierr)
     icontrol=1
100  if ( icontrol== 0 ) then
       STDERR 'iostat=',ierr
       stop 'init_var_bfm:I could not read bfm_save_nml'
     end if

     !---------------------------------------------
     ! Check variable to be saved and
     ! set the corresponding flag value in var_ids
     !---------------------------------------------
     STDERR "Variables selected for output"
     k=select_vars_from_list(0,var_save,NSAVE,var_names, &
                                          var_ids,var_ave,stPRFFluxE)
     STDERR "------------------------------------------------"
     STDERR "There are ",k, "sampled variables seleced for output"
     STDERR "------------------------------------------------"
     k=select_vars_from_list(1,ave_save,NSAVE,var_names, &
                                          var_ids,var_ave,stPRFFluxE)
     STDERR "------------------------------------------------"
     STDERR "There are ",k, "averaged variables seleced for output"
     STDERR "------------------------------------------------"
     LEVEL1 'nco_style=',nco_style_vars
     if (nco_style_vars) then
      LEVEL2 'All ( and ) in var. names are replace by _ (underscore)'
       do i=1,stPRFFluxE
        if (.not.( ((i.ge.stPelStateS).and.(i.le.stPelStateE)).or. &
             ((i.ge.stBenStateS).and.(i.le.stBenStateE)) ) ) then
          var_names(i)=set_nco_style(var_names(i))
        endif
       enddo
     endif
     end subroutine select_vars_for_output

     function select_vars_from_list( mode,save,nsave,names,ids,ave,n)
     use global_mem,only:LOGUNIT
     use string_functions, ONLY: getseq_number,empty

     implicit none
     integer                  ::select_vars_from_list
     integer,intent(in)       ::mode
     integer,intent(in)       ::nsave
     integer,intent(in)       ::n
     character(len=*),intent(IN),dimension(nsave) ::save
     character(len=*),intent(IN),dimension(n)     ::names
     integer,intent(INOUT),dimension(n)  ::ids
     logical,intent(INOUT),dimension(n)  ::ave
     
     integer      ::i,j,k

     k=0
     do i=1,nsave
       if (.NOT.empty(save(i))) then
         j=getseq_number(save(i),names,n,.TRUE.)
         if ( .NOT.empty(save(i)) .AND. j==0 ) then
            STDERR trim(save(i)), ':variable  does not exist!'
         elseif ( ids(j)==-2) then
           STDERR trim(save(i)), ': The variable ', &
           ' is used as pointer to an another variable'
           STDERR 'Therefor this variable can NOT be outputted'
           ids(j)=0
         else if ( ids(j)<0 ) then
             STDERR trim(save(i)),': Variable ', &
                 ' is already selected for output in var_save'
         else if ( j> 0 ) then
           ids(j)=-1
           if (mode==1)ave(j)=.true.
           STDERR trim(save(i))
           k=k+1
          end if
        end if
     end do
     select_vars_from_list=k
     if (.NOT.empty(save(NSAVE))) then
        STDERR 'Warning:',NSAVE-1,'or more variables are defined &
               & for output'
        STDERR 'Warning:Change parameter NSAVE in &
                 & init_var_bfm or limit number of var_save variables'
         STDERR 'Warning:Potential memory problem!!'
     endif
     end function select_vars_from_list

     subroutine unselect_var_from_output(text)
     use string_functions, ONLY: getseq_number,empty,set_nco_style
     implicit none
     character(len=*),intent(IN)    :: text
     integer                        ::j

     j=getseq_number(text,var_names,stPRFFluxE,.TRUE.)
     var_ids(j)=-2
     end subroutine unselect_var_from_output

     subroutine output_state_name(nr,mode,text)
       integer,intent(IN)             ::nr
       integer,intent(IN),optional    ::mode
       character(len=*)                ::text
       integer                        ::j
     
       if (present(mode)) j=mode
       text=trim(var_names(nr+j*stBenStateS-j))
     end subroutine output_state_name
     end module bfm_output


