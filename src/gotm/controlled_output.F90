#include "DEBUG.h"
#include "INCLUDE.h"

     module controlled_messages

     integer                        ::n=0
     integer,pointer                ::m_status
     integer,pointer                ::l
     !Controlled output---------------------------------------------
     character(len=80), dimension(:), pointer :: texts
     integer, dimension(:), pointer :: ntexts
     integer, dimension(:), pointer :: nrs
     integer, dimension(:), pointer :: empty_places
     integer, dimension(:), pointer :: actual_messages
     character(len=80), dimension(:,:,:),allocatable,target :: mm_texts
     integer, dimension(:,:,:), allocatable,target :: mm_ntexts
     integer, dimension(:,:,:), allocatable,target :: mm_nrs
     integer, dimension(:,:,:), allocatable,target :: mm_empty_places
     integer, dimension(:,:,:), allocatable,target :: mm_actual_messages
     integer,dimension(:,:,:),allocatable,target    ::mm_int
     integer,pointer                        :: nplaces,last_n
     integer,pointer                        :: lmessage
!
     integer,parameter    :: NSAVE=30
     character(len=64),dimension(NSAVE):: messtag
     character(len=80)                 :: var_output,text_output
     integer,dimension(NSAVE)          :: limit
     integer                           :: n_messtag,steps

      contains

     subroutine controlled_output(ltimes,text,tail,pel_id,ben_id,layers)
     use global_mem,only:LOGUNIT
     use bfm_output,ONLY: var_names,var_long, stPelStates,stBenStates
     implicit none
     integer,intent(IN)                   ::ltimes
     character(len=*),intent(IN)          ::text
     character(len=*),intent(IN),optional ::tail
     integer,intent(IN),optional          ::pel_id,ben_id
     integer,optional                     ::layers
     integer                              ::local_id,local_sys,j,k,l,m
     character(len=20)                    :: tlayers
     local_id=-1;local_sys=-1
     var_output='' ;tlayers=' '
     if (present(pel_id)) then
        local_sys=pel_id-1;local_id=stPelStates
     elseif (present(ben_id)) then
        local_sys=ben_id-1;local_id=stBenStates
     endif
      select case (local_id<0.or.ltimes==0)
        case (.true.) ; !call controlled_eval(ltimes,text)
        case (.false.)
          if (local_id>0)then ;var_output=var_names(local_sys+local_id)
                         else ;var_output=var_long(local_sys-local_id);endif
          j=len_trim(text);k=len_trim(var_output);l=len_trim(tail)
          if (present(layers)) then
          if (layers>0) then
            select case (ltimes>layers)
              case(.false.);write(tlayers,'('' in 1-'',i2,'' layers'')') layers
              case(.true.); write(tlayers,'('' in >'',i2,'' layers'')') layers
            end select
          endif
          endif
          m=len_trim(tlayers)
          if (j+k+l+m+2>80) then
            write(LOGUNIT,*) text,tail
            write(LOGUNIT,*) 'control_output :input too long<80'
            stop
          endif
          write(text_output,'(A,''-'',A,'' '',A,A)') &
            text(1:j),var_output(1:k),tail(1:l),tlayers(1:m)
          !call controlled_eval(ltimes,text_output)
      end select
      end subroutine controlled_output

      subroutine controlled_eval (ltimes,text,control)
      use global_mem,only:LOGUNIT
      use string_functions,only:index_trim,empty,getseq_number
      use BFM_ERROR_MSG,only:set_warning_for_getm,BFM_ERROR

      implicit none
      integer,intent(IN)         :: ltimes
      character(len=*),intent(IN)::text
      integer,intent(IN),optional::control

      logical                    ::new_item,llprint,ll_suppress,llk
      integer                    ::k,rc,j,i,is,it

      llprint=(ltimes>0)
      if (.not.present(control)) then
      ! call in model
        n=n+1
        ! in the first round through the model only the number
        ! of calls to controlled_output is counted!
        if (m_status>=1) then
          if (llprint) then
            new_item=(nrs(n)==0)
            if (.not.new_item ) new_item=(ntexts(nrs(n))==0)
            if (.not.new_item ) then
              k=nrs(n)
              new_item=(index_trim(texts(k),text).ne.1)
              if (.not.new_item) ntexts(k)=abs(ntexts(k))+1
              lmessage=lmessage+1
              actual_messages(lmessage)=k
            endif
            if (new_item.and.l<size(ntexts)) then
              select case (nplaces>0)
                case(.true.) ; i=empty_places(nplaces);nplaces=nplaces-1
                case(.false.)
                  l=l+1;i=l
                  if (l==size(ntexts)) write(LOGUNIT,*) 'Warning: Some&
                    & warning may be suppressed during this timestep'
              end select
              nrs(n)=i
              ntexts(i)=1
!             write(LOGUNIT,*) 'text=',text(1:len_trim(text))
              texts(i)=text(1:len_trim(text))
              lmessage=lmessage+1
              actual_messages(lmessage)=i
            endif
          endif
        endif
      else
!       write(LOGUNIT,*) "controlled_output control=",control
        ! this is called after call of model
        if (m_status==0) then
          m_status=-1
          !  m_status==-1: now first do calls
          !  to controlled-messages_pointer 0+1
        elseif (control==0) then
!         write(LOGUNIT,*) 'control==0', llprint,lmessage
!         write(LOGUNIT,*) 'actual_messages', actual_messages
          if ( llprint) then
            if (n_messtag==0) then
              j=0
              do i=1,NSAVE
                limit(i)=limit(i)*steps/100
                if (.not.empty(messtag(i))) then
                  n_messtag=n_messtag+1
                  if (j+1<i) then
                    messtag(j+1)=messtag(i);limit(j+1)=limit(i)
                  endif
                  j=j+1
                endif
              enddo
            endif
            if ( m_status.ne.1) return
            do i=1,n
              j=nrs(last_n+i);llk=.false.
              if (j> 0) then
                do rc=1,lmessage
                  llk=(.not.llk).and.(j.eq.actual_messages(rc))
                enddo
                if (.not.llk)then
                 lmessage=lmessage+1;actual_messages(lmessage)=j
                endif
              endif
            enddo
          endif
          if ( m_status.ne.1) return
          do k=1,lmessage
            j=actual_messages(k)
            i=ntexts(j)
            ! all messages are printed of llprint==.true.
            ! othere wise only messages are printed which doe not
            ! appear any more.
            if ((i<0.and.(.not.llprint)).or.llprint) then
              i=abs(i)
              if (.not.empty(texts(j))) then
!---------------------------------------------------------------------
! check if output must be supressed or not
                ll_suppress=(n_messtag>0)
                if  (ll_suppress) then
                  it=index(texts(j),' ')-1
                  ll_suppress=(it>0)
                  if (ll_suppress) then
                    is=getseq_number(texts(j)(1:it),messtag,n_messtag,.true.)
                    ll_suppress=(is>0)
                    if (ll_suppress) ll_suppress=(limit(is)>=i)
                  endif
!---------------------------------------------------------------------
                endif
                if ( .not.ll_suppress) then
                  if (i==1) then;write(LOGUNIT,'(A)') trim(texts(j))
                  else;write(LOGUNIT,'(A,''('',I3,''x)'')')trim(texts(j)),i
                  endif
                  call set_warning_for_getm
                endif
              endif
              texts(j)='';ntexts(j)=0
              if (nplaces<size(empty_places)) then
                 nplaces=nplaces+1;empty_places(nplaces)=j
              endif
            else
              ntexts(j)=-abs(ntexts(j))
            endif
          enddo
          ! all messages are written
          lmessage=0
          last_n=n
        elseif ( control==1.and.m_status==1) then
          select case (llprint)
            case(.true.);l=0;ntexts(:)=0;texts(:)='';nrs(:)=0;nplaces=0
          end select
          last_n=0;n=0
        elseif( control==2) then
!         write(LOGUNIT,'(2A,4I6)') &
!               trim(text),'n,nrs(n:n+1)',n,nrs(n:n+1),size(nrs)
        endif
      endif
      end subroutine controlled_eval

      subroutine controlled_output_point(mode,i,j)
      use global_mem,only:LOGUNIT
      use BFM_ERROR_MSG,only:BFM_ERROR
      implicit none
      integer,intent(in)::mode,i,j
      integer,parameter ::II_STATUS=1,II_L=2,II_NPLACES=3, &
                            II_LAST_n=4,II_LMESSAGE=5
      integer           ::i0,j0,k,rc
      select case (mode)
      ! allocate room for all pointer scalar integers
      case (-1)
        allocate(mm_int(1:i,1:j,1:II_LMESSAGE),stat=rc)
      ! before first call through bio-model initialize 1 values
        mm_int=0
        m_status=>mm_int(1,1,ii_STATUS)
        m_status=0
        mm_int(:,:,ii_STATUS)=m_status
      ! create space for all messages in all grid points
      case (0)
        if ( m_status.ne.-1 ) then
          write(LOGUNIT,*) 'call controlled_output_point(0,.,.) &
             & after frist call thougth the model'
           call BFM_ERROR('Stop controlled_output_point','')
        endif
        ! steps setting for 3d-Model....
        if (i >0) then
          steps=i
          write(LOGUNIT,*) 'controlled_output steps in outdelt=',steps
        endif
         write(LOGUNIT,*) 'controlled_output messsage in step',n
        mm_int(:,:,II_STATUS)=2
        i0=size(mm_int,dim=1);j0=size(mm_int,dim=2)
        n=n*steps
        k=max(5,(n+1)/5)
        allocate(mm_texts(i0,j0,1:k),stat=rc)
        allocate(mm_ntexts(i0,j0,1:k),stat=rc)
        allocate(mm_empty_places(i0,j0,1:k),stat=rc)
        allocate(mm_actual_messages(i0,j0,1:2*n),stat=rc)
        allocate(mm_nrs(i0,j0,1:n),stat=rc)
        write(LOGUNIT,*) 'controlled_output_point all vars are allocated'
        ! Set all m_status's for all gridpoints on 2
        !set pointer to grid point
        n=0
      case(1)
        m_status=>mm_int(i,j,II_STATUS)
        if (m_status==2) then
          m_status=1
          mm_ntexts(:,:,:)=0;mm_nrs(:,:,:)=0;mm_int(:,:,II_NPLACES)=0
          mm_empty_places(:,:,:)=0;mm_actual_messages(:,:,:)=0
          mm_int(:,:,II_L)=0;mm_int(:,:,II_LMESSAGE)=0;mm_int(:,:,II_LAST_N)=0
          mm_int(:,:,II_STATUS)=m_status;mm_texts(:,:,:)=''
          write(LOGUNIT,*) &
               'controlled_output_point all vars initialized '
        endif
        l=>mm_int(i,j,II_L)
        nplaces=>mm_int(i,j,II_NPLACES)
        last_n=>mm_int(i,j,II_LAST_N)
        lmessage=>mm_int(i,j,II_LMESSAGE)
        texts=>mm_texts(i,j,:)
        ntexts=>mm_ntexts(i,j,:)
        empty_places=>mm_empty_places(i,j,:)
        actual_messages=>mm_actual_messages(i,j,:)
        nrs=>mm_nrs(i,j,:)
        if (m_status.ne.1) then
           write(LOGUNIT,*) "i,j,m_status",i,j,m_status
           call BFM_ERROR('Stop controlled_output_point','')
        endif
      end select
      end subroutine controlled_output_point

      function controlled_output_status()
      use global_mem,only:LOGUNIT
      implicit none
      integer               ::controlled_output_status
      controlled_output_status=m_status
      end function controlled_output_status

     subroutine manage_controlled_output(namlst,steps0)
     use global_mem,only:LOGUNIT,NML_READ,error_msg_prn
     use string_functions, ONLY: getseq_number,empty

     implicit none
     integer,intent(IN)                ::namlst
     integer,intent(iN)                ::steps0
     !Maximum no controls which can be saved
     integer                           ::error,rc
!    integer,dimension(:,:),allocatable::tmp

     namelist /controlled_output_nml/ messtag, limit

!    allocate(tmp(1:size(limit),1:2),stat=rc)
!    if ( error /=0)call error_msg_prn(ALLOC,"contolled_output","tmp")

     ! steps setting for 1d-Model....
     if (steps0 >0) then
       steps=steps0
       write(LOGUNIT,*)  &
        'start  manage_controlled_output steps in outdelt=',steps
     endif

     messtag='';limit=0
     read(namlst,nml=controlled_output_nml,iostat=error,err=100)
!    deallocate(tmp,stat=rc)
     if (error==0) return
        write(LOGUNIT,*) 'manage_controlled_output iostat=',error
     return
100  call error_msg_prn(NML_READ,"controlled_output.f90", &
                                  "Controlled_output parameters")
     end subroutine manage_controlled_output



     function logic_test_vector(ll,mess,text)
      use global_mem,only:LOGUNIT
      implicit none
      character(len=*),intent(INOUT),optional ::mess
      character(len=*),intent(IN),optional    ::text
      logical,intent(IN)      ::ll(:)
      integer                 ::logic_test_vector
      integer                 ::i,n

      n=0
      do i=1,size(ll)
         if (ll(i)) n=n+1
      enddo
      if ( n>0.and.present(mess).and.present(text)) mess=text
      logic_test_vector=n
      end function logic_test_vector


      end module

!---------------------------------------------------------------------



