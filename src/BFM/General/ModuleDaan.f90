module mem_Daan

    use global_mem, only: RLEN, LOGUNIT, NMLUNIT, NML_OPEN, NML_READ, &
                          error_msg_prn

    implicit none
    
    public
    
    real(RLEN) :: seconds_between_changes=0.0
    real(RLEN) :: delta_p_smd=0.D0, p_smd_min=0.D0, p_smd_max=1.D10
!    real(RLEN) :: deltaN1p=0.D0, deltaN3n=0.D0, deltaN4n=0.D0, deltaN5s=0.D0

    public InitDaan

    contains
    
    subroutine InitDaan()

        namelist /Daan_parameters/ &
            seconds_between_changes, &
            delta_p_smd, p_smd_min, p_smd_max  !&
!            deltaN1p, deltaN3n, deltaN4n, deltaN5s

        write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
        write(LOGUNIT,*) "#  Reading Daan parameters.."
        open(NMLUNIT,file='Daan.nml',status='old',action='read',err=100)
        read(NMLUNIT,nml=Daan_parameters,err=101)
        close(NMLUNIT)
        write(LOGUNIT,*) "#  Namelist is:"
        write(LOGUNIT,nml=Daan_parameters)
        
        RETURN
        
100 call error_msg_prn(NML_OPEN,"ModuleDaan.f90","Daan.nml")
101 call error_msg_prn(NML_READ,"ModuleDaan.f90","Daan_parameters")
        
    end subroutine InitDaan

end module mem_Daan
