!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_lib
!
!
! DESCRIPTION:
!>  This module provide a library of post-process routines. This library
!! provide all that you need to build your own post-process.
!
!> @details
!!  This module provide all post-process implement in scales. It has to
!! be used to build your own post-process traitement. Note that for post-
!! performed during computation, you have to used the module
!! post_in_simu routine (another exe will be add later to perform post-
!! process separately from simulation).
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module energy_transfers

    USE wavenumber_tools 
    use datalayout
    use parallel_tools
    use fileio
    use data
    use transforms_tools
    use kin_energy_transfers_lib
    implicit none
    private


    ! ==== Public procedures ====
    public          :: compute_energy_transfers

contains


subroutine post_energy_transfers_init(Uin,Vin,Win,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank
    logical :: res
   IF(.NOT.wavesNumber_init(VelWN,spec_rank,Uin%Lx,Uin%Ly,Uin%Lz,Uin%nx,Uin%ny,Uin%nz)) RETURN

    IF (.NOT. initDataLayout("U_Spectral", &
            & Uk,(Uin%nx/2)+1,Uin%ny,Uin%nz,Uin%Lx,Uin%Ly,Uin%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for U velocities in fourier space.'
        RETURN
    ENDIF
    IF (.NOT. initDataLayout("V_Spectral", &
            & Vk,(Vin%nx/2)+1,Vin%ny,Vin%nz,Vin%Lx,Vin%Ly,Vin%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for V velocities in fourier space.'
        RETURN
    ENDIF
    IF (.NOT. initDataLayout("W_Spectral", &
            & Wk,(Win%nx/2)+1,Win%ny,Win%nz,Win%Lx,Win%Ly,Win%Lz,nbcpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
            & ': not enought memory for W velocities in fourier space.'
        RETURN
    ENDIF
   CALL ftran(Uin,Uk,res)
   IF (.NOT.res) RETURN
   CALL ftran(Vin,Vk,res)
   IF (.NOT.res) RETURN
   CALL ftran(Win,Wk,res)
   IF (.NOT.res) RETURN

end subroutine post_energy_transfers_init

! ===== Compute all the terms =====

!------------------------------------------------------------------------------
!> Compute and save mean values of all quantities in Ennergie équations .
!! @author Mouloud KESSAR
!!    @param[in]    spec_rank     = mpi rank inside spectral communicator
!!    @param[in]    Uin,Vin,Win   = Velocity fields ins physical space
!!    @param[in]    Bin           = magnetic Field
!!    @param[in]    WaveVel       = velocity wavenumber
!!    @param[in]    WaveB         = Bfield wavenumber
!!    @param[in]    nbcpus        = numbers of cpus
!!    @param[in]    simtime       = time in the simulation
!!    @param[in]    tstep         = step time in the simulation
!!    
!! @details
!!        This function compute the mean values of all quantities in Ennergie équations
!!    stores it, and print this vlaues at each iteration
!------------------------------------------------------------------------------
subroutine compute_energy_transfers(Uin,Vin,Win,&
                         & nbcpus,spec_rank)

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    INTEGER :: nbcpus,spec_rank



  call post_energy_transfers_init(Uin,Vin,Win,&
                         & nbcpus,spec_rank)
   call compute_spec_Ekin_transfers(VelWN,spec_rank)

end subroutine compute_energy_transfers

!------------------------------------------------------------------------------
!> Compute and save spectrum of kinetic energy.
!! @author Mouloud Kessar
!!    @param[in]    ite     = current iteration step
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of kinetic energy. It takes
!!    directly the field in Fourrier Space from the solver. This allows
!!    to avoid useless and expansive FFT transform, but it could only be
!!    used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
subroutine compute_spec_Ekin_transfers(wave,spec_rank) 
    ! Input/Output
    integer, intent(in)             :: spec_rank
    logical                         :: success
    TYPE(WaveNumbers),INTENT(IN) :: Wave
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win

    ! Others variables
    real(WP), dimension(:,:), allocatable     :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, ik      ! indices
!     real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
!     logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
!     if(present(time)) then
!         time_or_ite = time
!     else
!         time_or_ite = .false.
!     end if

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
        write(*,'(a)') '[warning] transfers function not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    !size_spec = size(VelWN%kx)*size(VelWN%ky)*size(VelWN%kz)
    size_spec = size(wave%kx) ! Only this because de-aliasing
    allocate(spectrum(size_spec,size_spec))
!     allocate(WN_norm(size_spec))
!     WN_norm = -1
    call compute_all_TXXnm(Uk,Vk,Wk,wave,size_spec,spectrum,spec_rank)
    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        do i=1,size_spec
        file_id = iopen()
!         if (.not. time_or_ite) then
            if (size_spec .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'transf_vel_for_k_',i,'.table'
            elseif (size_spec .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'transf_vel_for_k_',i,'.table'
            else
               write(file_name,'(a,i0,a)') 'transf_vel_for_k_',i,'.table'
            endif
!         else
!             if (max(n_iter, ite) .lt. 1000) then
!                write(file_name,'(a,i3.3,a)') 'spec_vel_at_',ite,'deltaT.table'
!             elseif (n_iter .lt. 1000000) then
!                write(file_name,'(a,i6.6,a)') 'spec_vel_at_',ite,'deltaT.table'
!             else
!                write(file_name,'(a,i0,a)') 'spec_vel_at_',ite,'deltaT.table'
!             endif
!         end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
        write(file_id,'(a)') '# Kinetic transfers'
        write(file_id,'(a)') '# wave_number and transfer '
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(i,ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
       end do
    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
200 write(6,'(a)') 'Unable to open file - velo "'//trim(file_name)//'"'
    return

end subroutine compute_spec_Ekin_transfers


end module energy_transfers
!> @}
