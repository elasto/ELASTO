!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_lib2D
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

module post_lib2D

    use precision_tools
    use datalayout
    use avs
    use toolbox 
    use stat_tools 
    use transforms_tools
    implicit none
    private

    ! ==== Public procedures ==== 
    public          :: compute_all_spectrum2D
    public          :: computeFieldPDF_2D
    public          :: computeScalarsPDF_2D
    public          :: compute_2D_AVG_quantities

    ! ==== Private variable ====
    logical, private:: init_mixedness = .true.

contains

! ===== Spectrum computations =====

!------------------------------------------------------------------------------
!> Compute and save spectrum of all field (scalars and velocity).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite             = current time iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    n_iter          = number of max iteration
!!    @param[in]    imhd            = logical equals to true if mhd is present
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of all scalars fields and of 
!!    kinetic energy. It takes directly the field in Fourrier Space from 
!!    the data module. This allows to avoid useless and expansive FFT transform,
!!    but it could only be used during a simulation (and not be performed 
!!    after).
!------------------------------------------------------------------------------
function compute_all_spectrum2D(ite, z_pos, spec_rank, n_iter,imhd, time) result(success)

    use datalayout

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    real(WP),dimension(:),intent(in):: z_pos
    integer, intent(in)             :: spec_rank
    logical, intent(in)             :: imhd
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    logical                         :: time_or_ite      ! output could be done at each N 
                                                        ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! Compute spectrum
    if (.not.compute_all_scalars_spectrum2D(ite, z_pos, spec_rank,n_iter,time_or_ite)) then
        write(*,'(a)') '[WARNING] Failed to compute 2D spectrum of scalars fiels'
        return
    else if (.not.compute_velo_spectrum2D(ite, z_pos, spec_rank,n_iter,time_or_ite)) then
        write(*,'(a)') '[WARNING] Failed to compute 2D spectrum of kinetic energy'
        return
   else if (imhd) then
       if (.not.compute_mhd_spectrum2D(ite, z_pos, spec_rank,n_iter, time_or_ite)) then
          write(*,'(a)') '[WARNING] Failed to compute spectrum of magnetic energy'
          return
       end if
    end if

    ! Return sucess
    success = .true.

end function compute_all_spectrum2D


!------------------------------------------------------------------------------
!> Compute and save spectrum of all scalars field.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite             = current time iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of all scalars fields. It 
!!    takes directly the field in Fourrier Space from the data module. This 
!!    allows to avoid useless and expansive FFT transform, but it could 
!!    only be used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
function compute_all_scalars_spectrum2D(ite, z_pos, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use parallel_tools
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    real(WP),dimension(:),intent(in):: z_pos
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:,:), allocatable   :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, sca, ik ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    character(len=15)                       :: name_short       ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.
    type(COMPLEX_DATA_LAYOUT)               :: ScalKbis         ! spectral transform of the scalar along X and Y only
    integer                                 :: indz             ! indice along z associate to z_pos


    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    ! -- Scalar solved with pseudo-spectral solver --
    if (allocated(ScalArrayk)) then
        do sca=1, size(ScalArrayk)
        if (res_t(sca)<1) then
            kc = 2.0_WP*acos(-1.0_WP)*(ScalArrayk(sca)%nx-1)/Uk%Lx
            size_spec = size(ScalWN(sca)%kx) ! Only this because de-aliasing
            ! -- Allocate temporaly fields --
            if (.not.(copystructonly(ScalArrayk(sca), Scalkbis))) write(*,'(a,i0,a)')   &
                &   '[error] 2D spectrum of Scal on process ',spec_rank,                &
                &   ': not enought memory to store spectral value along X and Y only'
            ! --FFT-back along Z to have spectral value anlong X and Y only --
            call btran_alongZ(ScalArrayK(sca), ScalKbis, success)
            if (.not. success) then
                write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
                &   ': failed in btran_alongZ for scal spec'
                stop
            end if
            ! -- Compute spectra --
            ! *kbis is alongZ, therefore, it contains all Z values !
            allocate(spectrum(size(z_pos),size_spec))
            allocate(WN_norm(size_spec))
            WN_norm = -1
            spectrum = 0
            do j =  ScalKbis%ymin,ScalKbis%ymax
                do i =  ScalKbis%xmin,ScalKbis%xmax
                    ! Compute norm of wave number
                    kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + &
                        & ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j)
                    kk = sqrt(real(kk))
                    if ((kk>eps).and.(kk<=kc)) then
                        ! Compute indice
                        ik = nint(kk/dk) + 1
                        do k = 1, size(z_pos)
                            ! Position above the middle
                            indz = nint((z_pos(k) + 1._WP)*ScalKbis%Nz/2.)
                            spectrum(k,ik) = spectrum(k,ik) &
                                & + ScalKbis%values(i,j,indz)*conjg(ScalKbis%values(i,j,indz))
                            ! Position above the middle
                            indz = nint((1._WP -z_pos(k)) * ScalKbis%Nz/2.)
                            spectrum(k,ik) = spectrum(k,ik) &
                                & + ScalKbis%values(i,j,indz)*conjg(ScalKbis%values(i,j,indz))
                        end do
                    end if
                end do
            end do
            ! Free memory
            call deleteDataLayout(ScalKbis)

            ! Sum values on all processes
            spectrum = DoLocalArray2DSum(spec_rank, spectrum)/2._WP
            ! Write output
            if (spec_rank==0) then
                do ik = 1, size_spec
                    WN_norm(ik) = dk*real(ik-1)
                end do
                do k = 1, size(z_pos)
                    file_id = iopen()
                    name_short = trim(adjustl(ScalArrayk(sca)%name))
                    if (.not. time_or_ite) then
                        if (n_iter .lt. 1000) then
                          write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'spec_',name_short, &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        elseif (n_iter .lt. 1000000) then
                          write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'spec_',name_short, &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        else
                          write(file_name,'(a,a,a,f4.2,a,i0,a)') 'spec_',name_short,   &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        endif
                    else
                        if (max(n_iter, ite) .lt. 1000) then
                          write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'spec_',name_short,  &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        elseif (n_iter .lt. 1000000) then
                          write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'spec_',name_short,  &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        else
                          write(file_name,'(a,a,a,f4.2,a,i0,a)') 'spec_',name_short,    &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        endif
                    end if
                    file_name = trim(adjustl(file_name))
                    open(unit=file_id, file=file_name, form="FORMATTED", ERR=100)
                    do ik = 1, size_spec
                        if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(k,ik)
                    end do
                    close(file_id)
                    file_id = iclose(file_id)
                end do
            end if
            ! Free memory for next scalar field
            deallocate(spectrum)
            deallocate(WN_norm)
        end if
        end do
    end if

    ! -- Scalar solved with mixed particle/spectral solver --
    if (allocated(Scal_partArrayk)) then
        do sca=1, size(Scal_partArrayk)
            kc = 2.0_WP*acos(-1.0_WP)*(Scal_partArrayk(sca)%nx-1)/Uk%Lx
            size_spec = size(Scal_partWN(sca)%kx)
            ! -- Allocate temporaly fields --
            if (.not.(copystructonly(Scal_partArrayk(sca), Scalkbis))) write(*,'(a,i0,a)')&
                &   '[error] 2D spectrum of Scal on process ',spec_rank,     &
                &   ': not enought memory to store spectral value along X and Y only'
            ! --FFT-back along Z to have spectral value anlong X and Y only --
            call btran_alongZ(Scal_partArrayK(sca), ScalKbis, success)
            if (.not. success) then
                write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
                &   ': failed in btran_alongZ for scal part'
                stop
            end if
            ! -- Compute spectra --
            ! *kbis is alongZ, therefore, it contains all Z values !
            allocate(spectrum(size(z_pos),size_spec))
            allocate(WN_norm(size_spec))
            WN_norm = -1
            spectrum = 0
            do j =  ScalKbis%ymin,ScalKbis%ymax
                do i =  ScalKbis%xmin,ScalKbis%xmax
                    ! Compute norm of wave number
                    kk = Scal_partWN(sca)%kx(i)*Scal_partWN(sca)%kx(i) + &
                        & Scal_partWN(sca)%ky(j)*Scal_partWN(sca)%ky(j)
                    kk = sqrt(real(kk))
                    if ((kk>eps).and.(kk<=kc)) then
                        ! Compute indice
                        ik = nint(kk/dk) + 1
                        do k = 1, size(z_pos)
                            ! Position above the middle
                            indz = nint((z_pos(k) + 1._WP)*ScalKbis%Nz/2.)
                            spectrum(k,ik) = spectrum(k,ik) &
                                & + ScalKbis%values(i,j,indz)*conjg(ScalKbis%values(i,j,indz))
                            ! Position above the middle
                            indz = nint((-z_pos(k) + 1._WP)*ScalKbis%Nz/2.)
                            spectrum(k,ik) = spectrum(k,ik) &
                                & + ScalKbis%values(i,j,indz)*conjg(ScalKbis%values(i,j,indz))
                        end do
                    end if
                end do
            end do
            ! Free memory
            call deleteDataLayout(ScalKbis)
            ! Sum values on all processes
            spectrum = DoLocalArray2DSum(spec_rank, spectrum)/2._WP
            ! Write output
            if (spec_rank==0) then
                do ik = 1, size_spec
                    WN_norm(ik) = dk*real(ik-1)
                end do
                do k = 1, size(z_pos)
                    file_id = iopen()
                    name_short = trim(adjustl(Scal_partArrayk(sca)%name))
                    if (.not. time_or_ite) then
                        if (n_iter .lt. 1000) then
                          write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'spec_',name_short, &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        elseif (n_iter .lt. 1000000) then
                          write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'spec_',name_short, &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        else
                          write(file_name,'(a,a,a,f4.2,a,i0,a)') 'spec_',name_short,   &
                            & '_z=',z_pos(k),'_at_',ite,'.table'
                        endif
                    else
                        if (max(n_iter, ite) .lt. 1000) then
                          write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'spec_',name_short,  &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        elseif (n_iter .lt. 1000000) then
                          write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'spec_',name_short,  &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        else
                          write(file_name,'(a,a,a,f4.2,a,i0,a)') 'spec_',name_short,    &
                            & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                        endif
                    end if
                    file_name = trim(adjustl(file_name))
                    open(unit=file_id, file=file_name, form="FORMATTED", ERR=150)
                    do ik = 1, size_spec
                        if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(k,ik)
                    end do
                    close(file_id)
                    file_id = iclose(file_id)
                end do
            end if
            ! Free memory for next scalar field
            deallocate(spectrum)
            deallocate(WN_norm)
        end do
    end if

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
100 write(6,'(a,i0,a, i0)') 'Unable to open file "'//trim(file_name)//'"for scal(',sca,') - rank ', spec_rank
    return
150 write(6,'(a,i0,a,i0)') 'Unable to open file "'//trim(file_name)//'"for sca_part(',sca,') - rank ', spec_rank
    return


end function compute_all_scalars_spectrum2D


!------------------------------------------------------------------------------
!> Compute and save spectrum of kinetic energy.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite     = current iteration step
!!    @param[in]    z_pos= wanted position along Z relative to LZ/2.
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending 
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum along X and Y of kinetic energy 
!!    in function of Z. Due to consideration about symmetry, it computes two
!!    spectrum and gives the average. It takes directly the field in Fourrier
!!    Space from the data module. This allows to avoid useless and expansive FFT:
!!    only a fft-backward along Z have to be performed and no communication
!!    is needed. The drawback is that it could only be used during a simulation 
!!    (and not be performed after).
!------------------------------------------------------------------------------
function compute_velo_spectrum2D(ite, z_pos, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    real(WP),dimension(:),intent(in):: z_pos
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:,:), allocatable   :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.
    type(COMPLEX_DATA_LAYOUT)               :: Ukbis,Vkbis,Wkbis! spectral transform of U,V,W along X and Y only
    integer                                 :: indz             ! indice along z associate to z_pos

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above X and Y direction =====
    if(Uk%Lx/=Uk%Ly) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    ! -- Allocate temporaly fields --
    if (.not.(copystructonly(Uk,Ukbis) .and. &
        &   copystructonly(Vk,Vkbis) .and. &
        &   copystructonly(Wk,Wkbis))) write(6,'(a,i0,a)')                  &
        &   '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': not enought memory to store spectral valuealong X and Y only'
    ! --FFT-back along Z to have spectral value anlong X and Y only --
    call btran_alongZ(Uk, Ukbis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for U'
        stop
    end if
    call btran_alongZ(Vk, Vkbis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for V'
        stop
    end if
    call btran_alongZ(Wk, Wkbis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for W'
        stop
    end if
    ! -- Compute spectra --
    ! *kbis is alongZ, therefore, it contains all Z values !
    size_spec = size(VelWN%kx) ! Only this because de-aliasing
    allocate(spectrum(size(z_pos),size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! -- Compute spectrum --
    do j =  UK%ymin,UK%ymax
        do i =  UK%xmin,UK%xmax
            ! Compute norm of wave number
            kk = VelWN%kx(i)*VelWN%kx(i) + &
                & VelWN%ky(j)*VelWN%ky(j)
            kk = sqrt(real(kk))
            if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                ik = nint(kk/dk) + 1
                do k = 1, size(z_pos)
                    ! And then update spectrum
                    ! Position above the middle
                    indz = nint((z_pos(k) + 1._WP)* Ukbis%Nz/2._WP)
                    spectrum(k,ik) = spectrum(k,ik) + Uk%values(i,j,indz)*conjg(Uk%values(i,j,indz)) &
                            & + Vk%values(i,j,indz)*conjg(Vk%values(i,j,indz)) &
                            & + Wk%values(i,j,indz)*conjg(Wk%values(i,j,indz))
                    ! Position above the middle
                    indz = nint((-z_pos(k) + 1._WP)*Ukbis%Nz/2._WP)
                    spectrum(k,ik) = spectrum(k,ik) + Uk%values(i,j,indz)*conjg(Uk%values(i,j,indz)) &
                            & + Vk%values(i,j,indz)*conjg(Vk%values(i,j,indz)) &
                            & + Wk%values(i,j,indz)*conjg(Wk%values(i,j,indz))
                    !WN_norm(ik) = dk*real(ik-1)
                end do
            end if
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoLocalArray2DSum(spec_rank, spectrum)/2._WP

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        do k = 1, size(z_pos)
            file_id = iopen()
            if (.not. time_or_ite) then
                if (n_iter .lt. 1000) then
                   write(file_name,'(a,f4.2,a,i3.3,a)') 'spec_vel_z=',z_pos(k),'_at_',ite,'.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,f4.2,a,i6.6,a)') 'spec_vel_z=',z_pos(k),'_at_',ite,'.table'
                else
                   write(file_name,'(a,f4.2,a,i0,a)') 'spec_vel_z=',z_pos(k),'_at_',ite,'.table'
                endif
            else
                if (max(n_iter, ite) .lt. 1000) then
                   write(file_name,'(a,f4.2,a,i3.3,a)') 'spec_vel_z=',z_pos(k),'_T_',ite,'deltaT.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,f4.2,a,i6.6,a)') 'spec_vel_z=',z_pos(k),'_T_',ite,'deltaT.table'
                else
                   write(file_name,'(a,f4.2,a,i0,a)') 'spec_vel_z=',z_pos(k),'_T_',ite,'deltaT.table'
                endif
            end if
            file_name = trim(adjustl(file_name))
            open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
            write(file_id,'(a)') '# Kinetic spectrum'
            write(file_id,'(a)') '# wave_number and spectrum '
            do ik = 1, size_spec
                if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(k,ik)
            end do
            close(file_id)
            file_id = iclose(file_id)
        end do
    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)
    call deleteDataLayout(Ukbis)
    call deleteDataLayout(Vkbis)
    call deleteDataLayout(Wkbis)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
200 write(6,'(a)') 'Unable to open file - velo "'//trim(file_name)//'"'
    return

end function compute_velo_spectrum2D


! !------------------------------------------------------------------------------
! !> Compute and save spectrum of magnetic energy.
! !! @author Jean-Baptiste Lagaert and Guillaume Balarac
! !!    @param[in]    ite     = current iteration step
! !!    @param[in]    n_iter  = max iteration step
! !!    @param[in]    time    = optional logical do adjust name if
! !!                              spectrum computation are done depending 
! !!                              time and not iterations.
! !!    @return       success = logical equals to true is everything is right.
! !! @details
! !!        This function compute the spectrum of magnetic energy. It takes 
! !!    directly the field in Fourrier Space from the data module. This allows 
! !!    to avoid useless and expansive FFT transform, but it could only be 
! !!    used during a simulation (and not be performed after).
! !------------------------------------------------------------------------------
function compute_mhd_spectrum2D(ite, z_pos, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    real(WP),dimension(:),intent(in):: z_pos
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:,:), allocatable   :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.
    type(COMPLEX_DATA_LAYOUT)               :: Bkxbis,Bkybis,Bkzbis! spectral transform of U,V,W along X and Y only
    integer                                 :: indz             ! indice along z associate to z_pos

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above X and Y direction =====
    if(Bk(1)%Lx/=Bk(1)%Ly) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Bk(1)%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(Bk(1)%nx-1)/Bk(1)%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    ! -- Allocate temporaly fields --
    if (.not.(copystructonly(Bk(1),Bkxbis) .and. &
        &   copystructonly(Bk(2),Bkybis) .and. &
        &   copystructonly(Bk(3),Bkzbis))) write(6,'(a,i0,a)')                  &
        &   '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': not enought memory to store spectral valuealong X and Y only'
    ! --FFT-back along Z to have spectral value anlong X and Y only --
    call btran_alongZ(Bk(1), Bkxbis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for U'
        stop
    end if
    call btran_alongZ(Bk(2), Bkybis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for V'
        stop
    end if
    call btran_alongZ(Bk(3), Bkzbis, success)
    if (.not. success) then
        write(*,'(a,i0,a)') '[error] compute2D spectrum on process ',spec_rank,      &
        &   ': failed in btran_alongZ for W'
        stop
    end if
    ! -- Compute spectra --
    ! *kbis is alongZ, therefore, it contains all Z values !
    size_spec = size(BfieldWN%kx) ! Only this because de-aliasing
    allocate(spectrum(size(z_pos),size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! -- Compute spectrum --
    do j =  Bk(1)%ymin,Bk(1)%ymax
        do i =  Bk(1)%xmin,Bk(1)%xmax
            ! Compute norm of wave number
            kk = BfieldWN%kx(i)*BfieldWN%kx(i) + &
                & BfieldWN%ky(j)*BfieldWN%ky(j)
            kk = sqrt(real(kk))
            if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                ik = nint(kk/dk) + 1
                do k = 1, size(z_pos)
                    ! And then update spectrum
                    ! Position above the middle
                    indz = nint((z_pos(k) + 1._WP)* Bkxbis%Nz/2._WP)
                    spectrum(k,ik) = spectrum(k,ik) + Bk(1)%values(i,j,indz)*conjg(Bk(1)%values(i,j,indz)) &
                            & + Bk(2)%values(i,j,indz)*conjg(Bk(2)%values(i,j,indz)) &
                            & + Bk(3)%values(i,j,indz)*conjg(Bk(3)%values(i,j,indz))
                    ! Position above the middle
                    indz = nint((-z_pos(k) + 1._WP)*Bkxbis%Nz/2._WP)
                    spectrum(k,ik) = spectrum(k,ik) + Bk(1)%values(i,j,indz)*conjg(Bk(1)%values(i,j,indz)) &
                            & + Bk(2)%values(i,j,indz)*conjg(Bk(2)%values(i,j,indz)) &
                            & + Bk(3)%values(i,j,indz)*conjg(Bk(3)%values(i,j,indz))
                    !WN_norm(ik) = dk*real(ik-1)
                end do
            end if
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoLocalArray2DSum(spec_rank, spectrum)/2._WP

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        do k = 1, size(z_pos)
            file_id = iopen()
            if (.not. time_or_ite) then
                if (n_iter .lt. 1000) then
                   write(file_name,'(a,f4.2,a,i3.3,a)') 'spec_mhd_z=',z_pos(k),'_at_',ite,'.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,f4.2,a,i6.6,a)') 'spec_mhd_z=',z_pos(k),'_at_',ite,'.table'
                else
                   write(file_name,'(a,f4.2,a,i0,a)') 'spec_mhd_z=',z_pos(k),'_at_',ite,'.table'
                endif
            else
                if (max(n_iter, ite) .lt. 1000) then
                   write(file_name,'(a,f4.2,a,i3.3,a)') 'spec_mhd_z=',z_pos(k),'_T_',ite,'deltaT.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,f4.2,a,i6.6,a)') 'spec_mhd_z=',z_pos(k),'_T_',ite,'deltaT.table'
                else
                   write(file_name,'(a,f4.2,a,i0,a)') 'spec_mhd_z=',z_pos(k),'_T_',ite,'deltaT.table'
                endif
            end if
            file_name = trim(adjustl(file_name))
            open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
            write(file_id,'(a)') '# magnetic spectrum'
            write(file_id,'(a)') '# wave_number and spectrum '
            do ik = 1, size_spec
                if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(k,ik)
            end do
            close(file_id)
            file_id = iclose(file_id)
        end do
    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)
    call deleteDataLayout(Bkxbis)
    call deleteDataLayout(Bkybis)
    call deleteDataLayout(Bkzbis)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
200 write(6,'(a)') 'Unable to open file - mhd "'//trim(file_name)//'"'
    return

end function compute_mhd_spectrum2D

! function compute_mhd_spectrum(ite, spec_rank,n_iter, time) result(success)
! 
!     use datalayout
!     use data 
!     use fileio
! 
!     ! Input/Output
!     integer, intent(in)             :: ite
!     integer, intent(in)             :: n_iter
!     integer, intent(in)             :: spec_rank
!     logical, intent(in), optional   :: time
!     logical                         :: success
!     ! Others variables
!     real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
!     real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
!     integer                                 :: size_spec        ! spectrum size
!     integer                                 :: i, j ,k, ik      ! indices
!     real(WP)                                :: kk               ! norm of wave number
!     real(WP)                                :: dk               ! space step in Fourrier space
!     real(WP)                                :: kc, eps          ! to filter spectrum
!     integer                                 :: file_id          ! id for file io
!     character(len=50)                       :: file_name        ! name of output file
!     logical                                 :: time_or_ite      ! output could be done at each N 
!                                                                 ! iterations or at each N*deltaT time.
! 
!     ! ===== Init =====
!     success = .false.
!     if(present(time)) then
!         time_or_ite = time
!     else
!         time_or_ite = .false.
!     end if
! 
!     ! ===== Check if spectral space step is the same above all direction =====
!     if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
!         write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
!         return
!     else
!         dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
!     end if
!     kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
!     eps = dk/2.0
! 
!     ! ===== Compute and write spectrum =====
!     !size_spec = size(VelWN%kx)*size(VelWN%ky)*size(VelWN%kz)
!     size_spec = size(BfieldWN%kx) ! Only this because de-aliasing
!     allocate(spectrum(size_spec))
!     allocate(WN_norm(size_spec))
!     WN_norm = -1
!     spectrum = 0
!     ! Compute
!     do k = Bk(1)%zmin,Bk(1)%zmax
!         do j =  Bk(1)%ymin,Bk(1)%ymax
!             do i =  Bk(1)%xmin,Bk(1)%xmax
!                 ! Compute norm of wave number
!                 kk = BfieldWN%kx(i)*BfieldWN%kx(i) + &
!                     & BfieldWN%ky(j)*BfieldWN%ky(j) + &
!                     & BfieldWN%kz(k)*BfieldWN%kz(k)
!                 kk = sqrt(real(kk))
!                 if ((kk>eps).and.(kk<=kc)) then
!                 ! Compute indice
!                     ik = nint(kk/dk) + 1
!                     ! And then update spectrum
!                     spectrum(ik) = spectrum(ik) + Bk(1)%values(i,j,k)*conjg(Bk(1)%values(i,j,k)) &
!                             & + Bk(2)%values(i,j,k)*conjg(Bk(2)%values(i,j,k)) &
!                             & + Bk(3)%values(i,j,k)*conjg(Bk(3)%values(i,j,k))
!                     !WN_norm(ik) = dk*real(ik-1)
!                 end if
!             end do
!         end do
!     end do
! 
!     ! ===== Sum values on all processes =====
!     spectrum = DoGlobalVectorSum(spec_rank, spectrum)
! 
!     ! ===== Write output =====
!     if (spec_rank==0) then
!         do ik = 1, size_spec
!             WN_norm(ik) = dk*real(ik-1)
!         end do
!         file_id = iopen()
!         if (.not. time_or_ite) then
!             if (n_iter .lt. 1000) then
!                write(file_name,'(a,i3.3,a)') 'spec_mhd_at_',ite,'.table'
!             elseif (n_iter .lt. 1000000) then
!                write(file_name,'(a,i6.6,a)') 'spec_mhd_at_',ite,'.table'
!             else
!                write(file_name,'(a,i0,a)') 'spec_mhd_at_',ite,'.table'
!             endif
!         else
!             if (max(n_iter, ite) .lt. 1000) then
!                write(file_name,'(a,i3.3,a)') 'spec_mhd_at_',ite,'deltaT.table'
!             elseif (max(n_iter, ite).lt. 1000000) then
!                write(file_name,'(a,i6.6,a)') 'spec_mhd_at_',ite,'deltaT.table'
!             else
!                write(file_name,'(a,i0,a)') 'spec_mhd_at_',ite,'deltaT.table'
!             endif
!         end if
!         file_name = trim(adjustl(file_name))
!         open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
!         write(file_id,'(a)') '# Magnetic spectrum'
!         write(file_id,'(a)') '# wave_number and spectrum '
!         do ik = 1, size_spec
!             if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
!         end do
!         close(file_id)
!         file_id = iclose(file_id)
!     end if
!     ! Free memory for next scalar field
!     deallocate(spectrum)
!     deallocate(WN_norm)
! 
!     ! ===== Return sucess =====
!     success = .true.
!     return
! 
!     ! ===== IO error =====
! 200 write(6,'(a)') 'Unable to open file - mhd "'//trim(file_name)//'"'
!     return
! 
! end function compute_mhd_spectrum


! Todo but not needed immediatly.
!!------------------------------------------------------------------------------
!!> Compute and save spectrum of a given field (scalars, not kinetic energy
!!! spectrum).
!!! @author Jean-Baptiste Lagaert
!!!    @param[in]    field           = field from wich spectrum has to been computed.
!!!    @param[in]    ite             = time iteration number
!!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!!    @param[in]    nbcpus          = number of processes
!!!    @return       success         = logical equals to true is everything is right.
!!! @details
!!!        This function compute the spectrum of a given field. It perform 
!!!    a FFT in order to obtain it in the Fourrier space (wich is required
!!!    to computing the spectrum). It could be perform "out a simulation"
!!!    (ie from a field write in an output).
!!------------------------------------------------------------------------------
!function compute_spectrum(field, ite, nbcpus, spec_rank) result(success)
!
!    use datalayout
!    use fileio
!    use parallel_tools
!    use wavenumber_tools 
!    use transforms_tools
!
!    ! Input/Output
!    type(REAL_DATA_LAYOUT), intent(in)      :: field
!    real(WP), intent(in)                    :: ite
!    integer, intent(in)                     :: spec_rank
!    integer, intent(in)                     :: nbcpus
!    logical                                 :: success
!    ! Others variables
!    type(COMPLEX_DATA_LAYOUT)               :: spec_field       ! to store spectral values of "field"
!    type(WaveNumbers)                       :: WN               ! associated wave numbers
!    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
!    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
!    integer                                 :: size_spec        ! spectrum size
!    integer                                 :: i, j ,k, ik      ! indices
!    integer                                 :: nk               ! indice for maximal wave number
!    real(WP)                                :: kk               ! norm of wave number
!    real(WP)                                :: dk               ! space step in Fourrier space
!    real(WP)                                :: kc, eps          ! to filter spectrum
!    integer                                 :: file_id          ! id for file io
!    character(len=50)                       :: file_name        ! name of output file
!    integer                                 :: ires             ! error code for allocation
!
!    ! Init
!    success = .false.
!
!    ! ===== Check if spectral space step is the same above all direction =====
!    if((field%Lx/=field%Ly).or.(field%Lx/=field%Lz)) then
!        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
!        return
!    else
!        dk=2.0_WP*acos(-1.0_WP)/field%Lx
!    end if
!    kc = 2.0_WP*acos(-1.0_WP)*(field%nx-1)/field%Lx
!    eps = dk/2.0
!
!
!    ! ===== Go in spectral space =====
!    ! Allocate storage in fourrier space
!    if (.not. initDataLayout(trim(field%name)//"_spectral", spec_field,(field%nx/2)+1, &
!            & field%ny,field%nz,field%Lx,field%Ly,field%Lz,nbcpus,spec_rank,alongZ)) then
!        write(6,'(a,i0,a,a)')'[ERROR] initSolver on process ',spec_rank,&
!            & ': not enought memory to transform field in fourier space in order ', &
!            & 'to compute spectrum.'
!        return
!    endif
!    ! Perform FFT
!    call ftran(field,spec_field,success)
!    if (.not. success) then
!        write(6,'(a)')'[ERROR] spectrum computation: cannot perform FFT'
!        return
!    end if
!    success = .false.
!
!    ! ===== Compute Wave Number =====
!    ! -- Allocate memory storage --
!    allocate(WN%kx(1:(field%nx/2)+1),WN%ky(1:field%ny), WN%kz(1:field%nz), stat=ires)
!    if (ires/=0) then
!     write(6,'(a,i0,a)')'[ERROR] on process ',spec_rank,&
!     & ': not enought memory for waves numbers.'
!     return
!    end if
!    ! -- Along X for Real <-> complex FFT --
!    do i = 1,(field%nx/2)+1
!       WN%kx(i)=real(i-1,WP)*dk
!    end do
!    ! along Y for complex <-> complex FFT
!    nk=(field%ny/2) + 1  
!    do i = 1,field%ny
!      if (i.gt.nk) then
!         WN%ky(i)=-real(field%ny+1-i,WP)*dk
!      else   
!         WN%ky(i)=real(i-1,WP)*dk
!      end if
!    end do
!    ! along Z for complex <-> complex FFT
!    nk=(field%nz/2) + 1  
!    do i = 1,field%nz
!      if (i.gt.nk) then
!         WN%kz(i)=-real(field%nz+1-i,WP)*dk
!      else   
!         WN%kz(i)=real(i-1,WP)*dk
!      end if
!    end do
!    ! -- Update Max --
!    WN%xkmax = WN%kx(field%nx/2+1)
!    WN%ykmax = WN%ky(field%ny/2+1)
!    WN%xkmax = WN%kz(field%nz/2+1)
!
!    ! ===== Compute spectrum =====
!    size_spec = size(WN%kx)
!    allocate(spectrum(size_spec))
!    allocate(WN_norm(size_spec))
!    WN_norm = -1
!    spectrum = 0
!    ! Compute
!    do k = spec_field%zmin,spec_field%zmax
!        do j =  spec_field%ymin,spec_field%ymax
!            do i =  spec_field%xmin,spec_field%xmax
!                ! Compute norm of wave number
!                kk = WN%kx(i)*WN%kx(i) + &
!                    & WN%ky(j)*WN%ky(j) + &
!                    & WN%kz(k)*WN%kz(k)
!                kk = sqrt(kk)
!                ! Compute indice
!                ik = nint(kk/dk) + 1
!                ! And then update spectrum
!                spectrum(ik) = spectrum(ik) &
!                    & + spec_field%values(i,j,k)*conjg(spec_field%values(i,j,k))
!            end do
!        end do
!    end do
!
!    ! ===== Sum values on all processes =====
!    spectrum = DoGlobalVectorSum(spec_rank, spectrum)
!
!    ! ===== Write output =====
!    if (spec_rank==0) then
!        do ik = 1, size_spec
!            WN_norm(ik) = dk*real(ik-1)
!        end do
!        file_id = iopen()
!        write(file_name,'(a,a,i0,a)') 'spec_',trim(field%name),ite,'.table'
!        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
!        do ik = 1, size_spec
!            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
!        end do
!        close(file_id)
!        file_id = iclose(file_id)
!    end if
!
!    ! ===== Free memory =====
!    deallocate(spectrum)
!    deallocate(WN_norm)
!    if (associated(WN%kx)) deallocate(WN%kx)
!    if (associated(WN%ky)) deallocate(WN%ky)
!    if (associated(WN%kz)) deallocate(WN%kz)
!    WN%kx=>  null()
!    WN%ky=>  null()
!    WN%kz=>  null()
!    call deleteDataLayout(spec_field)
!
!    ! ===== Return sucess =====
!    success = .true.
!    return
!
!    ! ===== IO error =====
!300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
!    return
!
!end function compute_spectrum


!------------------------------------------------------------------------------
!> Compute and print PDF (along X and Y, in function of Z) of a scalars field.
!! @author Antoine Vollant, LEGI
!!    @param[in]    ite         = current time iteration
!!    @param[in]    n_iter      = number max of iteration
!!    @param[in]    field       = field to compute PDF
!!    @param[in]    dime        = number of samples in PDF
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    time        = optional logical do adjust name if
!!                                  spectrum computation are done depending 
!!                                  time and not iterations.
!!    @return       success     = logical equals to true is everything is right.
!! @details
!!        This function compute the probability density function of a scalar field. 
!------------------------------------------------------------------------------
function computeFieldPDF_2D(ite,z_pos,field,dime,spec_rank,n_iter,minDelta,maxDelta,time) result(success)

    use datalayout
    use parallel_tools
    use fileio

    ! Input/Output
    TYPE(REAL_DATA_LAYOUT), intent(in)  :: field
    real(WP), dimension(:), intent(in)  :: z_pos
    integer, intent(in)                 :: ite
    integer, intent(in)                 :: n_iter
    integer, intent(in)                 :: dime 
    integer, intent(in)                 :: spec_rank
    real(WP),intent(in),optional        :: minDelta,maxDelta
    logical, intent(in), optional       :: time
    logical                             :: success
    ! Others variables
    real(WP), dimension(:,:), allocatable   :: yPdf  ! to save y coordinates of PDF 
    real(WP), dimension(:), allocatable     :: xPdf  ! to save x coordinates of PDF
    integer                                 :: i, j ,k ,ll ! indices
    real(WP)                                :: cpt  ! count 
    real(WP)                                :: mini,maxi 
    real(WP)                                :: delta ! size of one PDF's sample
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite     ! output could be done at each N 
                                                        ! iterations or at each N*deltaT time.
    integer                                 :: indz             ! indice along z associate to z_pos

    ! ===== Init =====
    success = .false.
    allocate ( yPdf(size(z_pos),dime))
    allocate(xPdf(dime))
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    if (present(minDelta) .and. present(maxDelta) ) then
        mini = minDelta
        maxi = maxDelta
    else
        mini = computeFieldMin( field , spec_rank )
        maxi = computeFieldMax( field , spec_rank )
    endif

    yPdf=0.0_WP
    cpt = 0.0_WP
    delta=(maxi-mini)/(dime-1.)
    do k = 1, size(z_pos)
        indz = nint((1._WP + z_pos(k))*field%Nz/2._WP)
        if ((indz>=field%zmin).and.(indz<=field%zmax)) then
            do j = field%ymin,field%ymax
                do i = field%xmin,field%xmax
                ll = computeIndice(mini , maxi , field%values(i,j,indz) , dime)
                yPdf(k,ll)=yPdf(k,ll)+1.       
                    cpt = cpt + 1.
                enddo
            enddo
        end if
        indz = nint((1._WP - z_pos(k))*field%Nz/2._WP)
        if ((indz>=field%zmin).and.(indz<=field%zmax)) then
            do j = field%ymin,field%ymax
                do i = field%xmin,field%xmax
                ll = computeIndice(mini , maxi , field%values(i,j,indz) , dime)
                yPdf(k,ll)=yPdf(k,ll)+1.       
                cpt = cpt + 1.
                enddo
            enddo
        end if
    enddo

        do ll = 1,dime
            xPdf(ll)=mini+(ll-0.5)*delta
        enddo
        yPdf= DoLocalArray2DSum(spec_rank,yPdf)/2._WP
        cpt= DoLocalSum(spec_rank,cpt)

    ! Write output
    if (spec_rank==0) then
        yPdf=yPdf/(real(cpt))
        do k = 1, size(z_pos)
            file_id = iopen()
            if (.not. time_or_ite) then
                if (n_iter .lt. 1000) then
                    write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_at_',ite,'.table'
                elseif (n_iter .lt. 1000000) then
                    write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_at_',ite,'.table'
                else
                    write(file_name,'(a,a,a,f4.2,a,i0,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_at_',ite,'.table'
                endif
            else
                if (max(n_iter, ite) .lt. 1000) then
                    write(file_name,'(a,a,a,f4.2,a,i3.3,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_at_',ite,'deltaT.table'
                elseif (max(n_iter,ite) .lt. 1000000) then
                    write(file_name,'(a,a,a,f4.2,a,i6.6,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                else
                    write(file_name,'(a,a,a,f4.2,a,i0,a)') 'pdf_',trim(adjustl(field%name)),&
                        & '_z=',z_pos(k),'_T_',ite,'deltaT.table'
                endif
            end if
            file_name = trim(adjustl(file_name))
            open(unit=file_id, file=file_name, form="FORMATTED")
            do i = 1,dime 
                write(file_id, '(e14.6,x,e14.6)') xPdf(i), yPdf(k,i)
            end do
            close(file_id)
            file_id = iclose(file_id)
        end do
    end if
    ! Free memory for next scalar field
    deallocate(yPdf)
    deallocate(xPdf)

    ! ===== Return sucess =====
    success = .true.
    return

end function computeFieldPDF_2D



!------------------------------------------------------------------------------
!> Compute and print PDF of a Scalars.
!! @author Antoine Vollant and Jean-Baptiste Lagaert (only small modifications), LEGI
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    scalarArray     = array of scalars
!!    @param[in]    nbscalar        = number of scalars in the array
!!    @param[in]    scalar_partArray= array of scalars solved with
!!                                      particle/spectral methods
!!    @param[in]    nbscalar_part   = number of scalars in the array - scalars
!!                                      solved with particle/spectral methods.
!!    @param[in]    dime            = number of samples in PDF
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the probability density function of all scalars fields. 
!------------------------------------------------------------------------------
function computeScalarsPDF_2D(ite,z_pos,scalArray,nbscal,scal_partArray, nbscal_part, &
    & dime,spec_rank,n_iter, time) result (success)

    use datalayout
    use data

    ! Input/Output
    type(REAL_DATA_LAYOUT), intent(in), dimension(:)    :: scalArray,scal_partArray
    real(WP), intent(in), dimension(:)    :: z_pos
    integer, intent(in) :: ite
    integer, intent(in) :: nbscal, nbscal_part
    integer, intent(in) :: dime 
    integer, intent(in) :: n_iter  
    integer, intent(in) :: spec_rank
    logical, intent(in), optional   :: time
    logical             :: success
    ! Others variables
    integer :: sca 
    logical :: time_or_ite

    ! Init
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! Compute pdf - for scalar solved with full spectral solver
    do sca=1,nbscal
        if (res_t(sca)<1) then
            if(.not.computeFieldPDF_2D(ite,z_pos,scalArray(sca),dime, &
                    & spec_rank,n_iter,time=time_or_ite)) then
                write(*,'(a,i3.3)') '[WARNING] Failed to compute PDF of scalars fiels ',sca
                return
            end if
        end if
    enddo
    ! Compute pdf - for scalar solved with full spectral solver
    do sca=1,nbscal_part
        if(.not.computeFieldPDF_2D(ite,z_pos,scal_partArray(sca),dime, &
                & spec_rank,n_iter, time=time_or_ite)) then
            write(*,'(a,i3.3)') '[WARNING] Failed to compute PDF of scalars_part fiels ',sca
            return
        end if
    enddo
    ! Return sucess
    success = .true.

end function computeScalarsPDF_2D


!> Performed some 2D post-traitement based on mean values. These post-traitement 
!! are actually well adapted for plane-jet and are gathered as some computed
!! values are utilized by several of them (and thus there computations are mutualized)
!! @autor Jean-Baptiste Lagaert, LEGI
!!    @param[in]    scalar          = array of scalars
!!    @param[in]    nbscalar        = number of scalars in the array
!!    @param[in]    scalar_part     = array of scalars solved with
!!                                      particle/spectral methods
!!    @param[in]    nbscalar_part   = number of scalars in the array - scalars
!!                                      solved with particle/spectral methods.
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along W
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!
!> @details
!!    The following post-processes are performed for each scalar field:
!!         - mean values are along X,Y in function of Z
!!         - RMS: <f'^2> = <ff> - <f>^2 (reuse mean value already computed)
!!         - Intensity of segregation: I = <f'^2>/[<f>(1-<f>)] or 0 if <f'2> vanishes
!!         - mixedness: Z(t) = (4/Lz) intregal belong Z of {<f>(z)[1-<f>(z)] dz}
!!    For the velocity, we compute the mean values and the RMS.
function compute_2D_AVG_quantities(Scal, nbscal, Scal_part, nbscal_part, &
            & U, V, W, ite, n_iter, spec_rank, nbcpus, time) result(success)

    use toolbox
    use fileio

    type(REAL_DATA_LAYOUT),intent(in),dimension(:)  :: Scal,Scal_part
    type(REAL_DATA_LAYOUT),intent(in)               :: U, V, W
    integer, intent(in)                             :: nbscal, nbscal_part
    integer, intent(in)                             :: ite
    integer, intent(in)                             :: n_iter
    integer, intent(in)                             :: spec_rank, nbcpus
    logical, intent(in), optional                   :: time
    logical                                         :: success

    ! Local variables
    type(REAL_DATA_LAYOUT)              :: tmp_data ! to save square field for RMS computation
    real(WP), dimension(:), allocatable :: avg,RMS,I! to store avg, RMS
    real(WP), dimension(:), allocatable :: Z        ! intensity of segregation and mixedness
    real(WP)                            :: dz       ! space step along Z
    integer                             :: file_id  ! to identify output files
    character(len=50)                   :: file_name! name of output file
    integer                             :: k, sca   ! for loops
    character(len=25)                   :: format_mix   ! format for writing mixedness of each scalar
    logical                             :: time_or_ite  ! output depend of time avancement or iterations ?

    ! ===== Initialisation =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if
    allocate(Z(nbscal+nbscal_part))
    file_name=" "
    
    ! ===== Initialisation for post-process on velocity =====
    allocate(avg(U%Nz))
    allocate(RMS(U%Nz))
    dz = U%Lz/U%Nz
    if (.not.(copystructonly(U, tmp_data))) write(*,'(a,i0,a)')     &
        &   '[error] 2D spectrum of Scal on process ',spec_rank,    &
        &   ': not enought memory to store U^2'

    ! ===== Velocity along X =====
    ! -- Mean values --
    if(.not.computeFieldAvg2D(U,avg,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute Avg2D for U'
    ! -- RMS --
    tmp_data%values = U%values*U%values
    if(.not.computeFieldAvg2D(tmp_data,RMS,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute RMS for U'
    RMS = RMS - avg*avg
    ! -- Write output --
    if (spec_rank==0) then
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_U_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_U_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'avg_U_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_U_T_',ite,'deltaT.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_U_T_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'avg_U_T_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=400)
        write(file_id,'(a)') '# Different post-process based on mean values'
        write(file_id,'(a)') '# position along Z / mean values / RMS'
        do k = 1, size(avg)
            write(file_id, '(3(e14.6,x))') k*dz, avg(k), RMS(k)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if


    ! ===== Velocity along Y =====
    ! -- Mean values --
    if(.not. computeFieldAvg2D(V,avg,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute Avg2D for V'
    ! -- RMS --
    tmp_data%values = V%values*V%values
    if(.not. computeFieldAvg2D(tmp_data,RMS,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute RMS for V'
    RMS = RMS - avg*avg
    ! -- Write output --
    if (spec_rank==0) then
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_V_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_V_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'avg_V_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_V_T_',ite,'deltaT.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_V_T_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'avg_V_T_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=300)
        write(file_id,'(a)') '# Different post-process based on mean values'
        write(file_id,'(a)') '# position along Z / mean values / RMS'
        do k = 1, size(avg)
            write(file_id, '(3(e14.6,x))') k*dz, avg(k), RMS(k)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Velocity along Z =====
    ! -- Mean values --
    if(.not. computeFieldAvg2D(W,avg,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute Avg2D for U'
    ! -- RMS --
    tmp_data%values = W%values*W%values
    if(.not. computeFieldAvg2D(tmp_data,RMS,spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute RMS for W'
    RMS = RMS - avg*avg
    ! -- Write output --
    if (spec_rank==0) then
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_W_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_W_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'avg_W_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'avg_W_T_',ite,'deltaT.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'avg_W_T_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'avg_W_T_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=300)
        write(file_id,'(a)') '# Different post-process based on mean values'
        write(file_id,'(a)') '# position along Z / mean values / RMS'
        do k = 1, size(avg)
            write(file_id, '(3(e14.6,x))') k*dz, avg(k), RMS(k)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Free memory =====
    deallocate(avg)
    deallocate(RMS)
    call deleteDataLayout(tmp_data)

    ! ===== Scalar solved with full pseudo-spectral solver =====
    do sca = 1, nbscal
        ! -- Some init --
        allocate(avg(Scal(sca)%Nz))
        allocate(RMS(Scal(sca)%Nz))
        allocate(I(Scal(sca)%Nz))
        dz = Scal(sca)%Lz/Scal(sca)%Nz
        if (.not.(copystructonly(Scal(sca), tmp_data))) write(*,'(a,i0,a)')     &
            &   '[error] 2D spectrum of Scal on process ',spec_rank,    &
            &   ': not enought memory to store ScalSpec^2'
        ! -- Mean values --
        if(.not. computeFieldAvg2D(Scal(sca),avg,spec_rank)) write(*,'(a,i0)')&
        & '[ERROR] Could not compute Avg2D for Scal_spec', sca
        ! -- RMS --
        tmp_data%values = Scal(sca)%values*Scal(sca)%values
        if(.not. computeFieldAvg2D(tmp_data,RMS,spec_rank)) write(*,'(a,i0)')&
        & '[ERROR] Could not compute RMS for Scal_spec', sca
        RMS = RMS - avg*avg
        ! -- Intensity of segregation --
        where(RMS==0)
            I = 0
        elsewhere
            I = RMS/(avg*(1._WP-avg))
        end where
        ! -- Mixedness --
        Z(sca) = 0
        do k = 2, size(avg)
            Z(sca) = Z(sca) + avg(k-1)*(1.0_WP-avg(k-1))/2._WP
            Z(sca) = Z(sca) + avg(k)*(1.0_WP-avg(k))
        end do
        Z(sca) = 4._WP*Z(sca)/(real(size(avg)-1,WP))
        ! -- Write output --
        if (spec_rank==0) then
            file_id = iopen()
            if (.not. time_or_ite) then
                if (n_iter .lt. 1000) then
                   write(file_name,'(a,a,a,i3.3,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_at_',ite,'.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,a,a,i6.6,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_at_',ite,'.table'
                else
                   write(file_name,'(a,a,a,i0,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_at_',ite,'.table'
                endif
            else
                if (max(n_iter, ite) .lt. 1000) then
                   write(file_name,'(a,a,a,i3.3,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_T_',ite,'deltaT.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,a,a,i6.6,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_T_',ite,'deltaT.table'
                else
                   write(file_name,'(a,a,a,i0,a)') 'avg_',trim(adjustl(Scal(sca)%name)),'_T_',ite,'deltaT.table'
                endif
            end if
            file_name = trim(adjustl(file_name))
            open(unit=file_id, file=file_name, form="FORMATTED", ERR=300)
            write(file_id,'(a)') '# Different post-process based on mean values'
            write(file_id,'(a)') '# position along Z / mean values / RMS'
            do k = 1, size(avg)
                write(file_id, '(4(e14.6,x))') k*dz, avg(k), RMS(k), I(k)
            end do
            close(file_id)
            file_id = iclose(file_id)
        end if
        ! -- Free memory --
        deallocate(avg)
        deallocate(RMS)
        deallocate(I)
        call deleteDataLayout(tmp_data)
    end do


    ! ===== Scalar solved with particle-spectral solver =====
    do sca = 1, nbscal_part
        ! -- Some init --
        allocate(avg(Scal_part(sca)%Nz))
        allocate(RMS(Scal_part(sca)%Nz))
        allocate(I(Scal_part(sca)%Nz))
        dz = Scal_part(sca)%Lz/Scal_part(sca)%Nz
        if (.not.(copystructonly(Scal_part(sca), tmp_data))) write(*,'(a,i0,a)')     &
            &   '[error] 2D spectrum of Scal on process ',spec_rank,    &
            &   ': not enought memory to store ScalPart^2'
        ! -- Mean values --
        if(.not. computeFieldAvg2D(Scal_part(sca),avg,spec_rank)) write(*,'(a,i0)')&
        & '[ERROR] Could not compute Avg2D for Scal_part', sca
        ! -- RMS --
        tmp_data%values = Scal_part(sca)%values*Scal_part(sca)%values
        if(.not. computeFieldAvg2D(tmp_data,RMS,spec_rank)) write(*,'(a,i0)')&
        & '[ERROR] Could not compute RMS for Scal_part', sca
        RMS = RMS - avg*avg
        ! -- Intensity of segregation --
        where(abs(RMS)<=1e-4)
            I = 0
        elsewhere
            I = RMS/(avg*(1._WP-avg))
        end where
        ! -- Mixedness --
        Z(sca+nbscal) = 0
        do k = 2, size(avg)
            Z(sca+nbscal) = Z(sca+nbscal) + avg(k-1)*(1.0_WP-avg(k-1))/2._WP
            Z(sca+nbscal) = Z(sca+nbscal) + avg(k)*(1.0_WP-avg(k))
        end do
        Z(sca+nbscal) = 4._WP*Z(sca+nbscal)/(real(size(avg)-1,WP))
        ! -- Write output --
        if (spec_rank==0) then
            file_id = iopen()
            if (.not. time_or_ite) then
                if (n_iter .lt. 1000) then
                   write(file_name,'(a,i0,a,i3.3,a)') 'avg_ScalPart_',sca,'_at_',ite,'.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,i0,a,i6.6,a)') 'avg_ScalPart_',sca,'_at_',ite,'.table'
                else
                   write(file_name,'(a,i0,a,i0,a)') 'avg_ScalPart_',sca,'_at_',ite,'.table'
                endif
            else
                if (max(n_iter, ite) .lt. 1000) then
                   write(file_name,'(a,i0,a,i3.3,a)') 'avg_ScalPart_',sca,'_T_',ite,'deltaT.table'
                elseif (n_iter .lt. 1000000) then
                   write(file_name,'(a,i0,a,i6.6,a)') 'avg_ScalPart_',sca,'_T_',ite,'deltaT.table'
                else
                   write(file_name,'(a,i0,a,i0,a)') 'avg_ScalPart_',sca,'_T_',ite,'deltaT.table'
                endif
            end if
            file_name = trim(adjustl(file_name))
            open(unit=file_id, file=file_name, form="FORMATTED", ERR=300)
            write(file_id,'(a)') '# Different post-process based on mean values'
            write(file_id,'(a)') '# position along Z / mean values / RMS'
            do k = 1, size(avg)
                write(file_id, '(4(e14.6,x))') k*dz, avg(k), RMS(k), I(k)
            end do
            close(file_id)
            file_id = iclose(file_id)
        end if
        ! -- Free memory --
        deallocate(avg)
        deallocate(RMS)
        deallocate(I)
        call deleteDataLayout(tmp_data)
    end do

    ! ===== Mixedness =====
    ! -- Write output --
    if (spec_rank==0) then
        file_id = iopen()
        if (.not. time_or_ite) then
            write(file_name,'(a,a,a)') 'mixedness_',trim(adjustl(Scal(1)%name)),'.table'
        else
            write(file_name,'(a,a,a)') 'mixedness_T_',trim(adjustl(Scal(1)%name)),'.table'
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", position="APPEND", ERR=350)
        if (init_mixedness) then
            write(file_id,'(a)') '# Mixedness for scalars (solved with particles method or not)'
            write(file_id,'(a,i0,a,i0,a)') '# ite  ', nbscal, ' scal_spec / ', nbscal_part, ' scal_part / '
            init_mixedness = .false.
        end if
        write(format_mix,'(a,i0,a)') '(i0,',size(Z),'(e14.6,x))'
        write(file_id, trim(format_mix)) ite, Z 
        close(file_id)
        file_id = iclose(file_id)
    end if
    ! -- Free memory --
    deallocate(Z)
    

    success = .true.
    return

    ! ===== IO error =====
400 write(6,'(a)') 'Unable to open file  bis - "'//trim(file_name)//'"'
300 write(6,'(a)') 'Unable to open file - "'//trim(file_name)//'"'
350 write(6,'(a)') 'Unable to open mixedness file - "'//trim(file_name)//'"'
    
end function compute_2D_AVG_quantities

end module post_lib2D
