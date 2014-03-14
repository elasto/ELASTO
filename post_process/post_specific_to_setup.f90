
module post_ref

    use precision_tools

    ! ==== Private variables ====
    ! -- mean values of 
    real(WP), private                                   :: mean_lambda, mean_eta, mean_etaB
    real(WP), private                                   :: mean_epsilon, mean_chi
    real(WP), private                                   :: mean_c_velo
    real(WP), dimension(:), allocatable, save, private  :: mean_c_OC, mean_c_B

contains

!------------------------------------------------------------------------------
!> Compute some quantities useful to plot theoritical decay of the spectrum of a
!! scalar with high schmidt number
!! @author Jean-Baptiste Lagaert
!!    @param[in]    U, V, W     = velocity field
!!    @param[in]    Scal        = scalar field
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @param[in]    nbcpus      = number of processes
!!    @return       success     = logical equals to true is everything is right.
!! @details
!!        This function compute chi, epsilon, eta_k and eta_B (more precisly the
!!    one associadted to the viscuous-convective range).
!------------------------------------------------------------------------------
function post_THI_spectrum_decay(U,V,W,Scal,spec_rank,nbcpus) result(success)

    use parser_tools
    use toolbox
    use stat_tools 
    use differential_tools 
    use datalayout
    use post_lib
    use post_velocity
    use post_scalar
    use transforms_tools

    ! Input/output
    type(real_data_layout), intent(in) :: U, V, W, Scal
!    type(real_data_layout), dimension(:), intent(in) :: Scal
    integer, intent(in) :: spec_rank, nbcpus
    logical             :: success

    ! Local variables
    type(complex_data_layout)           :: Uk, Vk, Wk
    type(real_data_layout),dimension(3) :: temp
    real(WP)                            :: E_kine, lambda, R_lambda,  eta_k, epsil, visco
    real(WP), dimension(3)              :: u_prime, di_u_i
    real(WP)                            :: eta_b, chi, schmidt_post, k_lambda
    type(wavenumbers)                   :: WN

    success = .false.

    call parser_read('Viscosity',visco)
    call parser_read('Schmidt number for scalar 1',schmidt_post)


    ! ===== Velocity =====
    ! -- Init --
    ! -- Velocity wave number --
    if(spec_rank==0) print*, 'start'
    success = initWN(WN,spec_rank,U%nx,U%ny,U%nz)
    success = computeWN(WN,spec_rank,U%Lx,U%Ly,U%Lz,U%nx,U%ny,U%nz)
    ! -- Get Uk, Vk and Wk --

    ! -- Spectral fields --
    !initialise COMPLEX_DATA_LAYOUT
    if(spec_rank==0) print*, 'Spec velo'
    if (.not. initDataLayout("U_Spectral", & 
            & Uk,(U%nx/2)+1,U%ny,U%nz, &
            & U%Lx,U%Ly,U%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for Uk : not enought memory!'
        return
    endif
    if (.not. initDataLayout("V_Spectral", & 
            & Vk,(V%nx/2)+1,V%ny,V%nz, &
            & V%Lx,V%Ly,V%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for Vk : not enought memory!'
        return
    endif
    if (.not. initDataLayout("W_Spectral", & 
            & Wk,(W%nx/2)+1,W%ny,W%nz, &
            & W%Lx,W%Ly,W%Lz,nbcpus,spec_rank,alongZ)) then
        write(6,'(a,i0,a)')'[ERROR] initDataLayout for Wk : not enought memory!'
        return
    endif

    ! -- compute FFT --
    call ftran(U,Uk,success)
    if ( .not. success) then
       write(6,'(a,i0,a)')'[ERROR] POs-THI - FFT for Uk: failed!'
       return
    endif
    call ftran(V,Vk,success)
    if ( .not. success) then
       write(6,'(a,i0,a)')'[ERROR] POs-THI - FFT for Vk: failed!'
       return
    endif
    call ftran(W,Wk,success)
    if ( .not. success) then
       write(6,'(a,i0,a)')'[ERROR] POs-THI - FFT for Wk: failed!'
       return
    endif

    ! -- Statistics --
    call compute_flow_stat(spec_rank, nbcpus, success,  & ! Working information
                        & U,V,W, Uk, Vk,Wk, WN, visco,  & ! Input
                        & E_kine, lambda, R_lambda,     & ! basic post-process
                        & u_prime, di_u_i,              & ! about isotropy
                        & epsil=epsil, eta=eta_k)         ! advanced post-process)
    k_lambda = 2.0_WP*acos(-1.0_WP)/lambda
    k_lambda = 1.0_WP/lambda
    ! -- Spectrum --
    !success = velo_spectrum_rescaled(Uk, Vk, Wk, WN, epsil, eta_k, nbcpus, visco, k_lambda, spec_rank)
    success = velo_spectrum_rescaled(Uk, Vk, Wk, WN, epsil, eta_k, k_lambda, spec_rank)

    ! -- Free memory --
    success = deleteWN(WN)
    call deleteDataLayout(Uk)
    call deleteDataLayout(Vk)
    call deleteDataLayout(Wk)

    ! ===== Scalar =====
    ! -- Init --
    success = copyStructOnly(Scal,temp(1))
    success = copyStructOnly(Scal,temp(2))
    success = copyStructOnly(Scal,temp(3))
    success = initWN(WN,spec_rank,Scal%nx,Scal%ny,Scal%nz)
    success = computeWN(WN,spec_rank,Scal%Lx,Scal%Ly,Scal%Lz,Scal%nx,Scal%ny,Scal%nz)
    ! -- Compute chi --
    success = computeFieldGradient(WN,Scal,temp(1),temp(2),temp(3),nbcpus,spec_rank)
    temp(1)%values = temp(1)%values**2 + temp(2)%values**2 + temp(3)%values**2
    chi = (visco/(schmidt_post*3._WP))*computeFieldAvg( temp(1) , spec_rank , nbcpus ) 

    ! -- Compute eta_k and eta_b --
    eta_b = eta_k/sqrt(schmidt_post)
    eta_B = eta_k 

    ! -- Screen output
    if (spec_rank==0) then
        write(*,'(a,e14.6)') 'lambda = ', lambda
        write(*,'(a,e14.6)') 'k_lambda = ', k_lambda
        write(*,'(a,e14.6)') 'chi = ', chi
        write(*,'(a,e14.6)') 'epsilon = ', epsil
        write(*,'(a,e14.6)') 'eta_k = ', eta_k
        write(*,'(a,e14.6)') 'eta_b = ', eta_b
    end if

    ! -- Compute rescaled scalar spectrum --
    success = compute_scalar_spectrum_rescaled(Scal, WN, chi, epsil, eta_k, eta_b, nbcpus, visco, k_lambda, spec_rank)

    ! -- Free memory
    call deleteDatalayout(temp(1))
    call deleteDatalayout(temp(2))
    call deleteDatalayout(temp(3))
    success = deleteWN(WN)

    success = .true.

end function post_THI_spectrum_decay



end module post_ref
