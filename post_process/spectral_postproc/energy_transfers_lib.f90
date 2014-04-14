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

module kin_energy_transfers_lib
    USE wavenumber_tools 
    use datalayout
    use parallel_tools
    use fileio

    implicit none
    private


    ! ==== Public procedures ====
    public :: compute_all_TXXnm
contains


subroutine compute_all_TXXnm(VecXk,VecYk,VecZk,wave,size_Matrix,spectrum,spec_rank)

    TYPE(WaveNumbers),INTENT(IN) :: Wave
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    INTEGER  :: spec_rank
    Integer :: size_Matrix,i
    INTEGER :: i_n,j_n,k_n
    INTEGER :: i_m,j_m,k_m
    INTEGER :: i_q,j_q,k_q

    INTEGER :: ik_n,ik_m

    REAL(WP), dimension (:,:) :: spectrum
    REAL(WP), dimension (:), allocatable :: transfer_buff
    real(WP)                             :: kk               ! norm of wave number
    real(WP)                             :: dk               ! space step in Fourrier space
    real(WP)                             :: kc, eps          ! to filter spectrum
    real(WP)                             :: kx_q,ky_q,kz_q

    spectrum = 0.0_WP
    ! ===== Check if spectral space step is the same above all direction =====
    if((VecXk%Lx/=VecXk%Ly).or.(VecXk%Lx/=VecXk%Lz)) then
        write(*,'(a)') '[warning] transfers function not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/VecXk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(VecXk%nx-1)/VecXk%Lx
    eps = dk/2.0

    do k_n = VecXk%zmin,VecXk%zmax
        do j_n =  VecXk%ymin,VecXk%ymax
            do i_n =  VecXk%xmin,VecXk%xmax
                ! Compute norm of wave number
                kk = wave%kx(i_n)*wave%kx(i_n) + &
                    & wave%ky(j_n)*wave%ky(j_n) + &
                    & wave%kz(k_n)*wave%kz(k_n)
                kk = sqrt(real(kk))
                if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                    ik_n = nint(kk/dk) + 1

                do k_m = VecXk%zmin,VecXk%zmax
                   do j_m =  VecXk%ymin,VecXk%ymax
                      do i_m =  VecXk%xmin,VecXk%xmax
                      ! Compute norm of wave number
                      kk = wave%kx(i_m)*wave%kx(i_m) + &
                         & wave%ky(j_m)*wave%ky(j_m) + &
                         & wave%kz(k_m)*wave%kz(k_m)
                      kk = sqrt(real(kk))
                      if ((kk>eps).and.(kk<=kc)) then
                        ! Compute indice
                        ik_m = nint(kk/dk) + 1
                      !!!!!!!!calcul du q
                      kx_q = -(wave%kx(i_m)+wave%kx(i_n))
                      ky_q = -(wave%ky(j_m)+wave%ky(j_n))
                      kz_q = -(wave%kz(k_m)+wave%kz(k_n))

                      i_q = nint(kx_q*VecXk%Lx/(2.0_WP*acos(-1.0_WP)) + 1)

                      IF ( ky_q .lt. 0.0_WP ) THEN
                         j_q = nint( +ky_q*VecXk%Ly/(2.0_WP*acos(-1.0_WP)) +VecXk%ny+1)
                      ELSE   
                         j_q = nint( +ky_q*VecXk%Ly/(2.0_WP*acos(-1.0_WP)) +1)
                      ENDIF

                      IF ( kz_q .lt. 0.0_WP ) THEN
                         k_q=nint( +kz_q*VecXk%Lz/(2.0_WP*acos(-1.0_WP)) +VecXk%nz+1)
                      ELSE   
                         k_q=nint( +kz_q*VecXk%Lz/(2.0_WP*acos(-1.0_WP)) +1)
                      ENDIF

                        ! And then update spectrum
                        spectrum(ik_n,ik_m) = spectrum(ik_n,ik_m) &
                                        & -aimag( &
                                        &  (  VecXk%values(i_n,j_n,k_n)*(VecXk%values(i_m,j_m,k_m)) &
                                        &  +  VecYk%values(i_n,j_n,k_n)*(VecYk%values(i_m,j_m,k_m)) &
                                        &  +  VecZk%values(i_n,j_n,k_n)*(VecZk%values(i_m,j_m,k_m)) )&
                                        & *( wave%kx(i_n)*VecXk%values(i_q,j_q,k_q) &
                                        &  + wave%ky(j_n)*VecYk%values(i_q,j_q,k_q) &
                                        &  + wave%kz(k_n)*VecZk%values(i_q,j_q,k_q))) &
                                        & -aimag( &
                                        &  ( VecXk%values(i_n,j_n,k_n)*(VecXk%values(i_q,j_q,k_q)) &
                                        &  + VecYk%values(i_n,j_n,k_n)*(VecYk%values(i_q,j_q,k_q)) &
                                        &  + VecZk%values(i_n,j_n,k_n)*VecZk%values(i_q,j_q,k_q))  &
                                        & *( wave%kx(i_n)*VecXk%values(i_m,j_m,k_m) &
                                        &  + wave%ky(j_n)*VecYk%values(i_m,j_m,k_m) &
                                        &  + wave%kz(k_n)*VecZk%values(i_m,j_m,k_m) )) 
                        !WN_norm(ik) = dk*real(ik-1)
                        end if
                      end do
                    end do
                end do

                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    allocate(transfer_buff(size_Matrix))
    do i=1,size_Matrix    
     transfer_buff(:)=spectrum(i,:)
     transfer_buff = DoGlobalVectorSum(spec_rank,transfer_buff)
     spectrum(i,:)=transfer_buff(:)
    enddo
    deallocate(transfer_buff)
    
end subroutine

end module kin_energy_transfers_lib
!> @}
