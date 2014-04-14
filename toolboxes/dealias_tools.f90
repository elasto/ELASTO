!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: dealias_tools
!
!! DESCRIPTION:
!> This module implements dealiasing methods.
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!>
!> @details
!! This module implements 3 choices for dealiasing: none, isotropic and anisotropic.
!! The last one is usefull when the number point in the 3 dimensions are very heterogeneous
!! to avoid loosing to much informations.
!! \image html dealias.png "Influence of dealiasing method on scalar min/max behavior."
!> @author
!> Guillaume Balarac et Patrick BEGOU, LEGI
!------------------------------------------------------------------------------
module dealias_tools

    use precision_tools
    use datalayout
    use wavenumber_tools
    use data 

    implicit none

    public dealiasIsotrope
    public dealiasAnisotrope
    public get_dealiasMethod

contains

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> Dealiasing. This routine should be used when nx,ny and nz dimensions
!> are    homogenous.
!> In input file set:
!> Dealiasing method            : isotrope
!> and/or
!> Dealiasing method for scalar x            : isotrope
!> @param [in,out] nlvalues longitudinal velocity (fourier space)
!> @param [in] wavenumb the WaveNumbers associated to nlvalues variable
!------------------------------------------------------------------------------
subroutine dealiasIsotrope(nlvalues,wavenumb)

    type(COMPLEX_DATA_LAYOUT), intent(inout):: nlvalues
    type(WaveNumbers), intent(in)           :: wavenumb

    integer :: i,j,k
    real(WP) :: kk, isotropcut

    isotropcut=sqrt(3.0)*min(wavenumb%xkcut,wavenumb%ykcut,wavenumb%zkcut)
    do k = nlvalues%zmin,nlvalues%zmax
        do j =    nlvalues%ymin,nlvalues%ymax
            do i =    nlvalues%xmin,nlvalues%xmax
                 kk = sqrt(wavenumb%kx(i)*wavenumb%kx(i)    + &
                        & wavenumb%ky(j)*wavenumb%ky(j)     + &
                        & wavenumb%kz(k)*wavenumb%kz(k))
                 if (kk>isotropcut) nlvalues%values(i,j,k)    = (0.0_WP,0.0_WP)
            enddo
        enddo
    enddo

end subroutine dealiasIsotrope


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> Dealiasing. This routine should be used when nx,ny and nz dimensions
!> are very heterogenous and where an isotropic dealiasing would remove
!> too much information.
!> In input file set:
!> Dealiasing method            : anisotrope
!> and/or
!> Dealiasing method for scalar x            : anisotrope
!> @param [in,out] nlvalues longitudinal velocity (fourier space)
!> @param [in] wavenumb the WaveNumbers associated to nlvalues variable
!------------------------------------------------------------------------------
subroutine dealiasAnisotrope(nlvalues,wavenumb)

    type(COMPLEX_DATA_LAYOUT), intent(inout):: nlvalues
    type(WaveNumbers),intent(in)            :: wavenumb

    integer :: i,j,k

    do k = nlvalues%zmin,nlvalues%zmax
        do j =    nlvalues%ymin,nlvalues%ymax
            do i =    nlvalues%xmin,nlvalues%xmax
                 if(    abs(wavenumb%kx(i))>wavenumb%xkcut .or. &
                     &  abs(wavenumb%ky(j))>wavenumb%ykcut .or. &
                     &  abs(wavenumb%kz(k))>wavenumb%zkcut )    &
                     &        nlvalues%values(i,j,k)    = (0.0_WP,0.0_WP)
            enddo
        enddo
    enddo

end subroutine dealiasAnisotrope


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This function get the dealiasing method to use. If called with
!> no parameter, it searchs for the methode for dealiasing velocities.
!> If called with an integer prameter its searchs for the methode for
!> dealiasing this scalar.
!> Default: no dealiasing.
!> @param [in] spec_rank my MPI processe rank in the range [0:nbcpus-1].
!> @param [in] sca the scalar number.
!> @return 0 for no dealiasing, 1 for isotropic dealiasing, 2 for anisotropic dealiasing.
!------------------------------------------------------------------------------
FUNCTION get_dealiasMethod(spec_rank,sca) result(res)

!    use param

    INTEGER, INTENT(IN)             :: spec_rank
    INTEGER, INTENT(IN), OPTIONAL   :: sca
    INTEGER :: res

    CHARACTER(LEN=32):: dealias
    CHARACTER(LEN=32):: request

    IF(PRESENT(sca)) THEN
        WRITE(request,'(a,i0)') 'Dealiasing method for scalar ',sca
        CALL parser_read(request,dealias)
    ELSE
        CALL parser_read('Dealiasing method',dealias)
    ENDIF

    IF (dealias .EQ. 'none') THEN
        res=0
    ELSEIF(dealias .EQ. 'isotrope') THEN
        res=1
    ELSEIF(dealias .EQ. 'anisotrope') THEN
        res=2
    ELSE
        IF(PRESENT(sca)) THEN
            IF(spec_rank.EQ.0) WRITE(6,'(a,i0,a)')'[WARNING] Unknown dealiasing method ['//TRIM(dealias) &
             & //'] for scalar N. ',sca,' : defaulting to no dealiasing.'
        ELSE
            IF(spec_rank.EQ.0) WRITE(6,'(a)')'[WARNING] Unknown dealiasing method ['//TRIM(dealias) &
             & //'] for velocities : defaulting to no dealiasing.'
        ENDIF
        res=0
    ENDIF

END FUNCTION get_dealiasMethod

end module dealias_tools
!> @}
