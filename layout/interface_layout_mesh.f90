!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!> @addtogroup cart_structure
!! @{

!-----------------------------------------------------------------------------
!
! MODULE: cart_mesh
!
!
! DESCRIPTION:
!>  This module provide a mesh structure. It is used for output and as future
!! base to deal with different scalar field computed with eventually different
!! resolutions.
!
!> @details
!!  This module provide structure to save mesh context associated to a field.
!! This allow to easily work with different resolutions and to know how
!! mesh interact with the mpi topology.
!! It provide the different tools to initialise the type to some default
!! value or to auto-complete it.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module interface_layout_mesh

    use precision_tools

    implicit none

    public

    ! ===== Public procedures =====
    !> Init a cartesian_mesh variable
    public      :: mesh_init

    ! ===== Private procedures =====
    ! Init a cartesian_mesh variable from a "real_data_layout" field.
    private     ::mesh_init_from_real_data_layout

    ! ===== Interface =====
    interface mesh_init
        module procedure mesh_init_from_real_data_layout
    end interface mesh_init


contains

!> Create the "cartesian_mesh" variable associated to a variable of type
!"real_data_layout".
!>    @param[out]   mesh    = variable of type cartesian_mesh where the data about mesh are save
!>    @param[in]    field   = variable of type real_data_layout
subroutine mesh_init_from_real_data_layout(mesh, field)

    use cart_mesh_tools
    use realdatalayout

    ! Input/Output
    type(cartesian_mesh), intent(out)       :: mesh
    type(real_data_layout), intent(in)      :: field
    ! Other local variables
    !integer                                 :: direction    ! integer matching to a direction (X, Y or Z)

    ! Range
    ! Absolute extend
    mesh%absolute_extend(1,:) = (/field%xmin, field%xmax/)
    mesh%absolute_extend(2,:) = (/field%ymin, field%ymax/)
    mesh%absolute_extend(3,:) = (/field%zmin, field%zmax/)
    ! Relative extend
    mesh%relative_extend(:,1) = 1
    mesh%relative_extend(:,2) = mesh%absolute_extend(:,2) - mesh%absolute_extend(:,1) + 1

    ! Number of mesh
    mesh%N = (/field%nx, field%ny, field%nz /)
    mesh%N_proc = mesh%relative_extend(:,2)

    ! Space step
    mesh%dx = (/field%Lx, field%Ly, field%Lz /)
    mesh%dx = mesh%dx/dble(mesh%N)

end subroutine mesh_init_from_real_data_layout

end module interface_layout_mesh
!> @}
