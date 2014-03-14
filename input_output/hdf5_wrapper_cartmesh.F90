!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: hdf5_wrapper_cartmesh
!
!
! DESCRIPTION:
!> This module provide all procedure needed to perform parallel i/o at
!! hdf5 format from field following a regular and cartesian 3D MPI-topologie
!
!! @details
!!        This module provides wrapper for hdf5 input/output defined in the
!!    hdf5_io module.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module hdf5_wrapper_cartmesh

    use cart_topology
    use precision_tools

#ifdef HDF5


    implicit none


    !> To write files
    interface hdf5_write_cartmesh
        module procedure hdf5_write_scalar_wrapper
    end interface hdf5_write_cartmesh

    !> To read files
    interface hdf5_read_cartmesh
        module procedure hdf5_read_scalar_wrapper
    end interface hdf5_read_cartmesh
!    interface hdf5_init_field
!        module procedure hdf5_init_field_basic, hdf5_init_field_iodata, &
!            & hdf5__init_field_realdatalayout
!    end interface hdf5_c_init_field
!
!    ! ===== Public procedures =====
!    !> Initialisation of the hdf5 parallel context
!    public                  :: hdf5_init
!    !> Generique procedure to write field in distribued vtk files
!    public                  :: hdf5_write
!    ! Initialisation of the mesh context
!    public                  :: hdf5_init_field
!
!    ! ===== Private procedures =====
!    !> Specific procedure to init output context for a specific field.
!    private                  :: hdf5_init_field_realdatalayout
!    !> Specific procedure to write field in distribued vtk files
!    private                 :: hdf5_scalar_realdatalayout

    ! ===== Private variables =====

    contains


!> Write a 3D field inside a hdf5 file
!!    @param[in]    field       = data to write
!!    @param[in]    fieldname   = name of data (e.g. "scalar_1")
!!    @param[in]    iteration   = iteration number (for file name generation)
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_write_scalar_wrapper(field, fieldname, iteration) result(success)

    use hdf5_io

    ! Input/output
    real(WP), dimension(:,:,:), intent(in)  :: field
    character(len=*), intent(in)            :: fieldname
    integer, intent(in)                     :: iteration
    logical                                 :: success

    ! Local variables
    character(len=str_medium)   :: filename
    integer, dimension(3,2)     :: mesh_info

    ! Init parameter
    write(filename,format_name) trim(fieldname), "_", iteration,".h5"
    ! First point of the local array is:
    mesh_info(:,2) = 1 + coord*mesh_sc%N_proc
    ! Array glocal size is:
    mesh_info(:,1) = mesh_sc%N

    success = hdf5_write(field, filename, mesh_info, fieldname)

end function hdf5_write_scalar_wrapper


!> Read file from an hdf5 file
!!    @param[in]    field       = data to write
!!    @param[in]    fieldname   = name of data (e.g. "scalar_1")
!!    @param[in]    iteration   = iteration number (for file name generation)
!!    @param[out]   success     = logical, true for success, false if an error occurs
subroutine hdf5_read_scalar_wrapper(field, filename, fieldname, success)

    use hdf5_io

    ! Input/output
    real(WP), dimension(:,:,:), intent(inout)   :: field
    character(len=*), intent(in)                :: filename
    character(len=*), intent(in)                :: fieldname
    logical, intent(out)                        :: success

    ! Local variables
    integer, dimension(3,2)             :: mesh_info


    success = .false.
    ! First point of the local array is:
    mesh_info(:,2) = 1 + coord*mesh_sc%N_proc
    ! Array glocal size is:
    mesh_info(:,1) = mesh_sc%N

    ! Read data from file
    call hdf5_read_scalar(field, filename, mesh_info, success, fieldname)

end subroutine hdf5_read_scalar_wrapper


#endif

end module hdf5_wrapper_cartmesh
!! @}
