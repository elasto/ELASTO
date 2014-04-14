!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: hdf5_wrapper_datalayout
!
!
! DESCRIPTION:
!> This module provide all procedure needed to perform parallel i/o at
!! hdf5 format from data stored with the "real_data_layout" type.
!
!! @details
!!        This module provides wrapper to use the basic hdf5 i/o to write
!!    and read data stored with the "real_datalayout" structure.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module hdf5_wrapper_datalayout

    use precision_tools
    use mpi
    use hdf5_io
    use datalayout

    implicit none


    !> To write files
    interface hdf5_write_datalayout
        module procedure hdf5_write_scalar_wrapper
    end interface hdf5_write_datalayout

    !> To read files
    interface hdf5_read_datalayout
        module procedure hdf5_read_scalar_wrapper
    end interface hdf5_read_datalayout
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
    !> Property list identifier for collective parallel I/O
    integer, private        :: plist_general_id

    contains


!> Write a 3D field inside a hdf5 file
!!    @param[in]    layout      = data to write
!!    @param[in]    iteration   = iteration number (for file name generation)
!!    @param[in]    cut         = to perform 2D cut, optional
!!    @param[in]    postname    = to add a postfix to the file name (for instance, for save, rediscretization, etc), optional
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_write_scalar_wrapper(layout, iteration, cut, postname) result(success)


    ! Input/output
    type(REAL_DATA_LAYOUT), intent(in)      :: layout
    integer, intent(in)                     :: iteration
    logical, optional, intent(in)           :: cut
    character(len=*), optional, intent(in)  :: postname
    logical                                 :: success

    ! Local variables
    character(len=str_medium)   :: filename
    integer, dimension(3,2)     :: mesh_info
    logical                     :: Dcut = .false.

#ifdef HDF5
    ! Init parameter
    write(filename,format_name) trim(layout%name), "_", iteration, ".h5"
    if (present(cut)) then
        if (cut) write(filename,format_name) trim(layout%name), "_2Dcut_", iteration, ".h5"
        Dcut = cut
    end if
    if (present(postname)) then
        write(filename,format_name) trim(layout%name), "_"//trim(postname)//"_", iteration, ".h5"
        !Dcut = .false.
    end if
    ! First point of the local array is:
    mesh_info(1,2) = layout%xmin
    mesh_info(2,2) = layout%ymin
    mesh_info(3,2) = layout%zmin
    ! Array glocal size is:
    mesh_info(1,1) = layout%nx
    mesh_info(2,1) = layout%ny
    mesh_info(3,1) = layout%nz

    success = hdf5_write(layout%values, filename, mesh_info, trim(layout%name), cut=Dcut)
#else
    success = .false.
#endif

end function hdf5_write_scalar_wrapper

!> Read file from an hdf5 file
!!    @param[out]   layout      = data to write
!!    @param[in]    filename    = name of the file to write
!!    @param[out]   success     = logical, true for success, false if an error occurs
!!    @param[in]    fieldname   = name of data (e.g. "scalar_1"), optional
subroutine hdf5_read_scalar_wrapper(layout, filename, success, fieldname, rank)

    ! Input/output
    type(REAL_DATA_LAYOUT), intent(inout)   :: layout
    character(len=*), intent(in)            :: filename
    logical, intent(out)                    :: success
    character(len=*), intent(in), optional  :: fieldname
    integer, optional, intent(in)           :: rank

    ! Local variables
    integer, dimension(3,2)             :: mesh_info


    success = .false.
#ifdef HDF5
    ! First point of the local array is:
    mesh_info(1,2) = layout%xmin
    mesh_info(2,2) = layout%ymin
    mesh_info(3,2) = layout%zmin
    ! Array glocal size is:
    mesh_info(1,1) = layout%nx
    mesh_info(2,1) = layout%ny
    mesh_info(3,1) = layout%nz

    ! Read data from file
    if(present(fieldname)) then
        call hdf5_read_scalar(layout%values, filename, mesh_info, success, &
            & dsetname=trim(fieldname), rank=rank)
    else
        call hdf5_read_scalar(layout%values, filename, mesh_info, success, &
            & dsetname = trim(layout%name), rank=rank)
    end if
#endif

end subroutine hdf5_read_scalar_wrapper


end module hdf5_wrapper_datalayout
!! @}
