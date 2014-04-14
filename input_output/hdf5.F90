!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: hdf5_io
!
!
! DESCRIPTION:
!> This module provide all procedure needed to perform parallel i/o at
!! hdf5 format.
!
!! @details
!!        This module developp tools to write or read data at hdf5 format into a MPI
!!    parallel context. Be aware that hdf5 files can be read directly by visit but
!!    requires in order to be read by paraview. Such a wrapper can be created, if
!!    needed, at a post-process stage and consists only of a xml description of
!!    the data as - for a 3D rectlinear grid - the type of grid (3D/2D structured),
!!    the number of points, their localisation (cell/node/etc centered), the grid
!!    spacing and origin.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module hdf5_io

    use precision_tools

#ifdef HDF5

    use hdf5

    implicit none


    !> To write files
    !interface hdf5_init_field
    !    module procedure hdf5_init_field_basic, hdf5_init_field_iodata, &
    !        & hdf5__init_field_realdatalayout
    !end interface hdf5_init_field

    interface hdf5_write
        module procedure hdf5_write_scalar, hdf5_write_vector
    end interface hdf5_write

    interface hdf5_read
        module procedure hdf5_read_scalar!, hdf5_read_vector
    end interface hdf5_read

    ! ===== Public procedures =====
    !> Initialisation of the hdf5 parallel context
    public                  :: hdf5_init
    !> Define output path
    public                  :: hdf5_set_path
    !> Generique procedure to write field in distribued vtk files
    public                  :: hdf5_write
    !> Procedure to close hdf5 context 
    public                  :: hdf5_close

    ! ===== Private procedures =====
    !> Specific procedure to write a simple field
    !private                 :: hdf5_write_scalar
    !> Specific procedure to write a 3D vector
    private                 :: hdf5_write_vector
    ! Initialisation of the mesh context
    private                 :: hdf5_init_subset_3D

    ! ===== Private variables =====
    !> Property list identifier for collective parallel I/O
    integer, private                :: plist_general_id
    !> write format for output name of *.hdf5
    character (len=100), protected  :: format_name
    !> path for output
    character (len=100), protected  :: path

    contains


!> Initialize the general context hdf5 context
!!    @param[in]    mpi_communicator    = mpi communicator
!!    @param[in]    nb_ite      = total number of iteration, in order to define i/o name format
!!    @return       success     = logical, true for success, false if an error occurs
!! @details
!!    This subroutine do the same job than the basic version but extract all information from the real_data_layout.
!!    Therefore, it is able to deal with subdomain of different size.
function hdf5_init(mpi_communicator, nb_ite) result(success)

    use mpi

    ! Input/output
    integer, intent(in)             :: mpi_communicator
    integer, intent(in), optional   :: nb_ite
    logical                         :: success

    ! Local variables
    integer                 :: size_ite             ! number of character used to write iteration numer
    integer                 :: error                ! error handler
    integer                 :: info                 ! for customisation of the mpi I/O (see, for instance, Babel's documentation on MPI I/O hints)

    ! Init
    success = .false.
    info = MPI_INFO_NULL

    ! Initialize HDF5 library and Fortran interfaces.
    call h5open_f(error)

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_general_id, error)
    call h5pset_fapl_mpio_f(plist_general_id, mpi_communicator, info, error)

    ! Define format of output name
    if (present(nb_ite)) then
        if (nb_ite<1000) then
            size_ite = 3
        elseif (nb_ite<1E6) then
            size_ite = 6
        else
            size_ite = 6
            do while(nb_ite>=10**size_ite)
                size_ite = size_ite + 1
            end do
        end if
    else
        size_ite = 6
    end if

    write(format_name,'(a,i0,a1,i0,a)') '(a,a,i',size_ite,'.',size_ite,',a)'
    path = trim("./")


    ! Check error
    if (error == 0) success = .true.

end function hdf5_init

!> Defined i/o path
!!    @param[in]    path        = required path
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_set_path(dir) result(success)

    character(len=*), intent(in)        :: dir
    logical                             :: success

    success = .false.

    path = trim(dir)
    success = .true.

end function hdf5_set_path


!> Write a 3D scalar field inside a hdf5 file
!!    @param[in]    field       = data to write
!!    @param[in]    filename    = name of the file to write
!!    @param[in]    mesh_info   = basic mesh info: size of local mesh and coordinate of the "first" local point inside the global grid
!!    @param[in]    dsetname    = name of data (e.g. "scalar_1"), optional
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_write_scalar(field, filename, mesh_info, dsetname, cut) result(success)

    ! Input/output
    real(WP), dimension(:,:,:), intent(in)  :: field
    character(len=*), intent(in)            :: filename
    integer, dimension(3,2), intent(in)     :: mesh_info
    character(len=*), optional, intent(in)  :: dsetname
    logical, optional, intent(in)           :: cut
    logical                                 :: success

    ! Local variables
    integer                         :: error        ! error flag
    logical                         :: cut_bis=.false.  ! to write 2D cut
    integer, dimension(3)           :: first_point  ! 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
    integer(hsize_t), dimension(3)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3)  :: global_size  ! size of the local field
    integer(HID_T)                  :: file_id      ! File identifier
    integer(HID_T)                  :: dset_id      ! Dataset identifier
    character(len=str_medium)      :: dsetname_bis ! Dataset identifier
    integer(HID_T)                  :: filespace    ! Dataspace identifier in file
    integer(HID_T)                  :: memspace     ! Dataspace identifier in memory
    integer(HID_T)                  :: plist_id     ! Property list identifier

    ! Init parameter
    success = .false.
    if (present(dsetname)) then
        dsetname_bis = trim(dsetname)
    else
        dsetname_bis = trim(filename)
    end if
    if(present(cut)) cut_bis = cut
    local_size  = shape(field)
    global_size = mesh_info(:,1)
    first_point = mesh_info(:,2)

    ! Create the file collectively.
    call h5fcreate_f(trim(path)//trim(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_general_id)

    ! Create the dataspace and the dataset
    if(cut_bis) then
        call hdf5_init_subset_2D(first_point, local_size, global_size, file_id, .true., dsetname_bis, &
            & dset_id, filespace, memspace, plist_id, success)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, global_size, error, &
            & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    else
        call hdf5_init_subset_3D(first_point, local_size, global_size, file_id, .true., dsetname_bis, &
            & dset_id, filespace, memspace, plist_id, success)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, global_size, error, &
            & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    end if

    ! write the dataset collectively.
    !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, global_size, error, &
    !    & file_space_id = filespace, mem_space_id = memspace)!, xfer_prp = plist_id)

    ! Free ressources
    success = hdf5_delete_subset_3D(filespace, memspace, dset_id, plist_id)
    call h5fclose_f(file_id, error)

end function hdf5_write_scalar


!> Write a 3D vector (3 component) inside a hdf5 file
!!    @param[in]    Vx          = first component of the vector
!!    @param[in]    Vy          = second component of the vector
!!    @param[in]    Vz          = third component of the vector
!!    @param[in]    filename    = name of the file to write
!!    @param[in]    dsetname    = name of each component to write (e.g. "U", "V", "W")
!!    @param[in]    mesh_info   = basic mesh info: size of local mesh and coordinate of the "first" local point inside the global grid
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_write_vector(Vx, Vy, Vz, filename, mesh_info, dsetname) result(success)

    ! Input/output
    real(WP), dimension(:,:,:), intent(in)  :: Vx, Vy, Vz
    character(len=*), intent(in)            :: filename
    character(len=*), dimension(3), intent(in)  :: dsetname
    integer, dimension(3,2), intent(in)     :: mesh_info
    logical                                 :: success

    ! Local variables
    integer                         :: error        ! error flag
    integer, dimension(3)           :: first_point  ! 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
    integer(hsize_t), dimension(3)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3)  :: global_size  ! size of the local field
    integer(HID_T)                  :: file_id      ! File identifier
    integer(HID_T)                  :: dset_id      ! Dataset identifier
    integer(HID_T)                  :: dset_bis_id  ! Dataset identifier
    integer(HID_T)                  :: filespace    ! Dataspace identifier in file
    integer(HID_T)                  :: memspace     ! Dataspace identifier in memory
    integer(HID_T)                  :: plist_id     ! Property list identifier

    ! Init parameter
    success = .false.
    local_size  = shape(Vx)
    global_size = mesh_info(:,1)
    first_point = mesh_info(:,2)

    ! Check size of each component (do they match?)
    if (.not.((size(Vx,1)==size(Vy,1)).and.(size(Vx,1)==size(Vz,1)))) return
    if (.not.((size(Vx,2)==size(Vy,2)).and.(size(Vx,2)==size(Vz,2)))) return
    if (.not.((size(Vx,3)==size(Vy,3)).and.(size(Vx,3)==size(Vz,3)))) return

    ! Create the file collectively.
    call h5fcreate_f(trim(path//filename), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_general_id)

    ! Create the dataspace and the dataset
    call hdf5_init_subset_3D(first_point, local_size, global_size, file_id, .true., &
        & dsetname(1), dset_id, filespace, memspace, plist_id, success)

    ! write the dataset collectively.
    ! Vx
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Vx, global_size, error, &
        & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, error)
    ! Vy
    call hdf5_add_dataset(local_size, global_size, file_id, .true., dsetname(2), dset_bis_id, plist_id, success)
    call h5dwrite_f(dset_bis_id, H5T_NATIVE_DOUBLE, Vy, global_size, error, &
        & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    call h5dclose_f(dset_bis_id, error)
    call h5pclose_f(plist_id, error)
    ! Vz
    call hdf5_add_dataset(local_size, global_size, file_id, .true., dsetname(3), dset_bis_id, plist_id, success)
    call h5dwrite_f(dset_bis_id, H5T_NATIVE_DOUBLE, Vz, global_size, error, &
        & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    call h5dclose_f(dset_bis_id, error)

    ! Free ressources
    success = hdf5_delete_subset_3D(filespace, memspace, dset_id, plist_id)
    call h5fclose_f(file_id, error)

end function hdf5_write_vector


!> Read file from an hdf5 file
!!    @param[in,out]    field       = data to write
!!    @param[in]        filename    = name of the file to write
!!    @param[in]        mesh_info   = basic mesh info: size of local mesh and coordinate of the "first" local point inside the global grid
!!    @param[out]       success     = logical, true for success, false if an error occurs
!!    @param[in]        dsetname    = name of data (e.g. "scalar_1"), optional
subroutine hdf5_read_scalar(field, filename, mesh_info, success, dsetname, rank)

    ! Input/output
    real(WP), dimension(:,:,:), intent(inout)  :: field
    character(len=*), intent(in)            :: filename
    integer, dimension(3,2), intent(in)     :: mesh_info
    logical, intent(out)                    :: success
    character(len=*), optional, intent(in)  :: dsetname
    integer, optional, intent(in)           :: rank

    ! Local variables
    integer                         :: error        ! error flag
    integer, dimension(3)           :: first_point  ! 3D absolute coordinate of the first point of the local field 
                                                    ! (ie the subset stored by the current mpi process)
    integer(hsize_t), dimension(3)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3)  :: global_size  ! size of the local field
    integer(HID_T)                  :: file_id      ! File identifier
    integer(HID_T)                  :: dset_id      ! Dataset identifier
    integer(HID_T)                  :: filespace    ! Dataspace identifier in file
    integer(HID_T)                  :: memspace     ! Dataspace identifier in memory
    integer(HID_T)                  :: plist_id     ! Property list identifier

    ! Init parameter
    local_size  = shape(field)
    global_size = mesh_info(:,1)    ! note that it is possible to get this (and the associate dataspace) directly from the file to read.
    first_point = mesh_info(:,2)

    ! Open the file collectively.
    call h5fopen_f(trim(path)//trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_general_id)
    if (present(rank).and.(error == -1)) then
        write(6,'(a,i0,a)')"ERROR on process ",rank,": cannot open file ["//trim(filename)//"]"
        return
    end if
    ! Create the dataspace and the dataset
    if (present(dsetname)) then
        call hdf5_init_subset_3D(first_point, local_size, global_size, file_id, .false., &
            & dsetname, dset_id, filespace, memspace, plist_id, success)
        if (present(rank).and.(.not. success)) then
            write(6,'(a,i0,a)')"ERROR on process ",rank,": cannot custom dataset ["//trim(dsetname)//"] in file ["//trim(filename)//"]"
            return
        end if
    else
        call hdf5_init_subset_3D(first_point, local_size, global_size, file_id, .false., &
            & filename, dset_id, filespace, memspace, plist_id, success)
        if (present(rank).and.(.not. success)) then
            write(6,'(a,i0,a)')"ERROR on process ",rank,": cannot default dataset ["//trim(filename)//"] in file ["//trim(filename)//"]"
            return
        end if
    end if

    ! Read data
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, local_size, error, &
        & mem_space_id = memspace,file_space_id = filespace, xfer_prp = plist_id)

    ! Free ressources
    success = hdf5_delete_subset_3D(filespace, memspace, dset_id, plist_id)
    call h5fclose_f(file_id, error)

end subroutine hdf5_read_scalar


!------------------------------------------------------------------------------
!> @author
!> Jean-Baptiste Lagaert
!
!> @details
!> This function get the size of the data field to read into an avs file.
!
!> @param[in]   spec_rank   = the rank number of the process in the range [0,ncpus - 1].
!> @param[in]   name        = the name of the output file.
!> @param[out]  mesh_size   = number of point along each direction
!> @param[out]  success     = .TRUE. if size is obtain
!------------------------------------------------------------------------------
subroutine hdf5_GetSizeFromFile(spec_rank,filename,dsetname,mesh_size, success)

    ! Input/Output
    integer, intent(in)                     :: spec_rank
    CHARACTER(len=*), intent(in)            :: filename
    character(len=*), intent(in)            :: dsetname 
    integer, dimension(3),intent(out)       :: mesh_size
    logical                                 :: success

    ! Local variables
    integer(HID_T)                  :: file_id      ! File identifier 
    integer(hid_t)                  :: filespace! Dataspace identifier 
    integer(HID_T)                  :: dset_id      ! Dataset identifier 
    integer(hsize_t), dimension(3)  :: dims     ! Array to store dimension sizes 
    integer(hsize_t), dimension(3)  :: maxdims  ! Array to store max dimension size
    integer                         :: hdferr   ! Error code = dataspace rank on success, -1 on failure

    ! Open the file collectively.
    call h5fopen_f(trim(path)//trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr, access_prp = plist_general_id)
    if (hdferr == -1) then
        write(6,'(a,i0,a)')"ERROR 'GetSizeFrom' on process ",spec_rank,": cannot open file ["//trim(filename)//"]"
        return
    end if

    ! Open dataset and get the dataspace
    call h5dopen_f(file_id, dsetname, dset_id, hdferr)
    if (hdferr == -1) then
        write(6,'(a,i0,a)')"ERROR 'GetSizeFrom' -  on process ",spec_rank,": cannot open data ["//trim(dsetname)//"] in file ["//trim(filename)//"]"
        return
    end if
    call h5dget_space_f(dset_id, filespace, hdferr)

    ! Obtain the dataspace dimension size
    call h5sget_simple_extent_dims_f(filespace, dims, maxdims, hdferr) 
    mesh_Size = dims

    ! Close all - free memory
    ! Close dataspaces.
    call h5sclose_f(filespace, hdferr)
    !if (hdferr /= 0) return
    ! Close the dataset.
    call h5dclose_f(dset_id, hdferr)
    !if (hdferr /= 0) return
    ! Close file
    call h5fclose_f(file_id, hdferr)

    !All is OK now!
    success=.true.

end subroutine hdf5_GetSizeFromFile


!> Initialise the parallel I/O context for a given geometry
!!    @param[in]    first_point = 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
!!    @param[in]    local_size  = size of the local field
!!    @param[in]    global_size = size of the total grid
!!    @param[in]    file_id     = File identifier
!!    @param[in]    new         = is the data set new or does it already exist? (true for output, false for input)
!!    @param[in]    dsetname    = dataset name (= name of the field to save in the dataset)
!!    @param[out]   dset_id     = Dataset identifier
!!    @param[out]   filespace   = Dataspace identifier in file
!!    @param[out]   memspace    = Dataspace identifier in memory
!!    @param[out]   plist_id    = Property list identifier
!!    @param[out]   success     = logical, true for success, false if an error occurs
subroutine hdf5_init_subset_3D(first_point, local_size, global_size, file_id, new, dsetname, dset_id, filespace, memspace, plist_id, success)

    ! Input/output
    integer, dimension(3), intent(in)           :: first_point  ! 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
    integer(hsize_t), dimension(3), intent(in)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3), intent(in)  :: global_size  ! size of the local field
    integer(HID_T), intent(in)                  :: file_id      ! File identifier
    logical, intent(in)                         :: new          ! is the data set new or does it already exist?
    character(len=*), intent(in)                :: dsetname     ! dataset name (= name of the field to save in the dataset)
    integer(HID_T), intent(out)                 :: dset_id      ! Dataset identifier
    integer(HID_T), intent(out)                 :: filespace    ! Dataspace identifier in file
    integer(HID_T), intent(out)                 :: memspace     ! Dataspace identifier in memory
    integer(HID_T), intent(out)                 :: plist_id     ! Property list identifier
    logical, intent(out)                        :: success

    ! Local variables
    integer                                 :: error        ! error flag
    integer                                 :: rank = 3     ! array rank
    integer(hsize_t),  dimension(3)         :: count        ! number of piece to write for each process
    integer(hsize_t),  dimension(3)         :: stride       ! space between each element
    integer(hsize_t),  dimension(3)         :: block        ! size of each piece
    integer(hsize_t),  dimension(3)         :: offset       ! position of the local first element inside the global 3D array

    success = .false.

    ! Create the data space for the  dataset.
    call h5screate_simple_f(rank, global_size, filespace, error)
    !if (error /= 0) return
    call h5screate_simple_f(rank, local_size, memspace, error)
    !if (error /= 0) return

    ! Create chunked dataset.
    if(new) then
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
        !if (error /= 0) return
        call h5pset_chunk_f(plist_id, rank, local_size, error)
        !if (error /= 0) return
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
            & dset_id, error, plist_id)
        !if (error /= 0) return
        call h5sclose_f(filespace, error)
        call h5pclose_f(plist_id, error)
        !if (error /= 0) return
    else
        !call h5dopen_f(file_id, dsetname, dset_id, error, plist_id)
        call h5dopen_f(file_id, dsetname, dset_id, error)
        !if (error /= 0) return
        call h5sclose_f(filespace, error)
    end if
    !if (error /= 0) return

    ! Each process defines dataset in memory and writes it to the hyperslab in
    ! the file
    stride   = 1
    count    = 1
    block    = local_size
    offset   = first_point - 1
    ! Select hyperslab in the file
    call h5dget_space_f(dset_id, filespace, error)
    !if (error /= 0) return
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, &
        & stride,block)
    !if (error /= 0) return

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    success = .true.

end subroutine hdf5_init_subset_3D


!> Initialise the parallel I/O context for a given geometry
!!    @param[in]    first_point = 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
!!    @param[in]    local_size  = size of the local field
!!    @param[in]    global_size = size of the total grid
!!    @param[in]    file_id     = File identifier
!!    @param[in]    new         = is the data set new or does it already exist? (true for output, false for input)
!!    @param[in]    dsetname    = dataset name (= name of the field to save in the dataset)
!!    @param[out]   dset_id     = Dataset identifier
!!    @param[out]   filespace   = Dataspace identifier in file
!!    @param[out]   memspace    = Dataspace identifier in memory
!!    @param[out]   plist_id    = Property list identifier
!!    @param[out]   success     = logical, true for success, false if an error occurs
subroutine hdf5_init_subset_2D(first_point, local_size, global_size, file_id, new, dsetname, dset_id, filespace, memspace, plist_id, success)

    ! Input/output
    integer, dimension(3), intent(in)           :: first_point  ! 3D absolute coordinate of the first point of the local field (ie the subset stored by the current mpi process)
    integer(hsize_t), dimension(3), intent(in)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3), intent(in)  :: global_size  ! size of the global field to write
    integer(HID_T), intent(in)                  :: file_id      ! File identifier
    logical, intent(in)                         :: new          ! is the data set new or does it already exist?
    character(len=*), intent(in)                :: dsetname     ! dataset name (= name of the field to save in the dataset)
    integer(HID_T), intent(out)                 :: dset_id      ! Dataset identifier
    integer(HID_T), intent(out)                 :: filespace    ! Dataspace identifier in file
    integer(HID_T), intent(out)                 :: memspace     ! Dataspace identifier in memory
    integer(HID_T), intent(out)                 :: plist_id     ! Property list identifier
    logical, intent(out)                        :: success

    ! Local variables
    integer                                 :: error        ! error flag
    integer                                 :: rank = 2     ! array rank
    integer(hsize_t),  dimension(3)         :: count        ! number of piece to write for each process
    integer(hsize_t),  dimension(3)         :: stride       ! space between each element
    integer(hsize_t),  dimension(3)         :: block        ! size of each piece
    integer(hsize_t),  dimension(3)         :: offset       ! position of the local first element inside the global 3D array
    integer(hsize_t),  dimension(2)         :: count_2D     ! number of piece to write for each process
    integer(hsize_t),  dimension(2)         :: stride_2D    ! space between each element
    integer(hsize_t),  dimension(2)         :: block_2D     ! size of each piece
    integer(hsize_t),  dimension(2)         :: offset_2D    ! position of the local first element inside the global 3D array

    success = .false.

    ! Create the data space for the  dataset.
    rank = 2
    call h5screate_simple_f(rank, global_size(2:3), filespace, error)
    !if (error /= 0) return
    rank = 3
    call h5screate_simple_f(rank, local_size, memspace, error)
    !if (error /= 0) return

    ! Create chunked dataset.
    if(new) then
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
        !if (error /= 0) return
        rank = 2
        call h5pset_chunk_f(plist_id, rank, local_size(2:3), error)
        !if (error /= 0) return
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
            & dset_id, error, plist_id)
        !if (error /= 0) return
        call h5sclose_f(filespace, error)
        call h5pclose_f(plist_id, error)
        !if (error /= 0) return
    else
        !call h5dopen_f(file_id, dsetname, dset_id, error, plist_id)
        call h5dopen_f(file_id, dsetname, dset_id, error)
        !if (error /= 0) return
        call h5sclose_f(filespace, error)
    end if
    !if (error /= 0) return

    ! Define wich part of the local data are inside the 2D cut to write
    ! Each process defines dataset in memory and writes it to the hyperslab in the file
    stride   = 1
    count    = 1
    block    = local_size
    offset   = 0
    offset(1)= block(1)/2
    if (offset(1)<0) offset(1) = 0
    block(1) = 1
    if (error /= 0) return
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, count, error, &
        & stride,block)
    !if (error /= 0) return

    ! Define wich part of the global data is dealed by the current process
    ! Each process defines dataset in memory and writes it to the hyperslab in the file
    stride_2D   = 1
    count_2D    = 1
    block_2D    = block(2:3)
    offset_2D   = first_point(2:3) - 1
    ! Select hyperslab in the file
    call h5dget_space_f(dset_id, filespace, error)
    !if (error /= 0) return
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_2D, count_2D, error, &
        & stride_2D,block_2D)
    !if (error /= 0) return

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    success = .true.

end subroutine hdf5_init_subset_2D


!> Add a dataset associated to a given dataspace
!!    @param[in]    local_size  = size of the local field
!!    @param[in]    global_size = size of the total grid
!!    @param[in]    file_id     = File identifier
!!    @param[in]    new         = is the data set new or does it already exist? (true for output, false for input)
!!    @param[in]    dsetname    = dataset name (= name of the field to save in the dataset)
!!    @param[out]   dset_id     = Dataset identifier
!!    @param[out]   plist_id    = Property list identifier
!!    @param[out]   success     = logical, true for success, false if an error occurs
subroutine hdf5_add_dataset(local_size, global_size, file_id, new, dsetname, dset_id, plist_id, success)

    ! Input/output
    integer(hsize_t), dimension(3), intent(in)  :: local_size   ! size of the local field
    integer(hsize_t), dimension(3), intent(in)  :: global_size  ! size of the global field
    integer(HID_T), intent(in)                  :: file_id      ! File identifier
    logical, intent(in)                         :: new          ! is the data set new or does it already exist?
    character(len=*), intent(in)                :: dsetname     ! dataset name (= name of the field to save in the dataset)
    integer(HID_T), intent(out)                 :: dset_id      ! Dataset identifier
    integer(HID_T), intent(out)                 :: plist_id     ! Property list identifier
    logical, intent(out)                        :: success

    ! Local variables
    integer                                 :: error        ! error flag
    integer                                 :: rank = 3     ! array rank
    integer(HID_T)                          :: filespace_loc! Dataspace identifier in file

    success = .false.

    ! Create the data space for the  dataset.
    call h5screate_simple_f(rank, global_size, filespace_loc, error)
    if (error /= 0) return

    ! Create chunked dataset.
    if(new) then
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
        !if (error /= 0) return
        call h5pset_chunk_f(plist_id, rank, local_size, error)
        !if (error /= 0) return
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace_loc, &
            & dset_id, error, plist_id)
        !if (error /= 0) return
        call h5sclose_f(filespace_loc, error)
        call h5pclose_f(plist_id, error)
        !if (error /= 0) return
    else
        !call h5dopen_f(file_id, dsetname, dset_id, error, plist_id)
        call h5dopen_f(file_id, dsetname, dset_id, error)
        !if (error /= 0) return
        call h5sclose_f(filespace_loc, error)
    end if
    !if (error /= 0) return

    success = .true.

end subroutine hdf5_add_dataset


!> Delete the parallel I/O context for a given geometry
!!    @param[out]   dset_id     = Dataset identifier
!!    @param[out]   filespace   = Dataspace identifier in file
!!    @param[out]   memspace    = Dataspace identifier in memory
!!    @param[out]   plist_id    = Property list identifier
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_delete_subset_3D(filespace, memspace, dset_id, plist_id) result(success)

    ! Input/output
    integer(HID_T), intent(in)      :: dset_id      ! Dataset identifier 
    integer(HID_T), intent(in)      :: filespace    ! Dataspace identifier in file 
    integer(HID_T), intent(in)      :: memspace     ! Dataspace identifier in memory
    integer(HID_T), intent(in)      :: plist_id     ! Property list identifier
    logical                         :: success

    ! Local variables
    integer                         :: error        ! error flag

    success = .false.

    ! Close dataspaces.
    call h5sclose_f(filespace, error)
    !if (error /= 0) return
    call h5sclose_f(memspace, error)
    !if (error /= 0) return
    ! Close the dataset.
    call h5dclose_f(dset_id, error)
    !if (error /= 0) return
    ! Close the property list.
    call h5pclose_f(plist_id, error)
    !if (error /= 0) return

    success = .true.

end function hdf5_delete_subset_3D



!> Delete the general context hdf5 context
!!    @return       success     = logical, true for success, false if an error occurs
function hdf5_close() result(success)

    ! Input/output
    logical                 :: success

    ! Local variables
    integer                 :: error                ! error flag

    success = .false.

    ! Close the property list.
    call h5pclose_f(plist_general_id, error)
    if (error /= 0) return

    ! Close FORTRAN interfaces and HDF5 library.
    CALL h5close_f(error)
    if (error /= 0) return

    success = .true.

end function hdf5_close


#endif

end module hdf5_io
!! @}
