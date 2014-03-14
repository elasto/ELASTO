!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: xdmf
!
!
! DESCRIPTION:
!> This module provide all procedure needed to add xdmf descriptor to h5 output
!! in order to visualise data containing into hdf5 files (*.h5).
!! hdf5 format.
!
!! @details
!!        This module developp tools to visualize any H5 files (hdf5 format) write
!!    by SCALES. It reach this purpose by creating an xdmf descriptor of the data.
!!        More precisly it means adding an *.xdmf file wich describes the (cartesian
!!    3D-)meshes and listing the data associated to each mesh. Note that data
!!    are not re-write and thus xmdf files is very short (few lines). These
!!    files are "human-readable" (few lines of xml).
!!        This way, data can be easily be visualize with classical tools (eg
!!    Paraview, VisIt).
!!        The add of descriptor can be done at post-process step on any single
!!    core. This step is very short and require few memory (even for big data) ; thus
!!    it could be performed sequentially on any small computer.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module xdmf

    use precision_tools

    implicit none

!    ! ===== Private attribut =====
!    integer,dimension(:,:),allocatable,save     :: meshes

    ! ===== Private member =====
    public      :: xdmf_write           ! write the wanted xdmf descriptor

    ! ===== Private member =====
    private     :: xdmf_init_file       ! create and init a new xdmf file
    private     :: xdmf_close_file      ! close a xdmf file
#ifdef HDF5
    private     :: xdmf_add_mesh_fields ! add all the fiels associate to a given mesh
#endif

contains

!> Write the wanted XDMF descriptor
!!    @param[in]    filename    = name of the xdmf file
!!    @return       success     = logical, true for success, false if an error occurs
subroutine xdmf_write(filename, file_id, success)

    use parser_tools

    ! Input/Output
    character(len=*), intent(in)    :: filename
    logical                         :: success

    ! Local variables
    integer                         :: file_id  ! file identifier
    integer                         :: mesh_id  ! mesh identifier (for do loop)
    integer                         :: mesh_nb  ! number of meshes
    integer                         :: field_id ! field identifier (for do loop)
    integer                         :: field_nb ! number of field for a given mesh
    integer                         :: ires     ! value for allocation success 
    character(len=str_medium), dimension(:,:), allocatable  :: h5_Sfield
    character(len=str_medium)       :: tag, tag_bis ! tag in input file

    success = .false.

#ifdef HDF5
    ! Create file and init it
    call xdmf_init_file(filename, file_id, success)

    ! Get number of meshes
    tag = "xdmf - number of mesh"
    if(parser_is_defined(tag)) then
        call parser_read(tag, mesh_nb)
    else
        write(6,'(a)')"ERROR for creating "//trim(filename)//": cannot get number of mesh from parser"
        return
    end if

    ! Add field for each mesh
    do mesh_id = 1, mesh_nb
        ! Get number of fields
        write(tag,'(a,i0)') 'xdmf - number of field for mesh ',mesh_id
        if(parser_is_defined(tag)) then
            call parser_read(tag, field_nb)
        else
            write(6,'(a,i0)')"ERROR for creating "//trim(filename)//": cannot get number of field from parser for mesh ", mesh_id
            return
        end if
        allocate(h5_Sfield(2,field_nb),stat=ires)
        if (ires .ne. 0) return
        ! Read fiels info (field name and file name)
        do field_id = 1, field_nb
            write(tag,'(2(a,i0),a)') 'xdmf - mesh ', mesh_id, ' field ', field_id, ' file'
            if(parser_is_defined(tag)) then
                call parser_read(tag, h5_Sfield(1,field_id))
            else
                write(6,'(3(a,i0))')"ERROR for creating "//trim(filename)//": cannot get file name for field ", &
                    field_id, " in mesh ", mesh_id
                return
            end if
            write(tag,'(2(a,i0),a)') 'xdmf - mesh ', mesh_id, ' field ', field_id, ' name'
            if(parser_is_defined(tag)) then
                call parser_read(tag, h5_Sfield(2,field_id))
            else
                write(6,'(3(a,i0))')"ERROR for creating "//trim(filename)//": cannot get field name for field ", &
                    field_id, " in mesh ", mesh_id
                return
            end if
        end do
        ! And write all in the XDMF file
        call  xdmf_add_mesh_fields(file_id, h5_Sfield, success)
        if (.not. success) then
            write(6,'(a,i0)')"ERROR for creating "//trim(filename)//": cannot not describe mesh number ", mesh_id
            return
        end if
        deallocate(h5_Sfield)
    end do

    ! Close all
    success = xdmf_close_file(file_id)
    if (.not. success) then
        write(6,'(a,i0)')"ERROR for "//trim(filename)//": cannot not close file environment."
        return
    end if

#else
    write(6,'(a)')"ERROR for writing XDMF descriptor for H5 file : not available H5 lib !!"

#endif

end subroutine xdmf_write


!> Create a new xdmf file and write the generic header of such a file
!!    @param[in]    filename    = name of the xdmf file
!!    @param[in]    file_id     = file identifier for future writes and to close it
!!    @return       success     = logical, true for success, false if an error occurs
subroutine xdmf_init_file(filename, file_id, success)

    use fileio

    character(len=*), intent(in)    :: filename
    integer, intent(out)            :: file_id
    logical                         :: success

    success = .false.

    ! Open file
    file_id = iopen()
    open(unit = file_id, file=trim(filename), form='formatted')

    write(file_id, '(a)') '<?xml version="1.0" ?>'
    write(file_id, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(file_id, '(a)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
    write(file_id, '(a)') '  <Domain>'

    success = .true.

end subroutine xdmf_init_file

#ifdef HDF5

!!    @param[in]    file_id     = file identifier for future writes and to close it
!!    @param[in]    h5_file     = list of hdf5 files with data at the same resolution
!!    @return       success     = logical, true for success, false if an error occurs
subroutine xdmf_add_mesh_fields(file_id, h5_Sfield, success)

    use hdf5
    use precision_tools

    !Input/Output
    integer, intent(in)                                     :: file_id
    logical, intent(out)                                    :: success
    character(len=str_medium), dimension(:,:), intent(in)   :: h5_Sfield

    ! Other local variables
    integer                                     :: hdferr       ! error flag
    integer                                     :: id_field     ! iterator on fields
    integer                                     :: mesh_rank    ! to know if fields are in 2D or 3D
    integer(hsize_t), dimension(:), allocatable :: mesh_size    ! number of cell in the mesh
    real(WP), dimension(:), allocatable         :: mesh_space_step  ! spatial step
    integer(hsize_t), dimension(:), allocatable :: maxdims      ! maximum size of H5 fields
    integer                                     :: rank         ! dimension of the current h5 field
    integer                                     :: ires         ! value for allocation success
    integer(hsize_t), dimension(:), allocatable :: dims         ! size of the current h5 field
    integer(HID_T)                              :: h5file_id    ! H5 file identifier for H5 routines
    integer(HID_T)                              :: dataset_id   ! Dataset identifier
    integer(HID_T)                              :: dataspace_id ! Dataspace identifier in file

    ! Init
    success = .false.

    ! ===== Get mesh size from the first field ====
    ! -- Get access to the field --
    ! Open field and dataset
    call h5fopen_f(trim(h5_Sfield(1,1)), H5F_ACC_RDONLY_F, h5file_id, hdferr)
    if (hdferr == -1) then
        write(6,'(a)')"ERROR for creating XDMF descriptor: cannot open file ["//trim(h5_Sfield(1,1))//"]"
        return
    end if
    call h5dopen_f(h5file_id, trim(h5_Sfield(2,1)), dataset_id, hdferr)
    if (hdferr == -1) then
        write(6,'(a)')"ERROR for creating XDMF descriptor: cannot open dataset ["//trim(h5_Sfield(2,1))//"]"
        return
    end if
    ! Get space
    call h5dget_space_f(dataset_id, dataspace_id, hdferr)
    call h5dclose_f(dataset_id, hdferr)
    call h5sget_simple_extent_ndims_f(dataspace_id, mesh_rank, hdferr)
    allocate(mesh_size(mesh_rank),stat=ires)
    if (ires .ne. 0) return
    allocate(maxdims(mesh_rank),stat=ires)
    if (ires .ne. 0) return
    allocate(dims(mesh_rank),stat=ires)
    if (ires .ne. 0) return
    call h5sget_simple_extent_dims_f(dataspace_id, mesh_size, maxdims, hdferr)
    call h5sclose_f(dataspace_id, hdferr)
     call h5fclose_f(h5file_id, hdferr)
 
     ! Get the space step if attribute is present
     ! Get time
     ! Close file
     ! Write header
     write(file_id,'(a)')          '    <Grid Name="Structured Grid" GridType="Uniform">'
     if (mesh_rank ==3) then
         write(file_id,'(a,3(i0,x),a)')'      <Topology TopologyType="3DCORECTMesh" NumberOfElements="', mesh_size,'"/>'   ! mesh = 3D cartesian one
         write(file_id,'(a)')          '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
         write(file_id,'(a)')          '          <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">'
         write(file_id,'(a)')          '            0 0 0'
         write(file_id,'(a)')          '          </DataItem>'
         write(file_id,'(a)')          '          <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">'
         write(file_id,'(12x,3(i0,x))')            mesh_space_step
     elseif(mesh_rank==2) then
         write(file_id,'(a,2(i0,x),a)')'      <Topology TopologyType="2DCORECTMesh" NumberOfElements="', mesh_size,'"/>'   ! mesh = 3D cartesian one
         write(file_id,'(a)')          '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
         write(file_id,'(a)')          '          <DataItem Dimensions="2 " NumberType="Float" Precision="4" Format="XML">'
         write(file_id,'(a)')          '            0 0'
         write(file_id,'(a)')          '          </DataItem>'
         write(file_id,'(a)')          '          <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">'
         write(file_id,'(12x,2(i0,x))')            mesh_space_step
     else
         write(*,'(a,a,i0,a)') ' [ERROR] In xdmf descriptor creation: '//trim(h5_Sfield(2,1))//' in ' ,&
             & trim(h5_Sfield(1,1))//' is a ',mesh_rank,'D array (and not a 2D or a 3D).'
         return
     end if
     write(file_id,'(a)')          '          </DataItem>'
     write(file_id,'(a)')          '     </Geometry>'
 
     ! ===== Then check size of the other fields and add them =====
     do id_field = 1, size(h5_Sfield,1)
         ! -- Open field and dataset --
         call h5fopen_f(trim(h5_Sfield(1,id_field)), H5F_ACC_RDONLY_F, h5file_id, hdferr)
         if (hdferr == -1) then
             write(6,'(a)')"ERROR for creating XDMF descriptor: cannot open file ["//trim(h5_Sfield(1,id_field))//"]"
             return
         end if
         call h5dopen_f(h5file_id, trim(h5_Sfield(2,id_field)), dataset_id, hdferr)
         if (hdferr == -1) then
             write(6,'(a)')"ERROR for creating XDMF descriptor: cannot open dataset ["//trim(h5_Sfield(2,id_field))//"]"
             return
         end if
 
         ! -- Check size --
         call h5dget_space_f(dataset_id, dataspace_id, hdferr)
         call h5dclose_f(dataset_id, hdferr)
         call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
         if (rank/=mesh_rank) then
             write(*,'(a,a,i0,a)') ' [ERROR] In xdmf descriptor creation: '//trim(h5_Sfield(2,id_field))//' in ' ,&
                 & trim(h5_Sfield(1,id_field))//' is not a ',mesh_rank,'D array.'
             return
         end if
         call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
         if (sum(abs(dims-mesh_size))>=1) then
             write(*,'(a,a,3(i0,x),a,3(i0,x))') ' [ERROR] In xdmf descriptor creation: '//trim(h5_Sfield(2,id_field))//' in ' ,&
                 & trim(h5_Sfield(1,id_field))//' size (',dims,') is not the same that current mesh (',mesh_size,').'
             return
         end if
 
         ! -- Close h5 file --
         call h5sclose_f(dataspace_id, hdferr)
         call h5fclose_f(h5file_id, hdferr)
 
         ! -- Write field info in the XDMF descriptor --
         write(file_id, '(a)') '      <Attribute Name="'//trim(h5_Sfield(2,id_field))//'" AttributeType="Scalar" Center="Node">'
         write(file_id, '(a,3(i0,x),a)') '        <DataItem Dimensions="', mesh_size, &
             & '" NumberType="Float" Precision="8" Format="HDF">'//trim(h5_Sfield(1,id_field))//':/'//trim(h5_Sfield(2,id_field))//'</DataItem>'
         write(file_id, '(a,3(i0,x),a)') '     </Attribute>'
     end do
 
     ! Free memory
     deallocate(mesh_size)
     deallocate(mesh_size)
     deallocate(mesh_size)
 
 
     success = .true.
 
 end subroutine xdmf_add_mesh_fields
 
#endif
 
 !> Close all xdmf environnement and then close the xmdf file
 !!    @param[in]    file_id     = file identifier for future writes and to close it
 !!    @return       success     = logical, true for success, false if an error occurs
 function xdmf_close_file(file_id) result(success)
 
     use fileio
 
     integer, intent(in)             :: file_id
     logical                         :: success
 
     success = .false.
 
     ! Close environement
     write(file_id, '(a)') '    </Grid>'
     write(file_id, '(a)') '  </Domain>'
     write(file_id, '(a)') '</Xdmf>'
 
     ! Close file
     close(file_id)
     if(iclose(file_id) == -1) return
 
     success = .true.
 
 end function xdmf_close_file
 
 
!function xdmf_init_in_simu(ite,U,V,W,Scal,ScalPart,Bx,By,Bz) result(success)
!
!    use datalayout
!
!    !I/O data
!    type(real_data_layout),intent(in),optional              :: U
!    type(real_data_layout),intent(in),optional              :: V
!    type(real_data_layout),intent(in),optional              :: W
!    type(real_data_layout),intent(in),dimension(:),optional :: Scal
!    type(real_data_layout),intent(in),dimension(:),optional :: ScalPart
!    type(real_data_layout),intent(in),optional              :: Bx
!    type(real_data_layout),intent(in),optional              :: By
!    type(real_data_layout),intent(in),optional              :: Bz
!    integer,intent(in)                                      :: ite
!    logical                                                 :: success
!
!    !local data
!    integer,dimension(:,:),allocatable                      :: temptab
!    integer,dimension(:,:),allocatable                      :: MeshLoc 
!    integer                                                 :: nbMeshCurrent
!    integer                                                 :: ires 
!    integer                                                 :: i,j,k,oldsize 
!    logical                                                 :: newVal 
!
!    success = .false.
!    nbMeshCurrent = 0
!    !for Velocity
!    if (present(U) .and. present(V) .and. present(W)) then
!      nbMeshCurrent = nbMeshCurrent + 1 
!      allocate(MeshLoc(nbMeshCurrent,3),stat=ires)
!      if (ires .eq. 0) return
!      MeshLoc(nbMeshCurrent,1)=U%nx
!      MeshLoc(nbMeshCurrent,2)=U%ny
!      MeshLoc(nbMeshCurrent,3)=U%nz
!    endif
!    !Same for Scal
!    if (present(Scal)) then
!      nbMeshCurrent = nbMeshCurrent + size(Scal,1) 
!      oldsize = 0
!      if ( allocated(MeshLoc)) then
!        oldsize = size(MeshLoc,1)
!        allocate(temptab(oldsize,3),stat=ires)
!        if (ires .eq. 0) return
!        temptab = MeshLoc
!        deallocate(MeshLoc)
!        allocate(MeshLoc(nbMeshCurrent,3),stat=ires)
!        if ( oldsize .ne. 0) then
!          do i=1,oldsize
!            do j=1,3
!              MeshLoc(i,j)=temptab(i,j)
!            enddo
!          enddo
!        endif
!        deallocate(temptab)
!      else
!        j=1
!        do i=oldsize+1,nbMeshCurrent
!          MeshLoc(i,1)=Scal(j)%nx
!          MeshLoc(i,2)=Scal(j)%ny
!          MeshLoc(i,3)=Scal(j)%nz
!          j=j+1
!        enddo
!      endif
!    endif
!    !Same for ScalPart
!    if (present(ScalPart)) then
!      nbMeshCurrent = nbMeshCurrent + size(ScalPart,1) 
!      oldsize = 0
!      if ( allocated(MeshLoc)) then
!        oldsize = size(MeshLoc,1)
!        allocate(temptab(oldsize,3),stat=ires)
!        if (ires .eq. 0) return
!        temptab = MeshLoc
!        deallocate(MeshLoc)
!        allocate(MeshLoc(nbMeshCurrent,3),stat=ires)
!        if ( oldsize .ne. 0) then
!          do i=1,oldsize
!            do j=1,3
!              MeshLoc(i,j)=temptab(i,j)
!            enddo
!          enddo
!        endif
!        deallocate(temptab)
!      else
!        j=1
!        do i=oldsize+1,nbMeshCurrent
!          MeshLoc(i,1)=ScalPart(j)%nx
!          MeshLoc(i,2)=ScalPart(j)%ny
!          MeshLoc(i,3)=ScalPart(j)%nz
!          j=j+1
!        enddo
!      endif
!    endif
!    !Same for B mais je ne sais pas si Mouloud le merite !!!!!! 
!    if (present(Bx) .and. present(By) .and. present(Bz) ) then
!      nbMeshCurrent = nbMeshCurrent + 1 
!      oldsize = 0
!      if ( allocated(MeshLoc)) then
!        oldsize = size(MeshLoc,1)
!        allocate(temptab(oldsize,3),stat=ires)
!        if (ires .eq. 0) return
!        temptab = MeshLoc
!        deallocate(MeshLoc)
!        allocate(MeshLoc(nbMeshCurrent,3),stat=ires)
!        if ( oldsize .ne. 0) then
!          do i=1,oldsize
!            do j=1,3
!              MeshLoc(i,j)=temptab(i,j)
!            enddo
!          enddo
!        endif
!        deallocate(temptab)
!      else
!        MeshLoc(nbMeshCurrent,1)=Bx%nx
!        MeshLoc(nbMeshCurrent,2)=Bx%ny
!        MeshLoc(nbMeshCurrent,3)=Bx%nz
!      endif
!    endif
!    
!    if (nbMeshCurrent .lt. 0) then    
!      j=1
!      allocate(tempTab(size(MeshLoc,1),size(MeshLoc,2)),stat=ires)
!      if (ires .ne. 0 ) return 
!      tempTab = 0 
!      do i=1,nbMeshCurrent 
!        if ( i.eq. 1 ) then
!          tempTab(j,:)=MeshLoc(i,:)
!        else
!          newVal = .false.
!          do k=1,j
!            if ( MeshLoc(i,1) .ne. tempTab(k,1) .and.&
!               & MeshLoc(i,2) .ne. tempTab(k,2) .and.&
!               & MeshLoc(i,3) .ne. tempTab(k,3)) then
!              newVal = .true.
!            endif
!          enddo
!          if ( newVal ) then
!            j = j + 1
!            tempTab(j,:) = MeshLoc(i,:)
!        endif
!      enddo
!      deallocate(MeshLoc)
!      allocate(meshes(j,3),stat=ires)
!      if (ires .ne. 0 ) return 
!      meshes(i,:)=tempTab(i,:)
!      deallocate(tempTab)
!    else
!      write (6,'(a)') '[WARNING] Aucun champs detecte par xdmf_init_in_simu'
!    endif
!    success=.true.
!
!end function xdmf_init_in_simu
!
! !> Write the wanted XDMF descriptor
! !!    @param[in]    filename    = name of the xdmf file
! !!    @return       success     = logical, true for success, false if an error occurs
! function xdmf_write_in_simu(U,V,W,Scal,ScalPart,B,ite,success)
! 
!     !I/O data
!     type(real_data_layout),intent(in),optional              :: U
!     type(real_data_layout),intent(in),optional              :: V
!     type(real_data_layout),intent(in),optional              :: W
!     type(real_data_layout),intent(in),dimension(:),optional :: Scal
!     type(real_data_layout),intent(in),dimension(:),optional :: ScalPart
!     type(real_data_layout),intent(in),optional              :: Bx
!     type(real_data_layout),intent(in),optional              :: By
!     type(real_data_layout),intent(in),optional              :: Bz
!     integer,intent(in)                                      :: ite
!     logical                                                 :: success
! 
!     success = .false.
! 
! #ifdef HDF5
! !nombre de grilles différentes sur tous les champs -> nbMesh
! 
! #else
!     write(6,'(a)')"ERROR for writing XDMF descriptor for H5 file : not available H5 lib !!"
! 
! #endif
! 
! end function xdmf_write_in_simu


end module xdmf
!! @}
