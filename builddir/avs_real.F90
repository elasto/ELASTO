!USEFORTEST io
!------------------------------------------------------------------------------
! MODULE avs_real
!
!> @author
!> Patrick BEGOU, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to realize all parallel I/O from datalayout to avs files.
!> This is a generic version wich should be processed to create a REAL version and a COMPLEX version.
!> DO NOT CHANGE this file: it is automagicaly created from avs_implement.Fortran
!> for REAL and COMPLEX data types. All changes must be done in avs_implement.Fortran
!> file only.
!------------------------------------------------------------------------------
MODULE avs_real


    IMPLICIT NONE


CONTAINS


!------------------------------------------------------------------------------
!> This function dump data into a binary files that could be read by avs.
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function dump in a binary file the values of the REAL_DATA_LAYOUT 
!> variable: val. If the file exist, it is overwriten. The file format is
!> stream (C like binary format). It contains 3 integers, for the full array
!> dimension, followed by the datas.
! 
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @param[in] name the name of the output file.
!> @param[in] val the REAL_DATA_LAYOUT variable to save in the file.
!> @return .TRUE. if file is successfuly written.
!------------------------------------------------------------------------------
FUNCTION real_WriteToFile(spec_rank,name,val)
!
! Parallel write of a 3D scalar array in one file named 'name'. All the processes write in the same file.
! U, V, W and SC are known as pointer to 3D array decalred in "data" module.
!-------------------------------------------------------------------------

    USE datalayout
    USE communicators_tools
    USE mpi

    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: spec_rank
    CHARACTER(LEN=*), INTENT(IN)      :: name
    TYPE(REAL_DATA_LAYOUT),INTENT(IN) :: val
    LOGICAL                           :: real_WriteToFile

     
    INTEGER :: gsizes(3)
    !> the size of the global array in the 3 dimensions.
    INTEGER :: lsizes(3)
    !> the size of the local array in the 3 dimensions.
    INTEGER :: offset(3)
    !> the offset of the local array in the global array for the 3 dimensions.
    
    INTEGER :: iunit, fileview, ierr
    
    !Offset in the file for the data: 3 x 4 bytes integers.
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp=12     
    !> the offset (in bytes) of the datas in the file to skip the header.
    
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    
    real_WriteToFile=.FALSE.
    
    ! Create the view for parallel I/O
    gsizes(1) = val%nx
    gsizes(2) = val%ny
    gsizes(3) = val%nz
    lsizes(1) = val%xmax - val%xmin + 1
    lsizes(2) = val%ymax - val%ymin + 1
    lsizes(3) = val%zmax - val%zmin + 1
    offset(1) = val%xmin - 1 
    offset(2) = val%ymin - 1
    offset(3) = val%zmin - 1
    CALL MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,offset,MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": unable to create subarray!"
      RETURN
    ENDIF
    
    CALL MPI_TYPE_COMMIT(fileview,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": subarray type commit fails!"
      RETURN
    ENDIF
    
    ! Open file and go to the right position corresponding to our local subdomain
    CALL MPI_FILE_OPEN(spec_communicator,trim(name),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,iunit,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": cannot open file ["//trim(name)//"]"
      RETURN
    ENDIF
    
    ! The processus of rank 0 create the header
    IF (spec_rank .EQ. 0) THEN
       CALL MPI_FILE_WRITE(iunit,gsizes,3,MPI_INTEGER,status,ierr)
       IF (ierr .NE. MPI_SUCCESS) THEN
         WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
                           &": cannot write dimensions in file ["//trim(name)//"]"
         RETURN
       ENDIF
    ENDIF   
    

    ! Positionnement
    CALL MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
                        &": cannot set view on file ["//trim(name)//"]"
      RETURN
    ENDIF
    
    ! Write data
    CALL MPI_FILE_WRITE_ALL(iunit,val%values,lsizes(1)*lsizes(2)*lsizes(3),MPI_REAL_WP,status,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
                        &": write error for file ["//trim(name)//"]"
      RETURN
    ENDIF

    ! Close file
    CALL MPI_FILE_CLOSE(iunit,ierr)
    IF (ierr .NE. MPI_SUCCESS) THEN
      WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
                        &": close error on file ["//trim(name)//"]"
      RETURN
    ENDIF

    !All is OK now!
    real_WriteToFile=.TRUE.
    RETURN

end function real_WriteToFile

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU and Jean-Baptiste Lagaert
!
!> @details
!> This function get the size of the data field to read into an avs file.
! 
!> @param[in]   spec_rank   = the rank number of the process in the range [0,ncpus - 1].
!> @param[in]   name        = the name of the output file.
!> @param[out]  mesh_size   = number of point along each direction
!> @param[out]  success     = .TRUE. if size is obtain
!------------------------------------------------------------------------------
subroutine real_GetSizeFromFile(spec_rank,name,mesh_size, success)
!
! Parallel read of a 3D scalar array in one file named 'name'. All the processes write in the same file.
! U, V, W and SC are pointers to 3D arrays defined in "data" module.
!-------------------------------------------------------------------------

    use communicators_tools
    use mpi

    ! Input/Output
    integer, intent(in)                 :: spec_rank
    CHARACTER(len=*), intent(in)        :: name
    integer, dimension(3),intent(out)   :: mesh_size
    logical                             :: success

    ! Local variables
    integer :: iunit, ierr
    integer, dimension(MPI_STATUS_SIZE) :: status
    
    ! ===== Initialisation =====
    success=.FALSE.
    
    ! ===== Read info about file data (mesh and co.) 
    ! Open file and go to the right position corresponding to our local subdomain
    call MPI_FILE_OPEN(spec_communicator,trim(name),MPI_MODE_RDONLY,MPI_INFO_NULL,iunit,ierr)
    if (ierr .NE. MPI_SUCCESS) then
        write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,": cannot open file ["//trim(name)//"]"
        return
    end if
    
    ! Reading header (which contians some mesh information)
    call MPI_FILE_READ(iunit,mesh_size,3,MPI_integer,status,ierr)
    if (ierr .NE. MPI_SUCCESS) then
        write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,&
                &": cannot read dimensions in file ["//trim(name)//"]"
        return
    end if

    ! ==== Close file and return result =====
    call MPI_FILE_CLOSE(iunit,ierr)
    if (ierr .NE. MPI_SUCCESS) then
      write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,&
                        &": close error on file ["//trim(name)//"]"
      return
    end if

    !All is OK now!
    success=.true.

end subroutine real_GetSizeFromFile

  
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU and Jean-Baptiste Lagaert
!
!> @details
!> This function dump in a binary file the values of the REAL_DATA_LAYOUT 
!> variable: val. If the file exist, it is overwriten. The file format is
!> stream (C like binary format). It contains 3 integers, for the full array
!> dismension, followed by the datas.
! 
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @param[in] name the nameof the output file.
!> @param[in] val the REAL_DATA_LAYOUT variable to save in the file.
!> @return .TRUE. if file is successfuly written.
!------------------------------------------------------------------------------
function real_ReadFromFile(spec_rank,name,val) result(success)
!
! Parallel read of a 3D scalar array in one file named 'name'. All the processes write in the same file.
! U, V, W and SC are pointers to 3D arrays defined in "data" module.
!-------------------------------------------------------------------------

    use datalayout
    use communicators_tools
    use mpi

    integer, intent(in)                  :: spec_rank
    CHARACTER(len=*), intent(in)         :: name
    type(REAL_DATA_LAYOUT),intent(inout) :: val
    logical                              :: success

    !> the size of the global array in the 3 dimensions.
    integer :: gsizes(3)
    !> the size of the local array in the 3 dimensions.
    integer :: lsizes(3)
    !> the offset of the local array in the global array for the 3 dimensions.
    integer :: offset(3)
    !> the buffer to read the  3 dimensions of the data in the file.
    integer, DIMENSION(3) ::nxyz
    integer :: iunit, fileview, ierr
    !Offset in the file for the data: 3 x 4 bytes integers.
    integer(KIND=MPI_OFFSET_KIND) :: disp=12     
    !> the offset (in bytes) of the datas in the file to skip the header.
    
    integer, dimension(MPI_STATUS_SIZE) :: status
    
    ! ===== Initialisation =====
    success=.FALSE.
    
    ! ===== Read info about file data (mesh and co.) 
    ! Open file and go to the right position corresponding to our local subdomain
    call MPI_FILE_OPEN(spec_communicator,trim(name),MPI_MODE_RDONLY,MPI_INFO_NULL,iunit,ierr)
    if (ierr .NE. MPI_SUCCESS) then
        write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,": cannot open file ["//trim(name)//"]"
        return
    end if
    
    ! Reading header (which contians some mesh information)
    call MPI_FILE_READ(iunit,nxyz,3,MPI_integer,status,ierr)
    if (ierr .NE. MPI_SUCCESS) then
        write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,&
                &": cannot read dimensions in file ["//trim(name)//"]"
        return
    end if
    
    ! ===== Read field =====
    ! == If the dimensions agree with val structure, just read data ... ==
    if (val%nx == nxyz(1) .OR. val%nz == nxyz(2) .OR. val%nz == nxyz(3)) then
        ! -- Create the view for parallel I/O --
        gsizes(1) = val%nx
        gsizes(2) = val%ny
        gsizes(3) = val%nz
        lsizes(1) = val%xmax - val%xmin + 1
        lsizes(2) = val%ymax - val%ymin + 1
        lsizes(3) = val%zmax - val%zmin + 1
        offset(1) = val%xmin - 1 
        offset(2) = val%ymin - 1
        offset(3) = val%zmin - 1
        call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,offset,MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview,ierr)
        if (ierr .NE. MPI_SUCCESS) then
            write(6,'(a,i0,a)')"ERROR real_ReadfromFile on process ",spec_rank,": unable to create subarray!"
            return
        end if
        call MPI_TYPE_COMMIT(fileview,ierr)
        if (ierr .NE. MPI_SUCCESS) then
            write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,": subarray type commit fails!"
            return
        end if
        ! -- Positionnement --
        call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
        if (ierr .NE. MPI_SUCCESS) then
          write(6,'(a,i0,a)')"ERROR real_ReadfromFile on process ",spec_rank,&
                            &": cannot set view on file ["//trim(name)//"]"
          return
        end if
        ! -- Read data --
        call MPI_FILE_READ_ALL(iunit,val%values,lsizes(1)*lsizes(2)*lsizes(3),MPI_REAL_WP,status,ierr)
        if (ierr .NE. MPI_SUCCESS) then
          write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,&
                            &": read error for file ["//trim(name)//"]"
          return
        end if

        ! == ... if it is not the case, the read interface does not do the job !! ==
    else
        if(spec_rank.eq.0) write(6,'(a)') '[ERROR] Read request on wrong mesh size for file '//trim(name)
    end if

    ! Close file
    call MPI_FILE_CLOSE(iunit,ierr)
    if (ierr .NE. MPI_SUCCESS) then
      write(6,'(a,i0,a)')"ERROR success on process ",spec_rank,&
                        &": close error on file ["//trim(name)//"]"
      return
    end if

    !All is OK now!
    success=.true.
    return

end function real_ReadFromFile


end module avs_real
