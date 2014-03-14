!USEFORTEST io
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: io_interface
!
!
! DESCRIPTION:
!> This module provide all procedure needed to io. It is designed to be used as
!! wrapper/interface to all the other io procedures: the other part of
!! the code are supposed - in future at least - only include this module.
!!
!! output at vtk format.
!
!! @details
!!        This module developp tools to write distribued output. This means that
!!    output of one field will be done in one file per (mpi) processus. These allows
!!    to write field computed with an high resolution without be limitated by the
!!    size of the file output and to avoid to big loading time during visualisation
!!    or file loading in order to initialize computation at a given setup.
!
!!         This first version provide only output tools. Some input procedures could
!!    be add in future works. The general context (number of field to save, physical
!!    dimension, number of processus, ...) is initiliazed by calling "parallel_io_init_all"
!!    the context specific to each field (mesh resolution, number of point, name of
!!    the ouput, information about time sequence, ...) is initialized or save by
!!    calling "parallel_io_init_field". After that, call "parallel_write" in order to
!!    create a new output of the field.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module io_interface

    use precision_tools

    implicit none

    ! ===========================================
    ! =====         Public procedures       =====
    ! ===========================================

    ! ===== Generic interfaces for input/output =====


    ! ===========================================
    ! =====         Private variables       =====
    ! ===========================================
    !> Init status
    logical, private                    :: init = .false.
    ! ==== I/O Format =====
    character(len=str_short), private   :: input_format = 'avs'
    character(len=str_short), private   :: output_format = 'avs'
    logical , private                   :: output_2Dslice= .false.
    !> write format for output name of *.hdf5
    character (len=100), private        :: avs_format_name

    contains

!> This subroutine init the input/output context
!! and avs files.
!! @param[in]   spec_rank   = the rank number of the process in the range [0,ncpus - 1].
!! @param[out]  success     = .TRUE. if size is obtain
function interface_io_init(spec_rank) result(success)

    use parser_tools
    use hdf5_io
    use communicators_tools

    integer, intent(in)     :: spec_rank
    logical                 :: success

    ! local variable
    integer     :: nb_ite   ! total number of iteration
    integer     :: size_ite ! number of character used to write iteration number

    success = .false.

    ! Set i/o format
    if(parser_is_defined("output format")) call parser_read("output format", output_format)
    if(parser_is_defined("output slice2D")) call parser_read("output slice2D", output_2Dslice)
    if(parser_is_defined("input format")) call parser_read("input format", input_format)

    ! Define output name format
    if(parser_is_defined("Simulation iterations")) then
        call parser_read("Simulation iterations", nb_ite)
    else
        nb_ite = 1
    end if

    ! Some specific intialisation depending on output format
    if ( trim(output_format) == "hdf5" .or. trim(input_format) == "hdf5" ) then
#ifdef HDF5
        if(.not. hdf5_init(spec_communicator, nb_ite)) then
            if(spec_rank==0) write(*,'(a,i0)') '[ERROR] Unable to init input/output on processus '&
                & , spec_rank
        end if
#else
        if(spec_rank==0) write(*,'(a,i0)') '[ERROR] Hdf5 output can not be performed : not hdf5 libraray detected !!'
        stop
#endif
    endif
    if ( trim(output_format) == "avs" .or. trim(input_format) == "avs" ) then
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
        write(avs_format_name,'(a,i0,a1,i0,a)') '(a,a,i',size_ite,'.',size_ite,')'
    endif

    init = .true.
    ! === Write io parameters ===
    if(spec_rank==0) then
        write(*,'(a,a)') '    [IO] output format is ',trim(output_format)
        write(*,'(a,L6)') '    [IO] output are 2D slices ',output_2Dslice
        write(*,'(a,a)') '    [IO] input  format is ',trim(input_format)
    end if

    success = .true.

end function interface_io_init

!> This subroutine prepare the income of a global input wrapper used for hdf5
!! and avs files.
!! @param[in]   spec_rank   = the rank number of the process in the range [0,ncpus - 1].
!! @param[in]   name        = the name of the output file.
!! @param[out]  success     = .TRUE. if size is obtain
!! @details
!! This first version just allow avs read. It allows to get the rediscretization
!! outside of the read it self (ie no recopy of interpolation in avs AND hdf5
!! read)
subroutine read_datalayout_scalar(spec_rank,name,val, success, fieldname)

    use avs_real
    use hdf5_wrapper_datalayout
    use datalayout
    use mpilayout_tools
    use rediscretization_tools

    ! Input/Output
    integer, intent(in)                     :: spec_rank
    CHARACTER(len=*), intent(in)            :: name
    type(REAL_DATA_LAYOUT),intent(inout)    :: val
    logical, intent(out)                    :: success
    CHARACTER(len=*), intent(in), optional  :: fieldname

    ! Local variables
    integer, dimension(3)               :: mesh_size    ! discretisation of the readed field
    type(REAL_DATA_LAYOUT)              :: read_tmp     ! field on "wrong" input discrectization
    type(COMPLEX_DATA_LAYOUT)           :: read_tmp_spec! field on "wrong" input discrectization
    integer                             :: nbcpus       ! number of cpus inside the spectral communicator
    logical                             :: res, existe
    character(len=str_medium)           :: fieldname_bis

    ! ===== Init =====
    success = .true.
    if (present(fieldname)) then
        fieldname_bis = trim(fieldname)
    else
        fieldname_bis = trim(val%name)
    end if

    if(.not.init) then
        if(.not.interface_io_init(spec_rank)) return
    end if

    ! ===== Check if file exist =====
    inquire( file=trim(adjustl(name)), exist=existe)
    if (.not. existe) then
      write(6,'(a,x,a,x,a,1x,i0)')'[ERROR] unable to find input file', &
        & trim(name), 'on rank', spec_rank
      return
    end if


    if(io_checkSizeFromFile(spec_rank,name, fieldname_bis, val, mesh_size)) then
        ! ===== Read data directly =====
        if(trim(input_format)=="hdf5") then
            call hdf5_read_scalar_wrapper(val, name, res, fieldname_bis, spec_rank)
        else
            res = real_ReadFromFile(spec_rank,name,val)
        end if
        if (.not. res) success = .false.
    else
        ! ===== Read data into a temp field and interpolate it to the right resolution =====
        if(spec_rank.eq.0) write(6,'(6(a,i0),a)')'[WARNING] success: Data in file '//trim(name)//&
           & ' [',mesh_size(1),',',mesh_size(2),',',mesh_size(3),                                      &
           & '] do not agree with the expected dimensions and will be interpolated'//   &
           & ' on a [',val%nx,',',val%ny,',',val%nz,'] grid.'
        ! -- Read data in a temp field --
        ! Init read_tmp
        nbcpus = getnbcpus()
        if (.not. initDataLayout(val%name,read_tmp,mesh_size(1),mesh_size(2),mesh_size(3),val%Lx,val%Ly,val%Lz,&
                & nbcpus,spec_rank)) then
            write(6,'(a,1x,a,1x,a,1x,i0)')'[ERROR] unable to create read temp field to interpolate file data for', &
                & trim(name),'on',spec_rank
            success = .false.
        end if
        ! Read data from file to read_tmp
        if(trim(input_format)=="hdf5") then
            call hdf5_read_scalar_wrapper(read_tmp, name, res, fieldname_bis, spec_rank)
        else
            res = real_ReadFromFile(spec_rank,name,read_tmp)
        end if
        if (.not. res) success = .false.
        ! -- Go in spectral space --
        ! Allocate storage for field in fourrier space
        if (.not. initDataLayout("read_tmp_spec", read_tmp_spec,(read_tmp%nx/2)+1, &
                & read_tmp%ny,read_tmp%nz, read_tmp%Lx,read_tmp%Ly,read_tmp%Lz, &
                & nbcpus,spec_rank,alongZ)) then
            write(6,'(a,x,a,x,i0)') "[ERROR] unable to create spectral read_tmp field to interpolate file data for",&
                & trim(name), 'on', spec_rank
            success = .false.
        end if
        ! Obtain spectral field
        call ftran(read_tmp, read_tmp_spec, success)
        if (.not.success) then
            write(6,'(a,x,i0)')'[ERROR] Unable to perform fft for tmp_read on rank', spec_rank
            success = .false.
        end if
        ! Real field is no more needed
        call deleteDataLayout(read_tmp)
        ! -- And perform the interpolation --
        if (.not. transdiscretization (read_tmp_spec, val, val)) then
            write(6,'(a,x,a,x,a,x,i0)')'[ERROR] Unable to perform spectral interpolation for file', &
            trim(name), 'on rank', spec_rank
            success = .false.
        end if
        ! -- Free memory --
        call deleteDataLayout(read_tmp_spec)
    end if  ! if test for mesh_size(field) = mesh_size(file)

end subroutine read_datalayout_scalar


!> Check if data present inside a file use the same discretisation (ie the same mesh) than a field
function io_checkSizeFromFile(spec_rank,filename,fieldname, layout, mesh_size) result(same_size)

    use datalayout
    use avs_real
    use hdf5_wrapper_datalayout

    ! Input/Output
    integer, intent(in)                     :: spec_rank
    CHARACTER(len=*), intent(in)            :: filename
    CHARACTER(len=*), intent(in)            :: fieldname
    type(REAL_DATA_LAYOUT),intent(inout)    :: layout
    integer, dimension(3), intent(out)      :: mesh_size    ! discretisation of the readed field
    logical                                 :: same_size

    ! Local variables
    logical                     :: success


    same_size=.false.

    if (modulo(layout%nx*layout%ny*layout%nz,getnbcpus()) .ne. 0 .and. trim(input_format)=="hdf5") then
      write (6,'(a)') '[INFO HDF5] number of grid points is not balanced on each processor: HDF5 not implemented !!!'
    endif
    if(trim(input_format)=="hdf5") then
#ifdef HDF5
        call hdf5_GetSizeFromFile(spec_rank,filename,fieldname,mesh_size, success)
#endif
    else
        call real_GetSizeFromFile(spec_rank,filename,mesh_size,success)
    end if

    if(success) same_size = (layout%nx == mesh_size(1)).and. &
                          & (layout%ny == mesh_size(2)).and. &
                          & (layout%nz == mesh_size(3))

end function io_checkSizeFromFile

!> Write a single scalar into a field
function write_datalayout_scalar(layout,iteration, spec_rank, postname) result(success)

    use datalayout
    use avs
    use hdf5_wrapper_datalayout

    ! Input/Output
    integer, intent(in)                     :: spec_rank
    integer, intent(in)                     :: iteration
    type(REAL_DATA_LAYOUT),intent(in)       :: layout
    logical                                 :: success
    character(len=*), optional, intent(in)  :: postname
    ! local variable
    CHARACTER(len=str_medium)                        :: outputname

    success = .false.

    if(.not.init) then
        if(.not.interface_io_init(spec_rank)) return
    end if

    if(present(postname)) then
        select case(trim(output_format))
        case("hdf5")
            if (modulo(layout%nx*layout%ny*layout%nz,getnbcpus()) .ne. 0 ) then 
              write (6,'(a,1x,a)') '[ERROR HDF5] number of grid points is not balanced on each processor for layout:',trim(adjustl(postname))
              return
            endif
            success = hdf5_write_scalar_wrapper(layout, iteration, postname=postname)
        case default ! avs
            write(outputname, avs_format_name) trim(layout%name), "_"//trim(postname)//"_", iteration
            success=dump_UVWP_ForAVS(outputname,spec_rank,layout%Lx,layout%Ly,layout%Lz,layout)
        end select
    else
        select case(trim(output_format))
        case("hdf5")
            if (modulo(layout%nx*layout%ny*layout%nz,getnbcpus()) .ne. 0 ) then 
              write (6,'(a,1x,a)') '[ERROR HDF5] number of grid points is not balanced on each processor for layout:',trim(adjustl(layout%name))
              return
            endif
            success = hdf5_write_scalar_wrapper(layout, iteration, cut=output_2Dslice)
        case default ! avs
            write(outputname, avs_format_name) trim(layout%name), "_", iteration
            success=dump_UVWP_ForAVS(outputname,spec_rank,layout%Lx,layout%Ly,layout%Lz,layout)
        end select
    end if


end function write_datalayout_scalar

!> This subroutine close the input/output context
!! @param[out]  success     = .TRUE. if size is obtain
function interface_io_close() result(success)

    use parser_tools
    use hdf5_io
    use communicators_tools

    logical                 :: success

    success = .false.

    ! If not initialize, then nothing to do.
    if(.not.init) then
        success = .true.
        return
    end if

    ! Some specific intialisation depending on output format
#ifdef HDF5
    if(trim(adjustl(output_format))=='hdf5' .or. trim(adjustl(input_format))=='hdf5') then
       success = hdf5_close()
    else
        success = .true.
    endif
#else
    success = .true.
#endif

    init = .false.

end function interface_io_close

end module io_interface
