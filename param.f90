!> @addtogroup simu_setup
!! @{
!------------------------------------------------------------------------------
!
! MODULE: param
!
!
! DESCRIPTION:
!> The module ``param'' gathers all the simulation parameters and setup.
!! @details
!!     This module gathers information about simulation setup and parameters
!! (like number of scalar, presence of magnetic field, ...) but also
!! about I/O and post-process.
!!     More precisly, these values/parameters have (eventually) default value
!! defines in this module and are update by using values defined in the
!! "input" file.
!!     All variables defined in this module are protected : they must be
!! visible by any other modules but must not be changed.
!!
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module param

    use precision_tools
    !use parser_tools, only : parser_read, parser_is_defined
    use parser_tools

    implicit none

    ! ##########################################
    ! ##########      Procedures      ##########
    ! ##########################################
    !
    ! ===== Public procedures =====
    public param_countScalars



contains

!> Initialize the simulation context (setup physical and numerical parameters)
!> @author  Jean-Baptiste Lagaert, LEGI
!! @details
!! Initialize the parser reseting all allocated datas.
!! Set to private because should only be called when reading a new file
!! by parser_parsefile
!!
!! @param[in]   input       = The name of the input file, optional.
!! @return      success     = .TRUE. if no error, else retruns .FALSE.
function param_init(input) result(success)

    use parser_tools

    logical                                 :: success
    character(len=*), intent(in), optional  :: input

    ! Init context
    success = .false.
    call parser_reset()

    ! Set default value
    if(.not. param_default()) then 
        write(6,'(a)')"Param : unable to init setup to default values"
        return
    end if

    ! Then update and complete the setup with the input file
    if(present(input)) then
        if(.not. parser_parsefile(input)) then 
            write(6,'(a)')"Parser : unable to update default setup from input file"
            return
        end if
    end if

    ! And check correctness.
    if(.not. param_check()) then 
        write(6,'(a)')"Param: unable to check if nothing is missing"
        return
    end if

    success = .true.

end function param_init

!> Set basic parameters to default values
!! @return  success = logical exit of the function (true if success, false if something go wrong)
function param_default() result(success)

    use parser_tools

    ! Input/Output variables
    logical     :: success

    ! Local variables

    ! ===== Initialisation =====
    success = .false.

    ! ===== General setup =====
    ! -- Physical parameters --
    if (.NOT. parser_newentry('Length of domain', '6.283185')) return
    if (.NOT. parser_newentry('Simulation name', 'Isotropic turbulence')) return
    ! -- Numerical setup --
    if (.NOT. parser_newentry('Default topology', '0')) return
    if (.NOT. parser_newentry('Simulation iterations', '1')) return
    if (.NOT. parser_newentry('Fixed time step', "-1.0")) return
    if (.NOT. parser_newentry('Maximum time step', "1000000.0")) return
    if (.NOT. parser_newentry('Convective time step number', "0.5")) return
    if (.NOT. parser_newentry('Viscous time step number', "1.0")) return
    if (.NOT. parser_newentry('Integration Method', "0")) return
! Integration Method
    ! ===== Parameters about the flow =====
    ! -- Physical parameters about the flow -----
    if (.NOT. parser_newentry('Viscosity', '1E-02')) return
    if (.NOT. parser_newentry('Forcing Method',"none")) return
    if (.NOT. parser_newentry('Power factor',"1.0")) return
    ! -- Numerical parameters about the flow --
    if (.NOT. parser_newentry('Number of points', '64')) return
    if (.NOT. parser_newentry('Dealiasing method', "none") ) return
    if (.NOT. parser_newentry('LES model for velocity', "none")) return
    if (.NOT. parser_newentry('Type of averaging for velocity', '0')) return
    if (.NOT. parser_newentry('Nullify first wave number', "no")) return
    if (.NOT. parser_newentry('Get pdf of sgs term for qdm',".false.")) return
    ! -- forcing parameters
    if (.NOT. parser_newentry('Forcing area', "1.0")) return
    if (.NOT. parser_newentry('Additional mode', "0")) return


    ! ===== Scalar parameters for scalars solved with pseudo-spectral methods =====
    if (.NOT. parser_newentry('Forcing scalar_spec', '0')) return

    ! ===== Scalar parameters for scalars solved with mixed particle/spectral methods =====
    if (.NOT. parser_newentry('Number of scalar for particle solver', '0')) return
    if (.NOT. parser_newentry('Constant for particle time step', '0.5')) return
    if (.NOT. parser_newentry('group_size', '5')) return
    if (.NOT. parser_newentry('Forcing scalar_part', '0')) return
    if (.NOT. parser_newentry('advection method', 'p_M6')) return
    if (.NOT. parser_newentry('advection interpol', 'none')) return

    ! ===== MHD =====
    if (.NOT. parser_newentry('Magnetic Prandl number', "0E+00")) return


    ! ===== Gravity =====
    if (.NOT. parser_newentry('X-gravity', "0.0")) return
    if (.NOT. parser_newentry('Y-gravity', "0.0")) return
    if (.NOT. parser_newentry('Z-gravity', "0.0")) return

    ! ===== Post-process =====
    if (.NOT. parser_newentry('Post-process frequency', "0")) return
    if (.NOT. parser_newentry('Show divergence', "0")) return
    if (.NOT. parser_newentry('Show memory', "0")) return
    if (.NOT. parser_newentry('Show basic info', "0")) return
    if (.NOT. parser_newentry('Show Magnetic growth rate', ".false.")) return
    if (.NOT. parser_newentry('Scalar - basic stats', "0")) return
    if (.NOT. parser_newentry('Compute scalar PDF','0')) return
    if (.NOT. parser_newentry('Compute spectrum - ite', "0")) return
    if (.NOT. parser_newentry('Compute spectrum - time', "0.0")) return
    if (.NOT. parser_newentry('Type of output', "avs")) return
    if (.NOT. parser_newentry('slice or 3D', "3D")) return
    if (.NOT. parser_newentry('Post-process out simulation','none')) return
    if (.NOT. parser_newentry('Save frequency','0')) return
    if (.NOT. parser_newentry('Write output - time','-1.0')) return
    if (.NOT. parser_newentry('Plane jet AVG - time','-1.0')) return
    if (.NOT. parser_newentry('Plane jet AVG - ite','-1')) return
    if (.NOT. parser_newentry('diff diff - ite','-1')) return
    if (.NOT. parser_newentry('diff diff - time','-1.0')) return
    if (.NOT. parser_newentry('diff diff scal1', '1')) return
    if (.NOT. parser_newentry('diff diff scal2', '1')) return
    if (.NOT. parser_newentry('mixing threshold', '0.02')) return
    if (.NOT. parser_newentry('Post Kolmo', '0')) return
    if (.NOT. parser_newentry('Bxext', '0E+00')) return
    if (.NOT. parser_newentry('Byext', '0E+00')) return
    if (.NOT. parser_newentry('Bzext', '0E+00')) return
    if (.NOT. parser_newentry('Number of points for rediscretization','0')) return
    if (.NOT. parser_newentry('a priori test in simu','0')) return
    if (.NOT. parser_newentry('filter size min for a priori test','8')) return
    if (.NOT. parser_newentry('filter size max for a priori test','8')) return
    if (.NOT. parser_newentry('type of filter for a priori test','Cutoff')) return
    if (.NOT. parser_newentry('Size of filter min','8'))return
    if (.NOT. parser_newentry('Size of filter max','8'))return
    if (.NOT. parser_newentry('Size of filter','8'))return
    if (.NOT. parser_newentry('MHD post process', "0")) return
    if (.NOT. parser_newentry('HD post process', "0")) return
    if (.NOT. parser_newentry('Scalar - moment stats', "0")) return
    if (.NOT. parser_newentry('anisotropic', "no")) return
    if (.NOT. parser_newentry('a priori test in simu for velocity', "0")) return
    if (.NOT. parser_newentry('No Energy', "no")) return

    success = .true.

end function param_default


!> Check correctness of setup and adjust it if needed.
!! @return  success = logical exit of the function (true if success, false if something go wrong)
function param_check() result(success)

    use parser_tools

    logical                                 :: success

    integer                     :: number_scal  ! number of scalar
    integer                     :: i            ! lood indice
    integer                     :: nx, ny,nz
    character(len=str_medium)   :: tag, tag_bis ! tag of the entry inside the input
    character(len=str_medium)   :: tag_value    ! value associated to a tag

    ! ===== Initialisation =====
    success = .false.

    ! ===== General context =====
    call parser_read('Simulation name', tag_value)
    if (trim(tag_value)=='Taylor Green') then
        if (.not. parser_is_defined('TG dimensions')) then
            if (.NOT. parser_newentry(tag,'1')) return
        end if
    else if (trim(tag_value)=='Isotropic turbulence') then
        tag = 'disc initialisation'
        if (.NOT. parser_is_defined(tag))  then
            if (.NOT. parser_newentry(tag,'0')) return
        end if
    else if (trim(tag_value)=='HIT_gotoh') then
        tag = 'k0'
        if (.NOT. parser_is_defined(tag))  then
          if (.NOT. parser_newentry(tag,'5')) return
        end if
        tag = 'kv'
        if (.NOT. parser_is_defined(tag))  then
          if (.NOT. parser_newentry(tag,'6')) return
        end if
        tag = 'ktheta'
        if (.NOT. parser_is_defined(tag))  then
          if (.NOT. parser_newentry(tag,'6')) return
        end if
        tag = 'theta_prime'
        if (.NOT. parser_is_defined(tag))  then
          if (.NOT. parser_newentry(tag,'1')) return
        end if
        tag = 'uprime'
        if (.NOT. parser_is_defined(tag))  then
          if (.NOT. parser_newentry(tag,'1')) return
        end if
    else if (trim(tag_value)=='Plane Jet') then
        tag = 'noise wave number ke'
        if (.NOT. parser_is_defined(tag))  then
            if (.not.GetNxNyNz('Number of points',nx,ny,nz)) return
            write(tag_value,'(f10.3)') (4./3.)*min(nx,ny,nz)/2.
            if (.NOT. parser_newentry(tag,tag_value)) return
        end if
    end if
    ! -- Group size for particle method --
    tag = 'advanced_group_size'
    tag_value = '.false.'
    if(parser_is_defined(tag)) then
        call parser_getsize(tag,i)
        if (i==3) then
            tag_value = '.true.'
        end if
    end if
    tag = 'advanced_gp_size'
    if (.NOT. parser_newentry(tag,tag_value)) return


    ! ===== Flow =====
    call parser_read('Forcing Method',tag_value)
    if (trim(tag_value)=='brandenburg') then
        if (.NOT. (&
        & parser_newentry('Helicoidal',"yes") .AND. &
        & parser_newentry('Peak forcing',"5.0_WP") .AND. &
        & parser_newentry('Forcing amplitude',"0.20_WP"))&
        & ) return
    end if
    if (trim(tag_value)=='alvelius') then
        if (.NOT. parser_is_defined('Helicoidal'))  then
            if (.NOT. parser_newentry('Helicoidal',"no")) return
        end if

        if (.NOT. parser_is_defined('No Energy'))  then
            if (.NOT. parser_newentry('No Energy',"no")) return
        end if

        if (.NOT. parser_is_defined('anisotropic'))  then
            if (.NOT. parser_newentry('anisotropic',"no")) return
        end if

    end if


    ! ===== Scalar simulated with pseudo-spectral methods =====
    number_scal=param_countScalars()
    if (number_scal>=1) then
        ! Residence term
        do i =1, number_scal
            write(tag,'(i0,a)') i, ' residence time for'
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'-1'))return
            end if
        end do
        ! Adjust forcing term
        call parser_read('Simulation name', tag_value)
        if (trim(tag_value)/='Isotropic turbulence') then
            do i =1, number_scal
                write(tag,'(a,i0)') 'Forcing for scalar_spec ',i
                if (.NOT.  parser_newentry(tag,'0')) return
            end do
        else
            do i =1, number_scal
                write(tag,'(a,i0)') 'Forcing for scalar_spec ',i
                if (.NOT. parser_is_defined(tag)) then
                    if (.NOT.  parser_newentry(tag,'0'))return
                end if
            end do
        end if
        do i =1, number_scal
            write(tag,'(a,i0)') 'Amplitude of spectral forcing for scalar_spec ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
        do i =1, number_scal
            write(tag,'(a,i0)') 'kmin of spectral forcing for scalar_spec ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
        do i =1, number_scal
            write(tag,'(a,i0)') 'kmax of spectral forcing for scalar_spec ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
        ! Adjust LES if needed
        write(tag,'(a)') 'Filter for velocity model'
        if (.NOT. parser_is_defined(tag)) then
            if (.NOT.  parser_newentry(tag,'0'))return
        end if
        do i =1, number_scal
            write(tag,'(a,i0)') 'LES model for scalar ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
            write(tag,'(a,i0)') 'Filter for scalar model ', i
            if (.NOT. parser_is_defined(tag))then
                if (.NOT.  parser_newentry(tag,'Cutoff'))return
            end if 
        end do
        ! Adjust Schmidt number if needed
        do i =1, number_scal
            write(tag,'(a,i0)') 'Schmidt number for scalar ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT. parser_newentry(tag,'0.7'))return
            end if 
        end do
        ! Adjust CoefRGM model (and coef in general)
        do i =1, number_scal
            write(tag,'(a,i0)') 'Model coefficient for scalar ',i
            if (.NOT. parser_is_defined(tag)) then
                write(tag_bis,'(a,i0)') 'LES model for scalar ',i
                call parser_read(tag_bis,tag_value)
                if ( trim(adjustl(tag_value)) .eq. 'scalesimilarity') then
                  if (.NOT. parser_newentry(tag,'0.2'))return
                else
                  if (.NOT. parser_newentry(tag,'0.08333334'))return
                endif 
            end if 
        end do
        ! Adjust type of averaging for dynamic procedure
        do i =1, number_scal
            write(tag,'(a,i0)') 'Type of averaging for Scalar ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT. parser_newentry(tag,'0'))return
            end if 
        end do
        ! Adjust dealiasing if needed
        do i =1, number_scal
            write(tag,'(a,i0)') 'Dealiasing method for scalar ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT. parser_newentry(tag,'none'))return
            end if 
        end do
    end if

    ! -- Scalar simulated with particle/spectral methods --
    call parser_read('Number of scalar for particle solver', number_scal)
    if (number_scal>=1) then
        ! Adjust mpi topology
        if (.NOT.   parser_newentry('Default topology', '3'))return
        ! Defined remesh formula if needed
        tag = 'advection method'
        if (.NOT. parser_is_defined(tag))then
            if (.NOT. parser_newentry(tag, 'p_O2'))return   
    end if 
        ! Adjust forcing term
        call parser_read('Simulation name', tag_value)
        if (trim(tag_value)/='Isotropic turbulence') then
            do i =1, number_scal
                write(tag,'(a,i0)') 'Forcing for scalar_part ',i
                if (.NOT. parser_newentry(tag,'0'))return
            end do
        else
            do i =1, number_scal
                write(tag,'(a,i0)') 'Forcing for scalar_part ',i
                if (.NOT. parser_is_defined(tag)) then
                    if (.NOT. parser_newentry(tag,'0'))return
                end if 
            end do
        end if
        ! Adjust Schmidt number if needed
        do i =1, number_scal
            write(tag,'(a,i0)') 'Schmidt number for scalar_part ',i
            if (.NOT. parser_is_defined(tag))then
                if (.NOT. parser_newentry(tag,'0.7'))return
            end if 
        end do
        do i =1, number_scal
            write(tag,'(a,i0)') 'Amplitude of spectral forcing for scalar_part ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
        do i =1, number_scal
            write(tag,'(a,i0)') 'kmin of spectral forcing for scalar_part ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
        do i =1, number_scal
            write(tag,'(a,i0)') 'kmax of spectral forcing for scalar_part ',i
            if (.NOT. parser_is_defined(tag)) then
                if (.NOT.  parser_newentry(tag,'0'))return
            end if
        end do
    end if

    ! ===== MHD =====
    call parser_read('Magnetic Prandl number', tag_value)
    if (trim(tag_value)=="0E+00_WP") then
        tag = 'Lorentz time step number'
        if (.NOT. parser_is_defined(tag)) then
            if (.NOT. parser_newentry(tag,'0.5'))return
        end if 
        tag = 'MHD filter size'
        if (.NOT. parser_is_defined(tag)) then
            if (.NOT. parser_newentry(tag,'4'))return
        end if 
        tag = 'MHD filter type'
        if (.NOT. parser_is_defined(tag)) then
        if (.NOT. parser_newentry(tag,'4'))return
    end if 
    end if

    ! ===== Post-process and others =====
    ! Output
    tag = 'Write frequency'
    call parser_read('Simulation iterations',tag_value)
    if (.NOT. parser_is_defined(tag))then
        if (.NOT. parser_newentry(tag,tag_value))return
    end if 
    ! Post in simu must be equal to the min of post-process traitment

    ! Post out simu
    call parser_read('Post-process out simulation',tag_value)
    if (trim(tag_value)=="yes") then
        tag = 'Post-processing'
        if (.NOT. parser_is_defined(tag)) then
            if (.NOT. parser_newentry(tag,"none"))return
        end if 
        call parser_read(tag,tag_value)
        if ((trim(tag_value)=='smagClarkOE').or. (trim(tag_value)=="diffSub")) then
            tag_bis = 'Size of filter'
            if (.NOT. parser_is_defined(tag))then
                if (.NOT. parser_newentry(tag,"1"))return
            end if 
        end if
        tag = 'Number of pool of file'
        if (.NOT. parser_is_defined(tag)) then
            if (.NOT. parser_newentry(tag,"1"))return
        end if 
    end if




    success = .true.

end function param_check

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This subroutine returns the number of scalar requested in the input file
!> We expect zero or more scalar definition.
!------------------------------------------------------------------------------

INTEGER FUNCTION param_countScalars()

  USE parser_tools
  
  IMPLICIT NONE
  CHARACTER(len=64):: request

  param_countScalars=0

  !first one with backward syntax compatibility
  if (parser_is_defined('Number of points for scalar') .OR. &
    & parser_is_defined('Number of points for scalar 1'))   then
        param_countScalars=1
  else
        return
  end if
  
  
  !other occurences
  DO
     WRITE(request,'(a,i0)') 'Number of points for scalar ',param_countScalars+1
     if (parser_is_defined(request)) then
        param_countScalars=param_countScalars+1
     else
        EXIT
     end if     
  ENDDO
    
  return
END FUNCTION param_countScalars

end module param
!! @}
