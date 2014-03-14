!USEFORTEST io
!USEFORTEST topo
!> @addtogroup output
!! @{

!------------------------------------------------------------------------------
!
! MODULE: parallel_out
!
!
! DESCRIPTION: 
!> This module provide all procedure needed to perform parallel and distribued
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
!!    dimension, number of processus, ...) is initiliazed by calling "vtkxml_plus_init_all"
!!    the context specific to each field (mesh resolution, number of point, name of
!!    the ouput, information about time sequence, ...) is initialized or save by
!!    calling "vtkxml_plus_init_field". After that, call "parallel_write" in order to
!!    create a new output of the field.
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module vtkxml_plus

    use precision_tools
    use vtkxml

    implicit none


    !> To write files
    interface vtkxml_plus_write
        !module procedure parallel_vect3D, parallel_scalar
        module procedure vtkxml_scalar, vtkxml_plus_scalar_realdatalayout
    end interface vtkxml_plus_write

    !> To init output context for a field
    interface vtkxml_plus_init_field
        module procedure vtkxml_init_field_basic, vtkxml_init_field_iodata, &
            & vtkxml_plus_init_field_realdatalayout
    end interface vtkxml_plus_init_field

    ! ===== Public procedures =====
    !> Generique procedure to write field in distribued vtk files
    public                  :: vtkxml_plus_write
    ! Initialisation of the mesh context
    public                  :: vtkxml_plus_init_field

    ! ===== Private procedures =====
    ! Specifique procedure to init output context for a specific field.
    private                  :: vtkxml_plus_init_field_realdatalayout
    ! Specifique procedure to write field in distribued vtk files
    private                 :: vtkxml_plus_scalar_realdatalayout

    contains


!> Initialize the context associated to a given field stored a real_data_layout variable.
!!    @param[out]   tag         = tag used to identifiate the field from the others one when using "vtkxml_write"
!!    @param[in]    realdata    = real_data_layout variable from which output will be generated.
!! @details
!!    This subroutine do the same job than the basic version but extract all information from the real_data_layout.
!!    Therefore, it is able to deal with subdomain of different size.
subroutine vtkxml_plus_init_field_realdatalayout(realdata, tag)

    use mpi
    use realdatalayout
    use interface_layout_mesh

    ! Input/Output
    type(real_data_layout), intent(in)          :: realdata
    integer, intent(out)                        :: tag
    ! Local variables
    integer                                     :: direction
    type(io_data)                               :: layout_data


    ! Save io_data for the field
    layout_data%f_name = realdata%name
    layout_data%iteration  = 0
    layout_data%iteration_slice = 0
    ! Generate a the right "cartesian_mesh" variable
    call mesh_init(layout_data%mesh, realdata)

    ! Compute piece extend
    allocate(layout_data%piece_extend(0:total_nb_proc-1,3,2))
    if (realdata%curlayout == alongX) then
        layout_data%piece_extend(:,1,1) = realdata%xmin
        layout_data%piece_extend(:,1,2) = realdata%xmax
        layout_data%piece_extend(:,2,:)=realdata%Y4layoutOnX(:,:)
        layout_data%piece_extend(:,3,:)=realdata%Z4layoutOnX(:,:)
        do direction = 2,3
            where (layout_data%piece_extend(:,direction,2) /= &
                    & layout_data%mesh%N(direction))
                layout_data%piece_extend(:,direction,2) = layout_data%piece_extend(:,direction,2) + 1
            end where
        end do
    else
        stop '[error] real_data_layout is not alongX - vtk xml output not yet implemented.'
    end if

    call vtkxml_init_field(layout_data, tag)

end subroutine vtkxml_plus_init_field_realdatalayout


! ===========================================================
! ==================== Private procedures ===================
! ===========================================================

!> Write an output at "vtk xml" format for an "one-component" field stored in a "real_data_layout" variables.
!!    @param[in,out]    tag     = tag associated to the field to save (ie matching indice in the table "field_data")
!!    @param[in]        field   = field at format real_data_layout
subroutine vtkxml_plus_scalar_realdatalayout(tag, field)

    use realdatalayout

    integer, intent(inout)                              :: tag
    type(real_data_layout), intent(in)                  :: field

    call vtkxml_scalar(tag, field%values, field%name)

end subroutine vtkxml_plus_scalar_realdatalayout


!> Write an output at "vtk xml" format for a slice (along X) of an "one-component" field stored in a "real_data_layout" variables.
!!    @param[in,out]    tag     = tag associated to the field to save (ie matching indice in the table "field_data")
!!    @param[in]        field   = field at format real_data_layout
subroutine vtkxml_plus_scalar_realdatalayout_slice(tag, field)

    use realdatalayout

    integer, intent(inout)                              :: tag
    type(real_data_layout), intent(in)                  :: field

    call vtkxml_scalar_slice(tag, field%values, field%name)

end subroutine vtkxml_plus_scalar_realdatalayout_slice

end module vtkxml_plus
!> @}
