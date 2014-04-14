!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!> @addtogroup layout 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: realdatalayout 
!
!> @author
!> Patrick Begou, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to provide the data structure used by codescalar.
!> It describes the blocs organisation on each MPI process. This is a generic version
!> wich should be processed to create a REAL version and a COMPLEX version.
!> DO NOT CHANGE this file: it is automagicaly created from implementdatalayout.Fortran
!> for REAL and COMPLEX data types. All changes must be done in implementdatalayout.Fortran
!> file only.
!------------------------------------------------------------------------------

module realVectordatalayout

    use maindatalayout
    use precision_tools
    use mpilayout_tools
    use parser_tools

    implicit none


    integer, private, parameter ::lennames=str_medium
  
    type REAL_VECTOR_DATA_LAYOUT
        ! -- Generic information --
        !< the name of the datas u, v, scalar 1....etc
        character(len=lennames)                       :: name
        !> The name of each component
        character(len=lennames), dimension(:),pointer :: names_comp => null()
        !> physical dimensions
        real(WP) ::lx,ly,lz
        !> number of components
        integer :: nb_components

        ! -- Grid information --
        integer :: nx         !< global grid size for the datas
        integer :: ny         !< global grid size for the datas
        integer :: nz         !< global grid size for the datas
        ! curlayout should be alongx or alongy or alongz
        integer :: curlayout  !< should be alongx or alongy or alongz
        integer :: xmin  !< start x position in the global array for the current layout
        integer :: xmax  !< final x position in the global array for the current layout
        integer :: ymin  !< start y position in the global array for the current layout
        integer :: ymax  !< final y position in the global array for the current layout
        integer :: zmin  !< start z position in the global array for the current layout
        integer :: zmax  !< final z position in the global array for the current layout
        !> a ncpus x 2 array: ymin and ymax values for each cpu when the layout is along x axis  
        integer, dimension(:,:),pointer :: y4layoutonx => null()
        !> a ncpus x 2 array: zmin and zmax values for each cpu when the layout is along x axis  
        integer, dimension(:,:),pointer :: z4layoutonx => null()
        !> a ncpus x 2 array: xmin and xmax values for each cpu when the layout is along y axis  
        integer, dimension(:,:),pointer :: x4layoutony => null()
        !> a ncpus x 2 array: zmin and zmax values for each cpu when the layout is along y axis  
        integer, dimension(:,:),pointer :: z4layoutony => null()
        !> a ncpus x 2 array: xmin and xmax values for each cpu when the layout is along z axis  
        integer, dimension(:,:),pointer :: x4layoutonz => null()
        !> a ncpus x 2 array: ymin and ymax values for each cpu when the layout is along z axis  
        integer, dimension(:,:),pointer :: y4layoutonz => null()

        ! -- Data --
        !> the datas (u, v....). the array size is (xmin:xmax,ymin:ymax,zmin:zmax,1:nb_components)
        REAL(WP), dimension(:,:,:,:), pointer :: values => null()
    
    end type REAL_VECTOR_DATA_LAYOUT
    
    
    interface assignment(=)
        module procedure real_vector_affectation
    end interface
    
    ! this internal subroutine inherited from maindatalayout is set to private
    ! it should not be seen out of this module.
!    private initlayoutarray
    private real_vector_affectation
  
contains

!------------------------------------------------------------------------------
!> vector data layout initialisation on each process
!! @author 
!! Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> @details
!! this subroutine allocate and organise the data for the local process, taking account
!! for the local distribution along on the the x, y or z axis. default is along x axis.
!
!> @param[in]       name            = the names_comp of this data (u, pressure, scalar_1...). it should not contains 
!!                                      spaces, tabs nor special characters as it can be used to create file names_comps. 
!> @param[in]       names_comp           = the name of each component, optional
!! @param[in,out]   val             = the variable of type vectro_data_layout to be instanciated
!! @param[in,out]   nb_components   = number of vector components
!! @param[in]       nx              = global domain discretization size along x axis
!! @param[in]       ny              = global domain discretization size along y axis
!! @param[in]       nz              = global domain discretization size along z axis
!! @param[in]       lx              = global domain lenght along x axis
!! @param[in]       ly              = global domain lenght along y axis
!! @param[in]       lz              = global domain lenght along z axis
!! @param[in]       nbcpus          = the total number of available cpus 
!! @param[in]       rank            = my rank in the cpus pool. range in  [0,nbcpus-1]
!! @param[in]       how             = the choosen distribution. this parameter is optional and defaults to alongx.
!!                                      acceptable values are alongx, alongy or alongz.
!! @return true if memory is successfuly allocated
!------------------------------------------------------------------------------
logical function real_vector_initdatalayout(name,names_comp,val,nb_components,nx,ny,nz,lx,ly,lz,nbcpus,rank,how)

    ! Input/Ouput
    character(len=*), intent(in)                :: name
    character(len=*), dimension(:), intent(in), optional  :: names_comp
    integer, intent(in)                         :: nx,ny,nz, nb_components
    real(wp), intent(in)                        :: lx,ly,lz
    type(REAL_VECTOR_DATA_LAYOUT),intent(inout) :: val
    integer, intent(in)                         :: nbcpus,rank
    integer, intent(in),optional                :: how   !default to alongx

    ! Local
    logical             :: res
    integer             :: ncpus1,ncpus2,ires
    integer             :: component    ! for loop!
    character(len=lennames)    :: request

    ! -- Init --
    real_vector_initdatalayout=.false.
    if (nb_components<1) then
        write(*,'(a,x,i0,x,a)') '[ERROR] number of component for a vector is', nb_components,'and must be >= 1'
        stop
    end if
    if (nb_components>size(names_comp)) then
        write(*,'(a)') '[ERROR] wrong number of names_comps for vector components'
        stop
    end if
    val%nb_components=nb_components
    
    ! If something is allocated, deallocate it
    if(associated(val%names_comp)) deallocate(val%names_comp)
    if(associated(val%values)) deallocate(val%values)
    if(associated(val%y4layoutonx)) deallocate(val%y4layoutonx)
    if(associated(val%z4layoutonx)) deallocate(val%z4layoutonx)
    if(associated(val%x4layoutony)) deallocate(val%x4layoutony)
    if(associated(val%z4layoutony)) deallocate(val%z4layoutony)
    if(associated(val%x4layoutonz)) deallocate(val%x4layoutonz)
    if(associated(val%y4layoutonz)) deallocate(val%y4layoutonz)

    ! Adjust names_comp
    allocate(val%names_comp(1:nb_components),stat=ires)
    if (ires .ne.0) then
        write(6,'(a,i0,a)')"error on process ",rank,": not enought memory for components names_comp !"
        return
    else
        do component = 1, nb_components
            write(request,'(a,i0)') trim(name), component
            val%names_comp = request
        end do
        if (present(names_comp)) then
            do component = 1, max(nb_components, size(names_comp))
                if (len_trim(names_comp(component)).lt.lennames) then
                    val%names_comp=names_comp(component)
                else
                    val%names_comp=names_comp(component)(1:lennames)
                endif
            end do
        end if
    endif

    val%nx=nx
    val%ny=ny
    val%nz=nz
        
    ! How many cpus in the two decomposed dimensions ?
    call computeprocsmatrix(nbcpus,ncpus1,ncpus2,rank,res)
    if (.not. res) return

    ! Get the data organisation on the cpus
    if (.not. ( &
        & initlayoutarray(val%y4layoutonx,val%z4layoutonx,ncpus1,ncpus2,ny,nz)&
        & .and. &
        & initlayoutarray(val%x4layoutony,val%z4layoutony,ncpus1,ncpus2,nx,nz)&
        & .and. &
        & initlayoutarray(val%x4layoutonz,val%y4layoutonz,ncpus1,ncpus2,nx,ny)&
        )) then
        write(6,'(a,i0,a)')"error on process ",rank,": not enought memory!"
        return
    endif

    ! Get the initial datalayout
    if (present(how)) then
        if (how.eq.alongx .or. how .eq. alongy .or. how .eq. alongz) then
            val%curlayout=how
        else
            write(6,'(a,i0,a,i0)')"error on process ",rank,": unknown layout ",how
            write(6,'(a)')        "      should be in [1,3]"
        endif
    else 
        val%curlayout=alongx  
    endif

    ! Allocate data array, we know curlayout is right so no default case.
    ! we decide do have realistic indexes, not sure it is a good idea!
    select case (val%curlayout)
    case (alongx)
        val%xmin=1
        val%xmax=nx
        val%ymin=val%y4layoutonx(rank+1,1)
        val%ymax=val%y4layoutonx(rank+1,2)
        val%zmin=val%z4layoutonx(rank+1,1)
        val%zmax=val%z4layoutonx(rank+1,2)
    case (alongy)
        val%xmin=val%x4layoutony(rank+1,1)
        val%xmax=val%x4layoutony(rank+1,2)
        val%ymin=1
        val%ymax=ny
        val%zmin=val%z4layoutony(rank+1,1)
        val%zmax=val%z4layoutony(rank+1,2)
    case (alongz)
        val%xmin=val%x4layoutonz(rank+1,1)
        val%xmax=val%x4layoutonz(rank+1,2)
        val%ymin=val%y4layoutonz(rank+1,1)
        val%ymax=val%y4layoutonz(rank+1,2)
        val%zmin=1
        val%zmax=nz
    end select

    val%lx=lx
    val%ly=ly
    val%lz=lz

    allocate(val%values(val%xmin:val%xmax ,val%ymin:val%ymax , val%zmin:val%zmax, 1:nb_components),stat=ires)
    if (ires .ne.0) then
        write(6,'(a,i0,a)')"error on process ",rank,": not enought memory for datas!"
        return
    else
        val%values=0.0_wp
        real_vector_initdatalayout=.true.
    endif

    return

end function real_vector_initdatalayout
    
    
!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> \brief
!> this routine should be called to free memory allocated in a data_layout variable.
!
!
!> @param[in,out] val is the variable of type REAL_VECTOR_DATA_LAYOUT to be freed
!------------------------------------------------------------------------------
subroutine real_vector_deletedatalayout(val)

    type(REAL_VECTOR_DATA_LAYOUT),intent(inout) :: val

    if(associated(val%names_comp))  deallocate(val%names_comp)
    if(associated(val%values))  deallocate(val%values)
    if(associated(val%y4layoutonx))  deallocate(val%y4layoutonx)
    if(associated(val%z4layoutonx))  deallocate(val%z4layoutonx)
    if(associated(val%x4layoutony))  deallocate(val%x4layoutony)  
    if(associated(val%z4layoutony))  deallocate(val%z4layoutony) 
    if(associated(val%x4layoutonz))  deallocate(val%x4layoutonz)
    if(associated(val%y4layoutonz))  deallocate(val%y4layoutonz)
    val%names_comp => null()
    val%values => null()
    val%y4layoutonx => null()
    val%z4layoutonx => null()
    val%x4layoutony => null()
    val%z4layoutony => null()
    val%x4layoutonz => null()
    val%y4layoutonz => null()
    val%xmin=0
    val%xmax=-1
    val%ymin=0
    val%ymax=-1
    val%zmin=0
    val%zmax=-1
    val%nb_components=0

end subroutine real_vector_deletedatalayout    

!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> \brief
!> print the current data layout on this node. for debug purpose only.
!> @param[in] val is the variable of type data_layout for wich the layout is printed
!> @param[in] me my rank in the cpus pool. range in  [0,nbcpus-1]
!------------------------------------------------------------------------------
subroutine real_vector_showlayout(val,me)

    implicit none
    integer, intent(in)::me
    type(REAL_VECTOR_DATA_LAYOUT),intent(in) :: val
    integer ::i
    write(me+10,'(a)')"===================================================================="
    write(me+10,'(a)')'i am a real data layout object called "'//trim(val%name)//'"'
    write(me+10,'(a)')"===================================================================="
    write(me+10,100) "global layout",1,val%nx,1,val%ny,1,val%nz
    if (val%curlayout.eq.alongx) write(me+10,'(a)') "organised along x"
    if (val%curlayout.eq.alongy) write(me+10,'(a)') "organised along y"
    if (val%curlayout.eq.alongz) write(me+10,'(a)') "organised along z"
    write(me+10,100) "local layout",val%xmin,val%xmax,val%ymin,val%ymax,val%zmin,val%zmax
    
    write(me+10,'(a)') "=========== y4layoutonx ==========="
    do i=1,size(val%y4layoutonx,1)
    write(me+10,101) "y, when alongx",i-1,val%y4layoutonx(i,1),val%y4layoutonx(i,2)
    enddo
    write(me+10,'(a)') "=========== z4layoutonx ==========="
    do i=1,size(val%z4layoutonx,1)
    write(me+10,101) "z, when alongx",i-1,val%z4layoutonx(i,1),val%z4layoutonx(i,2)
    enddo
    
    write(me+10,'(a)') "=========== x4layoutony ==========="
    do i=1,size(val%x4layoutony,1)
    write(me+10,101) "x, when alongy",i-1,val%x4layoutony(i,1),val%x4layoutony(i,2)
    enddo
    write(me+10,'(a)') "=========== z4layoutony ==========="
    do i=1,size(val%z4layoutony,1)
    write(me+10,101) "z, when alongy",i-1,val%z4layoutony(i,1),val%z4layoutony(i,2)
    enddo
    
    write(me+10,'(a)') "=========== x4layoutonz ==========="
    do i=1,size(val%x4layoutonz,1)
    write(me+10,101) "x, when alongz",i-1,val%x4layoutonz(i,1),val%x4layoutonz(i,2)
    enddo
    write(me+10,'(a)') "=========== y4layoutonz ==========="
    do i=1,size(val%y4layoutonz,1)
    write(me+10,101) "y, when alongz",i-1,val%y4layoutonz(i,1),val%y4layoutonz(i,2)
    enddo
    write(me+10,'(a)')"============================ end ==================================="
    
     
100 format(a," is [",i3,":",i3,",",i3,":",i3,",",i3,":",i3,"]")   
101 format(a," on cpu ",i0," is [",i3,":",i3,"]")   
 
end subroutine real_vector_showlayout

!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> \brief
!> check that the different REAL_VECTOR_DATA_LAYOUT argument have the same layout :
!> same global dimensions and same local dimensions.
!> @param[in] val is the variable of type REAL_VECTOR_DATA_LAYOUT for wich the layout is checked
!> @param[in] val1 is a variable of type REAL_VECTOR_DATA_LAYOUT to compare with
!> @param[in] val2 is an optional variable of type REAL_VECTOR_DATA_LAYOUT to compare with
!> @param[in] val3 is an optional variable of type REAL_VECTOR_DATA_LAYOUT to compare with
!> @param[in] val4 is an optional variable of type REAL_VECTOR_DATA_LAYOUT to compare with
!> @param[in] val5 is an optional variable of type REAL_VECTOR_DATA_LAYOUT to compare with
!> @return .true. if all the provided variable have the same layout
!------------------------------------------------------------------------------
recursive function real_vector_samelayout(val,val1,val2,val3,val4,val5) result(res)

    implicit none
    
    type(REAL_VECTOR_DATA_LAYOUT),intent(in) :: val, val1
    type(REAL_VECTOR_DATA_LAYOUT),intent(in), optional :: val2
    type(REAL_VECTOR_DATA_LAYOUT),intent(in), optional :: val3
    type(REAL_VECTOR_DATA_LAYOUT),intent(in), optional :: val4
    type(REAL_VECTOR_DATA_LAYOUT),intent(in), optional :: val5
    
    logical :: res
    
    if (present(val5)) then
        res=(real_vector_samelayout(val,val1) &
            &.and. &
        & real_vector_samelayout(val1,val2)&
            &.and. &
        & real_vector_samelayout(val2,val3)&
            &.and. &
        & real_vector_samelayout(val3,val4)&
            &.and. &
        & real_vector_samelayout(val4,val5))
    elseif (present(val4)) then
        res=(real_vector_samelayout(val,val1) &
            &.and. &
        & real_vector_samelayout(val1,val2)&
            &.and. &
        & real_vector_samelayout(val2,val3)&
            &.and. &
        & real_vector_samelayout(val3,val4))
    elseif (present(val3)) then
        res=(real_vector_samelayout(val,val1) &
            &.and. &
        & real_vector_samelayout(val1,val2)&
            &.and. &
        & real_vector_samelayout(val2,val3))
    elseif (present(val2)) then
        res=(real_vector_samelayout(val,val1) &
            &.and. &
            & real_vector_samelayout(val1,val2))
    else
        res= (val%nx .eq. val1%nx &
            & .and. &
         & val%ny .eq. val1%ny &
            & .and. &
         & val%nz .eq. val1%nz &
            & .and. &
         & val%xmin .eq. val1%xmin &
            & .and. &
         & val%xmax .eq. val1%xmax &
            & .and. &
         & val%ymin .eq. val1%ymin &
            & .and. &
         & val%ymax .eq. val1%ymax &
            & .and. &
         & val%zmin .eq. val1%zmin &
            & .and. &
         & val%zmax .eq. val1%zmax &
            & .and. &
         & val%lx .eq. val1%lx &
            & .and. &
         & val%ly .eq. val1%ly &
            & .and. &
         & val%lz .eq. val1%lz &
            & .and. &
         & val%nb_components .eq. val1%nb_components)
    end if
    
    return
    
end function real_vector_samelayout

!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> \brief
!> support for the = operator between two type REAL_VECTOR_DATA_LAYOUT:
!> same global dimensions and same local dimensions.
!> @param [in,out] left is the variable of type REAL_VECTOR_DATA_LAYOUT to set
!> @param [in] right is the variable of type REAL_VECTOR_DATA_LAYOUT to copy in left
!------------------------------------------------------------------------------
    subroutine real_vector_affectation(left, right)
    
        type (REAL_VECTOR_DATA_LAYOUT), intent(inout) :: left
        type (REAL_VECTOR_DATA_LAYOUT), intent(in) :: right

        integer :: ires, dim1, dim2
        ! to be shure that a syntax a=a will work we must use buffers
        !> the datas (u, v....). the array size is (xmin:xmax,ymin:ymax,zmin:zmax)
        REAL(WP), dimension(:,:,:,:), pointer :: values => null()     
        !> a ncpus x 2 array: ymin and ymax values for each cpu when the layout is along x axis  
        integer, dimension(:,:),pointer :: y4layoutonx => null()  
        !> a ncpus x 2 array: zmin and zmax values for each cpu when the layout is along x axis  
        integer, dimension(:,:),pointer :: z4layoutonx => null()    
        !> a ncpus x 2 array: xmin and xmax values for each cpu when the layout is along y axis  
        integer, dimension(:,:),pointer :: x4layoutony => null()    
        !> a ncpus x 2 array: zmin and zmax values for each cpu when the layout is along y axis  
        integer, dimension(:,:),pointer :: z4layoutony => null()    
        !> a ncpus x 2 array: xmin and xmax values for each cpu when the layout is along z axis  
        integer, dimension(:,:),pointer :: x4layoutonz => null()    
        !> a ncpus x 2 array: ymin and ymax values for each cpu when the layout is along z axis  
        integer, dimension(:,:),pointer :: y4layoutonz => null() 
    
        if (.not. associated(right%values)) then
          write(6,'(a,i0,a)')'[error] operator = for type (real_data_layout) : receive uninitialized data!'
          call mpi_finalize(ires)
          stop 'failed'
        endif  
    
        allocate(values(right%xmin:right%xmax, right%ymin:right%ymax, right%zmin:right%zmax,right%nb_components),stat=ires)
        if (ires .ne.0) then
          write(6,'(a,i0,a)')'[error] operator = for type (real_data_layout) : not enought memory!'
          call mpi_finalize(ires)
          stop 'failed'
        endif 
        values=right%values 
        
        dim1=size(right%y4layoutonx,1)
        dim2=size(right%y4layoutonx,2)
        allocate(y4layoutonx(dim1,dim2),z4layoutonx(dim1,dim2), &
               & x4layoutony(dim1,dim2),z4layoutony(dim1,dim2), &
           & x4layoutonz(dim1,dim2),y4layoutonz(dim1,dim2), &
           & stat=ires)
        if (ires .ne.0) then
          write(6,'(a,i0,a)')'[error] operator = for type (real_data_layout) : not enought memory!'
          call mpi_finalize(ires)
          stop 'failed'
        endif  
        y4layoutonx=right%y4layoutonx
        z4layoutonx=right%z4layoutonx
        x4layoutony=right%x4layoutony
        z4layoutony=right%z4layoutony
        x4layoutonz=right%x4layoutonz
        y4layoutonz=right%y4layoutonz
        
        ! now free memory of the left operand
        if (associated(left%names_comp))        deallocate(left%names_comp)
        if (associated(left%values))      deallocate(left%values)
        if (associated(left%y4layoutonx)) deallocate(left%y4layoutonx)
        if (associated(left%z4layoutonx)) deallocate(left%z4layoutonx)
        if (associated(left%x4layoutony)) deallocate(left%x4layoutony)
        if (associated(left%z4layoutony)) deallocate(left%z4layoutony)
        if (associated(left%x4layoutonz)) deallocate(left%x4layoutonz)
        if (associated(left%y4layoutonz)) deallocate(left%y4layoutonz)

        ! rebuild left operand
        left%names_comp           => right%names_comp
        left%nb_components  =  right%nb_components
        left%nx             =  right%nx
        left%ny             =  right%ny
        left%nz             =  right%nz
        left%curlayout      =  right%curlayout
        left%xmin           =  right%xmin
        left%xmax           =  right%xmax
        left%ymin           =  right%ymin
        left%ymax           =  right%ymax
        left%zmin           =  right%zmin
        left%zmax           =  right%zmax
        left%values         => values
        left%y4layoutonx    => y4layoutonx
        left%z4layoutonx    => z4layoutonx
        left%x4layoutony    => x4layoutony
        left%z4layoutony    => z4layoutony
        left%x4layoutonz    => x4layoutonz
        left%y4layoutonz    => y4layoutonz
        left%lx             =  right%lx
        left%ly             =  right%ly
        left%lz             =  right%lz

     end subroutine real_vector_affectation
 
!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou and Jean-Baptiste Lagaert (adaption for Vector from other dalalayout only)
!
!> \brief
!> support for the copy of the structure of type REAL_VECTOR_DATA_LAYOUT without
!> value initialisation (values are set to 0.0_wp). but all the architecture
!> and parallel information is duplicated.
!> @param[in] ori is the variable model of type REAL_VECTOR_DATA_LAYOUT 
!> @param[in,out] copy is the variable of type REAL_VECTOR_DATA_LAYOUT builded
!> @param[in] name_out is the names of the build variable instead of the name of the ori (optional)
!> @param[in] names_out is the names of the build variable instead of the names of the ori (optional)
!------------------------------------------------------------------------------
    function real_vector_copyStructOnly(ori, copy, name_out, names_compout) result(res)
    
        type (REAL_VECTOR_DATA_LAYOUT), intent(inout) :: copy
        type (REAL_VECTOR_DATA_LAYOUT), intent(in) :: ori
        character(len=*),intent(in),optional :: name_out
        character(len=*),dimension(:),intent(in),optional :: names_compout

        logical :: res
        integer :: ires, dim1, dim2, comp

        res=.false.

        if (.not. associated(ori%values)) then
          write(6,'(a,i0,a)')'[error]copystructonly for type (real_data_layout) : receive uninitialized data!'
          return
        endif  
        if (associated(copy%names_comp))  deallocate(copy%names_comp)
        if (associated(copy%values))      deallocate(copy%values)
        if (associated(copy%y4layoutonx)) deallocate(copy%y4layoutonx)
        if (associated(copy%z4layoutonx)) deallocate(copy%z4layoutonx)
        if (associated(copy%x4layoutony)) deallocate(copy%x4layoutony)
        if (associated(copy%z4layoutony)) deallocate(copy%z4layoutony)
        if (associated(copy%x4layoutonz)) deallocate(copy%x4layoutonz)
        if (associated(copy%y4layoutonz)) deallocate(copy%y4layoutonz)

        allocate(copy%values(ori%xmin:ori%xmax, ori%ymin:ori%ymax, ori%zmin:ori%zmax, ori%nb_components),stat=ires)
        if (ires .ne.0) then
            write(6,'(a,i0,a)')'[error] copystructonly for type (real_data_layout) : not enought memory!'
            return
        endif 
        copy%values=0.0_wp 

        copy%nb_components = ori%nb_components

        dim1=size(ori%y4layoutonx,1)
        dim2=size(ori%y4layoutonx,2)
        allocate(copy%y4layoutonx(dim1,dim2),copy%z4layoutonx(dim1,dim2), &
               & copy%x4layoutony(dim1,dim2),copy%z4layoutony(dim1,dim2), &
        & copy%x4layoutonz(dim1,dim2),copy%y4layoutonz(dim1,dim2), &
        & stat=ires)
        if (ires .ne.0) then
            write(6,'(a,i0,a)')'[error] copystructonly for type (real_data_layout) : not enought memory!'
            call mpi_finalize(ires)
            stop 'failed'
        endif  
        copy%y4layoutonx=ori%y4layoutonx
        copy%z4layoutonx=ori%z4layoutonx
        copy%x4layoutony=ori%x4layoutony
        copy%z4layoutony=ori%z4layoutony
        copy%x4layoutonz=ori%x4layoutonz
        copy%y4layoutonz=ori%y4layoutonz

        ! now free memory of the left operand

        ! rebuild copy operand
        if(present(name_out)) then
            copy%name     =trim(name_out)
        else
            copy%name     =ori%name
        endif
        if(present(names_compout)) then
            do comp = 1, max(size(names_compout), ori%nb_components)
                copy%names_comp(comp)     =trim(names_compout(comp))
            end do
        else
            copy%names_comp     =ori%names_comp
        endif
        copy%nx         =ori%nx
        copy%ny         =ori%ny
        copy%nz         =ori%nz
        copy%curlayout  =ori%curlayout
        copy%xmin       =ori%xmin
        copy%xmax       =ori%xmax
        copy%ymin       =ori%ymin
        copy%ymax       =ori%ymax
        copy%zmin       =ori%zmin
        copy%zmax       =ori%zmax
        copy%lx         =ori%lx
        copy%ly         =ori%ly
        copy%lz         =ori%lz
        res=.true.

    end function real_vector_copystructonly


!------------------------------------------------------------------------------
!> Manipulate DIRECTLY the value of one component at the REAL_DATA_LAYOUT format.
!! @author
!! Jean-Baptiste Lagaert
!
!> @details
!!     The returned values is REAL_DATA_LAYOUT whose each of is pointer component
!! are egal to the associated pointer from the vector. This mean that manipulate
!! the return argument mean manipulate the vector too. But you still can make a
!! copy (for instance copy=real_vector_manipulate component) by using the affection
!! "=". If "guru" argument is not defined or set to false, then, it will first delete
!! the datalayout from the output to avoid memoy leak. Be careful, if some of the output
!! component still point to part of the vector, delete the REAL_DATA_LAYOUT mean delete
!! part of the REAL_VECTOR_DATA_LAYOUT too.
!! "Basic" user are supposed to use each time a new output and to copy it.
!! To avoid memory leak, if the input "guru" argument is not given, the 
!! @param[in]       vector      = REAL_VECTOR_DATA_LAYOUT that you want manipulate.
!! @param[in]       comp        = indice of the component you want
!! @param[in,out]   comp_data   = REAL_DATA_LAYOUT whose datas are pointing to the chosen
!!                              component of the vector.
!! @param[in]       guru        = optional argument for advanced use.
!                                   in some case is it important to not deallocate all the pointer of
!!                                  comp_data at the initialisation. But in other general case, doing
!!                                  this deallocation avoid memory leak (and thus is done by default)/
!> @param[in,out] copy is the variable of type REAL_VECTOR_DATA_LAYOUT builded
!> @param[in] name_out is the names of the build variable instead of the name of the ori (optional)
!> @param[in] names_out is the names of the build variable instead of the names of the ori (optional)
!------------------------------------------------------------------------------
    subroutine real_vector_manipulate_comp(vector, comp, comp_data, guru)

    use realdatalayout

        type (REAL_VECTOR_DATA_LAYOUT), intent(in)  :: vector
        type (REAL_DATA_LAYOUT), intent(inout)      :: comp_data
        integer, intent(in)                         :: comp
        logical, intent(in), optional               :: guru

        ! Local
        logical                                     :: clean = .true.

        if (comp>vector%nb_components) then
            write(*,'(a,x,i0,x,a,x,i0)') &
                & '[ERROR] vector_manipulate_comp : indice of wanted component is', &
                & comp, 'and must be <= ', vector%nb_components
            stop
        end if

        if(present(guru)) clean = .not.(guru)
        if(clean) call real_deleteDatalayout(comp_data)

        comp_data%name          = vector%names_comp(comp)

        comp_data%Lx            = vector%Lx
        comp_data%Ly            = vector%Ly
        comp_data%Lz            = vector%Lz

        comp_data%nx            = vector%nx
        comp_data%ny            = vector%ny
        comp_data%nz            = vector%nz

        comp_data%curlayout     = vector%curlayout
        comp_data%xmin          = vector%xmin
        comp_data%xmax          = vector%xmax
        comp_data%ymin          = vector%ymin
        comp_data%ymax          = vector%ymax
        comp_data%zmin          = vector%zmin
        comp_data%zmax          = vector%zmax

        comp_data%Y4layoutOnX   => vector%Y4layoutOnX
        comp_data%Z4layoutOnX   => vector%Z4layoutOnX
        comp_data%X4layoutOnY   => vector%X4layoutOnY
        comp_data%Z4layoutOnY   => vector%Z4layoutOnY
        comp_data%X4layoutOnZ   => vector%X4layoutOnZ
        comp_data%Y4layoutOnZ   => vector%Y4layoutOnZ

        comp_data%values        => vector%values(:,:,:,comp)

    end subroutine real_vector_manipulate_comp

end module realVectordatalayout
!> @}
