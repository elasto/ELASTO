!------------------------------------------------------------------------------
!
! MODULE: data
!
!> @author
!! Patrick BEGOU and Jean-Baptiste LAGAERT, LEGI
!
! DESCRIPTION:
!> Ths module contains the basic field (velocity, scalar, mhd and associated field)
!! and the interface to allocate or dellocate all of them. Scalar may have different
!! resolutions
!------------------------------------------------------------------------------
MODULE data

    use datalayout
    use wavenumber_tools

    implicit none


! ###########################################
! #####     Functions and subroutines   #####
! ###########################################
    public data_init_luca
    public data_delete
    public data_read
    private dataReadSingleField

! ###########################################
! #####             Variables           #####
! ###########################################

    ! == Spectral velocity  ==
    !> Storage for the Spectral values of U
    type(COMPLEX_DATA_LAYOUT), public, save :: Uk
    !> Storage for the Spectral values of V
    type(COMPLEX_DATA_LAYOUT), public, save :: Vk
    !> Storage for the Spectral values of W
    type(COMPLEX_DATA_LAYOUT), public, save :: Wk
    ! == Other spectral fields ==
    !> Storage for the spectral values of each scalar
    type(COMPLEX_DATA_LAYOUT),allocatable, dimension(:), public :: ScalArrayK
    !> Storage for the spectral values of each scalar - scalar advected with particles method
    type(COMPLEX_DATA_LAYOUT),allocatable, dimension(:), public :: Scal_partArrayK
    !> Storage for the Spectral values of B
    type(COMPLEX_DATA_LAYOUT), allocatable, dimension(:), public :: Bk

    ! == Physical properties ==
    !> Physical size of the domain
    real(WP), protected                     :: Lx,Ly,Lz
    !> Storage for Schmidt numbers for scalars
    real(WP), dimension(:), pointer, public:: schmidt
    !> Storage for Schmidt numbers for scalars - scalar advected with particles method
    real(WP), dimension(:), pointer, public:: schmidt_part


    ! == Waves number ==
    !> Waves numbers for the velocities fields
    type(WaveNumbers), public, save::VelWN
    !> Waves number for the scalars fields
    type(WaveNumbers), dimension(:), pointer, public, save::ScalWN => NULL()
    !> Waves number for the scalars fields - scalar advected with particles method
    type(WaveNumbers), dimension(:), pointer, public, save::Scal_partWN => NULL()
    !> Waves number for the B-fields
    type(WaveNumbers), public, save ::BfieldWN

    !> Work arrays for compute non linear term
    type(COMPLEX_DATA_LAYOUT), public, save                     :: nl_x, nl_y, nl_z
    type(COMPLEX_DATA_LAYOUT), allocatable, dimension(:), public:: nl_Scalarray
    type(COMPLEX_DATA_LAYOUT), allocatable, dimension(:), public:: nl_B


    ! ===== Storage for specific quantity usefull only for specific setup =====
    !> Storage of setup about computing (or not) residence time
    integer, dimension(:), allocatable, public, save    :: res_t
    !> Mixing threshold for residence time
    real(WP), protected                                 :: mix_threshold

contains


!> @author
!> Luca Marradi
!
!>
!> @details
!> This subroutine initialise the datas storage for sigma
!
!> @param[in]     ncpus             = number of MPI processes to distribute data
!> @param[in]     spec_rank         = my MPI processe rank in the range
![0:ncpus-1]
!> @param[in,out] sigma                 = longitudinal velocity storage to
!                                     initialize
!> @return .TRUE. if initialization (data allocation) is successfull
!------------------------------------------------------------------------------
!=========================================================================================!
LOGICAL FUNCTION data_init_luca(ncpus,spec_rank,mult,sigma,sigma_old,Esigma,Lsigma,G,sigmak,n_activ,n_old,convol)
!=========================================================================================!
    USE datalayout
    USE param
    USE cart_topology

    IMPLICIT NONE

    !=== INPUT/OUTPUT DATA === 
    INTEGER, INTENT(IN)                         :: ncpus,spec_rank
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)        :: sigma,sigma_old,mult,n_old
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)        :: Esigma,Lsigma,n_activ,convol
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT)     :: G,sigmak 

    data_init_luca=.FALSE.

    IF (.NOT. GetLxLyLz("Length of domain",Lx,Ly,Lz)) THEN
        IF(spec_rank .EQ.0) &
        & WRITE(6,'(a)')"[ERROR] data_init_luca cannot get the dimensions of the domain "
        RETURN

    ENDIF
 
    !=== INITIALIZATION OF SIGMA SIGMAOLD GAMMADOT AND LSIGMA ===
    !-------------------------------
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "mult", 'Number of points',mult)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "sigma", 'Number of points',sigma)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "sigma_old", 'Number of points',sigma_old)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "Lsigma", 'Number of points',Lsigma)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "Esigma", 'Number of points',Esigma)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "n_activ", 'Number of points',n_activ)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "n_old", 'Number of points',n_old)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "convol", 'Number of points',convol)) RETURN

    !=== ALLOCATE SPACE FOR THE PROPAGATOR IN THE SPECTRAL SPACE ===
    IF(.NOT.initDataLayout("G",G,(sigma%nx/2)+1,sigma%ny,sigma%nz,sigma%Lx,&
           & sigma%Ly,sigma%Lz,ncpus,spec_rank,alongZ)) THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver on process ',spec_rank,&
           & ': not enought memory for the propagator G in the spectral space.'
        RETURN
    ENDIF

    !=== ALLOCATE SPACE FOR SIGMAK IN THE SPECTRAL SPACE ===
    IF (.NOT. initDataLayout("sigmak",sigmak,(sigma%nx/2)+1,sigma%ny,sigma%nz,& 
           & sigma%Lx,sigma%Ly,sigma%Lz,ncpus,spec_rank,alongZ))      THEN
        WRITE(6,'(a,i0,a)')'[ERROR] initSolver_luca on process ',spec_rank,&
           & ': not enought memory for sigmak velocities in fourier space.'
        RETURN
    ENDIF


    data_init_luca=.TRUE.
    
END FUNCTION data_init_luca


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This subroutine initialise the datas storage for U,V,W,P and
!> a variable number of scalars. Scalars solved with a full pseudo-spectral
!> method are stored in a dynamicaly allocated array ScalArray and can have
!> different resolution and scalar solved with mixed pseudo-spectral and
!> particle methods are stored in a dynamicaly allocated array called
!> Scal_partArray and shared all the same mesh resolution.
!
!> @param[in]     ncpus             = number of MPI processes to distribute data
!> @param[in]     spec_rank         = my MPI processe rank in the range [0:ncpus-1]
!> @param[in,out] U                 = longitudinal velocity storage to initialize
!> @param[in,out] V                 = vertical velocity  storage to initialize
!> @param[in,out] W                 = spanwise velocity  storage to initialize
!> @param[in,out] B                 = magnetic field (3D vector)
!> @param[in,out] P                 = Pressure storage to initialize (but is it usefull in the code ?)
!> @param[in,out] ScalArray         = a pointer to an array of scalars storage to initialize.
!> @param[in,out] Scal_partArray    = a pointer to an array of scalars (solved with particles solver for advection term) storage to initialize.
!> @param[out]    nbscal            = the number of scalar fields initialized.
!> @param[in]     nbscal_part       = the number of scalar fields initialized - for scalar advected with particles method
!> @param[in]     imhd              = logical set to true if mhd is active
!> The number of scalar created
!> depends of the information found in the input file
!> @return .TRUE. if initialization (data allocation) is successfull
!------------------------------------------------------------------------------
!=========================================================================================!
LOGICAL FUNCTION data_init(ncpus, spec_rank,U,V,W,B,ScalArray, Scal_partArray,nbscal, nbscal_part,imhd)
!=========================================================================================!
    USE datalayout
    USE param
    USE cart_topology

    IMPLICIT NONE

    ! Input/Output
    INTEGER, INTENT(IN)                                         :: ncpus, spec_rank
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT)                        :: U,V,W
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT),POINTER, DIMENSION(:)  :: ScalArray, Scal_partArray
    TYPE(REAL_DATA_LAYOUT),INTENT(INOUT),POINTER, DIMENSION(:)  :: B
!    TYPE(REAL_VECTOR_DATA_LAYOUT),INTENT(INOUT)                 :: Scal_partVector
    INTEGER, INTENT(OUT)                                        :: nbscal
    INTEGER, INTENT(IN)                                         :: nbscal_part
    LOGICAL, INTENT(IN)                                         :: imhd
    ! Local variables
    INTEGER                                 :: ires, size
    INTEGER                                 :: cpt, start
    CHARACTER(len=64)                       :: request
    CHARACTER(LEN=str_medium), DIMENSION(:), ALLOCATABLE :: name_of_field   ! to store field name

    data_init=.FALSE.

    IF (.NOT. GetLxLyLz("Length of domain",Lx,Ly,Lz)) THEN
        IF(spec_rank .EQ.0) &
        & WRITE(6,'(a)')"[ERROR] data_init cannot get the dimensions of the domain "
        RETURN
    ENDIF

    nbscal=param_countScalars()
    size = max(nbscal, nbscal_part, 3)
    ALLOCATE(name_of_field(size),stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,2(i0,a))')"[ERROR] data_init on process ",&
	&spec_rank,": not enought memory for name_of_field(",&
	&size,")!"
        RETURN
    ENDIF

    !-------------------------------
    !Initialization of U,V,W and P
    !-------------------------------
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "U", 'Number of points',U)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "V", 'Number of points',V)) RETURN
    IF (.not. data_init_Rsingle(ncpus, spec_rank, "W", 'Number of points',W)) RETURN

    !-------------------------------
    ! initialization of the array of scalars solved with only pseudo spectral
    ! method
    !-------------------------------
    ! How many scalars ?
    IF (ASSOCIATED(ScalArray)) DEALLOCATE(ScalArray)
    ScalArray=>null()
    IF (nbscal>0) THEN
        ALLOCATE(ScalArray(nbscal),stat=ires)
        IF (ires .NE.0) THEN
            WRITE(6,'(a,i0,a)')"[ERROR] on process ",spec_rank,": not enought memory!"
            RETURN
        ENDIF

        ! get the size and allocate
        start=1
        if (parser_is_defined('Number of points for scalar')) then
            !allocate first one (backward compatibility without numbering the scalar
            if (.not. data_init_Rsingle(ncpus, spec_rank, "ScalarS1", 'Number of points for scalar',ScalArray(1))) return
            start=2
        endif

        do cpt=start,nbscal
            write(request,'(a,i0)') 'Number of points for scalar ',cpt
            write(name_of_field(1),'(a,i0)') 'ScalarS',cpt
            if (.not. data_init_Rsingle(ncpus, spec_rank, name_of_field(1), request, ScalArray(cpt))) return
        end do
    end if

    !-------------------------------
    ! initialization of the array of scalars solved with
    ! particle method for advection and pseudo-spectral
    ! diffusion solver
    !-------------------------------
    ! How many scalars ?
    IF (ASSOCIATED(Scal_partArray)) DEALLOCATE(Scal_partArray)
    Scal_partArray=>null()
    IF (nbscal_part.GT.0) THEN
        ALLOCATE(Scal_partArray(nbscal_part),stat=ires)
        IF (ires .NE.0) THEN
           WRITE(6,'(a,i0,a)')"[ERROR] on process ",spec_rank,": not enought memory for Scal_part!"
           RETURN
        ENDIF

        ! get the size and allocate
        do cpt=1, nbscal_part
            write(name_of_field(cpt),'(a,i0)') 'ScalarP',cpt
        end do
        if (.not. data_init_Rarray(ncpus, spec_rank, name_of_field, &
            & 'Number of points for scalar_part' , Scal_partArray)) return

! XXX Todo - Trailer de la nouveauté à venir à la rentrée - XXX
!        IF (parser_is_defined('Part_Vector')) THEN
!            IF (.not. GetNxNyNz('Number of points for scalar_part',nx,ny,nz)) THEN
!                WRITE(6,'(a,i0,x,i0,x,i0)') '[ERROR] Resolution undefined for Scal_part'
!                RETURN
!            END IF
!            if(.not. real_vector_initDataLayout('ScalPartVect', Scal_partVector,nbscal_part,nx,ny,nz, &
!                & Lx,Ly,Lz, ncpus, spec_rank)) then
!                WRITE(6,'(a,i0)') '[ERROR] data_init Failed to init scalar_partVector on ', spec_rank
!                RETURN
!            ENDIF
!        ELSE
!            Scal_partVector%nb_components = 0
!        END IF
! XXX  End Todo
    ENDIF
    !-------------------------------
    ! initialization of the B-field components solved with only pseudo spectral
    ! method
    !-------------------------------
    IF (ASSOCIATED(B)) DEALLOCATE(B)
    B=>null()
    IF (imhd) THEN
        ALLOCATE(B(3),stat=ires)
        IF (ires .NE.0) THEN
            WRITE(6,'(a,i0,a)')"[ERROR] on process ",spec_rank,": not enought memory!"
            RETURN
        ENDIF
        ! get the size and allocate
        name_of_field(1:3) = (/"Bx", "By", "Bz"/)
        IF (.not. data_init_Rarray(ncpus, spec_rank, name_of_field(1:3), &
            & 'Number of points for B-Field' , B)) RETURN
    END IF

    CALL parser_read('mixing threshold', mix_threshold)
    DEALLOCATE(name_of_field)
    data_init=.TRUE.

END FUNCTION data_init


!> Allocate a single real field (with datalayout structure) automatically for the
!! righ mesh size
!! @param[in]       ncpus           = number of MPI processes to distribute datas
!! @param[in]       spec_rank       = my MPI processe rank in the range [0:ncpus-1]
!! @param[in]       name_of_field      = name of the field
!! @param[in]       size_tag        = parser tag to get the right mesh size
!! @param[in]       Lx              = physical size of the domain along X-axis
!! @param[in]       Ly              = physical size of the domain along Y-axis
!! @param[in]       Lz              = physical size of the domain along Z-axis
!! @param[in,out]   field           = real field to allocate
!! @return          success         = .true. if no error occurs
!! @author
!! Jean-Baptiste Lagaert, LEGI
function data_init_Rsingle(ncpus, spec_rank, name_of_field, size_tag, field) result(success)

    integer, intent(in)                     :: ncpus, spec_rank
    character(len=*), intent(in)            :: name_of_field, size_tag
    type(real_data_layout), intent(inout)   :: field
    logical                                 :: success

    ! Local variable
    integer         :: nx, ny, nz ! resolution

    ! init
    success = .false.

    ! get size
    if (.not. GetNxNyNz(trim(size_tag),nx,ny,nz)) then
        write(6,'(a,a,a,i0)')'[ERROR] unable to get resolution for ', name_of_field, ' on ', &
            & spec_rank ,'. Size tag is "', trim(size_tag), '.'
        return
    end if

    !=== ALLOCATE MEMORY FOR THE CHUNK FIELD (FREE DIRECTION IS X BY DEFAULT) ===
    !PAY ATTENTION IF THE FREE DIRECTION IS CHANGED LOOK AT THE FFT 
    if (.not. initDataLayout(trim(name_of_field),field,nx,ny,nz,Lx,Ly,Lz,ncpus,spec_rank)) THEN
        write(6,'(a,a,a,i0)')'[ERROR] data_init Failed to init ', name_of_field, ' on ', spec_rank
        return
    end if

    success = .true.

end function data_init_Rsingle


!> Allocate an array of real field (with datalayout structure) automatically for the
!! righ mesh size
!! @param[in]       ncpus           = number of MPI processes to distribute datas
!! @param[in]       spec_rank       = my MPI processe rank in the range [0:ncpus-1]
!! @param[in]       name_of_field      = name of the fields
!! @param[in]       size_tag        = parser tag to get the right mesh size
!! @param[in]       Lx              = physical size of the domain along X-axis
!! @param[in]       Ly              = physical size of the domain along Y-axis
!! @param[in]       Lz              = physical size of the domain along Z-axis
!! @param[in,out]   field_array     = array of real fields to allocate
!! @return          success         = .true. if no error occurs
!! @author
!! Jean-Baptiste Lagaert, LEGI
function data_init_Rarray(ncpus, spec_rank, name_of_field, size_tag, field_array) result(success)

    integer, intent(in)                                 :: ncpus, spec_rank
    character(len=*), intent(in)                        :: size_tag
    character(len=*), dimension(:), intent(in)          :: name_of_field
    type(real_data_layout), dimension(:), intent(inout) :: field_array
    logical                                             :: success

    ! Local variable
    integer         :: nx, ny, nz   ! resolution
    integer         :: nb_field     ! numlber of field to allocate
    integer         :: i            ! loop indice

    ! init
    success = .false.

    ! Check number of field to allocate
    nb_field = size(field_array)
    if(size(name_of_field)<nb_field) then
        write (6,'(4a,a,a)')'[ERROR] for allocate ', trim(name_of_field(1)), ' until ', &
            & trim(name_of_field(size(name_of_field))), ' not enough field name on rank ', &
            & spec_rank
        return
    end if

    ! get mesh size
    if (.not. GetNxNyNz(trim(size_tag),nx,ny,nz)) then
        write(6,'(a,a,a,i0)')'[ERROR] unable to get resolution for ', name_of_field, ' on ', &
            & spec_rank ,'. Size tag is "', trim(size_tag), '.'
        return
    end if

    ! create storage
    do i = 1, nb_field
        if (.not. initDataLayout(trim(name_of_field(i)),field_array(i),nx,ny,nz,Lx,Ly,Lz,ncpus,spec_rank)) then
            write(6,'(a,a,a,i0)')'[ERROR] data_init Failed to init ', trim(name_of_field(i)), ' on ', spec_rank
            return
        end if
    end do

    success = .true.

end function data_init_Rarray


!------------------------------------------------------------------------------
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!> Deallocate data field.
!!    @param[in,out]    ScalArray       = pointer to an array of scalars storage
!!    @param[in,out]    Scal_partArray  = pointer to an array of scalars (solved with particles solver for advection term) storage
!!    @param[in,out]    U               = longitudinal velocity
!!    @param[in,out]    V               = vertical velocity
!!    @param[in,out]    W               = spanwise velocity
!!    @param[in,out]    P               = pressure
!!    @return           sucess          = booleen equal to .TRUE. if deallocation is successfull
!!@details
!!    Datal field are save into array of REAL_DATA_LAYOUT. This type involve
!!    Pointers which have to be deallocate after the computation. This is done in
!!    this subroutine.
!------------------------------------------------------------------------------
function data_delete(ScalArray, Scal_PartArray, U, V, W, P) result(success)

    use datalayout

    type(REAL_DATA_LAYOUT), intent(inout),pointer, dimension(:) :: ScalArray, Scal_partArray
    type(REAL_DATA_LAYOUT), intent(inout)                       :: U, V, W, P
    logical                                                     :: success
    
    integer     :: i    ! counter for "do" loop

    success = .false.

    if (associated(ScalArray)) then
        do i = 1, size(ScalArray)
            call deleteDataLayout(ScalArray(i))
        end do
        deallocate(ScalArray)
    end if

    if (associated(Scal_partArray)) then
        do i = 1, size(Scal_partArray)
            call deleteDataLayout(Scal_partArray(i))
        end do
        deallocate(Scal_partArray)
    end if

    call deleteDataLayout(U)
    call deleteDataLayout(V)
    call deleteDataLayout(W)
    call deleteDataLayout(P)

    success = .true.

end function data_delete


!------------------------------------------------------------------------------
!> @author
!! Patrick Begou, LEGI
!
!!    @param[in,out]    ScalArray       = pointer to an array of scalars storage
!!    @param[in,out]    Scal_partArray  = pointer to an array of scalars (solved with particles solver for advection term) storage
!!    @param[in,out]    U               = longitudinal velocity
!!    @param[in,out]    V               = vertical velocity
!!    @param[in,out]    W               = spanwise velocity
!!    @param[in,out]    B               = magnetic field
!!    @param[in]        me              = rank of process in cpu pool
!!    @param[in]        nbscal_part     = number of scalar solved with particle method
!!    @param[in]        nbscal          = number of scalar
!!    @param[in]        iFic            = number of file to read (optional(
!!    @return           res             = booleen equal to .TRUE. if deallocation is successfull
!!@details
!------------------------------------------------------------------------------
!=========================================================================================!
FUNCTION data_read(U,V,W,B,imhd,ScalArray, Scal_partArray,nbscal, nbscal_part,me,iFic) result(res)
!=========================================================================================!

USE datalayout
USE parser_tools
USE parallel_tools

IMPLICIT NONE

  !I/O data
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT)                        :: U,V,W
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT), POINTER, DIMENSION(:) :: ScalArray
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT), POINTER, DIMENSION(:) :: Scal_partArray
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT), POINTER, DIMENSION(:) :: B
  INTEGER, INTENT(IN)                                          :: me,nbscal, nbscal_part
  INTEGER, INTENT(IN), OPTIONAL                                :: iFic
  LOGICAL                                                      :: res
  LOGICAL                                                      :: imhd

  !Local data
  CHARACTER(LEN=256)                                           :: filename
  INTEGER                                                      :: sca

  res=.FALSE.

  ! -- Velocity --
  IF(.NOT. dataReadSingleField(U,'U',me,filename,numFic=iFic) ) THEN
      WRITE(6,'(a,i0)') '[ERROR] Failed to read U from '//trim(filename)// ' on processus ', me
      RETURN
  END IF
  IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting U with datas from '//TRIM(filename)

  IF(.NOT. dataReadSingleField(V,'V',me,filename,numFic=iFic) ) THEN
      WRITE(6,'(a,i0)') '[ERROR] Failed to read V from '//trim(filename)// ' on processus ', me
      RETURN
  END IF
  IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting V with datas from '//TRIM(filename)

  IF(.NOT. dataReadSingleField(W,'W',me,filename,numFic=iFic) ) THEN
      WRITE(6,'(a,i0)') '[ERROR] Failed to read W from '//trim(filename)// ' on processus ', me
      RETURN
  END IF
  IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting W with datas from '//TRIM(filename)

  ! -- Magnetic field --
  IF (imhd) THEN
     IF(.NOT. dataReadSingleField(B(1),'Bx',me,filename,numFic=iFic) ) RETURN
     IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting Bx with datas from '//TRIM(filename)

     IF(.NOT. dataReadSingleField(B(2),'By',me,filename,numFic=iFic) ) RETURN
     IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting By with datas from '//TRIM(filename)

     IF(.NOT. dataReadSingleField(B(3),'Bz',me,filename,numFic=iFic) ) RETURN
     IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a)')'[INFO] Overwriting Bz with datas from '//TRIM(filename)
  ENDIF
  
  ! -- Scalar solve with pseudo-spectral method  --
  DO sca=1, nbscal
    IF(.NOT. dataReadSingleField(ScalArray(sca),'Scalar',me,filename,numSca=sca,numFic=iFic) ) RETURN
    IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a,i0,a)')'[INFO] Overwriting scalar n. ',sca,' with datas from ' &
                                                                 & //TRIM(filename)
  ENDDO

  ! -- Scalar solve with spectral/particle method  --
  DO sca=1, nbscal_part
    IF(.NOT. dataReadSingleField(Scal_partArray(sca),'Scalar_part',me,filename,numSca=sca,numFic=iFic) ) RETURN
    IF(me.EQ.0 .AND. LEN_TRIM(filename).GT.0) WRITE(6,'(a,i0,a)')'[INFO] Overwriting scalar_part n. ',sca,' with datas from '&
                                                                 &//TRIM(filename)
  ENDDO
  res=.TRUE.

END FUNCTION data_read


!------------------------------------------------------------------------------
!> @author
!! Antoine Vollant, LEGI
!
!!    @param[in,out]    field           = real data layout
!!    @param[in]i       nameField       = name of field in input file (Bx, U, Scalar ...)
!!    @param[in]        spec_rank       = number of processus in cpu pool 
!!    @param[in]        numSca          = the number of scalar to read (optional)
!!    @param[in]        numFic          = the number of field to read (optional)
!!    @param[in,out]    nameOfFile      = direction of file to read 
!!    @return           res             = booleen equal to .TRUE. if deallocation is successfull
!!@details
!------------------------------------------------------------------------------

FUNCTION dataReadSingleField(field,nameField,spec_rank,nameOfFile,numSca,numFic) result(res)

  USE datalayout
  USE parser_tools
  USE parallel_tools
  use io_interface

  !I/O data
  type(real_data_layout),intent(inout)    :: field
  character(len=*),intent(in)             :: nameField
  integer,intent(in)                      :: spec_rank
  integer,intent(in),optional             :: numSca
  integer,intent(in),optional             :: numFic
  character(len=*),intent(inout)          :: nameOfFile
  logical                                 :: res

  !Local data
  character(len=100)                      :: request, request2
  character(len=str_medium)               :: nameOfField

  res = .false.

  IF     (present(numSca) .AND. present(numFic) ) THEN
      WRITE(request,'(a,1x,i0,1x,a,1x,i0)') 'Data file',numFic,'to read for '//trim(nameField),numSca
      WRITE(request2,'(a,1x,i0,1x,a,1x,i0)') 'Dataset',numFic,'to read for '//trim(nameField),numSca
  ELSEIF (present(numSca) .AND. .NOT. present(numFic)) THEN
      WRITE(request,'(a,1x,i0)') 'Data file to read for '//trim(nameField),numSca
      WRITE(request2,'(a,1x,i0)') 'Dataset to read for '//trim(nameField),numSca
  ELSEIF (.NOT. present(numSca) .AND. present(numFic)) THEN
      WRITE(request,'(a,1x,i0,1x,a)') 'Data file',numFic,'to read for '//trim(nameField)
      WRITE(request2,'(a,1x,i0,1x,a)') 'Dataset',numFic,'to read for '//trim(nameField)
  ELSE
      WRITE(request,'(a)') 'Data file to read for '//trim(nameField)
      WRITE(request2,'(a)') 'Dataset to read for '//trim(nameField)
  ENDIF

  nameOfFile = " "
  nameOfField = " "
  IF(parser_is_defined(request)) THEN
    CALL parser_read(request,nameOfFile)
    IF (LEN_TRIM(nameOfFile).GT.0) THEN
        if(parser_is_defined(request2)) then
            call parser_read(request2,nameOfField)
            call read_datalayout_scalar(spec_rank,nameOfFile,field, res,nameOfField)
        else
            call read_datalayout_scalar(spec_rank,nameOfFile,field, res)
        end if
        IF (.NOT. res) THEN
            write (6,'(a,1x,a)') '[ERROR] indataReadSingleField for',trim(nameOfFile)
          RETURN
        ENDIF
    ENDIF
  ENDIF

  res = .true.

END FUNCTION dataReadSingleField

END MODULE data
