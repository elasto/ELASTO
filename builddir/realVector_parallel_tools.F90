!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST io
!------------------------------------------------------------------------------
! MODULE realVectorparallel
!
!> @author
!> Patrick BEGOU, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to realize all parallel communications for
!> changing data layout, for calculations and for parallel I/O. This is a generic version
!> wich should be processed to create a REAL version and a COMPLEX version.
!> DO NOT CHANGE this file: it is automagicaly created from implementparallel.Fortran
!> for REAL and COMPLEX data types. All changes must be done in implementparallel.Fortran
!> file only.
!------------------------------------------------------------------------------
MODULE realVectorparallel

    use realparallel

    IMPLICIT NONE

    !private subroutines & functions
    PRIVATE real_AlongXtoAlongY
    PRIVATE real_AlongYtoAlongX
    PRIVATE real_AlongYtoAlongZ
    PRIVATE real_AlongZtoAlongY

CONTAINS
  
!------------------------------------------------------------------------------
!> Surcharge of real_setAlongX for vector.
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the X axis and distributed on Y and Z among the processes. It use
!> the real_setAlongX routines to do its job.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION real_setAlongX_Vector(spec_rank,val)
    USE datalayout
    USE communicators_tools
    !USE realparallel
    USE realVectordatalayout

    TYPE(REAL_VECTOR_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)                         :: spec_rank
    LOGICAL                                     :: real_setAlongX_Vector

    TYPE(REAL_DATA_LAYOUT)  :: comp_data
    INTEGER                 :: comp
    

    real_setAlongX_Vector=.FALSE.

    IF (.NOT. CommInitialized()) THEN
       WRITE(6,'(a)')"Error real_setAlongX: Call initCommunicators first!"
       RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        real_setAlongX_Vector=.TRUE.
        RETURN
    CASE (alongY)
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, .true.)
            if(.not. real_AlongYtoAlongX(spec_rank,comp_data)) RETURN
        end do
    CASE (alongZ)
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, .true.)
            if(.not. real_AlongZtoAlongY(spec_rank,comp_data)) RETURN
            if(.not. real_AlongYtoAlongX(spec_rank,comp_data)) RETURN
        end do
    END SELECT

    val%xmin=1
    val%xmax=val%nx
    val%ymin=val%Y4layoutOnX(spec_rank+1,1)
    val%ymax=val%Y4layoutOnX(spec_rank+1,2)
    val%zmin=val%Z4layoutOnX(spec_rank+1,1)
    val%zmax=val%Z4layoutOnX(spec_rank+1,2)
    val%curlayout=alongX

    real_setAlongX_Vector=.TRUE.

END FUNCTION real_setAlongX_Vector

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the Y axis and distributed on X and Z among the processes.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION real_setAlongY_Vector(spec_rank,val)
    USE datalayout
    USE communicators_tools
    USE realVectordatalayout

    IMPLICIT NONE
    TYPE(REAL_VECTOR_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)                         :: spec_rank
    LOGICAL                                     :: real_setAlongY_Vector
    
    TYPE(REAL_DATA_LAYOUT)  :: comp_data
    INTEGER                 :: comp

    real_setAlongY_Vector=.FALSE.

    IF (.NOT. CommInitialized()) THEN
        WRITE(6,'(a)')"Error real_setAlongY: Call initCommunicators first!"
        RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, .true.)
            if(.not. real_AlongXtoAlongY(spec_rank,comp_data)) RETURN
        end do
    CASE (alongY)
        real_setAlongY_Vector=.TRUE.
        RETURN
    CASE (alongZ)
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, .true.)
            if(.not. real_AlongZtoAlongY(spec_rank,comp_data)) RETURN
        end do
    END SELECT

    val%xmin=val%X4layoutOnY(spec_rank+1,1)
    val%xmax=val%X4layoutOnY(spec_rank+1,2)
    val%ymin=1
    val%ymax=val%ny
    val%zmin=val%Z4layoutOnY(spec_rank+1,1)
    val%zmax=val%Z4layoutOnY(spec_rank+1,2)
    val%curlayout=alongY

    real_setAlongY_Vector=.TRUE.

END FUNCTION real_setAlongY_Vector


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the Y axis and distributed on X and Z among the processes.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION real_setAlongZ_Vector(spec_rank,val)
    USE datalayout
    USE communicators_tools
    USE realVectordatalayout

    IMPLICIT NONE
    TYPE(REAL_VECTOR_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)             :: spec_rank
    LOGICAL                         :: real_setAlongZ_Vector

    TYPE(REAL_DATA_LAYOUT)  :: comp_data
    INTEGER                 :: comp

    real_setAlongZ_Vector=.FALSE.

    IF (.NOT. CommInitialized()) THEN
        WRITE(6,'(a)')"Error real_setAlongZ: Call initCommunicators first!"
        RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, guru=.true.)
            if(.not. real_AlongXtoAlongY(spec_rank,comp_data)) RETURN
            if(.not. real_AlongYtoAlongZ(spec_rank,comp_data)) RETURN
        end do
    CASE (alongY)
        do comp = 1, val%nb_components
            call real_vector_manipulate_comp(val,comp, comp_data, guru=.true.)
            if(.not. real_AlongYtoAlongZ(spec_rank,comp_data)) RETURN
        end do
    CASE (alongZ)
        real_setAlongZ_Vector=.TRUE.
        RETURN
    END SELECT

    val%xmin=val%X4layoutOnY(spec_rank+1,1)
    val%xmax=val%X4layoutOnY(spec_rank+1,2)
    val%ymin=1
    val%ymax=val%ny
    val%zmin=val%Z4layoutOnY(spec_rank+1,1)
    val%zmax=val%Z4layoutOnY(spec_rank+1,2)
    val%curlayout=alongY

    real_setAlongZ_Vector=.TRUE.

END FUNCTION real_setAlongZ_Vector

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along X and are reorganized to be contiguous along Y 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
  FUNCTION real_AlongXtoAlongY(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(REAL_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             ::spec_rank
  LOGICAL                         ::real_AlongXtoAlongY
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  REAL(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  real_AlongXtoAlongY=.FALSE.
  !The idea is to create two arrays containing the dimensions of 
  !  - the subset of data to send to each processor (whatIsend(0:ncpus,6)
  !  - the subset of data I receive from each processor (whatIget(0:ncpus-1,6)
  ! The index correspond to the rank in the communicator
  ALLOCATE(whatIsend(0:SIZE(XY_Nodes)-1,6),whatIget(0:SIZE(XY_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR real_AlongXtoAlongY on process ",spec_rank,": not enought memory for internal arrays!"
     RETURN
  ENDIF
  !now fill in the arrays
  DO i=0,SIZE(XY_Nodes)-1
     whatIsend(i,1)=val%X4layoutOnY(XY_Nodes(i+1)+1,1) !xmin to send to rank 0 of local communicator
     whatIsend(i,2)=val%X4layoutOnY(XY_Nodes(i+1)+1,2) !xmax to send to rank 0 of local communicator
     whatIsend(i,3)=val%Y4layoutOnX(spec_rank+1,1)
     whatIsend(i,4)=val%Y4layoutOnX(spec_rank+1,2)
     whatIsend(i,5)=val%Z4layoutOnX(spec_rank+1,1)
     whatIsend(i,6)=val%Z4layoutOnX(spec_rank+1,2)
     
     whatIget(i,1)=val%X4layoutOnY(spec_rank+1,1)
     whatIget(i,2)=val%X4layoutOnY(spec_rank+1,2)
     whatIget(i,3)=val%Y4layoutOnX(XY_Nodes(i+1)+1,1)
     whatIget(i,4)=val%Y4layoutOnX(XY_Nodes(i+1)+1,2)
     whatIget(i,5)=val%Z4layoutOnX(spec_rank+1,1)
     whatIget(i,6)=val%Z4layoutOnX(spec_rank+1,2)
  ENDDO
#ifdef DEBUGMPI  
! this debug prints the content of the arrays
  DO i=0,SIZE(XY_Nodes)-1
     WRITE(10+spec_rank,100) spec_rank," sends to ",XY_Nodes(i+1),whatIsend(i,:)
  ENDDO
  DO i=0,SIZE(XY_Nodes)-1
     WRITE(10+spec_rank,100) spec_rank," gets from ",XY_Nodes(i+1),whatIget(i,:)
  ENDDO
100 format &
   &("process ",i0,a," process ",i0," [",2(i0,":",i0,","), 1(i0,":",i0,"]"))
#endif


!  !allocate receive buffer  along Y
  ALLOCATE(values(val%X4layoutOnY(spec_rank+1,1):val%X4layoutOnY(spec_rank+1,2), &
                & 1:val%ny, &
	        & val%Z4layoutOnY(spec_rank+1,1):val%Z4layoutOnY(spec_rank+1,2)), &
              & stat=res)
  IF(res.NE.0) THEN
       WRITE(6,'(a,i0,a)')"ERROR real_AlongXtoAlongY on process ",spec_rank,": not enought spec_rankmory for datas!"
       WRITE(6,101)  spec_rank, & 
                  & val%X4layoutOnY(spec_rank+1,1),val%X4layoutOnY(spec_rank+1,2), &
                  & 1,val%ny, &
		  & val%Z4layoutOnY(spec_rank+1,1),val%Z4layoutOnY(spec_rank+1,2), &
                  & " has FAILED"
101 format ("On process ",i0," allocation of values: [",2(i0,":",i0,","), &
    & 1(i0,":",i0,"]"),a)
       RETURN
  ENDIF

#ifdef DEBUGMPI  
       WRITE(10+spec_rank,101)  spec_rank, & 
                  & val%X4layoutOnY(spec_rank+1,1),val%X4layoutOnY(spec_rank+1,2), &
                  & 1,val%ny, &
		  & val%Z4layoutOnY(spec_rank+1,1),val%Z4layoutOnY(spec_rank+1,2), &
                  & " OK"
#endif

  ! now send & receive the datas
  DO i=0,SIZE(XY_Nodes)-1
    IF (spec_rank .EQ. XY_Nodes(i+1)) THEN
       !just copy
       values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &=&
       &val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))
    ELSE
       !send and receive data
       cptin=(whatIget(i,2)-whatIget(i,1)+1)*(whatIget(i,4)-whatIget(i,3)+1)*(whatIget(i,6)-whatIget(i,5)+1)
       cptout=(whatIsend(i,2)-whatIsend(i,1)+1)*(whatIsend(i,4)-whatIsend(i,3)+1)*(whatIsend(i,6)-whatIsend(i,5)+1)
#ifdef DEBUGMPI
       WRITE(10+spec_rank,110) spec_rank," sends  ",cptout, &
       & whatIsend(i,1),whatIsend(i,2),whatIsend(i,3),whatIsend(i,4),whatIsend(i,5),whatIsend(i,6),&
       & "to   ",XY_Nodes(i+1)
       WRITE(10+spec_rank,110) spec_rank," stores ",&
       & cptin,whatIget(i,1),whatIget(i,2),whatIget(i,3),whatIget(i,4),whatIget(i,5),whatIget(i,6),&
       & "from ",XY_Nodes(i+1)
110    FORMAT("process ",i0,a,i0,"values [",i3,":",i3,",",i3,":",i3,",",i3,":",i3,"] ", &
       & a,i0) 
#endif
       IF(.NOT. real_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_REAL_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_REAL_WP,i,MPI_ANY_TAG,XY_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure
  val%values = values
  DEALLOCATE(values)
  values=>null()
  val%xmin=val%X4layoutOnY(spec_rank+1,1)
  val%xmax=val%X4layoutOnY(spec_rank+1,2)
  val%ymin=1
  val%ymax=val%ny
  val%zmin=val%Z4layoutOnY(spec_rank+1,1)
  val%zmax=val%Z4layoutOnY(spec_rank+1,2)
  val%curlayout=alongY
  
  DEALLOCATE(whatIsend,whatIget)
  
  real_AlongXtoAlongY=.TRUE.
  RETURN
  
  END FUNCTION real_AlongXtoAlongY
  
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Y and are reorganized to be contiguous along X 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
! See real_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION real_AlongYtoAlongX(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(REAL_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: real_AlongYtoAlongX
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  REAL(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  real_AlongYtoAlongX=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(XY_Nodes)-1,6),whatIget(0:SIZE(XY_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR real_AlongYtoAlongX on process ",spec_rank,": not enought memory for internal arrays!"
     RETURN
  ENDIF
  
  !now fill in the arrays
  DO i=0,SIZE(XY_Nodes)-1
     whatIsend(i,1)=val%X4layoutOnY(spec_rank+1,1)
     whatIsend(i,2)=val%X4layoutOnY(spec_rank+1,2)
     whatIsend(i,3)=val%Y4layoutOnX(XY_Nodes(i+1)+1,1)
     whatIsend(i,4)=val%Y4layoutOnX(XY_Nodes(i+1)+1,2)
     whatIsend(i,5)=val%Z4layoutOnX(spec_rank+1,1)
     whatIsend(i,6)=val%Z4layoutOnX(spec_rank+1,2)
     
     whatIget(i,1)=val%X4layoutOnY(XY_Nodes(i+1)+1,1) 
     whatIget(i,2)=val%X4layoutOnY(XY_Nodes(i+1)+1,2)
     whatIget(i,3)=val%Y4layoutOnX(spec_rank+1,1)
     whatIget(i,4)=val%Y4layoutOnX(spec_rank+1,2)
     whatIget(i,5)=val%Z4layoutOnX(spec_rank+1,1)
     whatIget(i,6)=val%Z4layoutOnX(spec_rank+1,2)
  ENDDO

  !allocate receive buffer  along X
  ALLOCATE(values(1:val%nx, &
                & val%Y4layoutOnX(spec_rank+1,1):val%Y4layoutOnX(spec_rank+1,2), &
	        & val%Z4layoutOnX(spec_rank+1,1):val%Z4layoutOnX(spec_rank+1,2)), &
              & stat=res)
  IF(res.NE.0) THEN
       WRITE(6,'(a,i0,a)')"ERROR real_AlongYtoAlongX on process ",spec_rank,": not enought memory for datas!"
       WRITE(6,101)  spec_rank, & 
                  & val%X4layoutOnY(spec_rank+1,1),val%X4layoutOnY(spec_rank+1,2), &
                  & 1,val%ny, &
		  & val%Z4layoutOnY(spec_rank+1,1),val%Z4layoutOnY(spec_rank+1,2), &
                  & " has FAILED"
101 format ("On process ",i0," allocation of values: [",2(i0,":",i0,","), &
    & 1(i0,":",i0,"]"),a)
       RETURN
  ENDIF

  ! now send & receive the datas
  DO i=0,SIZE(XY_Nodes)-1
    IF (spec_rank .EQ. XY_Nodes(i+1)) THEN
       !just copy
       values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &=&
       &val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))
    ELSE
       !send and receive data
       cptin=(whatIget(i,2)-whatIget(i,1)+1)*(whatIget(i,4)-whatIget(i,3)+1)*(whatIget(i,6)-whatIget(i,5)+1)
       cptout=(whatIsend(i,2)-whatIsend(i,1)+1)*(whatIsend(i,4)-whatIsend(i,3)+1)*(whatIsend(i,6)-whatIsend(i,5)+1)

       IF (.NOT. real_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_REAL_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_REAL_WP,i,MPI_ANY_TAG,XY_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  val%values = values
  DEALLOCATE(values)
  values=>null()
  val%xmin=1
  val%xmax=val%nx
  val%ymin=val%Y4layoutOnX(spec_rank+1,1)
  val%ymax=val%Y4layoutOnX(spec_rank+1,2)
  val%zmin=val%Z4layoutOnX(spec_rank+1,1)
  val%zmax=val%Z4layoutOnX(spec_rank+1,2)
  val%curlayout=alongX
  
  DEALLOCATE(whatIsend,whatIget)
  
  real_AlongYtoAlongX=.TRUE.
  RETURN
  
  END FUNCTION real_AlongYtoAlongX

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Y and are reorganized to be contiguous along Z 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communicator
! See real_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION real_AlongYtoAlongZ(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(REAL_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: real_AlongYtoAlongZ
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  REAL(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  real_AlongYtoAlongZ=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(YZ_Nodes)-1,6),whatIget(0:SIZE(YZ_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR real_AlongYtoAlongZ on process ",spec_rank,": not enought memory for internal arrays!"
     RETURN
  ENDIF
  
  !now fill in the arrays
  DO i=0,SIZE(YZ_Nodes)-1
     whatIsend(i,1)=val%X4layoutOnY(spec_rank+1,1)
     whatIsend(i,2)=val%X4layoutOnY(spec_rank+1,2)
     whatIsend(i,3)=val%Y4layoutOnZ(YZ_Nodes(i+1)+1,1)
     whatIsend(i,4)=val%Y4layoutOnZ(YZ_Nodes(i+1)+1,2)
     whatIsend(i,5)=val%Z4layoutOnY(spec_rank+1,1)
     whatIsend(i,6)=val%Z4layoutOnY(spec_rank+1,2)
     
     whatIget(i,1)=val%X4layoutOnZ(spec_rank+1,1) 
     whatIget(i,2)=val%X4layoutOnZ(spec_rank+1,2)
     whatIget(i,3)=val%Y4layoutOnZ(spec_rank+1,1)
     whatIget(i,4)=val%Y4layoutOnZ(spec_rank+1,2)
     whatIget(i,5)=val%Z4layoutOnY(YZ_Nodes(i+1)+1,1)
     whatIget(i,6)=val%Z4layoutOnY(YZ_Nodes(i+1)+1,2)
  ENDDO

  !allocate receive buffer  along X
  ALLOCATE(values(val%X4layoutOnZ(spec_rank+1,1):val%X4layoutOnZ(spec_rank+1,2), &
	        & val%Y4layoutOnZ(spec_rank+1,1):val%Y4layoutOnZ(spec_rank+1,2), &
		& 1:val%nz), &
              & stat=res)
  IF(res.NE.0) THEN
       WRITE(6,'(a,i0,a)')"ERROR real_AlongYtoAlongZ on process ",spec_rank,": not enought memory for datas!"
       WRITE(6,101)  spec_rank, & 
                  & val%X4layoutOnZ(spec_rank+1,1),val%X4layoutOnZ(spec_rank+1,2), &
	          & val%Y4layoutOnZ(spec_rank+1,1),val%Y4layoutOnZ(spec_rank+1,2), &
		  & 1,val%nz, &
                  & " has FAILED"
101 format ("On process ",i0," allocation of values: [",2(i0,":",i0,","), &
    & 1(i0,":",i0,"]"),a)
       RETURN
  ENDIF

  ! now send & receive the datas
  DO i=0,SIZE(YZ_Nodes)-1
    IF (spec_rank .EQ. YZ_Nodes(i+1)) THEN
       !just copy
       values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &=&
       &val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))
    ELSE
       !send and receive data
       cptin=(whatIget(i,2)-whatIget(i,1)+1)*(whatIget(i,4)-whatIget(i,3)+1)*(whatIget(i,6)-whatIget(i,5)+1)
       cptout=(whatIsend(i,2)-whatIsend(i,1)+1)*(whatIsend(i,4)-whatIsend(i,3)+1)*(whatIsend(i,6)-whatIsend(i,5)+1)

       IF (.NOT. real_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_REAL_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_REAL_WP,i,MPI_ANY_TAG,YZ_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  val%values = values
  DEALLOCATE(values)
  values=>null()
  val%xmin=val%X4layoutOnZ(spec_rank+1,1)
  val%xmax=val%X4layoutOnZ(spec_rank+1,2)
  val%ymin=val%Y4layoutOnZ(spec_rank+1,1)
  val%ymax=val%Y4layoutOnZ(spec_rank+1,2)
  val%zmin=1
  val%zmax=val%nz
  val%curlayout=alongZ
  
  DEALLOCATE(whatIsend,whatIget)
  
  real_AlongYtoAlongZ=.TRUE.
  RETURN
  
  END FUNCTION real_AlongYtoAlongZ


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the REAL_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Z and are reorganized to be contiguous along Y 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the REAL_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communicator
! See real_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION real_AlongZtoAlongY(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(REAL_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: real_AlongZtoAlongY
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  REAL(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  real_AlongZtoAlongY=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(YZ_Nodes)-1,6),whatIget(0:SIZE(YZ_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR real_AlongZtoAlongY on process ",spec_rank,": not enought memory for internal arrays!"
     RETURN
  ENDIF
  
  !now fill in the arrays
  DO i=0,SIZE(YZ_Nodes)-1
     whatIget(i,1)=val%X4layoutOnY(spec_rank+1,1)
     whatIget(i,2)=val%X4layoutOnY(spec_rank+1,2)
     whatIget(i,3)=val%Y4layoutOnZ(YZ_Nodes(i+1)+1,1)
     whatIget(i,4)=val%Y4layoutOnZ(YZ_Nodes(i+1)+1,2)
     whatIget(i,5)=val%Z4layoutOnY(spec_rank+1,1)
     whatIget(i,6)=val%Z4layoutOnY(spec_rank+1,2)
     
     whatIsend(i,1)=val%X4layoutOnZ(spec_rank+1,1) 
     whatIsend(i,2)=val%X4layoutOnZ(spec_rank+1,2)
     whatIsend(i,3)=val%Y4layoutOnZ(spec_rank+1,1)
     whatIsend(i,4)=val%Y4layoutOnZ(spec_rank+1,2)
     whatIsend(i,5)=val%Z4layoutOnY(YZ_Nodes(i+1)+1,1)
     whatIsend(i,6)=val%Z4layoutOnY(YZ_Nodes(i+1)+1,2)
  ENDDO

  !allocate receive buffer  along X
  ALLOCATE(values(val%X4layoutOnY(spec_rank+1,1):val%X4layoutOnY(spec_rank+1,2), &
                & 1:val%ny, &
	        & val%Z4layoutOnY(spec_rank+1,1):val%Z4layoutOnY(spec_rank+1,2)), &
              & stat=res)
  IF(res.NE.0) THEN
       WRITE(6,'(a,i0,a)')"ERROR real_AlongZtoAlongY on process ",spec_rank,": not enought memory for datas!"
       WRITE(6,101)  spec_rank, & 
                  & val%X4layoutOnY(spec_rank+1,1),val%X4layoutOnY(spec_rank+1,2), &
                  & 1,val%ny, &
		  & val%Z4layoutOnY(spec_rank+1,1),val%Z4layoutOnY(spec_rank+1,2), &
                  & " has FAILED"
101 format ("On process ",i0," allocation of values: [",2(i0,":",i0,","), &
    & 1(i0,":",i0,"]"),a)
       RETURN
  ENDIF

  ! now send & receive the datas
  DO i=0,SIZE(YZ_Nodes)-1
    IF (spec_rank .EQ. YZ_Nodes(i+1)) THEN
       !just copy
       values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &=&
       &val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))
    ELSE
       !send and receive data
       cptin=(whatIget(i,2)-whatIget(i,1)+1)*(whatIget(i,4)-whatIget(i,3)+1)*(whatIget(i,6)-whatIget(i,5)+1)
       cptout=(whatIsend(i,2)-whatIsend(i,1)+1)*(whatIsend(i,4)-whatIsend(i,3)+1)*(whatIsend(i,6)-whatIsend(i,5)+1)

       IF (.NOT. real_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_REAL_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_REAL_WP,i,MPI_ANY_TAG,YZ_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  val%values = values
  DEALLOCATE(values)
  values=>null()
  val%xmin=val%X4layoutOnY(spec_rank+1,1)
  val%xmax=val%X4layoutOnY(spec_rank+1,2)
  val%ymin=1
  val%ymax=val%ny
  val%zmin=val%Z4layoutOnY(spec_rank+1,1)
  val%zmax=val%Z4layoutOnY(spec_rank+1,2)
  val%curlayout=alongY
  
  DEALLOCATE(whatIsend,whatIget)
  
  real_AlongZtoAlongY=.TRUE.
  RETURN
  
  END FUNCTION real_AlongZtoAlongY
  
  
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function dump in a binary file the values of the REAL_DATA_LAYOUT 
!!> variable: val. If the file exist, it is overwriten. The file format is
!!> stream (C like binary format). It contains 3 integers, for the full array
!!> dimension, followed by the datas.
!! 
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @param[in] name the name of the output file.
!!> @param[in] val the REAL_DATA_LAYOUT variable to save in the file.
!!> @return .TRUE. if file is successfuly written.
!!------------------------------------------------------------------------------
!FUNCTION real_WriteToFile(spec_rank,name,val)
!!
!! Parallel write of a 3D scalar array in one file named 'name'. All the processes write in the same file.
!! U, V, W and SC are known as pointer to 3D array decalred in "data" module.
!!-------------------------------------------------------------------------
!
!  USE datalayout
!  USE communicators_tools
!  USE mpi
!
!  IMPLICIT NONE
!!  INCLUDE 'mpif.h'
!  INTEGER, INTENT(IN)               :: spec_rank
!  CHARACTER(LEN=*), INTENT(IN)      :: name
!  TYPE(REAL_DATA_LAYOUT),INTENT(IN) :: val
!  LOGICAL                           :: real_WriteToFile
!
!   
!  INTEGER :: gsizes(3)
!  !> the size of the global array in the 3 dimensions.
!  INTEGER :: lsizes(3)
!  !> the size of the local array in the 3 dimensions.
!  INTEGER :: offset(3)
!  !> the offset of the local array in the global array for the 3 dimensions.
!  
!  INTEGER :: iunit, fileview, ierr
!  
!  !Offset in the file for the data: 3 x 4 bytes integers.
!  INTEGER(KIND=MPI_OFFSET_KIND) :: disp=12     
!  !> the offset (in bytes) of the datas in the file to skip the header.
!  
!  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
!  
!  real_WriteToFile=.FALSE.
!  
!  ! Create the view for parallel I/O
!  gsizes(1) = val%nx
!  gsizes(2) = val%ny
!  gsizes(3) = val%nz
!  lsizes(1) = val%xmax - val%xmin + 1
!  lsizes(2) = val%ymax - val%ymin + 1
!  lsizes(3) = val%zmax - val%zmin + 1
!  offset(1) = val%xmin - 1 
!  offset(2) = val%ymin - 1
!  offset(3) = val%zmin - 1
!  CALL MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,offset,MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": unable to create subarray!"
!    RETURN
!  ENDIF
!  
!  CALL MPI_TYPE_COMMIT(fileview,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": subarray type commit fails!"
!    RETURN
!  ENDIF
!  
!  ! Open file and go to the right position corresponding to our local subdomain
!  CALL MPI_FILE_OPEN(spec_communicator,trim(name),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,iunit,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,": cannot open file ["//trim(name)//"]"
!    RETURN
!  ENDIF
!  
!  ! The processus of rank 0 create the header
!  IF (spec_rank .EQ. 0) THEN
!     CALL MPI_FILE_WRITE(iunit,gsizes,3,MPI_INTEGER,status,ierr)
!     IF (ierr .NE. MPI_SUCCESS) THEN
!       WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
!                         &": cannot write dimensions in file ["//trim(name)//"]"
!       RETURN
!     ENDIF
!  ENDIF   
!  
!
!  ! Positionnement
!  CALL MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
!                      &": cannot set view on file ["//trim(name)//"]"
!    RETURN
!  ENDIF
!  
!  ! Write data
!  CALL MPI_FILE_WRITE_ALL(iunit,val%values,lsizes(1)*lsizes(2)*lsizes(3),MPI_REAL_WP,status,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
!                      &": write error for file ["//trim(name)//"]"
!    RETURN
!  ENDIF
!
!  ! Close file
!  CALL MPI_FILE_CLOSE(iunit,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_WriteToFile on process ",spec_rank,&
!                      &": close error on file ["//trim(name)//"]"
!    RETURN
!  ENDIF
!
!  !All is OK now!
!  real_WriteToFile=.TRUE.
!  RETURN
!
!END FUNCTION real_WriteToFile
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Antoine Vollant
!!
!!> @details
!!> This function makes a local sommation of the variable "value" on returns the
!!> result on node 0.
!! 
!!> @param[in] value the local variable to sum.
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @param[in] node the rank number of the root process in the range [0,ncpus - 1]
!!>            by default this value is at 0 (optional).
!!> @return the sommation of all the value variable of each node.
!!------------------------------------------------------------------------------
!FUNCTION real_DoLocalSum(spec_rank,value,node)
!!------------------------------------------------------------------------------
!
!  USE precision_tools
!  USE communicators_tools
!  USE mpi
!  
!  IMPLICIT NONE
!!  INCLUDE 'mpif.h'
!  
!  REAL(WP) ::real_DoLocalSum 
!  REAL(WP), INTENT(IN)          :: value
!  INTEGER, INTENT(IN)           :: spec_rank
!  INTEGER, INTENT(IN),OPTIONAL  :: node
!
!  REAL(WP) :: resultat
!  INTEGER  ::ierr
!  
!  IF ( PRESENT(node) ) THEN
!    CALL MPI_REDUCE(value,resultat,1,MPI_REAL_WP,MPI_SUM,node,spec_communicator,ierr)
!  ELSE
!    CALL MPI_REDUCE(value,resultat,1,MPI_REAL_WP,MPI_SUM,0,spec_communicator,ierr)
!  ENDIF
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_DoLocalSum on process ",spec_rank,&
!                      &": Unable to run local sum."
!    CALL MPI_FINALIZE(ierr)
!    STOP 'FAILED'
!  ENDIF
!  real_DoLocalSum=resultat
!  RETURN
!  
!END FUNCTION real_DoLocalSum
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function makes a global sommation of the variable "value" on returns the
!!> result on every node.
!! 
!!> @param[in] value the local variable to sum.
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @return the sommation of all the value variable of each node.
!!------------------------------------------------------------------------------
!FUNCTION real_DoGlobalSum(spec_rank,value)
!!------------------------------------------------------------------------------
!
!  USE precision_tools
!  USE communicators_tools
!  USE mpi
!  
!  IMPLICIT NONE
!!  INCLUDE 'mpif.h'
!  
!  REAL(WP) ::real_DoGlobalSum 
!  REAL(WP), INTENT(IN) :: value
!  INTEGER, INTENT(IN)  :: spec_rank
!
!  REAL(WP) :: resultat
!  INTEGER  ::ierr
!  
!  CALL MPI_ALLREDUCE(value,resultat,1,MPI_REAL_WP,MPI_SUM,spec_communicator,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalSum on process ",spec_rank,&
!                      &": Unable to run global sum."
!    CALL MPI_FINALIZE(ierr)
!    STOP 'FAILED'
!  ENDIF
!  real_DoGlobalSum=resultat
!  RETURN
!  
!END FUNCTION real_DoGlobalSum
!
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function makes a global sommation of the 1D array "vector" on returns the
!!> result on every node.
!! 
!!> @param[in] vector the local array to sum.
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @return the sommation of all the value of vector for each node.
!!------------------------------------------------------------------------------
!FUNCTION real_DoGlobalVectorSum(spec_rank,vector)
!!------------------------------------------------------------------------------
!
!  USE precision_tools
!  USE communicators_tools
!  USE mpi
!  
!  IMPLICIT NONE
!!  INCLUDE 'mpif.h'
!  
!  REAL(WP), INTENT(IN) :: vector(:)
!  REAL(WP) ::real_DoGlobalVectorSum (size(vector))
!  INTEGER, INTENT(IN)  :: spec_rank
!
!!  REAL(WP) :: resultat
!  INTEGER  ::ierr
!  
!  CALL MPI_ALLREDUCE(vector,real_DoGlobalVectorSum,size(vector),MPI_REAL_WP,MPI_SUM,spec_communicator,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalVectorSum on process ",spec_rank,&
!                      &": Unable to run global sum on the array."
!    CALL MPI_FINALIZE(ierr)
!    STOP 'FAILED'
!  ENDIF
!!  real_DoGlobalSum=resultat
!  RETURN
!  
!END FUNCTION real_DoGlobalVectorSum
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function makes a global max of the variable "value" and returns the
!!> result on every node.
!! 
!!> @param[in] value the local variable to test.
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @return the maximum of all the value variable over each node.
!!------------------------------------------------------------------------------
!FUNCTION real_DoGlobalMax(spec_rank,value)
!!------------------------------------------------------------------------------
!
!  USE precision_tools
!  USE communicators_tools
!  USE mpi
!  
!  IMPLICIT NONE
!  
!  REAL(WP) ::real_DoGlobalMax 
!  REAL(WP), INTENT(IN) :: value
!  INTEGER, INTENT(IN)  :: spec_rank
!
!  REAL(WP) :: resultat
!  INTEGER  ::ierr
!  
!  CALL MPI_ALLREDUCE(value,resultat,1,MPI_REAL_WP,MPI_MAX,spec_communicator,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalMax on process ",spec_rank,&
!                      &": Unable to run global Max."
!    CALL MPI_FINALIZE(ierr)
!    STOP 'FAILED'
!  ENDIF
!  real_DoGlobalMax=resultat
!  RETURN
!  
!END FUNCTION real_DoGlobalMax
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function makes a global max of the variable "value" and returns the
!!> result on every node.
!! 
!!> @param[in] value the local variable to test.
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @return the maximum of all the value variable over each node.
!!------------------------------------------------------------------------------
!FUNCTION real_DoGlobalMin(spec_rank,value)
!!------------------------------------------------------------------------------
!
!  USE precision_tools
!  USE communicators_tools
!  USE mpi
!  
!  IMPLICIT NONE
!  
!  REAL(WP) ::real_DoGlobalMin 
!  REAL(WP), INTENT(IN) :: value
!  INTEGER, INTENT(IN)  :: spec_rank
!
!  REAL(WP) :: resultat
!  INTEGER  ::ierr
!  
!  CALL MPI_ALLREDUCE(value,resultat,1,MPI_REAL_WP,MPI_MIN,spec_communicator,ierr)
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalMin on process ",spec_rank,&
!                      &": Unable to run global Min."
!    CALL MPI_FINALIZE(ierr)
!    STOP 'FAILED'
!  ENDIF
!  real_DoGlobalMin=resultat
!  RETURN
!  
!END FUNCTION real_DoGlobalMin
!
!!------------------------------------------------------------------------------
!!> @author 
!!> Patrick BEGOU
!!
!!> @details
!!> This function makes a global sommation of the 1D array "vector" on returns the
!!> result on every node.
!! 
!!> @param[in] val a REAL_DATA_LAYOUT from wich we want to extract a slice in the [Y,Z] plane
!!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!!> @param[in] nslice the slice position. defauklt value is 1 if not provided.
!!> @return a 2D [ny,nz] array containing the requested slice or exit.
!!------------------------------------------------------------------------------
!FUNCTION real_extractYZSlice(val,spec_rank,nslice) RESULT(slice)
!!------------------------------------------------------------------------------
!    USE precision_tools
!    USE datalayout
!    USE communicators_tools
!    USE mpi
!
!    IMPLICIT NONE
!    TYPE(REAL_DATA_LAYOUT),INTENT(IN) :: val
!    INTEGER, INTENT(IN)               :: spec_rank
!    INTEGER, INTENT(IN), OPTIONAL     :: nslice
!    REAL(WP), DIMENSION(:,:), ALLOCATABLE::buff
!    REAL(WP), DIMENSION(val%ny,val%nz):: slice
!
!    INTEGER:: curslice, i, ierr
!    !> Rank in spec_communicator when i is the rank in  YZ_COMMUNICATOR
!    INTEGER::rankInCW
!
!    IF(PRESENT(nslice)) THEN
!        IF(nslice.GE.1 .AND. nslice.LE.val%nx) THEN
!            curslice=nslice
!        ELSE
!            WRITE(6,'(3(a,i0),a)')'[ERROR] real_extractXZSlice: Invalid slice requested (',&
!                & nslice,') on ',spec_rank,' should be in [1:',val%nx,']'
!            CALL MPI_FINALIZE(ierr)
!        STOP 'FAILED'
!        ENDIF
!    ELSE
!        curslice=1
!    ENDIF
!
!    slice=0
!
!    IF(val%curlayout .EQ. alongZ) THEN
!        DO i=1,SIZE(YZ_Nodes,dim=1)
!            rankInCW =YZ_Nodes(i)
!            IF (val%X4layoutOnZ(rankInCW +1,1).LE.curslice .AND. &
!            &val%X4layoutOnZ(rankInCW +1,2).GE.curslice) THEN
!                ALLOCATE(buff(val%Y4layoutOnZ(rankInCW +1,1):val%Y4layoutOnZ(rankInCW+1,2),val%nz))
!                IF(spec_rank.EQ.rankInCW) &
!                    & buff=val%values(curslice,val%ymin:val%ymax,1:val%nz)
!                CALL MPI_BCAST( &
!                    & buff, SIZE(buff), MPI_REAL_WP, i-1, YZ_COMMUNICATOR, ierr)
!                IF(ierr.NE.MPI_SUCCESS) THEN
!                    WRITE(6,'(a,i0,a)')"ERROR real_extractYZSlice on process ",spec_rank,&
!                        &": MPIBCAST fails."
!                    CALL MPI_FINALIZE(ierr)
!                    STOP 'FAILED'
!                ENDIF
!
!                slice(val%Y4layoutOnZ(rankInCW +1,1):val%Y4layoutOnZ(rankInCW+1,2),:)=buff
!                DEALLOCATE(buff)
!            ENDIF
!        ENDDO
!    ELSE IF(val%curlayout .EQ. alongY) THEN
!        DO i=1,SIZE(YZ_Nodes,dim=1)
!            rankInCW =YZ_Nodes(i)
!            IF (val%X4layoutOnY(rankInCW +1,1).LE.curslice .AND. &
!                &val%X4layoutOnY(rankInCW +1,2).GE.curslice) THEN
!                ALLOCATE(buff(val%ny, val%Z4layoutOnY(rankInCW +1,1):val%Z4layoutOnY(rankInCW+1,2)))
!                IF(spec_rank.EQ.rankInCW) buff=val%values(curslice,1:val%ny,val%zmin:val%zmax)
!
!                CALL MPI_BCAST( &
!                    & buff, SIZE(buff), MPI_REAL_WP, i-1, YZ_COMMUNICATOR, ierr)
!                IF(ierr.NE.MPI_SUCCESS) THEN
!                    WRITE(6,'(a,i0,a)')"ERROR real_extractYZSlice on process ",spec_rank,&
!                        &": MPIBCAST fails."
!                    CALL MPI_FINALIZE(ierr)
!                    STOP 'FAILED'
!                ENDIF
!
!                slice(:,val%Z4layoutOnY(rankInCW +1,1):val%Z4layoutOnY(rankInCW+1,2))=buff
!                DEALLOCATE(buff)
!            ENDIF
!        ENDDO
!    ELSE IF(val%curlayout .EQ. alongX) THEN 
!        ! IF(spec_rank.EQ.0) WRITE(6,'(a)') &
!        !&'[WARNING] Calling real_extractXZSlice when layout along X is a very expensive choice'
!        DO i=1,SIZE(val%Y4layoutOnX,dim=1)
!            ALLOCATE(buff(val%Y4layoutOnX(i,1):val%Y4layoutOnX(i,2),&
!                &  val%Z4layoutOnX(i,1):val%Z4layoutOnX(i,2)))
!            IF(spec_rank.EQ.i-1) buff=val%values(curslice,val%ymin:val%ymax,val%zmin:val%zmax)
!
!            CALL MPI_BCAST( &
!                & buff, SIZE(buff), MPI_REAL_WP, i-1, spec_communicator, ierr)
!            IF(ierr.NE.MPI_SUCCESS) THEN
!                WRITE(6,'(a,i0,a)')"ERROR real_extractYZSlice on process ",spec_rank,&
!                    &": MPIBCAST fails."
!                CALL MPI_FINALIZE(ierr)
!                STOP 'FAILED'
!            ENDIF
!
!            slice(val%Y4layoutOnX(i,1):val%Y4layoutOnX(i,2),&
!                & val%Z4layoutOnX(i,1):val%Z4layoutOnX(i,2))=buff
!            DEALLOCATE(buff)
!        ENDDO
!    ELSE
!        WRITE(6,'(2(a,i0))')'[ERROR] real_extractXZSlice: Unknown layout for the data on ',&
!                       & spec_rank,' : ',val%curlayout
!        CALL MPI_FINALIZE(ierr)
!        STOP 'FAILED'
!    ENDIF
!
!    RETURN
!END FUNCTION real_extractYZSlice
!
!!------------------------------------------------------------------------------
!!> @author
!!> Patrick BEGOU
!!
!!> @details
!!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!!> data by hand. Arguments are similar to MPI_SENDRECV call in MPI subroutine.
!!
!!> @param[in] senddata a Fortran 3D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!!> @param[in] sendCount the number of item to send
!!> @param[in] sendType  the type of item to send
!!> @param[in] sendTo the rank of the process we send to.
!!> @param[in] sendTag the tag of the send communication
!!> @param[in,out] recvdata a Fortran 3D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!!> @param[in] recvCount the number of item to  receive.
!!> @param[in] recvType the type of item we receive.
!!> @param[in] recvType the type of item we receive.
!!> @param[in] recvFrom the rank of the process we receive from.
!!> @param[in] recvTag the tag of the receive communication
!!> @param[in] com the current communicator
!!> @param[in,out] status the MPI status of this message
!!> @param[in,out] res is set to MPI_SUCCESS if all is OK. Keeped for MPI_SENDRRECV analogy.
!!> @param[in] me the rank of this process in the current communicator
!!> @return  .TRUE. if data exchange is successfull.
!!------------------------------------------------------------------------------
!
!FUNCTION real_MPI_SENDRECV(senddata,sendCount,sendType,sendTo,sendTag,recvdata,recvCount,recvType,&
!                          & recvFrom,recvTag,com,status,res,me)
!
!  USE precision_tools
!  USE mpi
!  IMPLICIT NONE
!
!    INTEGER, INTENT(IN) :: sendCount,sendType,sendTo,sendTag
!    INTEGER, INTENT(IN) :: recvCount,recvType,recvFrom,recvTag
!    INTEGER, INTENT(IN) :: com, me
!    INTEGER, INTENT(INOUT) :: res, status(MPI_STATUS_SIZE)
!    REAL(WP), DIMENSION(:,:,:), INTENT(IN)    :: senddata
!    REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: recvdata
!    LOGICAL :: real_MPI_SENDRECV
!
!    REAL(WP), DIMENSION(:,:,:), ALLOCATABLE:: sendbuf,recvbuf
!    INTEGER::ires
!
!    real_MPI_SENDRECV=.FALSE.
!
!    ALLOCATE(sendbuf(SIZE(senddata,1),SIZE(senddata,2),SIZE(senddata,3)), &
!           & recvbuf(SIZE(recvdata,1),SIZE(recvdata,2),SIZE(recvdata,3)), &
!           & stat=ires)
!    IF (ires .NE.0) THEN
!        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
!                            &in real_MPI_SENDRECV!"
!        RETURN
!    ENDIF
!
!
!    sendbuf=senddata
!    recvbuf=0.0_WP
!
!    CALL MPI_SENDRECV(sendBuf,sendCount,sendType,sendTo,sendTag,recvBuf,recvCount,recvType,recvFrom,&
!                     &recvTag,com,status,res)
!
!    IF (res.NE.MPI_SUCCESS) THEN
!        WRITE(6,'(a,3(i0,a))')"[ERROR] on process ",me,": send/receive failure between",sendTo," and &
!                               &",recvFrom," in real_MPI_SENDRECV."
!        RETURN
!    ENDIF
!
!    recvdata=recvbuf
!    DEALLOCATE(sendbuf,recvbuf)
!    real_MPI_SENDRECV=.TRUE.
!    RETURN
!END FUNCTION real_MPI_SENDRECV


END MODULE realVectorparallel
