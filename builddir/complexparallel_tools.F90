!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!USEFORTEST io
!------------------------------------------------------------------------------
! MODULE cmplxparallel
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
MODULE cmplxparallel


    IMPLICIT NONE


    !private subroutines & functions
    PRIVATE cmplx_AlongXtoAlongY
    PRIVATE cmplx_AlongYtoAlongX
    PRIVATE cmplx_AlongYtoAlongZ
    PRIVATE cmplx_AlongZtoAlongY


CONTAINS
  
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the X axis and distributed on Y and Z among the processes.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION cmplx_setAlongX(spec_rank,val)
    USE datalayout
    USE communicators_tools
    IMPLICIT NONE
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)             :: spec_rank
    LOGICAL                         :: cmplx_setAlongX

    cmplx_setAlongX=.FALSE.

    IF (.NOT. CommInitialized()) THEN
       WRITE(6,'(a)')"Error cmplx_setAlongX: Call initCommunicators first!"
       RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        cmplx_setAlongX=.TRUE.
        RETURN
    CASE (alongY)
        cmplx_setAlongX=cmplx_AlongYtoAlongX(spec_rank,val) 
        RETURN
    CASE (alongZ)
        cmplx_setAlongX=cmplx_AlongZtoAlongY(spec_rank,val) 
        IF (cmplx_setAlongX) cmplx_setAlongX=cmplx_AlongYtoAlongX(spec_rank,val)
        RETURN
    END SELECT
    RETURN

END FUNCTION cmplx_setAlongX

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the Y axis and distributed on X and Z among the processes.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION cmplx_setAlongY(spec_rank,val)
    USE datalayout
    USE communicators_tools
    IMPLICIT NONE
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)             :: spec_rank
    LOGICAL                         :: cmplx_setAlongY

    cmplx_setAlongY=.FALSE.

    IF (.NOT. CommInitialized()) THEN
        WRITE(6,'(a)')"Error cmplx_setAlongY: Call initCommunicators first!"
        RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        cmplx_setAlongY=cmplx_AlongXtoAlongY(spec_rank,val) 
        RETURN
    CASE (alongY)
        cmplx_setAlongY=.TRUE.
        RETURN
    CASE (alongZ)
        cmplx_setAlongY=cmplx_AlongZtoAlongY(spec_rank,val) 
        RETURN
    END SELECT
    RETURN

END FUNCTION cmplx_setAlongY


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are set
!> contiguous along the Y axis and distributed on X and Z among the processes.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
FUNCTION cmplx_setAlongZ(spec_rank,val)
    USE datalayout
    USE communicators_tools
    IMPLICIT NONE
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)             :: spec_rank
    LOGICAL                         :: cmplx_setAlongZ

    cmplx_setAlongZ=.FALSE.

    IF (.NOT. CommInitialized()) THEN
        WRITE(6,'(a)')"Error cmplx_setAlongZ: Call initCommunicators first!"
        RETURN
    ENDIF

    SELECT CASE (val%curlayout)
    CASE (alongX)  
        cmplx_setAlongZ=cmplx_AlongXtoAlongY(spec_rank,val) 
        IF (cmplx_setAlongZ) cmplx_setAlongZ=cmplx_AlongYtoAlongZ(spec_rank,val)
        RETURN
    CASE (alongY)
        cmplx_setAlongZ=cmplx_AlongYtoAlongZ(spec_rank,val)
        RETURN
    CASE (alongZ)
        cmplx_setAlongZ=.TRUE.
        RETURN
    END SELECT
    RETURN

END FUNCTION cmplx_setAlongZ

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along X and are reorganized to be contiguous along Y 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
!------------------------------------------------------------------------------
  FUNCTION cmplx_AlongXtoAlongY(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             ::spec_rank
  LOGICAL                         ::cmplx_AlongXtoAlongY
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  COMPLEX(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  cmplx_AlongXtoAlongY=.FALSE.
  !The idea is to create two arrays containing the dimensions of 
  !  - the subset of data to send to each processor (whatIsend(0:ncpus,6)
  !  - the subset of data I receive from each processor (whatIget(0:ncpus-1,6)
  ! The index correspond to the rank in the communicator
  ALLOCATE(whatIsend(0:SIZE(XY_Nodes)-1,6),whatIget(0:SIZE(XY_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongXtoAlongY on process ",spec_rank,": not enought memory for internal arrays!"
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
       WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongXtoAlongY on process ",spec_rank,": not enought spec_rankmory for datas!"
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
       IF(.NOT. cmplx_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_COMPLEX_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_COMPLEX_WP,i,MPI_ANY_TAG,XY_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  DEALLOCATE(val%values)
  val%values=>values
  values=>null()
  val%xmin=val%X4layoutOnY(spec_rank+1,1)
  val%xmax=val%X4layoutOnY(spec_rank+1,2)
  val%ymin=1
  val%ymax=val%ny
  val%zmin=val%Z4layoutOnY(spec_rank+1,1)
  val%zmax=val%Z4layoutOnY(spec_rank+1,2)
  val%curlayout=alongY
  
  DEALLOCATE(whatIsend,whatIget)
  
  cmplx_AlongXtoAlongY=.TRUE.
  RETURN
  
  END FUNCTION cmplx_AlongXtoAlongY
  
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Y and are reorganized to be contiguous along X 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communitor
! See cmplx_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION cmplx_AlongYtoAlongX(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: cmplx_AlongYtoAlongX
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  COMPLEX(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  cmplx_AlongYtoAlongX=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(XY_Nodes)-1,6),whatIget(0:SIZE(XY_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongYtoAlongX on process ",spec_rank,": not enought memory for internal arrays!"
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
       WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongYtoAlongX on process ",spec_rank,": not enought memory for datas!"
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

       IF (.NOT. cmplx_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_COMPLEX_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_COMPLEX_WP,i,MPI_ANY_TAG,XY_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  DEALLOCATE(val%values)
  val%values=>values
  values=>null()
  val%xmin=1
  val%xmax=val%nx
  val%ymin=val%Y4layoutOnX(spec_rank+1,1)
  val%ymax=val%Y4layoutOnX(spec_rank+1,2)
  val%zmin=val%Z4layoutOnX(spec_rank+1,1)
  val%zmax=val%Z4layoutOnX(spec_rank+1,2)
  val%curlayout=alongX
  
  DEALLOCATE(whatIsend,whatIget)
  
  cmplx_AlongYtoAlongX=.TRUE.
  RETURN
  
  END FUNCTION cmplx_AlongYtoAlongX

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Y and are reorganized to be contiguous along Z 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communicator
! See cmplx_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION cmplx_AlongYtoAlongZ(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: cmplx_AlongYtoAlongZ
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  COMPLEX(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  cmplx_AlongYtoAlongZ=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(YZ_Nodes)-1,6),whatIget(0:SIZE(YZ_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongYtoAlongZ on process ",spec_rank,": not enought memory for internal arrays!"
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
       WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongYtoAlongZ on process ",spec_rank,": not enought memory for datas!"
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

       IF (.NOT. cmplx_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_COMPLEX_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_COMPLEX_WP,i,MPI_ANY_TAG,YZ_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  DEALLOCATE(val%values)
  val%values=>values
  values=>null()
  val%xmin=val%X4layoutOnZ(spec_rank+1,1)
  val%xmax=val%X4layoutOnZ(spec_rank+1,2)
  val%ymin=val%Y4layoutOnZ(spec_rank+1,1)
  val%ymax=val%Y4layoutOnZ(spec_rank+1,2)
  val%zmin=1
  val%zmax=val%nz
  val%curlayout=alongZ
  
  DEALLOCATE(whatIsend,whatIget)
  
  cmplx_AlongYtoAlongZ=.TRUE.
  RETURN
  
  END FUNCTION cmplx_AlongYtoAlongZ


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine changes the shape of the COMPLEX_DATA_LAYOUT variable: val. Datas are 
!> supposed to be contiguous along Z and are reorganized to be contiguous along Y 
!> on the current processor. This version uses specific datatype.
!> @param[in,out] val the COMPLEX_DATA_LAYOUT variable to redistribute.
!> @param[in] spec_rank my rank in spec_communicator
! See cmplx_AlongXtoAlongY() for more details and comments.
!------------------------------------------------------------------------------
  FUNCTION cmplx_AlongZtoAlongY(spec_rank,val)
  USE datalayout
  USE communicators_tools
  USE mpi
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
  INTEGER, INTENT(IN)             :: spec_rank
  LOGICAL                         :: cmplx_AlongZtoAlongY
  

  INTEGER                            :: res,i,cptin,cptout
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: whatIsend,whatIget
  COMPLEX(WP), DIMENSION(:,:,:),POINTER :: values => NULL()     
  INTEGER :: status(MPI_STATUS_SIZE)
  
  cmplx_AlongZtoAlongY=.FALSE.

  ALLOCATE(whatIsend(0:SIZE(YZ_Nodes)-1,6),whatIget(0:SIZE(YZ_Nodes)-1,6),stat=res)
  IF(res.NE.0) THEN
     WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongZtoAlongY on process ",spec_rank,": not enought memory for internal arrays!"
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
       WRITE(6,'(a,i0,a)')"ERROR cmplx_AlongZtoAlongY on process ",spec_rank,": not enought memory for datas!"
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

       IF (.NOT. cmplx_MPI_SENDRECV(&
       & val%values(whatIsend(i,1):whatIsend(i,2),whatIsend(i,3):whatIsend(i,4),whatIsend(i,5):whatIsend(i,6))&
       &,cptout,MPI_COMPLEX_WP,i,i &
       &,values(whatIget(i,1):whatIget(i,2),whatIget(i,3):whatIget(i,4),whatIget(i,5):whatIget(i,6)) &
       &,cptin,MPI_COMPLEX_WP,i,MPI_ANY_TAG,YZ_COMMUNICATOR,status,res,spec_rank)) RETURN
    ENDIF
  ENDDO

  !Update val structure pointing to this new array
  DEALLOCATE(val%values)
  val%values=>values
  values=>null()
  val%xmin=val%X4layoutOnY(spec_rank+1,1)
  val%xmax=val%X4layoutOnY(spec_rank+1,2)
  val%ymin=1
  val%ymax=val%ny
  val%zmin=val%Z4layoutOnY(spec_rank+1,1)
  val%zmax=val%Z4layoutOnY(spec_rank+1,2)
  val%curlayout=alongY
  
  DEALLOCATE(whatIsend,whatIget)
  
  cmplx_AlongZtoAlongY=.TRUE.
  RETURN
  
  END FUNCTION cmplx_AlongZtoAlongY
  
  
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!
!> @details
!> This function makes a local sommation of the variable "value" on returns the
!> result on node 0.
! 
!> @param[in] value the local variable to sum.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @param[in] node the rank number of the root process in the range [0,ncpus - 1]
!>            by default this value is at 0 (optional).
!> @return the sommation of all the value variable of each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoLocalSum(spec_rank,value,node)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  
  COMPLEX(WP) ::cmplx_DoLocalSum 
  COMPLEX(WP), INTENT(IN)          :: value
  INTEGER, INTENT(IN)           :: spec_rank
  INTEGER, INTENT(IN),OPTIONAL  :: node

  COMPLEX(WP) :: resultat
  INTEGER  :: ierr
 
  resultat = 0.0_WP
  IF ( PRESENT(node) ) THEN
    CALL MPI_REDUCE(value,resultat,1,MPI_COMPLEX_WP,MPI_SUM,node,spec_communicator,ierr)
  ELSE
    CALL MPI_REDUCE(value,resultat,1,MPI_COMPLEX_WP,MPI_SUM,0,spec_communicator,ierr)
  ENDIF
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoLocalSum on process ",spec_rank,&
                      &": Unable to run local sum."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  cmplx_DoLocalSum=resultat
  RETURN
  
END FUNCTION cmplx_DoLocalSum

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function makes a global sommation of the variable "value" on returns the
!> result on every node.
! 
!> @param[in] value the local variable to sum.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @return the sommation of all the value variable of each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoGlobalSum(spec_rank,value)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
  
  COMPLEX(WP) ::cmplx_DoGlobalSum 
  COMPLEX(WP), INTENT(IN) :: value
  INTEGER, INTENT(IN)  :: spec_rank

  COMPLEX(WP) :: resultat
  INTEGER  ::ierr
  
  CALL MPI_ALLREDUCE(value,resultat,1,MPI_COMPLEX_WP,MPI_SUM,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoGlobalSum on process ",spec_rank,&
                      &": Unable to run global sum."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  cmplx_DoGlobalSum=resultat
  RETURN
  
END FUNCTION cmplx_DoGlobalSum


!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function makes a global sommation of the 1D array "vector" on returns the
!> result on every node.
! 
!> @param[in] vector the local array to sum.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @return the sommation of all the value of vector for each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoGlobalVectorSum(spec_rank,vector)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  
  COMPLEX(WP), INTENT(IN) :: vector(:)
  COMPLEX(WP) ::cmplx_DoGlobalVectorSum (size(vector))
  INTEGER, INTENT(IN)  :: spec_rank

!  COMPLEX(WP) :: resultat
  INTEGER  ::ierr
  
  CALL MPI_ALLREDUCE(vector,cmplx_DoGlobalVectorSum,size(vector),MPI_COMPLEX_WP,MPI_SUM,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoGlobalVectorSum on process ",spec_rank,&
                      &": Unable to run global sum on the array."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
!  cmplx_DoGlobalSum=resultat
  RETURN
  
END FUNCTION cmplx_DoGlobalVectorSum


!------------------------------------------------------------------------------
!> @author
!> Jean-Baptiste Lagaert, Antoine Vollant
!
!> @details
!> This function makes a global sommation of the 1D array "vector" on returns the
!> result on nodes 0 (default) or in the node given in input.
!
!> @param[in] array2D the local 2D-array to sum.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @param[in] node the rank number of the root process in the range [0,ncpus - 1]
!>            by default this value is at 0 (optional).
!> @return the sommation of all the value of vector for each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoLocalArray2DSum(spec_rank,array2D,node) result(resultat)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi

  IMPLICIT NONE

  COMPLEX(WP), DIMENSION(:,:), INTENT(IN) :: array2D
  INTEGER, INTENT(IN)  :: spec_rank
  INTEGER, INTENT(IN), OPTIONAL :: node
  COMPLEX(WP), DIMENSION(size(array2D,1), size(array2D,2)) :: resultat
  INTEGER  ::ierr

  resultat = 1.0_WP
  IF ( PRESENT(node) ) THEN
    CALL MPI_REDUCE(array2D,resultat,size(array2D),MPI_COMPLEX_WP,MPI_SUM,node,spec_communicator,ierr)
  ELSE
    CALL MPI_REDUCE(array2D,resultat,size(array2D),MPI_COMPLEX_WP,MPI_SUM,0,spec_communicator,ierr)
  ENDIF
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoLocalArray2DSum on process ",spec_rank,&
                      &": Unable to run local sum."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF

END FUNCTION cmplx_DoLocalArray2DSum

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function makes a global max of the variable "value" and returns the
!> result on every node.
! 
!> @param[in] value the local variable to test.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @return the maximum of all the value variable over each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoGlobalMax(spec_rank,value)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
  
  COMPLEX(WP) ::cmplx_DoGlobalMax 
  COMPLEX(WP), INTENT(IN) :: value
  INTEGER, INTENT(IN)  :: spec_rank

  COMPLEX(WP) :: resultat
  INTEGER  ::ierr
  
  CALL MPI_ALLREDUCE(value,resultat,1,MPI_COMPLEX_WP,MPI_MAX,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoGlobalMax on process ",spec_rank,&
                      &": Unable to run global Max."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  cmplx_DoGlobalMax=resultat
  RETURN
  
END FUNCTION cmplx_DoGlobalMax

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function makes a global max of the variable "value" and returns the
!> result on every node.
! 
!> @param[in] value the local variable to test.
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @return the maximum of all the value variable over each node.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoGlobalMin(spec_rank,value)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
  
  COMPLEX(WP) ::cmplx_DoGlobalMin 
  COMPLEX(WP), INTENT(IN) :: value
  INTEGER, INTENT(IN)  :: spec_rank

  COMPLEX(WP) :: resultat
  INTEGER  ::ierr
  
  CALL MPI_ALLREDUCE(value,resultat,1,MPI_COMPLEX_WP,MPI_MIN,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR cmplx_DoGlobalMin on process ",spec_rank,&
                      &": Unable to run global Min."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  cmplx_DoGlobalMin=resultat
  RETURN
  
END FUNCTION cmplx_DoGlobalMin

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This function makes a global sommation of the 1D array "vector" on returns the
!> result on every node.
! 
!> @param[in] val a COMPLEX_DATA_LAYOUT from wich we want to extract a slice in the [Y,Z] plane
!> @param[in] spec_rank the rank number of the process in the range [0,ncpus - 1].
!> @param[in] nslice the slice position. defauklt value is 1 if not provided.
!> @return a 2D [ny,nz] array containing the requested slice or exit.
!------------------------------------------------------------------------------
FUNCTION cmplx_extractYZSlice(val,spec_rank,nslice) RESULT(slice)
!------------------------------------------------------------------------------
    USE precision_tools
    USE datalayout
    USE communicators_tools
    USE mpi

    IMPLICIT NONE
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN) :: val
    INTEGER, INTENT(IN)               :: spec_rank
    INTEGER, INTENT(IN), OPTIONAL     :: nslice
    COMPLEX(WP), DIMENSION(:,:), ALLOCATABLE::buff
    COMPLEX(WP), DIMENSION(val%ny,val%nz):: slice

    INTEGER:: curslice, i, ierr
    !> Rank in spec_communicator when i is the rank in  YZ_COMMUNICATOR
    INTEGER::rankInCW

    IF(PRESENT(nslice)) THEN
        IF(nslice.GE.1 .AND. nslice.LE.val%nx) THEN
            curslice=nslice
        ELSE
            WRITE(6,'(3(a,i0),a)')'[ERROR] cmplx_extractXZSlice: Invalid slice requested (',&
                & nslice,') on ',spec_rank,' should be in [1:',val%nx,']'
            CALL MPI_FINALIZE(ierr)
        STOP 'FAILED'
        ENDIF
    ELSE
        curslice=1
    ENDIF

    slice=0

    IF(val%curlayout .EQ. alongZ) THEN
        DO i=1,SIZE(YZ_Nodes,dim=1)
            rankInCW =YZ_Nodes(i)
            IF (val%X4layoutOnZ(rankInCW +1,1).LE.curslice .AND. &
            &val%X4layoutOnZ(rankInCW +1,2).GE.curslice) THEN
                ALLOCATE(buff(val%Y4layoutOnZ(rankInCW +1,1):val%Y4layoutOnZ(rankInCW+1,2),val%nz))
                IF(spec_rank.EQ.rankInCW) &
                    & buff=val%values(curslice,val%ymin:val%ymax,1:val%nz)
                CALL MPI_BCAST( &
                    & buff, SIZE(buff), MPI_COMPLEX_WP, i-1, YZ_COMMUNICATOR, ierr)
                IF(ierr.NE.MPI_SUCCESS) THEN
                    WRITE(6,'(a,i0,a)')"ERROR cmplx_extractYZSlice on process ",spec_rank,&
                        &": MPIBCAST fails."
                    CALL MPI_FINALIZE(ierr)
                    STOP 'FAILED'
                ENDIF

                slice(val%Y4layoutOnZ(rankInCW +1,1):val%Y4layoutOnZ(rankInCW+1,2),:)=buff
                DEALLOCATE(buff)
            ENDIF
        ENDDO
    ELSE IF(val%curlayout .EQ. alongY) THEN
        DO i=1,SIZE(YZ_Nodes,dim=1)
            rankInCW =YZ_Nodes(i)
            IF (val%X4layoutOnY(rankInCW +1,1).LE.curslice .AND. &
                &val%X4layoutOnY(rankInCW +1,2).GE.curslice) THEN
                ALLOCATE(buff(val%ny, val%Z4layoutOnY(rankInCW +1,1):val%Z4layoutOnY(rankInCW+1,2)))
                IF(spec_rank.EQ.rankInCW) buff=val%values(curslice,1:val%ny,val%zmin:val%zmax)

                CALL MPI_BCAST( &
                    & buff, SIZE(buff), MPI_COMPLEX_WP, i-1, YZ_COMMUNICATOR, ierr)
                IF(ierr.NE.MPI_SUCCESS) THEN
                    WRITE(6,'(a,i0,a)')"ERROR cmplx_extractYZSlice on process ",spec_rank,&
                        &": MPIBCAST fails."
                    CALL MPI_FINALIZE(ierr)
                    STOP 'FAILED'
                ENDIF

                slice(:,val%Z4layoutOnY(rankInCW +1,1):val%Z4layoutOnY(rankInCW+1,2))=buff
                DEALLOCATE(buff)
            ENDIF
        ENDDO
    ELSE IF(val%curlayout .EQ. alongX) THEN 
        ! IF(spec_rank.EQ.0) WRITE(6,'(a)') &
        !&'[WARNING] Calling cmplx_extractXZSlice when layout along X is a very expensive choice'
        DO i=1,SIZE(val%Y4layoutOnX,dim=1)
            ALLOCATE(buff(val%Y4layoutOnX(i,1):val%Y4layoutOnX(i,2),&
                &  val%Z4layoutOnX(i,1):val%Z4layoutOnX(i,2)))
            IF(spec_rank.EQ.i-1) buff=val%values(curslice,val%ymin:val%ymax,val%zmin:val%zmax)

            CALL MPI_BCAST( &
                & buff, SIZE(buff), MPI_COMPLEX_WP, i-1, spec_communicator, ierr)
            IF(ierr.NE.MPI_SUCCESS) THEN
                WRITE(6,'(a,i0,a)')"ERROR cmplx_extractYZSlice on process ",spec_rank,&
                    &": MPIBCAST fails."
                CALL MPI_FINALIZE(ierr)
                STOP 'FAILED'
            ENDIF

            slice(val%Y4layoutOnX(i,1):val%Y4layoutOnX(i,2),&
                & val%Z4layoutOnX(i,1):val%Z4layoutOnX(i,2))=buff
            DEALLOCATE(buff)
        ENDDO
    ELSE
        WRITE(6,'(2(a,i0))')'[ERROR] cmplx_extractXZSlice: Unknown layout for the data on ',&
                       & spec_rank,' : ',val%curlayout
        CALL MPI_FINALIZE(ierr)
        STOP 'FAILED'
    ENDIF

    RETURN
END FUNCTION cmplx_extractYZSlice

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_SENDRECV call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 3D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in] sendType  the type of item to send
!> @param[in] sendTo the rank of the process we send to.
!> @param[in] sendTag the tag of the send communication
!> @param[in,out] recvdata a Fortran 3D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive.
!> @param[in] recvType the type of item we receive.
!> @param[in] recvFrom the rank of the process we receive from.
!> @param[in] recvTag the tag of the receive communication
!> @param[in] com the current communicator
!> @param[in,out] status the MPI status of this message
!> @param[in,out] res is set to MPI_SUCCESS if all is OK. Keeped for MPI_SENDRRECV analogy.
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data exchange is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_MPI_SENDRECV(senddata,sendCount,sendType,sendTo,sendTag,recvdata,recvCount,recvType,&
                          & recvFrom,recvTag,com,status,res,me)

  USE precision_tools
  USE mpi
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount,sendType,sendTo,sendTag
    INTEGER, INTENT(IN) :: recvCount,recvType,recvFrom,recvTag
    INTEGER, INTENT(IN) :: com, me
    INTEGER, INTENT(INOUT) :: res, status(MPI_STATUS_SIZE)
    COMPLEX(WP), DIMENSION(:,:,:), INTENT(IN)    :: senddata
    COMPLEX(WP), DIMENSION(:,:,:), INTENT(INOUT) :: recvdata
    LOGICAL :: cmplx_MPI_SENDRECV

    COMPLEX(WP), DIMENSION(:,:,:), ALLOCATABLE:: sendbuf,recvbuf
    INTEGER::ires

    cmplx_MPI_SENDRECV=.FALSE.

    ALLOCATE(sendbuf(SIZE(senddata,1),SIZE(senddata,2),SIZE(senddata,3)), &
           & recvbuf(SIZE(recvdata,1),SIZE(recvdata,2),SIZE(recvdata,3)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_MPI_SENDRECV!"
        RETURN
    ENDIF


    sendbuf=senddata
    recvbuf=0.0_WP

    CALL MPI_SENDRECV(sendbuf,sendCount,sendType,sendTo,sendTag,recvBuf,recvCount,recvType,recvFrom,&
                     &recvTag,com,status,res)

    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,3(i0,a))')"[ERROR] on process ",me,": send/receive failure between",sendTo," and &
                               &",recvFrom," in cmplx_MPI_SENDRECV."
        RETURN
    ENDIF

    recvdata=recvbuf
    DEALLOCATE(sendbuf,recvbuf)
    cmplx_MPI_SENDRECV=.TRUE.
    RETURN
END FUNCTION cmplx_MPI_SENDRECV


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_SENDRECV call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in] sendTo the rank of the process we send to.
!> @param[in] sendTag the tag of the send communication
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive.
!> @param[in] recvFrom the rank of the process we receive from.
!> @param[in] recvTag the tag of the receive communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data exchange is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoSendRecvVector(senddata,sendCount,sendTo,sendTag,recvdata,recvCount,recvFrom,recvTag,me)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount,sendTo,sendTag
    INTEGER, INTENT(IN) :: recvCount,recvFrom,recvTag
    INTEGER, INTENT(IN) :: me
    INTEGER :: status(MPI_STATUS_SIZE)
    COMPLEX(WP), DIMENSION(:), INTENT(IN)    :: senddata
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT) :: recvdata
    LOGICAL :: cmplx_DoSendRecvVector

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: sendbuf,recvbuf
    INTEGER :: res
    INTEGER::ires

    cmplx_DoSendRecvVector=.FALSE.

    ALLOCATE(sendbuf(SIZE(senddata)), &
           & recvbuf(SIZE(recvdata)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoSendRecvVector!"
        RETURN
    ENDIF


    sendbuf=senddata
    recvbuf=0.0_WP

    CALL MPI_SENDRECV(sendbuf,sendCount,MPI_COMPLEX_WP,sendTo,sendTag,recvBuf,recvCount,MPI_COMPLEX_WP,recvFrom,&
                     &recvTag,spec_communicator,status,res)

    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,3(i0,a))')"[ERROR] on process ",me,": send/receive failure between",sendTo," and &
                               &",recvFrom," in cmplx_DoSendRecvVector."
        RETURN
    ENDIF

    recvdata=recvbuf
    DEALLOCATE(sendbuf,recvbuf)

    cmplx_DoSendRecvVector=.TRUE.
    RETURN

END FUNCTION cmplx_DoSendRecvVector

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_SEND call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran value of datas to send.
!> @param[in] sendType  the type of item to send
!> @param[in] sendTo the rank of the process we send to.
!> @param[in] sendTag the tag of the send communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data exchange is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoSendScal(senddata,sendTo,sendTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendTo,sendTag
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP), INTENT(IN)    :: senddata
    LOGICAL :: success 
    INTEGER :: res

    success=.FALSE.

    CALL MPI_SEND(senddata,1,MPI_COMPLEX_WP,sendTo,sendTag,spec_communicator,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoSendScal."
        RETURN
    ENDIF
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoSendScal

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_RECV call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module out of this module
!
!> @param[in,out] recvdata a Fortran value to store the received datas.
!> @param[in] recvType the type of item we receive.
!> @param[in] recvFrom the rank of the process we receive from.
!> @param[in] recvTag the tag of the receive communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data receive is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoRecvScal(recvdata,recvFrom,recvTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: recvFrom,recvTag
    INTEGER, INTENT(IN) :: me
    INTEGER :: status(MPI_STATUS_SIZE)
    COMPLEX(WP), INTENT(INOUT) :: recvdata
    LOGICAL :: success 
    INTEGER :: res

    success=.FALSE.

    CALL MPI_RECV(recvdata,1,MPI_COMPLEX_WP,recvFrom,recvTag,spec_communicator,status,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": receive failure in cmplx_DoRecvScal."
        RETURN
    ENDIF
    success=.TRUE.

    RETURN

END FUNCTION cmplx_DoRecvScal


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_SEND call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in] sendType  the type of item to send
!> @param[in] sendTo the rank of the process we send to.
!> @param[in] sendTag the tag of the send communication
!> @param[in] com the current communicator
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data exchange is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoSendVector(senddata,sendCount,sendTo,sendTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount,sendTo,sendTag
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP), DIMENSION(:), INTENT(IN)    :: senddata
    LOGICAL :: success 

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: sendbuf
    INTEGER :: res
    INTEGER::ires

    success=.FALSE.

    IF ( SIZE(senddata) .LT. sendCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": sendCount disagree with &
                            & SIZE(senddata) in cmplx_DoSendScal!"
        RETURN
    ENDIF
    ALLOCATE(sendbuf(SIZE(senddata)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoSendVector!"
        RETURN
    ENDIF
    sendbuf=senddata
    CALL MPI_SEND(sendbuf,sendCount,MPI_COMPLEX_WP,sendTo,sendTag,spec_communicator,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoSendVector."
        RETURN
    ENDIF
    DEALLOCATE(sendbuf)
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoSendVector

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper because MPI 2 do not knows about Fortran90 subarrays on some compilers (xlf...)
!> This will solve MPI wrong behavior on IBM Regatta with allocating buffer to store contiguous
!> data by hand. Arguments are similar to MPI_RECV call in MPI subroutine.
!> Be careful, there is an interface to use this function out of this module out of this module
!
!> @param[in,out] recvdata a Fortran array to store the received datas. Depending on compilers these data can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive.
!> @param[in] recvType the type of item we receive.
!> @param[in] recvFrom the rank of the process we receive from.
!> @param[in] recvTag the tag of the receive communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data receive is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoRecvVector(recvdata,recvCount,recvFrom,recvTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: recvCount,recvFrom,recvTag
    INTEGER, INTENT(IN) :: me
    INTEGER :: status(MPI_STATUS_SIZE)
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT) :: recvdata
    LOGICAL :: success 

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: recvbuf
    INTEGER :: res
    INTEGER ::ires

    success=.FALSE.

    IF ( SIZE(recvdata) .LT. recvCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": recvCount disagree with &
                            & SIZE(recvdata) in cmplx_DoRecvVector!"
        RETURN
    ENDIF
    ALLOCATE(recvbuf(SIZE(recvdata)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoRecvVector!"
        RETURN
    ENDIF
    recvbuf=0.0_WP
    CALL MPI_RECV(recvBuf,recvCount,MPI_COMPLEX_WP,recvFrom,recvTag,spec_communicator,status,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": receive failure in cmplx_DoRecvVector."
        RETURN
    ENDIF
    recvdata=recvbuf
    DEALLOCATE(recvbuf)
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoRecvVector 


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_SCATTER for COMPLEX values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive.
!> @param[in] root the rank of the process we receive from.
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data spread is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoScatterVector(senddata,sendCount,recvdata,recvCount,&
                          & root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount
    INTEGER, INTENT(IN) :: recvCount,root
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP), DIMENSION(:), INTENT(IN)    :: senddata
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT) :: recvdata
    LOGICAL :: success 

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: sendbuf,recvbuf
    INTEGER :: res
    INTEGER :: ires

    success=.FALSE.

    IF ( SIZE(senddata) .LT. sendCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": sendCount disagree with &
                            & SIZE(senddata) in cmplx_DoScatterVector!"
        RETURN
    ENDIF
    IF ( SIZE(recvdata) .LT. recvCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": recvCount disagree with &
                            & SIZE(recvdata) in cmplx_DoScatterVector!"
        RETURN
    ENDIF

    ALLOCATE(sendbuf(sendCount), &
           & recvbuf(recvCount), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoScatterVector!"
        RETURN
    ENDIF


    sendbuf=senddata
    recvbuf=0.0_WP

    CALL MPI_SCATTER(sendbuf,sendCount,MPI_COMPLEX_WP,recvbuf,recvCount,MPI_COMPLEX_WP,root,&
                    &spec_communicator,res)

    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoScatterVector."
        RETURN
    ENDIF
    recvdata=recvbuf
    DEALLOCATE(sendbuf,recvbuf)
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoScatterVector

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_GATHER for COMPLEX values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive equal to the number of data received per processus time the number of processus
!> @param[in] root the rank of the process which receive all datas 
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data reassemble is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoGatherVector(senddata,sendCount,recvdata,recvCount,&
                          & root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  USE mpilayout_tools , only : getnbcpus
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount
    INTEGER, INTENT(IN) :: recvCount,root
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP), DIMENSION(:), INTENT(IN)    :: senddata
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT) :: recvdata
    LOGICAL :: success 

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: sendbuf,recvbuf
    INTEGER :: ires,nbproc

    success=.FALSE.
    nbproc = getnbcpus()

    IF ( SIZE(senddata) .LT. sendCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": sendCount disagree with &
                            & SIZE(senddata) in cmplx_DoScatterVector!"
        RETURN
    ENDIF
    IF ( SIZE(recvdata) .LT. nbproc*recvCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": recvCount disagree with &
                            & SIZE(recvdata) in cmplx_DoScatterVector!"
        RETURN
    ENDIF

    ALLOCATE(sendbuf(sendCount), &
           & recvbuf(size(recvdata,1)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoGatherVector!"
        RETURN
    ENDIF

    sendbuf=senddata
    recvbuf=0.0_WP

    CALL MPI_GATHER(sendbuf,sendCount,MPI_COMPLEX_WP,recvbuf,recvCount,MPI_COMPLEX_WP,root,&
                    &spec_communicator,ires)

    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoScatterVector."
        RETURN
    ENDIF
    recvdata=recvbuf
    DEALLOCATE(sendbuf,recvbuf)

    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoGatherVector

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_GATHER for COMPLEX values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] sendCount the number of item to send
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] recvCount the number of item to  receive equal to the number of data received per processus time the number of processus
!> @param[in] root the rank of the process which receive all datas 
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data reassemble is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoBcastVector(sendrecvdata,sendrecvCount,root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  USE mpilayout_tools , only : getnbcpus
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendrecvCount
    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT) :: sendrecvdata
    LOGICAL :: success 

    COMPLEX(WP), DIMENSION(:), ALLOCATABLE:: sendrecvbuf
    INTEGER :: ires,nbproc

    success=.FALSE.
    nbproc = getnbcpus()

    IF ( SIZE(sendrecvdata) .LT. sendrecvCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": sendrecvCount disagree with &
                            & SIZE(sendrecvdata) in cmplx_DoScatterVector!"
        RETURN
    ENDIF

    ALLOCATE(sendrecvbuf(sendrecvCount), stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in cmplx_DoBcastVector!"
        RETURN
    ENDIF
    IF ( me .EQ. root) THEN
        sendrecvbuf=sendrecvdata
    ELSE
        sendrecvbuf=0.0_WP
    ENDIF

    CALL MPI_BCAST(sendrecvbuf,sendrecvCount,MPI_COMPLEX_WP,root,spec_communicator,ires)

    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoScatterVector."
        RETURN
    ENDIF
    sendrecvdata=sendrecvbuf
    DEALLOCATE(sendrecvbuf)

    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoBcastVector

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_BCAST for COMPLEX values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in,out] senrecvddata a Fortran 1D array of datas to send. Depending on compilers these data can have a non contiguous storage.
!> @param[in] root the rank of the process which receive all datas 
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data reassemble is successfull.
!------------------------------------------------------------------------------

FUNCTION cmplx_DoBcastScal(sendrecvdata,root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(IN) :: me
    COMPLEX(WP),INTENT(INOUT) :: sendrecvdata
    LOGICAL             :: success

    INTEGER :: ires

    success=.FALSE.

    CALL MPI_BCAST(sendrecvdata,1,MPI_COMPLEX_WP,root,spec_communicator,ires)

    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoScatterScalar."
        RETURN
    ENDIF
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoBcastScal

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_ALLGATHER for COMPLEX values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a COMPLEX value
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data reassemble is successfull.
!------------------------------------------------------------------------------
FUNCTION cmplx_DoAllGatherScal(senddata,recvdata,me) result(success)

    USE precision_tools
    USE mpi
    USE communicators_tools
    USE mpilayout_tools , only : getnbcpus
    IMPLICIT NONE

    COMPLEX(WP)              , INTENT(IN)       :: senddata
    COMPLEX(WP), DIMENSION(:), INTENT(INOUT)    :: recvdata
    INTEGER, INTENT(IN) :: me
    LOGICAL :: success 
    INTEGER :: ires

    success=.FALSE.

    IF ( SIZE(recvdata) .LT. getnbcpus()) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": recvCount disagree with &
                            & the number of processors in cmplx_DoAllGatherScal!"
        RETURN
    ENDIF
    CALL MPI_ALLGATHER (senddata,1,MPI_COMPLEX_WP,recvdata,1,MPI_COMPLEX_WP, spec_communicator ,ires)
    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in cmplx_DoAllGatherScal."
        RETURN
    ENDIF
    success=.TRUE.
    RETURN

END FUNCTION cmplx_DoAllGatherScal

END MODULE cmplxparallel


