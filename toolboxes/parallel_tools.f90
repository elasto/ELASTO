!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST avgcond 
!USEFORTEST io
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to realize all parallel communications for
!> changing data layout, for calculations and for parallel I/O
!> This is the generic interface of the underlying modules. The hierarchical organisation is shown in the 
!> figure. The generic implementation is located in the file implementparallel.Fortran.
!> It is processed to generate realparallel.F90 (module realparallel for REAL datas) and
!> cmplxparallel.F90 (module cmplxparallel for COMPLEX data).
!>\image html parallel.png "Hierarchical organisation of the parallel modules and files"
!------------------------------------------------------------------------------
MODULE parallel_tools

USE cmplxparallel
USE realparallel
USE cmplxVectorparallel
USE realVectorparallel

IMPLICIT NONE

PRIVATE real_setAlongX
PRIVATE cmplx_setAlongX
PRIVATE real_setAlongY
PRIVATE cmplx_setAlongY
PRIVATE real_setAlongZ
PRIVATE cmplx_setAlongZ

INTERFACE setAlongX
   module procedure real_setAlongX
   module procedure cmplx_setAlongX
   module procedure real_setAlongX_Vector
   module procedure cmplx_setAlongX_Vector
END INTERFACE setAlongX
  
INTERFACE setAlongY
   module procedure real_setAlongY
   module procedure cmplx_setAlongY
   module procedure real_setAlongY_Vector
   module procedure cmplx_setAlongY_Vector
END INTERFACE setAlongY
  
INTERFACE setAlongZ
   module procedure real_setAlongZ
   module procedure cmplx_setAlongZ
   module procedure real_setAlongZ_Vector
   module procedure cmplx_setAlongZ_Vector
END INTERFACE setAlongZ
 
INTERFACE DoLocalSum
   module procedure real_DoLocalSum
   module procedure cmplx_DoLocalSum
END INTERFACE DoLocalSum

INTERFACE DoGlobalSum
   module procedure real_DoGlobalSum
   module procedure cmplx_DoGlobalSum
   module procedure integer_DoGlobalSum 
   module procedure biginteger_DoGlobalSum 
END INTERFACE DoGlobalSum
 
INTERFACE DoGlobalVectorSum
   module procedure real_DoGlobalVectorSum
   module procedure cmplx_DoGlobalVectorSum
   module procedure integer_DoGlobalVectorSum
END INTERFACE DoGlobalVectorSum

INTERFACE DoLocalArray2DSum
    module procedure real_DoLocalArray2DSum, cmplx_DoLocalArray2DSum
END INTERFACE DoLocalArray2DSum

INTERFACE DoGlobalMax
   module procedure real_DoGlobalMax
   module procedure cmplx_DoGlobalMax
END INTERFACE DoGlobalMax

INTERFACE DoGlobalMin
   module procedure real_DoGlobalMin
   module procedure cmplx_DoGlobalMin
END INTERFACE DoGlobalMin
 
INTERFACE extractYZSlice
   module procedure real_extractYZSlice
   module procedure cmplx_extractYZSlice
END INTERFACE extractYZSlice

INTERFACE DoSendRecvVector
   MODULE PROCEDURE real_DoSendRecvVector
   MODULE PROCEDURE cmplx_DoSendRecvVector
   MODULE PROCEDURE integer_DoSendRecvVector
END INTERFACE DoSendRecvVector

INTERFACE DoSendScal
   MODULE PROCEDURE cmplx_DoSendScal
   MODULE PROCEDURE integer_DoSendScal
   MODULE PROCEDURE real_DoSendScal
END INTERFACE DoSendScal

INTERFACE DoRecvScal
   MODULE PROCEDURE real_DoRecvScal
   MODULE PROCEDURE integer_DoRecvScal
   MODULE PROCEDURE cmplx_DoRecvScal
END INTERFACE DoRecvScal

INTERFACE DoSendVector
   MODULE PROCEDURE cmplx_DoSendVector
   MODULE PROCEDURE real_DoSendVector
END INTERFACE DoSendVector

INTERFACE DoRecvVector
   MODULE PROCEDURE real_DoRecvVector
   MODULE PROCEDURE cmplx_DoRecvVector
END INTERFACE DoRecvVector

INTERFACE DoScatterVector
   MODULE PROCEDURE real_DoScatterVector
   MODULE PROCEDURE cmplx_DoScatterVector
END INTERFACE DoScatterVector

INTERFACE DoGatherVector
   MODULE PROCEDURE real_DoGatherVector
   MODULE PROCEDURE cmplx_DoGatherVector
END INTERFACE DoGatherVector

INTERFACE DoBcastVector
   MODULE PROCEDURE real_DoBcastVector
   MODULE PROCEDURE cmplx_DoBcastVector
   MODULE PROCEDURE integer_DoBcastVector
END INTERFACE DoBcastVector

INTERFACE DoBcastScal
   MODULE PROCEDURE real_DoBcastScal
   MODULE PROCEDURE cmplx_DoBcastScal
   MODULE PROCEDURE integer_DoBcastScal
END INTERFACE DoBcastScal

INTERFACE DoAllGatherScal
   MODULE PROCEDURE real_DoAllGatherScal
   MODULE PROCEDURE cmplx_DoAllGatherScal
   MODULE PROCEDURE integer_DoAllGatherScal
END INTERFACE DoAllGatherScal

PUBLIC :: parallel_error

CONTAINS

!> Printing the error message associated to a mpi-error.
SUBROUTINE parallel_error(code, procedure_name)

    !> @ author = Jean-Baptiste Lagaert, LEGI
    use precision_tools
    !USE mpi
    include 'mpif.h'

    !> @param[in]   code            = mpi error code
    !> @param[in]   procedure_name  = name of the procedure from where the error occurs
    integer, intent(in)                             :: code
    character, dimension(:), intent(in), optional   :: procedure_name

    character, dimension(MPI_MAX_ERROR_STRING)  :: message
    integer                                     :: message_lengh, ierr

    call mpi_error_string(code, message, message_lengh,ierr)

    if(present(procedure_name)) then
        print*, '[INFO] Error occurs in procedure ', procedure_name
    end if

    IF(ierr .NE.0) THEN
      WRITE(6,'(a)')'[ERROR] parallel_error: unable to provide any error message'
    ENDIF

    print*, '   -> error message : ', message

END SUBROUTINE parallel_error

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
 FUNCTION biginteger_DoGlobalSum(spec_rank,value) RESULT(resultat)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
  
  INTEGER(DI), INTENT(IN) :: value
  INTEGER, INTENT(IN)     :: spec_rank
  INTEGER(DI)             :: resultat
  INTEGER                 :: ierr
  
  !CALL MPI_ALLREDUCE(value,resultat,1,MPI_REAL_WP,MPI_SUM,spec_communicator,ierr)
  CALL MPI_ALLREDUCE(value,resultat,1,MPI_INTEGER8,MPI_SUM,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalSum on process ",spec_rank,&
                      &": Unable to run global sum."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  
 END FUNCTION biginteger_DoGlobalSum


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_SEND for integer values 
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a Fortran value of datas to send.
!> @param[in] sendType  the type of item to send
!> @param[in] sendTo the rank of the process we send to.
!> @param[in] sendTag the tag of the send communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data exchange is successfull.
!------------------------------------------------------------------------------

FUNCTION integer_DoSendScal(senddata,sendTo,sendTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    !I/O data
    INTEGER, INTENT(IN) :: sendTo,sendTag
    INTEGER, INTENT(IN) :: me
    INTEGER, INTENT(IN)    :: senddata
    LOGICAL :: success 
    !Local data
    INTEGER :: res

    success=.FALSE.

    CALL MPI_SEND(senddata,1,MPI_INTEGER,sendTo,sendTag,spec_communicator,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in integer_DoSendScal."
        RETURN
    ENDIF
    success=.TRUE.
    RETURN

END FUNCTION integer_DoSendScal


!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_SEND for integer values 
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in,out] recvdata a Fortran value to store the received datas.
!> @param[in] recvType the type of item we receive.
!> @param[in] recvFrom the rank of the process we receive from.
!> @param[in] recvTag the tag of the receive communication
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data receive is successfull.
!------------------------------------------------------------------------------

FUNCTION integer_DoRecvScal(recvdata,recvFrom,recvTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    !I/O data
    INTEGER, INTENT(IN) :: recvFrom,recvTag
    INTEGER, INTENT(IN) :: me
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER, INTENT(INOUT) :: recvdata
    LOGICAL :: success 
    !Local data
    INTEGER :: res

    success=.FALSE.

    CALL MPI_RECV(recvdata,1,MPI_INTEGER,recvFrom,recvTag,spec_communicator,status,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": receive failure in real_DoRecvScal."
        RETURN
    ENDIF
    success=.TRUE.

    RETURN

END FUNCTION integer_DoRecvScal


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
FUNCTION integer_DoGlobalVectorSum(spec_rank,vector)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
  
  INTEGER, INTENT(IN) :: vector(:)
  INTEGER :: integer_DoGlobalVectorSum (size(vector))
  INTEGER, INTENT(IN)  :: spec_rank

!  REAL(WP) :: resultat
  INTEGER  ::ierr
  
  CALL MPI_ALLREDUCE(vector,integer_DoGlobalVectorSum,size(vector),MPI_INTEGER,MPI_SUM,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalVectorSum on process ",spec_rank,&
                      &": Unable to run global sum on the array."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
!  real_DoGlobalSum=resultat
  RETURN
  
END FUNCTION integer_DoGlobalVectorSum

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

FUNCTION integer_DoSendRecvVector(senddata,sendCount,sendTo,sendTag,recvdata,recvCount,recvFrom,recvTag,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendCount,sendTo,sendTag
    INTEGER, INTENT(IN) :: recvCount,recvFrom,recvTag
    INTEGER, INTENT(IN) :: me
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER, DIMENSION(:), INTENT(IN)    :: senddata
    INTEGER, DIMENSION(:), INTENT(INOUT) :: recvdata
    LOGICAL :: success 

    INTEGER, DIMENSION(:), ALLOCATABLE:: sendbuf,recvbuf
    INTEGER :: res
    INTEGER :: ires

    success=.FALSE.

    ALLOCATE(sendbuf(SIZE(senddata)), &
           & recvbuf(SIZE(recvdata)), &
           & stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in real_DoSendRecvVector!"
        RETURN
    ENDIF
    sendbuf=senddata
    recvbuf=0.0_WP
    CALL MPI_SENDRECV(sendbuf,sendCount,MPI_INTEGER,sendTo,sendTag,recvBuf,recvCount,MPI_INTEGER,recvFrom,&
                     &recvTag,spec_communicator,status,res)
    IF (res.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,3(i0,a))')"[ERROR] on process ",me,": send/receive failure between",sendTo," and &
                               &",recvFrom," in real_DoSendRecvVector."
        RETURN
    ENDIF
    recvdata=recvbuf
    DEALLOCATE(sendbuf,recvbuf)

    success=.TRUE.
    RETURN

END FUNCTION integer_DoSendRecvVector

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
FUNCTION integer_DoGlobalSum(spec_rank,value)
!------------------------------------------------------------------------------

  USE precision_tools
  USE communicators_tools
  USE mpi
  
  IMPLICIT NONE
  
  INTEGER :: integer_DoGlobalSum 
  INTEGER, INTENT(IN) :: value
  INTEGER, INTENT(IN) :: spec_rank

  INTEGER :: resultat
  INTEGER :: ierr
  
  CALL MPI_ALLREDUCE(value,resultat,1,MPI_INTEGER,MPI_SUM,spec_communicator,ierr)
  IF (ierr .NE. MPI_SUCCESS) THEN
    WRITE(6,'(a,i0,a)')"ERROR real_DoGlobalSum on process ",spec_rank,&
                      &": Unable to run global sum."
    CALL MPI_FINALIZE(ierr)
    STOP 'FAILED'
  ENDIF
  integer_DoGlobalSum=resultat
  RETURN
  
END FUNCTION integer_DoGlobalSum

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_GATHER for REAL values
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

FUNCTION integer_DoBcastVector(sendrecvdata,sendrecvCount,root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  USE mpilayout_tools , only : getnbcpus
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: sendrecvCount
    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(IN) :: me
    INTEGER, DIMENSION(:), INTENT(INOUT) :: sendrecvdata
    LOGICAL :: success 

    INTEGER, DIMENSION(:), ALLOCATABLE:: sendrecvbuf
    INTEGER :: ires,nbproc

    success=.FALSE.
    nbproc = getnbcpus()

    IF ( SIZE(sendrecvdata) .LT. sendrecvCount) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": sendrecvCount disagree with &
                            & SIZE(sendrecvdata) in real_DoScatterVector!"
        RETURN
    ENDIF

    ALLOCATE(sendrecvbuf(sendrecvCount), stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": not enought memory for buffers &
                            &in real_DoBcastVector!"
        RETURN
    ENDIF
    IF ( me .EQ. root) THEN
        sendrecvbuf=sendrecvdata
    ELSE
        sendrecvbuf=0.0_WP
    ENDIF

    CALL MPI_BCAST(sendrecvbuf,sendrecvCount,MPI_INTEGER,root,spec_communicator,ires)

    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in real_DoScatterVector."
        RETURN
    ENDIF
    sendrecvdata=sendrecvbuf
    DEALLOCATE(sendrecvbuf)

    success=.TRUE.
    RETURN

END FUNCTION integer_DoBcastVector

FUNCTION integer_DoBcastScal(sendrecvdata,root,me) result(success)

  USE precision_tools
  USE mpi
  USE communicators_tools
  IMPLICIT NONE

    INTEGER, INTENT(IN) :: root
    INTEGER, INTENT(IN) :: me
    INTEGER, INTENT(INOUT) :: sendrecvdata
    LOGICAL :: success 

    INTEGER :: ires

    success=.FALSE.

    CALL MPI_BCAST(sendrecvdata,1,MPI_INTEGER,root,spec_communicator,ires)

    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in real_DoScatterVector."
        RETURN
    ENDIF

    success=.TRUE.
    RETURN

END FUNCTION integer_DoBcastScal

!------------------------------------------------------------------------------
!> @author
!> Patrick BEGOU modified by Antoine Vollant
!
!> @details
!> This function is a wrapper of MPI_ALLGATHER for REAL values
!> Be careful, there is an interface to use this function out of this module
!
!> @param[in] senddata a REAL value
!> @param[in,out] recvdata a Fortran 1D array to store the received datas. Depending on compilers this array can have a non contiguous storage.
!> @param[in] me the rank of this process in the current communicator
!> @return  .TRUE. if data reassemble is successfull.
!------------------------------------------------------------------------------
FUNCTION integer_DoAllGatherScal(senddata,recvdata,me) result(success)

    USE precision_tools
    USE mpi
    USE communicators_tools
    USE mpilayout_tools , only : getnbcpus
    IMPLICIT NONE

    INTEGER              , INTENT(IN)       :: senddata
    INTEGER, DIMENSION(:), INTENT(INOUT)    :: recvdata
    INTEGER, INTENT(IN) :: me
    LOGICAL :: success 
    INTEGER :: ires
    success=.FALSE.

    IF ( SIZE(recvdata) .LT. getnbcpus()) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": recvCount disagree with &
                            & the number of processors in integer_DoAllGatherScal!"
        RETURN
    ENDIF
    CALL MPI_ALLGATHER (senddata,1,MPI_INTEGER,recvdata,1,MPI_INTEGER, spec_communicator ,ires)
    IF (ires.NE.MPI_SUCCESS) THEN
        WRITE(6,'(a,i0,a)')"[ERROR] on process ",me,": send faillure in integer_DoAllGatherScal."
        RETURN
    ENDIF
    success=.TRUE.
    RETURN

END FUNCTION integer_DoAllGatherScal

END MODULE parallel_tools
!> @}
