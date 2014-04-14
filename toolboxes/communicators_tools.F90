!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST io
!USEFORTEST topo
!USEFORTEST avgcond 
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: parallel
!
!> @author
!> Patrick BEGOU, LEGI
!
! DESCRIPTION:
!> The aim of this module is to organise all parallel communications for
!> changing data layout, for calculations and for parallel I/O
!------------------------------------------------------------------------------
MODULE communicators_tools

  IMPLICIT NONE

  LOGICAL, PRIVATE:: initialized=.FALSE.
  ! communicators for collective communications
  !> the communicator used for alongX to alongY exchanges for this node
  INTEGER ::XY_COMMUNICATOR
  !> global communicator related to the datalayout used for fft and spectral solver
  INTEGER :: spec_communicator
  !> global rank related to spec_communicator
  INTEGER :: spec_rank_save
  !> the list of nodes in spec_communicatorfor XY_COMMUNICATOR communicator
  INTEGER, POINTER, DIMENSION(:)::XY_Nodes=>null()
  !> the communicator used for alongY to alongZ exchanges for this node
  INTEGER::YZ_COMMUNICATOR
  !> the list of nodes in spec_communicator for YZ_COMMUNICATOR communicator
  INTEGER, POINTER, DIMENSION(:)::YZ_nodes=>null()

  PUBLIC  initCommunicators
  PUBLIC  CommInitialized
  PRIVATE createCommXtoY
  PRIVATE createCommYtoZ

  CONTAINS

!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> Check that module is initialized.
!
!> @details
!> Check that communicators have been built and are rady for use.
!------------------------------------------------------------------------------
!===================================================
  LOGICAL FUNCTION CommInitialized()
!===================================================
  IMPLICIT NONE
  CommInitialized=initialized
  RETURN
  END FUNCTION CommInitialized
  
  
!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> 2D cpu distribution
!
!> @details
!> This subroutine initialized the (optional) mpi topology and
!> creates all communicator used in the code (including the current process 
!> communicators for all collective communications and those used by the
!> particles method).
!> The current process communicators for all collective communications 
!> are used to change the data layout from alongX to alongY or alongZ.
!> @param[in]   ncpus1      = the number of available cpus in direction 1 (Y or X)
!> @param[in]   ncpus2      = the number of available cpus in direction 2 (Z)
!> @param[in]   topo_dim    = dimension of the mpi topology
!> @param[out]  spec_rank   = the current process id inside the communicator used for FFT
!------------------------------------------------------------------------------
!===================================================
  SUBROUTINE initCommunicators(ncpus1, ncpus2, topo_dim, spec_rank)
!===================================================
  USE precision_tools
  use cart_topology

  IMPLICIT NONE
  INTEGER, INTENT(IN)   ::ncpus1, ncpus2, topo_dim
  INTEGER, INTENT(OUT)  :: spec_rank

  INTEGER:: ierr

  ! Initialize mpi topology
  call cart_create((/ncpus1, ncpus2/), ierr, spec_communicator, topo_dim)

  ! Compute rank inside the new communicator
  call mpi_comm_rank(spec_communicator, spec_rank, ierr)
  spec_rank_save = spec_rank


  !> Initialize communicators
  CALL createCommXtoY(ncpus1, ncpus2, spec_rank)
  CALL createCommYtoZ(ncpus1, ncpus2, spec_rank)

  !> Set the elemental data size for communicators
  call init_mpitype()

  initialized=.TRUE.
  RETURN
  END SUBROUTINE initCommunicators

!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> Set the elemental MPI data size
!
!> @details
!> This subroutine adjust the mpi size for real and complex in order to match
!! to the type defined in precision_tools.F90
!------------------------------------------------------------------------------
!===================================================
  subroutine init_mpitype()
!===================================================
  USE precision_tools
  use mpi
  implicit none

  INTEGER:: ierr
  INTEGER:: size_dp

  ! For real
  CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_dp,ierr)
  IF (WP.EQ.size_dp) THEN
      MPI_REAL_WP = MPI_DOUBLE_PRECISION
  ELSE
      CALL MPI_TYPE_SIZE(MPI_REAL,size_dp,ierr)
      IF (WP.EQ.size_dp) THEN
         MPI_REAL_WP = MPI_REAL
      ELSE
         STOP 'Unknown elemental FLOAT size for MPI'
      ENDIF
  ENDIF

  ! For complex
  CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,size_dp,ierr)
  IF ((2*WP).EQ.size_dp) THEN
      MPI_COMPLEX_WP = MPI_DOUBLE_COMPLEX
  ELSE
      CALL MPI_TYPE_SIZE(MPI_COMPLEX,size_dp,ierr)
      IF ((2*WP).EQ.size_dp) THEN
         MPI_COMPLEX_WP = MPI_COMPLEX
      ELSE
         STOP 'Unknown elemental COMPLEX size for MPI'
      ENDIF
  ENDIF
  END  SUBROUTINE init_mpitype

!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> 2D cpu distribution
!
!> @details
!> This subroutine creates the current process communicator for collective communications.
!> The communicator XY_COMMUNICATOR is used for changing the data layout from "along X"
!> to "along Y".
!> \image html xycommunicator.png "Data layout change with XY_communicator for 6 cpus"
!> @param[in] ncpus1 the number of available cpus in direction 1 (Y or X)
!> @param[in] ncpus2 the number of available cpus in direction 2 (Z)
!> @param[in] spec_rank the current process id
!------------------------------------------------------------------------------
!===================================================
SUBROUTINE createCommXtoY(ncpus1, ncpus2, spec_rank)
!===================================================
    USE mpi
    IMPLICIT NONE
    INTEGER, INTENT(IN)::ncpus1, ncpus2, spec_rank

#ifdef DEBUGMPI
    CHARACTER(LEN=64):: form1, form2
#endif

    INTEGER, DIMENSION(ncpus1) :: procArray
    INTEGER:: mpi_group_spec,ierr,newgroup,newcomm
    INTEGER::i,j,res

    DO i=0, ncpus2*ncpus1-1,ncpus1
        procArray=(/(j,j=i,i+ncpus1-1,1)/)
        ! Create communicators
        CALL MPI_COMM_GROUP(spec_communicator, mpi_group_spec,ierr)
        CALL MPI_GROUP_INCL(mpi_group_spec,ncpus1,procArray,newgroup,ierr)
        CALL MPI_COMM_CREATE(spec_communicator,newgroup,newcomm,ierr)

        IF (spec_rank .GE. i .AND.  spec_rank .LT. i+ncpus1) THEN
             ! recupération de la valeur du communicateurt ici
            XY_COMMUNICATOR=newcomm

            !keep trace of the list of nodes for my communicator
            IF(ASSOCIATED(XY_Nodes)) DEALLOCATE(XY_Nodes)
            ALLOCATE(XY_Nodes(ncpus1),stat=res)
            IF (res .NE.0) THEN
                WRITE(6,'(a,i0,a)')"ERROR on process ",spec_rank,": not enought memory!"
                STOP 'FAILED in createCommXtoY'
            ENDIF
            XY_Nodes=procArray
#ifdef DEBUGMPI
            WRITE(form2,'(i0)') ncpus1
            form1='("Process n. ",i2.2," [",'//trim(form2)//'(1x,i2),"] my XY communicator!")'
            WRITE(6,form1) spec_rank,procArray
        ELSE
            WRITE(form2,'(i0)') ncpus1
            form1='("Process n. ",i2.2," [",'//trim(form2)//'(1x,i2),"]")'
            WRITE(6,form1) spec_rank,procArray
#endif
        ENDIF
    ENDDO
END SUBROUTINE createCommXtoY

!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU
!
!> @details
!> This subroutine creates the current process communicator for collective communications.
!> The communicator YZ_COMMUNICATOR is used for changing the data layout from "along Y"
!> to "along Z".
!> \image html yzcommunicator.png "Data layout change with YZ_communicator for 6 cpus"
!> @param[in] ncpus1 the number of available cpus in direction 1 (X)
!> @param[in] ncpus2 the number of available cpus in direction 2 (Y or Z)
!> @param[in] spec_rank the current process id
!------------------------------------------------------------------------------
!===================================================
SUBROUTINE createCommYtoZ(ncpus1, ncpus2, spec_rank)
!===================================================
    USE mpi
    IMPLICIT NONE
    INTEGER, INTENT(IN)::ncpus1, ncpus2, spec_rank
     
#ifdef DEBUGMPI
    CHARACTER(LEN=64):: form1, form2  
#endif
     
    INTEGER, DIMENSION(ncpus2) :: procArray
    INTEGER:: mpi_group_spec,ierr,newgroup,newcomm
    INTEGER::i,j,res
    DO i=0, ncpus1-1
        procArray=(/(j,j=mod(i,ncpus1),ncpus1*ncpus2-1,ncpus1)/)
        !creation du communicateur ici
        CALL MPI_COMM_GROUP(spec_communicator, mpi_group_spec,ierr)
        CALL MPI_GROUP_INCL(mpi_group_spec,ncpus2,procArray,newgroup,ierr)
        CALL MPI_COMM_CREATE(spec_communicator,newgroup,newcomm,ierr)
        IF (any(procArray.EQ.spec_rank)) THEN
            ! Get the communicator value
            YZ_COMMUNICATOR=newcomm
            !keep trace of the list of nodes for my communicator
            IF(ASSOCIATED(YZ_Nodes)) DEALLOCATE(YZ_Nodes)
            ALLOCATE(YZ_Nodes(ncpus2),stat=res)
            IF (res .NE.0) THEN
                WRITE(6,'(a,i0,a)')"ERROR on process ",spec_rank,": not enought memory!"
                STOP 'FAILED in createCommYtoZ'
            ENDIF
            YZ_Nodes=procArray
#ifdef DEBUGMPI
            WRITE(form2,'(i0)') ncpus2
            form1='("Process n. ",i2.2," [",'//trim(form2)//'(1x,i2),"] my YZ communicator!")'
            WRITE(6,form1) spec_rank,procArray
        ELSE   
            WRITE(form2,'(i0)') ncpus2
            form1='("Process n. ",i2.2," [",'//trim(form2)//'(1x,i2),"]")'
            WRITE(6,form1) spec_rank,procArray
#endif
        ENDIF
    ENDDO
END SUBROUTINE createCommYtoZ

END MODULE communicators_tools
!> @}
