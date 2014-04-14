!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST avgcond 
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!> @addtogroup layout 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: maindatalayout
!
!> @author
!> Patrick Begou, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to provide global functionalities and parameters
!> to the underliying datalayout architecture
!------------------------------------------------------------------------------

MODULE maindatalayout

  USE precision_tools

  ! How is organized the bloc on the node
  INTEGER, PARAMETER:: alongX=1
  !> Datas are organised along X axis (no distribution between the processes for this dimension)
  INTEGER, PARAMETER:: alongY=2
  !> Datas are organised along Y axis (no distribution between the processes for this dimension)
  INTEGER, PARAMETER:: alongZ=3
  !> Datas are organised along Z axis (no distribution between the processes for this dimension)

CONTAINS
    
!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> \brief
!> Build the layout array for the datas distribution on all the cpu. Internal use only.
!
!> @details
!> This subroutine allocate and calculate the data layout along axis for each distribution scheme
!> on each cpu. It builds 2D integer arrays of size (ncpus,2) with first and final point along the choosen
!> axis for all the cpu. These arrays are included in each DATA_LAYOUT variable to be used for communication.
!> this routine is for internal use only in this module.
!
!> @param[out] layout1 A ncpus x 2 array:  start and final positions for each cpu when the layout is along the first axis. 
!> Datas are not distributed along the third axis.  
!> @param[out] layout2 A ncpus x 2 array:  start and final positions for each cpu when the layout is along the second axis  
!> Datas are not distributed along the third axis.  
!> @param[in] ncpus1 number of cpus along the first axis.
!> @param[in] ncpus2 number of cpus along the second axis.
!> @param[in] nbpt1 number of points along the first axis.
!> @param[in] nbpt2 number of points along the second axis.
!> @param[in] me my rank in the cpus pool. Range in  [0,nbcpus-1]. Unused but should be usefull if we decide to 
!> restrict the distribution information to the local process only (less memory used).
!> @return true if memory is successfuly allocated
!------------------------------------------------------------------------------
    ! IN THIS CASE LAYOUT IS SOMETHING LIKE X4layoutOnY Z4layoutY...
    FUNCTION initLayoutArray(layout1,layout2,ncpus1,ncpus2,nbpt1,nbpt2)

    USE mpilayout_tools
    IMPLICIT NONE
    LOGICAL ::initLayoutArray
    INTEGER, INTENT(IN)::ncpus1,ncpus2,nbpt1,nbpt2
    INTEGER, DIMENSION(:,:),POINTER ::layout1,layout2
    INTEGER, DIMENSION(:,:),ALLOCATABLE ::tab1,tab2
    
    INTEGER:: res,i,j
    
    initLayoutArray=.FALSE.
    !May be we do not need all this information, only the one for the corrent cpu.
    ALLOCATE(layout1(ncpus1*ncpus2,2),layout2(ncpus1*ncpus2,2),tab1(ncpus1,2),tab2(ncpus2,2),stat=res)
    IF (res.NE.0) THEN
       WRITE(6,'(a)')"ERROR in initLayoutArray: not enought memory!"
       RETURN
    ENDIF
    !RETURNS THE NUM OF POINT PER CPU ALONG A DIRECTION
    tab1=spreadOncpus(ncpus1,nbpt1)
    tab2=spreadOncpus(ncpus2,nbpt2)
    
    DO i=1,ncpus2
      DO j=1,ncpus1
         layout1((i-1)*ncpus1+j,:)=tab1(j,:)
         layout2((i-1)*ncpus1+j,:)=tab2(i,:)
      ENDDO
    ENDDO
    DEALLOCATE(tab1,tab2)
    
    initLayoutArray=.TRUE.
    
    RETURN
    END FUNCTION initLayoutArray

END MODULE maindatalayout
!> @}
