!USEFORTEST toolbox
!USEFORTEST postprocess
!USEFORTEST advec
!USEFORTEST io
!USEFORTEST topo
!USEFORTEST avgcond 
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: mpilayout_tools
!
!> @author
!> Patrick BEGOU, LEGI
!
! DESCRIPTION: 
!> The aim of this module is to compute the local layout of
!> the datas on each processor from the global data size and
!> the number of processors
!------------------------------------------------------------------------------
MODULE mpilayout_tools

   USE options

   IMPLICIT NONE

   !Internals subroutines & functions
   PRIVATE computeProcsDecomp
   PRIVATE computeProcsinXinY
   INTEGER, PRIVATE :: totalCpu=0
   INTEGER, PRIVATE :: world_rank=-1
   INTEGER, PRIVATE :: spec_rank=-1
   !> To save the decomposition
   INTEGER, DIMENSION(2), PRIVATE :: decomp_cpu
   !> To check if the decomposition has already been done
   LOGICAL, PRIVATE :: decomp_init

   
   CONTAINS

!------------------------------------------------------------------------------
!> @author
!> Patrick B.
!
!> \brief
!> @details
!> This returns the total number of cpus available.
!> @return the total cpu number
!------------------------------------------------------------------------------
   FUNCTION getnbcpus()
   IMPLICIT NONE
   INTEGER getnbcpus
   IF (totalCpu .EQ.0) THEN
      WRITE(6,'(a)') '[ERROR] Mpilayout module not initialized'
      STOP 'FAILED'
   ENDIF
   getnbcpus=totalCpu
   RETURN
   END FUNCTION getnbcpus

!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> @details
!> This returns the total number of cpus available.
!> @return the total cpu number
!------------------------------------------------------------------------------
   SUBROUTINE getdecompcpus(cpuY, cpuZ)
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: cpuY, cpuZ
   IF (totalCpu .EQ.0) THEN
      WRITE(6,'(a)') '[ERROR] Mpilayout module not initialized'
      STOP 'FAILED'
   ENDIF
   cpuY = decomp_cpu(1)
   cpuZ = decomp_cpu(2)
   RETURN
   END SUBROUTINE getdecompcpus
!------------------------------------------------------------------------------
!> This subroutine initialise the info about rank in the spectral communicator.
!! @author =  Jean-Baptiste Lagaert
!!    @param[in]    rank_spec   = rank inside the spectral communicator
!------------------------------------------------------------------------------
subroutine init_spec_rank(rank_spec)

    integer, intent(in) :: rank_spec

    IF (spec_rank /= -1) THEN
        WRITE(6,'(a)') '[WARNIG] spec_rank in mpilayout will be reinitialized.'
    ENDIF
    spec_rank = rank_spec

end subroutine init_spec_rank

!------------------------------------------------------------------------------
!> @author
!> Patrick B.
!
!> \brief
!> @details
!> This returns my rank in the spectral communicator.
!> @return the total cpu number
!------------------------------------------------------------------------------
   FUNCTION getnbmycpu()
   IMPLICIT NONE
   INTEGER getnbmycpu
   IF (spec_rank .EQ. -1) THEN
      WRITE(6,'(a)') '[ERROR] Mpilayout module not initialized'
      STOP 'FAILED'
   ENDIF
   getnbmycpu=spec_rank
   RETURN 
   END FUNCTION getnbmycpu
!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> @details
!> This returns my rank in the spectral communicator.
!> @return the total cpu number
!------------------------------------------------------------------------------
   FUNCTION getnbmycpu_world()
   IMPLICIT NONE
   INTEGER getnbmycpu_world
   IF (world_rank .EQ. -1) THEN
      WRITE(6,'(a)') '[ERROR] Mpilayout module not initialized'
      STOP 'FAILED'
   ENDIF
   getnbmycpu_world=world_rank
   RETURN 
   END FUNCTION getnbmycpu_world
!------------------------------------------------------------------------------
!> @author 
!> Patrick B.
!
!> \brief
!> 2D cpu distribution
!
!> @details
!> This subroutine cpus distribution for the two distributed dimensions
!> of the domain (the third direction is not distributed over the cpu but
!> fully located on one cpu.
!
!> @param[in] ncpus the number of available cpus
!> @param[out] ncpus1 the number of cpus in the first distributed dimension
!> @param[out] ncpus2 the number of cpus in the second distributed dimension
!> @param[in] me the rank number of the process in the range [0,ncpus - 1].
!> @param[out] res is set to true if the subroutine is successfull.
!------------------------------------------------------------------------------
SUBROUTINE computeProcsMatrix(ncpus,ncpus1,ncpus2,me,res)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::ncpus
    INTEGER, INTENT(OUT)::ncpus1,ncpus2
    INTEGER, INTENT(IN) ::me
    LOGICAL, INTENT(OUT)::res
    INTEGER::layout(2,4) ! the decomposition
    
    res=.FALSE.
    IF (.not. decomp_init) THEN 
          ! Compute the best decomposition
      
       IF (me.GE.0 .AND. me .LT. ncpus) THEN
          IF (totalCpu .EQ. 0) THEN
              totalCpu = ncpus
          ELSE IF (totalCpu .NE. ncpus) THEN
          WRITE(6,101)'total cpu number',totalCpu,ncpus
              RETURN
          ENDIF
          IF (world_rank .EQ. -1) THEN
              totalCpu = ncpus
          ELSE IF (world_rank.NE.me) THEN
              WRITE(6,101)'cpu id',world_rank,me
              RETURN
          ENDIF
          world_rank = me
          ELSE
              WRITE(6,100) ncpus,me
              RETURN
       ENDIF 
  
       ncpus1=GetNbCpus1()
       IF (ncpus1.GT.0) THEN
         ! distribution guided by user choice
         IF (mod (ncpus,ncpus1) .NE. 0) THEN
            IF (me.EQ.0) WRITE (6,'(3(a,i4),a)') &
            & "[ERROR ]The total number of cpus (",ncpus, &
            & ") cannot divides by your imposed distribution along X (",ncpus1,&
            & ") and ",mod (ncpus,ncpus1)," remains unused. It is not allowed." 
            ncpus1=0
            ncpus2=0
            RETURN
         ELSE
            ncpus2=ncpus/ncpus1
            res=.TRUE.
         ENDIF                 
       ELSE
          layout=RESHAPE((/2, 0, 3, 0, 5, 0, 7, 0/),(/2,4/))
          CALL computeProcsDecomp(ncpus, layout,res)
          IF (.NOT. res) return
          CALL computeProcsinXinY(layout,ncpus1,ncpus2)
        ENDIF
        decomp_init = .true.
        decomp_cpu = (/ncpus1, ncpus2 /)
        !WRITE(6,'("Matrix for ",i4," is [",i0," x ",i0,"]")') ncpus, ncpus1,ncpus2
    ELSE
        ! Check if nothing has changed
        IF (totalCpu .NE. ncpus) THEN
        WRITE(6,101)'total cpu number',totalCpu,ncpus
            RETURN
        ENDIF
        IF (spec_rank /= -1) THEN
            IF (spec_rank.NE.me) THEN
                WRITE(6,101)'cpu id',world_rank,me
                RETURN
            END IF
        END IF
        ncpus1 = decomp_cpu(1)
        ncpus2 = decomp_cpu(2)
        res = .true.
    END IF
100   FORMAT ('[ERROR] computeProcsMatrix: strange number of cpu [',i0, &
            &'] or cpu id [',i0,'] (0 < cpu id < maxcpu)')
101   FORMAT ('[ERROR] computeProcsMatrix: value of ',a, &
            & ' previously set to [',i0,&
            & '] and now set to [',i0,'] ?')
   END SUBROUTINE computeProcsMatrix

   
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!>
!>
!> This subroutine divide ncpus in a product of prime numbers
!> including [2,3,5,7].
!>
!> @param [in] ncpus the number of available cpus.
!> @param [in,out] layout a two columns array. The first column provides the prime numbers for the decomposition. 
!!>  The second columns contains the factors calculated by this subroutine. 
!> @param [out] res is set to true if the subroutine is successfull.
!------------------------------------------------------------------------------
SUBROUTINE computeProcsDecomp(ncpus, layout,res)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::ncpus
    INTEGER, INTENT(INOUT)::layout(:,:) ! the decomposition
    LOGICAL, INTENT(OUT)::res
    INTEGER::cpt, cptold, i

    res=.TRUE.
    cpt=ncpus
    DO WHILE (cpt.GT.1)
        cptold=cpt
        DO i=1,size(layout,2)
            IF (mod(cpt,layout(1,i)).EQ.0) THEN
                layout(2,i)=layout(2,i)+1
                cpt=cpt/layout(1,i)
                EXIT
            ENDIF
        ENDDO
        IF(cpt .EQ. cptold) THEN
            WRITE(6,'(a)') "The number of processor should divide by 2,3,5 or 7"
            WRITE(6,'(a,i0,a)') "The value of ",ncpus," does not."
            cpt=0
            res=.FALSE.
        ENDIF    
    ENDDO
END SUBROUTINE computeProcsDecomp

   
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!>
!>
!> This subroutine calculates the number of cpus in each direction
!>
!> @param [in] layout a two columns array. The first column provides the 
!> prime numbers for the decomposition and the second the occurences. 
!> @param [out] ncpus1 the number of available cpus in the X direction (layout
!> along Y or Z) or the Y direction when layout is along X.
!> @param [out] ncpus2 the number of available cpus in the Z direction (layout
!> along X or Y) or the Y direction when layout is along Z.
!>
!> \image html xycommunicator.png "Data layout change with XY_communicator for 6 cpus"
!> \image html yzcommunicator.png "Data layout change with YZ_communicator for 6 cpus"
!> 
!> For spectral domain organisation, ncpus1 is set lower or equal to ncpus2
!------------------------------------------------------------------------------
SUBROUTINE computeProcsinXinY(layout,ncpus1,ncpus2)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) ::layout(:,:) ! the decomposition
    INTEGER, INTENT(OUT):: ncpus1,ncpus2  

    INTEGER:: i
    LOGICAL:: ok
       
    ok=.TRUE.
    ncpus1=1
    ncpus2=1

    DO i=size(layout,2),1,-1
        DO WHILE (layout(2,i) .GT.0)
            IF(ok) THEN
                ncpus1=ncpus1*layout(1,i)
                IF (ncpus1.GT.ncpus2) ok=.NOT.ok
            ELSE
                ncpus2=ncpus2*layout(1,i)
                IF (ncpus2.GT.ncpus1) ok=.NOT.ok
            ENDIF
            layout(2,i)=layout(2,i)-1
        ENDDO
    ENDDO
    IF (ncpus1 .GT. ncpus2) THEN
       i=ncpus1
       ncpus1=ncpus2
       ncpus2=i
    ENDIF
   
END SUBROUTINE computeProcsinXinY
   
!------------------------------------------------------------------------------
!> @author 
!> Patrick BEGOU, LEGI
!
!>
!> @details
!> This subroutine give the number of points on each cpu 
!> for the current dimension.
!
!> @param [in] ncpus the number of available cpus.
!> @param [in] nbpt the number of points along the current dimension
!> @return a 2D array of start/end position for each cpu
!------------------------------------------------------------------------------
FUNCTION spreadOncpus(ncpus,nbpt)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  ::ncpus,nbpt

    INTEGER,DIMENSION(ncpus,2)::spreadOncpus

    INTEGER :: rab, supp, slice, i
    INTEGER :: debut, fin


    IF (nbpt .LT. ncpus) THEN
        WRITE(6,'(a,i0,a,i0,a)')"Less grid points (",nbpt,") than cpus (",ncpus,")"
        WRITE(6,'(a)') "This config is not allowed."
        STOP 'FAILED'
    ELSE
        rab=mod(nbpt,ncpus) ! si nbpts ne se divise pas par ncpus
        slice=nbpt/ncpus    ! taille moyenne par cpu sinon
        supp = 0            ! pour répartir le rab sur les cpus
      
        DO i=1, ncpus
            debut=(i-1)*slice+supp+1
            IF (rab .GT.0) THEN
                supp=supp+1
                rab=rab-1   
            ENDIF
            fin= i * slice + supp 
        
!        WRITE(6,'("cpu: ",i2," [",i0," , ",i0,"]")'),i,debut,fin 
            spreadOncpus(i,1)=debut
            spreadOncpus(i,2)=fin
        ENDDO
    ENDIF
END FUNCTION spreadOncpus
   
   
END MODULE mpilayout_tools
!> @}
