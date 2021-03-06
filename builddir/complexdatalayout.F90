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
! MODULE: cmplxdatalayout
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

MODULE cmplxdatalayout

USE maindatalayout
USE precision_tools
USE mpilayout_tools
USE parser_tools

IMPLICIT NONE


    INTEGER, PRIVATE, PARAMETER ::lenName=str_medium

    TYPE COMPLEX_DATA_LAYOUT
        CHARACTER(LEN=lenName) :: name !< The name of the datas U, V, scalar 1....etc
        INTEGER :: nx         !< Global grid size for the datas
        INTEGER :: ny         !< Global grid size for the datas
        INTEGER :: nz         !< Global grid size for the datas
        ! curlayout should be alongX or alongY or alongZ
        INTEGER :: curlayout  !< should be alongX or alongY or alongZ
        INTEGER :: xmin  !< start X position in the global array for the current layout
        INTEGER :: xmax  !< final X position in the global array for the current layout
        INTEGER :: ymin  !< start Y position in the global array for the current layout
        INTEGER :: ymax  !< final Y position in the global array for the current layout
        INTEGER :: zmin  !< start Z position in the global array for the current layout
        INTEGER :: zmax  !< final Z position in the global array for the current layout

        !> The datas (U, V....). The array size is (xmin:xmax,ymin:ymax,zmin:zmax)
        COMPLEX(WP), DIMENSION(:,:,:), POINTER :: values => NULL()
        !> A ncpus x 2 array: ymin and ymax values for each cpu when the layout is along X axis
        INTEGER, DIMENSION(:,:),POINTER :: Y4layoutOnX => NULL()
        !> A ncpus x 2 array: zmin and zmax values for each cpu when the layout is along X axis
        INTEGER, DIMENSION(:,:),POINTER :: Z4layoutOnX => NULL()
        !> A ncpus x 2 array: xmin and xmax values for each cpu when the layout is along Y axis
        INTEGER, DIMENSION(:,:),POINTER :: X4layoutOnY => NULL()
        !> A ncpus x 2 array: zmin and zmax values for each cpu when the layout is along Y axis
        INTEGER, DIMENSION(:,:),POINTER :: Z4layoutOnY => NULL()
        !> A ncpus x 2 array: xmin and xmax values for each cpu when the layout is along Z axis
        INTEGER, DIMENSION(:,:),POINTER :: X4layoutOnZ => NULL()
        !> A ncpus x 2 array: ymin and ymax values for each cpu when the layout is along Z axis
        INTEGER, DIMENSION(:,:),POINTER :: Y4layoutOnZ => NULL()
        !> Physical dimensions
        REAL(WP) ::Lx,Ly,Lz

    END TYPE COMPLEX_DATA_LAYOUT


    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE cmplx_affectation
    END INTERFACE

    ! this internal subroutine inherited from maindatalayout is set to PRIVATE
    ! it should not be seen out of this module.
    PRIVATE initLayoutArray
    PRIVATE cmplx_affectation

CONTAINS

!------------------------------------------------------------------------------
!> @author
!> Patrick Begou
!
!> \brief
!> Data layout initialisation on each process
!
!> @details
!> This subroutine allocate and organise the data for the local process, taking account
!> for the local distribution along on the X, Y or Z axis. Default is along X axis.
!
!> @param[in] name the name of this data (U, pressure, scalar_1...). It should not contains spaces, tabs nor special characters as it can be used to create file names. 
!> @param[in,out] val is the variable of type DATA_LAYOUT to be instanciated
!> @param[in] nx global domain discretization size along X axis
!> @param[in] ny global domain discretization size along Y axis
!> @param[in] nz global domain discretization size along Z axis
!> @param[in] Lx global domain lenght along X axis
!> @param[in] Ly global domain lenght along Y axis
!> @param[in] Lz global domain lenght along Z axis
!> @param[in] nbcpus the total number of available cpus 
!> @param[in] me my rank in the cpus pool. Range in  [0,nbcpus-1]
!> @param[in] how the choosen distribution. This parameter is optional and defaults to alongX. Acceptable values are alongX, alongY or alongZ.
!> @return true if memory is successfuly allocated
!------------------------------------------------------------------------------
LOGICAL FUNCTION cmplx_initDataLayout(name,val,nx,ny,nz,Lx,Ly,Lz,nbcpus,me,how)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER, INTENT(IN)           :: nx,ny,nz
    REAL(WP), INTENT(IN)          :: Lx,Ly,Lz
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
    INTEGER, INTENT(IN)           :: nbcpus,me
    INTEGER, INTENT(IN),OPTIONAL  :: how   !default to alongX

    LOGICAL ::res
    INTEGER :: ncpus1,ncpus2,ires

    cmplx_initDataLayout=.FALSE.
    
    !if something is allocated, deallocate it
    IF(ASSOCIATED(val%values)) DEALLOCATE(val%values)
    IF(ASSOCIATED(val%Y4layoutOnX)) DEALLOCATE(val%Y4layoutOnX)
    IF(ASSOCIATED(val%Z4layoutOnX)) DEALLOCATE(val%Z4layoutOnX)
    IF(ASSOCIATED(val%X4layoutOnY)) DEALLOCATE(val%X4layoutOnY)
    IF(ASSOCIATED(val%Z4layoutOnY)) DEALLOCATE(val%Z4layoutOnY)
    IF(ASSOCIATED(val%X4layoutOnZ)) DEALLOCATE(val%X4layoutOnZ)
    IF(ASSOCIATED(val%Y4layoutOnZ)) DEALLOCATE(val%Y4layoutOnZ)

    IF (LEN_TRIM(name).LT.lenName) THEN
        val%name=name
    ELSE
        val%name=name(1:lenName)
    ENDIF

    val%nx=nx
    val%ny=ny
    val%nz=nz
        
    ! how many cpus in the twho decomposed dimensions ?
    CALL computeProcsMatrix(nbcpus,ncpus1,ncpus2,me,res)
    IF (.NOT. res) return

    !get the data organisation on the cpus
    IF (.NOT. ( &
        & initLayoutArray(val%Y4layoutOnX,val%Z4layoutOnX,ncpus1,ncpus2,ny,nz)&
        & .AND. &
        & initLayoutArray(val%X4layoutOnY,val%Z4layoutOnY,ncpus1,ncpus2,nx,nz)&
        & .AND. &
        & initLayoutArray(val%X4layoutOnZ,val%Y4layoutOnZ,ncpus1,ncpus2,nx,ny)&
        )) THEN
        WRITE(6,'(a,i0,a)')"ERROR on process ",me,": not enought memory!"
        RETURN
    ENDIF

        ! get the initial datalayout
    IF (present(how)) THEN
        IF (how.EQ.alongX .OR. how .EQ. alongY .OR. how .EQ. alongZ) THEN
            val%curlayout=how
        ELSE
            WRITE(6,'(a,i0,a,i0)')"ERROR on process ",me,": unknown layout ",how
            WRITE(6,'(a)')        "      Should be in [1,3]"
        ENDIF
    ELSE
        !LUCA: was alongX 
        val%curlayout=alongX  
    ENDIF

    !Allocate data array, we know curlayout is right so no default case.
    !we decide do have realistic indexes, not sure it is a good idea!
    SELECT CASE (val%curlayout)
    CASE (alongX)
        val%xmin=1
        val%xmax=nx
        val%ymin=val%Y4layoutOnX(me+1,1)
        val%ymax=val%Y4layoutOnX(me+1,2)
        val%zmin=val%Z4layoutOnX(me+1,1)
        val%zmax=val%Z4layoutOnX(me+1,2)
    CASE (alongY)
        val%xmin=val%X4layoutOnY(me+1,1)
        val%xmax=val%X4layoutOnY(me+1,2)
        val%ymin=1
        val%ymax=ny
        val%zmin=val%Z4layoutOnY(me+1,1)
        val%zmax=val%Z4layoutOnY(me+1,2)
    CASE (alongZ)
        val%xmin=val%X4layoutOnZ(me+1,1)
        val%xmax=val%X4layoutOnZ(me+1,2)
        val%ymin=val%Y4layoutOnZ(me+1,1)
        val%ymax=val%Y4layoutOnZ(me+1,2)
        val%zmin=1
        val%zmax=nz
    END SELECT

    Val%Lx=Lx
    Val%Ly=Ly
    Val%Lz=Lz

    ALLOCATE(val%values(val%xmin:val%xmax ,val%ymin:val%ymax , val%zmin:val%zmax),stat=ires)
    IF (ires .NE.0) THEN
        WRITE(6,'(a,i0,a)')"ERROR on process ",me,": not enought memory for datas!"
        RETURN
    ELSE
        val%values=0.0_WP
        cmplx_initDataLayout=.TRUE.
    ENDIF
        RETURN

END FUNCTION cmplx_initDataLayout
    
    
!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> \brief
!> This routine should be called to free memory allocated in a DATA_LAYOUT variable.
!
!
!> @param[in,out] val is the variable of type DATA_LAYOUT to be freed
!------------------------------------------------------------------------------
    SUBROUTINE cmplx_deleteDataLayout(val)
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(INOUT) :: val
    IF(ASSOCIATED(val%values))  DEALLOCATE(val%values)
    IF(ASSOCIATED(val%Y4layoutOnX))  DEALLOCATE(val%Y4layoutOnX)
    IF(ASSOCIATED(val%Z4layoutOnX))  DEALLOCATE(val%Z4layoutOnX)
    IF(ASSOCIATED(val%X4layoutOnY))  DEALLOCATE(val%X4layoutOnY)  
    IF(ASSOCIATED(val%Z4layoutOnY))  DEALLOCATE(val%Z4layoutOnY) 
    IF(ASSOCIATED(val%X4layoutOnZ))  DEALLOCATE(val%X4layoutOnZ)
    IF(ASSOCIATED(val%Y4layoutOnZ))  DEALLOCATE(val%Y4layoutOnZ)
    val%values => NULL()
    val%Y4layoutOnX => NULL()
    val%Z4layoutOnX => NULL()
    val%X4layoutOnY => NULL()
    val%Z4layoutOnY => NULL()
    val%X4layoutOnZ => NULL()
    val%Y4layoutOnZ => NULL()
    val%xmin=0
    val%xmax=-1
    val%ymin=0
    val%ymax=-1
    val%zmin=0
    val%zmax=-1
    END SUBROUTINE cmplx_deleteDataLayout
    
!------------------------------------------------------------------------------
!> @author 
!> Patrick Begou
!
!> \brief
!> Print the current data layout on this node. For debug purpose only.
!> @param[in] val is the variable of type DATA_LAYOUT for wich the layout is printed
!> @param[in] me my rank in the cpus pool. Range in  [0,nbcpus-1]
!------------------------------------------------------------------------------
    SUBROUTINE cmplx_showlayout(val,me)

    IMPLICIT NONE
    INTEGER, INTENT(IN)::me
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN) :: val
    INTEGER ::i
    WRITE(me+10,'(a)')"===================================================================="
    WRITE(me+10,'(a)')'I am a COMPLEX data layout object called "'//TRIM(val%name)//'"'
    WRITE(me+10,'(a)')"===================================================================="
    WRITE(me+10,100) "global layout",1,val%nx,1,val%ny,1,val%nz
    IF (val%curlayout.EQ.alongX) WRITE(me+10,'(a)') "Organised along X"
    IF (val%curlayout.EQ.alongY) WRITE(me+10,'(a)') "Organised along Y"
    IF (val%curlayout.EQ.alongZ) WRITE(me+10,'(a)') "Organised along Z"
    WRITE(me+10,100) "local layout",val%xmin,val%xmax,val%ymin,val%ymax,val%zmin,val%zmax
    
    WRITE(me+10,'(a)') "=========== Y4layoutOnX ==========="
    DO i=1,SIZE(val%Y4layoutOnX,1)
    WRITE(me+10,101) "Y, when alongX",i-1,val%Y4layoutOnX(i,1),val%Y4layoutOnX(i,2)
    ENDDO
    WRITE(me+10,'(a)') "=========== Z4layoutOnX ==========="
    DO i=1,SIZE(val%Z4layoutOnX,1)
    WRITE(me+10,101) "Z, when alongX",i-1,val%Z4layoutOnX(i,1),val%Z4layoutOnX(i,2)
    ENDDO
    
    WRITE(me+10,'(a)') "=========== X4layoutonY ==========="
    DO i=1,SIZE(val%X4layoutonY,1)
    WRITE(me+10,101) "X, when alongY",i-1,val%X4layoutonY(i,1),val%X4layoutonY(i,2)
    ENDDO
    WRITE(me+10,'(a)') "=========== Z4layoutonY ==========="
    DO i=1,SIZE(val%Z4layoutonY,1)
    WRITE(me+10,101) "Z, when alongY",i-1,val%Z4layoutonY(i,1),val%Z4layoutonY(i,2)
    ENDDO
    
    WRITE(me+10,'(a)') "=========== X4layoutonZ ==========="
    DO i=1,SIZE(val%X4layoutonZ,1)
    WRITE(me+10,101) "X, when alongZ",i-1,val%X4layoutonZ(i,1),val%X4layoutonZ(i,2)
    ENDDO
    WRITE(me+10,'(a)') "=========== Y4layoutonZ ==========="
    DO i=1,SIZE(val%Y4layoutonZ,1)
    WRITE(me+10,101) "Y, when alongZ",i-1,val%Y4layoutonZ(i,1),val%Y4layoutonZ(i,2)
    ENDDO
    WRITE(me+10,'(a)')"============================ END ==================================="


100 FORMAT(a," is [",i3,":",i3,",",i3,":",i3,",",i3,":",i3,"]")
101 FORMAT(a," on cpu ",i0," is [",i3,":",i3,"]")


    END SUBROUTINE cmplx_showlayout

!------------------------------------------------------------------------------
!> @author
!> Patrick Begou
!
!> \brief
!> Check that the different COMPLEX_DATA_LAYOUT argument have the same layout :
!> same global dimensions and same local dimensions.
!> @param[in] val is the variable of type COMPLEX_DATA_LAYOUT for wich the layout is checked
!> @param[in] val1 is a variable of type COMPLEX_DATA_LAYOUT to compare with
!> @param[in] val2 is an optional variable of type COMPLEX_DATA_LAYOUT to compare with
!> @param[in] val3 is an optional variable of type COMPLEX_DATA_LAYOUT to compare with
!> @param[in] val4 is an optional variable of type COMPLEX_DATA_LAYOUT to compare with
!> @param[in] val5 is an optional variable of type COMPLEX_DATA_LAYOUT to compare with
!> @return .TRUE. if all the provided variable have the same layout
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION cmplx_samelayout(val,val1,val2,val3,val4,val5) result(res)

    IMPLICIT NONE

    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN) :: val, val1
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN), OPTIONAL :: val2
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN), OPTIONAL :: val3
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN), OPTIONAL :: val4
    TYPE(COMPLEX_DATA_LAYOUT),INTENT(IN), OPTIONAL :: val5

    LOGICAL :: res

    IF (PRESENT(val5)) THEN
      res=(cmplx_samelayout(val,val1) &
         &.AND. &
	 & cmplx_samelayout(val1,val2)&
         &.AND. &
	 & cmplx_samelayout(val2,val3)&
         &.AND. &
	 & cmplx_samelayout(val3,val4)&
         &.AND. &
	 & cmplx_samelayout(val4,val5)&
	 &)
    ELSEIF (PRESENT(val4)) THEN
      res=(cmplx_samelayout(val,val1) &
         &.AND. &
	 & cmplx_samelayout(val1,val2)&
         &.AND. &
	 & cmplx_samelayout(val2,val3)&
         &.AND. &
	 & cmplx_samelayout(val3,val4)&
	 &)
    ELSEIF (PRESENT(val3)) THEN
      res=(cmplx_samelayout(val,val1) &
         &.AND. &
	 & cmplx_samelayout(val1,val2)&
         &.AND. &
	 & cmplx_samelayout(val2,val3)&
	 &)
    ELSEIF (PRESENT(val2)) THEN
      res=(cmplx_samelayout(val,val1) &
         &.AND. &
	 & cmplx_samelayout(val1,val2)&
	 &)
    ELSE
      res= (&
         & val%nx .EQ. val1%nx &
	 & .AND. &
         & val%ny .EQ. val1%ny &
	 & .AND. &
         & val%nz .EQ. val1%nz &
	 & .AND. &
         & val%xmin .EQ. val1%xmin &
	 & .AND. &
         & val%xmax .EQ. val1%xmax &
	 & .AND. &
         & val%ymin .EQ. val1%ymin &
	 & .AND. &
         & val%ymax .EQ. val1%ymax &
	 & .AND. &
         & val%zmin .EQ. val1%zmin &
	 & .AND. &
         & val%zmax .EQ. val1%zmax &
	 & .AND. &
         & val%Lx .EQ. val1%Lx &
	 & .AND. &
         & val%Ly .EQ. val1%Ly &
	 & .AND. &
         & val%Lz .EQ. val1%Lz)
    ENDIF

    RETURN

    END FUNCTION cmplx_samelayout

!------------------------------------------------------------------------------
!> @author
!> Patrick Begou
!
!> \brief
!> Support for the = operator between two type COMPLEX_DATA_LAYOUT:
!> same global dimensions and same local dimensions.
!> @param [in,out] left is the variable of type COMPLEX_DATA_LAYOUT to set
!> @param [in] right is the variable of type COMPLEX_DATA_LAYOUT to copy in left
!------------------------------------------------------------------------------
    SUBROUTINE cmplx_affectation(left, right)

       IMPLICIT NONE
       TYPE (COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: left
       TYPE (COMPLEX_DATA_LAYOUT), INTENT(IN) :: right

       INTEGER :: ires, dim1, dim2
       ! to be shure that a syntax A=A will work we must use buffers
       !> The datas (U, V....). The array size is (xmin:xmax,ymin:ymax,zmin:zmax)
       COMPLEX(WP), DIMENSION(:,:,:), POINTER :: values => NULL()
       !> A ncpus x 2 array: ymin and ymax values for each cpu when the layout is along X axis
       INTEGER, DIMENSION(:,:),POINTER :: Y4layoutOnX => NULL()
       !> A ncpus x 2 array: zmin and zmax values for each cpu when the layout is along X axis
       INTEGER, DIMENSION(:,:),POINTER :: Z4layoutOnX => NULL()
       !> A ncpus x 2 array: xmin and xmax values for each cpu when the layout is along Y axis
       INTEGER, DIMENSION(:,:),POINTER :: X4layoutOnY => NULL()
       !> A ncpus x 2 array: zmin and zmax values for each cpu when the layout is along Y axis
       INTEGER, DIMENSION(:,:),POINTER :: Z4layoutOnY => NULL()
       !> A ncpus x 2 array: xmin and xmax values for each cpu when the layout is along Z axis
       INTEGER, DIMENSION(:,:),POINTER :: X4layoutOnZ => NULL()
       !> A ncpus x 2 array: ymin and ymax values for each cpu when the layout is along Z axis
       INTEGER, DIMENSION(:,:),POINTER :: Y4layoutOnZ => NULL()

       IF (.NOT. ASSOCIATED(right%values)) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] Operator = for TYPE (COMPLEX_DATA_LAYOUT) : receive uninitialized data!'
         CALL MPI_FINALIZE(ires)
         STOP 'FAILED'
       ENDIF

       ALLOCATE(values(right%xmin:right%xmax, right%ymin:right%ymax, right%zmin:right%zmax),stat=ires)
       IF (ires .NE.0) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] Operator = for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
         CALL MPI_FINALIZE(ires)
         STOP 'FAILED'
       ENDIF
       values=right%values

       dim1=size(right%Y4layoutOnX,1)
       dim2=size(right%Y4layoutOnX,2)
       ALLOCATE(Y4layoutOnX(dim1,dim2),Z4layoutOnX(dim1,dim2), &
              & X4layoutOnY(dim1,dim2),Z4layoutOnY(dim1,dim2), &
   	   & X4layoutOnZ(dim1,dim2),Y4layoutOnZ(dim1,dim2), &
   	   & stat=ires)
       IF (ires .NE.0) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] Operator = for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
         CALL MPI_FINALIZE(ires)
         STOP 'FAILED'
       ENDIF
       Y4layoutOnX=right%Y4layoutOnX
       Z4layoutOnX=right%Z4layoutOnX
       X4layoutOnY=right%X4layoutOnY
       Z4layoutOnY=right%Z4layoutOnY
       X4layoutOnZ=right%X4layoutOnZ
       Y4layoutOnZ=right%Y4layoutOnZ

       ! Now free memory of the left operand
       IF (ASSOCIATED(left%values))      DEALLOCATE(left%values)
       IF (ASSOCIATED(left%Y4layoutOnX)) DEALLOCATE(left%Y4layoutOnX)
       IF (ASSOCIATED(left%Z4layoutOnX)) DEALLOCATE(left%Z4layoutOnX)
       IF (ASSOCIATED(left%X4layoutOnY)) DEALLOCATE(left%X4layoutOnY)
       IF (ASSOCIATED(left%Z4layoutOnY)) DEALLOCATE(left%Z4layoutOnY)
       IF (ASSOCIATED(left%X4layoutOnZ)) DEALLOCATE(left%X4layoutOnZ)
       IF (ASSOCIATED(left%Y4layoutOnZ)) DEALLOCATE(left%Y4layoutOnZ)

       ! Rebuild left operand
       left%name       =right%name
       left%nx         =right%nx
       left%ny         =right%ny
       left%nz         =right%nz
       left%curlayout  =right%curlayout
       left%xmin       =right%xmin
       left%xmax       =right%xmax
       left%ymin       =right%ymin
       left%ymax       =right%ymax
       left%zmin       =right%zmin
       left%zmax       =right%zmax
       left%values      => values
       left%Y4layoutOnX => Y4layoutOnX
       left%Z4layoutOnX => Z4layoutOnX
       left%X4layoutOnY => X4layoutOnY
       left%Z4layoutOnY => Z4layoutOnY
       left%X4layoutOnZ => X4layoutOnZ
       left%Y4layoutOnZ => Y4layoutOnZ
       left%lx         =right%lx
       left%ly         =right%ly
       left%lz         =right%lz
    END SUBROUTINE cmplx_affectation

!------------------------------------------------------------------------------
!> @author
!> Patrick Begou
!
!> \brief
!> Support for the copy of the structure of type COMPLEX_DATA_LAYOUT without
!> value initialisation (values are set to 0.0_WP). But all the architecture
!> and parallel informaytion is duplicated.
!> @param[in] ori is the variable model of type COMPLEX_DATA_LAYOUT
!> @param[in,out] copy is the variable of type COMPLEX_DATA_LAYOUT builded
!> @param[in] nameOut is the name of the build variable instead of the name of the ori (optional)
!------------------------------------------------------------------------------
    FUNCTION cmplx_copyStructOnly(ori, copy, nameOut) result(res)

       IMPLICIT NONE
       TYPE (COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: copy
       TYPE (COMPLEX_DATA_LAYOUT), INTENT(IN) :: ori
       CHARACTER(len=*),INTENT(IN),OPTIONAL :: nameOut
       LOGICAL :: res

       INTEGER :: ires, dim1, dim2

       res=.FALSE.

       IF (.NOT. ASSOCIATED(ori%values)) THEN
         WRITE(6,'(a,i0,a)')'[ERROR]copyStructOnly for TYPE (COMPLEX_DATA_LAYOUT) : receive uninitialized data!'
         RETURN
       ENDIF
       IF (ASSOCIATED(copy%values))      DEALLOCATE(copy%values)
       IF (ASSOCIATED(copy%Y4layoutOnX)) DEALLOCATE(copy%Y4layoutOnX)
       IF (ASSOCIATED(copy%Z4layoutOnX)) DEALLOCATE(copy%Z4layoutOnX)
       IF (ASSOCIATED(copy%X4layoutOnY)) DEALLOCATE(copy%X4layoutOnY)
       IF (ASSOCIATED(copy%Z4layoutOnY)) DEALLOCATE(copy%Z4layoutOnY)
       IF (ASSOCIATED(copy%X4layoutOnZ)) DEALLOCATE(copy%X4layoutOnZ)
       IF (ASSOCIATED(copy%Y4layoutOnZ)) DEALLOCATE(copy%Y4layoutOnZ)

       ALLOCATE(copy%values(ori%xmin:ori%xmax, ori%ymin:ori%ymax, ori%zmin:ori%zmax),stat=ires)
       IF (ires .NE.0) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
         RETURN
       ENDIF
       copy%values=0.0_WP

       dim1=size(ori%Y4layoutOnX,1)
       dim2=size(ori%Y4layoutOnX,2)
       ALLOCATE(copy%Y4layoutOnX(dim1,dim2),copy%Z4layoutOnX(dim1,dim2), &
              & copy%X4layoutOnY(dim1,dim2),copy%Z4layoutOnY(dim1,dim2), &
   	   & copy%X4layoutOnZ(dim1,dim2),copy%Y4layoutOnZ(dim1,dim2), &
   	   & stat=ires)
       IF (ires .NE.0) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] copyStructOnly for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
         CALL MPI_FINALIZE(ires)
         STOP 'FAILED'
       ENDIF
       copy%Y4layoutOnX=ori%Y4layoutOnX
       copy%Z4layoutOnX=ori%Z4layoutOnX
       copy%X4layoutOnY=ori%X4layoutOnY
       copy%Z4layoutOnY=ori%Z4layoutOnY
       copy%X4layoutOnZ=ori%X4layoutOnZ
       copy%Y4layoutOnZ=ori%Y4layoutOnZ

       ! Now free memory of the left operand

       ! Rebuild copy operand
       if(present(nameOut)) then
         copy%name     =trim(nameOut)
       else
         copy%name     =ori%name
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
       res=.TRUE.
    END FUNCTION cmplx_copyStructOnly

!------------------------------------------------------------------------------
!! Permute two values of two (or three) datalayout which have the same structure
!! @author
!! Jean-Baptiste Lagaert
!> @details
!! for 2 fields: new 1 = old 2, new 2 = old 1
!! for 3 fields: new 1 = old 3, new 2 = old 1, new 3 = old 2
!! @param[in,out] first   is the first variable of type COMPLEX_DATA_LAYOUT
!! @param[in,out] second  is the second variable of type COMPLEX_DATA_LAYOUT
!! @param[in,out] third   is the second variable of type COMPLEX_DATA_LAYOUT
!! @return        res     = true if everything is ok.
!------------------------------------------------------------------------------
FUNCTION cmplx_permute(first,second,third) result(res)

  IMPLICIT NONE
  ! Input/Output
  TYPE (COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: first, second
  TYPE (COMPLEX_DATA_LAYOUT), INTENT(INOUT), optional :: third
  LOGICAL :: res
  ! Local
  COMPLEX(WP), DIMENSION(:,:,:), POINTER :: tmp => NULL()

  res = .false.

  if(present(third)) then
    if(.not. cmplx_samelayout(first,second, third)) then
      write(6,'(a)')'[ERROR] permute for TYPE (COMPLEX_DATA_LAYOUT) : not same layout!'
      return
    end if

    tmp => second%values
    second%values => first%values
    first%values => third%values
    third%values => tmp ! = old second
  else
    if(.not. cmplx_samelayout(first,second)) then
      write(6,'(a)')'[ERROR] permute for TYPE (COMPLEX_DATA_LAYOUT) : not same layout!'
      return
    end if

    tmp => second%values
    second%values => first%values
    first%values => tmp
  end if

  res = .true.

END FUNCTION cmplx_permute


!------------------------------------------------------------------------------
!! Gives the number of processus where global indices are located
!! @author
!! Antoine Vollant 
!> @details
!! @param[in] i   is the first indice in X direction
!! @param[in] j   is the second indice in Y direction
!! @param[in] k   is the third indice in the Z direction
!! @param[in] Layout is the sample of datalayout on which indices are located 
!! @return    proc_rank = number of processus where (i,j,k) are located
!------------------------------------------------------------------------------
FUNCTION cmplx_foundProcess(i,j,k,layout) result(proc_rank)

  IMPLICIT NONE
  ! Input/Output
  INTEGER , INTENT(IN)                :: i,j,k
  TYPE (COMPLEX_DATA_LAYOUT), INTENT(IN) :: layout
  INTEGER :: proc_rank 

  proc_rank = 1 

    SELECT CASE (layout%curlayout)
    CASE (alongX)
        DO proc_rank = 1 , getnbcpus()
          IF ( j .le. layout%Y4layoutOnX(proc_rank,2) .and. &
               k .le. layout%Z4layoutOnX(proc_rank,2)) EXIT
        ENDDO
    CASE (alongY)
        DO proc_rank = 1 , getnbcpus()
          IF ( i .le. layout%X4layoutOnY(proc_rank,2) .and. &
               k .le. layout%Z4layoutOnY(proc_rank,2)) EXIT
        ENDDO
    CASE (alongZ)
        DO proc_rank = 1 , getnbcpus()
          IF ( i .le. layout%X4layoutOnZ(proc_rank,2) .and. &
               j .le. layout%Y4layoutOnZ(proc_rank,2)) EXIT
        ENDDO
    END SELECT
  
  proc_rank = proc_rank - 1 

END FUNCTION cmplx_foundProcess


END MODULE cmplxdatalayout
!> @}
