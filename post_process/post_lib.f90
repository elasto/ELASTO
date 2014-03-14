!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_lib
!
!
! DESCRIPTION:
!>  This module provide a library of post-process routines. This library
!! provide all that you need to build your own post-process.
!
!> @details
!!  This module provide all post-process implement in scales. It has to
!! be used to build your own post-process traitement. Note that for post-
!! performed during computation, you have to used the module
!! post_in_simu routine (another exe will be add later to perform post-
!! process separately from simulation).
!
!> @author
!! Jean-Baptiste Lagaert, LEGI
!
!------------------------------------------------------------------------------

module post_lib

    use fileio
    use parallel_tools
    use datalayout
    use transforms_tools
    use wavenumber_tools
    use avs
    use subgrid_tools
    use vtkxml_plus
    use data
    use rediscretization_tools
    use toolbox
    use stat_tools
    use parser_tools
    use precision_tools
    use physical_values_tools

    implicit none
    private


    ! ==== Public procedures ====
    public          :: showScalarAvg
    public          :: showVeloMinMaxAvg
    public          :: showPostScalar
    public          :: paravieWrite,paravieWrite_luca
    public          :: avsWrite,avsWrite_luca
    public          :: dataWrite,dataWrite_luca
    public          :: postprocess_double
    public          :: differential_diffusion
    public          :: transDiscAvsWrite
    public          :: compute_all_spectrum
    public          :: compute_all_scalars_spectrum
    public          :: compute_velo_spectrum
    public          :: compute_mhd_spectrum
    public          :: spectrum_time_avg
    public          :: WriteQcOmega
    public          :: compute_Helicity_spectrum
    public          :: compute_Integral_scale
    public          :: compareVeloFields 
    type(REAL_DATA_LAYOUT), PROTECTED, SAVE :: InitScaField
    type(REAL_DATA_LAYOUT), PROTECTED, SAVE :: InitScaDNSField

contains



!------------------------------------------------------------------------------
! ===== Compute min/max/average =====
!> Show min, max and average of a scalar fields.
!! @author = Patrick BEGOU, LEGI
!!    @param[in]    scalarray   = arrays of scalars in physical space
!!    @param[in]    nbscal      = the number of scalars to process
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus      = number of process in the pool
!!    @param[in]    tag         = optional, name of the array of scalar 
!!                                  (spectral or particles scalars ?)
!------------------------------------------------------------------------------
subroutine showScalarAvg(scalarray,sime_time,nbscal,spec_rank,tag)

    !I/O data
    type(REAL_DATA_LAYOUT), intent(in), dimension(:)    :: scalarray
    integer, intent(in)                                 :: nbscal, spec_rank
    real(WP), intent(in)                                :: sime_time
    character (len=*), intent(in), optional             :: tag

    integer :: sca
    real(WP)::mini, maxi, avg, var

    do sca=1,nbscal

        mini = computeFieldMin( scalarray(sca) , spec_rank )
        maxi = computeFieldMax( scalarray(sca) , spec_rank )
        avg  = computeFieldAvg( scalarray(sca) , spec_rank )
        var  = computeFieldVar( scalarray(sca) , spec_rank )

        if (present(tag)) then
            if(spec_rank.EQ.0) write(6,'(a,i2.2,5(1x,a,1x,g15.8))')'    [INFO-STAT] '//trim(tag)//' n', &
                & sca,' time =',sime_time,' min= ',mini,' max= ',maxi,' average= ',avg,' variance= ',var
        else
            if(spec_rank.EQ.0) write(6,'(a,i2.2,4(1x,a,1x,g15.8))')'    [INFO-STAT] Scalar n. ', &
                & sca,' min= ',mini,' max= ',maxi,' average= ',avg,' variance= ',var
        end if
    end do

end subroutine showScalarAvg


!------------------------------------------------------------------------------
!> Show the equality between two fields 
!! @author = Antoine Vollant, LEGI 
!!    @param[in]    U1 =fisrt component of velocity 
!!    @param[in]    V1 =second component of velocity 
!!    @param[in]    W1 =third component of velocity 
!!    @param[in]    U2 =fisrt component of velocity 
!!    @param[in]    V2 =second component of velocity 
!!    @param[in]    W2 =third component of velocity 
!!    @param[in]    spec_rank = mpi rank inside the spectral communicator
!!    @param[inout] success = .true. or .false.
!------------------------------------------------------------------------------
subroutine compareVeloFields(U1,V1,W1,U2,V2,W2,spec_rank,success)


    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)    :: U1,V1,W1,U2,V2,W2 
    integer, intent(in)                   :: spec_rank
    logical                               :: success 

    !Local data
    type(REAL_DATA_LAYOUT)                :: Fieldfil1,Fieldfil2,Fieldfil3
    type(COMPLEX_DATA_LAYOUT)             :: tempK
    integer                               :: nbcpus
    logical                               :: res

    success = .false.
    nbcpus = getnbcpus()
    IF (.NOT. samelayout(U1,V1,W1)) THEN
         write(6,'(a,i0,a)')'[ERROR] in compareVeloFields first vector field has different components'
         return
    ENDIF
    IF (.NOT. samelayout(U2,V2,W2)) THEN
         write(6,'(a,i0,a)')'[ERROR] in compareVeloFields second vector field has different components'
         return
    ENDIF
    IF ((U1%nx .EQ. U2%nx) .AND. (U1%nx .EQ. U1%ny) .AND. (U1%nx .EQ. U1%nz) ) THEN
      IF (equalToVal(U1,U2) .AND. equalToVal(V1,V2) .AND. equalToVal(W1,W2)) THEN
         IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : BOTH VECTOR FIELDS ARE EQUAL'
      ELSE
         IF (.NOT. equalToVal(U1,U2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : U1 != U2' 
         ENDIF
         IF (.NOT. equalToVal(V1,V2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : V1 != V2' 
         ENDIF
         IF (.NOT. equalToVal(W1,W2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : W1 != W2' 
         ENDIF
      ENDIF
    ELSEIF (U1%nx .LT. U2%nx .AND. (U1%nx .EQ. U1%ny) .AND. (U1%nx .EQ. U1%nz) ) THEN 
      IF (.NOT. copyStructonly(U1,Fieldfil1)) RETURN
      IF (.NOT. copyStructonly(V1,Fieldfil2)) RETURN
      IF (.NOT. copyStructonly(W1,Fieldfil3)) RETURN
      IF (.NOT. initDataLayout("tempK", tempK,(U2%nx/2)+1,U2%ny,U2%nz, &
               & U2%Lx,U2%Ly,U2%Lz,nbcpus,spec_rank,alongZ)) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] in compareVeloFields : not enouhgt memory'
         RETURN
      ENDIF
      CALL ftran(U2,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,U1,Fieldfil1) ) RETURN 
      CALL ftran(V2,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,V1,Fieldfil2) ) RETURN 
      CALL ftran(W2,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,W1,Fieldfil3) ) RETURN 
      IF (equalToVal(U1,Fieldfil1) .AND. equalToVal(V1,Fieldfil2) .AND. equalToVal(W1,Fieldfil3)) THEN
         IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : BOTH VECTOR FIELDS ARE EQUAL'
      ELSE
         IF (.NOT. equalToVal(U1,Fieldfil1)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : U1 != U2transdiscretize' 
         ENDIF
         IF (.NOT. equalToVal(V1,Fieldfil2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : V1 != V2transdiscretize' 
         ENDIF
         IF (.NOT. equalToVal(W1,Fieldfil3)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : W1 != W2transdiscretize' 
         ENDIF
      ENDIF
      call deletedatalayout(tempk) 
      call deletedatalayout(Fieldfil1) 
      call deletedatalayout(Fieldfil2) 
      call deletedatalayout(Fieldfil3) 
    ELSEIF (U1%nx .GT. U2%nx .AND. (U1%nx .EQ. U1%ny) .AND. (U1%nx .EQ. U1%nz) ) then 
      IF (.NOT. copyStructonly(U2,Fieldfil1)) RETURN
      IF (.NOT. copyStructonly(V2,Fieldfil2)) RETURN
      IF (.NOT. copyStructonly(W2,Fieldfil3)) RETURN
      IF (.NOT. initDataLayout("tempK", tempK,(U1%nx/2)+1,U1%ny,U1%nz, &
               & U1%Lx,U1%Ly,U1%Lz,nbcpus,spec_rank,alongZ)) THEN
         WRITE(6,'(a,i0,a)')'[ERROR] in compareVeloFields : not enouhgt memory'
         RETURN
      ENDIF
      CALL ftran(U1,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,U2,Fieldfil1) ) RETURN 
      CALL ftran(V1,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,V2,Fieldfil2) ) RETURN 
      CALL ftran(W1,tempk,res)
      IF (.NOT. res) RETURN
      IF (.NOT. transdiscretization (tempK,W2,Fieldfil3) ) RETURN 
      IF (equalToVal(Fieldfil1,U2) .AND. equalToVal(Fieldfil2,V2) .AND. equalToVal(Fieldfil3,W2)) THEN
         IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : BOTH VECTOR FIELDS ARE EQUAL'
      ELSE
         IF (.NOT. equalToVal(Fieldfil1,U2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : U1transdiscretize != U2' 
         ENDIF
         IF (.NOT. equalToVal(Fieldfil2,V2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : V1transdiscretize != V2' 
         ENDIF
         IF (.NOT. equalToVal(Fieldfil3,W2)) THEN
           IF (spec_rank .eq. 0) write(6,'(a,i0,a)')'[INFO ] in compareVeloFields : W1transdiscretize != W2' 
         ENDIF
      ENDIF
      call deletedatalayout(tempk) 
      call deletedatalayout(Fieldfil1) 
      call deletedatalayout(Fieldfil2) 
      call deletedatalayout(Fieldfil3) 
    ELSE
         write(6,'(a,i0,a)')'[ERROR] in compareVeloFields, case not implemented'
         return
    ENDIF
    success = .true.

end subroutine compareVeloFields



!------------------------------------------------------------------------------
!> Show min, max, average and variance of a scalar fields with adding of filterd
!> quantities for DNS scalar resolution and comparison with LES scalar resolution
!> TAKE CARE IT IS ASSUMED THAT DNS RESOLUTION IS FOR SCALAR 1 AND LES FOR OTHER
!! @author = Patrick BEGOU, LEGI
!!    @param[in]    scalarray   = arrays of scalars in physical space
!!    @param[in]    nbscal      = the number of scalars to process
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus      = number of process in the pool
!!    @param[in]    tag         = optional, name of the array of scalar
!!                                  (spectral or particles scalars ?)
!------------------------------------------------------------------------------
subroutine showPostScalar(scalarray,nbscal,sim_time,ite,show_scalar,dime,n_iter,spec_rank,nbcpus)


    !I/O data
    type(REAL_DATA_LAYOUT), intent(in), dimension(:)    :: scalarray
    real(WP), intent(in)                                :: sim_time
    integer, intent(in)                                 :: nbscal, spec_rank, dime, ite, show_scalar, n_iter, nbcpus

    type(REAL_DATA_LAYOUT)                              :: scalfilDNS, scalerror
    type(COMPLEX_DATA_LAYOUT)                           :: tempK

    logical                                 :: res
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file


    integer :: sca
    real(WP)::mini, maxi, avg, var

    ! DNS scalar filtered on LES resolution
    ! tempK is DNS scalar in spectral space
    if (.not. initDataLayout("TempK", tempK,(ScalArray(1)%nx/2)+1,ScalArray(1)%ny,ScalArray(1)%nz, &
               & ScalArray(1)%Lx,ScalArray(1)%Ly,ScalArray(1)%Lz,nbcpus,spec_rank,alongZ)) then
         write(6,'(a,i0,a)')'[ERROR] initDataLayout for TYPE (COMPLEX_DATA_LAYOUT) : not enought memory!'
         return
    endif
    IF (.NOT. samelayout(ScalArray(1),ScalArray(2))) THEN
            IF (spec_rank .eq. 0) WRITE (6,'(a)') '    [INFO show scalar] First and second scalar are different.'
            IF (.NOT. copyStructOnly(ScalArray(2),scalfilDNS) ) THEN
                IF (spec_rank .eq. 0) WRITE (6,'(a)') '[ERROR] not enought memory!'
                RETURN
            ELSE
                call ftran(ScalArray(1),tempK,res)
                if (.not. res ) return
                IF (.NOT. transdiscretization (tempK,ScalArray(2),scalfilDNS) ) THEN
                  if (spec_rank .eq. 0) WRITE (6,'(a)') '[ERROR] in transdiscretization.' 
                  return
                ENDIF   
            ENDIF
    END IF
    call deletedatalayout(tempK)
    !IF (sim_time.eq.0.0) then
    IF (ite .eq. 0) then
       !IF (spec_rank .eq. 0) WRITE (6,'(a)') '[INFO show scalar] Initializing InitScaField and InitScaDNSField' 
       IF (.NOT. copyStructOnly(scalfilDNS, InitScaField) )THEN
         WRITE (6,'(a)') '[ERROR] in showPostScalar : not able to create InitScaField'
         RETURN
       ENDIF
       InitScaField%values = scalfilDNS%values
       IF (.NOT. copyStructOnly(ScalArray(1), InitScaDNSField) ) THEN
         WRITE (6,'(a)') '[ERROR] in showPostScalar : not able to create InitScaDNSField'
         RETURN
       ENDIF
       InitScaDNSField%values = ScalArray(1)%values
    END IF
    mini = computeFieldMin( InitScaField , spec_rank )
    maxi = computeFieldMax( InitScaField , spec_rank )
    avg  = computeFieldAvg( InitScaField , spec_rank )
    var  = computeFieldVar( InitScaField , spec_rank )
    if (spec_rank.EQ.0) then
       file_id = iopen()
       write(file_name,'(a)') 'InitfilDNS.stat'
       file_name = trim(adjustl(file_name))
       open(unit=file_id, file=file_name, form="FORMATTED",position="append")
          write(file_id,'(5(f25.10))' ) sim_time,mini, maxi, avg, var
       close(file_id)
       file_id = iclose(file_id)
    end if
    mini = computeFieldMin( InitScaDNSField , spec_rank )
    maxi = computeFieldMax( InitScaDNSField , spec_rank )
    avg  = computeFieldAvg( InitScaDNSField , spec_rank )
    var  = computeFieldVar( InitScaDNSField , spec_rank )
    if (spec_rank.EQ.0) then
       file_id = iopen()
       write(file_name,'(a)') 'InitDNS.stat'
       file_name = trim(adjustl(file_name))
       open(unit=file_id, file=file_name, form="FORMATTED",position="append")
          write(file_id,'(5(f25.10))' ) sim_time,mini, maxi, avg, var
       close(file_id)
       file_id = iclose(file_id)
    end if

    do sca=2,nbscal
         if (.not. samelayout(InitScaField,ScalArray(sca))) then
           write (6,'(a)') '[ERROR] InitScaField and one of the scalar are not the same'
           return
         endif
         res = computeCorr(var,InitScaField,ScalArray(sca),spec_rank)
         if (spec_rank.EQ.0) then
           file_id = iopen()
           write(file_name,'(a,i3.3,a)') 'TimeCorr_on_',sca,'.stat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
              write(file_id, '(3(f25.10))') sim_time, var
           close(file_id)
           file_id = iclose(file_id)
        end if
    end do

    res = computeCorr(var,InitScaField,ScalfilDNS,spec_rank)
    if (spec_rank.EQ.0) then
      file_id = iopen()
      write(file_name,'(a)') 'TimeCorr_filDNS.stat'
      file_name = trim(adjustl(file_name))
      open(unit=file_id, file=file_name, form="FORMATTED",position="append")
         write(file_id, '(3(f25.10))') sim_time, var
      close(file_id)
      file_id = iclose(file_id)
    end if

    res = computeCorr(var,InitScaDNSField,ScalArray(1),spec_rank)
    if (spec_rank.EQ.0) then
      file_id = iopen()
      write(file_name,'(a)') 'TimeCorr_DNS.stat'
      file_name = trim(adjustl(file_name))
      open(unit=file_id, file=file_name, form="FORMATTED",position="append")
         write(file_id, '(3(f25.10))') sim_time, var
      close(file_id)
      file_id = iclose(file_id)
    end if
    
    do sca=1,nbscal

        mini = computeFieldMin( scalarray(sca) , spec_rank )
        maxi = computeFieldMax( scalarray(sca) , spec_rank )
        avg  = computeFieldAvg( scalarray(sca) , spec_rank )
        var  = computeFieldVar( scalarray(sca) , spec_rank )
        
        if (spec_rank.EQ.0) then
           file_id = iopen()
           write(file_name,'(a,i3.3,a)') 'Scalar_',sca,'.stat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
              write(file_id, '(5(f25.10))') sim_time,mini, maxi, avg, var
           close(file_id)
           file_id = iclose(file_id)
        end if
        
    end do
    
    IF (.NOT. copyStructOnly(scalfilDNS, scalerror) ) RETURN
    do sca=2,nbscal
         if (.not. samelayout(ScalfilDNS,ScalArray(sca))) then
           write (6,'(a)') '[ERROR] ScalfilDNS and one of the scalar are not the same'
           return
         endif
         scalerror%values = (ScalfilDNS%values - ScalArray(sca)%values)**2.0
         avg = computeFieldAvg( scalerror , spec_rank )
         res = computeCorr(var,ScalfilDNS,ScalArray(sca),spec_rank)
         if (spec_rank.EQ.0) then
           file_id = iopen()
           write(file_name,'(a,i3.3,a)') 'ErCor_on_',sca,'.stat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
              write(file_id, '(3(f25.10))') sim_time, avg, var
           close(file_id)
           file_id = iclose(file_id)
        end if
    end do
    
    mini = computeFieldMin( ScalfilDNS , spec_rank )
    maxi = computeFieldMax( ScalfilDNS , spec_rank )
    avg  = computeFieldAvg( ScalfilDNS , spec_rank )
    var  = computeFieldVar( ScalfilDNS , spec_rank )
    if (spec_rank.EQ.0) then
       file_id = iopen()
       write(file_name,'(a)') 'Filt_DNS.stat'
       file_name = trim(adjustl(file_name))
       open(unit=file_id, file=file_name, form="FORMATTED",position="append")
          write(file_id,'(5(f25.10))' ) sim_time,mini, maxi, avg, var
       close(file_id)
       file_id = iclose(file_id)
    end if
    
    ! Show info of scalars 
    if (show_scalar>0 .and. mod(ite, show_scalar)== 0) then
        scalfilDNS%name = "Filt_DNS"
        if (.not.computeFieldPDF(ite,scalfilDNS,dime,spec_rank,n_iter)) then
           write(*,'(a)') '[WARNING] Failed to compute PDF of scalars fiels'
          return
        end if
    end if
    
    call deletedatalayout(scalfilDNS) 
    call deletedatalayout(scalerror) 

end subroutine showPostScalar


!-------------------------------------------------------------------------------------------
!> Save field at avs format in order to be able to restart a stopped simulation.
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    nbscal          = number scalar fields solved with pseudo-spectral method
!!    @param[in]    scalArray       = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    scalArray_part  = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    nbscal_part     = number scalar fields solved with pseudo-spectral/particle
!!                                  method
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along Z
!!    @param[in]    B               = magentic field (3D vector)
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    paraview_3D     = true for 3D output
!!    @param[in]    paraview_slice  = true for 2D slice of the 3D fields
!!    @return       success         = 1 if everything is right !
!! @author  Antoine Vollant
!-------------------------------------------------------------------------------------------

function paravieWrite(spec_rank,nbscal, scalArray, nbscal_part,scal_partArray, U, V, &
                      & W,B,imhd,tagScal,tagScalPart,tagVel,tagB, paraview_slice,    &
                      & paraview_3D) result(success)

   use vtkxml_plus

   !I/O data
   integer,intent(in)                                       :: spec_rank 
   integer,intent(in)                                       :: nbscal
   integer,intent(in)                                       :: nbscal_part
   type(real_data_layout),intent(in)                        :: U, V ,W
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:) :: B
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:) :: scalArray
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:) :: scal_partArray
! XXX Todo
!   type(real_vector_data_layout),intent(in)                :: scal_partVector
!   integer,dimension(:), intent(inout)                      :: tagScalPartVector
! XXX End Todo
   integer,dimension(:), intent(inout)                      :: tagScal
   integer,dimension(:), intent(inout)                      :: tagScalPart
   integer,dimension(3), intent(inout)                      :: tagVel
   integer,dimension(3), intent(inout)                      :: tagB
   logical,intent(in)                                       :: imhd
   logical,intent(in)                                       :: paraview_slice, paraview_3D
   logical                                                  :: success 


   !Local data
   integer                                                  :: i
! XXX Todo
!   type(real_data_layout)                                   :: ScalPart
! XXX End Todo


   success = .false.
   if(spec_rank==0) write(6,'(a)')'  [INFO-IO] Writing vtkxml output (paraview output)'

   if (nbscal .gt. 0) then
     do i=1,nbscal
        if(paraview_3D) then
            call vtkxml_write(tagScal(i), scalArray(i)%values , scalArray(i)%name )
            !if(spec_rank .eq. 0) write (6,'(a,1x,i0)') '  [INFO] Writing PARAVIEW field for scalar',i
        end if
        if(paraview_slice) call vtkxml_write_slice(tagScal(i), scalArray(i)%values, scalArray(i)%name )
     enddo
   endif
   if (nbscal_part .gt. 0) then
     do i=1,nbscal_part
        if(paraview_3D) then
            call vtkxml_write(tagScalPart(i), scal_partArray(i)%values , scal_partArray(i)%name )
            !if(spec_rank .eq. 0) write (6,'(a,1x,i0)') '  [INFO] Writing PARAVIEW field for scalar part',i
        end if
        if(paraview_slice) call vtkxml_write_slice(tagScalPart(i), scal_partArray(i)%values, scal_partArray(i)%name )

     enddo
   endif
   if (imhd) then
     do i=1,3
        if(paraview_3D) then
            call vtkxml_write(tagB(i), B(i)%values , B(i)%name )
            !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field B(',i,')'
        end if
        if(paraview_slice) call vtkxml_write_slice(tagB(i), B(i)%values , B(i)%name)
     enddo
   endif


   if(paraview_3D) then
        call vtkxml_write(tagVel(1), U%values , U%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for U'
        call vtkxml_write(tagVel(2), V%values , V%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for V'
        call vtkxml_write(tagVel(3), W%values , W%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for W'
   end if
   if(paraview_slice) then
        call vtkxml_write_slice(tagVel(1), U%values , U%name )
        call vtkxml_write_slice(tagVel(2), V%values , V%name )
        call vtkxml_write_slice(tagVel(3), W%values , W%name )
   end if

! XXX Todo
!   if (Scal_partVector%nb_components .gt. 0) then
!     do i=1, Scal_partVector%nb_components
!        call real_vector_manipulate_comp(Scal_partVector, i, ScalPart, guru=.true.)
!        if(paraview_3D) then
!            call vtkxml_write(tagScalPartVector(i), scalPart%values , scalPart%name)
!            !if(spec_rank .eq. 0) write (6,'(a,1x,i0)') '  [INFO] Writing PARAVIEW field for scalar part',i
!        end if
!        if(paraview_slice) call vtkxml_write_slice(tagScalPartVector(i), scalPart%values, scalPart%name )
!     enddo
!   endif
! XXX End Todo

   success = .true.

end function paravieWrite

!-------------------------------------------------------------------------------------------
!> Save field at avs format in order to be able to restart a stopped simulation.
!!  @author  
!!  Luca MARRADI LIPhy 2014
!!
!!  @param[in]    spec_rank       = rank inside the spectral communicator
!!  @param[in]    sigma           = stress tensor component XY
!!  @param[in]    paraview_3D     = true for 3D output
!!  @param[in]    paraview_slice  = true for 2D slice of the 3D fields
!!  @return       success         = 1 if everything is right !
!!
!!  @detalis
!> Save field in avs format in order to be able to restart at stopped simulation.
!-------------------------------------------------------------------------------------------

FUNCTION paravieWrite_luca(spec_rank,sigma,tagVel,paraview_slice,paraview_3D) result(success)

   USE vtkxml_plus

   !== INPUT/OUTPUT DATA ==
   integer,intent(in)                                       :: spec_rank 
   type(real_data_layout),intent(in)                        :: sigma
   integer,dimension(3), intent(inout)                      :: tagVel
   logical,intent(in)                                       :: paraview_slice, paraview_3D
   logical                                                  :: success 


   !== LOCAL DATA ==
   integer                                                  :: i

   success = .FALSE.
   IF(spec_rank==0) write(6,'(a)')'  [INFO-IO] Writing vtkxml output (paraview output)'

   !== PRINTS PARAVIEW SIGMA FILE ==
   ! Note tagVel represents the number of components of our stress tensor !!!
   if(paraview_3D) then
        call vtkxml_write(tagVel(1), sigma%values , sigma%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for U'
!!        call vtkxml_write(tagVel(2), V%values , V%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for V'
!!        call vtkxml_write(tagVel(3), W%values , W%name )
        !if(spec_rank .eq. 0) write (6,'(a,i0,a)') '  [INFO] Writing PARAVIEW field for W'
   end if
   if(paraview_slice) then
        call vtkxml_write_slice(tagVel(1), sigma%values , sigma%name )
!!        call vtkxml_write_slice(tagVel(2), V%values , V%name )
!!        call vtkxml_write_slice(tagVel(3), W%values , W%name )
   end if

   success = .TRUE.
END FUNCTION  paravieWrite_luca

!-------------------------------------------------------------------------------------------
!> Save field at avs format in order to be able to restart a stopped simulation.
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    scalArray       = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    scalArray_part  = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along Z
!!    @param[in]    B               = magentic field (3D vector)
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    iter            = number of the current iteration
!!    @param[in]    iterMax         = total number of iterations
!!    @param[in]    usename         = is an optional logical which allow to use the member "name" of
!!                                    data_layout structure in the name of the dump field (useful for transdiscretization) 
!!    @return       success         = 1 if everything is right !
!! @author  Antoine Vollant
!-------------------------------------------------------------------------------------------
 function avsWrite(spec_rank,scalArray, scal_partArray, U, V, &
                      & W,B,imhd,iter,iterMax,usename) result(success)


   !I/O data
   integer,intent(in)                                     :: spec_rank
   type(real_data_layout),intent(in)                   :: U, V ,W
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: B
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scalArray
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scal_partArray
   logical,intent(in)                                     :: imhd
   logical                                                :: success
   integer,intent(in)                                     :: iter
   integer,intent(in)                                     :: iterMax
   logical,intent(in),optional                            :: usename


   !Local data
   integer             :: j
   CHARACTER(LEN=64)   :: writecom
   CHARACTER(LEN=3)    :: scalcom
   logical             :: withname

   success = .false.
   if(spec_rank==0) write(6,'(a)')'  [INFO-IO] Writing AVS output'

      IF(iterMax.LT.1000) THEN
          WRITE(writecom,'(i3.3)') iter
      ELSEIF(iterMax.LT.1000000) THEN
          WRITE(writecom,'(i6.6)') iter
      ELSE
          WRITE(writecom,'(i0)') iter
      ENDIF

      if ( present(usename))then
        if (usename)then
          withname = .true.
        else
          withname = .false.
        endif
      else
        withname = .false.
      endif

      ! save U, V, W
      if ( withname ) then
        IF (.NOT. dump_UVWP_ForAVS(trim(U%name)//'_'//trim(writecom),spec_rank,U%lx,U%ly,U%lz,U,usename=.true.) .or. &
            .not. dump_UVWP_ForAVS(trim(V%name)//'_'//trim(writecom),spec_rank,V%lx,V%ly,V%lz,V,usename=.true.) .or. & 
            .not. dump_UVWP_ForAVS(trim(W%name)//'_'//trim(writecom),spec_rank,W%lx,W%ly,W%lz,W,usename=.true.)) THEN 
            WRITE (6,'(a)') 'FAILED in velocity field dump'
            RETURN
        ENDIF
      else
        IF (.NOT. dump_UVWP_ForAVS('Fields_at_'//trim(writecom),spec_rank,U%lx,U%ly,U%lz,U=U,V=V,W=W)) THEN 
            WRITE (6,'(a)') 'FAILED in velocity field dump'
           RETURN
        ENDIF
      endif
      ! save scalar - pseudo-spectral solver
      IF(ASSOCIATED(ScalArray)) THEN
          DO j=1,SIZE(ScalArray)
              WRITE(scalcom,'(i3.3)') j
              if ( withname ) then
                IF (.NOT. dump_UVWP_ForAVS(TRIM(scalArray(j)%name)//'_'//scalcom//'_at_'//TRIM(writecom), &
                   & spec_rank,scalArray(j)%lx,scalArray(j)%ly,scalArray(j)%lz,U=ScalArray(j),usename=.true.)) THEN
                    WRITE (6,'(a,i0)') 'FAILED in scalar field dump',j
                    RETURN
                ENDIF
              else
                IF (.NOT. dump_UVWP_ForAVS('Scalar_'//scalcom//'_at_'//TRIM(writecom), &
                   & spec_rank,scalArray(j)%lx,scalArray(j)%ly,scalArray(j)%lz,U=ScalArray(j))) THEN
                    WRITE (6,'(a,i0)') 'FAILED in scalar field dump',j
                    RETURN
                ENDIF
              endif
          ENDDO
      ENDIF  
      ! save scalar - spectral/particles solver
      IF(ASSOCIATED(Scal_partArray)) THEN 
          DO j=1,SIZE(Scal_partArray)
              WRITE(scalcom,'(i3.3)') j
              if ( withname ) then
                IF (.NOT. dump_UVWP_ForAVS(TRIM(Scal_partArray(j)%name)//'_'//scalcom//'_at_'//TRIM(writecom), &
                   & spec_rank,Scal_partArray(j)%lx,Scal_partArray(j)%ly,Scal_partArray(j)%lz,&
                   &U=Scal_partArray(j),usename=.true.)) THEN
                    WRITE (6,'(a,i0)') 'FAILED in scalar field dump',j
                    RETURN
                ENDIF
              else
                IF (.NOT. dump_UVWP_ForAVS('ScalarPart_'//scalcom//'_at_'//TRIM(writecom), &
                   & spec_rank,Scal_partArray(j)%lx,Scal_partArray(j)%ly,Scal_partArray(j)%lz,&
                   &U=Scal_partArray(j),usename=.true.)) THEN
                    WRITE (6,'(a,i0)') 'FAILED in scalar field dump',j
                    RETURN
                ENDIF
              endif
          ENDDO
      ENDIF  
      ! save mhd
      IF (imhd) THEN
          if (withname) then
            do j=1,3
              IF (.NOT. dump_UVWP_ForAVS(TRIM(B(j)%name)//'_'//TRIM(writecom), &
                 & spec_rank,B(j)%lx,B(j)%ly,B(j)%lz,U=B(j),usename=.true.)) THEN 
                   WRITE (6,'(a)') 'FAILED in B dump'
                   RETURN
              ENDIF
            enddo
          else
            IF (.NOT. dump_UVWP_ForAVS('Bx_at_'//TRIM(writecom), &
               & spec_rank,B(1)%lx,B(1)%ly,B(1)%lz,U=B(1))) THEN 
                 WRITE (6,'(a)') 'FAILED in Bx dump'
                 RETURN
            ENDIF
            IF (.NOT. dump_UVWP_ForAVS('By_at_'//TRIM(writecom), &
               & spec_rank,B(2)%lx,B(2)%ly,B(2)%lz,U=B(2))) THEN 
                 WRITE (6,'(a)') 'FAILED in By dump'
                 RETURN
            ENDIF
            IF (.NOT. dump_UVWP_ForAVS('Bz_at_'//TRIM(writecom), &
               & spec_rank,B(3)%lx,B(3)%ly,B(3)%lz,U=B(3))) THEN
                 WRITE (6,'(a)') 'FAILED in Bz dump'
                 RETURN
            ENDIF
          endif
      ENDIF

   success = .true.

 end function avsWrite

!-------------------------------------------------------------------------------------------
!> Save field at avs format in order to be able to restart a stopped simulation.
!!  @param[in]    spec_rank   = rank inside the spectral communicator
!!  @param[in]    sigma       = stress tensor 
!!  @param[in]    iter        = number of the current iteration
!!  @param[in]    iterMax     = total number of iterations
!!  @param[in]    usename     = is an optional logical which allow to use the member "name" of
!!                                data_layout structure in the name of the dump field (useful for transdiscretization) 
!!  @return       success     = 1 if everything is right !
!!@author  Luca MARRADI LIPhy 2014
!-------------------------------------------------------------------------------------------
 function avsWrite_luca(spec_rank,sigma,iter,iterMax,usename) result(success)


   !== INPUT/OUTPUT DATA ==
   IMPLICIT NONE

   integer,intent(in)                           :: spec_rank
   type(real_data_layout),intent(in)            :: sigma
   logical                                      :: success
   integer,intent(in)                           :: iter
   integer,intent(in)                           :: iterMax
   logical,intent(in),optional                  :: usename


   !== LOCAL DATA ==
   integer             :: j
   CHARACTER(LEN=64)   :: writecom
   CHARACTER(LEN=3)    :: scalcom
   logical             :: withname

   success = .FALSE.
   if(spec_rank==0) write(6,'(a)')'  [INFO-IO] Writing AVS output'

      IF(iterMax.LT.1000) THEN
          WRITE(writecom,'(i3.3)') iter
      ELSEIF(iterMax.LT.1000000) THEN
          WRITE(writecom,'(i6.6)') iter
      ELSE
          WRITE(writecom,'(i0)') iter
      ENDIF

      if ( present(usename))then
        if (usename)then
          withname = .TRUE.
        else
          withname = .FALSE.
        endif
      else
          withname = .FALSE.
      endif

      !== SAVE SIGMA AND ACTIVITY ==
      ! ACTIVITY SHOULD BE IMPLEMENTED
      IF( withname ) THEN
        if ( .NOT.dump_UVWP_ForAVS_luca(trim(sigma%name)//'_'//trim(writecom),spec_rank,sigma%lx,sigma%ly,sigma%lz,sigma,usename=.true.) ) then
            write (6,'(a)') 'FAILED in velocity field dump'
            return
        endif
      ELSE
        if ( .NOT.dump_UVWP_ForAVS_luca('Fields_at_'//trim(writecom),spec_rank,sigma%lx,sigma%ly,sigma%lz,sigma=sigma) ) then
            write (6,'(a)') 'FAILED in velocity field dump'
           return
        endif
      ENDIF
   success = .TRUE.
 END FUNCTION avsWrite_luca

!-------------------------------------------------------------------------------------------
!> Save field at avs format in order to be able to restart a stopped simulation.
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    scalArray       = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    scalArray_part  = array containing all the scalar fields (scalars solved
!!                                  with pseudo-spectral methods)
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along Z
!!    @param[in]    B               = magentic field (3D vector)
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    iter            = number of the current iteration
!!    @param[in]    iterMax         = total number of iterations
!!    @param[in]    usename         = is an optional logical which allow to use the member "name" of
!!                                    data_layout structure in the name of the dump field (useful for transdiscretization)
!!    @return       success         = 1 if everything is right !
!! @author  Antoine Vollant and Jean-Baptiste Lagaert
!-------------------------------------------------------------------------------------------
 function dataWrite(spec_rank,scalArray, scal_partArray, U, V, &
                      & W,B,imhd,iter,postname) result(success)

    use io_interface
    use mpilayout_tools , only : getnbcpus
    use parser_tools ,    only : GetNxNyNz
    use parser_tools ,    only : parser_read
    use parser_tools ,    only : parser_is_defined

    !I/O data
    integer,intent(in)                                          :: spec_rank
    type(real_data_layout),intent(in)                           :: U, V ,W
    type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: B
    type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scalArray
    type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scal_partArray
    logical,intent(in)                                          :: imhd
    logical                                                     :: success
    integer,intent(in)                                          :: iter
    character(len=*),intent(in),optional                        :: postname


    !Local data
    integer                :: j,nx,ny,nz,nbcpus
    type(REAL_DATA_LAYOUT) :: Sca_trans
    character(len=str_long):: request

    success = .false.
    nbcpus = getnbcpus()
    if(present(postname)) then
        if(spec_rank==0) write(6,'(a, i0)')'  [INFO-IO] Writing output '//trim(postname)// ' - indice = ', iter

        ! save U, V, W
        ! TODO : definir un write_vector et l'utiliser ! (en fait il faut surtout
        ! tester la version hdf5 et après y a plus qu'à)
        if(     .not. write_datalayout_scalar(U,iter,spec_rank, postname)    .or. &
                .not. write_datalayout_scalar(V,iter,spec_rank, postname)    .or. &
                .not. write_datalayout_scalar(W,iter,spec_rank, postname)    ) then
            write (6,'(a)') 'FAILED to write output for velocity field'
            return
        end if

        ! save scalar - pseudo-spectral solver
        IF(ASSOCIATED(ScalArray)) THEN
            DO j=1,SIZE(ScalArray)
                if(.not. write_datalayout_scalar(scalArray(j),iter,spec_rank, postname)) then
                    write (6,'(a,i0)') 'FAILED to write output for scalar ', j
                    return
                end if
            end do
        end if
        ! save scalar - spectral/particles solver
        IF(ASSOCIATED(Scal_partArray)) then
            DO j=1,SIZE(Scal_partArray)
                if(.not. write_datalayout_scalar(scal_partArray(j),iter,spec_rank, postname)) then
                   write (6,'(a,i0)') 'FAILED to write output for scalar_part ', j
                   return
                end if
            end do
        end if
        ! save mhd
        IF (imhd) then
              do j=1,3
                if(.not. write_datalayout_scalar(B(j),iter,spec_rank, postname)) then
                    write (6,'(a,i0)') 'FAILED to write output for B ', j
                    return
                end if
            end do
        end if
    else
        if(spec_rank==0) write(6,'(a)')'  [INFO-IO] Writing output (avs or hdf5)'

        ! save U, V, W
        ! TODO : definir un write_vector et l'utiliser ! (en fait il faut surtout
        ! tester la version hdf5 et après y a plus qu'à)
        if(     .not. write_datalayout_scalar(U,iter,spec_rank)    .or. &
                .not. write_datalayout_scalar(V,iter,spec_rank)    .or. &
                .not. write_datalayout_scalar(W,iter,spec_rank)    ) then
            write (6,'(a)') 'FAILED to write output for velocity field'
            return
        end if

        ! save scalar - pseudo-spectral solver
        IF(ASSOCIATED(ScalArray)) THEN
            DO j=1,SIZE(ScalArray)
                write (request,'(a,1x,i0)') 'Discretization to write Scalar',j
                if ( parser_is_defined(request)) then
                    if (.not. GetNxNyNz(request,nx,ny,nz)) return
                    if (.not. transdiscretizeField(ScalArray(j),Sca_trans,spec_rank,nbcpus,nx,ny,nz)) then
                        WRITE (6,'(a,1x,i0)') 'FAILED in transdiscretize Scalar',j 
                        RETURN
                    endif
                    Sca_trans%name=trim(adjustl(ScalArray(j)%name))//'_trans'
                    if(.not. write_datalayout_scalar(Sca_trans,iter,spec_rank)) then
                       write (6,'(a,i0)') 'FAILED to write output for scalar_part ', j
                       return
                    end if
                    call deletedatalayout(Sca_trans)
                else    
                    if(.not. write_datalayout_scalar(scalArray(j),iter,spec_rank)) then
                       write (6,'(a,i0)') 'FAILED to write output for scalar ', j
                       return
                    end if
                endif
            end do
        end if
        ! save scalar - spectral/particles solver
        IF(ASSOCIATED(Scal_partArray)) then
            DO j=1,SIZE(Scal_partArray)
                write (request,'(a,1x,i0)') 'Discretization to write Scalar_part',j
                if ( parser_is_defined(request)) then
                    if (.not. GetNxNyNz(request,nx,ny,nz)) return
                    if (.not. transdiscretizeField(scal_partArray(j),Sca_trans,spec_rank,nbcpus,nx,ny,nz)) then
                        WRITE (6,'(a,1x,i0)') 'FAILED in transdiscretize scal_partArray',j 
                        RETURN
                    endif
                    Sca_trans%name=trim(adjustl(scal_partArray(j)%name))//'_trans'
                    if(.not. write_datalayout_scalar(Sca_trans,iter,spec_rank)) then
                       write (6,'(a,i0)') 'FAILED to write output for scalar_part ', j
                       return
                    end if
                    call deletedatalayout(Sca_trans)
                else    
                    if(.not. write_datalayout_scalar(scal_partArray(j),iter,spec_rank)) then
                       write (6,'(a,i0)') 'FAILED to write output for scalar_part ', j
                       return
                    end if
                endif
            end do
        end if
        ! save mhd
        IF (imhd) then
              do j=1,3
                if(.not. write_datalayout_scalar(B(j),iter,spec_rank)) then
                    write (6,'(a,i0)') 'FAILED to write output for B ', j
                    return
                end if
            end do
        end if
    end if

    success = .true.

end function dataWrite

!-------------------------------------------------------------------------------------------
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    sigma           = stress tensor 
!!    @param[in]    iter            = number of the current iteration
!!    @param[in]    iterMax         = total number of iterations
!!    @param[in]    usename         = is an optional logical which allow to use the member "name" of
!!                                    data_layout structure in the name of the dump field (useful for transdiscretization)
!!    @return       success         = 1 if everything is right !
!! @author  Luca MARRADI LIPhy 2014
!! @details
!! Save field at avs format in order to be able to restart at stopped simulation
!-------------------------------------------------------------------------------------------
FUNCTION dataWrite_luca(spec_rank,sigma,iter,postname) result(success)

    use io_interface
    use mpilayout_tools , only : getnbcpus
    use parser_tools ,    only : GetNxNyNz
    use parser_tools ,    only : parser_read
    use parser_tools ,    only : parser_is_defined

    IMPLICIT NONE

    !== INPUT/OUTPUT DATA
    integer,intent(in)                                          :: spec_rank
    type(real_data_layout),intent(in)                           :: sigma
    logical                                                     :: success
    integer,intent(in)                                          :: iter
    character(len=*),intent(in),optional                        :: postname


    !LOCAL DATA
    integer                :: j,nx,ny,nz,nbcpus
    character(len=str_long):: request

    success = .FALSE.
    nbcpus = getnbcpus()
    IF(present(postname)) THEN
        IF(spec_rank==0) write(6,'(a, i0)')' [INFO-IO] Writing output '//trim(postname)// ' - indice = ', iter

        !== SAVE SIGMA AND ACTIVITY ==
        IF( .NOT. write_datalayout_scalar(sigma,iter,spec_rank, postname) ) THEN
            WRITE (6,'(a)') 'FAILED to write output for velocity field'
            RETURN
        ENDIF

    ELSE
       IF(spec_rank==0) write(6,'(a)')' [INFO-IO] Writing output (avs or hdf5)'

        !== SAVE SIGMA AND ACTIVITY == 
        ! tester la version hdf5 et après y a plus qu'�| )
        IF( .NOT.write_datalayout_scalar(sigma,iter,spec_rank) ) THEN
            WRITE (6,'(a)') 'FAILED to write output for velocity field'
            RETURN
        ENDIF

    ENDIF

    success=.TRUE.
    RETURN

END FUNCTION dataWrite_luca


!------------------------------------------------------------------------------
!> Show min and max and average of the velocity fields.
!! @author = Patrick BEGOU, LEGI
!!    @param[in]    U           = longitudinal velocity

!------------------------------------------------------------------------------
!> Show min and max and average of the velocity fields.
!! @author = Patrick BEGOU, LEGI
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @param[in]    nbcpus      = number of process in the pool
!------------------------------------------------------------------------------
subroutine showVeloMinMaxAvg(U,V,W,spec_rank, nbcpus, sim_time)


    type(REAL_DATA_LAYOUT), intent(in)  ::U,V,W
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus
    type(REAL_DATA_LAYOUT)              :: buf

    real(WP)    ::minu, maxu, avgu, avgu2
    real(WP)    ::minv, maxv, avgv, avgv2
    real(WP)    ::minw, maxw, avgw, avgw2
    real(WP)    :: sim_time
    real(WP)    :: eps = 1e-10

    maxu =computeFieldMax(U ,spec_rank)
    minu =computeFieldMin(U ,spec_rank)
    avgu = computeFieldAvg(U , spec_rank , nbcpus )
    buf = U
    buf%values = buf%values * buf%values
    avgu2 = computeFieldAvg(buf , spec_rank , nbcpus )

    maxv =computeFieldMax(V ,spec_rank)
    minv =computeFieldMin(V ,spec_rank)
    avgv = computeFieldAvg(V , spec_rank , nbcpus )
    buf = V
    buf%values = buf%values * buf%values
    avgv2 = computeFieldAvg(buf , spec_rank , nbcpus )

    maxw =computeFieldMax(W ,spec_rank)
    minw =computeFieldMin(W ,spec_rank)
    avgw = computeFieldAvg(W , spec_rank , nbcpus )
    buf = W
    buf%values = buf%values * buf%values
    avgw2 = computeFieldAvg(buf , spec_rank , nbcpus )

    call deletedatalayout(buf)

    if(spec_rank.EQ.0) then
        write(6,'(3(a,g15.8))')'    [INFO-STAT] U min= ',minu,' max= ',maxu, ' average= ', avgu
        write(6,'(3(a,g15.8))')'    [INFO-STAT] V min= ',minv,' max= ',maxv, ' average= ', avgv
        write(6,'(3(a,g15.8))')'    [INFO-STAT] W min= ',minw,' max= ',maxw, ' average= ', avgw

        if ( sim_time .lt. eps ) then
            open(10,file='vel_info.out',form='formatted')
            write(10,'(a)') '#time, umin, umax, vmin, vmax, wmin, wmax, meanu, meanv, meanw, meanu2, meanv2, meanw2'
        else
            open(10,file='vel_info.out',form='formatted',position='append')
        endif
        write(10,'(13(g14.6))') sim_time, minu, maxu, minv, maxv, minw, maxw, avgu, avgv, avgw, avgu2, avgv2, avgw2
        close(10)

    end if

end subroutine showVeloMinMaxAvg


!-------------------------------------------------------------------------------------------
!> Perform custom post-process and show custom information instead of the run
!! @author Guillaume Balarac 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    scalarray       = arrays of scalars in physical space
!!    @param[in]    nbscal          = the number of scalars to process
!!    @param[in]    scal_partArray  = arrays of scalars in physical space
!!    @param[in]    nbscal_part     = the number of scalars to process
!!    @param[in]    n_iter          = the number of max iteration
!!    @return       success         = logical equals to true is everything is right.
!-------------------------------------------------------------------------------------------
function postprocess_double(ncpus,spec_rank,U, V, W, B, imhd) result(success)


    ! Input/Output
    integer, intent(in)                                 :: ncpus
    integer, intent(in)                                 :: spec_rank
    type(REAL_DATA_LAYOUT), intent(inout)               ::U,V,W
    type(REAL_DATA_LAYOUT), intent(inout), dimension(:) ::B


    type(REAL_DATA_LAYOUT)                              ::Udb,Vdb,Wdb
    type(REAL_DATA_LAYOUT), dimension(3)                ::Bdb

    logical             :: success
    logical             :: imhd

    INTEGER             :: i, j, k
    INTEGER             :: nxdb,nydb,nzdb
    REAL(WP)            :: Lxdb,Lydb,Lzdb

    IF (ncpus.ne.1) THEN
       WRITE(6,*)'[ERROR] DOES NOT WORK IN PARELLEL'
       RETURN
       STOP
    ELSE

! To do a 2 times larger domain for Velocity
       nxdb = 2*U%nx
       nydb = 2*U%ny
       nzdb = 2*U%nz
       Lxdb = 2*U%Lx
       Lydb = 2*U%Ly
       Lzdb = 2*U%Lz
       IF (spec_rank.eq.0) WRITE(6,*)'Velocity: new nx, ny, nz:',nxdb,nydb,nzdb
       IF (spec_rank.eq.0) WRITE(6,*)'Velocity: new Lx, Ly, Lz:',Lxdb,Lydb,Lzdb
       IF (.NOT. initDataLayout("U",Udb,nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
          WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Udb on ', spec_rank
          RETURN
       ENDIF
       IF (.NOT. initDataLayout("V",Vdb,nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
          WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Vdb on ', spec_rank
          RETURN
       ENDIF
       IF (.NOT. initDataLayout("W",Wdb,nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
          WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Wdb on ', spec_rank
          RETURN
       ENDIF

       DO k = U%zmin, U%zmax
          DO j = U%ymin, U%ymax
             DO i = U%xmin, U%xmax
                Udb%values(i,j,k) = U%values(i,j,k)
                Vdb%values(i,j,k) = V%values(i,j,k)
                Wdb%values(i,j,k) = W%values(i,j,k)
                Udb%values(i+U%nx,j,k) = U%values(i,j,k)
                Vdb%values(i+U%nx,j,k) = V%values(i,j,k)
                Wdb%values(i+U%nx,j,k) = W%values(i,j,k)
                Udb%values(i,j+U%ny,k) = U%values(i,j,k)
                Vdb%values(i,j+U%ny,k) = V%values(i,j,k)
                Wdb%values(i,j+U%ny,k) = W%values(i,j,k)
                Udb%values(i,j,k+U%nz) = U%values(i,j,k)
                Vdb%values(i,j,k+U%nz) = V%values(i,j,k)
                Wdb%values(i,j,k+U%nz) = W%values(i,j,k)
                Udb%values(i+U%nx,j+U%ny,k) = U%values(i,j,k)
                Vdb%values(i+U%nx,j+U%ny,k) = V%values(i,j,k)
                Wdb%values(i+U%nx,j+U%ny,k) = W%values(i,j,k)
                Udb%values(i+U%nx,j,k+U%nz) = U%values(i,j,k)
                Vdb%values(i+U%nx,j,k+U%nz) = V%values(i,j,k)
                Wdb%values(i+U%nx,j,k+U%nz) = W%values(i,j,k)
                Udb%values(i,j+U%ny,k+U%nz) = U%values(i,j,k)
                Vdb%values(i,j+U%ny,k+U%nz) = V%values(i,j,k)
                Wdb%values(i,j+U%ny,k+U%nz) = W%values(i,j,k)
                Udb%values(i+U%nx,j+U%ny,k+U%nz) = U%values(i,j,k)
                Vdb%values(i+U%nx,j+U%ny,k+U%nz) = V%values(i,j,k)
                Wdb%values(i+U%nx,j+U%ny,k+U%nz) = W%values(i,j,k)
            ENDDO
         ENDDO
       ENDDO

       ! Dump new fields
       IF (.NOT. dump_UVWP_ForAVS('Larger-Fields',spec_rank,lxdb,lydb,lzdb,Udb,Vdb,Wdb)) STOP 'FAILED in avsdump'
       CALL deleteDataLayout(Udb)
       CALL deleteDataLayout(Vdb)
       CALL deleteDataLayout(Wdb)

! To do a 2 times larger domain for MHD 
       IF (imhd) THEN
          nxdb = 2*B(1)%nx
          nydb = 2*B(1)%ny
          nzdb = 2*B(1)%nz
          Lxdb = 2*B(1)%Lx
          Lydb = 2*B(1)%Ly
          Lzdb = 2*B(1)%Lz
          IF (spec_rank.eq.0) WRITE(6,*)'MHD: new nx, ny, nz:',nxdb,nydb,nzdb
          IF (spec_rank.eq.0) WRITE(6,*)'MHD: new Lx, Ly, Lz:',Lxdb,Lydb,Lzdb
          IF (.NOT. initDataLayout("Bx",Bdb(1),nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
             WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Bdb(1) on ', spec_rank
             RETURN
          ENDIF
          IF (.NOT. initDataLayout("By",Bdb(2),nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
             WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Bdb(2) on ', spec_rank
             RETURN
          ENDIF
          IF (.NOT. initDataLayout("Bz",Bdb(3),nxdb,nydb,nzdb,Lxdb,Lydb,Lzdb,ncpus,spec_rank)) THEN
             WRITE(6,'(a,i0)')'[ERROR] data_init Failed to init Bdb(3) on ', spec_rank
             RETURN
          ENDIF

          DO k = B(1)%zmin, B(1)%zmax
             DO j = B(1)%ymin, B(1)%ymax
                DO i = B(1)%xmin, B(1)%xmax
                   Bdb(1)%values(i,j,k) = B(1)%values(i,j,k)
                   Bdb(2)%values(i,j,k) = B(2)%values(i,j,k)
                   Bdb(3)%values(i,j,k) = B(3)%values(i,j,k)
                   Bdb(1)%values(i+B(1)%nx,j,k) = B(1)%values(i,j,k)
                   Bdb(2)%values(i+B(1)%nx,j,k) = B(2)%values(i,j,k)
                   Bdb(3)%values(i+B(1)%nx,j,k) = B(3)%values(i,j,k)
                   Bdb(1)%values(i,j+B(1)%ny,k) = B(1)%values(i,j,k)
                   Bdb(2)%values(i,j+B(1)%ny,k) = B(2)%values(i,j,k)
                   Bdb(3)%values(i,j+B(1)%ny,k) = B(3)%values(i,j,k)
                   Bdb(1)%values(i,j,k+B(1)%nz) = B(1)%values(i,j,k)
                   Bdb(2)%values(i,j,k+B(1)%nz) = B(2)%values(i,j,k)
                   Bdb(3)%values(i,j,k+B(1)%nz) = B(3)%values(i,j,k)
                   Bdb(1)%values(i+B(1)%nx,j+B(1)%ny,k) = B(1)%values(i,j,k)
                   Bdb(2)%values(i+B(1)%nx,j+B(1)%ny,k) = B(2)%values(i,j,k)
                   Bdb(3)%values(i+B(1)%nx,j+B(1)%ny,k) = B(3)%values(i,j,k)
                   Bdb(1)%values(i+B(1)%nx,j,k+B(1)%nz) = B(1)%values(i,j,k)
                   Bdb(2)%values(i+B(1)%nx,j,k+B(1)%nz) = B(2)%values(i,j,k)
                   Bdb(3)%values(i+B(1)%nx,j,k+B(1)%nz) = B(3)%values(i,j,k)
                   Bdb(1)%values(i,j+B(1)%ny,k+B(1)%nz) = B(1)%values(i,j,k)
                   Bdb(2)%values(i,j+B(1)%ny,k+B(1)%nz) = B(2)%values(i,j,k)
                   Bdb(3)%values(i,j+B(1)%ny,k+B(1)%nz) = B(3)%values(i,j,k)
                   Bdb(1)%values(i+B(1)%nx,j+B(1)%ny,k+B(1)%nz) = B(1)%values(i,j,k)
                   Bdb(2)%values(i+B(1)%nx,j+B(1)%ny,k+B(1)%nz) = B(2)%values(i,j,k)
                   Bdb(3)%values(i+B(1)%nx,j+B(1)%ny,k+B(1)%nz) = B(3)%values(i,j,k)
               ENDDO
            ENDDO
          ENDDO

          ! Dump new fields
          IF (.NOT. dump_UVWP_ForAVS('Larger-Bx-Fields',spec_rank,lxdb,lydb,lzdb,Bdb(1))) STOP 'FAILED in avsdump'
          IF (.NOT. dump_UVWP_ForAVS('Larger-By-Fields',spec_rank,lxdb,lydb,lzdb,Bdb(2))) STOP 'FAILED in avsdump'
          IF (.NOT. dump_UVWP_ForAVS('Larger-Bz-Fields',spec_rank,lxdb,lydb,lzdb,Bdb(3))) STOP 'FAILED in avsdump'
          CALL deleteDataLayout(Bdb(1))
          CALL deleteDataLayout(Bdb(2))
          CALL deleteDataLayout(Bdb(3))
       ENDIF

       success = .true.

     ENDIF

end function postprocess_double
 



!> Compute Q criterion and Vorticity vector and dump these quantities with paraview
!! @author Guillaume Balarac
!!    @return       success         = logical equals to true is everything is right.
function WriteQcOmega(U,Uk,Vk,Wk,VelWN,Qc,O1,O2,O3,tag_Qc,tag_O,spec_rank, iter) result(success)

    use io_interface


    ! Input/Output
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: Uk,Vk,Wk
    type(REAL_DATA_LAYOUT), intent(in)      :: U
    type(REAL_DATA_LAYOUT), intent(inout)   :: Qc,O1,O2,O3
    TYPE(WaveNumbers),INTENT(IN)            :: velWN
    integer, intent(inout)                  :: tag_Qc,tag_O(3)
    integer, intent(in)                     :: spec_rank
    integer, intent(in)                     :: iter
    logical                                 :: success

    ! Local varaibales
    type(REAL_DATA_LAYOUT)                  :: S11,S12,S13,S22,S23,S33
    character(len=str_short)                :: output_format, type_para

    success = .true.

    ! Compute Sij
    if (.not. copyStructOnly(U,S11) .or. &
       &.not. copyStructOnly(U,S12) .or. &
       &.not. copyStructOnly(U,S13) .or. &
       &.not. copyStructOnly(U,S22) .or. &
       &.not. copyStructOnly(U,S23) .or. &
       &.not. copyStructOnly(U,S33) )then
         success = .false.
         write(6,'(a,i0)')'[ERROR] in WriteQcOmega : not enought memory !',spec_rank
     endif
     call computeStressTensor(Uk,VK,WK,VelWN,S11,S12,S13,S22,S23,S33,success)

    ! Compute Oi
    if (.not. copyStructOnly(U,O1) .or. &
       &.not. copyStructOnly(U,O2) .or. &
       &.not. copyStructOnly(U,O3) )then
         success = .false.
         write(6,'(a,i0)')'[ERROR] in WriteQcOmega : not enought memory !',spec_rank
    endif
    call computeVorticity(UK,VK,WK,VelWN,O1,O2,O3,success)

    !compute Q criterion
    if (.not. copyStructOnly(U,Qc) ) then
         success = .false.
         write(6,'(a,i0)')'[ERROR] in WriteQcOmega : not enought memory !',spec_rank
    endif   
    Qc%values = 0.5_WP*( S11%values**2.0 + S22%values**2.0 + S33%values**2.0 + &
          &        2.0_WP * (S12%values**2.0 + S13%values**2.0 + S23%values**2.0) - &
          &        2.0_WP * (O1%values**2.0 + O2%values**2.0 + O3%values**2.0 ) )
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S22)
    call deleteDataLayout(S23)
    call deleteDataLayout(S33)
    Qc%name = 'Qcrit'
    O1%name = 'Ox'
    O2%name = 'Oy'
    O3%name = 'Oz'

    output_format = "paraview"
    if(parser_is_defined("output format")) call parser_read("output format", output_format)
    call parser_read('slice or 3D', type_para)

    if ((trim(output_format)=="paraview").or.(trim(type_para)=='3Dplus')) then
        call vtkxml_scalar(tag_Qc, Qc%values, Qc%name)
        call vtkxml_scalar(tag_O(1), O1%values, O1%name)
        call vtkxml_scalar(tag_O(2), O2%values, O2%name)
        call vtkxml_scalar(tag_O(3), O3%values, O3%name)
    else
        if(     .not. write_datalayout_scalar(O1,iter,spec_rank)    .or.    &
            &   .not. write_datalayout_scalar(O1,iter,spec_rank)    .or.    &
            &   .not. write_datalayout_scalar(O1,iter,spec_rank)) then
            write (6,'(a,i0)') 'FAILED to write output for vorticity'
            return
        end if
        O1%name = 'vorticity_name'
        O1%values = (O1%values**2) + (O2%values**2) + (O3%values**2)
        if(     .not. write_datalayout_scalar(O1,iter,spec_rank)    .or.    &
            &   .not. write_datalayout_scalar(Qc,iter,spec_rank)) then
            write (6,'(a,i0)') 'FAILED to write output for Q criterium and vorticityi norm'
            return
        end if
    end if
    call deleteDataLayout(Qc)
    call deleteDataLayout(O1)
    call deleteDataLayout(O2)
    call deleteDataLayout(O3)

 end function WriteQcOmega



!> Compute the quantity of differential diffusion and write a 2D slide and and
!! of its values
!! @autor  Jean-Baptiste Lagaert
!!    @param[in]    scalarArray     = array of scalars
!!    @param[in]    nbscalar        = number of scalars in the array
!!    @param[in]    scalar_partArray= array of scalars solved with
!!                                      particle/spectral methods
!!    @param[in]    nbscalar_part   = number of scalars in the array - scalars
!!                                      solved with particle/spectral methods.
!!    @param[in]    ite             = current time iteration
!!    @param[in]    n_iter          = number max of iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
 function differential_diffusion(scal_1, scal_2, nbscal, nbscal_part, ScalArray,   &
    &  Scal_partArray, diff_diff, mix_th, tag_diff_diff, spec_rank, &
    &  ite, n_iter, time) result(success)


    ! Input/Output
    integer, intent(in)                             :: scal_1, scal_2, nbscal, nbscal_part
    type(REAL_DATA_LAYOUT),dimension(:),intent(in)  :: scalArray, scal_partArray
    type(REAL_DATA_LAYOUT), intent(inout)           :: diff_diff
    real(WP)                                        :: mix_th
    integer                                         :: tag_diff_diff
    integer, intent(in)                             :: ite
    integer, intent(in)                             :: n_iter
    integer, intent(in)                             :: spec_rank
    logical, intent(in), optional                   :: time
    logical                                         :: success

    ! Local variables
    type(COMPLEX_DATA_LAYOUT)   :: tmp_cmplx
    integer                     :: sca, i, j, k
    real(WP), dimension(:), allocatable :: avg2D
    real(WP)                            :: dz       ! space step along Z
    logical                             :: time_or_ite  ! output depend of time avancement or iterations ?
    integer                             :: file_id  ! to identify output files
    character(len=50)                   :: file_name! name of output file

    ! ===== Initialisation =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! Be carful: for these setup, the second scalar is not equal to the computed
    ! scalar2 but is equal to (1-scalar2) (wich has the right init (1-scalar1)).

    ! Compute derivation from spectral space
    if(scal_1<=nbscal) then
        ! The scalar is a full-spectral one
        if(.not. copyStructOnly(ScalArray(scal_1), diff_diff)) return
        if(.not. copyStructOnly(ScalArrayK(scal_1), tmp_cmplx)) return
        do k = tmp_cmplx%zmin, tmp_cmplx%zmax
            do j = tmp_cmplx%ymin, tmp_cmplx%ymax
                do i = tmp_cmplx%xmin, tmp_cmplx%xmax
                    tmp_cmplx%values(i,j,k) = -(schmidt(scal_1)) *      &
                        & ((ScalWN(scal_1)%kx(i)*ScalWN(scal_1)%kx(i))+ &
                        & (ScalWN(scal_1)%ky(j)*ScalWN(scal_1)%ky(j)) + &
                        & (ScalWN(scal_1)%kz(k)*ScalWN(scal_1)%kz(k)))  &
                        & *ScalArrayK(scal_1)%values(i,j,k)
                end do
            end do
        end do
    else
        ! The scalar is a scalar solved with mixed particle/spectral method
        sca = scal_1 - nbscal
        if(.not. copyStructOnly(Scal_partArray(sca), diff_diff)) return
        if(.not. copyStructOnly(Scal_partArrayK(sca), tmp_cmplx)) return
        do k = tmp_cmplx%zmin, tmp_cmplx%zmax
            do j = tmp_cmplx%ymin, tmp_cmplx%ymax
                do i = tmp_cmplx%xmin, tmp_cmplx%xmax
                    tmp_cmplx%values(i,j,k) = -(schmidt_part(sca))* &
                        & ((Scal_partWN(sca)%kx(i)*Scal_partWN(sca)%kx(i))+ &
                        & (Scal_partWN(sca)%ky(j)*Scal_partWN(sca)%ky(j)) + &
                        & (Scal_partWN(sca)%kz(k)*Scal_partWN(sca)%kz(k)))  &
                        & *Scal_partArrayK(sca)%values(i,j,k)
                end do
            end do
        end do
    end if

    ! Go back to real space and on grid associate to the second scalar
    if(scal_2<=nbscal) then
        ! The second scalar is a full-spectral one
        if (samelayout(ScalArray(scal_2), diff_diff)) then
            ! Same resolution is used !!
            call btran(tmp_cmplx,diff_diff,success)
            where ((ScalArray(scal_2)%values<mix_th).or.(ScalArray(scal_2)%values>1-mix_th))
                diff_diff%values = 0
            end where
        else
            if(.not. copyStructOnly(ScalArray(scal_2), diff_diff)) return
            if (.not. transdiscretization (tmp_cmplx,ScalArray(scal_2),diff_diff)) return
            where ((ScalArray(scal_2)%values<mix_th).or.(ScalArray(scal_2)%values>1-mix_th))
                diff_diff%values = 0
            end where
        end if
    else
        ! Second scalar is solved with particle/spectral method
        sca = scal_2 - nbscal
        if (samelayout(Scal_partArray(sca), diff_diff)) then
            ! Same resolution is used !!
            call btran(tmp_cmplx,diff_diff,success)
            where ((Scal_partArray(sca)%values<mix_th).or.(Scal_partArray(sca)%values>1.-mix_th))
                diff_diff%values = 0
            end where
        else
            if(.not. copyStructOnly(Scal_partArray(sca), diff_diff)) return
            if (.not. transdiscretization (tmp_cmplx,Scal_partArray(sca),diff_diff)) return
            where ((Scal_partArray(sca)%values<mix_th).or.(Scal_partArray(sca)%values>1.-mix_th))
                diff_diff%values = 0
            end where
        end if
    end if
    call deleteDataLayout(tmp_cmplx)

    ! Computation done; now perform output: slice 2D and avg 2D
    ! Slice
    diff_diff%name = 'differential_diff'
    call vtkxml_write_slice(tag_diff_diff, diff_diff%values, diff_diff%name )
    ! avg 2D
    allocate(avg2D(diff_diff%Nz))
    if(.not. computeFieldAvg2D(diff_diff,avg2D,spec_rank)) then
        success = .false.
        return
    end if
    ! Write Avg 2D
    if (spec_rank==0) then
        dz = diff_diff%Lz/diff_diff%Nz
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'differentialDiff_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'differentialDiff_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'differentialDiff_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'differentialDiff_T_',ite,'deltaT.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'differentialDiff_T_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'differentialDiff_T_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=300)
        do k = 1, size(avg2D)
            write(file_id, '(2(e14.6,x))') k*dz, avg2D(k)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if


    call deleteDataLayout(diff_diff)
    deallocate(avg2D)

    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file - "'//trim(file_name)//'"'

 end function differential_diffusion


!-------------------------------------------------------------------------------------------
!> transdiscretize field at avs format in order to be able to restart a stopped simulation.
!!    @param[in]    spec_rank       = rank inside the spectral communicator
!!    @param[in]    scalArray       = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    scalArray_part  = array containing all the scalar fields (scalars solved 
!!                                  with pseudo-spectral methods)
!!    @param[in]    U               = velocity along X
!!    @param[in]    V               = velocity along Y
!!    @param[in]    W               = velocity along Z
!!    @param[in]    B               = magentic field (3D vector)
!!    @param[in]    imhd            = 1 if MHD is ON, 0 if it is OFF
!!    @param[in]    iter            = number of the current iteration
!!    @param[in]    iterMax         = total number of iterations
!!    @return       success         = .true. if everything is right !
!! @author  Antoine Vollant
!-------------------------------------------------------------------------------------------
 function transDiscAvsWrite(spec_rank,scalArray, scal_partArray, U, V, &
                      & W,B,imhd,iter,iterMax) result(success)

   use rediscretization_tools , only : transdiscretizeField

   !I/O data
   integer,intent(in)                                          :: spec_rank
   type(real_data_layout),intent(in)                           :: U, V ,W
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: B
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scalArray
   type(REAL_DATA_LAYOUT), intent(in),pointer, dimension(:)    :: scal_partArray
   logical,intent(in)                                          :: imhd
   logical                                                     :: success 
   integer,intent(in)                                          :: iter
   integer,intent(in)                                          :: iterMax


   !Local data
   integer                                 :: nx,ny,nz,j,nbcpus,ires 
   type(real_data_layout)                  :: U_trans, V_trans ,W_trans
   type(REAL_DATA_LAYOUT), dimension(:),pointer    :: B_trans
   type(REAL_DATA_LAYOUT), dimension(:),pointer    :: scalArray_trans
   type(REAL_DATA_LAYOUT), dimension(:),pointer    :: scal_partArray_trans
   
    success = .false.
    nbcpus=getnbcpus()
    if (.not. GetNxNyNz('Number of points for rediscretization',nx,ny,nz)) then
      WRITE (6,'(a)') 'Failed in GetNxNyNz' 
      RETURN
    endif
    ! transdiscretize velocity
    if (.not. transdiscretizeField(U,U_trans,spec_rank,nbcpus,nx,ny,nz)) then
      WRITE (6,'(a)') 'FAILED in transdiscretize U' 
      RETURN
    endif
    if (.not. transdiscretizeField(V,V_trans,spec_rank,nbcpus,nx,ny,nz)) then
      WRITE (6,'(a)') 'FAILED in transdiscretize V' 
      RETURN
    endif
    if (.not. transdiscretizeField(W,W_trans,spec_rank,nbcpus,nx,ny,nz)) then
      WRITE (6,'(a)') 'FAILED in transdiscretize W' 
      RETURN
    endif
    ! transdiscretize scalar - pseudo-spectral solver
    IF(ASSOCIATED(ScalArray)) THEN
        allocate(ScalArray_trans(size(ScalArray)),stat=ires)
        if ( ires .ne. 0 ) then
           WRITE (6,'(a)') 'FAILED in allocation ScalArray_trans'
           return
        endif
        DO j=1,SIZE(ScalArray)
          if (.not. transdiscretizeField(ScalArray(j),ScalArray_trans(j),spec_rank,nbcpus,nx,ny,nz)) then
            WRITE (6,'(a,i0)') 'FAILED in transdiscretize scalar field',j
            RETURN
          ENDIF
        ENDDO
    ENDIF  
    ! transdiscretize scalar - spectral/particles solver
    IF(ASSOCIATED(Scal_partArray)) THEN 
        allocate(Scal_partArray_trans(size(Scal_partArray)),stat=ires)
        if ( ires .ne. 0 ) then
           WRITE (6,'(a)') 'FAILED in allocation Scal_partArray_trans'
           return
        endif
        DO j=1,SIZE(Scal_partArray)
          if (.not. transdiscretizeField(Scal_partArray(j),Scal_partArray_trans(j),spec_rank,nbcpus,nx,ny,nz)) then
             WRITE (6,'(a,i0)') 'FAILED in transdiscretize scalar field',j
             RETURN
          ENDIF
        ENDDO
    ENDIF  
    ! transdiscretize mhd
    IF (imhd) THEN
      allocate(B_trans(3),stat=ires)
      if ( ires .ne. 0 ) then
         WRITE (6,'(a)') 'FAILED in allocation Scal_partArray_trans'
         return
      endif
      if (.not. transdiscretizeField(B(1),B_trans(1),spec_rank,nbcpus,nx,ny,nz) .or. &
         &.not. transdiscretizeField(B(2),B_trans(2),spec_rank,nbcpus,nx,ny,nz) .or. &
         &.not. transdiscretizeField(B(3),B_trans(3),spec_rank,nbcpus,nx,ny,nz)) then
        WRITE (6,'(a)') '[ERROR] FAILED in transdiscretize B'
        RETURN
      endif
    ENDIF

    if (.not. AvsWrite(spec_rank,ScalArray_trans, Scal_partArray_trans, U_trans, V_trans, &
                    & W_trans,B_trans,imhd,iter,iterMax,.true.)) then
      WRITE (6,'(a)') '[ERROR] in avs write'
      RETURN
    ENDIF

    if ( associated(B_trans) ) deallocate(B_trans)
    if ( associated(Scal_partArray_trans) ) deallocate(Scal_partArray_trans)
    if ( associated(ScalArray_trans) ) deallocate(ScalArray_trans)

    success = .true.

 end function transDiscAvsWrite



! ===== Spectrum computations =====

!------------------------------------------------------------------------------
!> Compute and save spectrum of all field (scalars and velocity).
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite             = current time iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    n_iter          = number of max iteration
!!    @param[in]    imhd            = logical equals to true if mhd is present
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of all scalars fields and of 
!!    kinetic energy. It takes directly the field in Fourrier Space from 
!!    the data. This allows to avoid useless and expansive FFT transform,
!!    but it could only be used during a simulation (and not be performed 
!!    after).
!------------------------------------------------------------------------------
function compute_all_spectrum(ite, spec_rank, n_iter,imhd, time) result(success)

    use datalayout

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
    logical, intent(in)             :: imhd
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    logical                         :: time_or_ite      ! output could be done at each N 
                                                        ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! Compute spectrum
    if (.not.compute_all_scalars_spectrum(ite, spec_rank,n_iter,time_or_ite)) then
        write(*,'(a)') '[WARNING] Failed to compute spectrum of scalars fiels'
        return
    else if (.not.compute_velo_spectrum(ite, spec_rank,n_iter,time_or_ite)) then
        write(*,'(a)') '[WARNING] Failed to compute spectrum of kinetic energy'
        return
    else if (imhd) then
        if (.not.compute_mhd_spectrum(ite, spec_rank,n_iter, time_or_ite)) then
           write(*,'(a)') '[WARNING] Failed to compute spectrum of kinetic energy'
           return
        end if
    end if

    ! Return sucess
    success = .true.

end function compute_all_spectrum


!------------------------------------------------------------------------------
!> Compute and save spectrum of all scalars field.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite             = current time iteration
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    time            = optional logical do adjust name if
!!                                      spectrum computation are done depending 
!!                                      time and not iterations.
!!    @return       success         = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of all scalars fields. It 
!!    takes directly the field in Fourrier Space from the data module. This 
!!    allows to avoid useless and expansive FFT transform, but it could 
!!    only be used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
function compute_all_scalars_spectrum(ite, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, sca, ik ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    character(len=15)                       :: name_short       ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    ! -- Scalar solved with pseudo-spectral solver --
    if (allocated(ScalArrayk)) then
        do sca=1, size(ScalArrayk)
            kc = 2.0_WP*acos(-1.0_WP)*(ScalArrayk(sca)%nx-1)/Uk%Lx
            size_spec = size(ScalWN(sca)%kx) ! Only this because de-aliasing
            allocate(spectrum(size_spec))
            allocate(WN_norm(size_spec))
            WN_norm = -1
            spectrum = 0
            ! Compute
            do k = ScalArrayK(sca)%zmin,ScalArrayK(sca)%zmax
                do j =  ScalArrayK(sca)%ymin,ScalArrayK(sca)%ymax
                    do i =  ScalArrayK(sca)%xmin,ScalArrayK(sca)%xmax
                        ! Compute norm of wave number
                        kk = ScalWN(sca)%kx(i)*ScalWN(sca)%kx(i) + &
                            & ScalWN(sca)%ky(j)*ScalWN(sca)%ky(j) + &
                            & ScalWN(sca)%kz(k)*ScalWN(sca)%kz(k)
                        kk = sqrt(real(kk))
                        if ((kk>eps).and.(kk<=kc)) then
                            ! Compute indice
                            ik = nint(kk/dk) + 1
                            ! And then update spectrum
                            !WN_norm(ik) = dk*real(ik-1)
                            spectrum(ik) = spectrum(ik) &
                                & + ScalArrayk(sca)%values(i,j,k)*conjg(ScalArrayk(sca)%values(i,j,k))
                        end if
                    end do
                end do
            end do
            ! Sum values on all processes
            spectrum = DoGlobalVectorSum(spec_rank, spectrum)
            ! Write output
            if (spec_rank==0) then
                do ik = 1, size_spec
                    WN_norm(ik) = dk*real(ik-1)
                end do
                file_id = iopen()
                name_short = trim(adjustl(ScalArrayk(sca)%name))
                if (.not. time_or_ite) then
                    if (n_iter .lt. 1000) then
                      write(file_name,'(a,a,a,i3.3,a)') 'spec_',name_short,'_at_',ite,'.table'
                    elseif (n_iter .lt. 1000000) then
                      write(file_name,'(a,a,a,i6.6,a)') 'spec_',name_short,'_at_',ite,'.table'
                    else
                      write(file_name,'(a,a,a,i0,a)') 'spec_',name_short,'_at_',ite,'.table'
                    endif
                else
                    if (max(n_iter, ite) .lt. 1000) then
                      write(file_name,'(a,a,a,i3.3,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    elseif (n_iter .lt. 1000000) then
                      write(file_name,'(a,a,a,i6.6,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    else
                      write(file_name,'(a,a,a,i0,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    endif
                end if
                file_name = trim(adjustl(file_name))
                open(unit=file_id, file=file_name, form="FORMATTED", ERR=100)
                do ik = 1, size_spec
                    if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
                end do
                close(file_id)
                file_id = iclose(file_id)
            end if
            ! Free memory for next scalar field
            deallocate(spectrum)
            deallocate(WN_norm)
        end do
    end if

    ! -- Scalar solved with mixed partcile/spectral solver --
    if (allocated(Scal_partArrayk)) then
        do sca=1, size(Scal_partArrayk)
            kc = 2.0_WP*acos(-1.0_WP)*(Scal_partArrayk(sca)%nx-1)/Uk%Lx
            size_spec = size(Scal_partWN(sca)%kx)
            allocate(spectrum(size_spec))
            allocate(WN_norm(size_spec))
            spectrum = 0
            WN_norm = -1
            ! Compute
            do k = Scal_partArrayK(sca)%zmin,Scal_partArrayK(sca)%zmax
                do j =  Scal_partArrayK(sca)%ymin,Scal_partArrayK(sca)%ymax
                    do i =  Scal_partArrayK(sca)%xmin,Scal_partArrayK(sca)%xmax
                        ! Compute norm of wave number
                        kk = Scal_partWN(sca)%kx(i)*Scal_partWN(sca)%kx(i) + &
                            & Scal_partWN(sca)%ky(j)*Scal_partWN(sca)%ky(j) + &
                            & Scal_partWN(sca)%kz(k)*Scal_partWN(sca)%kz(k)
                        kk = sqrt(kk)
                        if ((kk>eps).and.(kk<=kc)) then
                            ! Compute indice
                            ik = nint(kk/dk) + 1
                            ! And then update spectrum
                            !WN_norm(ik) = kk
                            spectrum(ik) = spectrum(ik) &
                                & + Scal_partArrayk(sca)%values(i,j,k)*conjg(Scal_partArrayk(sca)%values(i,j,k))
                        end if
                    end do
                end do
            end do
            ! Sum values on all processes
            spectrum = DoGlobalVectorSum(spec_rank, spectrum)
            ! Write output
            if (spec_rank==0) then
                do ik = 1, size_spec
                    WN_norm(ik) = dk*real(ik-1)
                end do
                file_id = iopen()
                name_short = trim(adjustl(Scal_partArrayk(sca)%name))
                if (.not. time_or_ite) then
                    if (n_iter .lt. 1000) then
                      write(file_name,'(a,a,a,i3.3,a)') 'spec_',name_short,'_at_',ite,'.table'
                    elseif (n_iter .lt. 1000000) then
                      write(file_name,'(a,a,a,i6.6,a)') 'spec_',name_short,'_at_',ite,'.table'
                    else
                      write(file_name,'(a,a,i0,a)') 'spec_',name_short,ite,'.table'
                    endif
                else
                    if (max(n_iter, ite) .lt. 1000) then
                      write(file_name,'(a,a,a,i3.3,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    elseif (n_iter .lt. 1000000) then
                      write(file_name,'(a,a,a,i6.6,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    else
                      write(file_name,'(a,a,a,i0,a)') 'spec_',name_short,'_at_',ite,'deltaT.table'
                    endif
                end if
                file_name = trim(adjustl(file_name))
                open(unit=file_id, file=file_name, form="FORMATTED", ERR=150)
                do ik = 1, size_spec
                    if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
                end do
                close(file_id)
                file_id = iclose(file_id)
            end if
            ! Free memory for next scalar field
            deallocate(spectrum)
            deallocate(WN_norm)
        end do
    end if

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
100 write(6,'(a,i0,a, i0)') 'Unable to open file "'//trim(file_name)//'"for scal(',sca,') - rank ', spec_rank
    return
150 write(6,'(a,i0,a,i0)') 'Unable to open file "'//trim(file_name)//'"for sca_part(',sca,') - rank ', spec_rank
    return


end function compute_all_scalars_spectrum


!------------------------------------------------------------------------------
!> Compute and save spectrum of kinetic energy.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite     = current iteration step
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of kinetic energy. It takes
!!    directly the field in Fourrier Space from the solver. This allows
!!    to avoid useless and expansive FFT transform, but it could only be
!!    used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
function compute_velo_spectrum(ite, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    !size_spec = size(VelWN%kx)*size(VelWN%ky)*size(VelWN%kz)
    size_spec = size(VelWN%kx) ! Only this because de-aliasing
    allocate(spectrum(size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! Compute
    do k = UK%zmin,UK%zmax
        do j =  UK%ymin,UK%ymax
            do i =  UK%xmin,UK%xmax
                ! Compute norm of wave number
                kk = VelWN%kx(i)*VelWN%kx(i) + &
                    & VelWN%ky(j)*VelWN%ky(j) + &
                    & VelWN%kz(k)*VelWN%kz(k)
                kk = sqrt(real(kk))
                if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                    ik = nint(kk/dk) + 1
                    ! And then update spectrum
                    spectrum(ik) = spectrum(ik) + Uk%values(i,j,k)*conjg(Uk%values(i,j,k)) &
                            & + Vk%values(i,j,k)*conjg(Vk%values(i,j,k)) &
                            & + Wk%values(i,j,k)*conjg(Wk%values(i,j,k))
                    !WN_norm(ik) = dk*real(ik-1)
                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
          WN_norm(ik) = dk*dble(ik-1)
        end do
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'spec_vel_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'spec_vel_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'spec_veli_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'spec_vel_at_',ite,'deltaT.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'spec_vel_at_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'spec_vel_at_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
        write(file_id,'(a)') '# Kinetic spectrum'
        write(file_id,'(a)') '# wave_number and spectrum '
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
200 write(6,'(a)') 'Unable to open file - velo "'//trim(file_name)//'"'
    return

end function compute_velo_spectrum


!------------------------------------------------------------------------------
!> Compute and save spectrum of magnetic energy.
!! @author Jean-Baptiste Lagaert and Guillaume Balarac
!!    @param[in]    ite     = current iteration step
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending 
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of magnetic energy. It takes 
!!    directly the field in Fourrier Space from the data module. This allows 
!!    to avoid useless and expansive FFT transform, but it could only be 
!!    used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
function compute_mhd_spectrum(ite, spec_rank,n_iter, time) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
    logical, intent(in), optional   :: time
    logical                         :: success
    ! Others variables
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
    if(present(time)) then
        time_or_ite = time
    else
        time_or_ite = .false.
    end if

    ! ===== Check if spectral space step is the same above all direction =====
    if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    !size_spec = size(VelWN%kx)*size(VelWN%ky)*size(VelWN%kz)
    size_spec = size(BfieldWN%kx) ! Only this because de-aliasing
    allocate(spectrum(size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! Compute
    do k = Bk(1)%zmin,Bk(1)%zmax
        do j =  Bk(1)%ymin,Bk(1)%ymax
            do i =  Bk(1)%xmin,Bk(1)%xmax
                ! Compute norm of wave number
                kk = BfieldWN%kx(i)*BfieldWN%kx(i) + &
                    & BfieldWN%ky(j)*BfieldWN%ky(j) + &
                    & BfieldWN%kz(k)*BfieldWN%kz(k)
                kk = sqrt(real(kk))
                if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                    ik = nint(kk/dk) + 1
                    ! And then update spectrum
                    spectrum(ik) = spectrum(ik) + Bk(1)%values(i,j,k)*conjg(Bk(1)%values(i,j,k)) &
                            & + Bk(2)%values(i,j,k)*conjg(Bk(2)%values(i,j,k)) &
                            & + Bk(3)%values(i,j,k)*conjg(Bk(3)%values(i,j,k))
                    !WN_norm(ik) = dk*real(ik-1)
                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        file_id = iopen()
        if (.not. time_or_ite) then
            if (n_iter .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'spec_mhd_at_',ite,'.table'
            elseif (n_iter .lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'spec_mhd_at_',ite,'.table'
            else
               write(file_name,'(a,i0,a)') 'spec_mhd_at_',ite,'.table'
            endif
        else
            if (max(n_iter, ite) .lt. 1000) then
               write(file_name,'(a,i3.3,a)') 'spec_mhd_at_',ite,'deltaT.table'
            elseif (max(n_iter, ite).lt. 1000000) then
               write(file_name,'(a,i6.6,a)') 'spec_mhd_at_',ite,'deltaT.table'
            else
               write(file_name,'(a,i0,a)') 'spec_mhd_at_',ite,'deltaT.table'
            endif
        end if
        file_name = trim(adjustl(file_name))
        open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
        write(file_id,'(a)') '# Magnetic spectrum'
        write(file_id,'(a)') '# wave_number and spectrum '
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(file_id, '(e14.6,x,e14.6)') WN_norm(ik), spectrum(ik)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
200 write(6,'(a)') 'Unable to open file - mhd "'//trim(file_name)//'"'
    return

end function compute_mhd_spectrum

!------------------------------------------------------------------------------
!> Read different spectrum file and compute the (time) average. In order to
!! smooth the noise.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    tag         = standart name of the file to read
!!    @param[in]    n_start     = first iteration
!!    @param[in]    n_end       = last iteration
!!    @param[in]    spec_rank   = mpi rank inside spectral communicator
!!    @return       success     = logical equals to true is everything is right.
!! @details
!!        This function compute the average of the spectrum for different file.
!!    Be careful of what you do if the original spectrum do not come from fields 
!!    with the same resolution.
!------------------------------------------------------------------------------
function spectrum_time_avg(spec_rank) result(success)

    use fileio
    use precision_tools
    use parser_tools

    ! Input/Output
    integer, intent(in)                     :: spec_rank
    logical                                 :: success
    ! Others variables
    integer, parameter :: line_length=4960
    character(len=50)                       :: tag              ! generic file name
    character(len=10)                       :: extention        ! file type
    integer                                 :: n_start, n_end   ! first and last file to read
    integer                                 :: n_inc            ! increment between to files
    integer                                 :: length           ! number of character to read the iteration indice
    character(len=line_length)              :: buffer           ! to count line number
    real(WP), dimension(:,:), allocatable   :: spectrum         ! to save spectrum
    real(WP), dimension(2)                  :: tmp_spectrum     ! to save current read
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    character(len=10)                       :: file_format      ! to write the file name at the write format
    integer                                 :: file_header      ! number of line of the header
    integer                                 :: file_ind         ! indice of current spectrum file
    integer                                 :: file_number      ! number of file readed
    integer                                 :: ierr       ! error code
    integer                                 :: iline, nlines    ! current line and total number of lines
    !integer                                 :: iunit,ierr,limiter,i,j

    ! Init
    success = .false.

    ! ===== Initialisation =====
    call parser_read('Spectrum files generic name', tag)
    extention = '.table'
    call parser_read('Spectrum files indice length', length)
    call parser_read('Spectrum files start at', n_start)
    call parser_read('Spectrum files end at', n_end)
    call parser_read('Spectrum files increment', n_inc)
    call parser_read('Spectrum files header', file_header)
    file_ind = n_start
    file_number = 1

    ! ===== Read the first file to know the lenght =====
    write(file_format,'(a,i0,a,i0,a)') '(a,i',length,'.',length,',a)'
    write(file_name, file_format) trim(tag), file_ind, trim(extention)
    file_id = iopen()
    open (file_id,file=file_name,form='formatted',status='old',err=300)

    ! Count the number of lines in the file
    ierr = 0
    nlines = 0
    do while (ierr .eq. 0)
       read(file_id,'(a)',iostat=ierr) buffer
       nlines = nlines + 1
    end do
    rewind(file_id)

    ! Allocate to the right size
    nlines = nlines -1
    allocate(spectrum(2,file_header+1:nlines))

    ! ===== Init wave number (and read the first spectrum) =====
    ! -- To check if fields seem to come from the same resolution --
    ierr = 0
    iline = 0
    ! Skip the header
    do while ((ierr .eq. 0).and.(iline<file_header))
        iline=iline+1
        read(file_id,'(a)',iostat=ierr) buffer
    end do
    ! And read data
    do while (ierr .eq. 0)
        iline=iline+1
        read(file_id, '(e14.6,x,e14.6)',iostat=ierr) tmp_spectrum(1),tmp_spectrum(2)
        if (ierr==0) then
            spectrum(:,iline) = tmp_spectrum
        else
            exit 
        end if
    enddo 

    ! Close file
    close(file_id)
    file_id = iclose(file_id)

    ! ===== Read other files and compute the average =====
    do while (file_ind<n_end)
        file_ind = file_ind + n_inc
        file_number = file_number + 1
        ! -- Read and update values --
        write(file_name, file_format) trim(tag), file_ind, trim(extention)
        open (file_id,file=file_name,form='formatted',status='old',err=300)
        ierr = 0
        iline = 0
        ! Skip the header
        do while ((ierr .eq. 0).and.(iline<file_header))
            iline=iline+1
            read(file_id,'(a)',iostat=ierr) buffer
        end do
        ! And read data
        do while (ierr .eq. 0)
            read(file_id, '(e14.6,x,e14.6)',iostat=ierr) tmp_spectrum(1),tmp_spectrum(2)
            if (ierr.ne.0) exit 
            iline = iline + 1
            if(iline<nlines) then
                if (abs(tmp_spectrum(1)-spectrum(1,iline))<0.00001) then
                    spectrum(2,iline) = ((file_number-1)*spectrum(2,iline) &
                        & + tmp_spectrum(2))/dble(file_number)
                else
                    write(*,'(a,a,a,i0,a)') '[ERROR] Spectrum avg : file ', trim(file_name),&
                        & ' seems not be computed with the right resolution: at line '      &
                        & , iline, 'WN norm is not right'
                    write(*,'(a,a,a,i0,a)') '[ERROR] Spectrum avg - WN read is: ', tmp_spectrum(1), &
                        & 'and must be ', spectrum(1,iline)
                    stop
                end if
            end if
        enddo 
        ! -- Check size --
        if (iline/= size(spectrum,2)+file_header) then
            write(*,'(a,a,a,i0,a,i0)') '[ERROR] Spectrum avg : file ',  &
                & trim(file_name), ' size is ', (iline), ' and must be ', size(spectrum, 2)+file_header
            stop
        end if
        ! -- Close file --
        close(file_id)
        file_id = iclose(file_id)
    end do


    ! ===== Write output =====
    if (spec_rank==0) then
        file_id = iopen()
        write(file_name, '(a,a,a)') trim(tag), '_avg_', trim(extention)
        open(unit=file_id, file=trim(adjustl(file_name)), form="FORMATTED", ERR=300)
        do iline = file_header+1, nlines
            write(file_id, '(e14.6,x,e14.6)') spectrum(1,iline), spectrum(2,iline)
        end do
        close(file_id)
        file_id = iclose(file_id)
    end if

    ! ===== Free memory =====
    deallocate(spectrum)

    ! ===== Return sucess =====
    success = .true.
    return

    ! ===== IO error =====
300 write(6,'(a)') 'Unable to open file "'//trim(file_name)//'" for spectrum'
    return

end function spectrum_time_avg

!------------------------------------------------------------------------------
!> Compute and save spectrum of magnetic energy.
!! @author Jean-Baptiste Lagaert and Guillaume Balarac
!!    @param[in]    ite     = current iteration step
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending 
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of magnetic energy. It takes 
!!    directly the field in Fourrier Space from the data module. This allows 
!!    to avoid useless and expansive FFT transform, but it could only be 
!!    used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
subroutine compute_Helicity_spectrum(VecXk,VecYk,VecZk,&
                                  &RVecXk,RVecYk,RVecZk,Wave, spec_rank) 

    use datalayout
!     use data
    use fileio

    ! Input/Output
!     integer, intent(in)             :: ite
!     integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
!     logical, intent(in), optional   :: time
!     logical                         :: success
    TYPE(COMPLEX_DATA_LAYOUT) :: VecXk,VecYk,VecZk,&
                                  &RVecXk,RVecYk,RVecZk
   TYPE(WaveNumbers), INTENT(IN) :: Wave

    ! Others variables
    real(WP), dimension(:), allocatable     :: spectrum         ! to save spectrum
    real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
    logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.


    ! ===== Check if spectral space step is the same above all direction =====
    if((VecXk%Lx/=VecXk%Ly).or.(VecXk%Lx/=VecXk%Lz)) then
        write(*,'(a)') '[warning] Spectrum not yet implemented if not(Lx=Ly=Lz)'
        return
    else
        dk=2.0_WP*acos(-1.0_WP)/VecXk%Lx
    end if
    kc = 2.0_WP*acos(-1.0_WP)*(VecXk%nx-1)/VecXk%Lx
    eps = dk/2.0

    ! ===== Compute and write spectrum =====
    !size_spec = size(VelWN%kx)*size(VelWN%ky)*size(VelWN%kz)
    size_spec = size(wave%kx) ! Only this because de-aliasing
    allocate(spectrum(size_spec))
    allocate(WN_norm(size_spec))
    WN_norm = -1
    spectrum = 0
    ! Compute
    do k = VecXk%zmin,VecXk%zmax
        do j =  VecXk%ymin,VecXk%ymax
            do i =  VecXk%xmin,VecXk%xmax
                ! Compute norm of wave number
                kk = wave%kx(i)*wave%kx(i) + &
                    & wave%ky(j)*wave%ky(j) + &
                    & wave%kz(k)*wave%kz(k)
                kk = sqrt(real(kk))
                if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
                    ik = nint(kk/dk) + 1
                    ! And then update spectrum
                    spectrum(ik) = spectrum(ik) + 0.5_WP*(VecXk%values(i,j,k)*conjg(RVecXk%values(i,j,k)) &
                            & + VecYk%values(i,j,k)*conjg(RVecYk%values(i,j,k)) &
                            & + VecZk%values(i,j,k)*conjg(RVecZk%values(i,j,k)) &
                            & + RVecXk%values(i,j,k)*conjg(VecXk%values(i,j,k)) &
                            & + RVecYk%values(i,j,k)*conjg(VecYk%values(i,j,k)) &
                            & + RVecZk%values(i,j,k)*conjg(VecZk%values(i,j,k)))
                    !WN_norm(ik) = dk*real(ik-1)
                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalVectorSum(spec_rank, spectrum)

    ! ===== Write output =====
    if (spec_rank==0) then
        do ik = 1, size_spec
            WN_norm(ik) = dk*real(ik-1)
        end do
        file_id = iopen()

       open(10,file='spec_Kin_Hel.table',form='formatted')
!        write(10,*) meanadvecQDM, meandiffvisQDM, meandissvisQDM,&
!                  & somme,meanEk
!        close(10)
 
! 
!                write(file_name,'(a)') 'spec_Kin_Hel.table'
!         file_name = trim(adjustl(file_name))
!         open(unit=file_id, file=file_name, form="FORMATTED", ERR=200)
!         write(file_id,'(a)') '# Kinetic Helicity spectrum'
!         write(file_id,'(a)') '# wave_number and spectrum '
        do ik = 1, size_spec
            if (WN_norm(ik)>=0) write(10, *) WN_norm(ik), spectrum(ik)
        end do
       close(10)
 
!         close(file_id)
!         file_id = iclose(file_id)

    end if
    ! Free memory for next scalar field
    deallocate(spectrum)
    deallocate(WN_norm)

    ! ===== Return sucess =====
!     success = .true.
!     return

    ! ===== IO error =====
! 200 write(6,'(a)') 'Unable to open file - Hkin "'//trim(file_name)//'"'
    return

end subroutine compute_Helicity_spectrum

!------------------------------------------------------------------------------
!> Compute and save spectrum of kinetic energy.
!! @author Jean-Baptiste Lagaert
!!    @param[in]    ite     = current iteration step
!!    @param[in]    n_iter  = max iteration step
!!    @param[in]    time    = optional logical do adjust name if
!!                              spectrum computation are done depending
!!                              time and not iterations.
!!    @return       success = logical equals to true is everything is right.
!! @details
!!        This function compute the spectrum of kinetic energy. It takes
!!    directly the field in Fourrier Space from the solver. This allows
!!    to avoid useless and expansive FFT transform, but it could only be
!!    used during a simulation (and not be performed after).
!------------------------------------------------------------------------------
function compute_Integral_scale(ite, spec_rank,n_iter, time,mu) result(success)

    use datalayout
    use data
    use fileio

    ! Input/Output
    integer, intent(in)             :: ite
    integer, intent(in)             :: n_iter
    integer, intent(in)             :: spec_rank
    real(wp), intent(in)   :: time
    logical                         :: success
    real(wp) :: ReLint,mu,Lint,Urms
    ! Others variables
    real(WP)     :: spectrum         ! to save integral of spectrum
    !real(WP), dimension(:), allocatable     :: WN_norm          ! to know the different wave number modules
    real(WP)     :: Lintspectrum         ! to save Lint 
    integer                                 :: size_spec        ! spectrum size
    integer                                 :: i, j ,k, ik      ! indices
    real(WP)                                :: kk               ! norm of wave number
    real(WP)                                :: dk               ! space step in Fourrier space
    real(WP)                                :: kc, eps          ! to filter spectrum
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file
!     logical                                 :: time_or_ite      ! output could be done at each N 
                                                                ! iterations or at each N*deltaT time.

    ! ===== Init =====
    success = .false.
!     if(present(time)) then
!         time_or_ite = time
!     else
!         time_or_ite = .false.
!     end if

    ! ===== Check if spectral space step is the same above all direction =====
!     if((Uk%Lx/=Uk%Ly).or.(Uk%Lx/=Uk%Lz)) then
!         write(*,'(a)') '[warning] Integral Scale not yet implemented if not(Lx=Ly=Lz)'
!         return
!     else
        dk=2.0_WP*acos(-1.0_WP)/Uk%Lx
!     end if
    kc = 2.0_WP*acos(-1.0_WP)*(Uk%nx-1)/Uk%Lx
    eps = dk/4.0

    ! ===== Compute and write spectrum =====

    !WN_norm = -1
    spectrum = 0
    Lintspectrum=0
    ! Compute
    do k = UK%zmin,UK%zmax
        do j =  UK%ymin,UK%ymax
            do i =  UK%xmin,UK%xmax
                ! Compute norm of wave number
                kk = VelWN%kx(i)*VelWN%kx(i) + &
                    & VelWN%ky(j)*VelWN%ky(j) + &
                    & VelWN%kz(k)*VelWN%kz(k)
                kk = sqrt(real(kk))
                if ((kk>eps).and.(kk<=kc)) then
                ! Compute indice
!                     ik = nint(kk/dk) + 1
                    ! And then update spectrum
                    spectrum = spectrum + Uk%values(i,j,k)*conjg(Uk%values(i,j,k)) &
                            & + Vk%values(i,j,k)*conjg(Vk%values(i,j,k)) &
                            & + Wk%values(i,j,k)*conjg(Wk%values(i,j,k))

                    Lintspectrum = Lintspectrum + (1/kk)*(Uk%values(i,j,k)*conjg(Uk%values(i,j,k)) &
                            & + Vk%values(i,j,k)*conjg(Vk%values(i,j,k)) &
                            & + Wk%values(i,j,k)*conjg(Wk%values(i,j,k)))

                    !WN_norm(ik) = dk*real(ik-1)
                end if
            end do
        end do
    end do

    ! ===== Sum values on all processes =====
    spectrum = DoGlobalSum(spec_rank, spectrum)
    Lintspectrum = DoGlobalSum(spec_rank, Lintspectrum)

    Lint= Lintspectrum*2*acos(-1.0_WP)/spectrum
    Urms= sqrt(2.0_WP*spectrum)
    ReLint=Urms*Lint/mu


        ! Rlambda
      if (spec_rank.EQ.0)then
        if ( time .lt. eps ) then
            open(10,file='Re_lint.out',form='formatted')
            write(10,'(a)') "#time, ReLint,Lint,Urms"
        else
            open(10,file='Re_lint.out',form='formatted',position='append')
        end if
             write(10,'(4(g15.8,1x))') time, ReLint,Lint,Urms
        close(10)
     endif
    success = .true.
    return


    ! ===== IO error =====
!200 write(6,'(a)') 'Unable to open file - velo "'//trim(file_name)//'"'
!    return

end function compute_Integral_scale

end module post_lib
!> @}
