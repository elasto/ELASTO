!> @addtogroup post_process
!! @{

!-------------------------------------------------------------------------------------------
! This module provides tool to post-process several input files with several scalars and mhd
! fields
!-------------------------------------------------------------------------------------------
module post_out_simu

  use conditional_mean_tools
  use param
  use precision_tools
  use datalayout
  use data
  use post_lib
  use post_scalar
  use post_velocity
  use wavenumber_tools
  use vtkxml_plus
  use cart_topology
  use post_ref
  use stat_tools
  use post_hd_out
  use post_visus
  use post_mhd_out
  use energy_transfers

  implicit none

  private

  !Number of files
  integer            :: iFile=0
  !Number of scalar
  integer            :: iScal=0
  !Number of scalar part
  integer            :: iScalPart=0
  !Flag which determine the post process
  integer            :: iPostprocess=0
  !> Type of paraview output
  logical            :: paraview_3D, paraview_slice
  !Min of test filter
  integer            :: iMinFilter
  !Max of test filter
  integer            :: iMaxFilter
  !Test filter type
  !character(len=64)  :: iFilterType
  integer            :: iFilterType

  character(len=64)  :: namepost
  integer,dimension(:),allocatable :: tagScal
  integer,dimension(:),allocatable :: tagScalPart
  integer,dimension(3)             :: tagVel
  integer,dimension(3)             :: tagB


  !Public routines
  public :: post_out_simu_init,post_out_simu_init_luca
  public :: post_out_simu_del
  public :: post_out

contains

!-------------------------------------------------------------------------------------------
!> Initialise the context of poost-processing out of the solver
!!
!! @author Antoine Volllant
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbscal          = the number of scalar
!!    @param[in]    nbscal_part     = the number of scalar solved with particle methods
!!    @return       res             = logical value
!-------------------------------------------------------------------------------------------
 function post_out_simu_init(spec_rank,nbscal_part,nbscal,scal_partArray,scalArray,U,V,W,B,imhd) result(res)

   use io_interface
   use vtkxml_plus
   use cart_topology

   !I/O
   type(REAL_DATA_LAYOUT), intent(in),pointer,dimension(:) :: scalArray, scal_partArray
   type(REAL_DATA_LAYOUT), intent(in)                      :: U,V,W
   type(REAL_DATA_LAYOUT),pointer, intent(in), dimension(:):: B
   integer,intent(in)                    :: spec_rank
   integer,intent(in)                    :: nbscal_part
   integer,intent(in)                    :: nbscal
   logical,intent(in)                    :: imhd
   logical                               :: res

   !Local data
   real(WP),dimension(3)  :: length_3D
   integer                :: nbField,i
   character(len=64)      :: request
   character (len=10)     :: type_para       ! For choosing if 3D, slice or both

   res = .false.

   length_3D(1)=U%Lx
   length_3D(2)=U%Ly
   length_3D(3)=U%Lz
   !Check if the input file has informations
   IF     (.NOT. parser_is_defined('Number of pool of file')) THEN
     WRITE(6,'(a)')'[Warning] Number of pool of file is not set: abort !!'
     RETURN
   ENDIF
   IF (.NOT. parser_is_defined('Post-processing')) THEN
     WRITE(6,'(a)')'[Warning] Post-processing is not set: abort !!'
     RETURN
   ENDIF

   !Informations about scalars and part_scalar
   iScal=nbscal
   iScalPart=nbscal_part

   !Informations about files to compute
   call parser_read('Number of pool of file',iFile)
   !Determine the good post processing
   call parser_read('Post-processing',namepost)
   !min value for filtering
   call parser_read('filter size min for a priori test',iMinFilter)
   !max value for filtering
   call parser_read('filter size max for a priori test',iMaxFilter)
   !type of filter for tests
   write(request,'(a)') 'type of filter for a priori test'
   if (.not. parseFilter(iFilterType,spec_rank,other=request) ) return 
   if (.not. interface_io_init(spec_rank)) return

   if (spec_rank .eq. 0 ) write(6,'(a,1x,a,1x,a)')'[INFO]',trim(namepost),'asked ...'
   if      ( namepost .eq. 'double') then
     iPostprocess = 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'momentScalarIncrement') then 
     iPostprocess = 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'avs') then 
     iPostprocess = 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'paraview') then
      call parser_read('slice or 3D', type_para)
      select case(trim(type_para))
      case('slice')
          paraview_slice = .true.
      case('3D')
          paraview_3D = .true.
      case default
          paraview_slice = .true.
          paraview_3D = .true.
      end select
     if (imhd) then
       nbField = nbscal+nbscal_part+6
     else
       nbField = nbscal+nbscal_part+3
     endif
     call vtkxml_init_all(nbField, nb_proc_dim , length_3D, spec_rank , coord)
     if (nbscal .gt. 0) then
        allocate(tagScal(nbscal))
        do i=1,nbscal
          call vtkxml_plus_init_field(scalArray(i), tagScal(i) )
        enddo
     endif
     if (nbscal_part .gt. 0) then 
        allocate(tagScalPart(nbscal_part))
        do i=1,nbscal_part
          call vtkxml_plus_init_field(scal_partArray(i), tagScalPart(i) )
        enddo
     endif
     if (imhd) then 
       do i=1,3
          call vtkxml_plus_init_field(B(i), tagB(i) )
        enddo
     endif
     call vtkxml_plus_init_field(U, tagVel(1) )
     call vtkxml_plus_init_field(V, tagVel(2) )
     call vtkxml_plus_init_field(W, tagVel(3) )
     iPostprocess = 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'rewrite_output') then
     iPostprocess = 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'dumpSGSFields') then
     iPostprocess = 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'testsikmoinskjsij') then
     iPostprocess = 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'dumpUScalaireotdTauijdxi') then 
     iPostprocess = 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'CTR2012_postprocess') then 
     iPostprocess = 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'CTR2012_postprocess_vel') then
     iPostprocess = 91
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'checkGMVel') then
     iPostprocess = 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'spectrum avg') then
     iPostprocess = 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioriRGMVel') then
     iPostprocess = 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'testTransfertNrjQDM') then
     iPostprocess = 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif (namepost .eq. 'post_high_schmidt_decay') then
     iPostprocess = 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif (namepost .eq. 'compareVelocityFields') then
     iPostprocess = 15
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestscalar') then
     iPostprocess = 27
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'transDiscretize') then
     iPostprocess = 28
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestgradient') then
     iPostprocess = 29
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestsmagorinskydyn') then
     iPostprocess = 30
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestclark') then
     iPostprocess = 31
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestsmagorinsky') then
     iPostprocess = 32
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestclarkfabre') then
     iPostprocess = 33
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestRGMmod') then
     iPostprocess = 34
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestRGM') then
     iPostprocess = 35
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postHDout') then
     iPostprocess = 36
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postvisuHkin') then
     iPostprocess = 37
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postMHDout') then
     iPostprocess = 38

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'spectral transfers') then
     iPostprocess = 39
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'real_interpol') then
     iPostprocess = 40

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   elseif  ( namepost .eq. 'none') then
     write(6,'(a)')'[Warning] You asked a post process and you did not set its type !!'
     iPostprocess = 100
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else
     write(6,'(a)')'[Warning] Unknow post-processing in post_out_simu_init: abort !!'
     return
   endif

   res = .true.

 end function post_out_simu_init

!-------------------------------------------------------------------------------------------
!> Initialise the context of poost-processing out of the solver
!!
!!  @author 
!!  Luca MARRADI LIPhy 2014
!!  @param[in]    spec_rank       = mpi rank inside spectral communicator
!!  @param[in]    sigma           = stress tensor 
!!  @return       res             = logical value
!-------------------------------------------------------------------------------------------
 FUNCTION post_out_simu_init_luca(spec_rank,sigma) result(res)

   use io_interface
   use vtkxml_plus
   use cart_topology

   !== INPUT/OUTPUT DATA ==
   type(REAL_DATA_LAYOUT), intent(in)    :: sigma
   integer,intent(in)                    :: spec_rank
   logical                               :: res

   !== LOCAL DATA ==
   real(WP),dimension(3)  :: length_3D
   integer                :: nbField,i
   character(len=64)      :: request
   character (len=10)     :: type_para       ! For choosing if 3D, slice or both

   res = .false.

   length_3D(1)=sigma%Lx
   length_3D(2)=sigma%Ly
   length_3D(3)=sigma%Lz

   !Check if the input file has informations
   IF     (.NOT. parser_is_defined('Number of pool of file')) THEN
     WRITE(6,'(a)')'[Warning] Number of pool of file is not set: abort !!'
     RETURN
   ENDIF

   IF (.NOT. parser_is_defined('Post-processing')) THEN
     WRITE(6,'(a)')'[Warning] Post-processing is not set: abort !!'
     RETURN
   ENDIF

   !Informations about files to compute
   call parser_read('Number of pool of file',iFile)
   !Determine the good post processing
   call parser_read('Post-processing',namepost)
   !min value for filtering
   call parser_read('filter size min for a priori test',iMinFilter)
   !max value for filtering
   call parser_read('filter size max for a priori test',iMaxFilter)
   !type of filter for tests
   write(request,'(a)') 'type of filter for a priori test'
   if (.not. parseFilter(iFilterType,spec_rank,other=request) ) return 
   if (.not. interface_io_init(spec_rank)) return

   if (spec_rank .eq. 0 ) write(6,'(a,1x,a,1x,a)')'[INFO]',trim(namepost),'asked ...'
   if      ( namepost .eq. 'double') then
     iPostprocess = 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'momentScalarIncrement') then 
     iPostprocess = 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'avs') then 
     iPostprocess = 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'paraview') then
      call parser_read('slice or 3D', type_para)
      select case(trim(type_para))
      case('slice')
          paraview_slice = .true.
      case('3D')
          paraview_3D = .true.
      case default
          paraview_slice = .true.
          paraview_3D = .true.
      end select

     call vtkxml_plus_init_field(sigma, tagVel(1) )
!!     call vtkxml_plus_init_field(V, tagVel(2) )
!!     call vtkxml_plus_init_field(W, tagVel(3) )
     iPostprocess = 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'rewrite_output') then
     iPostprocess = 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'dumpSGSFields') then
     iPostprocess = 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'testsikmoinskjsij') then
     iPostprocess = 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'dumpUScalaireotdTauijdxi') then 
     iPostprocess = 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'CTR2012_postprocess') then 
     iPostprocess = 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'CTR2012_postprocess_vel') then
     iPostprocess = 91
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'checkGMVel') then
     iPostprocess = 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'spectrum avg') then
     iPostprocess = 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioriRGMVel') then
     iPostprocess = 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'testTransfertNrjQDM') then
     iPostprocess = 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif (namepost .eq. 'post_high_schmidt_decay') then
     iPostprocess = 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif (namepost .eq. 'compareVelocityFields') then
     iPostprocess = 15
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestscalar') then
     iPostprocess = 27
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'transDiscretize') then
     iPostprocess = 28
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestgradient') then
     iPostprocess = 29
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestsmagorinskydyn') then
     iPostprocess = 30
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestclark') then
     iPostprocess = 31
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestsmagorinsky') then
     iPostprocess = 32
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestclarkfabre') then
     iPostprocess = 33
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestRGMmod') then
     iPostprocess = 34
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'aprioritestRGM') then
     iPostprocess = 35
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postHDout') then
     iPostprocess = 36
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postvisuHkin') then
     iPostprocess = 37
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'postMHDout') then
     iPostprocess = 38

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'spectral transfers') then
     iPostprocess = 39
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif  ( namepost .eq. 'real_interpol') then
     iPostprocess = 40

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   elseif  ( namepost .eq. 'none') then
     write(6,'(a)')'[Warning] You asked a post process and you did not set its type !!'
     iPostprocess = 100
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else
     write(6,'(a)')'[Warning] Unknow post-processing in post_out_simu_init: abort !!'
     return
   endif

   res = .TRUE.

 END FUNCTION post_out_simu_init_luca
!-------------------------------------------------------------------------------------------
!> Perform post-processing whitout the main solver. This function calls subroutines 
!!
!! @author Antoine Volllant 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    nbcpus          = number of process in the cpu pool 
!!    @param[in]    U               = Component of velocity 
!!    @param[in]    V               = Component of velocity 
!!    @param[in]    W               = Component of velocity
!!    @param[in]    B               = Array for magnetic fields components
!!    @param[in]    scalArray       = Array which contains scalar
!!    @param[in]    scal_partArray  = Array which contains scalar solved with particle methods
!!    @param[in]    imhd            = logical value for the existence of mhd 
!!    @param[in]    nbscal          = the number of scalar
!!    @param[in]    nbscal_part     = the number of scalar solved with particle methods
!!    @return       res             = logical value 
!-------------------------------------------------------------------------------------------
 function post_out (nbcpus ,spec_rank, nbscal, scalArray, nbscal_part,scal_partArray, U, V,&
                   & W,B,imhd) result(res)

    use system_tools
    use io_interface

   !I/O
   integer,intent(in)                                        :: nbcpus
   integer,intent(in)                                        :: spec_rank
   integer,intent(in)                                        :: nbscal
   integer,intent(in)                                        :: nbscal_part
   type(real_data_layout),intent(inout)                      :: U, V ,W
   type(REAL_DATA_LAYOUT),intent(inout),pointer,dimension(:) :: B
   type(REAL_DATA_LAYOUT),intent(inout),pointer,dimension(:) :: scalArray
   type(REAL_DATA_LAYOUT),intent(inout),pointer,dimension(:) :: scal_partArray
   logical,intent(in)                                        :: imhd
   logical                                                   :: res

   !Local data
   integer :: i,j
   type(real_data_layout)                                    :: Utemp, Vtemp ,Wtemp


   res = .false.
   if (iPostprocess .eq. 100 ) return
   DO i=1,iFile
     IF (.NOT. data_read(U,V,W,B,imhd,ScalArray, Scal_partArray,nbscal,&
        & nbscal_part,spec_rank,i)) THEN
        WRITE(6,'(a,1x,i0)')'[ERROR] Unable to read files',i
        RETURN
     ENDIF
     SELECT CASE (iPostprocess)
       CASE (1)
         IF (.NOT.  postprocess_double(nbcpus ,spec_rank, U, V ,W,B,imhd)) THEN
            WRITE(6,'(a)')'[ERROR] Unable to perform post-process to double length'
            RETURN
         ENDIF

       CASE (2)
         IF (nbscal.ge.1) THEN
         DO j=1,nbscal
           IF (.NOT. computeScalarMomentsStat(scalArray(j), i, spec_rank) ) THEN
              WRITE(6,'(a)')'[ERROR] Unable to perform post-process to double length'
              RETURN
           ENDIF
         ENDDO
         ENDIF
         IF (nbscal_part.ge.1) THEN
         DO j=1,nbscal_part
           IF (.NOT. computeScalarMomentsStat(scal_partArray(j), i, spec_rank) ) THEN
              WRITE(6,'(a)')'[ERROR] Unable to perform post-process to double length'
              RETURN
           ENDIF
         ENDDO
         ENDIF
       CASE (3)
         IF (.NOT.  avsWrite(spec_rank,scalArray, scal_partArray, U, V, &
            & W,B,imhd,i,iFile)) THEN
            WRITE(6,'(a)')'[ERROR] Unable to perform post-process avsWrite'
            RETURN
         ENDIF

       CASE (4)
         IF (.NOT. paravieWrite(spec_rank,nbscal, scalArray, nbscal_part,&
            &scal_partArray, U, V, W,B,imhd,tagScal,tagScalPart,tagVel,  &
            &tagB, paraview_slice, paraview_3D) )THEN
            WRITE (6,'(a)')'[ERROR] Failed in paraview dump'
            RETURN
         ENDIF

       CASE (5)
         IF (.NOT. dataWrite(spec_rank,scalArray, scal_partArray, U, V, &
            & W,B,imhd,i)) THEN
            WRITE(6,'(a)')'[ERROR] Unable to perform post-process dataWrite'
            RETURN
         ENDIF
         IF (.NOT. interface_io_close()) WRITE(6,'(a)')'[ERROR] Unable to exit hdf5 environnement'
       CASE (6)
!         IF (.NOT. dumpSGSFieldsForScalar(i,U,V,W,ScalArray(1),iMinFilter,iMaxFilter,iFilterType) ) THEN
!            WRITE(6,'(a)')'[ERROR] Unable to perform post-process dumpSGSFieldsForScalar'
!            RETURN
!         ENDIF
       CASE (7)
          call testDissipationduidxkdujdxk(spec_rank,res,U,V,W)
          if (.not. res) write (6,'(a)')'[ERROR] in testDissipationduidxkdujdxk'
       CASE (8)
!          call dumpUScalaireRotdTauijdxi(U,V,W,spec_rank) 
       CASE (9)
#ifdef POSTCTR
          if (.not. gradientmodOE(U,V,W,ScalArray(1),nbcpus,spec_rank)) then
            write (6,'(a)')'[ERROR] Failed in gradientmodOE'
            return
          endif
#else
          write (6,'(a)')'[ERROR] PLEASE use -DPOSTCTR in makefile'
          return
#endif
       CASE (91)
          if (.not. gradientvelmodOE(U,V,W,ScalArray(1),nbcpus,spec_rank)) then
            write (6,'(a)')'[ERROR] Failed in gradientmodOE'
            return
          endif
       CASE (10)
          call checkGradientModels(spec_rank,nbcpus,res,U,V,W)
          if (.not. res) write (6,'(a)')'[ERROR] in checkGradientModels'
       CASE (11)
          if (.not. spectrum_time_avg(spec_rank)) then
            write (6,'(a)')'[ERROR] Failed in spectrum_time_avg'
          else
            write (6,'(a)')'Success in spectrum_time_avg'
          end if
       CASE (12)
!          call aprioriRGMVel(i,spec_rank,nbcpus,res,U,V,W)
          call aprioriSikmoinsSkj(i,spec_rank,nbcpus,res,U,V,W)
          if (.not. res) write (6,'(a)')'[ERROR] in aprioriRGMVel'
          !VIDE
       CASE (13)
          call testTransfertNrjQDM(spec_rank,res,U,V,W) 
       CASE (14)
           if (.not. post_THI_spectrum_decay(U,V,W,ScalArray(1),spec_rank,nbcpus)) then
            write (6,'(a)')'[ERROR] in post_high_schmidt_decay'
           end if
       CASE (15)
           IF (iFile .ne. 2) THEN
            write (6,'(a)')'[ERROR] compareVelocityField : number of fields is not exact'
           ENDIF
           IF ( i .EQ. 1) THEN
             IF(.NOT. copyStructonly(U,Utemp)) RETURN
             IF(.NOT. copyStructonly(V,Vtemp)) RETURN
             IF(.NOT. copyStructonly(W,Wtemp)) RETURN
             Utemp%values = U%values
             Vtemp%values = V%values
             Wtemp%values = W%values
           ELSE
             CALL compareVeloFields(Utemp,Vtemp,Wtemp,U,V,W,spec_rank,res)
             IF (.NOT. res) write (6,'(a)')'[ERROR] compareVelocityField failed !'
            CALL deletedatalayout(Utemp)
             CALL deletedatalayout(Vtemp)
             CALL deletedatalayout(Wtemp)
           ENDIF
       CASE (16)
          !VIDE
       CASE (17)
          !VIDE
       CASE (18)
          !VIDE
       CASE (19)
          !VIDE
       CASE (20)
          !VIDE
       CASE (21)
          !VIDE
       CASE (22)
          !VIDE
       CASE (23)
          !VIDE
       CASE (24)
          !VIDE
       CASE (25)
          !VIDE
       CASE (26)
          !VIDE
       CASE (27)
           if (.not. apriori_ScalarStat(i,U,V,W,scalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType) ) then
             write (6,'(a)')'[ERROR] in computeScalarStat'
           endif

       CASE (28)
           if (.not. transDiscAvsWrite(spec_rank,scalArray,scal_partArray,U,V,W,B,imhd,i,iFile) ) then
             write (6,'(a)')'[ERROR] in computeScalarStat'
           endif

       CASE (29)
           if (.not. apriori_Gradient(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_Gradient'
           endif

       CASE (30)
           if (.not. apriori_SmagorinskyDyn(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_SmagDyn'
           endif

       CASE (31)
           if (.not. apriori_Clark(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_Clark'
           endif

       CASE (32)
           if (.not. apriori_Smagorinsky(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_Smag'
           endif

       CASE (33)
           if (.not. apriori_ClarkFabre(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_ClarkFabre'
           endif

       CASE (34)
           if (.not. apriori_RGM_analyt(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_RGM_mod'
           endif

       CASE (35)
           if (.not. apriori_RGM(i,U,V,W,ScalArray(1),spec_rank,iMinFilter,iMaxFilter,iFilterType)   ) then 
             write (6,'(a)')'[ERROR] in apriori_RGM'
           endif

       CASE (36)
             call computeallthetermsHDout(U,V,W,&
                         & nbcpus,spec_rank)
       CASE (37)
             call Helicity_Fieldsave(U,V,W,&
                         & nbcpus,spec_rank)
       CASE (38)
             call computeallthetermsMHDout(U,V,W,B,&
                         & nbcpus,spec_rank)
       CASE (39)
             call compute_energy_transfers(U,V,W,&
                         & nbcpus,spec_rank)
       CASE (40)
          if (spec_rank==0) then
            write (6,'(a)')' [STATUS] Start velo interpol'
          end if
          if (.not. velo_test_interpol(U, V, W, scal_partArray(1), nbcpus, spec_rank)) then
            write (6,'(a)')'[ERROR] Failed in test_interpol (real interpolation)'
            return
          endif
!        IF (.NOT. dataWrite(spec_rank,scal_partArray, scal_partArray, U, V, &
!           & W,B,imhd,i)) THEN
!           WRITE(6,'(a)')'[ERROR] Unable to perform post-process dataWrite'
!           RETURN
!        ENDIF
!       CASE (34)
!           if (.not. apriori_ScaleSimilarity(i,U,V,W,ScalArray(1),spec_rank,8,8,1)   ) then 
!             write (6,'(a)')'[ERROR] in computeGradientAprioriTest'
!           endif

       CASE DEFAULT
         WRITE(6,'(a)')'[ERROR] Unknow post-processing in post_out: abort !!'
         RETURN
     END SELECT
   ENDDO

   res = .true.

 end function post_out



!------------------------------------------------------------------------------
!> Delete post-process setup.
!! @author Antoine Vollant
!! @param[in]    nbscal = number of scalar solved with spectral method
!! @param[in]    nbscal_part = number of scalar solved with particles method
!! @return       success = logical equals to true is everything is right.
!------------------------------------------------------------------------------
 function post_out_simu_del(nbscal,nbscal_part) result(success)
    use io_interface        

    !I/O data
    logical :: success
    integer,intent(in)  :: nbscal,nbscal_part

    success = .false.
    if ( iPostprocess .eq. 4 )then
       call vtkxml_finish
       if (nbscal .gt. 0) deallocate(tagScal)
       if (nbscal_part .gt. 0) deallocate(tagScalPart)
    endif 
    IF (.NOT. interface_io_close()) WRITE(6,'(a)')'[ERROR] Unable to exit hdf5 environnement'
    success = .true.
        
 end function post_out_simu_del


end module post_out_simu
!> @}
