!> @addtogroup sub_gird_models
!! @{

!------------------------------------------------------------------------------
!
! MODULE: opt_estim_models 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to select the model for the right field
!> 
!------------------------------------------------------------------------------


module opt_estim_models 

 use datalayout
 use precision_tools 
 use scamodels
 use velmodels
 use wavenumber_tools
 use conditional_mean_tools 

 implicit none
 
 private

 public :: read_models_param
 public :: init_opt_estim_models
 public :: optimalEstimatorModelsSca
 public :: del_opt_estim_models
 public :: dxideltaSijminusdzdxj
 public :: dxideltaduidxjdzdxj
 public :: dxiduidxjdzdxj
 public :: dxisdzdxi
 public :: dxideltasdzdxi

 !Number of velocity model to use
 integer,dimension(:),allocatable                   :: numCondVarGlob
 !Number of filter type for velocity 
 integer,dimension(:),allocatable                   :: sizeEstimOptGlob
 !Type of averaging
 integer,dimension(:),allocatable                   :: typeAvgGlob
 !Local number of scalar
 integer,dimension(:),allocatable                   :: numScaOptEstimGlob
 !Number of scalar with EO model 
 integer                                            :: numScaTotGlob
 !Number of filter type for scalars 
 integer, dimension(:), allocatable                 :: numFilterGlob
 !Type of input file
 logical, dimension(:), allocatable                 :: fromANNGlob 
 !Path for optimal estimator file
 character(len=str_long),dimension(:),allocatable   :: nameEOFileGlob
 !Name of conditional variable 
 character(len=str_long),dimension(:,:),allocatable :: nameCondVarGlob

contains

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> @details
!> init global variable 
!------------------------------------------------------------------------------
 subroutine init_opt_estim_models()

   !Local data

   numScaTotGlob = 0

 end subroutine init_opt_estim_models 

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> @details
!> init global variable 
!------------------------------------------------------------------------------
 subroutine del_opt_estim_models()

   !Local data
   if (numScaTotGlob .gt. 0 ) then
     deallocate(numCondVarGlob)
     deallocate(sizeEstimOptGlob)
     deallocate(typeAvgGlob)
     deallocate(numScaOptEstimGlob)
     deallocate(numFilterGlob)
     deallocate(nameEOFileGlob)
     deallocate(nameCondVarGlob)
     deallocate(fromANNGlob)
   endif

 end subroutine del_opt_estim_models 

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> @details
!> init global variable 
!> @return  success
!------------------------------------------------------------------------------
 function incre_size_arrays(numCondVarLoc,nameCondVarLoc,nameEOFileLoc,numFilterLoc,numAvgLoc,nbLineLoc,numSca,fromANNLoc) result (success)

   !I/O data
   integer,intent(in)                        :: numCondVarLoc,numSca
   integer,intent(in)                        :: numFilterLoc,numAvgLoc,nbLineLoc
   logical,intent(in)                        :: fromANNLoc 
   character(len=*),intent(in)               :: nameEOFileLoc
   character(len=*),dimension(:),intent(in)  :: nameCondVarLoc
   logical                                   :: success

   !Local data
   integer,dimension(:),allocatable                   :: tempInt
   logical,dimension(:),allocatable                   :: tempLog
   character(len=str_long),dimension(:),allocatable   :: tempChar
   character(len=str_long),dimension(:,:),allocatable :: tempTabChar
   integer                                            :: ires,i,j

   success = .false.

   if ( numScaTotGlob .eq. 0 ) then
     allocate(numCondVarGlob(1),stat=ires)
     if (ires .ne. 0) return
     numCondVarGlob(1) = numCondVarLoc
     allocate(sizeEstimOptGlob(1),stat=ires)
     if (ires .ne. 0) return
     sizeEstimOptGlob(1) = nbLineLoc
     allocate(typeAvgGlob(1),stat=ires)
     if (ires .ne. 0) return
     typeAvgGlob(1) = numAvgLoc
     allocate(numFilterGlob(1),stat=ires)
     if (ires .ne. 0) return
     numFilterGlob(1) = numFilterLoc
     allocate(nameEOFileGlob(1),stat=ires)
     if (ires .ne. 0) return
     nameEOFileGlob(1) = nameEOFileLoc
     allocate(numScaOptEstimGlob(1),stat=ires)
     if (ires .ne. 0) return
     numScaOptEstimGlob(1) = numSca
     allocate(fromANNGlob(1),stat=ires)
     if (ires .ne. 0) return
     fromANNGlob(1) = fromANNLoc 
     allocate(nameCondVarGlob(1,numCondVarLoc),stat=ires)
     if (ires .ne. 0) return
     do i=1,numCondVarLoc
       nameCondVarGlob(1,i)=nameCondVarLoc(i)
     enddo
   else
     !fromANN
     allocate(tempLog(size(fromANNGlob,1)),stat=ires)
     if (ires .ne. 0) return
     tempLog = fromANNGlob 
     deallocate(fromANNGlob)
     allocate(fromANNGlob(size(tempLog,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       fromANNGlob(i) = tempLog(i)
     enddo
     fromANNGlob(numScaTotGlob+1) = fromANNLoc 
     deallocate(tempInt)
     !numCondVarGlob
     allocate(tempInt(size(numCondVarGlob,1)),stat=ires)
     if (ires .ne. 0) return
     tempInt = numCondVarGlob
     deallocate(numCondVarGlob)
     allocate(numCondVarGlob(size(tempInt,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       numCondVarGlob(i) = tempInt(i)
     enddo
     numCondVarGlob(numScaTotGlob+1) = numCondVarLoc
     deallocate(tempInt)
     !sizeEstimOptGlob
     allocate(tempInt(size(sizeEstimOptGlob,1)),stat=ires)
     if (ires .ne. 0) return
     tempInt = sizeEstimOptGlob 
     deallocate(sizeEstimOptGlob)
     allocate(sizeEstimOptGlob(size(tempInt,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       sizeEstimOptGlob(i) = tempInt(i)
     enddo
     sizeEstimOptGlob(numScaTotGlob+1) = nbLineLoc
     deallocate(tempInt)
     !typeAvgGlob
     allocate(tempInt(size(typeAvgGlob,1)),stat=ires)
     if (ires .ne. 0) return
     tempInt = typeAvgGlob
     deallocate(typeAvgGlob)
     allocate(typeAvgGlob(size(tempInt,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       typeAvgGlob(i) = tempInt(i)
     enddo
     typeAvgGlob(numScaTotGlob+1) = numAvgLoc
     deallocate(tempInt)
     !numFilterGlob
     allocate(tempInt(size(numFilterGlob,1)),stat=ires)
     if (ires .ne. 0) return
     tempInt = numFilterGlob 
     deallocate(numFilterGlob)
     allocate(numFilterGlob(size(tempInt,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       numFilterGlob(i) = tempInt(i)
     enddo
     numFilterGlob(numScaTotGlob+1) = numFilterLoc
     deallocate(tempInt)
     !nameEOFileGlob
     allocate(tempChar( size(nameEOFileGlob,1) ),stat=ires)
     if (ires .ne. 0) return
     tempChar = nameEOFileGlob 
     deallocate(nameEOFileGlob)
     allocate(nameEOFileGlob(size(tempChar,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       nameEOFileGlob(i) = tempChar(i)
     enddo
     nameEOFileGlob(numScaTotGlob+1) = nameEOFileLoc
     deallocate(tempChar)
     !nameCondVarGlob
     allocate(tempTabChar(size(nameCondVarGlob,1),size(nameCondVarGlob,2)),stat=ires)
     if (ires .ne. 0) return
     tempTabChar = nameCondVarGlob
     deallocate(nameCondVarGlob)
     allocate(nameCondVarGlob(size(tempTabChar,1)+1,size(tempTabChar,2)),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       do j=1,numCondVarLoc
         nameCondVarGlob(i,j) = tempTabChar(i,j)
       enddo
     enddo
     do i=1,numCondVarLoc
       nameCondVarGlob(numScaTotGlob+1,i) = nameCondVarLoc(i) 
     enddo
     deallocate(tempTabChar)
     !numScaOptEstimGlob
     allocate(tempInt( size(numScaOptEstimGlob,1) ),stat=ires)
     if (ires .ne. 0) return
     tempInt = numScaOptEstimGlob
     deallocate(numScaOptEstimGlob)
     allocate(numScaOptEstimGlob(size(tempInt,1)+1),stat=ires)
     if (ires .ne. 0) return
     do i=1,numScaTotGlob
       numScaOptEstimGlob(i) = tempInt(i)
     enddo
     numScaOptEstimGlob(numScaTotGlob+1) = numSca
     deallocate(tempInt)
   endif

   numScaTotGlob = numScaTotGlob + 1
   success = .true.

 end function incre_size_arrays

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> initialisation des tableaux recevant les estimateurs optimaux
!> selection des variables de conditionnement d'après l'input
!> @param [in]    numSca longitudinal velocity in physical space
!> @return        success : logical equal to .true. if ok and .false. otherwize
!------------------------------------------------------------------------------

 function read_models_param(numSca) result(success)

   use precision_tools
   use mpilayout_tools , only : getnbmycpu

   !I/O data
   integer,intent(in) :: numSca
   logical            :: success

   !Local data
   integer :: i,ires
   integer :: spec_rank 
   integer :: numCondVarLoc 
   integer :: numFilterLoc
   integer :: numAvgLoc
   integer :: nbLineLoc
   character(len=str_long),dimension(:),allocatable :: nameCondVarLoc
   character(len=str_long)    :: nameEOFileLoc
   character(len=str_medium)  :: request
   character(len=str_medium)  :: requestFile
   character(len=str_medium)  :: requestANN
   logical                    :: fromANNLoc


   success = .false.
   fromANNLoc = .false.
   spec_rank = getnbmycpu()
   !Number of conditionnal variables for 
   write(request,'(a,i0)') 'Number of conditionnal variables for scalar ',numSca
   call parser_read(request,numCondVarLoc)
   !INFO
   if (spec_rank .eq. 0 ) write (6,'(a,1x,a,1x,a,1x,i0)') '[INFO LES EO]',trim(adjustl(request)),'is',numCondVarLoc
   !Name of conditional variable
   allocate(nameCondVarLoc(numCondVarLoc),stat=ires)
   if (ires .ne. 0) return
   do i = 1 , numCondVarLoc
     write(request,'(a,1x,i0,1x,a,1x,i0)') 'Conditional value',i,'for scalar',numSca
     if (.not. parser_is_defined(request)) then
       write (6,'(a)') '[ERROR] A conditional value is not defined'
       return
     endif
     call parser_read(request,nameCondVarLoc(i))
     !INFO
     if (spec_rank .eq. 0 ) write (6,'(a,1x,a,1x,a,1x,a)') '[INFO LES EO]',trim(adjustl(request)),'is',trim(adjustl(nameCondVarLoc(i)))
   enddo
   !Filter for scalar
   if (.not. parseFilter(numFilterLoc,spec_rank,numSca)) return
   !Path of file of optimal estimator
   write(requestFile,'(a,1x,i0)') 'Path of file of optimal estimator for scalar',numSca
   write(requestANN,'(a,1x,i0)') 'Path of file of ANN optimal estimator for scalar',numSca
   if (parser_is_defined(trim(adjustl(requestFile)))) then
     fromANNLoc = .false.
     call parser_read(requestFile,nameEOFileLoc)
   endif
   !Path of file of optimal estimator build by ANN
   if (parser_is_defined(trim(adjustl(requestANN)))) then
     fromANNLoc = .true.
     call parser_read(requestANN,nameEOFileLoc)
   endif
   !INFO
   if (spec_rank .eq. 0 ) write (6,'(a,1x,a,1x,a,1x,a)') '[INFO LES EO]',trim(adjustl(request)),'is',trim(adjustl(nameEOFileLoc))
   !Type of averaging
   write(request,'(a,i0)') 'Type of averaging for Scalar ',numSca
   call parser_read(request,numAvgLoc)
   !INFO
   if (spec_rank .eq. 0 ) write (6,'(a,1x,a,1x,a,1x,i0)') '[INFO LES EO]',trim(adjustl(request)),'is',numAvgLoc
   !Nombre de points dans le fichier de l'estim opt
   if (.not. computeNbLine(nameEOFileLoc,nbLineLoc)) return
   !Remplissage des variables globales du module
   if (.not. incre_size_arrays(numCondVarLoc,nameCondVarLoc,nameEOFileLoc,numFilterLoc,numAvgLoc,nbLineLoc,numSca,fromANNLoc)) then
     return
   endif

   deallocate(nameCondVarLoc)
   success = .true.

 end function read_models_param

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> calcul de toutes les variables de conditionnement demandées dans l'input avec les grandeurs sous mailles
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    numSca 
!> @param [in]    condVarGlob is an array and contains conditionnal variable 
!> @param [inout] res logical
!> @return
!------------------------------------------------------------------------------

 function compute_conditionnal_var(Scalar,Uk,Vk,Wk,ScalK,ScalWN,numSca,condVarGlob) result(success)

   use mpilayout_tools , only : getnbcpus,getnbmycpu
   use datalayout
   use wavenumber_tools

   !I/O data
   type(real_data_layout),dimension(:),intent(inout)  :: condVarGlob
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   integer,intent(in)                   :: numSca
   logical                              :: success

   !Local data
   integer                                          :: nbcpus,var
   integer                                          :: spec_rank

   success = .false.
   nbcpus = getnbcpus()
   spec_rank = getnbmycpu()
   do var=1,numCondVarGlob(numSca)
!dxi sdzdxi
     if ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "sdzdxi" ) then
       if (.not. dxisdzdxi(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,nbcpus,condVarGlob(var)) ) return
!dxi deltasdzdxi
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "deltasdzdxi" ) then
        if (.not. dxideltasdzdxi(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,nbcpus,condVarGlob(var)) ) return
!dxi duidxjdzdxj
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "duidxjdzdxj" ) then
         if (.not. dxiduidxjdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,condVarGlob(var)) ) return
!dxi deltaduidxjdzdxj
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "deltaduidxjdzdxj" ) then
         if (.not. dxideltaduidxjdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,condVarGlob(var)) ) return
!!dxideltaSijminusdzdxi
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "deltaSijminusdzdxi" ) then
         if (.not. dxideltaSijminusdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,condVarGlob(var)) ) return
!varZ
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "varz" ) then
        condVarGlob(var)%values = computeFieldVar(Scalar,spec_rank)
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "varzxslice" ) then
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "varzyslice" ) then
     elseif ( trim(adjustl(nameCondVarGlob(numSca,var))) .eq. "varzzslice" ) then
     else
        if (spec_rank .eq. 0) write (6,'(a,1x,i0,1x,a,1x,i0,1x,a,a)') '[ERROR] in compute_conditionnal_var: unknown name&
                                              & for conditional variable',var,'for the',numSca,'th scalar&
                                              & computed with optimal estimator : ',trim(adjustl(nameCondVarGlob(numSca,var)))
        return
     endif
   enddo
   success = .true.

 end function compute_conditionnal_var


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> calcul de toutes les variables de conditionnement demandées dans l'input avec les grandeurs sous mailles
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [inout] divSgsOut
!> @return        res
!------------------------------------------------------------------------------
 function dxideltaSijminusdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,divSgsOut) result(res)

   use datalayout
   use scamodels , only : RGMSca_analyt
   use transforms_tools , only : btran

   !I/O data
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   type(real_data_layout),intent(inout) :: divSgsOut
   integer,intent(in)                   :: spec_rank
   logical                              :: res
   !Local data
   type(complex_data_layout)            :: tempK
   logical                              :: resLoc

   res = .false.
   if ( .not. copyStructOnly(ScalK,tempK) ) return
   call RGMSca_analyt(tempK,Scalar,Uk,Vk,Wk,ScalK,ScalWN,resLoc,spec_rank,1.0_WP)
   if (.not. resLoc) return
   call btran(tempK,divSgsOut,resLoc)
   if (.not. resLoc) return
   call deletedatalayout(tempK)
   res = .true.

 end function dxideltaSijminusdzdxj
!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> calcul de toutes les variables de conditionnement demandées dans l'input avec les grandeurs sous mailles
!>
!> @details
!> Sub-grid model based on NDCM one variable
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [inout] divSgsOut
!> @return        res
!------------------------------------------------------------------------------
 function dxideltaduidxjdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,divSgsOut) result(res)

   use datalayout
   use filtering_tools  , only : computeDelta
   use scamodels        , only : GradientSca
   use transforms_tools , only : btran

   !I/O data
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   type(real_data_layout),intent(inout) :: divSgsOut
   logical                              :: res
   !Local data
   type(complex_data_layout)            :: tempK
   logical                              :: resLoc

!        if ( .not. copyStructOnly(ScalK,tempK) ) then
!          write (6,'(a)') '[ERROR] in compute_conditionnal_var: not able to allocate memory'
!          return
!        endif
!        call GradientSca(tempK,Scalar,Uk,Vk,Wk,ScalK,Scalar,ScalWN,res,1,1.0_WP)
!        call btran(tempK,condVarGlob(var),res)
!        if (.not. res) return
!        condVarGlob(var)%values = computeDelta(Scalar)**2 * condVarGlob(var)%values
!        call deletedatalayout(tempK)
   res = .false.
   if ( .not. copyStructOnly(ScalK,tempK) ) return
   call GradientSca(tempK,Scalar,Uk,Vk,Wk,ScalK,Scalar,ScalWN,resLoc,1,1.0_WP)
   if (.not. resLoc) return
   call btran(tempK,divSgsOut,resLoc)
   if (.not. resLoc) return
   divSgsOut%values = computeDelta(Scalar)**2 * divSgsOut%values
   call deletedatalayout(tempK)
   res = .true.

 end function dxideltaduidxjdzdxj


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> 
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [inout] divSgsOut
!> @return        res
!------------------------------------------------------------------------------
 function dxiduidxjdzdxj(Scalar,Uk,Vk,Wk,ScalK,ScalWN,divSgsOut) result(res)

   use datalayout
   use scamodels        , only : GradientSca
   use transforms_tools , only : btran

   !I/O data
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   type(real_data_layout),intent(inout) :: divSgsOut
   logical                              :: res
   !Local data
   type(complex_data_layout)            :: tempK
   logical                              :: resLoc

!        if ( .not. copyStructOnly(ScalK,tempK) ) then
!          write (6,'(a)') '[ERROR] in compute_conditionnal_var: not able to allocate memory'
!          return
!        endif
!        call GradientSca(tempK,Scalar,Uk,Vk,Wk,ScalK,Scalar,ScalWN,res,1,1.0_WP)
!        call btran(tempK,condVarGlob(var),res)
!        if (.not. res) return
!        call deletedatalayout(tempK)
   res = .false.
   if ( .not. copyStructOnly(ScalK,tempK) ) return
   call GradientSca(tempK,Scalar,Uk,Vk,Wk,ScalK,Scalar,ScalWN,resLoc,1,1.0_WP)
   if (.not. resLoc) return
   call btran(tempK,divSgsOut,resLoc)
   if (.not. resLoc) return
   call deletedatalayout(tempK)
   res = .true.

 end function dxiduidxjdzdxj

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> 
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    numSca 
!> @param [in]    condVarGlob is an array and contains conditionnal variable 
!> @param [inout] res logical
!> @return
!------------------------------------------------------------------------------
 function dxisdzdxi(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,nbcpus,divSgsOut) result(res)

   use datalayout
   use subgrid_tools , only : computeSdZ
   use differential_tools 
   use toolbox

   !I/O data
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   type(real_data_layout),intent(inout) :: divSgsOut
   integer,intent(in)                   :: nbcpus
   integer,intent(in)                   :: spec_rank
   logical                              :: res

   !Local data
   type(real_data_layout),dimension(:),allocatable  :: temp
   logical                                          :: resLoc

   res = .false.
!        if ( .not. initWorkArray(Scalar,3,temp) ) then
!          write (6,'(a)') '[ERROR] in compute_conditionnal_var: not able to allocate memory'
!          return
!        endif
!        call computeSdZ(Scalar,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)
!        if (.not. res) return
!        if (.not. computeFieldDivergence(ScalWN,temp(1),temp(2),temp(3),condVarGlob(var),nbcpus,spec_rank)) return
!        res = deleteWorkArray(temp)
   if ( .not. initWorkArray(Scalar,3,temp) ) return
   call computeSdZ(Scalar,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),resLoc)
   if (.not. resLoc) return
   if (.not. computeFieldDivergence(ScalWN,temp(1),temp(2),temp(3),divSgsOut,nbcpus,spec_rank)) return
   if (.not. deleteWorkArray(temp)) return
   res = .true.

 end function dxisdzdxi


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> 
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    Scalar scalar in physical space 
!> @param [in]    Scalk scalar in spectral space 
!> @param [in]    ScalWN wave numbers for scalar 
!> @param [in]    spec_rank
!> @param [in]    nbcpus
!> @param [inout] divSgsOut
!> @return        res
!------------------------------------------------------------------------------

 function dxideltasdzdxi(Scalar,Uk,Vk,Wk,ScalK,ScalWN,spec_rank,nbcpus,divSgsOut) result(res)

   use datalayout
   use filtering_tools , only : computeDelta
   use subgrid_tools   , only : computeSdZ
   use differential_tools 
   use toolbox

   !I/O data
   type(complex_data_layout),intent(in) :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)         :: ScalWN
   type(real_data_layout),intent(in)    :: Scalar
   type(real_data_layout),intent(inout) :: divSgsOut
   integer,intent(in)                   :: nbcpus
   integer,intent(in)                   :: spec_rank
   logical                              :: res

   !Local data
   type(real_data_layout),dimension(:),allocatable  :: temp
   logical                                          :: resLoc

   res = .false.
!        if ( .not. initWorkArray(Scalar,3,temp) ) then
!          write (6,'(a)') '[ERROR] in compute_conditionnal_var: not able to allocate memory'
!          return
!        endif
!        call computeSdZ(Scalar,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),res)
!        if (.not. res) return
!        if (.not. computeFieldDivergence(ScalWN,temp(1),temp(2),temp(3),condVarGlob(var),nbcpus,spec_rank)) return
!        condVarGlob(var)%values = computeDelta(Scalar)**2 * condVarGlob(var)%values
!        res = deleteWorkArray(temp)
   if ( .not. initWorkArray(Scalar,3,temp) ) return
   call computeSdZ(Scalar,Uk,Vk,Wk,ScalK,ScalWN,temp(1),temp(2),temp(3),resLoc)
   if (.not. resLoc) return
   if (.not. computeFieldDivergence(ScalWN,temp(1),temp(2),temp(3),divSgsOut,nbcpus,spec_rank)) return
   divSgsOut%values = computeDelta(Scalar)**2 * divSgsOut%values
   if (.not. deleteWorkArray(temp)) return
   res = .true.

 end function dxideltasdzdxi

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> construction du terme sgs avec les variable calculées et l'estiamateur optimale chargé depuis le fichier
!> @param [in]    U  longitudinal velocity in physical space
!> @param [in]    V  velocity in physical space
!> @param [in]    W  velocity in physical space
!> @param [in]    Uk longitudinal velocity in spectral space
!> @param [in]    Vk velocity in spectral space
!> @param [in]    Wk velocity in spectral space
!> @param [in]    ScalIn scalar in physical space 
!> @param [in]    ScalkIn scalar in spectral space 
!> @param [in]    WaveNum wave numbers for scalar 
!> @param [in]    Filter is the type of test filter 
!> @param [in]    spec_rank is the rank of processus
!> @param [inout] res logical
!> @param [inout] div is the divergence of SGS model
!> @return
!------------------------------------------------------------------------------
 function compute_sgs_term(Scalar,dTidxiK,numSca,condVarGlob) result(success)

   use mpilayout_tools , only : getnbmycpu
   use fann

   !I/O data
   type(real_data_layout),intent(in)       :: Scalar
   type(complex_data_layout),intent(inout) :: dTidxiK
   integer,intent(in)                      :: numSca
   type(real_data_layout),dimension(:),intent(in)  :: condVarGlob
   logical                                 :: success

   !Local data
   real(WP),dimension(:,:),allocatable  :: tabEOFile
   type(real_data_layout)               :: dTidxi,dTidxiVerifTrans
   integer                              :: ires
   integer                              :: spec_rank
   logical                              :: res 

   success = .false.
   spec_rank = getnbmycpu()
   if ( .not. copyStructOnly(Scalar,dTidxi) .or. &
      & .not. copyStructOnly(Scalar,dTidxiVerifTrans) ) then
     write (6,'(a)') '[ERROR] not enought memory'
     return
   endif
   if (fromANNGlob(numSca)) then
     !Compute optimal estimator with ANN 
     if (.not. computeEOFromANN(condVarGlob,dTidxi,nameEOFileGlob(numSca))) return
   else
     !Compute optimal estimator with histogram technic
     if (spec_rank .eq. 0)  write (6,'(a)') '[INFO compute_sgs_term] allocating array...'
     allocate( tabEOFile(numCondVarGlob(numSca)+1,sizeEstimOptGlob(numSca)),stat = ires)
     if (ires .ne. 0) then
       write (6,'(a)') '[ERROR] not enought memory'
       return
     endif
     tabEOFile = 0.0_WP
     if (spec_rank .eq. 0)  write (6,'(a)') '[INFO compute_sgs_term] reading file...'
     if (.not. readModelOpt(trim(adjustl(nameEOFileGlob(numSca))),tabEOFile)) then
       write (6,'(a)') '[ERROR] Not able to read file readMod'
     endif
     if (spec_rank .eq. 0)  write (6,'(a)') '[INFO compute_sgs_term] building field...'
     if (.not. buildEOFieldFromEO(condVarGlob,tabEOFile,dTidxi)) then 
       write (6,'(a)') '[ERROR] buildEOFieldFromEO in compute_sgs_term'
       return
     endif
     deallocate(tabEOFile)
   endif
   call ftran(dTidxi,dTidxiK,res)
   if ( .not. res ) then 
     write (6,'(a)') '[ERROR] the ftran in compute_sgs_term failed because of a bad signal&
                     & from optimal estimator file'
     return
   endif
   !Verification de la FFT du signal de l'EO
   call btran(dTidxiK,dTidxiVerifTrans,res)
   if (.not. equalToVal(dTidxiVerifTrans,dTidxi)) then
     write (6,'(a)') '[WARNING] The signal build by optimal estimator is too bad so fft is not able to transform it correctly!'
     return
   endif 
   !free memory
   call deletedatalayout(dTidxi)
   call deletedatalayout(dTidxiVerifTrans)
  
   success = .true.

 end function compute_sgs_term

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> routine d'interface avec le module sub_grid_models
!> @param [in]    sca is the number of scalar in global enumeration 
!> @param [in]    Scalar i the scalar in physical space
!> @param [in]    dTidxiK is the sgs 
!> @param [inout] res logical
!------------------------------------------------------------------------------

 subroutine optimalEstimatorModelsSca(dTidxiK,Scalar,Uk,Vk,Wk,ScalK,ScalWN,numSca,res)

   use mpilayout_tools , only : getnbmycpu

   !I/O data
   type(complex_data_layout),intent(in)    :: Uk,Vk,Wk,ScalK 
   type(wavenumbers),intent(in)            :: ScalWN
   integer,intent(in)                      :: numSca
   type(real_data_layout),intent(in)       :: Scalar
   type(complex_data_layout),intent(inout) :: dTidxiK
   logical,intent(out)                     :: res

   !Local data
   integer                                 :: indice,i
   integer                                 :: spec_rank
   character(len=16)                       :: charIndice
   !Array for conditional variables
   type(real_data_layout),dimension(:),allocatable  :: condVarGlob

   res = .true.
   spec_rank = getnbmycpu()
   indice = findNumSca(numSca)
   if (.not. initWorkArray(Scalar,numCondVarGlob(indice),condVarGlob) ) then
     if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] in optimalEstimatorModelsSca: not able to allocate memory'
     res = .false.
   endif
   do i= 1 , numCondVarGlob(indice)
     write (charIndice,'(i2.2)') i
     condVarGlob(i)%name=trim(adjustl(Scalar%name))//'_VarCond_'//trim(adjustl(charIndice))
   enddo
   if (.not. compute_conditionnal_var(Scalar,Uk,Vk,Wk,ScalK,ScalWN,indice,condVarGlob) .or. &
       .not. compute_sgs_term(Scalar,dTidxiK,indice,condVarGlob) ) then
     if (spec_rank .eq. 0) write (6,'(a)') '[ERROR] in optimalEstimatorModelsSca: not able to compute sgs term'
     res = .false.
   endif
 
   if (.not. deleteWorkArray(condVarGlob)) then
     res = .false.
   endif

 end subroutine optimalEstimatorModelsSca

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> @param [in]    sca is the number of scalar in global enumeration 
!> @return        indice 
!------------------------------------------------------------------------------
 function findNumSca(sca) result(indice)

   !I/O data
   integer,intent(in)  :: sca
   integer             :: indice

   !Local data
   integer :: i
   indice = 0
   do i=1,numScaTotGlob
     if ( numScaOptEstimGlob(i) .eq. sca ) indice = i
   enddo

 end function findNumSca
 
end module opt_estim_models 
!> @}
