!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: physical_values_tools 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to provide elementaries piece of code which can be use elsewhere.
!> BEWARE there are no systematic tests in these routines in order to not degrade performances
!> 
!------------------------------------------------------------------------------


module physical_values_tools

 use precision_tools 
 use stat_tools
 use datalayout
 use wavenumber_tools
 use differential_tools
 
 implicit none

 public

    ! ==== Interface procedures ====

    interface computeEnergy
       module procedure computeEnergyArray
       module procedure computeEnergyScal
    end interface computeEnergy 

 contains


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> compute the term nu ( (du|idxjdu|idxj)^ - du|^idxj du|^idxk )
!> @details
!> Compute the 6 components of stress tensor 
!> @param [in] Ubar
!> @param [in] Ubark
!> @param [in] Vbark
!> @param [in] Wbark
!> @param [in] wave the WaveNumbers associated to inK variable
!> @param [in,out] out11 calculated term in physical space 
!> @param [in,out] out22 calculated term in physical space 
!> @param [in,out] out33 calculated term in physical space 
!> @param [in,out] out12 calculated term in physical space 
!> @param [in,out] out13 calculated term in physical space 
!> @param [in,out] out23 calculated term in physical space 
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
 subroutine computeViscousDissResolFiltered(ViscDissResFil,U,Ubark,Vbark,Wbark,wave,nu,filter,delta,success)

   use filtering_tools , only : computeFilter 

   !I/O
   type(real_data_layout),intent(in)     :: U
   type(complex_data_layout),intent(in)  :: Ubark,Vbark,Wbark
   type(WaveNumbers),intent(in)          :: wave
   integer,intent(in)                    :: filter
   real(wp),intent(in)                   :: delta,nu
   type(real_data_layout),intent(inout)  :: ViscDissResFil
   logical, intent(inout)                :: success  

   !Local data
   type(complex_data_layout)    :: t1k,t2k,t3k
   type(real_data_layout)       :: ddx,ddy,ddz
   type(real_data_layout)       :: temp
   logical                      :: res  


   success = .true.

   if (.not. copyStructOnly(Ubark,t1k) .or. &
      &.not. copyStructOnly(Ubark,t2k) .or. &
      &.not. copyStructOnly(Ubark,t3k) )then 
       success = .false.
       write(6,'(a)')'[ERROR] in computeViscousDissResolFiltered : not enought memory !'
   endif 

   if (.not. copyStructOnly(U,ddx) .or. &
      &.not. copyStructOnly(U,ddy) .or. &
      &.not. copyStructOnly(U,ddz) .or. &
      &.not. copyStructOnly(U,temp) )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeViscousDissResolFiltered : not enought memory !'
   endif 

   call computeGradientK(wave,Ubark,t1k,t2k,t3k)
   call btran(t1k,ddx,res) 
   if (.not. res) success = .false.
   call btran(t2k,ddy,res) 
   if (.not. res) success = .false.
   call btran(t3k,ddz,res) 
   if (.not. res) success = .false.
   ddx%values = ddx%values**2 + ddy%values**2 + ddz%values**2
   call computeFilter(Wave,delta,ddx,temp,filter)
   call computeFilter(wave,delta,t1k,ddx,filter)
   call computeFilter(wave,delta,t2k,ddy,filter)
   call computeFilter(wave,delta,t3k,ddz,filter)
   ViscDissResFil%values = temp%values - ddx%values**2 - ddy%values**2 - ddz%values**2

   call computeGradientK(wave,Vbark,t1k,t2k,t3k)
   call btran(t1k,ddx,res) 
   if (.not. res) success = .false.
   call btran(t2k,ddy,res) 
   if (.not. res) success = .false.
   call btran(t3k,ddz,res) 
   if (.not. res) success = .false.
   ddx%values = ddx%values**2 + ddy%values**2 + ddz%values**2
   call computeFilter(Wave,delta,ddx,temp,filter)
   call computeFilter(wave,delta,t1k,ddx,filter)
   call computeFilter(wave,delta,t2k,ddy,filter)
   call computeFilter(wave,delta,t3k,ddz,filter)
   ViscDissResFil%values = ViscDissResFil%values + temp%values - ddx%values**2 - ddy%values**2 - ddz%values**2

   call computeGradientK(wave,Vbark,t1k,t2k,t3k)
   call btran(t1k,ddx,res) 
   if (.not. res) success = .false.
   call btran(t2k,ddy,res) 
   if (.not. res) success = .false.
   call btran(t3k,ddz,res) 
   if (.not. res) success = .false.
   ddx%values = ddx%values**2 + ddy%values**2 + ddz%values**2
   call computeFilter(Wave,delta,ddx,temp,filter)
   call computeFilter(wave,delta,t1k,ddx,filter)
   call computeFilter(wave,delta,t2k,ddy,filter)
   call computeFilter(wave,delta,t3k,ddz,filter)
   ViscDissResFil%values = ViscDissResFil%values + temp%values - ddx%values**2 - ddy%values**2 - ddz%values**2
   ViscDissResFil%values = nu*ViscDissResFil%values

   call deleteDataLayout(t1k)
   call deleteDataLayout(t2k)
   call deleteDataLayout(t3k)
   call deleteDataLayout(ddx)
   call deleteDataLayout(ddy)
   call deleteDataLayout(ddz)
   call deleteDataLayout(temp)

 end subroutine computeViscousDissResolFiltered

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the 6 components of stress tensor 
!> @param [in] in1K first component in Fourier space
!> @param [in] in2K second component in Fourier space
!> @param [in] in3K third component in Fourier space
!> @param [in] wave the WaveNumbers associated to inK variable
!> @param [in,out] out11 calculated term in physical space 
!> @param [in,out] out22 calculated term in physical space 
!> @param [in,out] out33 calculated term in physical space 
!> @param [in,out] out12 calculated term in physical space 
!> @param [in,out] out13 calculated term in physical space 
!> @param [in,out] out23 calculated term in physical space 
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
 subroutine computeStressTensor(in1K,in2K,in3K,wave,out11,out12,out13, &
            & out22,out23,out33,success)

   implicit none

   !I/O
   type(complex_data_layout),intent(in)  :: in1K,in2K,in3K
   type(real_data_layout),intent(inout)  :: out11,out12,out13
   type(real_data_layout),intent(inout)  :: out22,out23,out33
   type(WaveNumbers),intent(in)          :: wave
   logical, intent(inout)                :: success  

   !Local data
   type(complex_data_layout)    :: dxin1K,dxin2K,dxin3K
   type(complex_data_layout)    :: dyin1K,dyin2K,dyin3K
   type(complex_data_layout)    :: dzin1K,dzin2K,dzin3K
   type(real_data_layout)       :: dyin1,dzin1,dxin2
   type(real_data_layout)       :: dzin2,dxin3,dyin3
   logical                      :: res  


   success = .true.

   if (.not. copyStructOnly(in1K,dxin1K) .or. &
      &.not. copyStructOnly(in1K,dyin1K) .or. &
      &.not. copyStructOnly(in1K,dzin1K) .or. &
      &.not. copyStructOnly(in1K,dxin2K) .or. &
      &.not. copyStructOnly(in1K,dyin2K) .or. &
      &.not. copyStructOnly(in1K,dzin2K) .or. &
      &.not. copyStructOnly(in1K,dxin3K) .or. &
      &.not. copyStructOnly(in1K,dyin3K) .or. &
      &.not. copyStructOnly(in1K,dzin3K) )then 
       success = .false.
       write(6,'(a)')'[ERROR] in computeStressTensor : not enought memory !'
   endif 

   if (.not. copyStructOnly(out11,dyin1) .or. &
      &.not. copyStructOnly(out11,dzin1) .or. &
      &.not. copyStructOnly(out11,dxin2) .or. &
      &.not. copyStructOnly(out11,dzin2) .or. &
      &.not. copyStructOnly(out11,dxin3) .or. &
      &.not. copyStructOnly(out11,dyin3) )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeStressTensor : not enought memory !'
   endif 

   call computeGradientK(wave,in1K,dxin1K,dyin1K,dzin1K)
   call computeGradientK(wave,in2K,dxin2K,dyin2K,dzin2K)
   call computeGradientK(wave,in3K,dxin3K,dyin3K,dzin3K) 

   call btran(dxin1K,out11,res)
   if ( .not. res) success = .false.
   call btran(dyin1K,dyin1,res)
   if ( .not. res) success = .false.
   call btran(dzin1K,dzin1,res)
   if ( .not. res) success = .false.
   call btran(dxin2K,dxin2,res)
   if ( .not. res) success = .false.
   call btran(dyin2K,out22,res)
   if ( .not. res) success = .false.
   call btran(dzin2K,dzin2,res)
   if ( .not. res) success = .false.
   call btran(dxin3K,dxin3,res)
   if ( .not. res) success = .false.
   call btran(dyin3K,dyin3,res)
   if ( .not. res) success = .false.
   call btran(dzin3K,out33,res)
   if ( .not. res) success = .false.

   out12%values = 0.5_WP*(dyin1%values+dxin2%values)
   out13%values = 0.5_WP*(dzin1%values+dxin3%values)
   out23%values = 0.5_WP*(dzin2%values+dyin3%values)

   call deleteDataLayout(dxin1K)
   call deleteDataLayout(dyin1K)
   call deleteDataLayout(dzin1K)
   call deleteDataLayout(dxin2K)
   call deleteDataLayout(dyin2K)
   call deleteDataLayout(dzin2K)
   call deleteDataLayout(dxin3K)
   call deleteDataLayout(dyin3K)
   call deleteDataLayout(dzin3K)
   call deleteDataLayout(dyin1)
   call deleteDataLayout(dzin1)
   call deleteDataLayout(dxin2)
   call deleteDataLayout(dzin2)
   call deleteDataLayout(dxin3)
   call deleteDataLayout(dyin3)

 end subroutine computeStressTensor


!------------------------------------------------------------------------------
!> Compute energy of a scalar field 
!! @author = Antoine Vollant, LEGI
!!    @param[in]    inField1    = first component of an array field
!!    @param[in]    inField2    = second component of an array field
!!    @param[in]    inField3    = third component of an array field
!!    @param[in]    inField3    = third component of an array field
!!    @param[in]    nbcpus      = number of process in the pool 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       energy 
!------------------------------------------------------------------------------
function computeEnergyScal(inField1,spec_rank,nbcpus) result(energy)


    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: inField1
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus
    real(WP)                            :: energy
 
    !Local data
    type(REAL_DATA_LAYOUT)              :: BufField

    if (.not. copyStructOnly(inField1,BufField) ) then
      write(6,'(a,i0,a)')'[ERROR] compute Energy, copy struct : failed!'
    endif

    BufField%values = inField1%values**2.0
    energy =computeFieldAvg( BufField , spec_rank , nbcpus )
    energy = 0.5000000_WP*energy

    call deleteDataLayout(BufField)

end function computeEnergyScal

!------------------------------------------------------------------------------
!> Compute energy of an array field 
!! @author = Antoine Vollant, LEGI
!!    @param[in]    inField1    = first component of an array field
!!    @param[in]    inField2    = second component of an array field
!!    @param[in]    inField3    = third component of an array field
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       energy 
!------------------------------------------------------------------------------
function computeEnergyArray(inField1,inField2,inField3,spec_rank,nbcpus) result(energy)


    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: inField1,inField2,inField3
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus 
    real(WP)                            :: energy
 
    !Local data
    type(REAL_DATA_LAYOUT)              :: BufField

    if (.not. copyStructOnly(inField1,BufField) ) then
      write(6,'(a,i0,a)')'[ERROR] compute Energy, copy struct : failed!'
    endif

    BufField%values = inField1%values**2.0 + inField2%values**2.0 + inField3%values**2.0
    energy =computeFieldAvg( BufField , spec_rank , nbcpus )
    energy = 0.5000000_WP*energy

    call deleteDataLayout(BufField)

end function computeEnergyArray 


!------------------------------------------------------------------------------
!> Compute palinstrophy of velocity 
!! @author = Antoine Vollant, LEGI
!!    @param[in]    VelWN       = velocity waveNumber
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    nbcpus      = number of processus in the pool 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       palinstrophy 
!------------------------------------------------------------------------------
function computePalinstrophy(VelWN,U,V,W,nbcpus,spec_rank) result(palinstrophy)

    use differential_tools
    
    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: U,V,W 
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus 
    type(waveNumbers),intent(in)        :: VelWN
    real(WP)                            :: palinstrophy
 
    !Local data
    type(REAL_DATA_LAYOUT)              :: bufField1,bufField2,bufField3
    type(REAL_DATA_LAYOUT)              :: bufField4,bufField5,bufField6

    if (.not. copyStructOnly(U,bufField1) .or. &
       &.not. copyStructOnly(U,bufField2) .or. &
       &.not. copyStructOnly(U,bufField3) .or. &
       &.not. copyStructOnly(U,bufField4) .or. &
       &.not. copyStructOnly(U,bufField5) .or. &
       &.not. copyStructOnly(U,bufField6) ) then
      write(6,'(a,i0,a)')'[ERROR] Palinstrophy, copyStruct : failed!'
    endif

    if (.not. computeFieldCurl(VelWN,U,V,W,bufField1,bufField2,bufField3,&
       &      nbcpus,spec_rank) ) then
       write(6,'(a,i0,a)')'[ERROR] Palinstrophy : computeFieldCurl failed!'
    endif
    if (.not. computeFieldCurl(VelWN,bufField1,bufField2,bufField3,&
       &      bufField4,bufField5,bufField6,nbcpus,spec_rank) ) then
       write(6,'(a,i0,a)')'[ERROR] Palinstrophy : computeFieldCurl failed!'
    endif

    bufField1%values = bufField4%values**2.0 + &
                     & bufField5%values**2.0 + &
                     & bufField6%values**2.0

    palinstrophy = computeFieldAvg( bufField1 , spec_rank , nbcpus)
    palinstrophy = 0.50000_WP*palinstrophy 

    call deleteDataLayout(bufField1)
    call deleteDataLayout(bufField2)
    call deleteDataLayout(bufField3)
    call deleteDataLayout(bufField4)
    call deleteDataLayout(bufField5)
    call deleteDataLayout(bufField6)

end function computePalinstrophy

!------------------------------------------------------------------------------
!> Compute enstrophy of velocity 
!! @author = Antoine Vollant, LEGI
!!    @param[in]    VelWN       = velocity waveNumber
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    nbcpus      = number of processus in the pool 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       enstrophy 
!------------------------------------------------------------------------------
function computeEnstrophy(VelWN,U,V,W,nbcpus,spec_rank) result(enstrophy)

    use differential_tools
    
    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: U,V,W 
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus 
    type(waveNumbers),intent(in)        :: VelWN
    real(WP)                            :: enstrophy 
 
    !Local data
    type(REAL_DATA_LAYOUT)              :: bufField1,bufField2,bufField3

    if (.not. copyStructOnly(U,bufField1) .or. &
       &.not. copyStructOnly(U,bufField2) .or. &
       &.not. copyStructOnly(U,bufField3) ) then
      write(6,'(a,i0,a)')'[ERROR] Enstrophy , copyStruct : failed!'
    endif

    if (.not. computeFieldCurl(VelWN,U,V,W,bufField1,bufField2,bufField3,&
       &      nbcpus,spec_rank) ) then
       write(6,'(a,i0,a)')'[ERROR] Enstrophy : computeFieldCurl failed!'
    endif

    bufField1%values = bufField1%values**2.0 + &
                     & bufField2%values**2.0 + &
                     & bufField3%values**2.0

    enstrophy = computeFieldAvg( bufField1 , spec_rank , nbcpus )
    enstrophy = 0.50000_WP*enstrophy 

    call deleteDataLayout(bufField1)
    call deleteDataLayout(bufField2)
    call deleteDataLayout(bufField3)

end function computeEnstrophy 


!------------------------------------------------------------------------------
!> Compute helicity of velocity 
!! @author = Antoine Vollant, LEGI
!!    @param[in]    VelWN       = velocity waveNumber
!!    @param[in]    U           = longitudinal velocity
!!    @param[in]    V           = vertical velocity
!!    @param[in]    W           = spanwise velocity
!!    @param[in]    nbcpus      = number of processus in the pool 
!!    @param[in]    spec_rank   = mpi rank inside the spectral communicator
!!    @return       helicity 
!------------------------------------------------------------------------------
function computeHelicity(VelWN,U,V,W,nbcpus,spec_rank) result(helicity)

    use differential_tools   
 
    !I/O data
    type(REAL_DATA_LAYOUT), intent(in)  :: U,V,W 
    integer, intent(in)                 :: spec_rank
    integer, intent(in)                 :: nbcpus 
    type(waveNumbers),intent(in)        :: VelWN
    real(WP)                            :: helicity 
 
    !Local data
    type(REAL_DATA_LAYOUT)              :: bufField1,bufField2,bufField3

    if (.not. copyStructOnly(U,bufField1) .or. &
       &.not. copyStructOnly(U,bufField2) .or. &
       &.not. copyStructOnly(U,bufField3) ) then
      write(6,'(a,i0,a)')'[ERROR] computeHelicity , copyStruct : failed!'
    endif

    if (.not. computeFieldCurl(VelWN,U,V,W,bufField1,bufField2,bufField3,&
       &      nbcpus,spec_rank) ) then
       write(6,'(a,i0,a)')'[ERROR] computeHelicity : computeFieldCurl failed!'
    endif

    bufField1%values = U%values*bufField1%values + &
                     & V%values*bufField2%values + &
                     & W%values*bufField3%values

    helicity = computeFieldAvg( bufField1 , spec_rank , nbcpus )
    helicity = 0.50000_WP*helicity 

    call deleteDataLayout(bufField1)
    call deleteDataLayout(bufField2)
    call deleteDataLayout(bufField3)

end function computeHelicity


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!> Compute the norm |sqrt(Sij Sij)| 
!> @param [in] in11 calculated term in physical space 
!> @param [in] in22 calculated term in physical space 
!> @param [in] in33 calculated term in physical space 
!> @param [in] in12 calculated term in physical space 
!> @param [in] in13 calculated term in physical space 
!> @param [in] in23 calculated term in physical space 
!> @param [out] norm scalar field which contains the norm
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize
!------------------------------------------------------------------------------
 subroutine computeTensorNorme1(in11,in12,in13,in22,in23,in33,norm,res)

   !I/O
   type(real_data_layout),intent(in)     :: in11,in12,in13
   type(real_data_layout),intent(in)     :: in22,in23,in33
   type(real_data_layout),intent(inout)  :: norm 
   logical, intent(inout)                :: res  
   
   res = .false.

   norm%values = sqrt(2.0_WP *( in11%values**2 + in22%values**2 + &
               & in33%values**2 + 2.0_WP*in12%values**2 + &
               & 2.0_WP*in13%values**2 + 2.0_WP*in23%values**2) )

   res = .true.

 end subroutine computeTensorNorme1

!------------------------------------------------------------------------------
!> @author 
!> Guillaume Balarac, LEGI
!>
!> @details
!> Compute the 3 components of Vorticty Vector 
!> @param [in] in1K first component in Fourier space
!> @param [in] in2K second component in Fourier space
!> @param [in] in3K third component in Fourier space
!> @param [in] wave the WaveNumbers associated to inK variable
!> @param [in,out] out1 calculated term in physical space -> W1 ou O23 ou -O32 
!> @param [in,out] out2 calculated term in physical space -> W2 ou O31 ou -O13
!> @param [in,out] out3 calculated term in physical space -> W3 ou O12 ou -O21
!> @param [out] success .TRUE. if all is succcessfull, .FALSE. otherwize.
!------------------------------------------------------------------------------
 subroutine computeVorticity(in1K,in2K,in3K,wave,out1,out2,out3,success)

   implicit none

   !I/O
   type(complex_data_layout),intent(in)  :: in1K,in2K,in3K
   type(real_data_layout),intent(inout)     :: out1,out2,out3
   type(WaveNumbers),intent(in)          :: wave
   logical, intent(inout)                :: success  

   !Local data
   type(complex_data_layout)    :: dxin1K,dxin2K,dxin3K
   type(complex_data_layout)    :: dyin1K,dyin2K,dyin3K
   type(complex_data_layout)    :: dzin1K,dzin2K,dzin3K
   type(real_data_layout)       :: dyin1,dzin1,dxin2
   type(real_data_layout)       :: dzin2,dxin3,dyin3
   logical                      :: res  


   success = .true.

   if (.not. copyStructOnly(in1K,dxin1K) .or. &
      &.not. copyStructOnly(in1K,dyin1K) .or. &
      &.not. copyStructOnly(in1K,dzin1K) .or. &
      &.not. copyStructOnly(in1K,dxin2K) .or. &
      &.not. copyStructOnly(in1K,dyin2K) .or. &
      &.not. copyStructOnly(in1K,dzin2K) .or. &
      &.not. copyStructOnly(in1K,dxin3K) .or. &
      &.not. copyStructOnly(in1K,dyin3K) .or. &
      &.not. copyStructOnly(in1K,dzin3K) )then 
       success = .false.
       write(6,'(a)')'[ERROR] in computeVorticty : not enought memory !'
   endif 

   if (.not. copyStructOnly(out1,dyin1) .or. &
      &.not. copyStructOnly(out1,dzin1) .or. &
      &.not. copyStructOnly(out1,dxin2) .or. &
      &.not. copyStructOnly(out1,dzin2) .or. &
      &.not. copyStructOnly(out1,dxin3) .or. &
      &.not. copyStructOnly(out1,dyin3) )then
       success = .false.
       write(6,'(a)')'[ERROR] in computeVorticity : not enought memory !'
   endif 

   call computeGradientK(wave,in1K,dxin1K,dyin1K,dzin1K)
   call computeGradientK(wave,in2K,dxin2K,dyin2K,dzin2K)
   call computeGradientK(wave,in3K,dxin3K,dyin3K,dzin3K) 

   call btran(dyin1K,dyin1,res)
   if ( .not. res) success = .false.
   call btran(dzin1K,dzin1,res)
   if ( .not. res) success = .false.
   call btran(dxin2K,dxin2,res)
   if ( .not. res) success = .false.
   call btran(dzin2K,dzin2,res)
   if ( .not. res) success = .false.
   call btran(dxin3K,dxin3,res)
   if ( .not. res) success = .false.
   call btran(dyin3K,dyin3,res)
   if ( .not. res) success = .false.

   out1%values = 0.5_WP*(dyin3%values-dzin2%values) !W1 ou O23 ou -O32
   out2%values = 0.5_WP*(dzin1%values-dxin3%values) !W2 ou O31 ou -O13
   out3%values = 0.5_WP*(dxin2%values-dyin1%values) !W3 ou O12 ou -O21

   call deleteDataLayout(dxin1K)
   call deleteDataLayout(dyin1K)
   call deleteDataLayout(dzin1K)
   call deleteDataLayout(dxin2K)
   call deleteDataLayout(dyin2K)
   call deleteDataLayout(dzin2K)
   call deleteDataLayout(dxin3K)
   call deleteDataLayout(dyin3K)
   call deleteDataLayout(dzin3K)
   call deleteDataLayout(dyin1)
   call deleteDataLayout(dzin1)
   call deleteDataLayout(dxin2)
   call deleteDataLayout(dzin2)
   call deleteDataLayout(dxin3)
   call deleteDataLayout(dyin3)

 end subroutine computeVorticity


!------------------------------------------------------------------------------
!> Compute Q criterion 
!! @author Guillaume Balarac 
!!    @param[in]     Uk
!!    @param[in]     Vk
!!    @param[in]     Wk
!!    @param[in]     VelWN
!!    @param[in,out] Qc
!!    @param[in,out] spec_rank
!!    @return        success         = logical equals to true is everything is right.
!------------------------------------------------------------------------------
function computeParaQc(U,Uk,Vk,Wk,VelWN,Qc,spec_rank) result(success)

    ! Input/Output
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk
    type(REAL_DATA_LAYOUT), intent(in)    :: U
    type(REAL_DATA_LAYOUT), intent(inout) :: Qc
    TYPE(WaveNumbers),INTENT(IN)          :: velWN
    integer, intent(in)                   :: spec_rank
    logical                               :: success

    ! Local varaibales
    type(REAL_DATA_LAYOUT)                :: S11,S12,S13,S22,S23,S33
    type(REAL_DATA_LAYOUT)                :: O1,O2,O3

    success = .true.
    if (.not. copyStructOnly(U,S11) .or. &
       &.not. copyStructOnly(U,S12) .or. &
       &.not. copyStructOnly(U,S13) .or. &
       &.not. copyStructOnly(U,S22) .or. &
       &.not. copyStructOnly(U,S23) .or. &
       &.not. copyStructOnly(U,S33) .or. &
       &.not. copyStructOnly(U,O1)  .or. &
       &.not. copyStructOnly(U,O2)  .or. &
       &.not. copyStructOnly(U,O3) )then
         success = .false.
         write(6,'(a,i0)')'[ERROR] in computeParaQc : not enought memory !',spec_rank
    endif
    ! Compute Sij
    call computeStressTensor(Uk,VK,WK,VelWN,S11,S12,S13,S22,S23,S33,success)
    ! Compute Oi
    call computeVorticity(UK,VK,WK,VelWN,O1,O2,O3,success)
    !compute Q criterion
    Qc%values = 0.5_WP*( S11%values**2.0 + S22%values**2.0 + S33%values**2.0 + &
          &        2.0_WP * (S12%values**2.0 + S13%values**2.0 + S23%values**2.0) - &
          &        2.0_WP * (O1%values**2.0 + O2%values**2.0 + O3%values**2.0 ) )
    Qc%name = 'Qcrit'
    call deleteDataLayout(S11)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S22)
    call deleteDataLayout(S23)
    call deleteDataLayout(S33)
    call deleteDataLayout(O1)
    call deleteDataLayout(O2)
    call deleteDataLayout(O3)

end function computeParaQc 

!------------------------------------------------------------------------------
!! @details
!> Compute angle between two vectors from arcsin [0;Pi/2] as arcsin is computed from |sin| 
!! @author Antoine Vollant 
!!    @param[in]     vecX1
!!    @param[in]     vecY1
!!    @param[in]     vecZ1
!!    @param[in]     vecX2
!!    @param[in]     vecY2
!!    @param[in]     vecZ2
!!    @return        arcsinus which is the angle in radian 
!------------------------------------------------------------------------------
function computeArcSinVectors(vecX1,vecY1,vecZ1,vecX2,vecY2,vecZ2) result(arcsinus)

  use precision_tools

  !I/O data
  real(WP),intent(in) :: vecX1,vecY1,vecZ1,vecX2,vecY2,vecZ2
  real(WP)            :: arcsinus

  !Local variables
  real(WP)            :: norm1,norm2,normProd
  real(WP)            :: prodVX,prodVY,prodVZ

  prodVX = vecY1*vecZ2-vecZ1*vecY2
  prodVY = vecZ1*vecX2-vecX1*vecZ2
  prodVZ = vecX1*vecY2-vecY1*vecX2
  normProd = sqrt(prodVX**2+prodVY**2+prodVZ**2)
  norm1 = sqrt(vecX1**2+vecY1**2+vecZ1**2)
  norm2 = sqrt(vecX2**2+vecY2**2+vecZ2**2)
  arcsinus = asin(normProd /( norm1 * norm2 ))

end function computeArcSinVectors


!------------------------------------------------------------------------------
!! @details
!> Compute angle between two vectors from arccosinus [0;Pi]
!! @author Antoine Vollant 
!!    @param[in]     vecX1
!!    @param[in]     vecY1
!!    @param[in]     vecZ1
!!    @param[in]     vecX2
!!    @param[in]     vecY2
!!    @param[in]     vecZ2
!!    @return        arccosinus 
!------------------------------------------------------------------------------
function computeArcCosVectors(vecX1,vecY1,vecZ1,vecX2,vecY2,vecZ2) result(arccosinus)

  use precision_tools

  !I/O data
  real(WP),intent(in) :: vecX1,vecY1,vecZ1,vecX2,vecY2,vecZ2
  real(WP)            :: arccosinus 

  !Local variables
  real(WP)            :: norm1,norm2,prodScal

  norm1 = sqrt(vecX1**2+vecY1**2+vecZ1**2)
  norm2 = sqrt(vecX2**2+vecY2**2+vecZ2**2)
  prodScal = vecX1*vecX2 + vecY1*vecY2 + vecZ1*vecZ2
  arccosinus = acos ( prodScal /( norm1 * norm2 ) )

end function

 !-------------------------------------------------------------------------------------------
 !! @details
 !> Compute the angle between two vectors and the ratio between their norms 
 !!
 !! @author Antoine Volllant 
 !!    @param[in]     vecRef
 !!    @param[in]     vecTest 
 !!    @param[inout]  norm datalayout which contains the ratio between norms
 !!    @param[inout]  angle datalyout which contains the angles between vectors
 !!    @param[in]     dumpFields logical value
 !!    @return        success = logical value 
 !-------------------------------------------------------------------------------------------
 function computeVectorSimilarity(vecRef,vecTest,norm,angle,dumpFields,ific) result (success)

   use avs , only : dump_UVWP_ForAVS

   !I/O data
   type(real_data_layout),dimension(:),intent(in) :: vecRef,vecTest 
   type(real_data_layout),intent(inout)        :: norm,angle 
   logical,intent(in)                             :: dumpFields 
   integer,intent(in),optional                    :: ific
   logical                                        :: success

   !Local data
   integer                                        :: i,j,k
   integer                                        :: spec_rank
   character(len=8)                               :: writecom

   success = .false.
   if ( (size(vecRef,1) .ne. 3) .or. (size(vecTest,1) .ne. 3) ) then
     write (6,'(a)') '[ERROR] in computeVectorSimilarity : parameters have wrong size'
     return
   endif
   if (.not. sameLayout(vecRef(1),vecTest(1),norm,angle) ) return
   do k=vecRef(1)%zmin,vecRef(1)%zmax
     do j=vecRef(1)%ymin,vecRef(1)%ymax
       do i=vecRef(1)%xmin,vecRef(1)%xmax
         angle%values(i,j,k) =  computeArcCosVectors( vecRef(1)%values(i,j,k) ,vecRef(2)%values(i,j,k), &
                                                    & vecRef(3)%values(i,j,k) ,vecTest(1)%values(i,j,k),&
                                                    & vecTest(2)%values(i,j,k),vecTest(3)%values(i,j,k) )
         angle%values(i,j,k) = angle%values(i,j,k) * 180.0_WP/ acos(-1.0_WP)
         norm%values(i,j,k)  = sqrt ( vecTest(1)%values(i,j,k)**2 + vecTest(2)%values(i,j,k)**2 + vecTest(3)%values(i,j,k)**2 ) 
         norm%values(i,j,k)  = norm%values(i,j,k) / sqrt ( vecRef(1)%values(i,j,k)**2 + vecRef(2)%values(i,j,k)**2 + vecRef(3)%values(i,j,k)**2 )
       enddo
     enddo
   enddo

   if (dumpFields) then
      spec_rank = getnbmycpu()
      if (spec_rank .eq. 0) write (6,'(a)') '[INFO] dumping angle and ratio norm for TiMod VS TiDNS in avs fields'
      if ( .not. present(ific)) return
      write (writecom,'(i6.6)') ific
      angle%name = 'angle'
      norm%name = 'normRatio' 
      if (.not. dump_UVWP_ForAVS(trim(angle%name)//'_'//trim(writecom),spec_rank,angle%lx,angle%ly,angle%lz,angle,usename=.true.)) return
      if (.not. dump_UVWP_ForAVS(trim(norm%name)//'_'//trim(writecom),spec_rank,norm%lx,norm%ly,norm%lz,norm,usename=.true.)) return
   endif
   success = .true.

 end function


!------------------------------------------------------------------------------
!> computes bernstein's polynome 
!!    @param[in,out]    m = number of points minus one in the caracteritic polygone 
!!    @param[in]        i = current indice in cumputation
!!    @param[in]        param = center coordinate (optional, if not precise, then the domain center)
!!    @return           val 
!------------------------------------------------------------------------------
function computeBernsteinPol( m , i , param ) result(val)
 
    !I/O data
    real(WP),intent(in)                   :: param
    integer,intent(in)                    :: m 
    integer,intent(in)                    :: i
    real(WP)                              :: val

    if ( param .gt. 1.0_WP ) write (6,'(a)') '[ERROR] param .gt. 1.0_WP in computeBernsteinPol'

    val = factoriel(m) / (factoriel(i) * factoriel(m-i))
    val = val * param ** i
    val = val * (1.0_WP-param)**(m-i)

end function computeBernsteinPol


!------------------------------------------------------------------------------
!> computes factoriel of a integer 
!!    @param[in]        val to compute factoriel 
!!    @return           factoriel val!
!------------------------------------------------------------------------------
function factoriel(val) result(fac)

    !I/O data
    integer,intent(in) :: val
    integer            :: fac

    !Local data
    integer            :: i

    fac=1
    if ( val .gt. 0) then
      fac = 1
      do i=1,val
        fac = fac * i
      enddo
    endif

end function factoriel


!------------------------------------------------------------------------------
!> computes levi-Civita indice
!!    @param[in]        i first indice
!!    @param[in]        j second indice
!!    @param[in]        k third indice
!!    @return           res equal to 1 , -1 or 0 according to (i,j,k) in [1,2,3]
!------------------------------------------------------------------------------
function leviCivita(i,j,k) result(res)

  !I/O data
  integer,intent(in) :: i,j,k
  integer            :: res

  if     ( i == 1 .and. j == 2 .and. k == 3 .or. &
           i == 2 .and. j == 3 .and. k == 1 .or. &
           i == 3 .and. j == 1 .and. k == 2 ) then
      res = 1
  elseif ( i == 1 .and. j == 3 .and. k == 2 .or. &
           i == 2 .and. j == 1 .and. k == 3 .or. &
           i == 3 .and. j == 2 .and. k == 1 ) then
      res = -1
  else
      res = 0
  endif

end function leviCivita

end module physical_values_tools
!> @}
