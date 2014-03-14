module Eqdm_lorentz_models

    USE wavenumber_tools
    USE transforms_tools
    use toolbox
    use post_lib
    use subgrid_tools
    use datalayout
    use differential_tools 

implicit none




 contains
! 
! subroutine InducMod(Bin,nlB,WaveB,InduCcase,delta)
! 
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN), DIMENSION(:) :: B 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT), DIMENSION(:) :: nlB
!   INTEGER, INTENT(IN)  :: spec_rank,InduCcase
!   
! 
!   call InducSmago(Bin,WaveB,delta,T1Ind,T2Ind,T3Ind)
!   
!   nlB(1)%values = nlB(1)%values - T1Ind%values
!   nlB(2)%values = nlB(2)%values - T2Ind%values
!   nlB(3)%values = nlB(3)%values - T3Ind%values
! 
! end subroutine InducMod

! subroutine Velmod(U,V,W,Uk,Vk,Wk,Bin,Bk,nlUx,nlUy,nlUz,WaveB,delta)
! 
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN), DIMENSION(:) :: B 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT), DIMENSION(:) :: nlB
!   INTEGER, INTENT(IN)  :: spec_rank
!   
!   call Velunik()
! end subroutine Velmod
! 
! subroutine Velbb(U,V,W,Uk,Vk,Wk,Bin,Bk,WaveB,delta,T1,T2,T3)
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN), DIMENSION(:) :: B 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN), DIMENSION(:) :: Bk
!   TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
!   TYPE(REAL_DATA_LAYOUT), DIMENSION(:) :: Bf 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3 
!   INTEGER, INTENT(IN)  :: spec_rank
!   TYPE(REAL_DATA_LAYOUT) :: t12,t13,t21,t23,t31,t32
!   logical :: res
!   REAL(WP) :: SmagCst
!   
! end subroutine Velbb



 subroutine magneticEqdmHamba(Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  deltx,T11,T22,T33,T12,T13,T23)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
!   type(REAL_DATA_LAYOUT), intent(inout)          :: rot1,rot2,rot3  
  real(WP), intent(in)                        :: deltx,C_Lambda,C_mu
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: J23,J12,J13  
  type(real_data_layout)          :: Jdens,normS

  type(real_data_layout)          :: S11,S22,S33 
  type(real_data_layout)          :: S23,S12,S13  
 
  type(real_data_layout)          :: T11,T22,T33 
  type(real_data_layout)          :: T23,T12,T13  
 

 res = .true.

  if (.not. copyStructOnly(T11,J12) .or. &
     &.not. copyStructOnly(T11,J13) .or. &
     &.not. copyStructOnly(T11,J23) .or. &
     &.not. copyStructOnly(T11,Jdens) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(T11,S11) .or. &
     &.not. copyStructOnly(T11,S22) .or. &
     &.not. copyStructOnly(T11,S33) .or. &
     &.not. copyStructOnly(T11,S12) .or. &
     &.not. copyStructOnly(T11,S13) .or. &
     &.not. copyStructOnly(T11,S23)  )then

       res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
  call computeTensorJdens1(J12,J13,J23,Jdens,res)

  call computeStressTensor(Uk,VK,WK,wave,S11,S12,S13,S22,S23,S33,res)
!   call computeTensorNorme1(S11,S12,S13,S22,S23,S33,normS,res)

  T11%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S11%values
  T22%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S22%values
  T33%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S33%values
  T12%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S12%values
  T13%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S13%values
  T23%values= -C_mu*deltx**2 *sqrt(C_lambda)*Jdens%values*S23%values

    call deleteDataLayout(Jdens)
    call deleteDataLayout(J12)
    call deleteDataLayout(J13)
    call deleteDataLayout(J23)

!     call deleteDataLayout(normS)
    call deleteDataLayout(S11)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)

 end subroutine magneticEqdmHamba


 subroutine HambaYoshizawaNSequation(Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  deltx,T11,T22,T33,T12,T13,T23)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
!   type(REAL_DATA_LAYOUT), intent(inout)          :: rot1,rot2,rot3  
  real(WP), intent(in)                        :: deltx,C_Lambda,C_mu
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: J23,J12,J13  
  type(real_data_layout)          :: Jdens,normS

  type(real_data_layout)          :: S11,S22,S33 
  type(real_data_layout)          :: S23,S12,S13  
 
  type(real_data_layout)          :: T11,T22,T33 
  type(real_data_layout)          :: T23,T12,T13  
 

 res = .true.

  if (.not. copyStructOnly(T11,J12) .or. &
     &.not. copyStructOnly(T11,J13) .or. &
     &.not. copyStructOnly(T11,J23) .or. &
     &.not. copyStructOnly(T11,Jdens) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(T11,S11) .or. &
     &.not. copyStructOnly(T11,S22) .or. &
     &.not. copyStructOnly(T11,S33) .or. &
     &.not. copyStructOnly(T11,S12) .or. &
     &.not. copyStructOnly(T11,S13) .or. &
     &.not. copyStructOnly(T11,S23) .or. &
     &.not. copyStructOnly(T11,normS) )then

       res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
  call computeTensorJdens1(J12,J13,J23,Jdens,res)

  call computeStressTensor(Uk,VK,WK,wave,S11,S12,S13,S22,S23,S33,res)
  call computeTensorNorme1(S11,S12,S13,S22,S23,S33,normS,res)

  T11%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S11%values
  T22%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S22%values
  T33%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S33%values
  T12%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S12%values
  T13%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S13%values
  T23%values= -C_mu*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*S23%values

    call deleteDataLayout(Jdens)
    call deleteDataLayout(J12)
    call deleteDataLayout(J13)
    call deleteDataLayout(J23)

    call deleteDataLayout(normS)
    call deleteDataLayout(S11)
    call deleteDataLayout(S22)
    call deleteDataLayout(S33)
    call deleteDataLayout(S12)
    call deleteDataLayout(S13)
    call deleteDataLayout(S23)

 end subroutine HambaYoshizawaNSequation


subroutine DivHambaYoshizawaNSequation(Uin,Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  deltx,T1,T2,T3)

implicit none

  TYPE(Complex_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk,Bxk,Byk,Bzk
  INTEGER :: spec_rank
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  REAL(WP) :: deltx,C_Lambda,C_mu
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(REAL_DATA_LAYOUT) :: t11,t12,t13,t22,t23,t33,Uin


  IF ((.NOT. copyStructOnly(Uin,t11)) .OR. &
     &(.NOT. copyStructOnly(Uin,t12)) .OR. &
     &(.NOT. copyStructOnly(Uin,t13)) .OR. &
     &(.NOT. copyStructOnly(Uin,t22)) .OR. &
     &(.NOT. copyStructOnly(Uin,t23)) .OR. &
     &(.NOT. copyStructOnly(Uin,t33))    ) RETURN 

  call HambaYoshizawaNSequation(Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  deltx,T11,T22,T33,T12,T13,T23)

  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,wave,res)

    CALL deleteDataLayout(t11)
    CALL deleteDataLayout(t12)
    CALL deleteDataLayout(t13)
    CALL deleteDataLayout(t22)
    CALL deleteDataLayout(t23)
    CALL deleteDataLayout(t33)

end subroutine DivHambaYoshizawaNSequation

 subroutine MullerCarratiNSequation(Uk,Vk,Wk,Bxk,Byk,Bzk,wave,spec_rank,res, &
                          &  deltx,T23,T12,T13,T11,T22,T33)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
!   type(REAL_DATA_LAYOUT)   :: rot1,rot2,rot3  
  type(REAL_DATA_LAYOUT), intent(inout)   :: T23,T12,T13,T11,T22,T33  
  


!   type(REAL_DATA_LAYOUT)                  :: Bx  

  real(WP), intent(in)                    :: deltx
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: normSuSb

  type(real_data_layout)          :: SU11,SU22,SU33 
  type(real_data_layout)          :: SU23,SU12,SU13  
 
  type(real_data_layout)          :: SB11,SB22,SB33 
  type(real_data_layout)          :: SB23,SB12,SB13  
 

 res = .true.

  if (.not. copyStructOnly(T12,SU12) .or. &
     &.not. copyStructOnly(T12,SU13) .or. &
     &.not. copyStructOnly(T12,SU23) .or. &
     &.not. copyStructOnly(T12,SU11) .or. &
     &.not. copyStructOnly(T12,SU22) .or. &
     &.not. copyStructOnly(T12,SU33)  )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiNS : not enought memory !',spec_rank
   endif


  if (.not. copyStructOnly(T12,SB12) .or. &
     &.not. copyStructOnly(T12,SB13) .or. &
     &.not. copyStructOnly(T12,SB23) .or. &
     &.not. copyStructOnly(T12,SB11) .or. &
     &.not. copyStructOnly(T12,normSuSb) .or. & 
     &.not. copyStructOnly(T12,SB22) .or. &
     &.not. copyStructOnly(T12,SB33)  )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiNS : not enought memory !',spec_rank
   endif


  call computeStressTensor(Uk,VK,WK,wave,SU11,SU12,SU13,SU22,SU23,SU33,res)
  call computeStressTensor(Bxk,ByK,BzK,wave,SB11,SB12,SB13,SB22,SB23,SB33,res)

  normSuSb%values=sqrt(abs(SU11%values*SB11%values + SU22%values*SB22%values + SU33%values*SB33%values &
                         & +2.0*(SU12%values*SB12%values + SU13%values*SB13%values + SU23%values*SB23%values) ))

  T12%values= -0.015*deltx**2 *normSuSb%values*SU12%values 
  T13%values= -0.015*deltx**2 *normSuSb%values*SU13%values 
  T23%values= -0.015*deltx**2 *normSuSb%values*SU23%values 

  T11%values= -0.015*deltx**2 *normSuSb%values*SU11%values 
  T22%values= -0.015*deltx**2 *normSuSb%values*SU22%values 
  T33%values= -0.015*deltx**2 *normSuSb%values*SU33%values 

    call deleteDataLayout(SB11)
    call deleteDataLayout(SB22)
    call deleteDataLayout(SB33)
    call deleteDataLayout(SB12)
    call deleteDataLayout(SB13)
    call deleteDataLayout(SB23)

    call deleteDataLayout(SU11)
    call deleteDataLayout(SU22)
    call deleteDataLayout(SU33)
    call deleteDataLayout(SU12)
    call deleteDataLayout(SU13)
    call deleteDataLayout(SU23)

    call deleteDataLayout(normSuSb)

 end subroutine MullerCarratiNSequation


subroutine DivMullerCarratiNSequation(Uin,Uk,Vk,Wk,Bxk,Byk,Bzk,wave,spec_rank,res, &
                          &  deltx,T1,T2,T3)

implicit none

  TYPE(Complex_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk,Bxk,Byk,Bzk
  INTEGER :: spec_rank
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: T1,T2,T3
  REAL(WP) :: deltx,C_Lambda,C_mu
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(REAL_DATA_LAYOUT) :: t11,t12,t13,t22,t23,t33,Uin


  IF ((.NOT. copyStructOnly(Uin,t11)) .OR. &
     &(.NOT. copyStructOnly(Uin,t12)) .OR. &
     &(.NOT. copyStructOnly(Uin,t13)) .OR. &
     &(.NOT. copyStructOnly(Uin,t22)) .OR. &
     &(.NOT. copyStructOnly(Uin,t23)) .OR. &
     &(.NOT. copyStructOnly(Uin,t33))    ) RETURN 

  call MullerCarratiNSequation(Uk,Vk,Wk,Bxk,Byk,Bzk,wave,spec_rank,res, &
                          &  deltx,T23,T12,T13,T11,T22,T33)

  call ComputeDivSymTijVec1Vec2(T1,T2,T3,T12,T13,T23,T11,T22,T33,wave,res)

    CALL deleteDataLayout(t11)
    CALL deleteDataLayout(t12)
    CALL deleteDataLayout(t13)
    CALL deleteDataLayout(t22)
    CALL deleteDataLayout(t23)
    CALL deleteDataLayout(t33)

end subroutine DivMullerCarratiNSequation
end module Eqdm_lorentz_models

