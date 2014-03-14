module mhd_model

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


 subroutine DynsmagorinskyMHDTij(T12,T13,T23,U,V,W,Uk,Vk,Wk,Bxk,Byk,Bzk,Bx,By,Bz,wave,spec_rank,res, &
                          &  numVelFilter,deltx)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
!   type(complex_data_layout),intent(inout) :: div1,div2,div3
  type(real_data_layout),intent(in)       :: U,V,W
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
  integer, intent(in),optional            :: numVelFilter
  type(real_data_layout)          :: T23,T12,T13  

  !Local data
  type(real_data_layout)          :: J23,J12,J13  
  type(real_data_layout)          :: J23f,J12f,J13f  
  type(real_data_layout)          :: Uf, Vf, Wf
  type(real_data_layout)          :: normJ12f,normJ13f,normJ23f  
  type(real_data_layout)          :: coefff
  type(complex_data_layout),dimension(:),allocatable :: tempK, tempKf
  type(real_data_layout)          :: Jdens, Jdensf
  real(WP)                        :: deltx, num, den
  TYPE(REAL_DATA_LAYOUT) :: Wtab1,Wtab2
  TYPE(REAL_DATA_LAYOUT) :: L12,L13,L23

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bx,By,Bz 
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk

  TYPE(REAL_DATA_LAYOUT) :: Bxf,Byf,Bzf 
  TYPE(COMPLEX_DATA_LAYOUT) :: Ukf,Vkf,Wkf
  TYPE(COMPLEX_DATA_LAYOUT) :: Bxkf,Bykf,Bzkf

 res = .true.

! call cpu_time(tps1)
  if (.not. copyStructOnly(U,J12) .or. &
     &.not. copyStructOnly(U,J13) .or. &
     &.not. copyStructOnly(U,J23) .or. &
     &.not. copyStructOnly(U,coefff) .or. &
     &.not. copyStructOnly(U,Jdens) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in smagorinsky for vel : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
  call computeTensorJdens1(J12,J13,J23,Jdens,res)


     ! Dynamic model --- Formulation
     ! Leonard Terms computation -> (u_i u_j)f - uf_i uf_j

  IF ((.NOT. copyStructOnly(Bxk,Ukf)) .OR. &
     &(.NOT. copyStructOnly(Bxk,Vkf)) .OR. &
     &(.NOT. copyStructOnly(Bxk,Wkf)) .OR. &
     &(.NOT. copyStructOnly(Bxk,Bxkf)) .OR. &
     &(.NOT. copyStructOnly(Bxk,Bykf)) .OR. &
     &(.NOT. copyStructOnly(Bxk,Bzkf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Bxf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Byf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Bzf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Uf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Vf)) .OR. &
     &(.NOT. copyStructOnly(Bx,Wf)) .OR. &
     &(.NOT. copyStructOnly(Bx,L12)) .OR. &
     &(.NOT. copyStructOnly(Bx,L13)) .OR. &
     &(.NOT. copyStructOnly(Bx,L23)) .OR. &
     &(.NOT. copyStructOnly(Bx,Wtab1)) .OR. &
     &(.NOT. copyStructOnly(Bx,Wtab2)) ) then

       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for vel : not enought memory !',spec_rank
     endif
     res = initWorkArray(Uk,3,tempK)
     res = initWorkArray(Uk,3,tempKf)

     ! uf_i
     call computeFilter(wave,2.0*deltx,Uk,Ukf,numVelFilter) 
     call computeFilter(wave,2.0*deltx,Vk,Vkf,numVelFilter) 
     call computeFilter(wave,2.0*deltx,Wk,Wkf,numVelFilter) 
     call btran(Ukf,Uf,res)
     call btran(Vkf,Vf,res)
     call btran(Wkf,Wf,res)

     call computeFilter(wave,2.0*deltx,Bxk,Bxkf,numVelFilter) 
     call computeFilter(wave,2.0*deltx,Byk,Bykf,numVelFilter) 
     call computeFilter(wave,2.0*deltx,Bzk,Bzkf,numVelFilter)
     call btran(Bxkf,Bxf,res)
     call btran(Bykf,Byf,res)
     call btran(Bzkf,Bzf,res)

     CALL computeLijMHD(U,V,W,Uf,Vf,Wf,Bx,By,Bz,Bxf,Byf,Bzf&
                     &,Bxk,L12,L13,L23,deltx,wave,numVelFilter)

     ! Compute M_ij
     ! ( norme * Sij )f
     if (.not. copyStructOnly(U,normJ12f) .or. &
        &.not. copyStructOnly(U,normJ13f) .or. &
        &.not. copyStructOnly(U,normJ23f)  ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for Induc 2 : not enought memory !',spec_rank
     endif
     normJ12f%values = Jdens%values * J12%values
     normJ13f%values = Jdens%values * J13%values
     normJ23f%values = Jdens%values * J23%values

     call ftran(normJ12f,tempK(1),res)
     call ftran(normJ13f,tempK(2),res)
     call ftran(normJ23f,tempK(3),res)

     call computeFilter(wave,2.0*deltx,tempK(1),tempKf(1),numVelFilter)
     call computeFilter(wave,2.0*deltx,tempK(2),tempKf(2),numVelFilter) 
     call computeFilter(wave,2.0*deltx,tempK(3),tempKf(3),numVelFilter) 
     call btran(tempKf(1),normJ12f,res) 
     call btran(tempKf(2),normJ13f,res)
     call btran(tempKf(3),normJ23f,res) 

     ! normef * S11f
     if (.not. copyStructOnly(U,J12f) .or. &
        &.not. copyStructOnly(U,J13f) .or. &
        &.not. copyStructOnly(U,J23f) .or. &
        &.not. copyStructOnly(U,Jdensf) ) then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in dyn smagorinsky for induc 2 : not enought memory !',spec_rank
     endif

     call computeJijTensor(BxKf,ByKf,BzKf,wave,J12f,J13f,J23f, res)
     call computeTensorJdens1(J12f,J13f,J23f,Jdensf,res)

     normJ12f%values = 2.0*(2.0 * deltx)**2.0 * Jdensf%values * J12f%values - 2.0*( deltx )**2.0 * normJ12f%values ! = M12
     normJ13f%values = 2.0*(2.0 * deltx)**2.0 * Jdensf%values * J13f%values - 2.0*( deltx )**2.0 * normJ13f%values ! = M13
     normJ23f%values = 2.0*(2.0 * deltx)**2.0 * Jdensf%values * J23f%values - 2.0*( deltx )**2.0 * normJ23f%values ! = M23

     ! L_ij * M_ij
     Wtab1%values  =    (L12%values*normJ12f%values + L13%values*normj13f%values + L23%values*normJ23f%values )
     Wtab2%values  =    (normJ12f%values**2.0 + normJ13f%values**2.0 + normJ23f%values**2.0 )
     
        num = computeFieldAvg( Wtab1 , spec_rank )
        den = computeFieldAvg( Wtab2 , spec_rank )
        coefff%values = num / den
        if (spec_rank.eq.0) print*,'DynCoef induc : ',num / den

     
    call deleteDataLayout(Uf)
    call deleteDataLayout(Vf)
    call deleteDataLayout(Wf)
    call deleteDataLayout(L12)
    call deleteDataLayout(L13)
    call deleteDataLayout(L23)
    call deleteDataLayout(normJ12f)
    call deleteDataLayout(normJ13f)
    call deleteDataLayout(normJ23f)
    call deleteDataLayout(J12f)
    call deleteDataLayout(J13f)
    call deleteDataLayout(J23f)
    call deleteDataLayout(Jdensf)


    CALL deleteDataLayout(Ukf)
    CALL deleteDataLayout(Vkf)
    CALL deleteDataLayout(Wkf)

    CALL deleteDataLayout(Bxkf)
    CALL deleteDataLayout(Bykf)
    CALL deleteDataLayout(Bzkf)
    CALL deleteDataLayout(Bxf)
    CALL deleteDataLayout(Byf)
    CALL deleteDataLayout(Bzf)
    CALL deleteDataLayout(Wtab1)
    CALL deleteDataLayout(Wtab2)

    res = deleteWorkArray(tempK)
    res = deleteWorkArray(tempKf)

  
  T12%values = 2.0*coefff%values * (deltx)**2.0 * Jdens%values * J12%values
  T13%values = 2.0*coefff%values * (deltx)**2.0 * Jdens%values * J13%values
  T23%values = 2.0*coefff%values * (deltx)**2.0 * Jdens%values * J23%values


    call deleteDataLayout(coefff)
    call deleteDataLayout(Jdens)
    call deleteDataLayout(J12)
    call deleteDataLayout(J13)
    call deleteDataLayout(J23)

 end subroutine DynsmagorinskyMHDTij

 subroutine smagorinskyMHDTij(T12,T13,T23,Bxk,Byk,Bzk,wave,spec_rank,res, &
                          &  deltx)


  !I/O
!   type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
  type(real_data_layout), intent(inout)          :: T23,T12,T13  
  real(WP), intent(in)                        :: deltx
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: J23,J12,J13  
  type(real_data_layout)          :: Jdens


 res = .true.

  if (.not. copyStructOnly(T12,J12) .or. &
     &.not. copyStructOnly(T12,J13) .or. &
     &.not. copyStructOnly(T12,J23) .or. &
     &.not. copyStructOnly(T12,Jdens) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in smagorinsky for induction : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
  call computeTensorJdens1(J12,J13,J23,Jdens,res)


  T12%values = 2.0*0.91_WP * (deltx)**2.0 * Jdens%values * J12%values
  T13%values = 2.0*0.91_WP * (deltx)**2.0 * Jdens%values * J13%values
  T23%values = 2.0*0.91_WP * (deltx)**2.0 * Jdens%values * J23%values


    call deleteDataLayout(Jdens)
    call deleteDataLayout(J12)
    call deleteDataLayout(J13)
    call deleteDataLayout(J23)

 end subroutine smagorinskyMHDTij



!!!!!!!!!!simil!!!!!!!!
subroutine SimilInduc(t12,t13,t23,U,V,W,Uk,Vk,Wk,Bx,By,Bz,Bxk,Byk,Bzk,Wave,delta,filter,spec_rank)


implicit none


  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bx,By,Bz
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk
  TYPE(REAL_DATA_LAYOUT) :: Uf,Vf,Wf
  TYPE(REAL_DATA_LAYOUT) :: Bxf,Byf,Bzf 
  INTEGER :: spec_rank,nbcpus,filter
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: t12,t13,t23
  REAL(WP) :: delta
  REAL(WP) :: delta2,Ct
  logical :: res
  TYPE(WaveNumbers), INTENT(IN) :: Wave
  TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2
  TYPE(REAL_DATA_LAYOUT) :: Wtab1
  nbcpus = getnbcpus()
  delta2=2*delta
  Ct=0.91_WP

  IF ( (.NOT. copyStructOnly(U,Uf)) .OR. &
     & (.NOT. copyStructOnly(U,Vf)) .OR. &
     & (.NOT. copyStructOnly(U,Wf)) .OR. &
     & (.NOT. copyStructOnly(U,Wtab1)) .OR. &
     & (.NOT. copyStructOnly(Uk,tabk1)) .OR. &
     & (.NOT. copyStructOnly(Uk,tabk2)) .OR. &
     & (.NOT. copyStructOnly(Bx,Bxf)) .OR. &
     & (.NOT. copyStructOnly(By,Byf)) .OR. &
     & (.NOT. copyStructOnly(Bz,Bzf))    ) RETURN 

    CALL computeFilter(wave,delta2,Uk,tabk2,filter)
    call btran(tabk2,Uf,res)

    CALL computeFilter(wave,delta2,Vk,tabk2,filter)
    call btran(tabk2,Vf,res)

    CALL computeFilter(wave,delta2,Wk,tabk2,filter)
    call btran(tabk2,Wf,res)

    CALL computeFilter(wave,delta2,Bxk,tabk2,filter)
    call btran(tabk2,Bxf,res)

    CALL computeFilter(wave,delta2,Byk,tabk2,filter)
    call btran(tabk2,Byf,res)

    CALL computeFilter(wave,delta2,Bzk,tabk2,filter)
    call btran(tabk2,Bzf,res)

    !calcul de T12

    Wtab1%values = Bx%values*V%values - U%values*By%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t12%values = Ct*(Wtab1%values - Bxf%values*Vf%values + Uf%values*Byf%values) 

   !calcul de T13

    wtab1%values = Bx%values*W%values - U%values*Bz%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t13%values = Ct*(Wtab1%values - Bxf%values*Wf%values + Uf%values*Bzf%values) 

    !calcul de T23

    wtab1%values= By%values*W%values - W%values*Bz%values
    call ftran(wtab1,tabk1,res)
    CALL computeFilter(wave,delta2,tabk1,tabk2,filter)
    call btran(tabk2,Wtab1,res)
    
    t23%values = Ct*(Wtab1%values - Byf%values*Wf%values + Vf%values*Bzf%values) 
    
    
    CALL deleteDataLayout(Uf)
    CALL deleteDataLayout(Vf)
    CALL deleteDataLayout(Wf)
    CALL deleteDataLayout(Bxf)
    CALL deleteDataLayout(Byf)
    CALL deleteDataLayout(Bzf)
    CALL deleteDataLayout(tabk1)
    CALL deleteDataLayout(tabk2)
    CALL deleteDataLayout(Wtab1)


end subroutine SimilInduc

!!!!!!!!!!gradient!!!!!!!!
subroutine InducGradient(U,V,W,Bx,By,Bz,Wave,delta,t12,t13,t23,filter,spec_rank )
  implicit none

  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: U,V,W
  TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bx,By,Bz 
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Uk,Vk,Wk
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Bxk,Byk,Bzk
  TYPE(REAL_DATA_LAYOUT) :: dfibdx, dfibdy, dfibdz
  TYPE(REAL_DATA_LAYOUT) :: dfiudx, dfiudy, dfiudz
  TYPE(REAL_DATA_LAYOUT) :: dfjbdx, dfjbdy, dfjbdz
  TYPE(REAL_DATA_LAYOUT) :: dfjudx, dfjudy, dfjudz
  INTEGER :: spec_rank
  TYPE(REAL_DATA_LAYOUT) :: t12,t13,t23
  logical :: res
  REAL(WP) :: delta
  Integer :: filter,nbcpus
  TYPE(WaveNumbers), INTENT(IN) :: Wave

  nbcpus = getnbcpus()

  IF ((.NOT. copyStructOnly(Bx,dfibdx)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfibdy)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfibdz)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjbdx)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjbdy)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjbdz)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfiudx)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfiudy)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfiudz)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjudx)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjudy)) .OR. &
     &(.NOT. copyStructOnly(Bx,dfjudz))    ) RETURN 

    res = computeFieldGradient(wave, Bx,dfibdx, dfibdy, dfibdz,nbcpus,spec_rank)
    res = computeFieldGradient(wave, By,dfjbdx, dfjbdy, dfjbdz,nbcpus,spec_rank)
    res = computeFieldGradient(wave, U, dfiudx, dfiudy, dfiudz,nbcpus,spec_rank)
    res = computeFieldGradient(wave, V, dfjudx, dfjudy, dfjudz,nbcpus,spec_rank)

    t12%values = (delta*delta/12)*( dfibdx%values*dfjudx%values + dfibdy%values*dfjudy%values + dfibdz%values*dfjudz%values &
                                & - dfjbdx%values*dfiudx%values - dfjbdy%values*dfiudy%values - dfjbdz%values*dfiudz%values )

    res = computeFieldGradient(wave,Bz,dfjbdx, dfjbdy, dfjbdz,nbcpus,spec_rank)
    res = computeFieldGradient(wave, W, dfjudx, dfjudy, dfjudz,nbcpus,spec_rank)


    t13%values = (delta*delta/12)*( dfibdx%values*dfjudx%values + dfibdy%values*dfjudy%values + dfibdz%values*dfjudz%values &
                                & - dfjbdx%values*dfiudx%values - dfjbdy%values*dfiudy%values - dfjbdz%values*dfiudz%values )

    res = computeFieldGradient(wave,By,dfibdx, dfibdy, dfibdz,nbcpus,spec_rank)
    res = computeFieldGradient(wave,V, dfiudx, dfiudy, dfiudz,nbcpus,spec_rank)


    t23%values = (delta*delta/12)*( dfibdx%values*dfjudx%values + dfibdy%values*dfjudy%values + dfibdz%values*dfjudz%values &
                                & - dfjbdx%values*dfiudx%values - dfjbdy%values*dfiudy%values - dfjbdz%values*dfiudz%values )

    CALL deleteDataLayout(dfibdx)
    CALL deleteDataLayout(dfibdy)
    CALL deleteDataLayout(dfibdz)
    CALL deleteDataLayout(dfjbdx)
    CALL deleteDataLayout(dfjbdy)
    CALL deleteDataLayout(dfjbdz)

    CALL deleteDataLayout(dfiudx)
    CALL deleteDataLayout(dfiudy)
    CALL deleteDataLayout(dfiudz)
    CALL deleteDataLayout(dfjudx)
    CALL deleteDataLayout(dfjudy)
    CALL deleteDataLayout(dfjudz)

end subroutine InducGradient

!!!!!!!!!!! Tools !!!!!!!!!!
! subroutine computeMijInducSmag(Bxinf,Byinf,Bzinf,Bxkf,Bykf,Bzkf,t12,t13,t23,wave,&
!                               &M12,M13,M23,delta,filter)
! implicit none
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: t12,t13,t23 
!   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Bxinf,Byinf,Bzinf
!   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: M12,M13,M23
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: Bxkf,Bykf,Bzkf
! 
!   TYPE(WAVENUMBERs), INTENT(IN) :: wave
!   TYPE(REAL_DATA_LAYOUT) :: tab,tabf
!   TYPE(COMPLEX_DATA_LAYOUT) :: tabk,tabkf
!   TYPE(REAL_DATA_LAYOUT) :: t12f,t13f,t23f 
!   TYPE(REAL_DATA_LAYOUT) :: tab1,tab2,tab3
!   REAL(WP), INTENT(IN) :: delta
!   logical :: res
!   integer :: filter,spec_rank
! 
!   IF ((.NOT. copyStructOnly(Bxinf,t12f)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,t13f)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,t23f)) .OR. &
!      &(.NOT. copyStructOnly(Bxkf,tabk)) .OR. &
!      &(.NOT. copyStructOnly(Bxkf,tabkf)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,tab)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,tabf)) ) RETURN 
! 
!      CALL inducSmagotij(Bxkf,Bykf,Bzkf,Bxinf,Byinf,Bzinf,Wave,2.0*delta,t12f,t13f,t23f,spec_rank )
! 
!      CALL ftran(t12,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M12%values = tabf%values - t12f%values
! 
!      CALL ftran(t13,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M13%values = tabf%values - t13f%values
! 
!      CALL ftran(t23,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M23%values = tabf%values - t23f%values
! 
!   call deleteDataLayout(tab)
!   call deleteDataLayout(tabf)
!   call deleteDataLayout(tabk)
!   call deleteDataLayout(tabkf)
! 
!   call deleteDataLayout(t12f)
!   call deleteDataLayout(t13f)
!   call deleteDataLayout(t23f)
! 
! end subroutine computeMijInducSmag

! subroutine computeMijInducGrad(Bxinf,Byinf,Bzinf,Bxkf,Bykf,Bzkf,t12,t13,t23,wave,&
!                               &M12,M13,M23,delta,filter)
! implicit none
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Bxin,Byin,Bzin
!   TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: t12,t13,t23 
!   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: Bxinf,Byinf,Bzinf
!   TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: M12,M13,M23
!   TYPE(COMPLEX_DATA_LAYOUT), INTENT(INOUT) :: Bxkf,Bykf,Bzkf
! 
!   TYPE(WAVENUMBERs), INTENT(IN) :: wave
!   TYPE(REAL_DATA_LAYOUT) :: tab,tabf
!   TYPE(COMPLEX_DATA_LAYOUT) :: tabk,tabkf
!   TYPE(REAL_DATA_LAYOUT) :: t12f,t13f,t23f 
!   TYPE(REAL_DATA_LAYOUT) :: tab1,tab2,tab3
!   REAL(WP), INTENT(IN) :: delta
!   logical :: res
!   integer :: filter,spec_rank
! 
!   IF ((.NOT. copyStructOnly(Bxinf,t12f)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,t13f)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,t23f)) .OR. &
!      &(.NOT. copyStructOnly(Bxkf,tabk)) .OR. &
!      &(.NOT. copyStructOnly(Bxkf,tabkf)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,tab)) .OR. &
!      &(.NOT. copyStructOnly(Bxinf,tabf)) ) RETURN 
! 
!      CALL inducSmagotij(Bxkf,Bykf,Bzkf,Bxinf,Byinf,Bzinf,Wave,2.0*delta,t12f,t13f,t23f,spec_rank )
! 
!      CALL ftran(t12,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M12%values = tabf%values - t12f%values
! 
!      CALL ftran(t13,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M13%values = tabf%values - t13f%values
! 
!      CALL ftran(t23,tabk,res)
!      CALL computeFilter(wave,2.0*delta,tabk,tabkf,Filter) 
!      CALL btran(tabkf,tabf,res)
!      M23%values = tabf%values - t23f%values
! 
!   call deleteDataLayout(tab)
!   call deleteDataLayout(tabf)
!   call deleteDataLayout(tabk)
!   call deleteDataLayout(tabkf)
! 
!   call deleteDataLayout(t12f)
!   call deleteDataLayout(t13f)
!   call deleteDataLayout(t23f)
! 
! end subroutine computeMijInducGrad
! 

 subroutine HambaYoshizawaInduc(Bx,Uk,Vk,Wk,Bxk,Byk,Bzk,C_Lambda,C_mu,wave,spec_rank,res, &
                          &  deltx,rot1k,rot2k,rot3k)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
!   type(REAL_DATA_LAYOUT)   :: rot1,rot2,rot3  
  type(REAL_DATA_LAYOUT)                  :: Wtab1,Wtab2,Wtab3  
  type(COMPLEX_DATA_LAYOUT)               :: Wtab1k,Wtab2k,Wtab3k  
  type(COMPLEX_DATA_LAYOUT), intent(inout):: rot1k,rot2k,rot3k  

  type(REAL_DATA_LAYOUT)                  :: Jx,Jy,Jz  

  type(REAL_DATA_LAYOUT)                  :: Bx  

  real(WP), intent(in)                    :: deltx,C_Lambda,C_mu
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: J23,J12,J13  
  type(real_data_layout)          :: Jdens,normS

  type(real_data_layout)          :: S11,S22,S33 
  type(real_data_layout)          :: S23,S12,S13  
 


 res = .true.

  if (.not. copyStructOnly(Bx,J12) .or. &
     &.not. copyStructOnly(Bx,J13) .or. &
     &.not. copyStructOnly(Bx,J23) .or. &
     &.not. copyStructOnly(Bx,Jdens) )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif
  if (.not. copyStructOnly(Bx,Wtab1) .or. &
     &.not. copyStructOnly(Bx,Wtab2) .or. &
     &.not. copyStructOnly(Bx,Wtab3) )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif
  if (.not. copyStructOnly(Uk,Wtab1k) .or. &
     &.not. copyStructOnly(Uk,Wtab2k) .or. &
     &.not. copyStructOnly(Uk,Wtab3k) )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif
  if (.not. copyStructOnly(Bx,S11) .or. &
     &.not. copyStructOnly(Bx,S22) .or. &
     &.not. copyStructOnly(Bx,S33) .or. &
     &.not. copyStructOnly(Bx,S12) .or. &
     &.not. copyStructOnly(Bx,S13) .or. &
     &.not. copyStructOnly(Bx,S23) .or. &
     &.not. copyStructOnly(Bx,normS) )then

       res = .false.
       write(6,'(a,i0)')'[ERROR] in HambaYoshizawaNSequation : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
  call computeTensorJdens1(J12,J13,J23,Jdens,res)

  call computeStressTensor(Uk,VK,WK,wave,S11,S12,S13,S22,S23,S33,res)
  call computeTensorNorme1(S11,S12,S13,S22,S23,S33,normS,res)

     call computeCurlK(wave,Bxk,Byk,Bzk,Wtab1k,Wtab2k,Wtab3k)

     call btran(Wtab1k,Jx,res)
     call btran(Wtab2k,Jy,res)
     call btran(Wtab3k,Jz,res)

  Wtab1%values= -C_lambda*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*Jx%values 
  Wtab2%values= -C_lambda*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*Jy%values
  Wtab3%values= -C_lambda*deltx**2 *sqrt((0.5_WP*C_mu*(normS%values**2) + C_lambda*Jdens%values**2))*Jz%values

     call ftran(Wtab1,Wtab1k,res)
     call ftran(Wtab2,Wtab2k,res)
     call ftran(Wtab3,Wtab3k,res)

     call computeCurlK(wave,Wtab1k,Wtab2k,Wtab3k,rot1k,rot2k,rot3k)

!      call btran(rot1k,rot1,res)
!      call btran(rot2k,rot2,res)
!      call btran(rot3k,rot3,res)


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

    call deleteDataLayout(Wtab1)
    call deleteDataLayout(Wtab2)
    call deleteDataLayout(Wtab3)

    call deleteDataLayout(Wtab1k)
    call deleteDataLayout(Wtab2k)
    call deleteDataLayout(Wtab3k)

 end subroutine HambaYoshizawaInduc

 subroutine MullerCarratiInduc(Uk,Vk,Wk,Bxk,Byk,Bzk,wave,spec_rank,res, &
                          &  deltx,T23,T12,T13)


  !I/O
  type(complex_data_layout),intent(in)    :: Uk,Vk,Wk 
  type(WaveNumbers), intent(in)           :: wave
  integer, intent(in)                     :: spec_rank
  logical, intent(inout)                  :: res
!   type(REAL_DATA_LAYOUT)   :: rot1,rot2,rot3  
  type(REAL_DATA_LAYOUT)                  :: Wtab1,Wtab2,Wtab3  
  type(COMPLEX_DATA_LAYOUT)               :: Wtab1k,Wtab2k,Wtab3k  
  type(REAL_DATA_LAYOUT), intent(inout)   :: T23,T12,T13  

  type(REAL_DATA_LAYOUT)                  :: Jx,Jy,Jz  
  type(REAL_DATA_LAYOUT)                  :: Vortx,Vorty,Vortz  

!   type(REAL_DATA_LAYOUT)                  :: Bx  

  real(WP), intent(in)                    :: deltx
  TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)   :: Bxk,Byk,Bzk


  !Local data
  type(real_data_layout)          :: J23,J12,J13  
!   type(real_data_layout)          :: Jdens,normS

!   type(real_data_layout)          :: S11,S22,S33 
!   type(real_data_layout)          :: S23,S12,S13  
 


 res = .true.

  if (.not. copyStructOnly(T12,J12) .or. &
     &.not. copyStructOnly(T12,J13) .or. &
     &.not. copyStructOnly(T12,J23)  )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiInduc : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(T12,Wtab1) .or. &
     &.not. copyStructOnly(T12,Wtab2) .or. &
     &.not. copyStructOnly(T12,Wtab3) )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiInduc : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(Uk,Wtab1k) .or. &
     &.not. copyStructOnly(Uk,Wtab2k) .or. &
     &.not. copyStructOnly(Uk,Wtab3k) )then
     res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiInduc : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(T12,Jx) .or. &
     &.not. copyStructOnly(T12,Jy) .or. &
     &.not. copyStructOnly(T12,Jz) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiInduc : not enought memory !',spec_rank
   endif

  if (.not. copyStructOnly(T12,Vortx) .or. &
     &.not. copyStructOnly(T12,Vorty) .or. &
     &.not. copyStructOnly(T12,Vortz) )then
       res = .false.
       write(6,'(a,i0)')'[ERROR] in MullerCarratiInduc : not enought memory !',spec_rank
   endif

  call computeJijTensor(BxK,ByK,BzK,wave,J12,J13,J23, res)
!   call computeTensorJdens1(J12,J13,J23,Jdens,res)

!   call computeStressTensor(Uk,VK,WK,wave,S11,S12,S13,S22,S23,S33,res)
!   call computeTensorNorme1(S11,S12,S13,S22,S23,S33,normS,res)

     call computeCurlK(wave,Bxk,Byk,Bzk,Wtab1k,Wtab2k,Wtab3k)

     call btran(Wtab1k,Jx,res)
     call btran(Wtab2k,Jy,res)
     call btran(Wtab3k,Jz,res)

     call computeCurlK(wave,Uk,Vk,Wk,Wtab1k,Wtab2k,Wtab3k)

     call btran(Wtab1k,Vortx,res)
     call btran(Wtab2k,Vorty,res)
     call btran(Wtab3k,Vortz,res)

  T12%values= -0.015*deltx**2 *sign(1.0_WP,(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*&
                           &sqrt(abs(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*J12%values 
  T13%values= -0.015*deltx**2 *sign(1.0_WP,(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*&
                           &sqrt(abs(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*J13%values
  T23%values= -0.015*deltx**2 *sign(1.0_WP,(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*&
                           &sqrt(abs(Jx%values*Vortx%values + Jy%values*Vorty%values + Jz%values*Vortz%values))*J23%values

!     call deleteDataLayout(Jdens)
    call deleteDataLayout(J12)
    call deleteDataLayout(J13)
    call deleteDataLayout(J23)

    call deleteDataLayout(Wtab1)
    call deleteDataLayout(Wtab2)
    call deleteDataLayout(Wtab3)

    call deleteDataLayout(Wtab1k)
    call deleteDataLayout(Wtab2k)
    call deleteDataLayout(Wtab3k)

    call deleteDataLayout(Jx)
    call deleteDataLayout(Jy)
    call deleteDataLayout(Jz)

    call deleteDataLayout(Vortx)
    call deleteDataLayout(Vorty)
    call deleteDataLayout(Vortz)

 end subroutine MullerCarratiInduc

end module

