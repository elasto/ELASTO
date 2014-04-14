!> @addtogroup post_process
!! @{

!------------------------------------------------------------------------------

!
! MODULE: post_mhd
!
!
! DESCRIPTION: 
!>  This module provide a library of post-process routines, specifics to kolmogorov flow simulations. This library
!! provide all that you need to build your own post-process.
!
!> @details 
!! 
!> @author
!! Mouloud KESSAR, LEGI
!
!------------------------------------------------------------------------------

module post_kolmo

    use datalayout
    use toolbox
    use stat_tools 
    use differential_tools 
    USE wavenumber_tools 
    use transforms_tools
    use post_lib2D
    implicit none
    private

   ! ==== Public procedures ==== 
    public          :: post_kolmo_init
    public          :: compute_hydro_stat
    public          :: compute_Mhd_stat
    ! ======= Mean terms ========
    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanU, meanUV,meanDU,meanDU2
    integer, save                                          :: save_freq       ! frequency for saving statistics

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanB, meanBiBj,meanDB,meanDB2

    real(WP), PROTECTED, ALLOCATABLE, dimension(:,:), SAVE :: meanBJ, meanUW,meanUiBj

contains

subroutine post_kolmo_init(U,imhd)

    implicit none

    TYPE(REAL_DATA_LAYOUT), intent(in)    :: U
    logical :: imhd

    CALL parser_read('Save frequency',save_freq)

    allocate(meanU(1:U%nz,3))
    allocate(meanUV(1:U%nz,6))
    allocate(meanDU(1:U%nz,9))
    allocate(meanDU2(1:U%nz,9))

    meanU(:,:) = 0.0_WP
    meanUV(:,:) = 0.0_WP
    meanDU(:,:) = 0.0_WP
    meanDU2(:,:) = 0.0_WP

   if (imhd) Then
    allocate(meanB(1:U%nz,3))
    allocate(meanBiBj(1:U%nz,6))
    allocate(meanDB(1:U%nz,9))
    allocate(meanDB2(1:U%nz,9))

    allocate(meanBJ(1:U%nz,3))
    allocate(meanUW(1:U%nz,3))
    allocate(meanUiBj(1:U%nz,6))

   endif
end subroutine post_kolmo_init

subroutine compute_hydro_stat(Uin,Vin,Win,Ukin,Vkin,Wkin,WaveVel,spec_rank,simtime,tstep,ite, save_freq)

    use fileio

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: Uin,Vin,Win
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: Ukin,Vkin,Wkin
    TYPE(WaveNumbers),INTENT(IN) :: WaveVel
    INTEGER, INTENT(IN) :: spec_rank,ite, save_freq
    REAL(WP), INTENT(IN) :: simtime,tstep

    INTEGER :: k
    REAL(WP) :: timesimu, dt
    REAL(WP),allocatable, dimension (:,:) :: buffer_meanU, buffer_meanUV,buffer_meanDU, buffer_meanDU2
    integer                                 :: file_id          ! id for file io
    character(len=50)                       :: file_name        ! name of output file


    dt=tstep
    timesimu=simtime

    ! Compute Mean veolcities component
    ! spatial averaging

    allocate(buffer_meanU(1:Uin%nz,3))
    allocate(buffer_meanUV(1:Uin%nz,6))
    allocate(buffer_meanDU(1:Uin%nz,9))
    allocate(buffer_meanDU2(1:Uin%nz,9))

    call compute_Vec_stat(Uin,Vin,Win,Ukin,Vkin,Wkin,WaveVel,spec_rank,buffer_meanU, buffer_meanUV &
                          &, buffer_meanDU, buffer_meanDU2 )

    if (simtime.eq.0) then
        meanUV(:,:) = buffer_meanUV(:,:)
        meanU(:,:) = buffer_meanU(:,:)
        meanDU(:,:) = buffer_meanDU(:,:)
        meanDU2(:,:) = buffer_meanDU2(:,:)
    else
        meanUV(:,:) = (meanUV(:,:)*(simtime-dt) + buffer_meanUV(:,:)*dt)/simtime
        meanU(:,:) = (meanU(:,:)*(simtime-dt) + buffer_meanU(:,:)*dt)/simtime
        meanDU(:,:) = (meanDU(:,:)*(simtime-dt) + buffer_meanDU(:,:)*dt)/simtime
        meanDU2(:,:) = (meanDU2(:,:)*(simtime-dt) + buffer_meanDU2(:,:)*dt)/simtime
    end if




    if((save_freq>0).and.(mod(ite, save_freq)==0)) then
       if (spec_rank.eq.0) then

          file_id = iopen()
           write(file_name,'(a)') 'U_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(5(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanU(k,1), meanU(k,2), meanU(k,3), simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'UU_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(8(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanUV(k,1), meanUV(k,2), meanUV(k,3),&
                                               & meanUV(k,4), meanUV(k,5), meanUV(k,6),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'DU_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(11(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanDU(k,1), meanDU(k,2), meanDU(k,3),&
                                               & meanDU(k,4), meanDU(k,5), meanDU(k,6),meanDU(k,7), &
                                               & meanDU(k,8), meanDU(k,9),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)


          file_id = iopen()
           write(file_name,'(a)') 'DU2_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(11(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanDU2(k,1), meanDU2(k,2), meanDU2(k,3),&
                                               & meanDU2(k,4), meanDU2(k,5), meanDU2(k,6),meanDU2(k,7), &
                                               & meanDU2(k,8), meanDU2(k,9),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

       end if
    end if

    deallocate(buffer_meanUV)
    deallocate(buffer_meanU)

    deallocate(buffer_meanDU)
    deallocate(buffer_meanDU2)

end subroutine compute_hydro_stat

subroutine compute_Vec_stat(VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,spec_rank,&
                          & buffer_meanVi, buffer_meanViVj, buffer_meanDjVi, buffer_meanDjVi2 )

!     use fileio

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    TYPE(WaveNumbers),INTENT(IN) :: Wave
    INTEGER, INTENT(IN) :: spec_rank

    REAL(WP),dimension (:,:) :: buffer_meanVi, buffer_meanViVj
    REAL(WP),dimension (:,:) :: buffer_meanDjVi, buffer_meanDjVi2
    TYPE(REAL_DATA_LAYOUT) :: ViVj,Wtab
    TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3
    logical :: res1

    IF ( .NOT. copyStructOnly(VecX,ViVj) ) RETURN

    ! Compute Mean veolcities component
    ! spatial averaging
    if(.not. computeFieldAvg2D(VecX,buffer_meanVi(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for Uin'
    if(.not. computeFieldAvg2D(VecY,buffer_meanVi(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for Vin'
    if(.not. computeFieldAvg2D(VecZ,buffer_meanVi(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for Win'


    ViVj%values = VecX%values*VecX%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UU'

    ViVj%values = VecX%values*VecY%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UV'

    ViVj%values = VecX%values*VecZ%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UW'

    ViVj%values = VecY%values*VecY%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,4),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for VV'

    ViVj%values = VecY%values*VecZ%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,5),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for VW'

    ViVj%values = VecZ%values*VecZ%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanViVj(:,6),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for WW'



     IF ((.NOT. copyStructOnly(VecX,Wtab)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk1)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk2)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk3)) ) RETURN

    !!!!!DVecXDj et (DVecXDj)2
    call computeGradientK(wave,VecXk,tabk1,tabk2,tabk3)

    call btran(tabk1,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDx'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDx2'

    call btran(tabk2,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDy'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDy2'

    call btran(tabk3,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDz'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecXDz2'

    !!!!!DVecYDj et (DVecYDj)2
    CALL computeGradientK(wave,VecYk,tabk1,tabk2,tabk3)

    call btran(tabk1,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,4),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDx'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,4),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDx2'

    call btran(tabk2,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,5),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDy'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,5),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDy2'

    call btran(tabk3,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,6),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDz'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,6),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecYDz2'

    !!!!!DVecZDj et (DVecZDj)2
    CALL computeGradientK(wave,VecZk,tabk1,tabk2,tabk3)

    call btran(tabk1,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,7),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDx'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,7),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDx2'

    call btran(tabk2,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,8),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDy'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,8),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDy2'

    call btran(tabk3,Wtab,res1)
    if(.not. computeFieldAvg2D(Wtab,buffer_meanDjVi(:,9),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDz'

    ViVj%values = Wtab%values*Wtab%values
    if(.not. computeFieldAvg2D(ViVj,buffer_meanDjVi2(:,9),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for DVecZDz2'


    CALL deletedatalayout(ViVj)
    CALL deletedatalayout(tabk1)
    CALL deletedatalayout(tabk2)
    CALL deletedatalayout(tabk3)
    CALL deletedatalayout(Wtab)


end subroutine compute_Vec_stat

subroutine compute_Mhd_stat(Uin,Vin,Win,Ukin,Vkin,Wkin,Bin,Bkin,wave,spec_rank,simtime,tstep,ite, save_freq)

    use fileio

   implicit none
   integer, intent(in)    :: spec_rank,ite,save_freq
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
   TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN)          :: Ukin,Vkin,Wkin
   TYPE(REAL_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bin
   TYPE(COMPLEX_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bkin
   TYPE(WaveNumbers), INTENT(IN) :: wave

   INTEGER :: k
   REAL(WP), INTENT(IN) :: simtime,tstep
   REAL(WP),allocatable, dimension (:,:) :: buffer_meanBi, buffer_meanBiBj,buffer_meanDjBi, buffer_meanDjBi2
   REAL(WP),allocatable, dimension (:,:) :: buffer_meanBJ, buffer_meanUiBj,buffer_meanUW
   integer                                 :: file_id          ! id for file io
   character(len=50)                       :: file_name        ! name of output file

    allocate(buffer_meanBi(1:Uin%nz,3))
    allocate(buffer_meanBiBj(1:Uin%nz,6))
    allocate(buffer_meanDjBi(1:Uin%nz,9))
    allocate(buffer_meanDjBi2(1:Uin%nz,9))

    allocate(buffer_meanBJ(1:Uin%nz,3))
    allocate(buffer_meanUW(1:Uin%nz,3))
    allocate(buffer_meanUiBj(1:Uin%nz,6))

   !!Termes liées purement au champ de vitesse
   call compute_hydro_stat(Uin,Vin,Win,Ukin,Vkin,Wkin,wave,spec_rank,simtime,tstep,ite,save_freq)

   call compute_ArotA_stat(Uin,Vin,Win,Ukin,Vkin,Wkin,Wave,spec_rank,&
                          & buffer_meanUW )

   !!!!! termes liés au champ B
   call compute_Vec_stat(Bin(1),Bin(2),Bin(3),Bkin(1),Bkin(2),Bkin(3),Wave,spec_rank,&
                    & buffer_meanBi, buffer_meanBiBj, buffer_meanDjBi, buffer_meanDjBi2 )

   call compute_ArotA_stat(Bin(1),Bin(2),Bin(3),Bkin(1),Bkin(2),Bkin(3),Wave,spec_rank,&
                          & buffer_meanBJ )
   !!! termes croisés
   call compute_UB_stat(Uin,Vin,Win,Bin,spec_rank,buffer_meanUiBj)

    if (simtime.eq.0) then
        meanBiBj(:,:) = buffer_meanBiBj(:,:)
        meanB(:,:) = buffer_meanBi(:,:)
        meanDB(:,:) = buffer_meanDjBi(:,:)
        meanDB2(:,:) = buffer_meanDjBi2(:,:)

        meanBJ(:,:) = buffer_meanBJ(:,:)
        meanUiBj(:,:) = buffer_meanUiBj(:,:)
        meanUW(:,:) = buffer_meanUW(:,:)

    else
        meanBiBj(:,:) = (meanBiBj(:,:)*(simtime-tstep) + buffer_meanBiBj(:,:)*tstep)/simtime
        meanB(:,:) = (meanB(:,:)*(simtime-tstep) + buffer_meanBi(:,:)*tstep)/simtime
        meanDB(:,:) = (meanDB(:,:)*(simtime-tstep) + buffer_meanDjBi(:,:)*tstep)/simtime
        meanDB2(:,:) = (meanDB2(:,:)*(simtime-tstep) + buffer_meanDjBi2(:,:)*tstep)/simtime

        meanBJ(:,:) = (meanBJ(:,:)*(simtime-tstep) + buffer_meanBJ(:,:)*tstep)/simtime
        meanUiBj(:,:) = (meanUiBj(:,:)*(simtime-tstep) + buffer_meanUiBj(:,:)*tstep)/simtime
        meanUW(:,:) = (meanUW(:,:)*(simtime-tstep) + buffer_meanUW(:,:)*tstep)/simtime

    end if



    if((save_freq>0).and.(mod(ite, save_freq)==0)) then
       if (spec_rank.eq.0) then

          file_id = iopen()
           write(file_name,'(a)') 'B_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(5(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanB(k,1), meanB(k,2), meanB(k,3), simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'BB_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(8(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanUV(k,1), meanBiBj(k,2), meanBiBj(k,3),&
                                               & meanBiBj(k,4), meanBiBj(k,5), meanBiBj(k,6),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'DB_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(11(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanDB(k,1), meanDB(k,2), meanDB(k,3),&
                                               & meanDB(k,4), meanDB(k,5), meanDB(k,6),meanDB(k,7), &
                                               & meanDB(k,8), meanDB(k,9),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)


          file_id = iopen()
           write(file_name,'(a)') 'DU2_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(11(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanDB2(k,1), meanDB2(k,2), meanDB2(k,3),&
                                               & meanDB2(k,4), meanDB2(k,5), meanDB2(k,6),meanDB2(k,7), &
                                               & meanDB2(k,8), meanDB2(k,9),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'BJ_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(5(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanBJ(k,1), meanBJ(k,2), meanBJ(k,3),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'UW_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(5(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanUW(k,1), meanUW(k,2), meanUW(k,3),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

          file_id = iopen()
           write(file_name,'(a)') 'UB_MEAN.dat'
           file_name = trim(adjustl(file_name))
           open(unit=file_id, file=file_name, form="FORMATTED",position="append")
               do k = 1,Uin%nz
                  write(file_id,'(8(f25.10))' ) (k-1)*Uin%Lz/Uin%nz,meanUiBj(k,1), meanUiBj(k,2), meanUiBj(k,3),&
                                              & meanUiBj(k,4), meanUiBj(k,5), meanUiBj(k,6),simtime
               end do
           close(file_id)
           file_id = iclose(file_id)

       end if
    end if


    deallocate(buffer_meanBi)
    deallocate(buffer_meanBiBj)
    deallocate(buffer_meanDjBi)
    deallocate(buffer_meanDjBi2)

    deallocate(buffer_meanBJ)
    deallocate(buffer_meanUW)
    deallocate(buffer_meanUiBj)


end subroutine compute_Mhd_stat

subroutine compute_UB_stat(Uin,Vin,Win,Bin,spec_rank,buffer_meanUiBj)

   implicit none
   integer, intent(in)    :: spec_rank
   TYPE(REAL_DATA_LAYOUT), INTENT(IN)          :: Uin,Vin,Win
   TYPE(REAL_DATA_LAYOUT),dimension(:), INTENT(IN) :: Bin(:)

   REAL(WP),dimension (:,:) :: buffer_meanUiBj
   TYPE(REAL_DATA_LAYOUT) :: UiBj

    IF ( .NOT. copyStructOnly(Uin,UiBj) ) RETURN

    ! Compute Mean veolcities component

    ! spatial averaging

    UiBj%values = Uin%values*Bin(2)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UxBy'

    UiBj%values = Uin%values*Bin(3)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UxBz'

    UiBj%values = Vin%values*Bin(1)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UyBx'

    UiBj%values = Vin%values*Bin(3)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,4),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UyBz'

    UiBj%values = Win%values*Bin(1)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,5),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UzBx'

    UiBj%values = Win%values*Bin(2)%values
    if(.not. computeFieldAvg2D(UiBj,buffer_meanUiBj(:,6),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for UzBy'

    CALL deletedatalayout(UiBj)
end subroutine compute_UB_stat

subroutine compute_ArotA_stat(VecX,VecY,VecZ,VecXk,VecYk,VecZk,Wave,spec_rank,&
                          & buffer_meanArotA )

!     use fileio

    implicit none

    TYPE(REAL_DATA_LAYOUT), INTENT(IN) :: VecX,VecY,VecZ
    TYPE(COMPLEX_DATA_LAYOUT), INTENT(IN) :: VecXk,VecYk,VecZk
    TYPE(WaveNumbers),INTENT(IN) :: Wave
    Integer :: spec_rank
    REAL(WP),dimension (:,:) :: buffer_meanArotA
    TYPE(REAL_DATA_LAYOUT) :: Wtab1,Wtab2,Wtab3,Wtab
    TYPE(COMPLEX_DATA_LAYOUT) :: tabk1,tabk2,tabk3
    logical :: res1

     IF ((.NOT. copyStructOnly(VecX,Wtab1)) .OR. &
        &(.NOT. copyStructOnly(VecX,Wtab2)) .OR. &
        &(.NOT. copyStructOnly(VecX,Wtab3)) .OR. &
        &(.NOT. copyStructOnly(VecX,Wtab)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk1)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk2)) .OR. &
        &(.NOT. copyStructOnly(VecXk,tabk3)) ) RETURN

    CALL computeCurlKmain(wave,VecXk,VecYk,VecZk,tabk1,tabk2,tabk3)
    CALL btran(tabk1,Wtab1,res1)
    CALL btran(tabk2,Wtab2,res1)
    CALL btran(tabk3,Wtab3,res1)

    Wtab%values=VecX%values*Wtab1%values
        if(.not. computeFieldAvg2D(Wtab,buffer_meanArotA(:,1),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for AxrotAx'

    Wtab%values=VecY%values*Wtab2%values
        if(.not. computeFieldAvg2D(Wtab,buffer_meanArotA(:,2),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for AyrotAy'

    Wtab%values=VecZ%values*Wtab3%values
        if(.not. computeFieldAvg2D(Wtab,buffer_meanArotA(:,3),spec_rank)) write(*,'(a)')&
        & '[ERROR] Could not compute AVG for AzrotAz'

    CALL deletedatalayout(tabk1)
    CALL deletedatalayout(tabk2)
    CALL deletedatalayout(tabk3)
    CALL deletedatalayout(Wtab)
    CALL deletedatalayout(Wtab1)
    CALL deletedatalayout(Wtab2)
    CALL deletedatalayout(Wtab3)

end subroutine compute_ArotA_stat

end module post_kolmo
!> @}

