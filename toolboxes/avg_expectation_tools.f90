!USEFORTEST avgcond
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: conditional_mean 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> The aim of this module is to provide elementaries piece of code which can be use elsewhere.
!> BEWARE there are no systematic tests in these routines in order to not degrade performances
!> 
!------------------------------------------------------------------------------


module conditional_mean_tools

   use datalayout
   use parallel_tools
   use precision_tools
   use stat_tools 
   use mpilayout_tools

   IMPLICIT NONE

    PRIVATE
    
! Structure de données contenant un estimateur optimal  
!    type estimOp
!      type(varCond),dimension(:)       :: condVar 
!      real(WP),dimension(:)            :: estimOptDisc 
!      type(real_data_layout)           :: fieldOpt
!      integer(DI),dimension(:,:)       :: pdfCond
!      integer(DI),dimension(:)         :: pdfField
!    end type  estimOp

! Structure de données décrivant une variable conditionnelle
    TYPE varCond
     integer                                  :: typeDisc !1 pour cst !2 pour variable
     integer                                  :: bin
     integer(DI),dimension(:),allocatable     :: nbValues 
     real(WP)                                 :: minCond, maxCond
     real(WP)                                 :: delta
     real(WP)                                 :: deltaTot
     real(WP) , dimension(:,:,:), allocatable :: valCond
     real(WP) , dimension(:),allocatable      :: valCondDisc
    END TYPE varCond

! Variable servant de test d'initialisation
    logical                               :: init = .false.
    logical                               :: initPdf = .false.
! Variable servant à connaitre l'état de l'intervalle forcé(min/max) par rapport à
! l'intervalle réel de la variable de conditionnement  
    logical                               :: interSousDim = .false.
! Variable contenant le nombre de variables conditionnelles
    integer                               :: nbCondTot
! Indices des noeuds min/max locaux au processus
    integer                               :: xminLoc,xmaxLoc
    integer                               :: yminLoc,ymaxLoc
    integer                               :: zminLoc,zmaxLoc
! Nombre de noeuds locaux au processus
    integer                               :: nxLoc,nyLoc,nzLoc
! Nombre de noeuds totale
    integer                               :: nxGlob,nyGlob,nzGlob
! Variable utilisée lors de l'initialisation des variables
    integer                               :: nbCondCurrent
! Variable contenant le nombre d'intervalles sur lesquels on va
! calculer les moyennes conditionnelles
    integer                               :: allBin
! Vecteur contenant les coefficients pour la decomposition (resp.
! composition) de l'espace de dimension nbCondTot à 1 (resp. 1 à
! nbCondTot)
! integer, dimension (:), allocatable ,save :: tabExposant
    integer, dimension (:), allocatable   :: tabExposant
! Tableau du champ à moyenner
    real(WP),dimension(:,:,:),allocatable :: fieldIn
! Tableau du champ moyenné
    real(WP),dimension(:,:,:),allocatable :: fieldOut
! Tableau des valeurs du champ de sortie
    real(WP) , dimension(:) , allocatable :: valuesFieldOut
! Tableau du nombre des valeurs du champ de sortie
   ! real(WP) , dimension(:) , allocatable :: nbValuesFieldOut
    integer(DI) , dimension(:) , allocatable :: nbValuesFieldOut

! Structure des données contenant les variables de conditionnement
    type(varCond),dimension(:),allocatable:: condAll

    PUBLIC :: init_meanCond         !initialiser le module
    PUBLIC :: finalize_meanCond     !detruire le module
    PUBLIC :: setConds              !donner les variables de conditionnement
    PUBLIC :: setField              !donner le champs a moyenner
    PUBLIC :: compute_meanCond      !calcul de la moyenne conditionnelle
    PUBLIC :: computePdfConds       !calcul les pdf des variables conditionnelles
    PUBLIC :: showFieldDisc         !afficher les valeurs du champ moyenné 
    PUBLIC :: showConnectTable      !afficher la table de connexion entre la nième val du champ moyenné
                                    ! et les numéros des bin de chaque variable conditionnelle 
    PUBLIC :: showPdfCond           !afficher les pdf des variables conditionnelles
    PUBLIC :: showCondDisc          !afficher les variables conditionnelles discrétisées
    PUBLIC :: showPdfConjointe      !afficher la pdfconjointe de toutes les variables conditionnelles
    PUBLIC :: getMeanField          !recuperer le champ moyenne sur les noeuds du domaine de calcul
    PUBLIC :: getFieldDisc          !récupérer les valeurs du champs moyenne discretise
    PUBLIC :: getCondDisc           !récupérer les valeurs des variables conditionnelles
    PUBLIC :: getConnectTable       !récupérer la table de connexion
    PUBLIC :: getPdfConjointe       !récupérer la pdf conjointe
    PUBLIC :: getPdfCond            !récupérer la pdf de chaque variable conditionnelle
    PUBLIC :: computeEONVar 
    PUBLIC :: computeSlideAvg
    PUBLIC :: readModelOpt 
    PUBLIC :: computeNbLine 
    PUBLIC :: buildEOFieldFromEO
    PUBLIC :: computeEoObjectives

  CONTAINS

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Initialize the conditional mean context
!> @details
!>
!> @param [in] nbCondExt is the number of variable 
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------
  function init_meanCond( nbCondExt,spec_rank ) result (success)
    !I/O
    integer , intent (in)                      :: nbCondExt
    integer , intent (in)                      :: spec_rank
    logical                                    :: success

    !Local data
    integer :: check 

    success = .false.

!    if ( spec_rank .eq. 0) write(6,'(a)')'[INFO AVG EXPECTATION] -> init_meanCond'
    if ( nbCondExt .lt. 1) then
      write(6,'(a)') '[ERROR] Le nombre de condition doit etre superieur a 1.'
      return
    else
      nbCondTot = nbCondExt
      allocate (condAll( nbCondTot ) , stat=check)
      if  ( check .ne. 0 ) then
        write(6,'(a)') '[ERROR] Not enought memory'
        return
      endif
      allocate (tabExposant( nbCondTot ), stat=check)
      if (check .ne. 0) then
        write(6,'(a)') '[ERROR] Not enought memory'
        return
      endif
      init = .true. 
      nbCondCurrent = 0
    endif

    success = .true.

  end function init_meanCond


!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Delete the conditional mean context
!------------------------------------------------------------------------------
  subroutine finalize_meanCond()

    !Local data
    integer               :: i

!    if ( getnbmycpu() .eq. 0) write(6,'(a)')'[INFO AVG EXPECTATION] -> finalizing...'
    init = .false. 
    initPdf = .false.
    interSousDim = .false.
    nbCondTot = 0
    nbCondCurrent = 1 
    do i=1,nbCondTot
      deallocate (condAll(i)%nbValues )
      deallocate (condAll(i)%valCond )
      deallocate (condAll(i)%valCondDisc ) 
    enddo
    deallocate (condAll)
    deallocate (fieldIn)
    deallocate (fieldOut)
    deallocate (valuesFieldOut)
    deallocate (nbValuesFieldOut)
    deallocate (tabExposant)

  end subroutine finalize_meanCond

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!>
!> @details
!>
!> @param [in] field champs scalaire a moyenner  
!> @param [in] nbVar nombre de variables conditionnelles  
!> @param [in] Var vecteur contenant les variables conditionnelles 
!> @param [in] 
!> @param [in,out] ScalArray a pointer to an array of scalars. 
!> @param [in] nbcpus the total amount of processes (or cpus)
!> @param [in] me my MPI processe rank in the range [0:nbcpus-1]
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------

  function setConds( condext , typeDisc , spec_rank , bin ) result(success)

    use parallel_tools , only : DoGlobalMax
    use parallel_tools , only : DoGlobalMin

    !I/O data
    type(real_data_layout) , intent (in)  :: condext
    integer , intent (in)                 :: typeDisc
    integer , intent (in)                 :: spec_rank 
    integer , intent (in)                 :: bin
    logical                               :: success 

    !Local data
    integer                               :: i,j,k
    integer                               :: check 

    success = .false. 

    nbCondCurrent = nbCondCurrent + 1
    if ( condext%xmin .ne. xminLoc .or. &
       & condext%xmax .ne. xmaxLoc .or. &
       & condext%ymin .ne. yminLoc .or. &
       & condext%ymax .ne. ymaxLoc .or. &
       & condext%zmin .ne. zminLoc .or. &
       & condext%zmax .ne. zmaxLoc ) then
       write(6,'(a)')'[ERROR] in getMeanField, fields are not the same'
       return
    endif
!    if (spec_rank .eq. 0) write(6,'(a,1x,i3)')'[INFO AVG EXPECTATION] -> setConds numero:',nbCondCurrent
    if (.not. init) then
       write(6,'(a)') '[ERROR] meanCond pas initialise'
       return
    elseif ( nbCondCurrent .gt. nbCondTot ) then
       write(6,'(a)') '[ERROR] le nombre de condition depasse celui initialise'
       return
    elseif ( bin .lt. 1 ) then
       write(6,'(a)') '[ERROR] le nombre de bin est incorrect' 
       return
    else
       allocate(condAll(nbCondCurrent)%valCond(xminLoc:xmaxLoc,yminLoc:ymaxLoc,&
               &zminLoc:zmaxLoc),stat=check)
       if ( check .ne. 0 )then
         write(6,'(a)') '[ERROR] Not enought memory'
         return
       endif
    endif 
!print *,'condext%values ( condext%xmax , condext%ymax , condext%zmax )',&
!       &condext%values ( condext%xmax , condext%ymax , condext%zmax )
    do k = condext%zmin,condext%zmax
      do j = condext%ymin,condext%ymax 
        do i = condext%xmin,condext%xmax 
          condAll( nbCondCurrent )%valCond( i , j , k ) = condext%values ( i , j , k )
        enddo
      enddo
    enddo
    condAll( nbCondCurrent )%minCond=doglobalmin(spec_rank,minval(condext%values))
    condAll( nbCondCurrent )%maxCond=doglobalmax(spec_rank,maxval(condext%values))
    if ( equalToVal(condAll( nbCondCurrent )%minCond , condAll( nbCondCurrent )%maxCond )) then
         write(6,'(a,1x,i0,1x,a)') '[ERROR] Conditional variable',nbCondCurrent,' is a uniform field'
         return
    endif
    if ( typeDisc .eq. 1) then
       condAll( nbCondCurrent )%bin = bin
       condAll( nbCondCurrent )%typeDisc = typeDisc
if (bin.lt. 1) print *,'ERROR 1'
       allocate( condAll(nbCondCurrent)%nbValues( bin ), stat=check )
       allocate( condAll(nbCondCurrent)%valCondDisc( bin ), stat=check )
       condAll(nbCondCurrent)%nbValues = 0_DI
       if ( check .ne. 0) then
         write(6,'(a,i0,a)') '[ERROR] Not able to allocate nbValues at condVal ',nbCondCurrent,'&
                             & Not enought memory'
         return
       endif
    elseif ( typeDisc .eq. 2) then
       condAll( nbCondCurrent )%bin = bin
       condAll( nbCondCurrent )%typeDisc = typeDisc
if (bin.lt. 1) print *,'ERROR 2'
       allocate( condAll(nbCondCurrent)%nbValues( bin ), stat=check )
       allocate( condAll(nbCondCurrent)%valCondDisc( bin ), stat=check )
       condAll(nbCondCurrent)%nbValues = 0_DI
       if ( check .ne. 0) then
         write(6,'(a,i0,a)') '[ERROR] Not able to allocate nbValues at condVal ',nbCondCurrent,'&
                             & Not enought memory'
         return
       endif
    else
       if (spec_rank .eq. 0) then
          write(6,'(a)')'[ERROR] type de discretisation non valide'
          return
       endif
    endif
    success = .true. 

  end function setConds

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Initialize the field to estimate
!>
!> @details
!> @param [in] field champs scalaire a moyenner  
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------

  function setField( field ) result (success)

    !I/O data
    type(real_data_layout), intent (in)                     :: field
    logical                                                 :: success

    !Local data
    integer :: i,j,k
    integer :: check 

    success = .false.

    if (.not. init ) then
      write(6,'(a)')'[ERROR] meanCond pas initialise'
      return
    else
      xminLoc=field%xmin
      xmaxLoc=field%xmax
      yminLoc=field%ymin
      ymaxLoc=field%ymax
      zminLoc=field%zmin
      zmaxLoc=field%zmax
      nxGlob=field%nx
      nyGlob=field%ny
      nzGlob=field%nz
      nxLoc=xmaxLoc-xminLoc+1
      nyLoc=ymaxLoc-yminLoc+1
      nzLoc=zmaxLoc-zminLoc+1
      allocate(fieldIn(xminLoc:xmaxLoc,yminLoc:ymaxLoc,zminLoc:zmaxLoc),stat=check)
      if (check.ne.0) then
        write(6,'(a)') '[ERROR] Not enought memory'
        return
      endif 
      allocate(fieldOut(xminLoc:xmaxLoc,yminLoc:ymaxLoc,zminLoc:zmaxLoc),stat=check)
      if (check.ne.0) then
        write(6,'(a)') '[ERROR] Not enought memory'
        return
      endif 
      do k = zminLoc,zmaxLoc 
        do j = yminLoc,ymaxLoc 
          do i = xminLoc,xmaxLoc
            fieldIn( i , j , k ) = field%values( i , j , k )
          enddo
        enddo
      enddo
    endif
    success = .true.

  end function setField

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Compute the conditionnal mean
!> @details
!>  @param[in] spec_rank
!> @return .TRUE. if initialization is successfull. 
!------------------------------------------------------------------------------

!  function compute_meanCond(spec_rank,minGlob,maxGlob,rangOfValue) result (success)
  function compute_meanCond(spec_rank,minGlob,maxGlob) result (success)

    !I/O data
    integer, intent(in)                       :: spec_rank 
    real(WP),intent(in),dimension(:),optional :: minGlob
    real(WP),intent(in),dimension(:),optional :: maxGlob
    logical                                   :: success
!    logical,optional                          :: rangOfValue 

    !Local data
    integer              :: i
    integer              :: check 
    
    success = .false.

!    if (spec_rank .eq. 0) write(6,'(a)')'[INFO AVG EXPECTATION] -> compute_meanCond'
    if (.not. init ) then
       if (spec_rank .eq. 0 ) write(6,'(a)')'[ERROR] -> meanCond pas initialise'
       return
    elseif (nbCondCurrent .ne. nbCondTot  ) then
       if (spec_rank .eq. 0 ) write(6,'(a)')'[ERROR] -> Il manque une variable de conditionnement'
       return
    else
      if ( present(minGlob) .neqv. present(maxGlob)) then
         write(6,'(a)')'[ERROR] An optinal argument is missing: abort !'
         return
      endif
      if ( present(minGlob) .and. present(maxGlob)) then
        if ( size(minGlob) .ne. nbCondTot .and. size(maxGlob) .ne. nbCondTot) then
          write(6,'(a)')'[ERROR] Values Min/Max for conditional variable&
                            & disagree with number of variables: abort !'
          return
        endif  
      endif
      do i=1,nbCondTot
        if ( present(minGlob) .and. present(maxGlob) ) then
          call computeDiscCond( condAll( i ) , spec_rank ,minGlob(i), maxGlob(i))
        else
          call computeDiscCond( condAll( i ) , spec_rank )
        endif
      enddo
      call computeAllBin()
      ! Initialisation du tableau contenant les valeurs de l'EO
      allocate (valuesFieldOut(allBin),stat=check)
      if (check .ne. 0 ) then
        write(6,'(a)')'[ERROR] Not enought memory for valuesFieldOut: abort!' 
        return
      endif
      valuesFieldOut=0.0_WP
      ! Initialisation du tableau contenant le nombre de valeurs dans chaque intervalle
      allocate (nbValuesFieldOut(allBin),stat=check)
      if (check .ne. 0 ) then
        write(6,'(a)')'[ERROR] Not enought memory for nbvaluesFieldOut: abort!' 
        return
      endif
      nbValuesFieldOut=0_DI
      call computeTabExposant()
      call meanCond(spec_rank)
      !call meanCond(spec_rank,present(rangOfValue))
    endif
    success = .true.

  end function compute_meanCond


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> To get the Mean field at the end of computational 
!> @details
!>  @param[in] ChampOut is the real_data_layout where to write the mean field 
!>  @param[in] spec_rank
!------------------------------------------------------------------------------
  function getMeanField(ChampOut,spec_rank) result (success)
    
    !I/O data
    type(real_data_layout), intent(inout) :: ChampOut
    integer ,intent(in)                   :: spec_rank
    logical                               :: success

    !Local data
    integer :: i,j,k

    success = .false.

    if ( ChampOut%xmin .ne. xminLoc .or. &
       & ChampOut%xmax .ne. xmaxLoc .or. &
       & ChampOut%ymin .ne. yminLoc .or. &
       & ChampOut%ymax .ne. ymaxLoc .or. &
       & ChampOut%zmin .ne. zminLoc .or. &
       & ChampOut%zmax .ne. zmaxLoc ) then
       write(6,'(a)')'[ERROR] in getMeanField, fields are not the same'
       return
    endif  

!    if (spec_rank .eq. 0) write(6,'(a)')'[INFO AVG EXPECTATION]  get_meanField'
    do k = zminLoc,zmaxLoc 
      do j = yminLoc,ymaxLoc 
        do i = xminLoc,xmaxLoc
          ChampOut%values( i , j , k ) = fieldOut ( i , j , k )
        enddo
      enddo
    enddo

    success = .true.

  end function getMeanField

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> To get the Mean field at the end of computational 
!> @details
!>  @param[inout] ChampDiscOut contains the values of the conditional mean field 
!>  @param[inout] indiceOut contains index of conditional mean field i.e. :
!>                the value of ChampDiscOut for the nth conditional value at its 
!>                kth bin is : ChampDiscOut( indiceOut( k , n ) )
!------------------------------------------------------------------------------
  subroutine showFieldDisc(spec_rank)

    !Local data
    integer,intent(in)                :: spec_rank 
    integer                           :: i 
                  
    do i = 1,allBin
       if (spec_rank .eq. 0 ) write(6,'(a,2x,i0,5x,g20.10)') '[INFO AVG EXPECTATION] FieldDisc',i,valuesFieldOut(i)
    enddo

  end subroutine showFieldDisc

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> To get the connection between the ith value of meanField and the jth bin of
!> kth conditional value 
!> @details
!>  @param[in] spec_rank is the rank of processus in the pool of cpu 
!> 
!------------------------------------------------------------------------------
  subroutine showConnectTable(spec_rank)

    !Local data
    integer,intent(in)                :: spec_rank 
    integer                           :: i 
    integer                           :: maxbin,sizeBin 
    integer                           :: sizeAllBin 
    character(len=30)                 :: formatType                  
                  
    sizeBin = 1
    maxBin = maxval(condAll(:)%bin) 
    do while (maxBin/(10**sizeBin) .ge. 1)
        sizeBin = sizeBin + 1 
    end do
    sizeAllBin = 1
    do while (allBin/(10**sizeAllBin) .ge. 1)
        sizeAllBin = sizeAllBin + 1 
    end do
    write(formatType,'(a,i0,a,i0,a,i0,a)')'(a,2x,i',sizeAllBin,',',nbCondTot,'(2x,i',sizeBin,'))'
!    if (spec_rank .eq. 0 ) print *,formatType
    do i = 1,allBin
       if (spec_rank .eq. 0 ) write(6,formatType) '[INFO AVG EXPECTATION] ConnectTable',i,indice_jj1d_to_jj(i)
    enddo

  end subroutine showConnectTable



!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> To get the connection between the ith value of meanField and the jth bin of
!> kth conditional value 
!> @details
!>  @param[inout]  ConnectTable contains in the first column the ith value of meanField.
!>                 The kth+1 column at the ith row contains the number of bin of the kth 
!>                 conditionnal value concerned by the ith values of meanField
!>  @param[inout]  dim1 number of row of ConnectTable 
!>  @param[inout]  dim2 number of column of ConnectTable 
!>  @return        success
!------------------------------------------------------------------------------
  function getConnectTable(ConnectTable,dim1,dim2) result(success)

    !I/O data
    integer,intent(inout),dimension(:,:),allocatable  :: ConnectTable 
    integer,intent(out)                               :: dim1,dim2 
    logical                                           :: success

    !Local data
    integer                               :: i,j
    integer,dimension(nbCondTot)          :: temp 

    success=.false.

    dim1 = allBin
    dim2 = nbCondTot
    if ( size(ConnectTable(:,1)) .ne. dim1 .or. size(ConnectTable(1,:)) .ne. dim2 ) then
      write (6,'(a)') '[ERROR] Size of ConnectTable: Wrong dimensions !' 
      return
    endif
    do i = 1,allBin
      temp = indice_jj1d_to_jj( i )
      do j = 1,nbCondTot
        ConnectTable( i , j ) = temp ( j )
      enddo
    enddo

    success=.true.

  end function getConnectTable


!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> To get the Mean field at the end of computational 
!> @details
!>  @param[inout] ChampDiscOut contains the values of the conditional mean field 
!>  @param[inout] dim1 number of values of ChampDiscOut
!------------------------------------------------------------------------------
  function getFieldDisc( ChampDiscOut,dim1 ) result(success)

    !I/O data
    real(WP) , intent(inout), dimension(:)  :: ChampDiscOut
    integer , intent(out)                   :: dim1 
    logical                                 :: success 

    !Local data
    integer                                 :: i
  
    success = .false.   

    dim1 = allBin
    if ( size(ChampDiscOut(:)) .ne. dim1 ) then
      write (6,'(a)') '[ERROR] Size of ChampDiscOut: Wrong dimensions !' 
      return
    endif
    do i = 1,allBin
       ChampDiscOut( i ) = valuesFieldOut ( i )
    enddo

    success = .true.   

  end function getFieldDisc

!------------------------------------------------------------------------------
  function getCondDisc( CondDisc , num , dim1 ) result(success)

    !I/O data
    !real(WP) , intent(inout), dimension(:),allocatable    :: CondDisc
    real(WP) , intent(inout), dimension(:)       :: CondDisc
    integer , intent(in)                         :: num 
    integer , intent(out)                      :: dim1 
    logical                                      :: success 

    !Local data
    integer                                      :: i
  
    success = .false.   

    if ( num .gt. nbCondTot) then
      write (6,'(a)') '[ERROR] Conditional value called do not exist : abort !' 
      return
    endif
    dim1 = condAll(num)%bin
    if ( size(CondDisc(:)) .ne. dim1 ) then
      write (6,'(a)') '[ERROR] Size of CondDisc : Wrong dimensions !' 
      return
    endif
    if ( condAll(num)%typeDisc .eq. 1 ) then
      do i = 1,condAll(num)%bin
         CondDisc( i ) = condAll(num)%valCondDisc(i)
      enddo
    elseif ( condAll(num)%typeDisc .eq. 2 ) then
      do i = 1,condAll(num)%bin
         CondDisc( i ) = condAll(num)%valCondDisc(i)
      enddo
    endif 

    success = .true.   

  end function getCondDisc


!------------------------------------------------------------------------------
  subroutine showCondDisc( spec_rank , num ) 

    !I/O data
    integer , intent(in)                                  :: num 
    integer , intent(in)                                  :: spec_rank 

    !Local data
    integer                                               :: i
  
    if ( num .gt. nbCondTot) then
      write (6,'(a)') '[ERROR] Conditional value called do not exist : abort !' 
    endif
    if (condAll(num)%typeDisc .eq. 1 ) then
      do i = 1,condAll(num)%bin
         if (spec_rank .eq. 0 ) write(6,'(a,2x,i0,2x,i0,2x,g20.10)') '[INFO AVG EXPECTATION] CondVar', &
                    & num,i,condAll(num)%minCond + (i-0.500_WP) * condAll(num)%delta 
      enddo
    else
      write (6,'(a)') '[ERROR] a coder'
    endif

  end subroutine showCondDisc


!------------------------------------------------------------------------------
  function getPdfCond( PdfCond , num , dim1 ) result(success)

    !I/O data
   ! real(WP) , intent(inout), dimension(:),allocatable    :: PdfCond
   ! real(WP) , intent(inout), dimension(:)   :: PdfCond
    integer(DI) , intent(inout), dimension(:)  :: PdfCond
    integer , intent(in)                       :: num 
    integer , intent(out)                      :: dim1 
    logical                                    :: success 

    !Local data
    integer                                    :: i
  
    success = .false.   

    if ( num .gt. nbCondTot) then
      write (6,'(a)') '[ERROR] Conditional value called do not exist : abort !' 
      return
    endif
    dim1 = condAll(num)%bin
    if ( size(PdfCond(:)) .ne. dim1 ) then
      write (6,'(a)') '[ERROR] Size of PdfCond : Wrong dimensions !' 
      return
    endif
    do i = 1,condAll(num)%bin
       PdfCond( i ) = condAll(num)%nbValues(i) 
    enddo

    success = .true.   

  end function getPdfCond


!------------------------------------------------------------------------------
  subroutine showPdfCond(spec_rank,num)


    !I/O data
    integer , intent(in)                                  :: num 
    integer , intent(in)                                  :: spec_rank 

    !Local data
    integer                                               :: i
 
    if (.not. initPdf ) then
      write (6,'(a)') '[ERROR] You have to compute pdf before'
    endif
    if ( num .gt. nbCondTot) then
      write (6,'(a)') '[ERROR] Conditional value called do not exist : abort !' 
    endif
    do i = 1,condAll(num)%bin
      if (spec_rank .eq. 0 ) write(6,'(a,2x,i0,2x,i0,2x,i0)') '[INFO AVG EXPECTATION] pdf CondVar', &
                             & num,i,condAll(num)%nbValues(i)
    enddo

  end subroutine showPdfCond

!------------------------------------------------------------------------------
  subroutine showPdfConjointe(spec_rank)

    !I/O data
    integer , intent(in)  :: spec_rank 

    !Local data
    integer               :: i
    integer(DI)           :: maxPdf,sizePdf 
    integer               :: sizeAllBin 
    character(len=30)     :: formatType
  
    maxPdf = maxval(nbValuesFieldOut)
    sizePdf = 1_DI
    sizeAllBin = 1
    do while (maxPdf/(10**sizePdf) .ge. 1)
        sizePdf = sizePdf + 1 
    end do
    do while (allbin/(10**sizeAllBin) .ge. 1)
        sizeAllBin = sizeAllBin + 1 
    end do
    write(formatType,'(a,i0,a,i0,a)')'(a,2x,i',sizeAllBin,',2x,i',sizePdf,')'
    do i = 1,allbin
      if (spec_rank .eq. 0 ) write(6,formatType) '[INFO AVG EXPECATION] Pdf Conjointe', &
                             & i,nbValuesFieldOut( i )
    enddo

  end subroutine showPdfConjointe

!------------------------------------------------------------------------------
  function getPdfConjointe( PdfConjointe , dim1 ) result(success)

    !I/O data
    integer(DI) , intent(inout), dimension(:)   :: PdfConjointe
    integer , intent(out)                       :: dim1 
    logical                                     :: success 

    !Local data
    integer                                     :: i
  
    success = .false.   

    dim1 = allBin 
    if ( size(PdfConjointe(:)) .ne. dim1 ) then
      write (6,'(a)') '[ERROR] Size of PdfConjointe : Wrong dimensions !' 
      return
    endif
    do i = 1,allBin
!       PdfConjointe( i ) = idint(aint (nbValuesFieldOut(i))) 
       PdfConjointe( i ) = nbValuesFieldOut(i) 
    enddo

    success = .true.   

  end function getPdfConjointe

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Compute pdf of conditional variables
!> @details
!>  @param[in] cond is the condition 
!------------------------------------------------------------------------------
  subroutine computePdfConds()

    !Local Data
    integer                       :: i,j 
    integer,dimension(nbCondTot)  :: tab_indice

    do i=1,allBin
      tab_indice = indice_jj1d_to_jj( i )
      do j=1,nbCondTot
if (tab_indice(j) .lt. 1) print *,'ERROR 3'
        condAll(j)%nbValues( tab_indice(j) ) = condAll(j)%nbValues( tab_indice(j) )&
                                             &  + 1_DI 
      enddo
    enddo
    initPdf = .true.

  end subroutine computePdfConds

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                        PRIVATE ROUTINES                                                                  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Compute the discretization of conditional variables 
!> @details
!>  @param[in] cond is the condition 
!------------------------------------------------------------------------------
  subroutine computeDiscCond( cond , spec_rank ,minGlob, maxGlob)
    
    !I/O
    type(varCond), intent(inout) :: cond
    integer, intent(in)          :: spec_rank 
    real(WP),intent(in),optional :: minGlob,maxGlob

    !Local Data
    integer  :: i

    if (cond%typeDisc .eq. 1) then
      if ( present(minGlob) .and. present(maxGlob) ) then
        call deltaCondCst( cond , spec_rank , minGlob, maxGlob )
      else
        call deltaCondCst ( cond , spec_rank )
      endif
!      if (spec_rank .eq. 0) write(6,'(5(a,1x,g15.5,1x))')'[INFO AVG EXPECTATIONS] Discretization type: ',cond%typeDisc,&
!                         'sampling: ',cond%bin,'min of condition: ',cond%minCond,&
!                         &'max of condition: ',cond%maxCond,'delta: ',cond%delta
      cond%deltaTot = cond%maxCond - cond%minCond
      do i=1,cond%bin
        cond%valCondDisc(i) = cond%minCond + ( i - 0.5_WP ) * cond%deltaTot / real(cond%bin,WP)
      enddo
    elseif (cond%typeDisc .eq. 2) then
      if ( present(minGlob) .and. present(maxGlob) ) then
        call deltaCondVariable( cond , spec_rank , minGlob, maxGlob )
      else
        call deltaCondVariable( cond , spec_rank )
      endif
!      if (spec_rank .eq. 0) write(6,'(5(a,1x,g15.5,1x))')'[INFO AVG EXPECTATIONS] Discretization type: ',cond%typeDisc,&
!                         'sampling: ',cond%bin,'min of condition: ',cond%minCond,&
!                         &'max of condition: ',cond%maxCond,'delta: ',cond%delta
    endif

  end subroutine computeDiscCond

!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant, LEGI
!> Computation of the conditional mean 
!> @details
!>  @param[in] ChampDiscOut
!>  @param[in] indiceOut 
!------------------------------------------------------------------------------
!  subroutine meanCond( spec_rank , rangOfValue )
  subroutine meanCond( spec_rank)

    !I/O
    integer, intent(in)                    :: spec_rank 

    !Local data
    integer                                :: i,j,k,l,jj1d
    integer                                :: check
    integer,dimension(:,:,:),allocatable   :: indic
    integer , dimension( nbCondTot )       :: jj
    real(WP),dimension(:), allocatable     :: valueloc
    integer(DI),dimension(:), allocatable  :: nombreloc
!    logical                                :: rangOfValue
!    logical                                :: outofbound


    allocate( valueloc(allBin) , stat=check)
    if ( check .ne. 0) then
      write (6,'(a)') '[ERROR] Not able to allocate valueloc: Not enough memory' 
    endif
    allocate( nombreloc(allBin) , stat=check)
    if ( check .ne. 0) then
      write (6,'(a)') '[ERROR] Not able to allocate nombreloc: Not enough memory' 
    endif
    allocate( indic(xminLoc:xmaxLoc,yminLoc:ymaxLoc,zminLoc:zmaxLoc) , stat=check)
    if ( check .ne. 0) then
      write (6,'(a)') '[ERROR] Not able to allocate indic: Not enough memory' 
    endif

    valueloc= 0.0_WP
    nombreloc= 0_DI

    if (spec_rank .eq. 0) then
!      write (6,'(a)') '[INFO AVG EXPECTATION] Computing conditional mean ...' 
!      if (interSousDim) then
!        write (6,'(a)') '[INFO AVG EXPECTATION] attention les valeurs hors bornes ne sont pas prise en compte !' 
!      endif
    endif
    if (.not. interSousDim ) then
    do k = zminLoc,zmaxLoc 
      do j = yminLoc,ymaxLoc 
        do i = xminLoc,xmaxLoc
          do l = 1,nbCondTot
            jj(l) = get_indice( condAll( l ) ,i,j,k)
!            if ( spec_rank .eq. 0 .and. i.eq.xminLoc .and. j.eq.yminLoc .and. k.eq.zminLoc) write (6,'(a,1x,i0)') &
!                          '[INFO AVG EXPECTATION] getting indices with typeDisc',condAll(l)%typeDisc 
            if (jj(l).lt. 1) print *,'ERROR 4'
            condAll(l)%nbValues( jj(l) ) = condAll(l)%nbValues( jj(l) ) + 1_DI
          enddo
          !transformer le vect des indices jj en un indice jj1d
          jj1d = indice_jj_to_jj1d( jj )
          !sommation des valeurs du champ objectif en vue de calculer la moyenne
          valueloc( jj1d ) = valueloc( jj1d ) + fieldIn( i,j,k )
          !incrémentation du compteur de moyenne
          nombreloc( jj1d )=nombreloc( jj1d ) + 1_DI
          !sauvegarde du numéro du bin concerné par le noeud (i,j,k)
          indic(i,j,k) = jj1d
        enddo
      enddo
    enddo
    else
    do k = zminLoc,zmaxLoc 
      do j = yminLoc,ymaxLoc 
        do i = xminLoc,xmaxLoc
!           !calcul qui n'intègre pas les valeurs hors borne dans les bin min et max
!           if (rangOfValue) then
!            outofbound = .false.
!            do l = 1,nbCondTot
!              jj(l) = get_indice( condAll( l ) ,i,j,k)
!              if ( (jj(l) .lt. 1) .or. (jj(l) .gt. condAll( l )%bin)) then 
!                outofbound = .true.
!              endif
!            enddo
!            if (.not. outofbound) then
!              do l = 1,nbCondTot
!                jj(l) = get_indice( condAll( l ) ,i,j,k)
!                condAll(l)%nbValues( jj(l) ) = condAll(l)%nbValues( jj(l) ) + 1
!              enddo
!              jj1d = indice_jj_to_jj1d( jj )
!              valueloc( jj1d ) = valueloc( jj1d ) + fieldIn( i,j,k )
!              nombreloc( jj1d )=nombreloc( jj1d ) + 1_DI
!              indic(i,j,k) = jj1d
!            endif
!          else
            do l = 1,nbCondTot
              jj(l) = get_indice( condAll( l ) ,i,j,k)
              if ( jj(l) .lt. 1 ) jj(l) = 1  
              if ( jj(l) .gt. condAll( l )%bin ) jj(l) = condAll(l)%bin
if (jj(l).lt. 1) print *,'ERROR 5'
              condAll(l)%nbValues( jj(l) ) = condAll(l)%nbValues( jj(l) ) + 1
            enddo
            jj1d = indice_jj_to_jj1d( jj )
            valueloc( jj1d ) = valueloc( jj1d ) + fieldIn( i,j,k )
            nombreloc( jj1d )=nombreloc( jj1d ) + 1_DI
            indic(i,j,k) = jj1d
!          endif
        enddo
      enddo
    enddo
    endif

    if (spec_rank .eq. 0) then
!      write (6,'(a)') '[INFO AVG EXPECTATION] Computing averages over all bin...' 
    endif
    do i=1,allBin
      !Concatenation des resultats de tous les procs
      valuesFieldOut(i)= DoGlobalSum(spec_rank,valueloc(i))
      nbValuesFieldOut(i)= DoGlobalSum(spec_rank,nombreloc(i))
      !Calcul des moyennes
      if ( nbValuesFieldOut(i) .ne. 0_DI ) then
        valuesFieldOut(i) = valuesFieldOut(i) / real(nbValuesFieldOut(i),WP)
      endif
      !if (rank .eq. 0) write(6,'(a,1x,f10.3)')'[INFO AVG EXPECTATION]->valuesFieldOut(i):',valuesFieldOut(i)
    enddo
    !Attribution des valeurs de moyenne cond aux points de l'espace
    do k = zminLoc,zmaxLoc 
      do j = yminLoc,ymaxLoc 
        do i = xminLoc,xmaxLoc
          fieldOut(i,j,k) = valuesFieldOut( indic(i,j,k) )
        enddo
      enddo
    enddo
    deallocate( valueloc ) 
    deallocate( nombreloc )
    deallocate( indic )

  end subroutine meanCond

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "computeAllBin"
!>
!> @return  Calcul le nombre total d'intervalle d'echantillonnage et le stock
!>          dans la variable globale allBin
!!------------------------------------------------------------------------------
  subroutine computeAllBin()

    !Local Data      
    integer :: temp
    integer :: i

    temp = 1
    do i=1,nbCondTot
      temp  = temp * condAll(i)%bin
    enddo
    allBin = temp
    !if (rank .eq. 0 ) write(6,'(a8,1x,i15)')'allBin =',allBin

  end subroutine computeAllBin


!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "indice_jj_to_jj1d"
!>
!> @param [in]      jj :Vecteur de taille nbCondTot. Les indices de ce
!>            vecteur representent les coordonnees cartesiennes dans l'espace 
!>            d'ordre nbConTot de l'intervalle sur lequel est moyenne le champs
!> @return      le cardinal de l'intervalle considere dans une dimension
!!------------------------------------------------------------------------------
  function indice_jj_to_jj1d( jj ) result(indice)

    !I/O
    integer,dimension( nbCondTot ), intent(in) :: jj
    integer                                    :: indice

    !Local data
    integer                                     :: i

    indice = 0 
    do i =1,nbCondTot
      indice = indice + ( jj( i ) - 1 ) * tabExposant( i )
    enddo
    indice = indice + 1

  end function indice_jj_to_jj1d

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "indice_jj1d_to_jj"
!>
!> @param [in]      valjj1d : Le cardinal de l'intervalle considere dans une dimension
!> @return      Vecteur de taille nbCondTot. Les indices de ce
!>            vecteur representent les coordonnees cartesiennes dans l'espace 
!>            d'ordre nbConTot de l'intervalle sur lequel est moyenne le champ
!!------------------------------------------------------------------------------
  function indice_jj1d_to_jj( valjj1d ) result(tab_indice)

    !I/O
    integer ,intent(in)           :: valjj1d
    integer ,dimension(nbCondTot) :: tab_indice

    !Local Data
    integer                       :: i
    integer                       :: tempCible

    tab_indice = 0
    tempCible = valjj1d-1
 
    do i = 0,nbCondTot-1
      tab_indice( nbCondTot-i ) = tempCible / tabExposant( nbCondTot-i )
      tempCible = tempCible - tab_indice( nbCondTot-i ) * tabExposant( nbCondTot-i )
      tab_indice( nbCondTot-i ) = tab_indice( nbCondTot-i ) + 1
    enddo

  end function indice_jj1d_to_jj

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "computeTabExposant"
!>
!> @return  Calculs les coefficients servant a la decomposition/composition des
!>          intervalles de moyennes conditionnelles et les stocke dans le vecteur
!>          global tabExposant
!!------------------------------------------------------------------------------
  subroutine computeTabExposant()

    !Local data
    integer :: i

    tabExposant( 1 ) = 1
!    write(6,'(a12,i1,a2,i5)'),'tabExposant(',1,')=',tabExposant( 1 )
    do i=2,nbCondTot
      tabExposant( i ) = tabExposant( i - 1 ) * condAll( i-1 )%bin
!      write(6,'(a12,i1,a2,i5)'),'tabExposant(',i,')=',tabExposant( i )
    enddo

  end subroutine computeTabExposant

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "get_indice"
!>
!> @param   cond : structure de type varCond contenant une variable de conditionnement
!> @param   i : indice spatial
!> @param   j : indice spatial
!> @param   k : indice spatial
!> @return  l'indice de l'intervalle dans lequel se trouve la variable de conditionnement
!>          au point de l'espace (i,j,k) 
!!------------------------------------------------------------------------------
  function get_indice(cond,i,j,k) result(indice)

   use sortfind_tools , only : seekDichoto

    !I/O data
    type( varCond ), intent(in) :: cond
    integer, intent(in)         :: i,j,k
    integer                     :: indice

    indice = -1
    if (cond%typeDisc .eq. 1) then
      indice = indiceCondCst(cond%minCond , cond%valCond(i,j,k), cond%deltaTot ,cond%bin )
      !indice = seekDichoto(cond%valCond(i,j,k),cond%valCondDisc,1,cond%bin)
    elseif (cond%typeDisc .eq. 2) then
      indice = seekDichoto(cond%valCond(i,j,k),cond%valCondDisc,1,cond%bin)
    else 
      
    endif
    
  end function get_indice

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "indiceCondCst"
!>
!> @param [in]  phiMin : La plus petite valeur de la variable de conditionnement
!> @param [in]  phiCible : valeur dont on cherche la position dans le partionnement
!> @param [in]  deltaTot : plage de variation de la variable de conditionnement 
!> @param [in]  bin : nombre d'intervalles
!> @return  indice de l'intervalle dans lequel se trouve la variable de conditionnement
!>          dans le cas d'un decoupage constant de la variable de conditionnnement 
!!------------------------------------------------------------------------------
  function indiceCondCst ( phiMin , phiCible , deltaTot , bin ) result(indice)

    !I/O data
    real(WP),intent(in) :: phiMin
    real(WP),intent(in) :: phiCible
    real(WP),intent(in) :: deltaTot
    integer,intent(in)  :: bin
    integer             :: indice
    !Local data
    real(WP)            :: tx,jjfloat

!print *,'[INFO AVG EXPECTATION] phiMin=',phiMin
!print *,'[INFO AVG EXPECTATION] phiCible=',phiCible

    tx=phiCible-phiMin
    tx=tx/deltaTot
    jjfloat=tx*bin
    indice = ceiling( jjfloat )
if (indice .lt. 1 ) print *,'WARNING in indiceCondCst phiMin->',phiMin,' phiCible->',phiCible,' deltaTot->',deltaTot,' bin->',bin

  end function indiceCondCst


!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "deltaCondCst"
!>
!> @param [in]  cond : Valeurs de la variable de conditionnement 
!> @param [in]  bin : Nombre d'intervalles d'echantillonnement
!> @param [out]  minCond : La plus petite valeur de la variable de conditionnement 
!> @param [out]  maxCond : La plus grande valeur de la variable de conditionnement 
!> @param [out]  delta : intervalle elementaire d'echantillonnement
!!------------------------------------------------------------------------------
  subroutine deltaCondVariable( cond , spec_rank , minGlob, maxGlob )
     
    use sortfind_tools , only : parallelSort
    use parallel_tools , only : DoBcastScal
    use mpilayout_tools , only : getnbcpus

    !I/O
    type(varCond),  intent(inout)       :: cond
    integer , intent(in)                :: spec_rank
    real(WP) , intent(in),optional      :: minGlob
    real(WP) , intent(in),optional      :: maxGlob

    real(WP),dimension(:),allocatable   :: tabToSort,tabSorted
    integer      :: GridPts,processor,ptsParEchan,indiceLoc
    integer      :: i,j,k,l,ires
    !logical      :: existe
    !character(len=64)  :: file_name,numEtape
    
    if (present(minGlob) .and. present (maxGlob)) then
    endif

    gridPts=nxLoc*nyLoc*nzLoc
    allocate(tabToSort(gridPts),tabSorted(gridPts),stat=ires)
    if (ires.ne.0) return 
    l = 1
    do k = zminLoc,zmaxLoc
      do j = yminLoc,ymaxLoc
        do i = xminLoc,xmaxLoc
          tabToSort(l)=cond%valCond(i,j,k)
          l = l + 1
        enddo
      enddo
    enddo
    if (spec_rank .eq. 0) write(6,'(a)') '[INFO AVG EXPECTATION] Parallel sort in progress...'
    if (.not. parallelSort(tabToSort,tabSorted,spec_rank)) then 
      write(6,'(a)') '[ERROR AVG EXPECTATION] parallel sort failed' 
      return
    endif
    ptsParEchan=floor(real(nxGlob,WP)*real(nyGlob,WP)*real(nzGlob,WP)/cond%bin)
    do i=1,cond%bin
      processor = 0
      if (ptsParEchan*i-ptsParEchan/2 .gt. gridPts) then
        indiceLoc = ptsParEchan*i-ptsParEchan/2-gridPts
        processor = processor + 1
        do while (indiceLoc .gt. gridPts)
          processor = processor + 1
          indiceLoc = indiceLoc - gridPts
        enddo
      else
        indiceLoc = ptsParEchan*i-ptsParEchan/2
      endif
      if (indiceLoc==0) indiceLoc=1
      if (spec_rank .eq. processor) cond%valCondDisc(i)=tabSorted(indiceLoc)
      if (.not. DoBcastScal(cond%valCondDisc(i),processor,spec_rank)) return
!      if (spec_rank .eq. 0) then 
!         inquire( file='param_indices.table', exist=existe)
!         open(unit=18, file='param_indices.table', form="FORMATTED", position='append')
!         if (.not. existe) write(18,'(a)')'#cond%bin,i,indiceLoc,processor,ptsParEchan,ptsParEchan*i-ptsParEchan/2'
!         write (18,'(6(1x,i0),1x,g15.8)') cond%bin,i,indiceLoc,processor,ptsParEchan,ptsParEchan*i-ptsParEchan/2,cond%valCondDisc(i)
!        close(18)
!      endif
    enddo
    deallocate(tabSorted)
    deallocate(tabToSort)
    

  end subroutine deltaCondVariable

!!------------------------------------------------------------------------------
!> @author 
!> Antoine Vollant
!>
!> @details de "deltaCondCst"
!>
!> @param [in]  cond : Valeurs de la variable de conditionnement 
!> @param [in]  bin : Nombre d'intervalles d'echantillonnement
!> @param [out]  minCond : La plus petite valeur de la variable de conditionnement 
!> @param [out]  maxCond : La plus grande valeur de la variable de conditionnement 
!> @param [out]  delta : intervalle elementaire d'echantillonnement
!!------------------------------------------------------------------------------
  subroutine deltaCondCst( cond , spec_rank , minGlob, maxGlob )

    !I/O
    type(varCond),  intent(inout)       :: cond
    integer , intent(in)                :: spec_rank
    real(WP) , intent(in),optional      :: minGlob
    real(WP) , intent(in),optional      :: maxGlob

    !Local data
    real(WP)                            :: correction
    real(WP)                            :: tolerance !defini le demi taux d'agrandissement de l'intervalle
    real(WP)                            :: minCondLoc
    real(WP)                            :: maxCondLoc
    real(WP)                            :: temp

    tolerance = 0.0001_WP

    !Recherche du min et du max sur tous les procs
    maxCondLoc = DoGlobalMax(spec_rank,maxval(cond%valCond))
    minCondLoc = DoGlobalMin(spec_rank,minval(cond%valCond))

    if ( present (minGlob) ) then
      if ( minGlob .gt. minCondLoc .or. maxGlob .lt. maxCondLoc ) then
        if (spec_rank .eq. 0) write(6,'(a)')'[WARNING] Beware the min/max choosen is &
        &between the bounds of conditional variable'
        interSousDim = .true.
      endif
      maxCondLoc = maxGlob 
      minCondLoc = minGlob
    endif

    !Calcul de l'amplitude de la condition Cond avec une correction
    temp = maxCondLoc-minCondLoc
    correction = tolerance*temp
    temp = temp+2.0_WP*correction
    cond%deltaTot = temp 
    !Correction des condMin et condMax
    cond%maxCond=maxCondLoc+correction
    cond%minCond=minCondLoc-correction
    !cond%delta=(maxCondLoc-minCondLoc)/REAL(cond%bin,WP)
    cond%delta=(cond%maxCond-cond%minCond)/REAL(cond%bin,WP)

  end subroutine deltaCondCst

!-------------------------------------------------------------------------------------------
!> Compute conditional expectations 
!! @author Antoine Volllant 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    Var             = array of real_data_layout conditonal variables 
!!    @param[in]    TypeDisc        = 1 for constant decomposition of conditional variables 
!!    @param[in]    bin             = values for number of samples. if TypeDisc = 1 it's the number of intervals 
!!    @param[in]    minVar          = force to set the min of a conditional variable 
!!    @param[in]    maxVar          = force to set the max of a conditional variable 
!!    @param[inout] estimOpt        = contains the tabulation of optimal estimator (optional)
!!                                    size of array are estimOpt(nbCond+1,bin) 
!!    @param[inout] fieldOpt        = real_data_layout with optimal estimator applied to the field (optional)
!!    @param[in]    pdfCond         = gives the pdf of each conditional variable 
!!    @param[in]    pdfField        = gives the pdf of optimal estiamtor field 
!!    @return       success         = logical value 
!-------------------------------------------------------------------------------------------
! function computeEONVar( spec_rank,Var,Field,typeDisc,bin,minVar,maxVar,estimOpt,fieldOpt,&
!                       & pdfCond,pdfField,rangeOfVal ) result(success)
 function computeEONVar( spec_rank,Var,Field,typeDisc,bin,minVar,maxVar,estimOpt,fieldOpt,&
                       & pdfCond,pdfField ) result(success)

    !I/O data
    type(real_data_layout),intent(in),dimension(:)    :: Var
    type(real_data_layout),intent(in)                 :: Field
    integer,intent(in)                                :: spec_rank
    integer,intent(in),dimension(:)                   :: bin
    integer,intent(in),dimension(:)                   :: typeDisc 
    real(WP),intent(in),dimension(:),optional         :: minVar,maxVar
    real(WP),intent(inout),dimension(:,:),optional    :: estimOpt 
    type(real_data_layout),intent(inout),optional     :: fieldOpt
    integer(DI),intent(inout),dimension(:,:),optional :: pdfCond
    integer(DI),intent(inout),dimension(:),optional   :: pdfField
!    logical,optional                                  :: rangeOfVal
    logical                                           :: success

    !Local data 
    integer,dimension(:,:),allocatable     :: ConnectTable
    real(WP),dimension(:,:),allocatable    :: condDisc
    real(WP),dimension(:),allocatable      :: valDiscField
    integer                                :: nbcond
    integer                                :: dim1,dim2
    integer                                :: i,j 
    integer                                :: binTot,check

    success = .false.
    !Verification de la taillle des tableaux
    nbcond = size(Var)

    if ( nbcond .lt. 1 ) then
      write (6,'(a)')'[ERROR] The number of conditional variables is 0 : abort!' 
      return
    endif 
    if (present(minVar) .neqv. present(maxVar)) then
      write (6,'(a)')'[ERROR] Min or Max value is missing : abort!' 
      return
    endif 
    if ( size(bin) .ne. nbcond ) then
      write (6,'(a)')'[ERROR] Wrong size of vect bin : abort!' 
      return
    endif 
    if ( size(typeDisc) .ne. nbcond ) then
      write (6,'(a)')'[ERROR] Wrong size of vect typeDisc : abort!' 
      return
    endif 
    if ( present(minVar) .and. present(maxVar) ) then
       if (size(minVar) .ne. nbcond .and. size(maxVar) .ne. nbcond ) then
          write (6,'(a)')'[ERROR] Wrong size of vect min or max : abort!' 
          return
       endif 
    endif 
    if ( present(estimOpt) ) then 
       if (size(estimOpt,1) .ne. nbcond+1 ) then
          write (6,'(a)')'[ERROR] Wrong size of vec estimOpt(:,1) : abort!' 
          return
       end if
    endif 
    if ( present(pdfCond)) then
       if (size(pdfCond(:,1)) .ne. nbcond ) then
          write (6,'(a)')'[ERROR] Wrong size of vect pdfCond(:,1) : abort!' 
          return
      endif
    endif
    binTot = 1
    do i= 1,nbcond
      binTot = binTot*bin(i)
      if ( present(pdfCond) ) then
         if ( (size(pdfCond(i,:)) .lt. bin(i)) ) then
            write (6,'(a,i0,a)')'[ERROR] Wrong size of vect pdfCond(',i,',:) : abort!' 
            return
         endif
      endif
    enddo
    if ( present(estimOpt)) then
       if (size(estimOpt,2) .ne. binTot ) then
          write (6,'(a)')'[ERROR] Wrong size of vec estimOpt(1,:) : abort!' 
          return
       endif 
    endif 
    if ( present(pdfField)) then 
       if (size(pdfField) .ne. binTot ) then
          write (6,'(a)')'[ERROR] Wrong size of vect pdfField : abort!' 
          return
       endif
    endif


    !allocation des tableaux temporaires
    allocate(ConnectTable(binTot,nbcond),stat=check)
    if ( check .ne. 0 ) then
      write (6,'(a)')'[ERROR] Not able to allocate ConnectTable : abort!'
      return
    endif
    allocate(condDisc( nbcond ,maxval(bin)),stat=check)
    if ( check .ne. 0 ) then
      write (6,'(a)')'[ERROR] Not able to allocate condDisc : abort!'
      return
    endif
    allocate(valDiscField(binTot),stat=check)
    if ( check .ne. 0 ) then
      write (6,'(a)')'[ERROR] Not able to allocate valDiscField : abort!'
      return
    endif
    !Init Moyenne conditionnelle
    if (.not. init_meanCond( nbcond ,spec_rank )) then
      write(6,'(a)')'[ERROR] Not able to initialize meanCond : abort !!'   
      return
    endif
    if (.not. setField( Field ) ) then
      write(6,'(a)')'[ERROR] Not able to setFields : abort !!'   
      return
    endif
    do i=1,nbcond
      if (.not. setConds( Var(i) , typeDisc(i) , spec_rank , bin(i) ) ) then
        write(6,'(a)')'[ERROR] Not able to setConds 1 : abort !!'   
        return
      endif
    enddo

    !Calcul de la moyenne conditionnelle
    if ( present(minVar) .and. present(maxVar)) then
!      if (present (rangeOfVal)) then
!        if (.not. compute_meanCond(spec_rank,minVar,maxVar,rangeOfVal)) then
!          write(6,'(a)')'[ERROR] Not able to initialize conds : abort !!'   
!          return
!        endif
!      else
        if (.not. compute_meanCond(spec_rank,minVar,maxVar)) then
          write(6,'(a)')'[ERROR] Not able to initialize conds : abort !!'   
          return
        endif
!      endif
    else
      if (.not. compute_meanCond(spec_rank)) then
        write(6,'(a)')'[ERROR] Not able to initialize conds : abort !!'   
        return
      endif
    endif

    !Récupération de la tabulation de l'EO
    if ( present(estimOpt) ) then
      if (.not. getConnectTable(ConnectTable,dim1,dim2)) then
        write(6,'(a)')'[ERROR] Not able to initialize conds : abort !'   
        return
      endif
      do i=1,nbcond
        if (.not. getCondDisc( condDisc(i,:),i,dim1 )) then
          write(6,'(a,i0,a)')'[ERROR] Not able to get conds ',i,' : abort !'   
          return
        endif
      enddo
      if (.not. getFieldDisc( valDiscField,dim1 )) then
        write(6,'(a)')'[ERROR] Not able to initialize conds : abort !'   
        return
      endif
      do i=1,binTot
        do j=1,nbcond
          estimOpt(j,i) = condDisc(j,ConnectTable(i,j))
        enddo
        estimOpt(nbcond+1,i) = valDiscField(i)
      enddo
    endif


    !Récupération l'EO appliqué au domaine de calcul
    if ( present(fieldOpt) ) then
      if (.not. getMeanField(fieldOpt,spec_rank)) then
        write(6,'(a)')'[ERROR] Not able to get fielOpt : abort !'   
        return
      endif
    endif
    !Récupération des pdf de chacune des variables conditionnelles 
    if ( present(pdfcond) ) then
      call computePdfConds
      do i=1,nbcond
        if (.not. getPdfCond(pdfcond(i,:) , i , dim1 )) then
          write(6,'(a,i0,a)')'[ERROR] Not able to get pdfcond ',i,' : abort !'   
          return
        endif
      enddo
    endif


    !Récupération de la pdf conjointe de toutes les variables
    if ( present(pdfField) ) then
      if (.not. getPdfConjointe( pdfField , dim1 )) then
        write(6,'(a)')'[ERROR] Not able to get pdfField : abort !'   
        return
      endif
    endif

    deallocate(ConnectTable)
    deallocate(condDisc)
    deallocate(valDiscField)
    call finalize_meanCond


    success = .true.

 end function computeEONVar

!-------------------------------------------------------------------------------------------
!> Update mean over several iteration 
!! @author Antoine Volllant 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[inout] valTot          = array of mean values until t-1 at the input and 
!!                                    mean values until t at output
!!    @param[inout] nbvalTot        = array of number of values to compute mean until t-1 at 
!!                                    the input  and number of values until t at output
!!    @param[in]    valLoc          = array of mean values at t 
!!    @param[in]    nbvalLoc        = array of number of values to compute mean at t 
!!    @param[in]    bin             = size of arrays 
!!    @return       success         = logical value 
!-------------------------------------------------------------------------------------------
 function computeSlideAvg(valTot,nbValtot,valLoc,nbValLoc,bin) result(success)


   !I/O data
   real(WP),intent(inout), dimension(:)   :: valTot
   integer(DI),intent(inout), dimension(:):: nbValtot
   real(WP),intent(in), dimension(:)      :: valLoc 
   integer(DI),intent(in), dimension(:)   :: nbValLoc
   integer,intent(in)                     :: bin
   logical                                :: success

   !Local data
   integer                                :: j

   success = .false.

   if ( size(valTot(:)) .lt. bin ) then
     write (6,'(a)') '[ERROR] Size of valTot : Wrong dimensions !' 
     return
   endif
   if ( size(nbValTot(:)) .lt. bin ) then
     write (6,'(a)') '[ERROR] Size of nbValTot : Wrong dimensions !' 
     return
   endif
   if ( size(valLoc(:)) .lt. bin ) then
     write (6,'(a)') '[ERROR] Size of valLoc : Wrong dimensions !' 
     return
   endif
   if ( size(nbValLoc(:)) .lt. bin ) then
     write (6,'(a)') '[ERROR] Size of nbvalLoc : Wrong dimensions !' 
     return
   endif

   do j=1,bin
     if (nbValLoc(j) .ne. 0 ) then
       valTot(j)   = real(nbValtot(j),WP)*valTot(j)+real(nbValLoc(j),WP)*valLoc(j)
       nbValTot(j) = nbValTot(j) + nbValLoc(j)
       valTot(j)   = valTot(j) /real(nbValTot(j),WP)
     endif 
   enddo

   success = .true.

 end function computeSlideAvg 

!-------------------------------------------------------------------------------------------
!> Read an ASCII file which contaains an optimal estimator 
!! @author Antoine Volllant 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[inout] valTot          = array of mean values until t-1 at the input and 
!!                                    mean values until t at output
!!    @param[inout] nbvalTot        = array of number of values to compute mean until t-1 at 
!!                                    the input  and number of values until t at output
!!    @param[in]    valLoc          = array of mean values at t 
!!    @param[in]    nbvalLoc        = array of number of values to compute mean at t 
!!    @param[in]    bin             = size of arrays 
!!    @return       success         = logical value 
!-------------------------------------------------------------------------------------------
 function computeNbLine(ficName,nbLine) result(success)

    !USE fileio

    !I/O data
    character(len=*), intent(in)   :: ficName 
    integer,intent(inout)          :: nbLine 
    logical                        :: success 

    !Local Data
    integer             :: nlines,iunit,ierr
    character(len=4960) :: buffer
    logical             :: check     

    success=.false.
    iunit=15
    ! Open the file
    ierr = 0
    inquire( file=ficName, exist=check)
    if ( check ) then
      open (iunit,file=ficName,form='formatted',status='old',iostat=ierr)
      if (ierr .NE. 0) then 
          write(6,'(a)')"Parser : unable to find file ["//TRIM(ficName)//"]"
          return
      endif 
    else
        write(6,'(a)')"the file "//TRIM(ficName)//" do not exist: abort!"
        return
    endif 

    ! Count the number of lines in the file
    ierr = 0 
    nlines = 0 
    do while(ierr .eq. 0)
       read(iunit,'(a)',iostat=ierr) buffer
       nlines = nlines + 1 
    enddo 
    !write (6,'(a,1x,i0)') 'FIN DU COMPTAGE',nlines-1
    nbLine = nlines-1
    close(iunit)

    success=.true.

 end function computeNbLine 

 
!-------------------------------------------------------------------------------------------
!> Read an ASCII file which contains an optimal estimator 
!! @author Antoine Volllant 
!!    @param[in]    spec_rank       = mpi rank inside spectral communicator
!!    @param[in]    ficName         = array of mean values until t-1 at the input and 
!!                                    mean values until t at output
!!    @param[inout] modelValues     = value table filled by the routine. its size are :
!!                                    number of lines equal to number of colums in file and the
!!                                    number of colums equal to the number of lines in the file
!!                                    number of lines = Nb conditional variable + 2
!!                                    number of columns = Nb of samples
!!    @return       success         = logical value 
!-------------------------------------------------------------------------------------------
 function readModelOpt(ficName,modelValues) result(success)

    !USE fileio

    !I/O data
    real(WP), dimension(:,:),intent(inout)  :: modelValues 
    character(len=*), intent(in)            :: ficName 
!    integer,intent(inout)                   :: modelNbVal
    logical                                 :: success 

    !Local Data
    integer             :: iunit,ierr,i
    logical             :: check     
    integer             :: spec_rank,sizevar
    integer             :: sizevalues

    spec_rank=getnbmycpu()
    sizevar=size(modelValues,1)
    sizevalues=size(modelValues,2)
    success=.false.
    iunit=66
    ! Open the file
    ierr = 0
    inquire( file=ficName, exist=check)
    if ( check ) then
       open (iunit,file=ficName,form='formatted',status='old',&
            &action='read',position='rewind',iostat=ierr)
      if (ierr .NE. 0) then 
          write(6,'(a)')"Parser : unable to find EO model file ["//TRIM(ficName)//"]"
          return
      endif 
    else
        write(6,'(a)')"the file "//TRIM(ficName)//" do not exist: abort!"
        return
    endif 
    i=1
    ierr = 0 
    do while (ierr .eq. 0 .and. i .le. sizevalues)
      read(iunit,*,iostat=ierr) modelValues(:,i)
      if (ierr.ne.0) then
        write(6,'(a,1x,i0)') '[ERROR] Not able to read File. Stop at line:',i
        return 
      endif
      i = i + 1
    enddo
    ! Close de file
    close(iunit)
!    if (spec_rank .eq. 0 ) then
!      print *,'sizevar=',sizevar,'sizevalues=',sizevalues
!      !do i= 1 , sizevalues , 1000
!      do i= 1 , sizevalues
!        print *,'[FICHIER IN]',modelValues(1,i),modelValues(2,i),modelValues(3,i)
!      enddo
!    endif

    success=.true.

 end function readModelOpt



!-------------------------------------------------------------------------------------------
!> Re-build Optimal estimator from an extra-file 
!! files 
!! @author Antoine Volllant 
!!    @param[in]    inVar  is an array of real_data_layout to look for in tabEO. The size of
!!                  array is equal at the number of conditional variable 
!!    @param[in]    tabEO  contains the optimal estimator. first dimension contains input variables + output and
!!                  second dimension contains number of points 
!!    @param[in]    EOField is the optimal estimator applied to a real_data_layout 
!!    @return       .true. if everything is OK and .false. otherwise
!-------------------------------------------------------------------------------------------
 function buildEOFieldFromEO(inVarCond,tabEO,EOFieldOut) result(success)

   use sortfind_tools , only : seekDichoto
   use parallel_tools , only : DoLocalSum
   use mpilayout_tools , only : getnbmycpu

  !I/O data
  type(real_data_layout),intent(inout)             :: EOFieldOut
  type(real_data_layout),dimension(:),intent(in)   :: inVarCond
  real(WP),dimension(:,:),intent(in)               :: tabEO
  logical                                          :: success

  !Local data
  integer            :: nbVar,dimTabVarCond
  integer            :: indice,i,j,k,l 
  integer            :: indiceMin,indiceMax
  integer            :: spec_rank,ligne,var
  integer            :: out_of_bounds
  integer            :: out_of_boundsGlob
  logical            :: verbose,stopLoop

  success =.false.
  out_of_bounds = 0
  nbVar = size(inVarCond)
  dimTabVarCond = size(tabEO,2)
  spec_rank=getnbmycpu()
  verbose = .false.
  
  if (verbose) then
    if (spec_rank .eq. 0 ) write (6,'(a,i0,a,i0)') '[INFO AVG EXPECTATION] nbVar= ',nbVar,' dimTabVarCond= ',dimTabVarCond
    if ( nbVar .ne. size(tabEO,1)-1 ) then
      write (6,'(a)') '[ERROR] on the number of conditional variable'  
      return
    endif
  endif
  if ( .not. samelayout(EOFieldOut,inVarCond(1))) then
    write (6,'(a)') '[ERROR] real_data_layout do not match in function buildEOFieldFromEO'
    return
  endif
!Verification de l'ordre du fichier a explorer
  ligne = 2
  stopLoop = .false.
  do while ( ligne .eq. dimTabVarCond .and. .not. stopLoop )
    if ( tabEO( nbVar , ligne-1) .gt. tabEO( nbVar , ligne )) then
      stopLoop = .true.
      if (spec_rank .eq. 0) write (6,'(a,1x,i0,1x,a,1x,i0)') '[WARNING] optimal estimator is not sorted at var',nbVar,'at line',ligne
    endif
    if ( nbVar .ge. 2 ) then
      do var = nbVar-1,-1,1
        if ( (tabEO(var,ligne-1) .gt. tabEO(var,ligne)) .and. &
           & equalToVal(tabEO(var+1,ligne-1),tabEO(var+1,ligne)) ) then
           stopLoop = .true.
           if (spec_rank .eq. 0) write (6,'(a,1x,i0,1x,a,1x,i0)') '[WARNING] optimal estimator is not sorted at var',var,'at line',ligne
        endif
      enddo
    endif
    ligne = ligne + 1
  enddo
  if (stopLoop) then
    write (6,'(a)') '[ERROR] in buildEOFieldFromEO: file is not correctly ordered!'
    return
!  else
!    if (spec_rank .eq. 0) write (6,'(a)') '[INFO buildEOFieldFromEO] File is correct!' 
  endif
  !Recherche de points
  if ( nbVar .eq. 1 ) then
    if (spec_rank .eq. 0 ) write (6,'(a)') '[INFO buildEOFieldFromEO] Look for values in file with 1 variable' 
    do k=EOFieldOut%zmin,EOFieldOut%zmax
      do j=EOFieldOut%ymin,EOFieldOut%ymax
        do i=EOFieldOut%xmin,EOFieldOut%xmax
          indice = seekDichoto(inVarCond(1)%values(i,j,k),tabEO(1,:),1,dimTabVarCond )
          !!!!!!!!!!!!!!!!!!!!!
          !Verif sortie indice
          if ( indice .eq. -1 ) then
            write (6,'(a)') '[ERROR] seekDichoto 1 returns a bad index in buildEOFieldFromEO: abort !'
            return
          endif
          if ( indice .eq. 1 .or. indice .eq. dimTabVarCond ) then
!            if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
!                            &first conditional var=',inVarCond(1)%values(i,j,k),' min found = ',tabEO(1,1),' max found = ',tabEO(1,dimTabVarCond) 
             out_of_bounds = out_of_bounds + 1
          endif
          !!!!!!!!!!!!!!!!!!!!!
          EOFieldOut%values(i,j,k) = tabEO( nbVar+1 , indice)
          if (verbose) then
            if ( i.eq.EOFieldOut%xmax/2 .and. j.eq.EOFieldOut%ymax/2 .and. k.eq.EOFieldOut%zmax/2 .and. &
               & spec_rank .eq. 0 ) then
              write (6,'(a,i0,a,i0)') '[INFO AVG EXPECTATION] nbVar= ',nbVar,' dimTabVarCond= ',dimTabVarCond
              write (6,'(a,1x,f15.7)') '[INFO AVG EXPECTATION] LOOK FOR :inVarCond(1)=',inVarCond(1)%values(i,j,k)
              write (6,'(a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] BEFORE : tabEO(1,',indice-1,')=',tabEO(1,indice-1)
              write (6,'(a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] MATCH : tabEO(1,',indice,')=',tabEO(1,indice)
              write (6,'(a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] AFTER : tabEO(1,',indice+1,')=',tabEO(1,indice+1)
              write (6,'(a,1x,f15.7)') '[INFO AVG EXPECTATION] dTidxi FOUND=',tabEO(2,indice)
            endif
          endif
        enddo
      enddo
    enddo
  endif
  if (nbVar .eq. 2 ) then
    if (spec_rank .eq. 0 ) write (6,'(a)') '[INFO buildEOFieldFromEO] Look for values in file with 2 variables' 
    do k=EOFieldOut%zmin,EOFieldOut%zmax
      do j=EOFieldOut%ymin,EOFieldOut%ymax
        do i=EOFieldOut%xmin,EOFieldOut%xmax
          indice = seekDichoto(inVarCond(2)%values(i,j,k),tabEO(2,:),1,dimTabVarCond,.true.)
          !!!!!!!!!!!!!!!!!!!!!
          !Verif sortie indice
          if ( indice .eq. -1 ) then
            write (6,'(a)') '[ERROR] seekDichoto 2 returns a bad index in buildEOFieldFromEO: abort !'
            return
          endif
          if ( indice .eq. 1 .or. indice .eq. dimTabVarCond ) then
!            if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
!                            &first conditional var=',inVarCond(2)%values(i,j,k),' min found = ',tabEO(2,1),' max found = ',tabEO(2,dimTabVarCond) 
            out_of_bounds = out_of_bounds + 1
          endif
          !!!!!!!!!!!!!!!!!!!!!
          indiceMin = indice
          if ( indiceMin .gt. 1 ) then
            do while ( equalToVal(tabEO(2,indiceMin),tabEO(2,indiceMin-1) ) .and. &
                   & indiceMin .gt. 2 )
              indiceMin = indiceMin -1 
            enddo
          endif
          indiceMax = indice
          if (indiceMax .lt. dimTabVarCond ) then
            do while ( equalToVal(tabEO(2,indiceMax),tabEO(2,indiceMax+1) ) .and. &
                   & indiceMax .lt. dimTabVarCond-1 )
              indiceMax = indiceMax +1 
            enddo
          endif
          indice = seekDichoto(inVarCond(1)%values(i,j,k),tabEO(1,:),indiceMin,indiceMax)
          !!!!!!!!!!!!!!!!!!!!!
          !Verif sortie indice
          if ( indice .eq. -1 ) then
            write (6,'(a)') '[ERROR] seekDichoto 3 returns a bad index in buildEOFieldFromEO: abort !'
            return
          endif
          if ( indice .eq. indiceMin .or. indice .eq. indiceMax) then
!            if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
!                            &first conditional var=',inVarCond(1)%values(i,j,k),' min found = ',tabEO(1,indiceMin),' max found = ',tabEO(1,indiceMax) 
            out_of_bounds = out_of_bounds + 1
          endif
          !!!!!!!!!!!!!!!!!!!!!
          EOFieldOut%values(i,j,k) = tabEO( 3 , indice)
          if (verbose) then
             if ( i.eq.EOFieldOut%xmax/2 .and. j.eq.EOFieldOut%ymax/2 .and. k.eq.EOFieldOut%zmax/2 .and. &
                & spec_rank .eq. 0 ) then
               do l = 1,nbVar
                 write (6,'(a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] LOOK FOR :inVarCond(',l,')=',inVarCond(l)%values(i,j,k)
               enddo
               do l = 1,nbVar
                if (indice .gt. 0) then
                 write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] BEFORE : tabEO(',l,',',indice-1,')=',tabEO(l,indice-1)
                endif
               enddo
               do l = 1,nbVar
                 write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] MATCH : tabEO(',l,',',indice,')=',tabEO(l,indice)
               enddo
               do l = 1,nbVar
                if (indice .lt. size(tabEO,2)) then
                  write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] AFTER : tabEO(',l,',',indice+1,')=',tabEO(l,indice+1)
                endif
               enddo
               write (6,'(a,1x,f15.7)') '[INFO AVG EXPECTATION] dTidxi FOUND=',tabEO(3,indice)
               write (6,'(a)') '-'
             endif
          endif
        enddo
      enddo
    enddo
  endif
  if (nbVar .gt. 2) then
    if (spec_rank .eq. 0 ) write (6,'(a,1x,i0,1x,a)') '[INFO buildEOFieldFromEO] Look for values in file with',nbVar,'variables' 
    do k=EOFieldOut%zmin,EOFieldOut%zmax
      do j=EOFieldOut%ymin,EOFieldOut%ymax
        do i=EOFieldOut%xmin,EOFieldOut%xmax
          indice = seekDichoto(inVarCond(nbVar)%values(i,j,k),tabEO(nbVar,:),1,dimTabVarCond)
          !!!!!!!!!!!!!!!!!!!!!
          !Verif sortie indice
          if ( indice .eq. -1 ) then
            write (6,'(a)') '[ERROR] seekDichoto 2 returns a bad index in buildEOFieldFromEO: abort !'
            return
          endif
          if ( indice .eq. 1 .or. indice .eq. dimTabVarCond ) then
            if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
                            &first conditional var=',inVarCond(nbVar)%values(i,j,k),' min found = ',tabEO(nbVar,1),&
                            &' max found = ',tabEO(nbVar,dimTabVarCond) 
            out_of_bounds = out_of_bounds + 1
          endif
          !!!!!!!!!!!!!!!!!!!!!
          indiceMin = indice
          if ( indiceMin .gt. 1 ) then
            do while ( equalToVal(tabEO(nbVar,indiceMin),tabEO(nbVar,indiceMin-1) ) .and. &
                   & indiceMin .gt. 2 )
              indiceMin = indiceMin -1 
            enddo
          endif
          indiceMax = indice
          if (indiceMax .lt. dimTabVarCond ) then
            do while ( equalToVal(tabEO(nbVar,indiceMax),tabEO(nbVar,indiceMax+1) ) .and. &
                   & indiceMax .lt. dimTabVarCond-1 )
              indiceMax = indiceMax +1 
            enddo
          endif
          do l=nbVar-1,2,-1
            indice = seekDichoto(inVarCond(l)%values(i,j,k),tabEO(l,:),indiceMin,indiceMax)
            !!!!!!!!!!!!!!!!!!!!!
            !Verif sortie indice
            if ( indice .eq. -1 ) then
              write (6,'(a)') '[ERROR] seekDichoto 3 returns a bad index in buildEOFieldFromEO: abort !'
              return
            endif
            if ( indice .eq. indiceMin .or. indice .eq. indiceMax) then
              if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
                                     &first conditional var=',inVarCond(l)%values(i,j,k),' min found = ',tabEO(l,indiceMin),&
                                     &' max found = ',tabEO(l,indiceMax) 
              out_of_bounds = out_of_bounds + 1
            endif
            !!!!!!!!!!!!!!!!!!!!!
            indiceMin = indice
            if ( indiceMin .gt. 1 ) then
              do while ( equalToVal(tabEO(nbVar,indiceMin),tabEO(l,indiceMin-1) ) .and. &
                    & indiceMin .gt. 2 )
                indiceMin = indiceMin -1 
              enddo
            endif
            indiceMax = indice
            if (indiceMax .lt. dimTabVarCond ) then
              do while ( equalToVal(tabEO(nbVar,indiceMax),tabEO(l,indiceMax+1) ) .and. &
                     & indiceMax .lt. dimTabVarCond-1 )
                indiceMax = indiceMax +1 
              enddo
            endif
          enddo
          indice = seekDichoto(inVarCond(1)%values(i,j,k),tabEO(1,:),indiceMin,indiceMax)
          !!!!!!!!!!!!!!!!!!!!!
          !Verif sortie indice
          if ( indice .eq. -1 ) then
            write (6,'(a)') '[ERROR] seekDichoto 4 returns a bad index in buildEOFieldFromEO: abort !'
            return
          endif
          if ( indice .eq. indiceMin .or. indice .eq. indiceMax) then
            !if (spec_rank .eq. 0 ) write (6,'(3(a,f12.3,a))') '[WARNING] in buildEOFieldFromEO: one of the bound is reached for &
            !                &first conditional var=',inVarCond(1)%values(i,j,k),' min found = ',tabEO(1,indiceMin),' max found = ',tabEO(1,indiceMax) 
            out_of_bounds = out_of_bounds + 1
          endif
          !!!!!!!!!!!!!!!!!!!!!
          EOFieldOut%values(i,j,k) = tabEO( nbVar+1 , indice)
          if (verbose) then
             if ( i.eq.EOFieldOut%xmax/2 .and. j.eq.EOFieldOut%ymax/2 .and. k.eq.EOFieldOut%zmax/2 .and. &
                & spec_rank .eq. 0 ) then
               do l = 1,nbVar
                 write (6,'(a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] LOOK FOR :inVarCond(',l,')=',inVarCond(l)%values(i,j,k)
               enddo
               do l = 1,nbVar
                if (indice .gt. 0) then
                 write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] BEFORE : tabEO(',l,',',indice-1,')=',tabEO(l,indice-1)
                endif
               enddo
               do l = 1,nbVar
                 write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] MATCH : tabEO(',l,',',indice,')=',tabEO(l,indice)
               enddo
               do l = 1,nbVar
                if (indice .lt. size(tabEO,2)) then
                  write (6,'(a,i0,a,i0,a,1x,f15.7)') '[INFO AVG EXPECTATION] AFTER : tabEO(',l,',',indice+1,')=',tabEO(l,indice+1)
                endif
               enddo
               write (6,'(a,1x,f15.7)') '[INFO AVG EXPECTATION] dTidxi FOUND=',tabEO(nbVar+1,indice)
               write (6,'(a)') '-'
             endif
          endif
        enddo
      enddo
    enddo
  endif
  out_of_boundsGlob = DoGlobalSum(spec_rank,out_of_bounds)
  if (spec_rank .eq. 0 ) write (6,'(a,1x,f7.3,1x,a)') '[INFO buildEOFieldFromEO] there are',real(out_of_boundsGlob,WP)&
                                                       &/real(EOFieldOut%nx*EOFieldOut%ny*EOFieldOut%nz*nbVar,WP)*100.0_WP,&
                                '% of points out of bounds' 

  success =.true.

 end function buildEOFieldFromEO

!-------------------------------------------------------------------------------------------
!> Compute the SGS term and the viscous term of scalar equation and it dumps results in ASCII
!! files 
!! @author Antoine Volllant 
!!    @param[in]    dTidxiDNS = divergence of SGS flux from DNS filtered     
!!    @param[in]    phi  = array of datalayout which contains the mean conditional variable 
!!    @param[in]    tabEO = optimal estimator
!!    @param[in]    spec_rank = number of process in cpu pool
!!    @param[in]    adimVar = logical value to ask mse adimensioned by variance
!!    @param[in]   Scal            = Scalar to compute 
!!    @param[in]    msedTidxi = mean square error between dTidxi from DNS and dTidxi from optimal estimator 
!!    @param[in]    mseZdTidxi = mean square error between ZdTidxi from DNS and ZdTidxi from optimal estimator 
!!    @param[in]    moyZdTidxiDNS = average of ZdTidxi from DNS 
!!    @param[in]    moyZdTidxiEO = average of ZdTidxi from optimal estimator 
!!    @return       res             = logical value 
!-------------------------------------------------------------------------------------------
 function computeEoObjectives( dTidxiDNS,phi,Scal,tabEO,adimVar,spec_rank,&
                             &msedTidxi,moyZdTidxiDNS,moyZdTidxiEo,&
                             &mseZdTidxi ) result(success)

  use toolbox

  !I/O data
  type(real_data_layout),intent(in)              :: dTidxiDNS,Scal
  type(real_data_layout),dimension(:),intent(in) :: phi 
  integer,intent(in)                             :: spec_rank
  real(WP),dimension(:,:),intent(in)             :: tabEO
  logical,intent(in)                             :: adimVar 
  real(WP),intent(inout)                           :: msedTidxi
  real(WP),intent(inout)                           :: mseZdTidxi
  real(WP),intent(inout)                           :: moyZdTidxiDNS
  real(WP),intent(inout)                           :: moyZdTidxiEo 
  logical                                        :: success

  !Local data
  type(real_data_layout),dimension(:),allocatable :: temp

   success = .false.

   if (.not.  initWorkArray(dTidxiDNS,3,temp)) then
     write (6,'(a)') '[ERROR] allocation of memory failed in computeEoObjectives : abort !'
     return
   endif
   if (.not. buildEOFieldFromEO(phi,tabEO,temp(3)) ) then
     write (6,'(a)') '[ERROR] buildEOFieldFromEO failed in computeEoObjectives : abort !'
     return
   endif
   temp(1)%values = temp(3)%values * Scal%values
   temp(2)%values = dTidxiDNS%values * Scal%values
   msedTidxi = computeMSE(dTidxiDNS,temp(3),spec_rank,adimVar)
   mseZdTidxi = computeMSE(temp(2),temp(1),spec_rank,adimVar)
   moyZdTidxiDNS = computeFieldAvg(temp(2),spec_rank)
   moyZdTidxiEo = computeFieldAvg(temp(1),spec_rank)
   if (.not. deleteWorkArray(temp)) then
     write (6,'(a)') '[ERROR] Memory leak in computeEoObjectives'
   endif
   success = .true.

 end function computeEoObjectives

END MODULE conditional_mean_tools
!> @}
