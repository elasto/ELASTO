!USEFORTEST toolbox
!USEFORTEST avgcond 
!USEFORTEST postprocess
!> @addtogroup toolbox 
!! @{
!------------------------------------------------------------------------------
!
! MODULE: sort_find_tools 
!
!> @author
!> Antoine Vollant
!
! DESCRIPTION: 
!> 
!------------------------------------------------------------------------------


module sortfind_tools 

  use precision_tools

  implicit none

  public

 interface parallelSort
   module procedure parallelSort_field
   module procedure parallelSort_tab
 end interface parallelSort 

contains


function parallelSort_field(fieldIn,tabOut,me) result(success)

  use datalayout

  type(real_data_layout),intent(in)         :: fieldIn
  real(WP),dimension(:),intent(inout)       :: tabout
  integer,intent(in)                        :: me
  logical                                   :: success 
  real(WP),dimension(:),allocatable         :: tabin
  integer                                   :: i,j,k,l,ires,sizeLocField

  success = .false.
  sizeLocField=(fieldIn%zmax-fieldIn%zmin+1)*(fieldIn%ymax-fieldIn%ymin+1)*(fieldIn%xmax-fieldIn%xmin+1)
  allocate(tabIn(sizeLocField),stat=ires)
  if (ires .ne. 0) return
  l=1
  do k=fieldIn%zmin,fieldIn%zmax
    do j=fieldIn%ymin,fieldIn%ymax
      do i=fieldIn%xmin,fieldIn%xmax
        tabIn(l) = fieldIn%values(i,j,k)
        l = l + 1
      enddo
    enddo
  enddo
  if (.not. parallelSort(tabIn,tabOut,me)) return
  deallocate(tabIn)
  success = .true.

end function parallelSort_field

!---------------------------------------------------------------------------
!> @details 
!> parallel sorting 
!! @author  Antoine Vollant, LEGI
!! @param[in]       fieldIn contains datas to sort over all processors 
!! @param[inout]    fieldOutSorted contains values in fieldIn sorted in i,j,k increasing 
!! @param[in]       sortedArray contains sorted values of fielIn in increasing index
!! @param[in]       my_rank is the number of processor in cpu pool 
!! @return          success equal to .true. if everything is OK and .false. otherwise
!---------------------------------------------------------------------------

function parallelSort_tab(tabIn,tabOut,me) result(success)

  use stat_tools , only : equalToVal 
  use parallel_tools
  use datalayout
  use mpilayout_tools , only : getnbcpus
  use stat_tools , only : equalToVal 

  !I/O data
  real(WP),dimension(:),intent(in)    :: tabIn    !tableau à trier
  real(WP),dimension(:),intent(inout) :: tabOut   !tableau trié
  integer,intent(in)     :: me       !rang processus
  logical                :: success

  !Local data
  integer                             :: firstMe,lastMe
  real(WP),dimension(2)               :: extremBef
  integer,dimension(:),allocatable    :: tabWorkGroup,tabCom
  real(WP),dimension(:),allocatable   :: tabSend,tabRecv
  real(WP) :: valMinDown,valMaxDown,valMinUp,valMaxUp
  integer  :: nbValUp=0,nbValDown=0,iMinDown=0,iMaxDown=0,iMinUp=0,iMaxUp=0
  integer  :: hierarch,nbValEchange,cas,comLoc,ires
  integer  :: nbProc,nbParGroupe,etape,groupWork,sortAll
  logical  :: fusion

  success = .false.

  nbProc=getnbcpus()
  nbParGroupe=size(tabIn,1)
  firstMe=1
  lastMe=nbParGroupe

  fusion=.false.
  if ( equalToVal( real(nint(log(real(nbParGroupe))/log(2.0)))*log(2.0) , log(real(nbParGroupe)) ) ) then 
    if (me .eq. 0) write (6,'(a)') '[INFO // SORT] fusion sort...'
    fusion=.true.
  else
    if (me .eq. 0) write (6,'(a)') '[INFO // SORT] quick sort...'
  endif
  if ( mod(nbParGroupe,2) .ne. 0 ) then
    write (6,'(a)') '[ERROR // SORT] number of values per processus are not pair'
    return 
  endif
  !tableaux
  if (nbParGroupe .lt. 2 ) return
  allocate(tabWorkGroup(nbProc),stat=ires)
  if ( ires .ne. 0 ) return
  allocate(tabCom(nbProc),stat=ires)
  if ( ires .ne. 0 ) return
  tabOut  = tabIn
  !initiallisation
  etape=0
  Cas = 0
  sortAll = 0
  do while ( sortAll .ne. nbProc )
      sortAll = 0
      !creation des groupes
      call formationGroup(etape,me,nbProc,groupWork,hierarch) 
      if (.not. DoAllGatherScal(groupWork,tabWorkGroup,me)) return
      call formationCom(etape,me,nbProc,groupWork,hierarch,comLoc)
      if (.not. DoAllGatherScal(comLoc,tabCom,me)) return
      !tri sequentiel
      if ( fusion ) then 
        if (.not. seqFusionSort(tabOut, tabOut , nbParGroupe )) return 
      else
        if (.not. seqQuickSort(tabOut ,firstMe ,lastMe)) return 
      endif
      !selection du cas à envisager entre deux processus
      if (groupWork .ne. 0) then
          if ( hierarch .eq. 0 ) then
              extremBef(1)=tabOut(firstMe)
              extremBef(2)=tabOut(lastMe)
              if (.not. DoSendVector( extremBef,2,tabCom(me+1) ,groupWork,me)) return 
              if (.not. DoRecvScal( Cas,tabCom(me+1),groupWork,me)) return 
          else if ( hierarch .eq. 1 ) then 
              if (.not. DoRecvVector( extremBef,2,tabCom(me+1) ,groupWork,me)) return 
              if ( (extremBef(2) .lt. tabOut(firstMe)) .or. equalToVal(extremBef(2),tabOut(firstMe)) ) then 
                Cas=1
              else if ( extremBef(1) .gt. tabOut(lastMe) ) then 
                Cas=5
              else if ( extremBef(1) .lt. tabOut(firstMe) .and. &
                      & extremBef(2) .gt. tabOut(firstMe) .and. &
                      & extremBef(2) .lt. tabOut(lastMe) ) then 
                Cas=2
              else if ( extremBef(1) .lt. tabOut(firstMe) .and. &
                      & extremBef(2) .gt. tabOut(lastMe) ) then 
                Cas=3
              else if ( extremBef(1) .gt. tabOut(firstMe).and. &
                      & extremBef(2) .gt. tabOut(lastMe) .and. &
                      & extremBef(1) .lt. tabOut(lastMe) ) then 
                Cas=5
              else if ( extremBef(1) .gt. tabOut(firstMe) .and. & 
                      & extremBef(2) .lt. tabOut(lastMe) ) then 
                Cas=3
              else if ( equalToVal(extremBef(1),tabOut(firstMe)) .and. &
                      & extremBef(2) .lt. tabOut(lastMe) ) then 
                Cas=3
              else if ( (extremBef(1) .gt. tabOut(firstMe)) .and. &
                      & equalToVal(extremBef(2),tabOut(lastMe)) ) then 
                Cas=3
              else if ( equalToVal(extremBef(1),tabOut(firstMe)) .and. &
                      & equalToVal(extremBef(2),tabOut(lastMe)) ) then 
                Cas=3
              endif
              if (.not. DoSendScal( Cas,tabCom(me+1),groupWork,me)) return 
          endif
      !echange des données selon le cas
          select case (Cas)
          case(1)
              sortAll=1
          case(2)
               allocate(tabSend(3),tabRecv(3),stat=ires)
               if (ires .ne. 0) return
               if ( hierarch .eq. 0 ) then
                  iMaxDown=lastMe
                  valMaxDown=tabOut(iMaxDown)
                  if (.not. DoRecvScal( valMinUp,tabCom(me+1),groupWork,me)) return 
                  iMinDown = seekDichoto( valMinUp , tabOut , firstMe, lastMe )
                  nbValDown = iMaxDown - iMinDown + 1
                  tabSend(1)=nbValDown
                  tabSend(2)=iMaxDown
                  tabSend(3)=iMinDown
                  if (.not. DoSendVector( tabSend,3,tabCom(me+1),groupWork,me)) return 
                  if (.not. DoRecvVector( tabRecv,3,tabCom(me+1),groupWork,me)) return 
                  nbValUp=tabRecv(1)
                  iMaxUp= tabRecv(2)
                  iMinUp= tabRecv(3)
              else if ( hierarch .eq. 1 ) then
                  iMinUp=firstMe
                  valMinUp=tabOut(iMinUp)
                  valMaxDown=extremBef(2)
                  if (.not. DoSendScal( valMinUp,tabCom(me+1),groupWork,me)) return 
                  iMaxUp = seekDichoto( valMaxDown , tabOut, firstMe , lastMe)
                  nbValUp = iMaxUp - iMinUp + 1
                  tabSend(1)=nbValUp
                  tabSend(2)=iMaxUp
                  tabSend(3)=iMinUp
                  if (.not. DoRecvVector( tabRecv,3,tabCom(me+1),groupWork,me)) return
                  if (.not. DoSendVector( tabSend,3,tabCom(me+1),groupWork,me)) return
                  nbValDown=tabRecv(1)
                  iMaxDown= tabRecv(2)
                  iMinDown= tabRecv(3)
              endif
              deallocate(tabRecv,tabSend)
              if ( nbValUp .gt. nbValDown ) then
                  nbValUp = nbValDown
                  iMaxUp = iMinUp + nbValDown - 1
              elseif ( nbValUp .lt. nbValDown ) then
                  nbValDown=nbValUp
                  iMinDown = iMaxDown - nbValUp + 1
              endif
              if ( hierarch .eq. 0 ) then
                  valMinDown=tabOut(iMinDown)
                  if (.not. DoSendScal( valMinDown,tabCom(me+1),groupWork,me)) return
                  if (.not. DoRecvScal( valMaxUp,tabCom(me+1),groupWork,me  )) return
              else if ( hierarch .eq. 1 ) then
                  valMaxUp=tabOut(iMaxUp)
                  if (.not. DoRecvScal( valMinDown,tabCom(me+1),groupWork,me)) return
                  if (.not. DoSendScal( valMaxUp,tabCom(me+1),groupWork,me)  ) return
              endif
              do while ( valMaxUp .gt. valMinDown  )
                  iMinDown = iMinDown + 1
                  iMaxUp = iMaxUp - 1
                  if ( hierarch .eq. 0 ) then
                      valMinDown=tabOut(iMinDown)
                      if (.not. DoSendScal( valMinDown,tabCom(me+1),groupWork,me)) return
                      if (.not. DoRecvScal( valMaxUp,tabCom(me+1),groupWork,me)  ) return
                  else if ( hierarch .eq. 1 ) then
                      valMaxUp=tabOut(iMaxUp)
                      if (.not. DoRecvScal( valMinDown,tabCom(me+1),groupWork,me)) return 
                      if (.not. DoSendScal( valMaxUp,tabCom(me+1),groupWork,me)  ) return
                  endif
              end do
              nbValEchange = iMaxDown - iMinDown + 1
              allocate(tabSend(nbValEchange),tabRecv(nbValEchange),stat=ires)
              if (ires .ne. 0 ) return
              if (hierarch .eq. 0) tabSend = tabOut(iMinDown:iMaxDown)
              if (hierarch .eq. 1) tabSend = tabOut(iMinUp:iMaxUp)
              if (.not. DoSendRecvVector(tabSend,nbValEchange,tabCom(me+1),groupWork,&
                             &tabRecv, nbValEchange,tabCom(me+1),groupWork,me)) return 
              if (hierarch .eq. 0) then
                  tabOut(iMinDown:nbParGroupe)=tabRecv(1:nbValEchange)
              else if (hierarch .eq. 1) then
                  tabOut(iMinUp:iMaxUp)=tabRecv(1:nbValEchange)
              endif
              deallocate(tabSend)
              deallocate(tabRecv)
          case(3)
              allocate(tabSend(1:nbParGroupe/2),tabRecv(1:nbParGroupe/2),stat=ires)
              if (ires .ne. 0 ) return
              if ( hierarch .eq. 1 ) then
                  tabSend(1:nbParGroupe/2) = tabOut(1:nbParGroupe/2)
              else if ( hierarch .eq. 0 ) then
                  tabSend(1:nbParGroupe/2) = tabOut(nbParGroupe/2+1:lastMe)
              endif
              if (.not. DoSendRecvVector(tabSend,nbParGroupe/2,tabCom(me+1),groupWork,&
                         &tabRecv,nbParGroupe/2,tabCom(me+1),groupWork,me)) return 
              if (hierarch .eq. 0) then
                  tabOut(nbParGroupe/2+1:nbParGroupe) = tabRecv(1:nbParGroupe/2)
              else if (hierarch .eq. 1) then
                  tabOut(1:nbParGroupe/2)=tabRecv(1:nbParGroupe/2)
              endif
              deallocate(tabSend)
              deallocate(tabRecv)
          case(5)
              if (.not. DoSendRecvVector(tabOut,nbParGroupe,tabCom(me+1),groupWork,&
                      &tabOut,nbParGroupe,tabCom(me+1),groupWork,me)) return 
          end select
      endif
      !somme du sortAll
      sortAll=DoGlobalSum(me,sortAll)
      etape = etape + 1
!      if (me.eq.0) write(6,'(3(a,1x,i0,1x))') '[INFO // SORT] Step',etape,'->',sortAll,'/',nbProc
   enddo
   !tri sequentiel
   if (fusion) then
     if (.not. seqFusionSort(tabOut , tabOut , nbParGroupe )) return 
   else
     if (.not. seqQuickSort(tabOut ,firstMe ,lastMe)) return 
   endif
   deallocate(tabWorkGroup)
   deallocate(tabCom)

   success = .true.

 end function parallelSort_tab


 
!---------------------------------------------------------------------------
!> @details 
!> sequential sorting 
!! @author  Antoine Vollant, LEGI
!! @param[in]       tabIn contains datas to sort over one processor 
!! @param[inout]    tabOut contains values of tabIn sorted in index increasing 
!! @param[in]       sizeTab is the size of the array tabIn 
!! @return          success equal to .true. if everything is OK and .false. otherwise
!---------------------------------------------------------------------------
 function seqFusionSort( tabIn, tabOut , sizeTab) result(success)

 
   !I/O data
   real(WP),intent(in),dimension(:) :: tabIn 
   real(WP),intent(inout),dimension(:) :: tabOut
   integer,intent(in) :: sizeTab 
   logical            :: success
 
   !Local data
   integer :: i,j,k,l,m,n,offSet
   integer :: tailleUplet,nbUplet,reste
   integer :: tailleUpletDown,tailleUpletUp
   integer :: up,down
   integer :: entier
   real(WP),dimension(:),allocatable :: tabTemp
   real(WP) :: reel,diff
   integer,dimension(:,:),allocatable :: tabTailleUplet
   integer,dimension(:,:,:),allocatable :: tabTailleUpletPrec
   integer,dimension(:),allocatable :: tabNbUplet
   integer :: ires
 
   success = .true.

   allocate(tabTemp(sizeTab),stat=ires)
   if (  ires .ne. 0 ) then
     write (6,'(a)') '[ERROR] cannot allocate memory'
     return
   endif
   tabTemp=0.0_WP
 ! copie des valeurs
   tabOut=tabIn
   reel = log(real(sizeTab,WP))/log(2.0_WP) 
   entier = floor(reel) 
   diff = reel - real(entier,WP)
  
 ! test si la taille du tableau est une puissance de 2    
   if ( diff .eq. 0. ) then
     i=1
     tailleUplet=2
     do while ( tailleUplet .le. sizeTab )
       nbUplet=sizeTab/tailleUplet
       tailleUpletDown=tailleUplet/2
       tailleUpletUp=tailleUplet/2
       do j=0,nbUplet-1
         m=0
         l=0
         k=0
         do while( k .lt. tailleUplet )
           up=m+j*tailleUplet+tailleUpletDown
           down=l+j*tailleUplet
           if ( (m .lt. tailleUpletUp) .and. (l .lt. tailleUpletDown)  ) then
             if ( tabOut(down+1) .gt. tabOut(up+1) ) then
               tabTemp(k+1) = tabOut(up+1)
               m=m+1
             else
               tabTemp(k+1) = tabOut(down+1)
               l=l+1
             endif
             k=k+1
           endif
           if ( m .ge. tailleUpletUp .and. l .lt. tailleUpletDown  ) then
             tabTemp(k+1) = tabOut(down+1)
             l=l+1
             k=k+1
           endif
           if ( (m .lt. tailleUpletUp) .and. (l .ge. tailleUpletDown)  ) then
             tabTemp(k+1) = tabOut(up+1)
             m=m+1
             k=k+1
           endif
         enddo
         do l=0,tailleUplet-1
           tabOut(j*tailleUplet+l+1)=tabTemp(l+1)
         enddo
       enddo
       i=i+1
       tailleUplet=2**i
     enddo
   else
   ! Si la taille du tableau n'est pas une puissance de 2
     !tabNbUplet = calloc ( entier+1 , sizeof( unsigned int *) )
     allocate(tabNbUplet(entier+1),stat=ires)
     if (  ires .ne. 0 ) then
       write (6,'(a)') '[ERROR] cannot allocate memory'
       return
     endif
     allocate(tabTailleUplet(entier+1,sizeTab/2+1),stat=ires)
     if (  ires .ne. 0 ) then
       write (6,'(a)') '[ERROR] cannot allocate memory'
       return
     endif
     allocate(tabTailleUpletPrec(entier+1,sizeTab/4+1,2),stat=ires)
     if (  ires .ne. 0 ) then
       write (6,'(a)') '[ERROR] cannot allocate memory'
       return
     endif
     tabNbUplet(1)=sizeTab/2
     reste=sizeTab-tabNbUplet(1)*2
     do j=0,tabNbUplet(1)-1
       tabTailleUplet(1,j+1)=2
     enddo
     if (reste > 0) then
       tabTailleUplet(1,tabNbUplet(1)) = 1
       tabNbUplet(1)= tabNbUplet(1) +1
     endif
     i=0
     do while ( tabNbUplet(i+1) > 1 )
       if ( mod(tabNbUplet(i+1),2) .ne. 0 ) then
         do j=-1,tabNbUplet(i+1)/2
           tabTailleUplet(i+1+1,j+1) = tabTailleUplet(i+1,2*j+1) + tabTailleUplet(i+1,2*j+1+1)
           tabTailleUpletPrec(i+1,j+1,1) = tabTailleUplet(i+1,2*j+1) 
           tabTailleUpletPrec(i+1,j+1,2) = tabTailleUplet(i+1,2*j+1+1) 
         enddo
         tabNbUplet(i+1+1) = tabNbUplet(i+1)/2
       else
         if ( tabTailleUplet(i+1,1) .gt. tabTailleUplet(i+1,tabNbUplet(i+1)-1+1 ) ) then
           do j=tabNbUplet(i-1)/2,1,-1
             tabTailleUplet(i+1+1,j+1) = tabTailleUplet(i+1,2*j-1+1) + tabTailleUplet(i+1,2*j+1)
             tabTailleUpletPrec(i+1,j+1,2) = tabTailleUplet(i+1,2*j+1) 
             tabTailleUpletPrec(i+1,j+1,1) = tabTailleUplet(i+1,2*j-1+1) 
           enddo
           tabTailleUplet( i+1+1,1) = tabTailleUplet( i+1,1)
           tabTailleUpletPrec(i+1,1,1) = tabTailleUplet(i+1,1) 
           tabTailleUpletPrec(i+1,1,2) = 0 
           tabNbUplet( i+1+1 ) = tabNbUplet( i+1 )/2+1
         else
           do j=0,tabNbUplet(i+1)/2-1
             tabTailleUplet(i+1+1,j+1) = tabTailleUplet(i+1,2*j+1) + tabTailleUplet(i+1,2*j+1+1)
             tabTailleUpletPrec(i+1,j+1,1) = tabTailleUplet(i+1,2*j+1) 
             tabTailleUpletPrec(i+1,j+1,2) = tabTailleUplet(i+1,2*j+1+1) 
           enddo
           tabTailleUplet( i+1+1,tabNbUplet(i+1)/2+1 ) = tabTailleUplet( i+1,tabNbUplet(i+1)-1+1 )
           tabTailleUpletPrec(i+1,tabNbUplet(i+1)/2+1,1) = tabTailleUplet( i+1,tabNbUplet(i+1)-1+1 )  
           tabTailleUpletPrec(i+1,tabNbUplet(i+1)/2+1,2) = 0 
           tabNbUplet( i+1+1 ) = tabNbUplet( i+1 )/2+1
         endif
       endif
       i=i+1
     enddo
     do i=0,entier+1-1
       do j=0,tabNbUplet(i+1)-1
         m=0
         l=0
         k=0
         tailleUplet = tabTailleUplet(i,j)
         if ( i .eq. 0) then
           tailleUpletDown = tailleUplet/2 
           tailleUpletUp = tailleUplet - tailleUplet/2 
         else
           tailleUpletDown = tabTailleUpletPrec(i-1+1,+1,1) 
           tailleUpletUp = tabTailleUpletPrec(i-1+1,j+1,2) 
           if ( tailleUpletUp .eq. 0 ) then
             k = tailleUplet
           endif
         endif
         if ( k .lt. tailleUplet ) then
           do while ( k .lt. tailleUplet )
             if (k .eq. 0) then
               offSet=0
               do n=0,j-1
                 offSet= offSet + tabTailleUplet(i+1,n+1)
               enddo
             endif
             up=m+offSet+tailleUpletDown
             down=l+offSet
             if ( (m .lt. tailleUpletUp) .and. (l .lt. tailleUpletDown) ) then
               if ( tabOut(down+1) .gt. tabOut(up+1) ) then
                 tabTemp(k+1) = tabOut(up+1)
                 m=m+1
               else
                 tabTemp(k+1) = tabOut(down+1)
                 l=l+1
               endif
               k=k+1
             endif
             if ( (m .ge. tailleUpletUp) .and. (l < tailleUpletDown)  ) then
               tabTemp(k+1) = tabOut(down+1)
               l=l+1
               k=k+1
             endif
             if ( (m .lt. tailleUpletUp) .and. (l .ge. tailleUpletDown)  ) then
               tabTemp(k+1) = tabOut(up+1)
               m=m+1
               k=k+1
             endif
           enddo
           do l=0,tailleUplet-1
             tabOut(offSet+l+1)=tabTemp(l+1)
           enddo
         endif
       enddo
     enddo
     deallocate(tabNbUplet)
     deallocate(tabTailleUplet)
     deallocate(tabTailleUpletPrec)
   endif
   deallocate(tabTemp)

   success = .true.


 end function seqFusionSort


!---------------------------------------------------------------------------
!> @details 
!> sequential sorting 
!! @author  Antoine Vollant, LEGI
!! @param[in]       tabIn contains datas to sort over one processor 
!! @param[inout]    tabOut contains values of tabIn sorted in index increasing 
!! @param[in]       sizeTab is the size of the array tabIn 
!! @return          success equal to .true. if everything is OK and .false. otherwise
!---------------------------------------------------------------------------
 function seqInsertionSort(tabIn,sizeTab,inf,sup) result(success)

  use mpilayout_tools

  !I/O data
  real(WP),dimension(:),intent(inout)           :: tabIn
  integer,intent(in),optional                   :: sizeTab
  integer,intent(in),optional                   :: inf,sup 
  logical                                       :: success

  !Local data
  integer :: i,indiceSup,infLoc,isave

  success = .false.

  if (present(inf) .and. present(sup)) then
    i=inf
    indiceSup=sup
  elseif (present(sizeTab))then
    i=1
    indiceSup=sizeTab
  else
    return
  endif
  infLoc=i
  do while ( i .le. indiceSup-1)
    isave = i
    do while ( (i.gt.infLoc) .and. (tabIn(i).gt.tabIn(i+1)) )
      call swapVal(tabIn(i+1),tabIn(i))
      i = i - 1
    enddo
    if ( i.eq.infLoc .and. tabIn(i).gt.tabIn(i+1)) then
      call swapVal(tabIn(i),tabIn(i+1))
    endif
    i = isave
    i=i+1
  enddo

  success = .true.

 end function seqInsertionSort



!---------------------------------------------------------------------------
!> @details 
!> Quick sort
!! @author  Antoine Vollant, LEGI
!! @param[in]       tabIn contains datas to sort over one processor 
!! @param[inout]    tabOut contains values of tabIn sorted in index increasing 
!! @param[in]       sizeTab is the size of the array tabIn 
!! @return          success equal to .true. if everything is OK and .false. otherwise
!---------------------------------------------------------------------------
 recursive function seqQuickSort(tab,firstIndice,lastIndice)  result(success)

   use stat_tools , only : equalToVal

   !I/O data
   REAL(WP), DIMENSION(:), INTENT(INOUT) :: tab
   integer,intent(in)                    :: firstIndice,lastIndice 
   logical                               :: success 
  
   !Local data
   REAL(WP)     :: pivot 
   INTEGER      :: i,j

   success = .false.

   if( firstIndice-lastIndice+1 .gt. size(tab,1)) return

   if (firstIndice-lastIndice .le. 14) then
     if (.not. seqInsertionSort(tab,inf=firstIndice,sup=lastIndice)) return
   else
     pivot=tab( (lastIndice+firstIndice)/2 )
     i = firstIndice
     j = lastIndice
     do while(i .lt. j)
       do while ( tab(i) .lt. pivot .or. equalToVal(pivot,tab(i)) )
         i=i+1
       enddo
       do while ( tab(j) .gt. pivot .or. equalToVal(pivot,tab(j)) )
         j=j-1
       enddo
       call swapVal(tab(i),tab(j),i.lt.j)
     enddo
     if (.not. seqQuickSort(tab,firstIndice,j         )) return
     if (.not. seqQuickSort(tab,j+1        ,lastIndice)) return
   endif

   success = .true.

 end function seqQuickSort

!-------------------------------------------------------------------------------------------
! @details Création des groupes de communication entre les processus
!          aux differentes etapes
! @param[in]  etape     = cycle de triage
! @param[in]  my_rank   = rang du processus dans le groupe de cpu
! @param[in]  nbProcess = nombre total de processus
! @param[out] groupWork = numéro du groupe dont le processus fait parti 
! @param[out] hierarch = dans chaque groupe, le processus de rang n a
!                        la valeur 0 et le processus de rang n+1 a la
!                        valeur 1 
!-------------------------------------------------------------------------------------------
 subroutine formationGroup(etape,my_rank,nbProcess,groupWork,hierarch)

  !I/O Data
  integer,intent(in)    :: etape,my_rank,nbProcess
  integer,intent(inout) :: groupWork,hierarch

  !ETAPE PAIRE
  if ( mod(etape,2) .eq. 0) then
    groupWork = (my_rank/2)+1
    hierarch = mod((my_rank+2),2)
    !NOMBRE de processus IMPAIR
    if ( mod(nbProcess,2) .ne. 0 ) then
      if (my_rank == nbProcess-1) then
        groupWork = 0
      endif
    endif
  endif
  !ETAPE IMPAIRE
  if ( mod(etape,2) .ne. 0) then 
    !Cas particulier de 2 processus
    if ( nbProcess .eq. 2) then
      groupWork= 1
      hierarch = my_rank
    else
      groupWork = (my_rank+1)/2
      hierarch = mod((my_rank+1),2)
      if ( my_rank .eq. 0 ) then
        groupWork = 0
      endif
      !NOMBRE de processus PAIR
      if ( (mod(nbProcess,2) .eq. 0) .and. (my_rank .ne. 0)) then
        if (nbProcess .eq. 2 ) then
          groupWork= 1
          hierarch = my_rank
        endif
        if (my_rank .eq. nbProcess-1) then
          groupWork = 0
        endif
      endif
      !NOMBRE de processus IMPAIR
      if ( (mod(nbProcess,2) .ne. 0) .and. (my_rank .ne. 0) ) then
      endif
    endif
  endif

 end subroutine formationGroup

!-------------------------------------------------------------------------------------------
! @details Donne le numéro du processus avec lequel communiquer 
! @param[in]  etape     = cycle de triage
! @param[in]  my_rank   = rang du processus dans le groupe de cpu
! @param[in]  nbProcess = nombre total de processus
! @param[out] groupWork = numéro du groupe dont le processus fait parti 
! @param[out] hierarch  = dans chaque groupe, le processus de rang n a
!                        la valeur 0 et le processus de rang n+1 a la
!                        valeur 1 
! @param[out] ComLoc    = numero du processus avec lequel communiquer
!-------------------------------------------------------------------------------------------
 subroutine formationCom(etape,my_rank,nbProcess,groupWork,hierarch,ComLoc)

  !I/O Data
  integer,intent(in)    :: etape,my_rank,nbProcess
  integer,intent(inout) :: groupWork,ComLoc,hierarch

  !ETAPE PAIRE
  if ( mod(etape,2) .eq. 0) then
    ComLoc = 2*(groupWork)-(hierarch+1)
    !NOMBRE de processus IMPAIR
    if ( mod(nbProcess,2) .ne. 0 ) then
      if (my_rank .eq. nbProcess-1) ComLoc = -1
    endif
  endif
  !ETAPE IMPAIRE
  if ( mod(etape,2) .ne. 0) then
    ComLoc = 2*groupWork-hierarch
    !Si le rang du processus est 0
    if ( my_rank .eq. 0 ) ComLoc = -1
    !NOMBRE de processus PAIR
    if ( mod(nbProcess,2) .eq. 0 ) then
      if ( (my_rank .eq. nbProcess-1) .and. (nbProcess .gt. 2) ) then
        ComLoc = -1
      endif
      if ( nbProcess .eq. 2) ComLoc = 1-hierarch
    endif
  endif

 end subroutine formationCom


!---------------------------------------------------------------------------
!> @details 
!> find the index of the nearest value of targetPoint in the list Values 
!! @author  Antoine Vollant, LEGI
!! @param[in]    targetPoint is the objective value
!! @param[in]    Values is the list of values where to find
!! @param[in]    minIndice is the smaller index in value's list to take account
!! @param[in]    maxIndice is the last index to take into account
!! @return       index of the nearest value of targetPoint in Value
!---------------------------------------------------------------------------


function seekDichoto(targetPoint,Values,minIndice,maxIndice,view) result(indice)
   
  !I/O data
  real(WP), intent(in)               :: targetPoint
  real(WP),dimension(:), intent(in)  :: Values 
  integer,intent(in)                 :: minIndice,maxIndice
  integer                            :: indice 
  logical,optional                   :: view

  !Local data
  integer :: indiceDown,indiceUp,indiceMid

  indice = -1
  if ( maxIndice .gt. size(Values)  .or.  minIndice .lt. 1 ) then
    write(6,'(a)')'[ERROR] borneIndice is not in the range of Values: abort !'
    write(6,'(a,i0,1x,a,i0)')'size(Values)= ',size(Values),'minIndice= ',minIndice
    return
  endif
 
!  if (present(view)) then
!    if (view) then
!      print *,'minIndice=',minIndice
!      print *,'maxIndice=',maxIndice
!      print *,'Values(minIndice)',Values(minIndice) 
!      print *,'targetPoint',targetPoint
!      print *,'Values(maxIndice)',Values(maxIndice)
!    endif
!  endif
 
  if ( targetPoint .le. Values(minIndice) ) then
    indice = minIndice 
!  endif
  elseif ( targetPoint .ge. Values(maxIndice) ) then
    indice = maxIndice
!  endif
  elseif ( targetPoint .gt. Values(minIndice) .and. &
     & targetPoint .lt. Values(maxIndice) ) then
    indiceDown = minIndice
    indiceUp = maxIndice 
    do while ( indiceDown + 1 .ne. indiceUp )
      if ( mod(indiceUp-indiceDown,2) .eq. 0) then
        indiceMid = (indiceUp+indiceDown)/2
      else
        indiceMid = (indiceUp+indiceDown-1)/2
      endif
      if ( Values(indiceMid) .gt. targetPoint ) then
        indiceUp=indiceMid
      else
        indiceDown=indiceMid
      endif
    enddo
    if ( (targetPoint-Values(indiceDown))**2 .gt. &
       & (targetPoint-Values(indiceUp))**2 ) then
      indice = indiceUp
    else
      indice = indiceDown
    endif
  else
   ! if ( isnan(targetPoint) ) then
   !   write(6,'(a)') 'targetPoint is a NaN in seekDichoto!'
   ! else
      write(6,'(a,f25.5,1x,a,f25.5,1x,a,f25.5)')'targetPoint= ',targetPoint,&
                             &'Values(minIndice)= ',Values(minIndice),&
                             &'Values(maxIndice)= ',Values(maxIndice)
   ! endif
  endif

end function seekDichoto


!---------------------------------------------------------------------------
!> @details 
!> find the index of the nearest value of targetPoint in the list Values 
!! @author  Antoine Vollant, LEGI
!! @param[in]    targetPoint is the objective value
!! @param[in]    Values is the list of values where to find
!! @param[in]    firstIndice is the smaller index in value's list to take account
!! @param[in]    borneIndice is the last index to take into account
!! @param[out]   sens equal 1 (resp. -1) if the seek is from low (resp. high)
!!               index to high (resp. low) index 
!! @return       index of the nearest value of targetPoint in Value. If there is
!!               a mistake it return -1
!---------------------------------------------------------------------------
function seekSame(Values,firstIndice,sens) result(indice)

  !I/O data
  real(WP),dimension(:),intent(in) :: Values
  real(WP)                         :: valToTest
  real(WP)                         :: ValTop
  real(WP)                         :: ValDown
  integer,intent(in)               :: firstIndice
  integer,intent(in)               :: sens
  integer                          :: indice

  !Local data
  real(WP)     :: eps
  logical      :: stopLoop

  indice = -1
  if ( sens .ne. -1 .and. sens .ne. 1 ) then
    write(6,'()')'[ERROR] sens do not have the right value: abort !'
    return
  endif
  if ( firstIndice .gt. size(Values)  .or.  firstIndice .lt. 1 ) then
    write(6,'()')'[ERROR] firstIndice is not in the range of Values: abort !'
    return
  endif

  eps = 1.0E-8_WP
  stopLoop = .false.
  indice = firstIndice
  if(sens .eq. -1 .and. firstIndice .eq. 1) then
    indice = 1
  elseif ( sens .eq. 1 .and. firstIndice .eq. size(Values)) then
    indice = size(Values)
  else 
  valToTest =  Values(indice)
  valTop    =  Values(indice+sens)+eps 
  valDown   =  Values(indice+sens)-eps 
  if ( valToTest .lt. valTop .and. &
     & valToTest .gt. valDown ) then
    do while ( valToTest .lt. valTop .and. &
             & valToTest .gt. valDown .and. &
             & .not. stopLoop )
      valToTest =  Values(indice)
      valTop    =  Values(indice+sens)+eps 
      valDown   =  Values(indice+sens)-eps 
      indice = indice + sens
      if ( indice .eq. 1 .or. indice .eq. size(Values) ) then
        stopLoop=.true.
      endif 
    enddo
  endif
  endif

end function seekSame


!---------------------------------------------------------------------------
!> @details 
!> find the index of the nearest value of targetPoint in the list Values 
!! @author  Antoine Vollant, LEGI
!! @param[in]    targetPoint is the objective value
!! @param[in]    Values is the list of values where to find
!! @param[in]    firstIndice is the smaller index in value's list to take account
!! @param[in]    borneIndice is the last index to take into account
!! @param[out]   sens equal 1 (resp. -1) if the seek is from low (resp. high)
!!               index to high (resp. low) index 
!! @return       index of the nearest value of targetPoint in Value. If there is
!!               a mistake it return -1
!---------------------------------------------------------------------------
function seekIter(targetPoint,Values,firstIndice,borneIndice,sens) result(indice)


  !I/O data
  real(WP), intent(in)             :: targetPoint
  real(WP),dimension(:),intent(in) :: Values
  integer,intent(in)               :: firstIndice
  integer,intent(in)               :: borneIndice
  integer,intent(in)               :: sens
  integer                          :: indice
  real(WP)                         :: valToTest
  real(WP)                         :: ValTop
  real(WP)                         :: ValDown

  !Local data
  real(WP)     :: distBef,distAft
  real(WP)     :: eps
  logical      :: stopLoop,gap
  integer      :: indiceSave=0

  indice = -1
  if ( sens .eq. -1 .and. firstIndice .lt. borneIndice ) then
    write(6,'()')'[ERROR] borneIndice > firstIndice : abort !' 
    return
  endif
  if ( sens .eq. 1 .and. firstIndice .gt. borneIndice ) then
    write(6,'()')'[ERROR] borneIndice < firstIndice : abort !' 
    return
  endif
  if ( sens .ne. -1 .and. sens .ne. 1 ) then
    write(6,'()')'[ERROR] sens do not have the right value: abort !'
    return
  endif
  if ( borneIndice .gt. size(Values)  .or.  borneIndice .lt. 1 ) then
    write(6,'()')'[ERROR] borneIndice is not in the range of Values: abort !'
    return
  endif
  if ( firstIndice .gt. size(Values)  .or.  firstIndice .lt. 1 ) then
    write(6,'()')'[ERROR] firstIndice is not in the range of Values: abort !'
    return
  endif

  eps = 1.0E-10_WP
  stopLoop = .false.
  gap = .false.
  indice = firstIndice
  if ( indice .ne. borneIndice ) then
    distBef=( targetPoint-Values(indice) )**2
    distAft=( targetPoint-Values(indice+sens) )**2
    do while ( .not. stopLoop )
      if ( indice .eq. borneIndice-sens ) stopLoop=.true.
      if ( distBef .lt. distAft+eps .and. &
         & distBef .gt. distAft-eps  ) then
        indiceSave = indice
        valToTest = Values(indice)
        valTop    = Values(indice+sens)+eps 
        valDown   = Values(indice+sens)-eps 
        do while ( valToTest .lt. ValTop .and. &
                 & valToTest .gt. ValDown .and. &
                 & .not. stopLoop )
          valToTest = Values(indice)
          valTop    = Values(indice+sens)+eps 
          valDown   = Values(indice+sens)-eps 
          indice = indice + sens
          if ( indice .eq. borneIndice-sens ) then
            stopLoop=.true.
          endif 
        enddo
!        if (.not. stopLoop ) then
          distBef=( targetPoint-Values(indiceSave) )**2
          distAft=( targetPoint-Values(indice+sens) )**2
!        endif
        if (.not. gap) gap =.true.
      endif
      if ( distBef .gt. distAft ) then
        indice = indice + sens
        distBef = distAft
        distAft = (targetPoint-Values(indice))**2
        gap =.false.
      endif
      if ( distBef .lt. distAft ) then
        stopLoop=.true.
        if ( gap ) indice = indiceSave
      endif
    enddo
  else
    indice = borneIndice
  endif

end function seekIter

!---------------------------------------------------------------------------
!> @details 
!> exchange two real values 
!! @author  Antoine Vollant, LEGI
!! @param[inout]    a value to change with b 
!! @param[inout]    b value to change with a 
!! @param[in]       test is a logical. If equal to .false. nothing happen. If equal to .true.
!!                  a and b are exchanged 
!---------------------------------------------------------------------------
subroutine swapVal(a,b,test)

  use precision_tools

  !I/O data
  real(WP),intent(inout)      :: a,b
  logical,intent(in),optional :: test
  !Local data
  real(WP) :: temp

  if ( present(test) ) then
    if (test) then
      temp = a
      a = b
      b = temp 
    endif
  else
    temp = a
    a = b
    b = temp 
  endif

end subroutine swapVal


subroutine testSort(tabStart,tabCurrent,my_rank,tag)

  use parallel_tools
  use datalayout
  use mpilayout_tools , only : getnbcpus
  use stat_tools , only : equalToVal 
  
  !I/O data
  real(WP),dimension(:),intent(in)  :: tabStart,tabCurrent  
  integer,intent(in)                :: my_rank
  character(len=64),intent(in)      :: tag
  !Local data
  integer                           :: i,j,nbproc,nbParGroupe
  integer                           :: erreurConserv,erreurOrdre
  integer                           :: ires,nbTot
  logical                           :: found
  real(WP),dimension(:),allocatable :: tabStartTot
  real(WP),dimension(:),allocatable :: tabCurrentTot

  nbproc=getnbcpus()
  nbParGroupe=size(tabStart,1)
  nbTot=nbproc*nbParGroupe
  allocate(tabStartTot(nbTot),stat=ires)
  if (ires .ne. 0) print *,'[ERROR] in testSort' 
  allocate(tabCurrentTot(nbTot),stat=ires)
  if (ires .ne. 0) print *,'[ERROR] in testSort' 
  if (.not. DoGatherVector(tabStart  ,nbParGroupe,tabStartTot  ,nbParGroupe,0,my_rank)) return
  if (.not. DoGatherVector(tabCurrent,nbParGroupe,tabCurrentTot,nbParGroupe,0,my_rank)) return
  if (my_rank .eq. 0) then
    erreurOrdre   = 0
    erreurConserv = 0
    do i=1,nbTot
      j=1
      found=.false.
      do while ( (j .le. nbTot) .and. (.not. found) )
        if ( equalToval(tabCurrentTot(j),tabStartTot(i)) ) found = .true.
        j=j+1
      enddo
      if ( .not. found ) erreurConserv = erreurConserv + 1 
      if ( i .lt. nbTot) then
        if ( tabCurrentTot(i+1) .lt. tabCurrentTot(i) )then
          erreurOrdre = erreurOrdre + 1 
        endif
      endif
    enddo
    write (6,'(a,a64,a,1x,i0,1x,a,1x,i0)') '[INFO testSort] ',trim(adjustl(tag)),' erreurOrdre',erreurOrdre,'erreurConserv',erreurConserv
  endif
  deallocate(tabStartTot)
  deallocate(tabCurrentTot)

end subroutine testSort

END MODULE sortfind_tools 
!> @}
