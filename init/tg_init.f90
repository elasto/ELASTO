MODULE tailorgreen

CONTAINS

SUBROUTINE TailorGreen_init(U, V, W)

  USE datalayout
  USE parser_tools


  IMPLICIT NONE

  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: U
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: V
  TYPE(REAL_DATA_LAYOUT), INTENT(INOUT) :: W
  
  INTEGER :: i,j,k,tgdim,tgdir
  
  REAL(WP):: Lx, Ly, Lz
  REAL(WP):: dx, dy, dz
  REAL(WP):: kx, ky, kz 
  REAL(WP):: yj,zk,xi
  REAL(WP):: twopi, pi
  REAL(WP):: Umag
  

  CALL parser_read('TG dimensions',tgdim)
  IF (.NOT. GetLxLyLz("Length of domain",Lx,Ly,Lz)) THEN
    STOP 'FAILED in TailorGreen_init'
  ENDIF
  
  IF (.NOT. samelayout(U,V,W)) THEN
    STOP 'FAILED layout error in TailorGreen_init'
  ENDIF

  twopi=2.0_WP*acos(-1.0)
  pi=acos(-1.0)

  kx = twopi/Lx 
  ky = twopi/Ly 
  kz = twopi/Lz 

  IF (tgdim.eq.1) THEN
     CALL parser_read('TG direction',tgdir)
     IF (tgdir.eq.1) THEN
        DO k = U%zmin,U%zmax
           zk = Lz/real(U%nz,WP) * real(k-1,WP)
           DO j = U%ymin,U%ymax
              yj = Ly/real(U%ny,WP) * real(j-1,WP)
              U%values(:,j,k) = 0.0_WP
              V%values(:,j,k) =  kz/ky*cos(yj*ky)*sin(zk*kz) 
              W%values(:,j,k) = -sin(yj*ky)*cos(zk*kz)
           ENDDO
        ENDDO
     ELSE IF (tgdir.eq.2) THEN
        DO k = U%zmin,U%zmax
           zk = Lz/real(U%nz,WP) * real(k-1,WP)
           DO i = U%xmin,U%xmax
              xi = Lx/real(U%nx,WP) * real(i-1,WP)
              U%values(i,:,k) =  kz/kx*cos(xi*kx)*sin(zk*kz)
              V%values(i,:,k) =  0.0_WP
              W%values(i,:,k) = -sin(xi*kx)*cos(zk*kz)
           ENDDO
        ENDDO
     ELSE IF (tgdir.eq.3) THEN
        DO j = U%ymin,U%ymax
           yj = Ly/real(U%ny,WP) * real(j-1,WP)
           DO i = U%xmin,U%xmax
              xi = Lx/real(U%nx,WP) * real(i-1,WP)
              U%values(i,j,:) =  ky/kx*cos(xi*kx)*sin(yj*ky)
              V%values(i,j,:) = -sin(xi*kx)*cos(yj*ky)  
              W%values(i,j,:) =  0.0_WP 
           ENDDO
        ENDDO
     END IF
  ELSE IF (tgdim.eq.2) THEN
     dx = Lx/real(U%nx,WP)
     dy = Ly/real(U%ny,WP)
     dz = Lz/real(U%nz,WP)
     CALL parser_read('TG direction',tgdir)
     DO k = U%zmin,U%zmax
        zk = real(k-1,WP)*dz
        DO j = U%ymin,U%ymax
           yj = real(j-1,WP)*dy
           DO i = U%xmin,U%xmax
              xi = real(i-1,WP)*dx
              U%values(i,j,k) =  cos((yj-pi)*ky)*cos(zk*kz)
              V%values(i,j,k) =  cos((xi-pi)*kx)*cos(zk*kz)
              IF (tgdir.eq.3) THEN
                 W%values(i,j,k) =  0.0_WP
              ELSE 
                 W%values(i,j,k) =  cos((xi-pi)*kx)*cos((yj-pi)*ky)
              END IF
           ENDDO
        ENDDO
     ENDDO
  ELSE IF (tgdim.eq.3) THEN
     dx = Lx/real(U%nx,WP)
     dy = Ly/real(U%ny,WP)
     dz = Lz/real(U%nz,WP)
     DO k = U%zmin,U%zmax
        zk = real(k-1,WP)*dz
        DO j = U%ymin,U%ymax
           yj = real(j-1,WP)*dy
           DO i = U%xmin,U%xmax
              xi = real(i-1,WP)*dx
              U%values(i,j,k) =  +   cos(xi)*sin(yj)*cos(zk)
              V%values(i,j,k) =  +   sin(xi)*cos(yj)*cos(zk)
              W%values(i,j,k) =  2.0*sin(xi)*sin(yj)*sin(zk)
           ENDDO
        ENDDO
     ENDDO
  ELSE 
     dx = Lx/real(U%nx,WP)
     dy = Ly/real(U%ny,WP)
     !dz = Lz/real(U%nz,WP)
     DO k = U%zmin,U%zmax
        !zk = real(k-1,WP)*dz
        DO j = U%ymin,U%ymax
           yj = real(j-1,WP)*dy
           DO i = U%xmin,U%xmax
              xi = real(i-1,WP)*dx
              U%values(i,j,:) =  cos(xi)*sin(yj)
              V%values(i,j,:) = -sin(xi)*cos(yj)
              W%values(i,j,:) =  sin(xi)*cos(yj)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  CALL parser_read('TG velocity magnitude',Umag)
  U%values = Umag*U%values
  V%values = Umag*V%values
  W%values = Umag*W%values
  
END SUBROUTINE TailorGreen_init

END MODULE tailorgreen
