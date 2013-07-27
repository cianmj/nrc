SUBROUTINE CALCRATE(n2,dn,epsilon,c,Hl,Ham,dr,dp,da,db,dc)

USE FNVAR

IMPLICIT NONE

INTERFACE
 SUBROUTINE CHECK(n,MAT,char)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  DOUBLE COMPLEX, INTENT(in) :: MAT(0:n,0:n)
  CHARACTER(4),INTENT(in) :: char
 END SUBROUTINE

 SUBROUTINE SVD(loop,dn,dn1,DrLC,DpLC,DaLC,DbLC,HlC,dr,dp,da,db)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: loop, dn, dn1
  DOUBLE COMPLEX, INTENT(in),DIMENSION(0:dn) :: DrLC,DpLC,DaLC,DbLC,HlC
  DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
 END SUBROUTINE

 SUBROUTINE CALCQUAD(dr,dp,da,db)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
 END SUBROUTINE

 SUBROUTINE SINGLE(dr,dp,da,db)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
 END SUBROUTINE

 SUBROUTINE NUMRND(value)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(inout) :: value
 END SUBROUTINE
END INTERFACE

!! Variables passed in/out of subroutine
INTEGER, INTENT(in) :: n2,dn
DOUBLE PRECISION, INTENT(in) :: epsilon
DOUBLE COMPLEX, INTENT(in) :: c(0:n),Hl(0:dn,0:n),Ham(0:n,0:n)
DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
DOUBLE COMPLEX, INTENT(out) :: dc(0:n)
!
INTEGER :: k,l,lim2,dn1
DOUBLE PRECISION :: test,test1
DOUBLE PRECISION, SAVE :: Rsqsum(0:2*dn+1), TsumRsq
DOUBLE COMPLEX :: cmplx1, cmplx2
DOUBLE COMPLEX, DIMENSION(0:n2,0:n) :: Drkl, Dpkl, Dakl, Dbkl
DOUBLE COMPLEX, DIMENSION(0:dn) :: DrLC,DpLC,DaLC,DbLC,HlC
DOUBLE COMPLEX, DIMENSION(0:n) :: expnH,expnR,expnP,expnA,expnB
DOUBLE COMPLEX, DIMENSION(0:n,0:n) :: tpHam, tpR ,tpP, tpA, tpB
CHARACTER(4) :: Drmat,Dpmat,Damat,Dbmat
LOGICAL :: lstop
!
! Intrinsic functions
REAL, INTRINSIC :: REAL
DOUBLE PRECISION, INTRINSIC :: DBLE, IMAG, ABS
DOUBLE COMPLEX, INTRINSIC :: CONJG, MATMUL, SQRT
!
! External function
DOUBLE PRECISION :: factrl

lim2 = 8 !	! Order for the matrix exponential expansion

!!! Character definition
Drmat = 'Drkl' ; Dpmat = 'Dpkl' ; Damat = 'Dakl' ; Dbmat = 'Dbkl'

!!! Default values
dr = 0.0d0; dp = 0.0d0; da = 0.0d0; db = 0.0d0

!!! Matrices needed for propagation of coefficients !!!

!WRITE(*,*) 'Drkl matrix'
DO k=0,n2
 DO l=0,n
  IF (k.eq.l) THEN
   Drkl(k,l) = -(0.0d0,1.0d0)*p/hbar
  ELSE IF (k.eq.(l-1)) THEN
   Drkl(k,l) = SQRT(DBLE(k)+1.0d0)*((-a-(0.0d0,1.0d0)*b)/SQRT(a))
  ELSE IF (k.eq.(l+1)) THEN
   Drkl(k,l) = SQRT(DBLE(k))*( (a-(0.0d0,1.0d0)*b)/SQRT(a))
  ELSE
   Drkl(k,l) = 0.0d0
  ENDIF
 ENDDO
ENDDO

!WRITE(*,*) 'Dpkl matrix'
DO k=0,n2
 DO l=0,n
  IF (k.eq.(l-1)) THEN
   Dpkl(k,l) = (0.0d0,1.0d0)*SQRT((DBLE(k)+1.0d0)/a) / (2.0d0*hbar)
  ELSE IF (k.eq.(l+1)) THEN
   Dpkl(k,l) = (0.0d0,1.0d0)*SQRT(DBLE(k)/a)/(2.0d0*hbar)
  ELSE
   Dpkl(k,l) = 0.0d0
  ENDIF
 ENDDO
ENDDO

!WRITE(*,*) 'Dakl matrix'
DO k=0,n2
 DO l=0,n
  IF (k.eq.(l-2)) THEN
   Dakl(k,l) = SQRT((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))/(4.0d0*a)
  ELSE IF (k.eq.(l+2)) THEN
   Dakl(k,l) = -SQRT(DBLE(k)*(DBLE(k)-1.0d0))/(4.0d0*a)
  ELSE
   Dakl(k,l) = 0.0d0
  ENDIF
 ENDDO
ENDDO

!WRITE(*,*) 'Dbkl matrix'
DO k=0,n2
 DO l=0,n
  IF (k.eq.l) THEN
   Dbkl(k,l) = (0.0d0,1.0d0)*(DBLE(k)+0.5d0)/(2.0d0*a)
  ELSE IF (k.eq.(l-2)) THEN
   Dbkl(k,l) = (0.0d0,1.0d0)*SQRT((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))/(4.0d0*a)
  ELSE IF (k.eq.(l+2)) THEN
   Dbkl(k,l) = (0.0d0,1.0d0)*SQRT(DBLE(k)*(DBLE(k)-1.0d0))/(4.0d0*a)
  ELSE
   Dbkl(k,l) = 0.0d0
  ENDIF
 ENDDO
ENDDO

!! Check that the above matrices are correct
!CALL CHECK(n,Drkl(0:n,0:n),Drmat)
!CALL CHECK(n,Dpkl(0:n,0:n),Dpmat)
!CALL CHECK(n,Dakl(0:n,0:n),Damat)
!CALL CHECK(n,Dbkl(0:n,0:n),Dbmat)
!
!!! For the propagation constraint equations

 DrLC(0:dn) = MATMUL(Drkl(n+1:n2,:),c(:))
 DpLC(0:dn) = MATMUL(Dpkl(n+1:n2,:),c(:))
 DaLC(0:dn) = MATMUL(Dakl(n+1:n2,:),c(:))
 DbLC(0:dn) = MATMUL(Dbkl(n+1:n2,:),c(:))
 HlC(0:dn) = -(0.0d0,1.0d0)*MATMUL(Hl(:,:),c(:))/hbar

!!! Solve system of equations for the unknown time derivatives !!!
 dn1 = dn + 1   ! need to change to correspond with subr.

 CALL SVD(loop,dn,dn1,DrLC,DpLC,DaLC,DbLC,HlC,dr,dp,da,db)
! WRITE(*,*) dr,dp,da,db
! read(*,*)
!!
!! For quadratic potential: --
! CALL CALCQUAD(dr,dp,da,db)
! WRITE(*,*) dr,dp,da,db
! read(*,*)
!!
!! For single wave packet in a harmonic potential: --
! CALL SINGLE(dr,dp,da,db)
! WRITE(*,*) dr,dp,da,db
! read(*,*)

cmplx1 = DCMPLX(dr,dp); cmplx2 = DCMPLX(da,db)
CALL NUMRND(cmplx1); CALL NUMRND(cmplx2);
dr = DBLE(cmplx1); dp = DIMAG(cmplx1)
da = DBLE(cmplx2); db = DIMAG(cmplx2)

!! Verification of the solutions outputed from the subroutine

IF (loop.lt.2) THEN
 Rsqsum(:)=0.0d0; TsumRsq = 0.0d0
ENDIF
lstop = .false.

IF (loop.ne.1) THEN  !!!------ skip the test section
 DO k=0,dn
  test=DBLE(DrLC(k))*dr+DBLE(DpLC(k))*dp+DBLE(DaLC(k))*da+DBLE(DbLC(k))*db
  test1=DIMAG(DrLC(k))*dr+DIMAG(DpLC(k))*dp+DIMAG(DaLC(k))*da+DIMAG(DbLC(k))*db
  IF ((ABS(test-DBLE(HlC(k))).gt.epsilon).OR. &
   & (ABS(test1-DIMAG(HlC(k))).gt.epsilon)) THEN
   WRITE(*,*) 'Program stopped on loop #',loop
   WRITE(*,*) 'Solutions to linslv have lost precision - ',epsilon
   WRITE(*,*) 'Residuals for ',k,' := '
   WRITE(*,*) abs(test-DBLE(HlC(k))),abs(test1-DIMAG(HlC(k)))
   lstop = .true.
  ENDIF

  IF ((ABS(test-DBLE(HlC(k)))**2*dt.gt.1E-50).OR. & 
    & (ABS(test1-DIMAG(HlC(k)))**2*dt.gt.1E-50)) THEN
   Rsqsum(2*k) = Rsqsum(2*k) + (test-DBLE(HlC(k)))**2*dt
   Rsqsum(2*k+1) = Rsqsum(2*k+1) + (test1-DIMAG(HlC(k)))**2*dt
   TsumRsq = TsumRsq + Rsqsum(2*k) + Rsqsum(2*k+1)
  ENDIF
 ENDDO
ENDIF

IF ((TsumRsq/(loop*dt)).gt.epsilon**2) THEN
 WRITE(*,*) 'Residuals have become too large : -', TsumRsq/(loop*dt) 
ENDIF

IF (lstop) STOP !	! Stop program because of residue size


! Matrix exponential expansion for propagation equations
!-!-! (lim2)th order exp. ... for coefficient propagation

!write(*,*) 'ham'
!write(*,*) ham(:,:)

expnH(:) = (0.0d0,0.0d0)
expnR(:) = (0.0d0,0.0d0); expnP(:) = (0.0d0,0.0d0)
expnA(:) = (0.0d0,0.0d0); expnB(:) = (0.0d0,0.0d0)
DO k=1,lim2
 tpHam(:,:)=Ham(:,:)
 tpR(:,:)=Drkl(0:n,0:n); tpP(:,:)=Dpkl(0:n,0:n)
 tpA(:,:)=Dakl(0:n,0:n); tpB(:,:)=Dbkl(0:n,0:n)
 DO l=1,k-1
  tpHam = MATMUL(tpHam,Ham)
  tpR = MATMUL(tpR,Drkl(0:n,:)); tpP = MATMUL(tpP,Dpkl(0:n,:))
  tpA = MATMUL(tpA,Dakl(0:n,:)); tpB = MATMUL(tpB,Dbkl(0:n,:))
 ENDDO
 expnH(:) = expnH(:) + (-(0.0d0,1.0d0)/hbar)**k * dt**(k-1) * &
  & MATMUL(tpHam,c(:)) / factrl(k)
 expnR(:) = expnR(:) + (dr)**k * dt**(k-1) * MATMUL(tpR,c(:))/factrl(k)
 expnP(:) = expnP(:) + (dp)**k * dt**(k-1) * MATMUL(tpP,c(:))/factrl(k)
 expnA(:) = expnA(:) + (da)**k * dt**(k-1) * MATMUL(tpA,c(:))/factrl(k)
 expnB(:) = expnB(:) + (db)**k * dt**(k-1) * MATMUL(tpB,c(:))/factrl(k)
ENDDO

!write(*,*) 'c'
!write(*,*) c(:)
!write(*,*) 'tpham'
!write(*,*) tpham(:,:)
!write(*,*) 'matmul'
!write(*,*) MATMUL(tpHam,c(:))

DO k=0,n
 CALL NUMRND(expnR(k)); CALL NUMRND(expnP(k)); CALL NUMRND(expnA(k))
 CALL NUMRND(expnB(k)); CALL NUMRND(expnH(k))
ENDDO

!write(*,*) 'r'
!write(*,*) expnR(:)
!write(*,*) 'p'
!write(*,*) expnP(:)
!write(*,*) 'a'
!write(*,*) expnA(:)
!write(*,*) 'b'
!write(*,*) expnB(:)
!write(*,*) 'H'
!write(*,*) expnH(:)
!read(*,*)

!!! Coefficient equation !!!

dc(:) = 0.0d0
dc(:) = expnH(:) - expnR(:) - expnP(:) - expnA(:) - expnB(:)

!!
!dc(:) = -((0.0d0,1.0d0)/hbar)*MATMUL(Ham,c) - dr*MATMUL(Drkl(0:n,:),c(:)) - &
!& dp*MATMUL(Dpkl(0:n,:),c(:)) - da*MATMUL(Dakl(0:n,:),c(:)) - &
!& db*MATMUL(Dbkl(0:n,:),c(:))
!!!

RETURN

END SUBROUTINE

!
!!!!!!!!!!!!!!
!
SUBROUTINE CHECK(n,MAT,char)
IMPLICIT NONE

INTEGER :: k,l
INTEGER, INTENT(in) :: n
DOUBLE COMPLEX, INTENT(in) :: MAT(0:n,0:n)
CHARACTER(4),INTENT(in) :: char

DO k=0,n
 DO l=0,n
  IF ( ABS((MAT(k,l)+CONJG(MAT(l,k)))).gt.1E-10) THEN
   WRITE(*,*) k,l
   WRITE(*,*) MAT(k,l),CONJG(MAT(l,k))
   WRITE(*,*) REAL(MAT(k,l)+CONJG(MAT(l,k)))
   WRITE(*,*) char,' matrix is incorrect'
   STOP
  ENDIF
 ENDDO
ENDDO

RETURN
END SUBROUTINE
!
!!!!!!!!!!!!!!
!
SUBROUTINE NUMRND(value)
IMPLICIT NONE

INTEGER :: vk
DOUBLE PRECISION, PARAMETER :: round = 1E-12
DOUBLE PRECISION :: real, imag, nk
DOUBLE COMPLEX, INTENT(inout) :: value
!
INTEGER, INTRINSIC :: NINT

real = DBLE(value)
imag = DIMAG(value)

IF (ABS(real).lt.round) THEN
 vk = NINT(real)
 nk = real - vk
 nk = NINT(nk*round)/round
 real = vk + nk
 value = DCMPLX(real,imag)
ENDIF

IF (ABS(imag).lt.round) THEN
 vk = NINT(imag)
 nk = imag - vk
 nk = NINT(nk*round)/round
 imag = vk + nk
 value = DCMPLX(real,imag)
ENDIF

END SUBROUTINE
