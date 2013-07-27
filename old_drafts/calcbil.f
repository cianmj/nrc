SUBROUTINE CALCBIL(calc,Amtrx,Ham,dr,dp,da,db,c,dc)

USE FNVAR

IMPLICIT NONE

INTERFACE
 SUBROUTINE CHECK(MAT,char)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(in) :: MAT(0:n,0:n)
  CHARACTER(4),INTENT(in) :: char
 END SUBROUTINE

 SUBROUTINE NUMRND(value)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(inout) :: value
 END SUBROUTINE
END INTERFACE

!! Variables passed in/out of subroutine
LOGICAL :: calc
DOUBLE PRECISION, INTENT(in), OPTIONAL :: dr,dp,da,db
DOUBLE PRECISION, INTENT(in) :: Amtrx(0:n,0:n,0:lim)
DOUBLE COMPLEX, INTENT(in), OPTIONAL :: c(0:n)
DOUBLE COMPLEX, INTENT(out) :: Ham(0:n,0:n)
DOUBLE COMPLEX, INTENT(out), OPTIONAL :: dc(0:n)
!
INTEGER :: i,k,l,lim2
DOUBLE PRECISION :: Potvar(0:lim)
DOUBLE COMPLEX, DIMENSION(0:n,0:n) :: KE, PE, Drkl, Dpkl, Dakl, Dbkl
DOUBLE COMPLEX, DIMENSION(0:n,0:n) :: Qmtrx, Qtp, idty, expnQ
CHARACTER(4) :: Drmat,Dpmat,Damat,Dbmat
!
! Intrinsic functions
REAL, INTRINSIC :: REAL
DOUBLE PRECISION, INTRINSIC :: DBLE, IMAG, ABS
DOUBLE COMPLEX, INTRINSIC :: CONJG, MATMUL, SQRT
!
! External function
DOUBLE PRECISION, EXTERNAL :: factrl, DPoten


lim2 = 6 !	! Order for the matrix exponential expansion

!!! Character definition
Drmat = 'Drkl' ; Dpmat = 'Dpkl' ; Damat = 'Dakl' ; Dbmat = 'Dbkl'

DO k=0,n
 DO l=0,n
  IF (k.eq.l) THEN
   KE(k,l) = -(p**2)/(hbar**2)-((b**2+a**2)/a)*(2.0d0*DBLE(k)+1.0d0)

  ELSE IF (k.eq.(l+1)) THEN
   KE(k,l) = (2.0d0*p/hbar)*SQRT(DBLE(k))*(-b-(0.0d0,1.0d0)*a)/SQRT(a)

  ELSE IF (k.eq.(l-1)) THEN
   KE(k,l) = (2.0d0*p/hbar)*SQRT(DBLE(k)+1.0d0)*(-b+(0.0d0,1.0d0)*a)/SQRT(a)

  ELSE IF (k.eq.(l+2)) THEN
   KE(k,l) = SQRT(DBLE(k)*(DBLE(k)-1.0d0))*(((a**2-b**2)/a) - &
 &  2.0d0*(0.0d0,1.0d0)*b)

  ELSE IF (k.eq.(l-2)) THEN
   KE(k,l) = SQRT((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))*(((a**2-b**2)/a) + &
 &  2.0d0*(0.0d0,1.0d0)*b)

  ELSE
   KE(k,l) = 0.0d0

  ENDIF
   KE(k,l) = -(hbar**2/(2.0d0*mass))*KE(k,l)

 ENDDO
ENDDO

!!! Potential Energy matrix

DO i=0,lim
 Potvar(i) = DPoten(i)
ENDDO

! WRITE(*,*) 'Potential energy matrix'
DO k=0,n
 DO l=0,n
  PE(k,l) = 0.0d0
  DO i=0,lim
   IF ((MODULO((k+l+i),2).eq.0).AND.(l.ge.(k-i)).AND.(l.le.(k+i))) THEN
     PE(k,l) = PE(k,l)+Amtrx(k,l,i)*Potvar(i)/((2.0d0*a)**(DBLE(i)/2.0d0))
   ELSE
    CYCLE
   ENDIF
  ENDDO
 ENDDO
ENDDO

!!! Hamlitonian matrix / check whether it is Hermitian !!!

! WRITE(*,*) 'Hamiltonian matrix'
DO k=0,n
 DO l=0,n
  Ham(k,l) = KE(k,l) + PE(k,l)
 ENDDO
ENDDO

! Check whether matrix is Hermitian
DO k=0,n
 DO l=0,k
  IF ((ABS(REAL(Ham(k,l)-CONJG(Ham(l,k)))).gt.epsilon).OR. &
   & (ABS(AIMAG(Ham(k,l)-CONJG(Ham(l,k)))).gt.epsilon)) THEN
    WRITE(UNIT=31,FMT='(A)') 'Hamiltonian matrix is not Hermitian'
    WRITE(UNIT=31,FMT='(A4,I8)') 'loop',loop
    STOP
  ENDIF
 ENDDO
ENDDO

!!!!!!!!!!!!!!!!Find coefficient propagation!!!!!!!!!!!!!!!!!!!!!

IF (calc) THEN

 !!! Matrices needed for propagation of coefficients !!!

 DO k=0,n
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

 DO k=0,n
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

 DO k=0,n
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

 DO k=0,n
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
 CALL CHECK(Drkl(0:n,0:n),Drmat)
 CALL CHECK(Dpkl(0:n,0:n),Dpmat)
 CALL CHECK(Dakl(0:n,0:n),Damat)
 CALL CHECK(Dbkl(0:n,0:n),Dbmat)


 ! Matrix exponential expansion for propagation equations
 !-!-! (lim2)th order exp. ... for coefficient propagation

 Qmtrx(:,:) = - (0.0d0,1.0d0)*Ham(0:n,0:n) - dr*Drkl(0:n,0:n) - &
	& dp*Dpkl(0:n,0:n) - da*Dakl(0:n,0:n) - db*Dbkl(0:n,0:n)

 Qmtrx = Qmtrx*dt/2.0d0

 expnQ(:,:) = (0.0d0,0.0d0); Qtp(:,:) = (0.0d0,0.0d0)
 DO k=0,n
  Qtp(k,k) = (1.0d0,0.0d0)
 ENDDO
 DO k = 0,lim2
  expnQ(:,:) = expnQ(:,:) + Qtp(:,:) / factrl(k)
  Qtp(:,:) = MATMUL(Qtp,Qmtrx)
 ENDDO

 DO k = 0,n
  DO l = 0,n
   CALL NUMRND(expnQ(k,l))
  ENDDO
 ENDDO


 !!! coefficient equation !!!

 idty(:,:) = (0.0d0,0.0d0)
 DO k=0,n
  idty(k,k) = (1.0d0,0.0d0)
 ENDDO

 expnQ(:,:) = expnQ(:,:) - idty(:,:)
 dc(:) = MATMUL(expnQ(:,:),c(:)) / (dt/2.0d0)

 DO k=0,n
  CALL NUMRND(dc(k))
 ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN

END SUBROUTINE

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE CHECK(MAT,char)
USE FNVAR
IMPLICIT NONE

INTEGER :: k,l
DOUBLE COMPLEX, INTENT(in) :: MAT(0:n,0:n)
CHARACTER(4),INTENT(in) :: char

DO k=0,n
 DO l=0,n
  IF ( ABS((MAT(k,l)+CONJG(MAT(l,k)))).gt.epsilon) THEN
   WRITE(*,*) k,l
   WRITE(*,*) MAT(k,l),CONJG(MAT(l,k))
   WRITE(*,*) REAL(MAT(k,l)+CONJG(MAT(l,k)))
   WRITE(*,*) char,' matrix is incorrect'
   WRITE(*,*) 'Program stopped :'
   WRITE(*,*) loop,r,p,a,b
   STOP
  ENDIF

  CALL NUMRND(MAT(k,l))

 ENDDO
ENDDO

RETURN
END SUBROUTINE
!
!!!!!!!!!!!!!!
!
SUBROUTINE NUMRND(value)
USE FNVAR
IMPLICIT NONE

INTEGER :: vk
DOUBLE PRECISION :: real, imag, nk
DOUBLE COMPLEX, INTENT(inout) :: value
LOGICAL :: rounded
LOGICAL, SAVE :: wrt = .true.
!
INTEGER, INTRINSIC :: NINT
DOUBLE PRECISION, INTRINSIC :: ABS

rounded = .false.

real = DBLE(value)
imag = DIMAG(value)

IF (ABS(real).lt.epsilon) THEN
 vk = NINT(real)
 nk = real - vk
 nk = NINT(nk*epsilon)/epsilon
 real = vk + nk
 value = DCMPLX(real,imag)
 rounded = .true.
ENDIF

IF (ABS(imag).lt.epsilon) THEN
 vk = NINT(imag)
 nk = imag - vk
 nk = NINT(nk*epsilon)/epsilon
 imag = vk + nk
 value = DCMPLX(real,imag)
 rounded = .true.
ENDIF

IF ((rounded).AND.(wrt)) THEN
 WRITE(*,*) 'Numbers have been rounded in calcbil.'
 wrt = .false.
ENDIF

END SUBROUTINE
!
!!!!!!!!!!!!!!
