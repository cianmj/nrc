PROGRAM NRGMIN

USE FNVAR
USE VARSAVE

IMPLICIT NONE

INTERFACE
 SUBROUTINE Hamdiag(Ham,opt,nrg,c)
  USE FNVAR
  IMPLICIT NONE
  LOGICAL :: opt
  DOUBLE COMPLEX, intent(in) :: Ham(0:n,0:n)
  DOUBLE COMPLEX, intent(out), OPTIONAL :: c(0:n)
  DOUBLE PRECISION, intent(out) :: nrg(0:n)
 END SUBROUTINE
END INTERFACE

! Input parameters
DOUBLE PRECISION :: v,step,stepdef
!
INTEGER :: k,l,i,j
INTEGER :: s,err,err1,max
DOUBLE PRECISION :: DrE,DpE,DaE,DbE, randnum
DOUBLE PRECISION :: rksum,pksum,pksum2,rxvl,pxvl,Eexp
!
LOGICAL :: opt
DOUBLE PRECISION, ALLOCATABLE :: cH(:,:),Amtrx(:,:,:),nrg(:),potvar(:)
DOUBLE COMPLEX, ALLOCATABLE :: c(:)
DOUBLE COMPLEX, ALLOCATABLE :: KE(:,:),PE(:,:),Ham(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: DpT(:,:),DaT(:,:),DaV(:,:),DbT(:,:)
!
!! Gradient method
DOUBLE PRECISION :: unit
DOUBLE PRECISION, DIMENSION(4) :: grad,ugrad

! External functions
DOUBLE PRECISION, EXTERNAL :: factrl,dbfactrl
DOUBLE PRECISION, EXTERNAL :: DPoten, Potent
!
! Intrinsic functions
REAL, INTRINSIC :: REAL
DOUBLE PRECISION, INTRINSIC :: DBLE,DOT_PRODUCT,AIMAG
DOUBLE COMPLEX, INTRINSIC :: CONJG
!
! Variables used passed through subroutines
DOUBLE PRECISION, ALLOCATABLE :: rvec(:),Er(:,:), Vr(:)

!!! Universal constants
Pi = 4.0d0*atan2(1.0d0,1.0d0)

!! Parameter specification section

step = 0.01d0

hm = lim + 1!         ! hm must stay below 21 (-roots of Hermite max.)
!                       ! -odd number is better-
p = mass*v
lim1 = lim + 1


!! Allocation section

Allocate(c(0:n),KE(0:n,0:n),PE(0:n,0:n),Ham(0:n,0:n),cH(0:lim,0:lim), &
 & Amtrx(0:n,0:n,0:lim),DpT(0:n,0:n),DaT(0:n,0:n),DaV(0:n,0:n), &
 & nrg(0:n),Vr(0:lim1),rvec(0:lim1),Er(0:lim1,0:lim1), DbT(0:n,0:n), &
 & Vrsave(0:lim1),Potvar(0:lim1), STAT = err)

IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


!! Creating file to save data

OPEN (UNIT=11,FILE='pstn.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err)
OPEN (UNIT=12,FILE='nrg.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err1)
IF ((err.ne.0).OR.(err1.ne.0)) THEN
WRITE(*,*) 'Error opening file'
 STOP
ENDIF

WRITE(UNIT=12,FMT='(A14)') 'Energy levels:'

WRITE(*,*) '-Initial parameters-'
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameter',a,' Phase ', b
WRITE(*,*) '# basis functions:',n, '   PE approx.:',lim


!! Recursive term for Potential energy matrix

cH(0,0) = 1.0d0
DO i=1,lim
 cH(i,0) = 0.0d0
 DO l=1,i/2
  cH(i,0) = cH(i,0) + (-1.0d0)**(l-1)*2.0d0**(l-2) &
             *(dbfactrl(2*l-1)/DBLE(l)) * cH(i-1,2*l-1)
 ENDDO
 cH(i,0) = DBLE(i)*cH(i,0)

 DO j=1,i
  cH(i,j) = (DBLE(i)/(2.0d0*DBLE(j)))*cH(i-1,j-1)
 ENDDO
ENDDO


!! Values of the Hermite integrals

DO k=0,n
 DO l=0,n
  DO i=0,lim
   Amtrx(k,l,i) = 0.0d0
   DO j=0,i
    IF (MODULO((k+l+j),2).eq.0) THEN
     IF ((abs(k-l).le.j).AND.(abs(j-k).le.l).AND.(abs(l-j).le.k)) THEN
      s = (k+l+j)/2
      Amtrx(k,l,i) = Amtrx(k,l,i) + cH(i,j)*(2.0d0**s*sqrt(Pi)*factrl(k) &
         *factrl(l)*factrl(j))/(factrl(s-k)*factrl(s-l)*factrl(s-j)) 
     ENDIF
    ENDIF
   ENDDO
   Amtrx(k,l,i) = 1.0d0 / (factrl(i)*sqrt(2.0d0**(k+l)*factrl(k)*factrl(l) &
         *Pi)) * Amtrx(k,l,i)
  ENDDO
 ENDDO
ENDDO


!!! Initialization of loop !!!

DO
loop = loop + 1

! WRITE(*,*) 'Kinetic energy matrix'
DO k=0,n
 DO l=0,n
  IF (k.eq.l) THEN
   KE(k,l) = -(p**2)/(hbar**2)-((b**2+a**2)/a)*(2.0d0*DBLE(k)+1.0d0)

  ELSE IF (k.eq.(l+1)) THEN
   KE(k,l) = (2.0d0*p/hbar)*sqrt(DBLE(k))*(-b-(0,1)*a)/sqrt(a)

  ELSE IF (k.eq.(l-1)) THEN
   KE(k,l) = (2.0d0*p/hbar)*sqrt(DBLE(k)+1.0d0)*(-b+(0,1)*a)/sqrt(a)

  ELSE IF (k.eq.(l+2)) THEN
   KE(k,l) = sqrt(DBLE(k)*(DBLE(k)-1.0d0))*(((a**2-b**2)/a) - &
 &  2.0d0*(0,1)*b)

  ELSE IF (k.eq.(l-2)) THEN
   KE(k,l) = sqrt((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))*(((a**2-b**2)/a) + &
 &  2.0d0*(0,1)*b)

  ELSE
   KE(k,l) = 0.0d0

  ENDIF
   KE(k,l) = -(hbar**2/(2.0d0*mass))*KE(k,l)

 ENDDO
ENDDO

!!! Potential Energy matrix

DO i=0,lim+1
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

DO k=0,n
 DO l=0,k
  IF (Ham(k,l).ne.CONJG(Ham(l,k))) THEN
   WRITE(*,*) 'Hamiltonian is not Hermitian'
   WRITE(*,*) k,l,Ham(k,l),CONJG(Ham(l,k))
   WRITE(*,*) 'loop = ', loop
   write(*,*) r,p,a,b
   STOP
  ENDIF
 ENDDO
ENDDO

!!! Diagonalization of Hamiltonian matrix to find initial coeff's
IF (loop.eq.1) THEN

 DO i =0,n
  c(i) = (0.0d0,0.0d0)
 ENDDO
 c(0) = (1.0d0,0.0d0)
ENDIF


!!! Normalization of the coefficients

 CALL NORM(c)

opt = .false.
CALL Hamdiag(Ham,opt,nrg,c)

!!! Expectation values for position, momentum and energy

rksum = 0.0d0; pksum = 0.0d0; pksum2 = 0.0d0
DO l=0,n-1
 rksum = rksum + sqrt(DBLE(l)+1.0d0)*DBLE(CONJG(c(l))*c(l+1))
 pksum = pksum + sqrt(DBLE(l)+1.0d0)*DBLE(CONJG(c(l+1))*c(l))
 pksum2 = pksum2 + sqrt(DBLE(l)+1.0d0)*AIMAG(CONJG(c(l))*c(l+1))
ENDDO

rxvl = r + rksum/sqrt(a)
pxvl = p + 2.0d0*b*pksum/sqrt(a) + 2.0d0*sqrt(a)*pksum2


!!! Derivative of Energy wrt R

DrE = 0.0d0
DO k=0,n
 DO l=0,n
  DO i=0,lim
   DrE = DrE + (CONJG(c(k))*c(l))*(Amtrx(k,l,i)*Potvar(i+1) / &
    ((2.0d0*a)**(DBLE(i)/2.0d0)))
  ENDDO
 ENDDO
ENDDO

!!! Derivative of Energy wrt p

DpE = 0.0d0
DO k=0,n
 DO l=0,n
  IF (k.eq.l) THEN
   DpT(k,l) = p / mass  
  ELSE IF (k.eq.(l-1)) THEN
   DpT(k,l) = -(1.0d0/mass)*sqrt(DBLE(k)+1.0d0)*(-b+(0.0d0,1.0d0)*a)/sqrt(a)
  ELSE IF (k.eq.(l+1)) THEN
   DpT(k,l) = -(1.0d0/mass)*sqrt(DBLE(k))*(-b-(0.0d0,1.0d0)*a)/sqrt(a)
  ELSE
   DpT(k,l) = 0.0d0
  ENDIF
  DpE = DpE + (CONJG(c(k))*c(l))*DpT(k,l)
 ENDDO
ENDDO


!!! Derivative of Energy wrt a

DaE = 0.0d0
DO k=0,n
 DO l=0,n
  IF (k.eq.l) THEN
   DaT(k,l) = (2.0d0*DBLE(k)+1.0d0)*(1.0d0 - b**2/a**2)/(2.0d0*mass)
  ELSE IF (k.eq.(l-1)) THEN
   DaT(k,l) = -p*SQRT(DBLE(k)+1.0d0)*(b/a-(0.0d0,1.0d0))/(2.0d0*mass*SQRT(a))
  ELSE IF (k.eq.(l+1)) THEN
   DaT(k,l) = -p*SQRT(DBLE(k))*(b/a + (0,1)) / (2.0d0*mass*SQRT(a))
  ELSE IF (k.eq.(l-2)) THEN
   DaT(k,l) = -sqrt((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))*(1.0d0+b**2/a**2)/(2.0d0*mass)
  ELSE IF (k.eq.(l+2)) THEN
   DaT(k,l) = -sqrt((DBLE(k)-1.0d0)*DBLE(k))*(1.0d0+b**2/a**2)/(2.0d0*mass)
  ELSE
   DaT(k,l) = 0.0d0
  ENDIF

  DaV(k,l) = 0.0d0
  DO i=1,lim
   DaV(k,l) = DaV(k,l)-(Amtrx(k,l,i)*DBLE(i)*Potvar(i)) / &
    ((2.0d0*a)**((DBLE(i)+2.0d0)/2.0d0))
  ENDDO

  DaE = DaE + CONJG(c(k))*c(l)*(DaT(k,l)+DaV(k,l))
 ENDDO
ENDDO


!!! Derivative of Energy wrt b

DbE = 0.0d0
DO k=0,n
 DO l=0,n
  IF (k.eq.l) THEN
   DbT(k,l) = (2.0d0*DBLE(k) + 1.0d0)*b / (mass*a)
  ELSE IF (k.eq.(l-1)) THEN
   DbT(k,l) = p*SQRT(DBLE(k)+1.0d0) / (mass*SQRT(a))
  ELSE IF (k.eq.(l+1)) THEN
   DbT(k,l) = p*SQRT(DBLE(k)) / (mass*SQRT(a))
  ELSE IF (k.eq.(l-2)) THEN
   DbT(k,l) = -sqrt((DBLE(k)+1.0d0)*(DBLE(k)+2.0d0))* &
		(2.0d0*(0.0d0,1.0d0)-2.0d0*b/a) / (2.0d0*mass)
  ELSE IF (k.eq.(l+2)) THEN
   DbT(k,l) = sqrt(DBLE(k)*(DBLE(k)-1.0d0)) * &
	& (2.0d0*(0.0d0,1.0d0)+2.0d0*b/a) / (2.0d0*mass)
  ELSE
   DbT(k,l) = 0.0d0
  ENDIF

  DbE = DbE + CONJG(c(k))*c(l)*DbT(k,l)
 ENDDO
ENDDO


!!! Gradient method - move variable in direction of greatest descent

grad = (/DBLE(DrE),DBLE(DpE),DBLE(DaE),DBLE(DbE)/)

unit = sqrt(DOT_PRODUCT(grad,grad))
ugrad = grad / unit

max = 0
DO
10 max = max + 1
 IF (max.gt.3) EXIT
 DO j=1,4-max
  IF (abs(ugrad(max)).lt.abs(ugrad(max+j))) THEN
   goto 10
  ELSE
   CYCLE
  ENDIF
 ENDDO
 EXIT
ENDDO

CALL RANDOM_NUMBER(randnum)

 IF ((ABS(DrE)+ABS(DpE)+ABS(DaE)+ABS(DbE)).lt.step) THEN
  step = randnum*step
 ENDIF

! IF (MODULO(loop,50000).eq.0) step = 0.01d0*step

!SELECT CASE(max) ! Positive sign .... or negative !
!CASE(1)
 r = r - step*ugrad(1)
!CASE(2)
 p = p - step*ugrad(2)
!CASE(3)
 a = a - step*ugrad(3)
!CASE(4)
 b = b - step*ugrad(4)
!END SELECT


Eexp = 0.0d0
DO i=0,n
  Eexp = Eexp + CONJG(c(i))*c(i)*nrg(i)
ENDDO


IF (MODULO(loop,1001).eq.0) THEN
 WRITE(UNIT=11,FMT='(A2,I8,A2)') '**',loop,'**'
 WRITE(UNIT=11,FMT='(A4,F14.10)') 'Pos:',r
 WRITE(UNIT=11,FMT='(A6,F14.10)') 'Momen:',p
 WRITE(UNIT=11,FMT='(A6,F14.10)') 'Width:',a
 WRITE(UNIT=11,FMT='(A6,F14.10)') 'b-para:',b
 WRITE(UNIT=11,FMT='(A4,F16.12)') 'DrE=',REAL(DrE)
 WRITE(UNIT=11,FMT='(A4,F16.12)') 'DpE=',REAL(DpE)
 WRITE(UNIT=11,FMT='(A4,F16.12)') 'DaE=',REAL(DaE)
 WRITE(UNIT=11,FMT='(A4,F16.12)') 'DbE=',REAL(DbE)
 WRITE(UNIT=11,FMT='(F16.12)') ABS(DrE)+ABS(DpE)+ABS(DaE)+ABS(DbE)
 WRITE(UNIT=11,FMT='(A5,F16.12)') 'step=',step
! WRITE(UNIT=11,FMT='(A5,F16.12)') 'rxvl=',rxvl
! WRITE(UNIT=11,FMT='(A5,F16.12)') 'pxvl=',pxvl

 WRITE(UNIT=12,FMT='(A2,I8,A2)') '**',loop,'**'
 WRITE(UNIT=12,FMT='(A5,F14.10)') 'step:',step
 DO i=0,n
  WRITE(UNIT=12,FMT='(I2,A3,F12.7,F12.7)') i,' = ',nrg(i), &
	& Eexp
 ENDDO

ENDIF

!IF (MODULO(loop,1001).eq.0) THEN
! write(*,*) '**',loop,'**'
! write(*,*) nrg(:)*627.51d0,Eexp*627.51d0
!ENDIF

IF (((ABS(DrE)+ABS(DpE)+ABS(DaE)+ABS(DbE)).lt.1E-12).OR.(loop.gt.1E7)) EXIT


!! Reset all the variables
DO k=0,n
 DO l=0,n
  KE(k,l)=(0.0d0,0.0d0);PE(k,l)=(0.0d0,0.0d0);Ham(k,l)=(0.0d0,0.0d0)
  DpT(k,l)=(0.0d0,0.0d0);DaT(k,l)=(0.0d0,0.0d0);DaV(k,l)=(0.0d0,0.0d0)
 ENDDO
ENDDO


ENDDO !!!!!!!

CLOSE(11); CLOSE(12)

WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameter',a,' Phase ', b
WRITE(*,*) 'loop',loop
WRITE(*,*) 'END PROGRAM'

END PROGRAM
