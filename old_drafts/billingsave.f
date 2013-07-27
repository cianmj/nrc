PROGRAM BILLING

USE FNVAR
USE VARSAVE
USE ROOT

IMPLICIT NONE

INTERFACE
 SUBROUTINE Hamdiag(Ham,opt,nrg,c)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE COMPLEX, intent(in) :: Ham(0:n,0:n)
  LOGICAL, intent(in) :: opt
  DOUBLE PRECISION, intent(out) :: nrg(0:n)
  DOUBLE COMPLEX, intent(out), OPTIONAL :: c(0:n)
 END SUBROUTINE

 SUBROUTINE Hzeros(m,zm,Em)
  IMPLICIT NONE
  INTEGER, intent(in) :: m
  DOUBLE PRECISION, intent(out) :: zm(0:m),Em(0:m,0:m)
 END SUBROUTINE

 SUBROUTINE NORM(c)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(inout) :: c(0:n)
 END SUBROUTINE

 SUBROUTINE EVALWF(c,Eexp)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: Eexp
  DOUBLE COMPLEX, INTENT(in) :: c(0:n)
 END SUBROUTINE

END INTERFACE

! Input parameters
INTEGER :: m
DOUBLE PRECISION :: v, time
!
INTEGER :: k,l,i,j,pp
INTEGER :: s,err,err1,err2,limit
DOUBLE PRECISION :: acc,r1,v1,p1,a1,b1,rtemp,QC1,QC2
DOUBLE COMPLEX :: w,w1
LOGICAL :: opt,quantum
!
DOUBLE PRECISION, ALLOCATABLE :: cH(:,:),Amtrx(:,:,:),nrg(:),Potvar(:)
DOUBLE COMPLEX, ALLOCATABLE :: c(:),dc(:), Wa(:)
DOUBLE COMPLEX, ALLOCATABLE :: KE(:,:),PE(:,:),Ham(:,:)
DOUBLE PRECISION, ALLOCATABLE :: EPoten(:), WInt(:,:),Wnn(:,:)
DOUBLE PRECISION, ALLOCATABLE :: zm(:),Em(:,:), WTm(:)
DOUBLE PRECISION, ALLOCATABLE :: Mint(:,:,:),Mnki(:,:,:),Snm(:,:),Ssum(:,:)
!
DOUBLE PRECISION :: DPr,DDPr,DPtemp,xx
DOUBLE PRECISION :: rsum,psum,psum2,rexp,pexp,Eexp
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl, dbfactrl, DPoten, Potential
!
! Intrinsic functions
REAL, INTRINSIC :: REAL, ABS
DOUBLE PRECISION, INTRINSIC :: DBLE, AIMAG
DOUBLE COMPLEX, INTRINSIC :: CONJG, CMPLX, SQRT
!
! Variable passed through other subroutines
DOUBLE PRECISION, ALLOCATABLE :: rvec(:),Er(:,:),Vr(:)

!!! Universal constants
Pi = 4.0d0*atan2(1.0d0,1.0d0)

!!! Parameter specification section

quantum = .false. !		! include quantum corrections ?
!quantum = .true.
!limit = 10 !			! order of the "	"

hm = lim + 1!		! hm must stay below 21 (-roots of Hermite max.)
m = n + 1

p = mass*v
lim1 = lim + 1

!! Allocation section

Allocate(cH(0:lim,0:lim),Amtrx(0:n,0:n,0:lim),c(0:n), Wnn(0:n,0:n), &
 & Ham(0:n,0:n),dc(0:n),nrg(0:n),EPoten(0:m-1),WInt(0:n,0:n), & 
 & rvec(0:lim1),Er(0:lim1,0:lim1),Vr(0:lim1),Vrsave(0:lim1), &
 & Potvar(0:lim),Wa(0:n),zm(0:m),Em(0:m,0:m),WTm(0:m), &
 & Mint(0:n,0:n,0:2*n), Mnki(0:n,0:n,0:2*n),Snm(0:n,0:n),Ssum(0:n,0:n), &
 & STAT = err)
IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


!! Creating file to save data

OPEN (UNIT=21,FILE='databil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err)
OPEN (UNIT=22,FILE='plotbil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err1)
OPEN (UNIT=23,FILE='densitybil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err2)
OPEN (UNIT=24,FILE='potential.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err2)
IF ((err.ne.0).OR.(err1.ne.0).OR.(err2.ne.0)) THEN
WRITE(*,*) 'Error opening file'
 STOP
ENDIF

xx = -10.01d0
Do i = 1,2000
 xx = xx + 0.01d0
 WRITE(UNIT=24,FMT='(F16.5,F16.8)') xx,Potential(xx)
ENDDO

WRITE(UNIT=21,FMT='(A)') 'Expectation values: pos,mom,nrg'
WRITE(UNIT=22,FMT='(A)') 'Values of pos,mom,width,phase'

WRITE(*,*) '-Initial parameters-'
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameters:',a,b
WRITE(*,*) '# basis functions:',n+1,' n=',n,'   PE approx.:',lim


!!! Recursive term for Potential energy matrix

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


!!! Values of the Hermite integrals

DO k=0,n
 DO l=0,n
  DO i=0,lim
   Amtrx(k,l,i) = 0.0d0
   DO j=0,i
    IF (MODULO((k+l+j),2).eq.0) THEN
     IF ((abs(k-l).le.j).AND.(abs(j-k).le.l).AND.(abs(l-j).le.k)) THEN
      s = (k+l+j)/2
      Amtrx(k,l,i) = Amtrx(k,l,i) + cH(i,j)*(2.0d0**s*SQRT(Pi)*factrl(k) &
         *factrl(l)*factrl(j))/(factrl(s-k)*factrl(s-l)*factrl(s-j)) 
     ENDIF
    ENDIF
   ENDDO
   Amtrx(k,l,i) = 1.0d0 / (factrl(i)*SQRT(2.0d0**(k+l)*factrl(k)*factrl(l) &
         *Pi)) * Amtrx(k,l,i)
  ENDDO
 ENDDO
ENDDO

DEALLOCATE (cH)

!!! Finding roots of mth order Hermite polynomials and corrpdn. weights

CALL Hzeros(m,zm,Em)

 DO j=0,m-1
  WTm(j)= 2.0d0**(m-1) * factrl(m) * SQRT(Pi) / (DBLE(m)**2 * (Em(m-1,j))**2)
  IF (Em(m-1,j).eq.0) then
   WRITE(*,*) 'Error finding weighting function'
   STOP
  ENDIF
 ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Beginning of main iterative of loop !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO
loop = loop + 1
time = DBLE(loop)*dt

Allocate(KE(0:n,0:n),PE(0:n,0:n), STAT = err)
IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF

p = mass*v
w = CMPLX(b,a)

!!! Kinetic Energy matrix

! WRITE(*,*) 'Kinetic energy matrix'
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
    WRITE(UNIT=21,FMT='(A)') 'Hamiltonian matrix is not Hermitian'
    WRITE(UNIT=21,FMT='(A4,I8)') 'loop',loop
    STOP
  ENDIF
 ENDDO
ENDDO


!! Diagonalization of Hamiltonian matrix to find initial coefficients

IF (loop.eq.1) THEN

DO k=0,n
 c(k) = (0.0d0,0.0d0)
! c(k) = (1.0d0,0.0d0)
ENDDO
c(0) = (1.0d0,0.0d0)

 opt = .false.
 IF (opt) WRITE(*,*) 'Finding initial coefficients.'
 CALL Hamdiag(Ham,opt,nrg,c)

 CALL NORM(c)

ENDIF

DEALLOCATE(KE,PE)


!!! Quantum corrections

IF (quantum) THEN
 IF (loop.lt.2) WRITE(*,*) 'Including quantum corrections.'
 
 Mint(:,:,:)= 0.0d0; Mnki(:,:,:) = 0.0d0
 DO i=0,n
  DO j=0,n
   DO l = 0,2*n
    DO k=0,m-1
     Mint(i,j,l)=Mint(i,j,l)+WTm(k)*Em(i,k)*Em(j,k) * &
	& (SQRT(hbar/(2.0d0*a))*zm(k))**l
    ENDDO
   Mnki(i,j,l)=(1.0d0/(factrl(i)*factrl(j)*2.0d0**(i+j)))*Mint(i,j,l)
! /SQRT(Pi) !
   ENDDO
  ENDDO
 ENDDO

 Snm(:,:) = 0.0d0; Ssum(:,:) = 0.0d0
 DO i=0,n
  DO j=0,n
   DO k=0,n
    DO l=0,n
     DO pp=0,n
       Ssum(i,j) = Ssum(i,j) + ( Mnki(k,pp,i)*Mnki(pp,l,j)-Mnki(k,l,i+j) )
     ENDDO
     Snm(i,j) = Snm(i,j) + ( CONJG(c(k))*c(l) * Ssum(i,j) )
    ENDDO
   ENDDO
  ENDDO
 ENDDO

 QC1 = 0.0d0; QC2 = 0.0d0
 DO i = 3,limit
  QC1 = QC1 + (DPoten(i)/factrl(i))*(Snm(2,2)*Snm(1,i)-Snm(1,2)*Snm(2,i)) / &
		& (Snm(1,1)*Snm(2,2)-Snm(1,2)*Snm(2,1))
  QC2 = QC2 + (DPoten(i)/factrl(i))*(Snm(1,1)*Snm(2,i)-Snm(2,1)*Snm(1,i)) / &
                & (Snm(1,1)*Snm(2,2)-Snm(1,2)*Snm(2,1))
 ENDDO

ELSE

 QC1 = 0.0d0
 QC2 = 0.0d0

ENDIF


write(*,*) 'billingbck'
write(*,*) 'time = ',time
write(*,*) r,p,a,b
write(*,*) c(:)
read(*,*)


!!! Calculation of propagation equations

If (loop.eq.1) THEN
 DPr = DPoten(1) + QC1
 DDPr = DPoten(2) + QC2
ENDIF

acc = - DPr/mass
v1 = v + acc*dt/2.0d0
r1 = r + v1*dt/2.0d0
r = r1 + v1*dt/2.0d0
DPr = DPoten(1) + QC1
acc = -DPr/mass
v = v1 + acc*dt/2.0d0

w1 = w - DCMPLX(2.0d0*(w*w)/mass + 0.5d0*DDPr)*(dt/2.0d0)
DDPr = DPoten(2) + QC2
w = w1 - DCMPLX(2.0d0*(w1*w1)/mass + 0.5d0*DDPr)*(dt/2.0d0)

CALL NUMRND(w)

p1 = mass*v1
a1 = AIMAG(w1); b1 = DBLE(w1)
a = AIMAG(w); b = DBLE(w)


!!! Calculation of the effective potential (with 1/2 time-step variables)

rtemp = r
DO i=0,m-1
 r = SQRT(hbar/(2.0d0*a1))*zm(i) + r1
 DPtemp = DPoten(0)
 r = r1
 EPoten(i)=DPtemp-DPoten(0)-(DPoten(1)+QC1)*(SQRT(hbar/(2.0d0*a1))*zm(i)) &
  & - 0.5d0*(DPoten(2)+QC2)*(SQRT(hbar/(2.0d0*a1))*zm(i))**2
ENDDO
 r = rtemp
       
WInt(:,:)= 0.0d0; Wnn(:,:)= 0.0d0
DO i=0,n
 DO j=0,n
  DO k=0,m-1
   WInt(i,j)= WInt(i,j) + WTm(k)*Em(i,k)*Em(j,k)*EPoten(k)
  ENDDO
  Wnn(i,j)= (1.0d0/(factrl(i)*factrl(j)*2.0d0**(i+j)*SQRT(Pi))) * &
	& SQRT(hbar/(2.0d0*a1))*WInt(i,j)
 ENDDO
ENDDO

!! Calculation of the summation term of the effective potential
!! over the coefficients; as seen in the coeff. propa'n. equation !

Wa(:) = (0.0d0,0.0d0)
DO i = 0,n
  DO k = 0,n
    Wa(i) = Wa(i) + Wnn(i,k)*c(k)
  ENDDO
ENDDO

!!! Propagation of coefficients

DO i=0,n
c(i) = c(i)+((0.0d0,-1.0d0)/hbar)*(Wa(i) + c(i)*(p1**2/(2.0d0*mass) &
 & + (2.0d0*DBLE(i)+1.0d0)*hbar*a1/mass))*dt
ENDDO

DO i = 0,n
 CALL NUMRND(c(i))
ENDDO

a = AIMAG(w); b = DBLE(w)
p = mass*v

CALL NORM(c)

!!! Expectation values for position, momentum and energy

rsum = 0.0d0; psum = 0.0d0; psum2 = 0.0d0
DO l=0,n-1
 rsum = rsum + SQRT(DBLE(l)+1.0d0)*DBLE(CONJG(c(l))*c(l+1))
 psum = psum + SQRT(DBLE(l)+1.0d0)*DBLE(CONJG(c(l+1))*c(l))
 psum2 = psum2 + SQRT(DBLE(l)+1.0d0)*AIMAG(CONJG(c(l))*c(l+1))
ENDDO

rexp = r + rsum/SQRT(a)
pexp = p + 2.0d0*b*psum/SQRT(a) + 2.0d0*SQRT(a)*psum2

Eexp = 0.0d0
DO i=0,n
 DO j=0,n
   Eexp = Eexp + CONJG(c(i))*c(j)*Ham(i,j)
 ENDDO
ENDDO


IF (MODULO(INT(time/dt),10).eq.0) THEN
 WRITE(UNIT=21,FMT='(F12.7,F12.7,F12.7,F12.7)') time,rexp,pexp,Eexp
 WRITE(UNIT=22,FMT='(F12.7,F12.7,F12.7,F12.7,F12.7)') time,r,p,a,b

 CALL EVALWF(c,Eexp)
ENDIF

IF ((r.lt.0.0).OR.(a.lt.0.0)) THEN
 WRITE(*,*) 'r or a less than 0.0'
 WRITE(*,*) 'r=',r,'  a=',a
 EXIT
ENDIF
IF (loop.gt.3E4) EXIT


!! Reset variables to zero ::

DO k=0,n
 DO l=0,n
   Ham(k,l) = 0.0d0
 ENDDO
ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO !	! ! !Ending of propagation loop ! ! ! !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CLOSE(21); CLOSE(22); CLOSE(23); CLOSE(24)

WRITE(*,*) 'loop',loop
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameters:',a,b
WRITE(*,*) 'Coefficients'
WRITE(*,*) (k,c(k),k=0,n)

WRITE(*,*) 'Program endded'

END PROGRAM
!
!!!!!!!!!!!
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

