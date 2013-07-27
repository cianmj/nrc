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

 SUBROUTINE CALCBIL(calc,Amtrx,Ham,dr,dp,da,db,c,dc)
  USE FNVAR
  IMPLICIT NONE
  LOGICAL :: calc
  DOUBLE PRECISION, INTENT(in), OPTIONAL :: dr,dp,da,db
  DOUBLE PRECISION, INTENT(in) :: Amtrx(0:n,0:n,0:lim)
  DOUBLE COMPLEX, INTENT(in), OPTIONAL :: c(0:n)
  DOUBLE COMPLEX, INTENT(out), OPTIONAL :: dc(0:n)
  DOUBLE COMPLEX, INTENT(out) :: Ham(0:n,0:n)
 END SUBROUTINE

 SUBROUTINE NUMRND(value)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(inout) :: value
 END SUBROUTINE

 SUBROUTINE EVALWF(prog,c,Eexp)
  USE FNVAR
  IMPLICIT NONE
  INTEGER, INTENT(in) :: prog
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
INTEGER, PARAMETER :: prog = 33
DOUBLE PRECISION :: acc,r0,v0,p0,a0,b0,r1,v1,rtemp,QC1,QC2, coeffsum
DOUBLE PRECISION :: dr, dp, da, db
DOUBLE COMPLEX :: w, w1, w0, dw
LOGICAL :: opt,quantum,calc, wrt = .true.
!
DOUBLE PRECISION, ALLOCATABLE :: cH(:,:),Amtrx(:,:,:),nrg(:)
DOUBLE COMPLEX, ALLOCATABLE :: c(:), dc(:), c0(:)
DOUBLE COMPLEX, ALLOCATABLE :: Ham(:,:)
DOUBLE PRECISION, ALLOCATABLE :: zm(:),Em(:,:), WTm(:)
DOUBLE PRECISION, ALLOCATABLE :: Mint(:,:,:),Mnki(:,:,:),Snm(:,:),Ssum(:,:)
!
DOUBLE PRECISION :: DPr, DPr1, DDPr, DDPr1, DPtemp
DOUBLE PRECISION :: xx, psave(1:4), dsave(1:4)
DOUBLE PRECISION :: rsum, psum, psum2, rexp, pexp, Eexp
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl, dbfactrl, DPoten, Potential
!
! Intrinsic functions
REAL, INTRINSIC :: REAL, ABS
DOUBLE PRECISION, INTRINSIC :: DBLE, DIMAG
DOUBLE COMPLEX, INTRINSIC :: CONJG, DCMPLX, SQRT
!
! Variable passed through other subroutines
DOUBLE PRECISION, ALLOCATABLE :: rvec(:),Er(:,:),Vr(:)

!!! Universal constants
Pi = 4.0d0*atan2(1.0d0,1.0d0)

!!! Parameter specification section

hm = lim + 1!		! hm must stay below 21 (-roots of Hermite max.)
!m = MAX(n+1,9)
m = 7

p = mass*v
lim1 = lim + 1

quantum = .false. !		! include quantum corrections ?
!quantum = .true.
limit = lim1 !			! order of the "	"

!! Allocation section

Allocate(cH(0:lim,0:lim), Amtrx(0:n,0:n,0:lim), c(0:n), c0(0:n), &
 & Ham(0:n,0:n), dc(0:n), nrg(0:n), rvec(0:lim1), Er(0:lim1,0:lim1), &
 & Vr(0:lim1), Vrsave(0:lim1), zm(0:m), Em(0:m,0:m), WTm(0:m), &
 & Mint(0:n,0:n,0:limit+2),Mnki(0:n,0:n,0:limit+2),Snm(1:2,1:limit), &
 & Ssum(1:2,1:limit), STAT = err)
IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


!! Creating file to save data

OPEN (UNIT=31,FILE='databil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err)
OPEN (UNIT=32,FILE='plotbil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err1)
OPEN (UNIT=33,FILE='densitybil.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err2)
OPEN (UNIT=34,FILE='potential.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err2)
IF ((err.ne.0).OR.(err1.ne.0).OR.(err2.ne.0)) THEN
WRITE(*,*) 'Error opening file'
 STOP
ENDIF

!
xx = -10.01d0
Do i = 1,20000
 xx = xx + 0.01d0
 WRITE(UNIT=34,FMT='(F16.2,F16.5)') xx,Potential(xx)
ENDDO
!

WRITE(UNIT=31,FMT='(A)') 'Expectation values: pos,mom,nrg'
WRITE(UNIT=32,FMT='(A)') 'Values of pos,mom,width,phase'

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

p = mass*v
w = DCMPLX(b,a)


!! Diagonalization of Hamiltonian matrix to find initial coefficients

IF (loop.eq.1) THEN

DO k=0,n
 c(k) = (0.0d0,0.0d0)
! c(k) = (1.0d0,0.0d0)
ENDDO
c(:) = (1.0d0,0.0d0)

 calc = .false.
 CALL CALCBIL(calc,Amtrx,Ham)

 opt = .false.
 IF (opt) WRITE(*,*) 'Finding initial coefficients.'
 CALL Hamdiag(Ham,opt,nrg,c)

 CALL NORM(c)

ENDIF


!!! Quantum corrections

IF (quantum) THEN
 IF (loop.eq.1) WRITE(*,*) 'Including quantum corrections.'
 
 Mint(:,:,:)= 0.0d0; Mnki(:,:,:) = 0.0d0
 DO i=0,n
  DO j=0,n
   DO l = 0,limit+2
    DO k=0,m-1
     Mint(i,j,l)=Mint(i,j,l)+WTm(k)*Em(i,k)*Em(j,k) * &
	& (SQRT(hbar/(2.0d0*a))*zm(k))**l
    ENDDO
   Mnki(i,j,l)=(1.0d0/(factrl(i)*factrl(j)*2.0d0**(i+j)*SQRT(Pi) )) * & 
		& Mint(i,j,l)*SQRT(2.0d0*a/hbar)
   ENDDO
  ENDDO
 ENDDO

 Snm(:,:) = 0.0d0; Ssum(:,:) = 0.0d0
 DO i=1,2
  DO j=1,limit
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

!write(*,*) 'Snm',Snm(:,:)

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


!!! Calculation of propagation equations
calc = .true.

r0 = r; v0 = v; p0 = p; a0 = a; b0 = 0; w0 = w; c0(:) = c(:)

DPr = DPoten(1) + QC1
DDPr = DPoten(2) + QC2

dr = v
dp = -DPr
dw = - DCMPLX(2.0d0*(w*w)/mass + 0.5d0*DDPr)
da = DIMAG(dw); db = DBLE(dw)

CALL CALCBIL(calc,Amtrx,Ham,dr,dp,da,db,c,dc)

acc = dp/mass
r = r + dr*dt/2.0d0
v = v + acc*dt/2.0d0
p = mass*v
w = w + dw*dt/2.0d0
a = DIMAG(w); b = DBLE(w)

c(:) = c(:) + dc(:)*dt/2.0d0

CALL NORM(c)

!write(*,*) 'half-time step'
!write(*,*) r,p,a,b
!write(*,*) c(:)

! - !

DPr = DPoten(1) + QC1
DDPr = DPoten(2) + QC2

dr = v
dp = -DPr
dw = - DCMPLX(2.0d0*(w*w)/mass + 0.5d0*DDPr)
da = DIMAG(dw); db = DBLE(dw)

CALL CALCBIL(calc,Amtrx,Ham,dr,dp,da,db,c,dc)

acc = dp/mass
r = r0 + dr*dt
v = v0 + acc*dt
p = mass*v
w = w0 + dw*dt
a = DIMAG(w); b = DBLE(w)

c(:) = c0(:) + dc(:)*dt


DO i = 0,n
 CALL NUMRND(c(i))
ENDDO
CALL NORM(c)

! Check for normalization
 coeffsum = 0.0d0
 DO i=0,n
  coeffsum = coeffsum + CONJG(c(i))*c(i)
 ENDDO
 IF ((ABS(coeffsum-1.0d0).gt.epsilon)) WRITE(*,*) 'Coefficients not normalized'



!write(*,*) 'billing'
!write(*,*) 'time = ',time
!write(*,*) r,p,a,b
!write(*,*) c(:)
!write(*,*) DPoten(1), QC1
!write(*,*) DPoten(2), QC2
!read(*,*)

!!! Expectation values for position, momentum and energy

rsum = 0.0d0; psum = 0.0d0; psum2 = 0.0d0
DO l=0,n-1
 rsum = rsum + SQRT(DBLE(l)+1.0d0)*DBLE(CONJG(c(l))*c(l+1))
 psum = psum + SQRT(DBLE(l)+1.0d0)*DBLE(CONJG(c(l+1))*c(l))
 psum2 = psum2 + SQRT(DBLE(l)+1.0d0)*DIMAG(CONJG(c(l))*c(l+1))
ENDDO

rexp = r + rsum/SQRT(a)
pexp = p + 2.0d0*b*psum/SQRT(a) + 2.0d0*SQRT(a)*psum2

Eexp = 0.0d0
DO i=0,n
 DO j=0,n
   Eexp = Eexp + CONJG(c(i))*c(j)*Ham(i,j)
 ENDDO
ENDDO


IF (MODULO(INT(time),10).eq.1) wrt = .true.

IF ((MODULO(INT(time),10).eq.0).AND.(wrt)) THEN
 WRITE(UNIT=31,FMT='(F12.7,F12.7,F12.7,F12.7)') time,rexp,pexp,Eexp
 WRITE(UNIT=32,FMT='(F12.7,F12.7,F12.7,F12.7,F12.7)') time,r,p,a,b
 CALL EVALWF(prog,c,Eexp)
 wrt = .false.
ENDIF

!IF (r.lt.0.0) THEN
! WRITE(*,*) 'r less than 0.0'
! WRITE(*,*) 'r=',r
! EXIT
!ENDIF

IF (time.gt.runtime) EXIT


!! Reset variables to zero ::

DO k=0,n
 DO l=0,n
   Ham(k,l) = 0.0d0
 ENDDO
ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO !	! ! !Ending of propagation loop ! ! ! !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CLOSE(31); CLOSE(32); CLOSE(33); CLOSE(34)

WRITE(*,*) 'loop',loop
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameters:',a,b
WRITE(*,*) 'Coefficients'
WRITE(*,*) (k,c(k),k=0,n)

WRITE(*,*) 'Program endded'

END PROGRAM
