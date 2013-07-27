PROGRAM PROPGTN

USE FNVAR
USE VARSAVE

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

 SUBROUTINE CONDCHECK(Ham)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(in) :: Ham(0:n,0:n)
 END SUBROUTINE

 SUBROUTINE NORM(c)
  USE FNVAR
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(inout) :: c(0:n)
 END SUBROUTINE

 SUBROUTINE CALCRATE(n2,dn,c,Hl,Ham,dr,dp,da,db,dc,expnQ)
  USE FNVAR
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n2,dn
  DOUBLE COMPLEX, INTENT(in) :: c(0:n),Hl(0:dn,0:n),Ham(0:n,0:n)
  DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
  DOUBLE COMPLEX, INTENT(out) :: dc(0:n), expnQ(0:n,0:n)
 END SUBROUTINE

 SUBROUTINE EVALWF(prog,c,Eexp)
  USE FNVAR
  IMPLICIT NONE
  INTEGER, INTENT(in) :: prog
  DOUBLE PRECISION, INTENT(in) :: Eexp
  DOUBLE COMPLEX, INTENT(in) :: c(0:n)
 END SUBROUTINE

END INTERFACE

INTEGER :: k,l,i,j
INTEGER :: s,err,err1,err2,err3,n2,dn
INTEGER, PARAMETER :: prog = 23
DOUBLE PRECISION :: v, time, coeffsum
DOUBLE PRECISION :: dr,dp,da,db,xx
DOUBLE PRECISION :: test,test1,r0,p0,a0,b0,r1,p1,a1,b1
DOUBLE PRECISION :: rsum,psum,psum2,rexp,pexp,Eexp
DOUBLE PRECISION, ALLOCATABLE :: cH(:,:),Amtrx(:,:,:),nrg(:),Potvar(:)
DOUBLE COMPLEX, ALLOCATABLE :: c(:),dc(:),c0(:),c1(:),expnQ(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: KE(:,:),PE(:,:),Ham(:,:),Hl(:,:)
LOGICAL :: opt, wrt = .true.
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl,dbfactrl
DOUBLE PRECISION, EXTERNAL :: DPoten,Potential
!
! Intrinsic functions
REAL, INTRINSIC :: REAL,ABS
DOUBLE PRECISION, INTRINSIC :: DBLE, AIMAG
DOUBLE COMPLEX, INTRINSIC :: CONJG, SQRT, MATMUL
!
! Variable passed through other subroutines
DOUBLE PRECISION, ALLOCATABLE :: AA(:,:),BB(:,:),rvec(:),Er(:,:),Vr(:)
DOUBLE COMPLEX, ALLOCATABLE :: Drkl(:,:),Dpkl(:,:),Dakl(:,:),Dbkl(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: DrLC(:),DpLC(:),DaLC(:),DbLC(:),HlC(:)
DOUBLE COMPLEX, ALLOCATABLE :: HH(:,:)

!!! Universal constants
Pi = 4d0*atan2(1d0,1d0)

!!! Parameter specification section

! lim = expansion order for PE ( must be EVEN !)

p = mass*v
lim1 = lim + 1
hm = lim + 1!		! hm must stay below 21 (-roots of Hermite max.)

n2 = n + 2  ! this will give 4 equations with 4 unknowns !
dn = n2 - n - 1

!! Allocation section

Allocate(cH(0:lim,0:lim),Amtrx(0:n2,0:n,0:lim),c(0:n),expnQ(0:n,0:n), &
 & Ham(0:n,0:n),Hl(0:dn,0:n),Drkl(0:n2,0:n),Dpkl(0:n2,0:n), &
 & Dakl(0:n2,0:n),Dbkl(0:n2,0:n),dc(0:n),nrg(0:n),c0(0:n),c1(0:n), & 
 & DrLC(0:dn),DpLC(0:dn),DaLC(0:dn),DbLC(0:dn),HlC(0:dn),AA(2*dn,4), &
 & BB(2*dn,1),rvec(0:lim1),Er(0:lim1,0:lim1),Vr(0:lim1),Vrsave(0:lim1), &
 & HH(n,n),Potvar(0:lim), STAT = err)

IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


!! Creating file to save data

OPEN (UNIT=21,FILE='data.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
OPEN (UNIT=22,FILE='plot.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=err1)
OPEN (UNIT=23,FILE='density.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=err2)
OPEN (UNIT=24,FILE='potential.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=err3)
IF ((err.ne.0).OR.(err1.ne.0).OR.(err2.ne.0).OR.(err3.ne.0)) THEN
WRITE(*,*) 'Error opening file'
 STOP
ENDIF


WRITE(UNIT=21,FMT='(A)') 'Expectation values: pos,mom,nrg'
WRITE(UNIT=22,FMT='(A)') 'Values of pos,mom,width,phase'

WRITE(*,*) '-Initial parameters-'
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameters:',a,b
WRITE(*,*) '# basis functions:',n+1,' n=',n,'   PE approx.:',lim


!!! For plotting the potential - potential.txt

!xx = -0.01d0
!Do i = 1,2000

xx = -10.01d0
Do i = 1,20000
 xx = xx + 0.01d0
 WRITE(UNIT=24,FMT='(F16.2,F16.5)') xx,Potential(xx)
ENDDO

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

DO k=0,n2
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


!!! Set default variable to zero

loop = 0
r0 = r ; p0 = p; a0 = a; b0 = b
r1 = r ; p1 = p; a1 = a; b1 = b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Beginning of main iterative loop !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO
loop = loop + 1

time = DBLE(loop)*dt/2.0d0

Allocate(KE(0:n2,0:n),PE(0:n2,0:n), STAT = err)
IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF

!!! Kinetic Energy matrix

DO k=0,n2
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

DO k=0,n2
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

DO k=0,n
 DO l=0,n
  Ham(k,l) = KE(k,l) + PE(k,l)
 ENDDO
ENDDO

! Check whether Hamiltonian matrix is Hermitian
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

! Check condition number of Hamiltonian matrix
 CALL CONDCHECK(Ham)


!! Hl part of the Hamiltonian matrix (not Hermitian)

DO k=n+1,n2
 DO l=0,n
  Hl(k-n-1,l) = KE(k,l) + PE(k,l)
 ENDDO
ENDDO

DEALLOCATE(KE,PE)


!! Diagonalization of Hamiltonian matrix to find initial coefficients

IF (loop.eq.1) THEN

 DO k=0,n
  c(k) = (0.0d0,0.0d0)
!  c(k) = (1.0d0,0.0d0)
 ENDDO
 c(:) = (1.0d0,0.0d0)

 opt = .false.
! opt = .true.
 IF (opt) WRITE(*,*) 'Will find the initial coefficients'
 CALL Hamdiag(Ham,opt,nrg,c)

  CALL NORM(c)

  c0(:) = c(:); c1(:) = c(:)
ENDIF


!!! Calculation of propagation equations

CALL CALCRATE(n2,dn,c,Hl,Ham,dr,dp,da,db,dc,expnQ)

!write(*,*) r,p,a,b
!write(*,*) dr,dp,da,db
!write(*,*) dc(:)
!read(*,*)

!!! Propagation of the variables and coefficients (Leap-Frog Algorithm)

IF (loop.eq.1) THEN !	! Move half time-step (initialization)
 r1 = r0 + dr*dt/2.0d0
 p1 = p0 + dp*dt/2.0d0
 a1 = a0 + da*dt/2.0d0
 b1 = b0 + db*dt/2.0d0

 c1(:) = c0(:) + dc(:)*dt/2.0d0
! c1(:) = (0.0d0,0.0d0) !	! both methods verified to be the same !
! c1(:) = MATMUL(expnQ,c0)

 CALL NORM(c1)

 r = r1; p = p1; a = a1; b = b1; c(:) = c1(:)

ELSE

 IF (MODULO(loop,2).eq.0) THEN
  r0 = r0 + dr*dt
  p0 = p0 + dp*dt
  a0 = a0 + da*dt
  b0 = b0 + db*dt

 c0(:) = c0(:) + dc(:)*dt
! c0(:) = (0.0d0,0.0d0)
! c0(:) = MATMUL(expnQ,c0)

  CALL NORM(c0)

  r = r0; p = p0; a = a0; b = b0; c(:) = c0(:)

 ELSE

  r1 = r1 + dr*dt
  p1 = p1 + dp*dt
  a1 = a1 + da*dt
  b1 = b1 + db*dt

  c1(:) = c1(:) + dc(:)*dt
! c1(:) = (0.0d0,0.0d0)
! c1(:) = MATMUL(expnQ,c1)

  CALL NORM(c1)

  r = r1; p = p1; a = a1; b = b1; c(:) = c1(:)

 ENDIF

ENDIF


! Check normalization
 coeffsum = 0.0d0
 DO i=0,n
  coeffsum = coeffsum + CONJG(c(i))*c(i)
 ENDDO
 IF ((ABS(coeffsum-1.0d0).gt.epsilon)) WRITE(*,*) 'Coefficients not normalized'


!write(*,*) 'propgtn'
!write(*,*) 'time = ',time
!write(*,*) r,p,a,b
!write(*,*) c(:)
!read(*,*)


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

!IF (loop.lt.0) THEN
! WRITE(*,FMT='(A6,A6,A6)') ' loop ','  au  ',' kcal '
! WRITE(*,FMT='(I2,F12.7,F12.7)') loop, Eexp, Eexp*627.51
! WRITE(*,*) 'Harmonic energy: hbar*w =',SQRT(0.3684173348d0/mass)
!ENDIF


!!! Evaluating the total wavefunction over a range of positions

IF (MODULO(INT(time),10).eq.1) wrt = .true.

IF ((MODULO(INT(time),10).eq.0).AND.(wrt)) THEN
 WRITE(UNIT=21,FMT='(F12.7,F12.7,F12.7,F12.8)') time,rexp,pexp,Eexp
 WRITE(UNIT=22,FMT='(F12.7,F12.7,F12.7,F12.7,F12.7)') time,r0,p0,a0,b0

 CALL EVALWF(prog,c,Eexp)
 wrt = .false.
ENDIF


!IF (r0.lt.0.0d0) THEN
!IF ((r0.lt.0.0).OR.(a0.lt.0.0)) THEN
! WRITE(*,*) 'Program stopped at loop#: ',loop
! WRITE(*,*) 'r0 less than 0.0'
! WRITE(*,*) 'r0=',r0
! WRITE(*,*) 'r0 or a0 less than 0.0'
! WRITE(*,*) 'r0=',r0,'  a0=',a0
! CALL EVALWF(prog,c,Eexp)
! STOP
!ENDIF

IF (time.gt.runtime) EXIT


!! Reset variables to zero ::

dr = 0.0d0; dp = 0.0d0; da = 0.0d0; db = 0.0d0; dc(:) = 0.0d0
DO k=0,n2
 DO l=0,n
  IF (k.le.n) THEN
   Ham(k,l) = 0.0d0
  ELSE
   Hl(k-n-1,l) = 0.0d0
  ENDIF
 ENDDO
ENDDO


ENDDO !	! ! ! ! ! ! Ending of propagation loop


CLOSE(21); CLOSE(22); CLOSE(23); CLOSE(24)

WRITE(*,*) 'loop',loop
WRITE(*,*) 'Position',r0,'Momemtum',p0
WRITE(*,*) 'Width Parameters:',a0,b0
WRITE(*,*) 'Coefficients'
WRITE(*,*) (k,c(k),k=0,n)

WRITE(*,*) 'Program Finished'


END PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
