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
END INTERFACE

! Input parameters
INTEGER :: m
DOUBLE PRECISION :: v,p,b,dt
DOUBLE PRECISION :: hbar,Pi,epsilon
! 
INTEGER :: k,l,i,j
INTEGER :: s,err,err1
DOUBLE PRECISION :: acc,v1,rtemp
DOUBLE COMPLEX :: w,w1
LOGICAL :: opt
!
DOUBLE PRECISION, ALLOCATABLE :: cH(:,:),Amtrx(:,:,:),nrg(:),Potvar(:)
DOUBLE COMPLEX, ALLOCATABLE :: c(:),dc(:)
DOUBLE COMPLEX, ALLOCATABLE :: KE(:,:),PE(:,:),Ham(:,:)
DOUBLE PRECISION, ALLOCATABLE :: EPoten(:), WInt(:,:),Wnn(:,:),rv2(:)
DOUBLE PRECISION, ALLOCATABLE :: Wa(:), Wa1(:), zm(:),Em(:,:), WTm(:)
!
DOUBLE PRECISION :: coeffsum,norm,DPr,DDPr,DPtemp
DOUBLE PRECISION :: rsum,psum,psum2,rexp,pexp,Eexp
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl,dbfactrl
DOUBLE PRECISION, EXTERNAL :: DPoten,Potential
!
! Intrinsic functions
REAL, INTRINSIC :: REAL,ABS
DOUBLE PRECISION, INTRINSIC :: DBLE, DOT_PRODUCT, AIMAG
DOUBLE COMPLEX, INTRINSIC :: CONJG, CMPLX
!
! Variable passed through other subroutines
DOUBLE PRECISION, ALLOCATABLE :: rvec(:),Er(:,:),Vr(:)

!!! Universal constants
hbar = 1.0d0
Pi = 4d0*atan2(1d0,1d0)

!!! Parameter specification section

! WRITE(*,*) 'Input the number of desired basis functions'
! READ(*,*) n + 1
n = 10

mass = 2.0d0*(1822.88853d0)   ! -- used by Serguei --

r = 1.5d0 !			! initial position
v = 0.000001d0 !		! initial velocity
a = 1.0d0 !			! width parameter
b = 0.5d0 !			! phase of wave function

lim = 8 !		! expansion order for PE ( must be EVEN !)
dt = 0.1d0 !		! propagation time step
epsilon = 1E-10 !	! accuracy limit --

hm = lim + 1!		! hm must stay below 21 (-roots of Hermite max.)
m = n + 1

p = mass*v
lim1 = lim + 1

!! Allocation section

Allocate(cH(0:lim,0:lim),Amtrx(0:n,0:n,0:lim),c(0:n), Wnn(0:n,0:n), &
 & Ham(0:n,0:n),dc(0:n),nrg(0:n),EPoten(0:m),WInt(0:n,0:n), & 
 & rvec(0:lim1),Er(0:lim1,0:lim1),Vr(0:lim1),Vrsave(0:lim1),rv2(0:m), &
 & Potvar(0:lim),Wa(0:n),Wa1(0:n),zm(0:m),Em(0:m,0:m),WTm(0:m), STAT = err)
IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


!! Creating file to save data

OPEN (UNIT=21,FILE='data.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err)
OPEN (UNIT=22,FILE='plot.txt',STATUS='REPLACE',ACTION='WRITE', &
 & IOSTAT=err1)
IF ((err.ne.0).OR.(err1.ne.0)) THEN
WRITE(*,*) 'Error opening file'
 STOP
ENDIF


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

DEALLOCATE (cH)

!!! Finding roots of mth order Hermite polynomials and corrpdn. weights

CALL Hzeros(m,zm,Em)

 DO j=0,m-1
  WTm(j)= 2.0d0**(m-1) * factrl(m) * sqrt(Pi) / (DBLE(m)**2 * (Em(m-1,j))**2)
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

! c(0) = (1.0d0,0.0d0)
! DO k=1,n
!  c(k) = (0.0d0,0.0d0)
! ENDDO

!
! WRITE(*,*) 'Finding initial coefficients...'
 opt = .true.
 CALL Hamdiag(Ham,opt,nrg,c)
!
!  WRITE(*,*) (c(k),k=0,n)

ENDIF

DEALLOCATE(KE,PE)


!!! Calculation of propagation equations

If (loop.eq.1) THEN
 DPr = DPoten(1)
 DDPr = DPoten(2)
ENDIF

acc = - DPr/mass
v1 = v + acc*dt/2.0d0
r = r + v1*dt
DPr = DPoten(1)
acc = -DPr/mass
v = v1 + acc*dt/2.0d0

w1 = w - (2.0d0*(w*w)/mass + 0.5d0*DDPr)*(dt/2.0d0)
DDPr = DPoten(2)
w = w - (2.0d0*(w1*w1)/mass + 0.5d0*DDPr)*(dt/2.0d0)

a = AIMAG(w1)
b = DBLE(w1)


!!! Calculation of the effective potential
! - ? I am using a half time step of 'a' with a full step of 'r' ?

rtemp = r
DO i=0,m
 rv2(i) = sqrt(hbar/(2.0d0*a))*zm(i) + r
 r = rv2(i)
 DPtemp = DPoten(0)
 r = rtemp
 EPoten(i) = DPtemp-DPoten(0)-DPoten(1)*(sqrt(hbar/(2.0d0*a))*zm(i)) &
  & - 0.5d0*DPoten(2)*(sqrt(hbar/(2.0d0*a))*zm(i))**2
ENDDO
       
DO i=0,n
 DO j=0,n
  WInt(i,j)=WTm(0)*Em(i,0)*Em(j,0)*EPoten(0)
  DO k=1,m-1
   WInt(i,j)= WInt(i,j)+WTm(k)*Em(i,k)*Em(j,k)*EPoten(k)
  ENDDO
  Wnn(i,j)= (1.0d0/(sqrt(factrl(i)*factrl(j)*2.0d0**(i+j)*Pi)))*WInt(i,j) !&
!!     &*sqrt(hbar/(2.0*a))
! write(*,*) "Wnn",i,j,Wnn(i,j)
 ENDDO
ENDDO

!! Calculation of the summation term of the effective potential
!! over the coefficients; as seen in the coeff. propa'n. equation !

!        WRITE(*,*) "coeff:",(i,c(i),i=0,n)
         DO i=0,n
           Wa1(i)=Wnn(i,0)*c(0)
          IF (n.ne.1) THEN
           DO k=1,n-1
             Wa(i)=Wa1(i)+ Wnn(i,k)*c(k)
             Wa1(i)=Wa(i)
           ENDDO
          ELSE
           Wa(i) = Wa1(i)
          ENDIF
         ENDDO


!!! Propagation of coefficients

DO i=0,n
c(i) = c(i)+((0.0d0,-1.0d0)/hbar)*(Wa(i) + c(i)*((mass*v1**2)/2.0d0 &
 & + (2.0d0*i+1.0d0)*hbar*a/mass))*dt
ENDDO


a = AIMAG(w)
b = DBLE(w)


!!! Normalization of the coefficients

coeffsum = 0.0d0
DO i=0,n
 coeffsum = coeffsum + CONJG(c(i))*c(i)
ENDDO
!write(*,*) coeffsum
!read(*,*)

IF ((ABS(coeffsum-1.0).gt.epsilon)) THEN
! WRITE(*,*) 'Coefficients have not been normalized:',coeffsum
 norm = 1.0d0/sqrt(coeffsum)
 DO i=0,n
  c(i) = norm*c(i)
 ENDDO
ENDIF


!!! Expectation values for position, momentum and energy

rsum = 0.0d0; psum = 0.0d0; psum2 = 0.0d0
DO l=0,n-1
 rsum = rsum + sqrt(DBLE(l)+1.0d0)*DBLE(CONJG(c(l))*c(l+1))
 psum = psum + sqrt(DBLE(l)+1.0d0)*DBLE(CONJG(c(l+1))*c(l))
 psum2 = psum2 + sqrt(DBLE(l)+1.0d0)*AIMAG(CONJG(c(l))*c(l+1))
ENDDO

rexp = r + rsum/sqrt(a)
pexp = p + 2.0d0*b*psum/sqrt(a) + 2.0d0*sqrt(a)*psum2


!! Energy of system

Eexp = 0.0d0
DO i=0,n
 DO j=0,n
   Eexp = Eexp + CONJG(c(i))*c(j)*Ham(i,j)
 ENDDO
ENDDO


IF (MODULO(loop,101).eq.0) THEN

! WRITE(*,*) '***',loop,'***'
! WRITE(*,*) 'New parameters:'
! WRITE(*,*) 'Position',r,'Momemtum',p
! WRITE(*,*) 'Width Parameters:',a,b
! WRITE(*,*) (c(i),i=0,n)
! WRITE(*,*) 'Time derivatives'
! WRITE(*,*) (dc(i),i=0,n)
! READ(*,*)

 WRITE(UNIT=21,FMT='(I7,F12.7,F12.7,F12.7)') loop,rexp,pexp,Eexp*627.51
 WRITE(UNIT=22,FMT='(I7,F12.7,F12.7,F12.7,F12.7)') loop,r,p,a,b

ENDIF


IF ((r.lt.0.0).OR.(a.lt.0.0)) THEN
 WRITE(*,*) 'r or a less than 0.0'
 WRITE(*,*) 'r=',r,'  a=',a
 EXIT
ENDIF
IF (loop.gt.1E6) EXIT


!! Reset variables to zero ::

DO k=0,n
 DO l=0,n
   Ham(k,l) = 0.0d0
 ENDDO
ENDDO


ENDDO !	! ! ! ! ! ! Ending of propagation loop


CLOSE(21)
CLOSE(22)

WRITE(*,*) 'loop',loop
WRITE(*,*) 'Position',r,'Momemtum',p
WRITE(*,*) 'Width Parameters:',a,b
WRITE(*,*) 'Coefficients'
WRITE(*,*) (k,c(k),k=0,n)

WRITE(*,*) 'Program endded'


END PROGRAM
