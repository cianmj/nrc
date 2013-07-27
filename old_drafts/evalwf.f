SUBROUTINE EVALWF(prog,c,Eexp)

USE FNVAR

IMPLICIT NONE

INTEGER, INTENT(in) :: prog
DOUBLE PRECISION, INTENT(in) :: Eexp
DOUBLE COMPLEX, INTENT(in) :: c(0:n)
!
INTEGER :: i,k,int
DOUBLE PRECISION :: x, rx(0:n), Hx(0:n,0:n), density, time, Ek
DOUBLE COMPLEX :: psi(0:n), wfn
!
DOUBLE PRECISION, EXTERNAL :: factrl,Potential
DOUBLE COMPLEX, EXTERNAL :: EXPN
DOUBLE COMPLEX, INTRINSIC :: EXP, CONJG


WRITE(UNIT=prog,FMT='(A,I8)') '#', loop

time = (dt)*DBLE(loop)
Ek = 0.5d0*mass*v**2

x = - 4.001d0

DO i = 1,10000
!DO i = 1,5000
 x = x + 0.002d0

 rx(0) = SQRT(2.0d0*a)*(x-r)
 IF (n.ge.1) THEN
  DO k=1,n
   rx(k) = 0.0d0
  ENDDO
 ENDIF
 CALL hermite(n,rx,Hx)

 wfn = (0.0d0,0.0d0)
 DO k = 0,n
  psi(k) = 1.0d0/SQRT((2.0d0**k)*factrl(k))*(2.0d0*a/Pi)**(0.25d0)* &
	& Hx(k,0) * EXPN(x)
  wfn = wfn + c(k)*psi(k)
 ENDDO

density = wfn * CONJG(wfn)

! - morse -
!WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density*0.02d0 + Ek
! - quartic anharmonic -
!WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density + Eexp
! - harmonic -
!WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density*0.25d0 + Eexp
! - serguei -
!WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density*0.00075d0 + Eexp
! - erkart -
!WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density*0.02d0 + Ek
! - hard barrier !
WRITE(UNIT=prog,FMT='(F12.5,5X,F15.8)') x, density*0.3 + Eexp

ENDDO

WRITE(UNIT=prog,FMT='(A1)')
WRITE(UNIT=prog,FMT='(A1)')

!write(*,*) 'Finished creating wfn. density file'


END SUBROUTINE



FUNCTION EXPN(x)
 USE FNVAR
 IMPLICIT NONE

 DOUBLE PRECISION :: x
 DOUBLE COMPLEX :: EXPN
!
 DOUBLE COMPLEX, INTRINSIC :: EXP

 EXPN= EXP((0.0,1.0)*p*(x-r)/hbar)* EXP((-a + (0.0,1.0)*b)*((x-r)**2))

END FUNCTION EXPN
