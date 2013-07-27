SUBROUTINE DAPPRX(Vr)

USE FNVAR
USE ROOT

IMPLICIT NONE

INTERFACE
 SUBROUTINE Hzeros(hm,z,E)
  IMPLICIT NONE
  INTEGER, intent(in) :: hm
  DOUBLE PRECISION, intent(out) :: z(0:hm),E(0:hm,0:hm)
 END SUBROUTINE

 SUBROUTINE hermite(num,zero,Eval)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: num
  DOUBLE PRECISION, INTENT(in) :: zero(0:num)
  DOUBLE PRECISION, INTENT(out) :: Eval(0:num,0:num)
 END SUBROUTINE
END INTERFACE

INTEGER :: j,i,k
DOUBLE PRECISION :: Potvar2(0:hm-1)
DOUBLE PRECISION :: whv(0:lim1),rvec(0:lim1),Er(0:lim1,0:lim1)
DOUBLE PRECISION, intent(out) :: Vr(0:lim1)
!
!INTEGER :: lim,hm !		! variable from module
!DOUBLE PRECISION :: r,a
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl,dbfactrl
DOUBLE PRECISION, EXTERNAL :: Potential
!
!Intrinsic function
REAL, INTRINSIC :: REAL

!!-calculate once at beginning of program-!!

IF ((count.eq.0).AND.(loop.lt.2)) THEN

 ALLOCATE (z(0:hm),E(0:hm,0:hm),WT(0:hm-1))

 !!! Zeros of (hm)th order Hermite / Evaluation of all Hermite

 Call Hzeros(hm,z,E)
! WRITE(*,*) (z(i),i=0,hm-1)
 

 !!! Weighting function (for roots of (lim+1)th order Hermite)

 !write(*,*) 'weighting..'
 DO j=0,hm-1
  WT(j)= 2.0**(hm-1) * factrl(hm) * sqrt(Pi) &
 &  /((hm)**2 * (E(hm-1,j))**2)
 
  IF (E(hm-1,j).eq.0) then
   WRITE(*,*) 'Error finding weighting function'
   STOP
  ENDIF
! write(*,*) j,WT(j)
 ENDDO

 DO i=0,30
  factvar(i) = factrl(i)
 ENDDO

 count = count + 1
ENDIF

!!-calculate below after each iteration-!!

!!! Expansion coefficients of the potential (for Hermite basis)

DO j=0,hm-1
 Potvar2(j) = Potential(sqrt(hbar/(2.0d0*a))*z(j)+r)
ENDDO

DO i=0,lim1
 whv(i) = 0.0d0
 DO j=0,hm-1
   whv(i) = whv(i) + WT(j)*E(i,j)*Potvar2(j)
 ENDDO
 whv(i) = (1.0d0/(sqrt(Pi)*(2.0d0**i)*factvar(i))) * whv(i)
ENDDO

!!! Evaluate ith order Hermite at position r (first term in rvec)
!! - ? I believe that I am expanding the potential about r = 0 !!

DO i=0,lim1
 rvec(i) = 0.0d0
ENDDO

Call hermite(lim1,rvec,Er)


!!! Approximate values of the kth derivatives of the potential

DO k=0,lim1
 Vr(k) = 0.0d0
 DO i=k,lim1
   Vr(k) = Vr(k) +whv(i)*(factvar(i)/factvar(i-k))* Er(i-k,0)
 ENDDO
 Vr(k) = (sqrt(2.0d0*a/hbar)*2.0d0)**k * Vr(k)
ENDDO


END SUBROUTINE
