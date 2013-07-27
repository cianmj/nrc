SUBROUTINE NORM(c)

USE FNVAR

IMPLICIT NONE

INTEGER :: i
DOUBLE COMPLEX, INTENT(inout) :: c(0:n)
DOUBLE PRECISION :: coeffsum, nrm
DOUBLE PRECISION, INTRINSIC :: SQRT, ABS
DOUBLE COMPLEX, INTRINSIC :: DCMPLX

coeffsum = 0.0d0
DO i=0,n
 coeffsum = coeffsum + CONJG(c(i))*c(i)
ENDDO

IF ((ABS(coeffsum-1.0d0).gt.epsilon)) THEN
! WRITE(*,*) 'Coefficients are being normalized: ',coeffsum
 nrm = 1.0d0/SQRT(coeffsum)
 DO i=0,n
  c(i) = DCMPLX(nrm*DBLE(c(i)),nrm*DIMAG(c(i)))
 ENDDO
ENDIF

END SUBROUTINE

