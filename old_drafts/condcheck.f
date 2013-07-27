SUBROUTINE CONDCHECK(Ham)

USE FNVAR

IMPLICIT NONE

INTERFACE
 SUBROUTINE ZGESVD(JOBU,JOBVT,M,NN,HH,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,RWORK,INFO)
  CHARACTER :: JOBU, JOBVT
  INTEGER  :: INFO, LDA, LDU, LDVT, LWORK, M, NN
  DOUBLE PRECISION :: RWORK(5*NN), S(NN)
  DOUBLE COMPLEX :: HH(LDA,NN),U(LDU,NN),VT(LDVT,NN),WORK(LWORK)
 END SUBROUTINE
END INTERFACE

INTEGER :: k,l
INTEGER, SAVE :: save
DOUBLE COMPLEX, INTENT(in) :: Ham(0:n,0:n)
DOUBLE PRECISION :: cutoff
!
!! Variables passed through ZGESVD
CHARACTER ::JOBU, JOBVT
INTEGER :: INFO, LDA, LDU, LDVT, LWORK, M, NN
DOUBLE PRECISION :: RWORK(1:5*(n+1)),S(1:n+1)
DOUBLE COMPLEX :: HH(1:n+1,1:n+1),U(1:n+1,1:n+1), VT(1:n+1,1:n+1)
DOUBLE COMPLEX, ALLOCATABLE :: WORK(:)

cutoff = 500

IF (loop.lt.2) THEN
 WRITE(*,*) 'The condition number of the Hamiltonian will be calculated'
 LWORK = 5000
 ALLOCATE (WORK(LWORK))
ELSE
 LWORK = SAVE
 ALLOCATE (WORK(LWORK))
ENDIF

JOBU = 'N'
JOBVT = 'N'

M = n + 1 !		! # rows of matrix A
NN = M !		! # columns   "   "
LDA = M
LDU = M
LDVT = M

!!! Construct matrix

DO k = 0,n
 DO l = 0,n
  HH(k+1,l+1) = Ham(k,l)
 ENDDO
ENDDO

CALL ZGESVD(JOBU,JOBVT,M,NN,HH,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,RWORK,INFO)

IF (INFO.ne.0) THEN
 WRITE(*,*) 'Error with the subroutine ZGESVD'
 WRITE(*,*) 'INFO',INFO
ENDIF

!WRITE(*,*) 'Condition number is given by difference between &
! & the largest and smallest singular values '
IF ((MODULO(loop,101).eq.0).AND.(ABS(S(NN)-S(1)).gt.cutoff)) THEN
 WRITE(UNIT=21,FMT='(A,F20.18)') 'Condition number : ',ABS(S(NN)-S(1))
ENDIF

SAVE = WORK(1)
DEALLOCATE(WORK)

RETURN

END SUBROUTINE
