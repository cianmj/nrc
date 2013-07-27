       SUBROUTINE Hamdiag(Ham,opt,nrg,c0)

	USE FNVAR

       IMPLICIT NONE

	INTERFACE
	 SUBROUTINE NORM(c0)
	  USE FNVAR
	  IMPLICIT NONE
	  DOUBLE COMPLEX, INTENT(inout) :: c0(0:n)
	 END SUBROUTINE
	END INTERFACE

       INTEGER :: i,j,num,min,n1
       DOUBLE COMPLEX, intent(in) :: Ham(0:n,0:n)
       LOGICAL :: opt
       DOUBLE PRECISION, intent(out) :: nrg(0:n)
       DOUBLE COMPLEX, intent(out), OPTIONAL :: c0(0:n)
!
       CHARACTER :: JOBVL, JOBVR
       INTEGER :: INFO, LDA, LDVL, LDVR, LWORK
       DOUBLE PRECISION :: lambda(n+1)
       DOUBLE COMPLEX :: AA(n+1,n+1), VL(n+1,n+1),VR(n+1,n+1),W(n+1)
       DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
       DOUBLE COMPLEX, ALLOCATABLE :: WORK(:)

       JOBVL= 'N'      !do not compute the left generalized eigenvectors
       JOBVR= 'V'      !compute the right generalized eigenvectors

       n1 = n+1

       ALLOCATE(RWORK(2*n1),WORK(10*n1))

       DO i=0,n
        DO j=0,n
          AA(i+1,j+1) = Ham(i,j)
        ENDDO
       ENDDO

       LDA = n1

       LDVL = n1 !		! Left eigenvectors not computed
       LDVR = n1
       LWORK = 10*n1


      CALL ZGEEV( JOBVL, JOBVR, n1, AA, LDA, W, VL, LDVL, VR, LDVR, &
     &                  WORK, LWORK, RWORK, INFO )

        IF (INFO.ne.0) THEN
          WRITE(*,*) 'PROMBLEM WITH ZGEEV'
          WRITE(*,*) 'INFO',INFO
          STOP
        ENDIF

!       DO i=1,n1
!        WRITE(*,*) 'Eigenvalue i',i,W(i)
!       ENDDO


! Picks out the eigenvector corresponding to the smallest eigenvalue !

       num = 0
       DO
10      num = num+1
        IF (num.gt.n1) EXIT
        DO j=1,n1-num
         IF (REAL(W(num)).gt.REAL(W(num+j))) goto10
        ENDDO
        EXIT
       ENDDO


       DO i=0,n
        nrg(i) = W(i+1)
       ENDDO


IF (opt) THEN

 DO i=0,n
  c0(i) = VR(i+1,num) !		! vector stored in columns..
 ENDDO

!!! Normalization of the coefficients

 CALL NORM(c0)

ENDIF

       END SUBROUTINE
