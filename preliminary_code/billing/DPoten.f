
      FUNCTION DPoten(i)

       USE FNVAR
       USE VARSAVE

       IMPLICIT NONE

       INTERFACE
        SUBROUTINE DAPPRX(Vr)
         USE FNVAR
         USE ROOT
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(out) :: Vr(0:lim1)
        END SUBROUTINE
       END INTERFACE

       INTEGER,INTENT(in) :: i
       INTEGER :: j
       DOUBLE PRECISION :: DPoten, AngPot
       DOUBLE PRECISION :: Vr(0:lim1)
!
       IF (REAL(r).ne.REAL(rsave)) THEN
        CALL DAPPRX(Vr)
        rsave = r
        DO j=0,lim1
         Vrsave(j) = Vr(j)
        ENDDO
       ELSE
         Vr(i) = Vrsave(i)
       ENDIF

       DPoten = Vr(i) + AngPot(i)

      END FUNCTION DPoten

!ccccccccccccccc

       FUNCTION factrl(a)

       IMPLICIT NONE

       INTEGER :: a,i
       DOUBLE PRECISION :: factrl

       IF (a.lt.0) THEN
        WRITE(*,*) 'Factorial of negative number cannot be taken'
        STOP
       ELSE IF (a.eq.0) THEN
        factrl = 1.0d0
       ELSE
        factrl = 1.0d0
        DO i=1,a
          factrl=factrl*DBLE(i)
        ENDDO
       ENDIF

       END FUNCTION


!ccccccccccccccc

       FUNCTION dbfactrl(a)

       IMPLICIT NONE

       INTEGER :: a,i
       DOUBLE PRECISION :: dbfactrl

       IF (a.lt.0) THEN
        WRITE(*,*) 'dbFactorial of negative number cannot be taken'
        STOP
       ELSE IF (a.eq.0) THEN
        dbfactrl = 1.0d0
       ELSE
        dbfactrl = 1.0d0
        IF (MODULO(a,2).eq.0) THEN
         DO i=2,a,2
          dbfactrl=dbfactrl*DBLE(i)
         ENDDO
        ELSE
         DO i=1,a,2
          dbfactrl=dbfactrl*DBLE(i)
         ENDDO
        ENDIF
       ENDIF

       END FUNCTION

