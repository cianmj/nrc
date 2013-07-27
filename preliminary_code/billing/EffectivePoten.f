       SUBROUTINE EffectivePoten(Wnn,n)

	USE FNVAR
	USE ROOT

       IMPLICIT NONE

       INTEGER :: i,j,k,factrl,m
       INTEGER, intent(in) :: n
       DOUBLE PRECISION :: Pi, hbar=1.0d0
       DOUBLE PRECISION :: WT(0:100)
       DOUBLE PRECISION :: WInt(0:100,0:100)
       DOUBLE PRECISION, intent(out) :: Wnn(0:100,0:100)
       DOUBLE PRECISION :: EPoten(0:100)
       DOUBLE COMPLEX :: WInt1(0:100,0:100)
!
	DOUBLE PRECISION, EXTERNAL :: DPoten

	m = n + 1
       Pi=4d0*atan2(1d0,1d0)

ccccc Integration of the effective potential w, with the
ccccc  ith and jth order Hermite basis set, over all x
ccccc The number of Hermite function used is n-1, evaluated at
ccccc the grid points of an mth order H. polynomial !

         DO j=0,m
           WT(j)= 2.0**(m-1) * factrl(m) * sqrt(Pi)
     &                /(m**2 * (E(m-1,j))**2)

           IF (E(m-1,j).eq.0) then
             WRITE(*,*) 'Error finding weighting function'
             STOP
           ENDIF
         ENDDO

        DO i=0,m
          EPoten(i) = Poten(sqrt(hbar/(2.0*a))*z(i)+x) - Poten(0) &
     &              - DPoten(1)*(sqrt(hbar/(2.0*a))*z(i)) & 
     &              - 0.5d0*DPoten(2)*(sqrt(hbar/(2.0*a))*z(i))**2
        ENDDO

       DO i=0,n
         DO j=0,n
           WInt1(i,j)=WT(0)*E(i,0)*E(j,0)*EPoten(0)
           DO k=1,m-1
             WInt(i,j)= WInt1(i,j)+WT(k)*E(i,k)*E(j,k)*EPoten(k)
             WInt1(i,j)=WInt(i,j)
           ENDDO
            Wnn(i,j)= (1.0/(sqrt(factrl(i)*factrl(j)*2.0**(i+j)*Pi)))*
c!!!!     &              sqrt(hbar/(2.0*a))*
     &   WInt(i,j)
c           write(*,*) "Wnn",i,j,Wnn(i,j)
         ENDDO
       ENDDO

ccccc The Eff. Potential integral has been found, this can now be
ccccc used to find the initial coefficient or be used in the
ccccc coefficient equation, in the Verlet algorithm !


       END SUBROUTINE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       FUNCTION factrl(a)

       IMPLICIT NONE

       INTEGER :: a,i,factrl

       factrl = 1
       DO i=1,a
         factrl=factrl*i
       ENDDO

       END FUNCTION
   
