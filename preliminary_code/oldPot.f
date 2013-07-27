
       FUNCTION Potent(i,x0)

       USE FNVAR

       IMPLICIT NONE

       INTEGER :: i
       INTEGER, SAVE :: neg = 0 !		! default value
       DOUBLE PRECISION :: Potent,au,l,x,Ang
       DOUBLE PRECISION,intent(in) :: x0

       l = 0.0d0 !			! for angular momentum term

!  CONVERSION FACTORS
       au = 1.0d0/627.51d0 !		! kcal/mol to Hartrees
       Ang = 0.529177249d0
       x = Ang*x0 !		! au to Angstrons

     IF ((x.le.0.0d0).AND.(neg.ne.1)) THEN
      WRITE(*,*) 'Negative input into the Potential'
      neg = 1
     ENDIF


!
!ccccc INPUT THE POTENTIAL ! - serguei v2(r)
!
!!!!!
       IF (i==0) THEN
         Potent = (0.5d0*35.80125d0*(x-0.74144d0)**2)/27.212d0
       ELSE IF (i==1) THEN
         Potent = (35.80125d0*(x-0.74144d0)/27.212d0) *Ang**i
       ELSE IF (i==2) THEN
         Potent = (35.80125d0/27.212d0) * Ang**i
       ELSE
         Potent = 0.0d0
       END IF

  IF (i.eq.1000) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (i.eq.0) THEN
        Potent = (0.6706d0*x**6 +0.5367d0*x**4 +1.370d0*x*x -1.999d0)*au &
 & + l*(l+1.0d0)/(2.0d0*mass*x0**2)
       ELSE IF (i.eq.1) THEN
        Potent = (4.0236d0*x**5+2.1468d0*x**3+2.74d0*x)*au*Ang**i &
 & - l*(l+1.0d0)/(mass*x0**3)
       ELSE IF (i.eq.2) THEN
        Potent = (20.118d0*x**4 + 6.4404d0*x**2 + 2.74d0)*au*Ang**i &
 & + 3.0d0*l*(l+1.0d0)/(mass*x0**4)
       ELSE IF (i.eq.3) THEN
        Potent = (80.472d0*x**3 + 12.8808d0*x)*au*Ang**i &
 & - 12.0d0*l*(l+1.0d0)/(mass*x0**5)
       ELSE IF (i.eq.4) THEN
        Potent = (241.416d0*x**2 + 12.8808d0)*au*Ang**i &
 & + 60.0d0*l*(l+1.0d0)/(mass*x0**6)
       ELSE IF (i.eq.5) THEN
        Potent = (482.832d0*x)*au*Ang**i &
 & - 360.0d0*l*(l+1.0d0)/(mass*x0**7)
       ELSE IF (i.eq.6) THEN
        Potent = (482.832d0)*au*Ang**i &
 & + 2520.0d0*l*(l+1.0d0)/(mass*x0**8)
       ELSE
        WRITE(*,*) 'Potential expansion order too high'
	WRITE(*,*) 'i=',i
!	WRITE(*,*) 'loop=',loop
        Potent = 0.0d0
       ENDIF

!

       IF (i.eq.0) THEN
        Potent = (0.373d0*x**6 -0.140d0*x**4 -1.210d0*x*x -1.84d0)*au &
 & + l*(l+1.0d0)/(2.0d0*mass*x0**2)
!
       ELSE IF (i.eq.1) THEN
        Potent = (2.238d0*x**5 - 0.56d0*x**3 -2.42d0*x)*au*Ang**i &
 & - l*(l+1.0d0)/(mass*x0**3)
!
       ELSE IF (i.eq.2) THEN
        Potent = (11.19d0*x**4 - 1.68d0*x**2 -2.42d0)*au*Ang**i &
 & + 3.0d0*l*(l+1.0d0)/(mass*x0**4)
!
       ELSE IF (i.eq.3) THEN
        Potent = (44.76d0*x**3 - 3.36d0*x)*au*Ang**i &
 & - 12.0d0*l*(l+1.0d0)/(mass*x0**5)
!
       ELSE IF (i.eq.4) THEN
        Potent = (134.28d0*x**2 - 3.36d0)*au*Ang**i &
 & + 60.0d0*l*(l+1.0d0)/(mass*x0**6)
!
       ELSE IF (i.eq.5) THEN
        Potent = (268.56d0*x)*au*Ang**i &
 & - 360.0d0*l*(l+1.0d0)/(mass*x0**7)
!
       ELSE IF (i.eq.6) THEN
        Potent = (268.56d0)*au*Ang**i &
 & + 2520.0d0*l*(l+1.0d0)/(mass*x0**8)
!
       ELSE
!
        WRITE(*,*) 'Potential expansion order too high'
	WRITE(*,*) 'i=',i
!	WRITE(*,*) 'loop=',loop
        Potent = 0.0d0
        READ(*,*)
!
       ENDIF

  ENDIF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    END FUNCTION Potent
