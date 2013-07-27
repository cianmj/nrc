
       FUNCTION AngPot(i)

       USE FNVAR

       IMPLICIT NONE

       INTEGER,INTENT(in) :: i
       INTEGER, SAVE :: neg = 0 !		! default value
       DOUBLE PRECISION :: AngPot,l,x0
       DOUBLE PRECISION :: mass

       l = 0.0d0 !			! for angular momentum term
       x0 = r

     IF ((r.le.0.0d0).AND.(neg.ne.1)) THEN
      WRITE(*,*) '*'
      neg = 1
     ENDIF

       IF (i.eq.0) THEN
        AngPot = l*(l+1.0d0)/(2.0d0*mass*x0**2)
       ELSE IF (i.eq.1) THEN
        AngPot = - l*(l+1.0d0)/(mass*x0**3)
       ELSE IF (i.eq.2) THEN
        AngPot = 3.0d0*l*(l+1.0d0)/(mass*x0**4)
       ELSE IF (i.eq.3) THEN
        AngPot = - 12.0d0*l*(l+1.0d0)/(mass*x0**5)
       ELSE IF (i.eq.4) THEN
        AngPot = 60.0d0*l*(l+1.0d0)/(mass*x0**6)
       ELSE IF (i.eq.5) THEN
        AngPot = - 360.0d0*l*(l+1.0d0)/(mass*x0**7)
       ELSE IF (i.eq.6) THEN
        AngPot = + 2520.0d0*l*(l+1.0d0)/(mass*x0**8)
       ELSE IF (i.eq.7) THEN
        AngPot = + 20160.0d0*l*(l+1.0d0)/(mass*x0**9)
       ELSE IF (i.eq.8) THEN
        AngPot = + 181440.0d0*l*(l+1.0d0)/(mass*x0**10)
       ELSE IF (i.eq.9) THEN
        AngPot = + 1814400.0d0*l*(l+1.0d0)/(mass*x0**11)
       ELSE IF (i.eq.10) THEN
        AngPot = + 19958400d0*l*(l+1.0d0)/(mass*x0**12)
       ELSE IF (i.eq.11) THEN
        AngPot = + 239500800d0*l*(l+1.0d0)/(mass*x0**13)
       ELSE IF (i.eq.12) THEN
        AngPot = + 3113510400d0*l*(l+1.0d0)/(mass*x0**14)
       ELSE IF (i.eq.13) THEN
        AngPot = + 43589145600d0*l*(l+1.0d0)/(mass*x0**15)
       ELSE IF (i.eq.14) THEN
        AngPot = + 653837184000d0*l*(l+1.0d0)/(mass*x0**16)
       ELSE
        WRITE(*,*) 'Potential expansion order too high'
        WRITE(*,*) 'Program stopped in AngPot function'
	STOP
       ENDIF

       END FUNCTION AngPot
