SUBROUTINE SINGLE(dr,dp,da,db)

USE FNVAR

IMPLICIT NONE

DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
!
DOUBLE PRECISION, EXTERNAL :: DPoten

dr = p / mass
dp = - DPoten(1) - DPoten(3)/(8.0d0*a) - DPoten(5)/(128.0d0*a**2)
da = -4.0d0*a*b/mass
db = 2.0d0*(a**2-b**2)/mass - 0.5d0*DPoten(2) - DPoten(4)/(16.0d0*a)

RETURN

END SUBROUTINE
