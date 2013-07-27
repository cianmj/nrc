SUBROUTINE CALCQUAD(dr,dp,da,db)

USE FNVAR

IMPLICIT NONE

DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
!
DOUBLE PRECISION, EXTERNAL :: DPoten

dr = p / mass
dp = - DPoten(1)
da = - 4.0d0*a*b/mass
db = (2.0d0*(a**2-b**2)/mass) - 0.5d0*DPoten(2)

END SUBROUTINE

