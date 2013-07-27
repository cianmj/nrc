
       FUNCTION Potential(x0)

       IMPLICIT NONE

       DOUBLE PRECISION :: Potential,au,x,Ang,k
       DOUBLE PRECISION :: k
       DOUBLE PRECISION :: De,Re,ae
       DOUBLE PRECISION,intent(in) :: x0
       DOUBLE PRECISION :: hbar = 1.0d0
       DOUBLE PRECISION, INTRINSIC :: EXP

       De = 0.1756d0 !          ! Everything is in Atomic Units (au)
       Re = 1.41d0 !    ! H2 -Re = 1.41d0  ... it was set to 2.41
       ae = 0.40551d0

!  CONVERSION FACTORS
       au = 1.0d0/627.51d0 !            ! kcal/mol to Hartrees
       Ang = 0.529177249d0
       x = Ang * x0 !           ! au to Angstrons

!
!ccccc INPUT THE POTENTIAL
!
 ! - serguei v(r) !
!
	Potential = (0.373d0*x**6 -0.140d0*x**4 -1.210d0*x*x -1.84d0)*au
!
!	Potential = (0.6706d0*x**6+0.5367d0*x**4+1.370d0*x*x-1.999d0)*au 
!
 ! quadratic potential !
!
!       Potential = (0.25d0 * x**2)*au
!
 ! morse potential !
!
!       Potential = De*(1.0d0 - EXP(-ae*(x-Re)))**2
!
!        Potential = 0.1756d0*(1.0d0 - EXP(-0.40551d0*(x-1.41d0)))**2
!
	END FUNCTION Potential
