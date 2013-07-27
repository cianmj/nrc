
       FUNCTION Potential(x0)

       IMPLICIT NONE

       DOUBLE PRECISION :: Potential,au,x,Ang,k
       DOUBLE PRECISION :: k
       DOUBLE PRECISION :: De,Re,ae
       DOUBLE PRECISION,intent(in) :: x0
       DOUBLE PRECISION :: hbar = 1.0d0
       DOUBLE PRECISION, INTRINSIC :: EXP, SQRT

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
!	Potential = (0.373d0*x**6 -0.140d0*x**4 -1.210d0*x*x -1.84d0)*au
!
!	Potential = (0.6706d0*x**6+0.5367d0*x**4+1.370d0*x*x-1.999d0)*au 
!
 ! quartic anharmonic oscillator
!
	Potential = 0.5d0*SQRT(2.0d0)*(x-0.74144d0)**2 + &
 		& 0.1d0*(x-0.74144d0)**3 + 0.1d0*(x-0.74144d0)**4
!
 ! quadratic potential !
!
!       Potential = (100.0d0*((x-1.0d0)**2))*au
!
!     Potential = (0.5d0*35.80125d0*(x-0.74144d0)**2)/27.212d0
!
!     Potential = (0.5d0*35.80125d0*(x-0.74144d0)**2)/27.212d0 + 0.075*x0
!
 ! morse potential !
!
!       Potential = De*(1.0d0 - EXP(-ae*(x0-Re)))**2
!
!       Potential = 0.1d0*(1.0d0 - EXP(-1.0d0*(x-0.74144d0)))**2
!
!       Potential = 0.1745d0*(1.0d0 - EXP(-1.9045d0*(x-0.74144d0)))**2
!
 ! Eckart barrier !
!
!	Potential = 0.0273386d0*(2.0d0/(EXP(0.4d0*(x0-9.0d0)) + &
!		& EXP(-0.4d0*(x0-9.0d0))))**2
!
!	Potential = 0.0273386d0*(2.0d0/(EXP(0.4d0*(x0-6.0d0)) + &
!		& EXP(-0.4d0*(x0-6.0d0))))**2
!			!! SECH(0.4d0*(x0-6.0d0)))**2 !!
!
 ! Hard Barrier !
!
!        Potential = 1.0d0 / sqrt( 1.0d0 + exp(-(x0-4.0d0)/0.005d0) )
!
!
	END FUNCTION Potential
