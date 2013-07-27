 module settings
   use numerov
   implicit none
   public
   real(rk) :: lengthScale ! multiplying length in internal units by this
                           ! should give Angstroms
   real(rk) :: energyScale ! multiplying energy in internal units by this
                           ! should give kcal/mol
 end module settings
!
 function potential(x) result(v)
   use settings
   implicit none

   real(rk), intent(in)  :: x
   real(rk)              :: v
   real(rk)              :: r

   r = x * lengthScale
   v = 0.372596_rk * r**6 - 0.13986_rk * r**4 - 1.21002_rk * r**2 - 1.84_rk
!   v = 0.6706_rk * r**6 + 0.5367_rk * r**4 + 1.370_rk * r**2 - 1.999_rk
!
! ! Quadratic Potentian
!    v = 0.25_rk * r**2
!
! ! Morse Potential
!   v = 0.1756d0*(1.0d0 - EXP(-0.40551d0*(r-1.41d0)))**2
!!

   v = v / energyScale
 end function potential
!
 program h2S
   use settings
   implicit none
   character(len=60), parameter :: system    = "H2 inside S cage" ! What this is?
   character(len=60), parameter :: prefix    = "h2-S"         ! Prefix for any output file names
   integer, parameter        :: npoints   = 1000              ! Number of points to use in Numerov's method
   integer, parameter        :: intervals = 1000              ! Number of intervals, used to subdivide energy range
   integer, parameter        :: minl      =    0              ! Smallest angular momentum to use
   integer, parameter        :: maxl      =    4              ! Lagest angular momentum to use
   integer, parameter        :: maxeig    =   10              ! Maximum number of solutions to look for, for a given L
   real(rk), parameter       :: rcage     =  2.0_rk           ! Cage radius, Angstrom
   real(rk), parameter       :: mass      =  2.0_rk           ! Mass of the guest, AMU
   real(rk), parameter       :: emin      = -20.0_rk           ! Minimum energy, kcal/mol
   real(rk), parameter       :: emax      = 20.0_rk           ! Maximum energy, kcal/mol

   INCLUDE 'process.f90'

 end program h2S
