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
   v = 0.6706_rk * r**6 + 0.5367_rk * r**4 + 1.370_rk * r**2 - 1.999_rk
   v = v / energyScale
 end function potential
!
 program he4c60
   use settings
   implicit none
   character(len=60), parameter :: system    = "He-4 inside C60" ! What this is?
   character(len=60), parameter :: prefix    = "he4-c60"         ! Prefix for any output file names
   integer, parameter        :: npoints   = 1000              ! Number of points to use in Numerov's method
   integer, parameter        :: intervals = 1000              ! Number of intervals, used to subdivide energy range
   integer, parameter        :: minl      =    0              ! Smallest angular momentum to use
   integer, parameter        :: maxl      =   30              ! Lagest angular momentum to use
   integer, parameter        :: maxeig    =   30              ! Maximum number of solutions to look for, for a given L
   real(rk), parameter       :: rcage     =  1.5_rk           ! Cage radius, Angstrom
   real(rk), parameter       :: mass      =  4.0_rk           ! Mass of the guest, AMU
   real(rk), parameter       :: emin      = -2.0_rk           ! Minimum energy, kcal/mol
   real(rk), parameter       :: emax      = 10.0_rk           ! Maximum energy, kcal/mol

   INCLUDE 'process.f90'

 end program he4c60
