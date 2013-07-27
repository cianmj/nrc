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
   v = v / energyScale
 end function potential
!
 program d2S
   use settings
   implicit none
   character(len=60), parameter :: system    = "HD inside S cage" ! What this is?
   character(len=60), parameter :: prefix    = "hd-S"         ! Prefix for any output file names
   integer, parameter        :: npoints   = 1000              ! Number of points to use in Numerov's method
   integer, parameter        :: intervals =10000              ! Number of intervals, used to subdivide energy range
   integer, parameter        :: minl      =   33              ! Smallest angular momentum to use
   integer, parameter        :: maxl      =   34              ! Lagest angular momentum to use
   integer, parameter        :: maxeig    =   30              ! Maximum number of solutions to look for, for a given L
   real(rk), parameter       :: rcage     =  2.0_rk           ! Cage radius, Angstrom
   real(rk), parameter       :: mass      =  3.0_rk           ! Mass of the guest, AMU
   real(rk), parameter       :: emin      =  7.5_rk           ! Minimum energy, kcal/mol
   real(rk), parameter       :: emax      =  9.0_rk           ! Maximum energy, kcal/mol

   INCLUDE 'process.f90'

 end program d2S
