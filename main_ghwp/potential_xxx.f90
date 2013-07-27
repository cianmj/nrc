 !
 ! 1D external potentials for ghwp evolution.
 ! Last modified: August 17, 2004
 !
 !
 function v(x)
   use ghwp1D_options
   real(rk), intent(in)  :: x
   !
   real(rk)              :: v
   real(rk)              :: r
   real(rk)              :: energyScale
   real(rk), parameter   :: mass            = 2.*1822.88853_rk
   real(rk), parameter   :: Bohr2Angstrom   = 0.529177249_rk
   real(rk), parameter   :: Hartree2KcalMol = 627.51_rk
   real(rk), parameter   :: Emass2Kg        = 9.1093897e-31_rk
   real(rk), parameter   :: AMU2Kg          = 1.6605402e-27_rk
   real(rk), parameter   :: AMU2Emass       = AMU2Kg/Emass2Kg
   !
   energyScale = Hartree2KcalMol
   !
   r = x * Bohr2Angstrom
   !
   select case (ghwp1D_potential)
     case default
       print *, ' Potential ', potential, ' is not recognized'
       !
     case ('harmonic')
       v = (0.5_rk*35.80125_rk*(r-0.74144_rk)**2)/27.212_rk
       !
     case ('quartic')
       v = 0.5_rk*sqrt(2._rk)*(x-0.74144_rk)**2 + &
               & 0.1_rk*(x-0.74144_rk)**3 + 0.1_rk*(x-0.74144_rk)**4
       !
     case ('hard') 
       !           ! -Fermi function- !
       v  = 1.0_rk / sqrt( 1.0_rk + exp(-(x-4._rk)/ghwp1D_potential_fermi_width) )
       !
     case ('eckart')
       v = 0.0273386_rk*( 2._rk/( exp(0.4*(x0-6._rk)) + exp(-0.4*(x0-6._rk)) ) )**2
       !
   end select
   !
 end function v
!
! v = 0.373_rk * r**6 - 0.14_rk * r**4 - 1.21_rk * r**2 - 1.84_rk
! v = 0.6706_rk * r**6 + 0.5367_rk * r**4 + 1.370_rk * r**2 - 1.999_rk
!
!   v = v / energyScale
