 module ghwp1D_options
   integer, parameter :: ik = kind(1)
   integer, parameter :: rk = kind(1d0)
   !
   !  Expansion order, used in the evaluation of matrix exponential
   !  for time evolution of the w.f. coefficients. Order 1 corresponds
   !  to simple leap-frog approach.
   !
   integer(ik) :: ghwp1D_propagation_order           = 20
   !
   !  Number of basis functions which are being implicitly suppressed
   !  during propagation.
   !
   integer(ik) :: ghwp1D_implicit_functions          = 2
   !
   !  Steepness of the Fermi weighting function for the contributions
   !  from the implicitly propagated basis functions. 
   !
   !  ghwp1D_weight can be either 'fermi' or 'step'
   !
   character(len=10) :: ghwp1D_weight                = 'fermi'
   real(rk)          :: ghwp1D_implicit_fermi_shift  = 0.0_rk
   real(rk)          :: ghwp1D_implicit_fermi_width  = 0.05_rk
   !
   !  Expansion order for the external potential. Must be even to
   !  prevent infinite spreading of the wavepacket.
   !
   integer(ik) :: ghwp1D_potential_order             = 8
   ! 
   !  Selection of the potential for the 1D system.
   !  ghwp1D_potential can be : 'harmonic', 'quartic', 'hard' or 'eckart'
   !    For the hard barrier potential.. we use the Fermi function with 
   !    width given by 'ghwp1D_potential_fermi_width'.
   !
   character(len=8)  :: ghwp1D_potential             = 'quartic'
   real(rk)          :: ghwp1D_potential_fermi_width = 0.005_rk
   !
   !
 end module ghwp1D_options
