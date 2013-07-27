module numerov
!
!  Integration of 1-dimentional Schroedinger equation, with two-sided
!  boundary conditions, using Numerov's method.
!
   implicit none
   private
   public rk, numerovSolve
!
!
!
   logical, parameter    :: verbose   = .false.
   integer, parameter    :: rk        = selected_real_kind(14)
   real(rk)              :: eps          ! A value about halfway down the dynamic
                                         ! range of the rk-kind
   real(rk)              :: bigval       ! A value about halfway up the dynamic 
                                         ! range of rk-kind
   integer               :: npoints      ! Number of grid points
   integer               :: lstitch      ! The leftmost grid point to consider for stitching
   integer               :: rstitch      ! The rightmost grid point to consider for stitching
   integer               :: lmatch       ! The leftmost point to check for W.F. matching
   integer               :: rmatch       ! The rightmost point to check for W.F. mathing
   integer               :: lx           ! L value
   real(rk)              :: h            ! Integration step
   real(rk)              :: energy       ! Current guess for the energy. This value is set
                                         ! by every call to numerovGoal
   real(rk)              :: goal         ! Current value for the goal function. This value is
                                         ! set by every call to numerovGoal
   real(rk), allocatable :: potential(:) ! Potential on grid, including
                                         ! centrifugal barrier
   real(rk), allocatable :: psi(:,:)     ! Solution wavefunctions on grid
                                         ! psi(:,1) <- consensus wavefunction
                                         ! psi(:,2) <- left wavefunction
                                         ! psi(:,3) <- right wavefunction

contains
!
!  numerovSolve     finds a numerical solution of radial Schroedinger 
!                   equation for a given potential, angular momentum
!                   and energy range. The units must be choosen such
!                   that the prefactor in front of the kinetic energy
!                   operator is 1 (i.e. h-bar**2/(2*m) = 1).
!
   subroutine numerovSolve(lval,rmax,npt,vfun,emin,emax,eval,efun,error,nroots)
      integer, intent(in)   :: lval        ! Angular momentum
      real(rk), intent(in)  :: rmax        ! Maximum radius. Beyound this point,
                                           ! potential is assumed to be infinite
      integer, intent(in)   :: npt         ! Number of integration points
      real(rk), external    :: vfun        ! Function, returning value of the potential
                                           ! for a given coordinate. The potential is
                                           ! assumed to be regular, so that it is never
                                           ! evaluated at r=0.
      real(rk), intent(in)  :: emin        ! Lower bound of the energy range
      real(rk), intent(in)  :: emax        ! Upper bound of the energy range
      real(rk), intent(out) :: eval        ! Eigenvalue of the solution
      real(rk), intent(out) :: efun(2,npt) ! Coordinates and values of the eigenfuntion
      real(rk), intent(out) :: error       ! Residual error
      integer,  intent(out) :: nroots      ! Number of roots of the wavefuntion (not counting the
                                           ! boundary conditions)
!
!    Initializaton
!
!     eps        = tiny(1._rk)**0.25_rk
!     bigval     = huge(1._rk)**0.25_rk
      eps        = 1e-30_rk
      bigval     = 1e+30_rk
      lx         = lval
      npoints    = npt
!
      lstitch    = (1*npoints)/3
      rstitch    = (2*npoints)/3
!
!    This definition is more accurate in locating a root
!
!     lmatch     = (1+lstitch)/2
!     rmatch     = (rstitch+npoints)/2
!
!    ... but this one is much more likely to -find- a root in the first place!
!
      lmatch     = lstitch
      rmatch     = rstitch
!
      h          = rmax / (npoints-1)
      allocate (potential(npoints),psi(npoints,3))
      call numerovPotential(vfun,efun,rmax)
!
!
!    Try to find the minimum by repeated bisections. The location of the
!    minimum is returned in (energy,goal) global variables.
!
      call numerovFindMinimum(emin,emax)
!
!    Copy out the final result
!
      efun(2,1:npoints) = psi(1:npoints,1)
      eval              = energy
      error             = goal
      nroots            = countRoots(psi(1:npoints,1))
!
!    Clean up and leave
!
      deallocate (potential,psi)
   end subroutine numerovSolve
!
!  countRoots figures out how many times the wavefuntion touches or crosses zero
!
   integer function countRoots(wf)
      real(rk), intent(in)   :: wf(npoints)
      character(len=5)       :: status
      integer                :: ipt
      logical                :: zero, plus, minus

      countRoots = -1        ! Account for boundary conditions
      status     = 'ZERO'
      do ipt=1,npoints
         zero  = abs(wf(ipt))<= 2*3e-2_rk
         plus  =     wf(ipt) >    3e-2_rk
         minus =     wf(ipt) <   -3e-2_rk
         select case (status)
            case ('ZERO')
               if (.not.zero) then
                  if (plus)  status = 'PLUS'
                  if (minus) status = 'MINUS'
               end if
            case ('PLUS')
               if (.not.plus) then
                  countRoots = countRoots + 1 
                  if (minus) status = 'MINUS'
                  if (zero)  status = 'ZERO'
               end if
            case ('MINUS')
               if (.not.minus) then
                  countRoots = countRoots + 1 
                  if (plus)  status = 'PLUS' 
                  if (zero)  status = 'ZERO'
               end if
         end select
      end do

   end function countRoots
!
!  numerovFindMinimum located the minimum of the goal function by repeated
!                     bisection of the interval. The energy is optimized to
!                     machine precision.
!
   subroutine numerovFindMinimum(iel,ier)
      real(rk), intent(in) :: iel, ier
      real(rk)             :: el, em, er, gl, gm, gr

      el = iel
      er = ier
      gl = numerovGoal(el)
      gr = numerovGoal(er)

      do while(abs(er-el)>max(1e-9_rk,spacing(abs(er)+abs(el))))
         em = (el+er)/2
         gm = numerovGoal(em)
         if (gl.lt.gr) then ! Take (el,em)
            er = em
            gr = gm
         else               ! Take (em,er)
            el = em
            gl = gm
         end if
      end do

      gm = numerovGoal((el+er)/2)
   end subroutine numerovFindMinimum
!
!  numerovPotential evaluates the potential (including the centrifugal barrier)
!                   on grid points, and stores coordinates of grid points for 
!                   the benefit of the calling routine
!
   subroutine numerovPotential(vfun,efun,rmax)
      real(rk), external      :: vfun
      real(rk), intent(inout) :: efun(2,1:npoints)
      real(rk), intent(in)    :: rmax
      integer                 :: ipt

      do ipt=1,npoints
         efun(1,ipt) = h * (ipt-1) 
         if (efun(1,ipt) < 0    ) efun(1,ipt) = 0
         if (efun(1,ipt) > rmax ) efun(1,ipt) = rmax
         potential(ipt) = 0
         if (ipt>1) then ! Centrifugal barrier is trouble at r=0!
            potential(ipt) = vfun(efun(1,ipt))
            potential(ipt) = potential(ipt) + (lx*(lx+1))/efun(1,ipt)**2
         end if
      end do
      if (verbose) then
         write(6,"(' Potential (not energy-shifted) :')")
         write(6,"((2f25.10))") (efun(1,ipt), potential(ipt), ipt=1,npoints)
      end if
   end subroutine numerovPotential
!
!  numerovGoal      calculates goal function, which is minimized for given
!                   energy value
!
   real(rk) function numerovGoal(e)
      real(rk), intent(in) :: e

      energy      = e
      call numerovIntegrate
      call numerovStitch
      numerovGoal = numerovError()
      goal        = numerovGoal
      if (verbose) write(6,"('    for e = ',f15.8,' goal function is ',f15.8)") e, goal
   end function numerovGoal
!
!  numerovError     calculates the mismatch between the left and right 
!                   solutions, by integrating the square of the difference
!                   over the stitching interval
!
   real(rk) function numerovError
      numerovError = h * sum((psi(lmatch:rmatch,2)-psi(lmatch:rmatch,3))**2)
   end function numerovError
!
!  numerovStitch    integrates Schroedinger equation from both sides, then
!                   stitches the solutions at a midpoint, such that the
!                   resulting "consensus" solution is propetly normalized
!                   and continution (but not necessarily smooth).
!
   subroutine numerovStitch
      real(rk) :: leftScale, rightScale
      real(rk) :: leftMax,   rightMax
      real(rk) :: wgt, maxwgt, overlap
      integer  :: stitchPoint, ipt
!
!    Find the "consensus" solution. To do this, find a point in the middle
!    of the interval where both left and right solutions are reasonably large,
!    and scale the smaller solution until wavefunctions match at this point
!
      leftMax     = maxval(abs(psi(lstitch:rstitch,2)))
      rightMax    = maxval(abs(psi(lstitch:rstitch,3)))
      stitchPoint = lstitch
      maxwgt      = 0
      do ipt=lstitch, rstitch
         wgt = (psi(ipt,2)/leftMax)**2 + (psi(ipt,3)/rightMax)**2
         if (wgt.gt.maxwgt) then
            maxwgt      = wgt
            stitchPoint = ipt
         end if
      end do
      wgt        = maxval(abs(psi(stitchPoint,2:3)))
      leftScale  = sign(wgt/psi(stitchPoint,2),psi(stitchPoint,2))
      rightScale = sign(wgt/psi(stitchPoint,3),psi(stitchPoint,3))
!
!    Ensure normalization of the wavefunction. 
!
      overlap    =           leftScale **2 * sum(psi(1:stitchPoint,        2)**2)
      overlap    = overlap + rightScale**2 * sum(psi(stitchPoint+1:npoints,3)**2)
      overlap    = sqrt(h*overlap)
      leftScale  = leftScale  / overlap
      rightScale = rightScale / overlap
!
!    Rescale the left and right solutions, so that they match the
!    consensus solution.
!
      psi(1:npoints,2) = psi(1:npoints,2) * leftScale
      psi(1:npoints,3) = psi(1:npoints,3) * rightScale
      psi(1:stitchPoint,        1) = psi(1:stitchPoint,        2)
      psi(stitchPoint+1:npoints,1) = psi(stitchPoint+1:npoints,3)
   end subroutine numerovStitch
!
!  numerovIntegrate calculates solutions of Schroedinger equation for 
!                   the current energy guess, starting from left and
!                   right constraints
!
   subroutine numerovIntegrate
!
!    First, integrate from the left. The initial guess is determined by
!    the asymptotic behaviour of the u(r) = r**(lx+1))
!
      psi(1,2) = 0
      psi(2,2) = h**(lx+1)
      call numerovIntegrateSide(psi(1:npoints,2),3,npoints,1)
      psi(2,2) = sum(psi(1:3:2,2))/2
!
!    Now, integrate from the right. The initial guess is -almost- zero
!    Once integration is done, we set it to exactly zero, to avoid 
!    trouble in rescaling.
!
      psi(npoints,  3) = 0
      psi(npoints-1,3) = eps
      call numerovIntegrateSide(psi(1:npoints,3),npoints-2,1,-1)
      psi(npoints-1,3) = sum(psi(npoints-2:npoints:2,3))/2
   end subroutine numerovIntegrate
!
!  numerovIntegrateSide integrates Schroedinger equation, given two starting 
!                       points, either at the left, or at the start of the
!                       integration interval.
!
   subroutine numerovIntegrateSide(u,first,last,step)
      real(rk), intent(inout) :: u(npoints)  ! Wavefunction
      integer, intent(in) :: first           ! First point to generate solutions for
                                             ! wf(first-step) and wf(first-2*step)
                                             ! contain the initial guess
      integer, intent(in) :: last            ! Last point to generate solutions for
      integer, intent(in) :: step            ! Step - either +1 (left to right) 
                                             !            or -1 (right to left)
      real(rk)            :: w(npoints)      ! Temporary integration variable

      integer             :: i, m1, m2

      do i=first-2*step,first-step
        w(i) = u(i) * u2w(i)
      end do

      do i=first,last,step
        m1  = i -   step
        m2  = i - 2*step
        w(i) = w(m1)
        w(i) = w(i) + (w(m1)-w(m2))
        w(i) = w(i) + h**2 * u(m1) * (potential(m1) - energy)
        if (abs(w(i))<=eps   ) w(i) = sign(eps,   w(i))
        if (abs(w(i))>=bigval) w(i) = sign(bigval,w(i))
        u(i) = w(i)/u2w(i)
      end do

      contains 
         real(rk) function u2w(j)
           integer, intent(in) :: j
           u2w = 1._rk - h**2 * (potential(j)-energy)
         end function u2w
   end subroutine  numerovIntegrateSide

end module numerov
