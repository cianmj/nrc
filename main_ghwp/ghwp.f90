 ! 
 !  Time propagation of a 1D Gauss-Hermite wavepacket
 !
 !  Last modified: August 17, 2004
 !
 !  Written by: S.P and C.M.-J.
 !
 !  The basis function form is:
 !
 !                      1          / 2 A \1/4   /     1/2       \
 !    \phi_n(x) = ---------------- | --- |    H | (2A)   (x-R) )|  
 !                (2**n n!)**(1/2) \ Pi  /     n\               /
 !
 !                   /           \    /              2 \
 !                exp| i P*(x-R) | exp| (-A+iB) (x-R)  |
 !                   \           /    \                /
 !
 module ghwp1D
   use ghwp1D_options
   use ghwp1D_pe
   use lapack
   use ghwp1D_properties
   implicit none
   private
   public ghwp1D_propagate, ghwp1D_propagation_order, ghwp1D_implicit_functions, &
          ghwp_coefficient_normalization
   public ghwp1D_property
   !
   type ghwp1D_property
     real(rk) :: r, p  ! Expectation values of r and p
     real(rk) :: ekin  ! Kinetic energy
     real(rk) :: epot  ! Potential energy
   end type ghwp1D_property
   !
   contains
     !
     subroutine ghwp1D_propagate(propagator,nmax,mass,v,tstep,rpab,c,prop)
       character(len=10), intent(in) :: propagator 
                                             ! Propagator to use, either
                                             ! 'optimal' or 'classical'
       integer(ik), intent(in)      :: nmax  ! Highest order of the basis functions
                                             ! in GHWP expansion. Number of basis
                                             ! functions is (ghwp_nmax+1)
       real(rk), intent(in)         :: mass  ! Particle mass, in atomic (electron) units
       real(rk), external           :: v     ! Function of one real argument, defining
                                             ! the external potential
       real(rk), intent(in)         :: tstep ! Time step, in atomic (electron) units
       type(ghwp1D_property), intent(out) :: prop  ! Properties of the wavefuntion
       !
       !  Wavefunction parameters and coefficients.
       !  The last index corresponds to the time step.
       !    1 = time (ghwp_time) (ie current time step)
       !    0 = time (ghwp_time-ghwp_tstep) (ie previous time step)
       !
       real(rk), intent(inout)    :: rpab  (4,0:1)
       complex(rk), intent(inout) :: c(0:nmax,0:1)
       !
       real(rk)    :: drpab(4)          ! Time derivatives of w.f. parameters
       complex(rk) :: dh(0:nmax,0:nmax) ! Time evolution operator of w.f. coefficients
       complex(rk) :: dc(0:nmax)        ! Time difference of the w.f. coefficients
       !
       real(rk)    :: coeff_sum         ! normalization coefficient for c(:)
       !
       !  Calculate time derivatives of the w.f. parameters, and w.f.
       !  evolution operator (which coincides with the Hamiltonian if
       !  w.f. parameters are constants)
       !
       call ghwp_parameter_derivatives(propagator,nmax,mass,v, &
                                       rpab(:,1),c(:,1),drpab(:),dh(:,:), &
                                       prop%ekin,prop%epot)
! write(*,"(4(f16.8,2x))") rpab(:,1)
! write(*,"(4(f16.8,2x)/)") drpab(:)
! write(*,"(4(f12.8,1x,f12.8,4x))") c(:,1)
       !
       !  Calculate change in the wavefunction coefficients with time.
       !  We have to calculate the change from (t-tstep) to (t+tstep),
       !  hence the time difference is 2*tstep.
       !
       call ghwp_wavefunction_increment(nmax,2*tstep,dh,c(:,0),dc)
       !
       !  Propagate coefficients and parameters in time
       !
       drpab(:) = 2*tstep*drpab(:) + rpab(:,0)
       dc   (:) = 2*tstep*dc   (:) + c   (:,0)
       !
       !  Shift wf&paramters from time slot t to time slot (t-tstep) 0
       !  and put new w.f. and parameters in time slot (t) 1
       !
       rpab(:,0) = rpab(:,1)
       c   (:,0) = c   (:,1)
       rpab(:,1) = drpab(:)
       c   (:,1) = dc   (:)
       !
       call ghwp_coefficient_normalization(nmax,c(:,1))
       !
       ! Calculation of the expectation values of the position and momentum
       !
       call ghwp_expectation_values(nmax,rpab,c(0:nmax,1),prop%r,prop%p)
       !
!write(*,*) prop%r,prop%p,prop%ekin,prop%epot
!read(*,*)
       !
     end subroutine ghwp1D_propagate
     !
     subroutine ghwp_wavefunction_increment(nmax,dt,dh,c0,dc)
       integer(ik), intent(in)  :: nmax      ! Highest order of the basis functions
       real(rk), intent(in)     :: dt        ! Change in time
       complex(rk), intent(in)  :: dh(:,:)   ! Time evolution operator 
       complex(rk), intent(in)  :: c0(:)     ! Original w.f.
       complex(rk), intent(out) :: dc(:)     ! Change in the w.f.
       !
       integer(ik) :: n
       complex(rk) :: exph(0:nmax,0:nmax)  ! current approximation to matrix 
                                           ! exponential of h, minus 1
       complex(rk) :: hn  (0:nmax,0:nmax)  ! running power of h (h**n)
       complex(rk) :: identity(0:nmax,0:nmax) ! identity matrix
       !
       exph  = 0
       hn    = dh*dt
       exponentiate: do n=1,ghwp1D_propagation_order
         exph = exph + hn
         if (n>=ghwp1D_propagation_order) exit exponentiate
         hn = (dt/n)*matmul(hn,dh)
       end do exponentiate
       !
       identity(:,:) = (0._rk,0_rk)
       do n=1,nmax
         identity(n,n) = (1._rk,0._rk)
       end do
       !
       hn = hn - identity
       dc(:) = matmul(hn(:,:),c0(:)) / dt
       !
     end subroutine ghwp_wavefunction_increment
     !
     subroutine ghwp_parameter_derivatives(propagator,nmax,mass,v,rpab,c,drpab,dh, &
                                           ekin, epot)
       character(len=*), intent(in) :: propagator
       integer(ik), intent(in)  :: nmax      ! Highest order of the basis functions
       real(rk), intent(in)     :: mass      ! Particle mass
       real(rk), external       :: v         ! External potential
       complex(rk), intent(in)  :: c(0:nmax) ! W.f. coefficients 
       real(rk), intent(in)     :: rpab(4)   ! W.f. parameters
       real(rk), intent(out)    :: drpab(4)  ! Time derivatives of the w.f. parameters
       complex(rk), intent(out) :: dh(:,:)   ! Time evolution operator for w.f. 
                                             ! coefficients.
       real(rk), intent(out)    :: ekin      ! Kinetic energy
       real(rk), intent(out)    :: epot      ! Potential energy
       !
       integer(ik)              :: ncheck    ! Highest order of the illicit
                                             ! basis functions
       integer(ik)              :: alloc     ! Allocation status
       integer(ik)              :: i
       integer(ik)              :: j
       !
       !  Matrix elements of the operators, (0:ncheck,0:nmax)
       !
       complex(rk), allocatable :: frpab(:,:,:) ! The last index correspond to
                                                ! matrix elements of the operators:
                                                ! 1: \frac{\partial}{\partial r}
                                                ! 2: \frac{\partial}{\partial p}
                                                ! 3: \frac{\partial}{\partial a}
                                                ! 4: \frac{\partial}{\partial b}
       complex(rk), allocatable :: h     (:,:)  ! Hamiltonian matrix elements
       !
       ncheck = nmax + ghwp1D_implicit_functions
       !
       allocate (frpab(0:ncheck,0:nmax,4),h(0:ncheck,0:nmax), stat=alloc)
       if (alloc/=0) then
         stop 'ghwp_parameter_derivatives - memory allocation'
       end if
       !
       call ghwp_grad_r   (nmax,ncheck,     rpab,frpab(:,:,1))
       call ghwp_grad_p   (nmax,ncheck,     rpab,frpab(:,:,2))
       call ghwp_grad_a   (nmax,ncheck,     rpab,frpab(:,:,3))
       call ghwp_grad_b   (nmax,ncheck,     rpab,frpab(:,:,4))
       !
       h = 0
       call ghwp_kinetic  (nmax,ncheck,mass,rpab,h)
       !       
       ! Calculation of the expectation value for the kinetic energy
       !
       ekin = dot_product(conjg(c),matmul(h(0:nmax,0:nmax),c))
       !
       call ghwp_potential(nmax,ncheck,v,   rpab,h)
       !
       ! Calculation of the expectation value for the potential energy
       !
       epot = dot_product(conjg(c),matmul(h(0:nmax,0:nmax),c)) - ekin
       !
       select case (propagator) 
         case default
           print *, ' Propagator ', propagator, ' is not recognized'
         case ('optimal')
           call ghwp_satisfy_constraints(nmax,ncheck,frpab,h,c,drpab)
         case ('classical')
           call ghwp_classical(rpab,mass,v,drpab)
       end select
       !
       dh(:,:) = -(0._rk,1._rk)*h(0:nmax,:)
       hamiltonian_corrections: do i=1,4
         dh(:,:) = dh(:,:) - drpab(i)*frpab(0:nmax,:,i)
       end do hamiltonian_corrections
       !
     end subroutine ghwp_parameter_derivatives
     !
     subroutine ghwp_grad_a(nmax,ncheck,rpab,da)
       integer(ik), intent(in)  :: nmax    ! Highest order of the explicit 
                                           ! basis functions
       integer(ik), intent(in)  :: ncheck  ! Highest order of the implicit 
                                           ! basis functions
       real(rk), intent(in)     :: rpab(4)
       complex(rk), intent(out) :: da(0:ncheck,0:nmax) 
                                           ! Matrix elements of the operator 
       !                                   ! \frac{\partial}{\partial a}
       integer(ik)              :: i, j
       real(rk)                 :: a
       !
       a = rpab(3)
       !
       do i=0,ncheck ! <i|
         do j=0,nmax  ! |j>
           if (i==(j-2)) then
             da(i,j) =  sqrt((j)  *real(j-1,rk))/(4._rk*a)
           else if (i==(j+2)) then
             da(i,j) = -sqrt((j+1)*real(j+2,rk))/(4._rk*a)
           else
             da(i,j) = 0._rk
           end if
         end do
       end do
       !
     end subroutine ghwp_grad_a
     !
     subroutine ghwp_grad_b(nmax,ncheck,rpab,db)
       integer(ik), intent(in)  :: nmax    ! Highest order of the explicit 
                                           ! basis functions
       integer(ik), intent(in)  :: ncheck  ! Highest order of the implicit 
                                           ! basis functions
       real(rk), intent(in)     :: rpab(4)
       complex(rk), intent(out) :: db(0:ncheck,0:nmax) 
                                           ! Matrix elements of the operator 
       !                                   ! \frac{\partial}{\partial b}
       integer(ik)              :: i, j
       real(rk)                 :: a
       !
       a = rpab(3)
       !
       do i=0,ncheck ! <i|
         do j=0,nmax  ! |j>
           if (i==j) then
             db(i,j) = cmplx(0._rk,(j+0.5_rk)/(2._rk*a),kind=rk)
           else if (i==j-2) then
             db(i,j) = cmplx(0._rk,sqrt( j   *real(j-1,rk))/(4._rk*a),kind=rk)
           else if (i==j+2) then
             db(i,j) = cmplx(0._rk,sqrt((j+1)*real(j+2,rk))/(4._rk*a),kind=rk)
           else
             db(i,j) = 0._rk
           end if
         end do
       end do
       !
     end subroutine ghwp_grad_b
     !
     subroutine ghwp_grad_r(nmax,ncheck,rpab,dr)
       integer(ik), intent(in)  :: nmax    ! Highest order of the explicit 
                                           ! basis functions
       integer(ik), intent(in)  :: ncheck  ! Highest order of the implicit 
                                           ! basis functions
       real(rk), intent(in)     :: rpab(4)
       complex(rk), intent(out) :: dr(0:ncheck,0:nmax) ! Matrix elements of the operator
       !                                               ! \frac{\partial}{\partial r}
       integer(ik)              :: i, j
       real(rk)                 :: a, b, p
       !
       a = rpab(3) ; b = rpab(4) ; p = rpab(2)
       !
       do i=0,ncheck ! <i|
         do j=0,nmax  ! |j>
           if (i==j) then
             dr(i,j) = - cmplx(0._rk,p,kind=rk)
           else if (i==j+1) then 
             dr(i,j) = -sqrt(real(j+1,rk))*cmplx(-sqrt(a),b/sqrt(a),kind=rk)
           else if (i==j-1) then
             dr(i,j) = -sqrt(real(j,  rk))*cmplx( sqrt(a),b/sqrt(a),kind=rk)
           else
             dr(i,j) = 0._rk
           end if
         end do
       end do
       !
     end subroutine ghwp_grad_r
     !
     subroutine ghwp_grad_p(nmax,ncheck,rpab,dp)
       integer(ik), intent(in)  :: nmax    ! Highest order of the explicit
                                           ! basis functions
       integer(ik), intent(in)  :: ncheck  ! Highest order of the implicit
                                           ! basis functions
       real(rk), intent(in)     :: rpab(4)
       complex(rk), intent(out) :: dp(0:ncheck,0:nmax) ! Matrix elements of the operator
       !                                               ! \frac{\partial}{\partial p}
       integer(ik)              :: i, j
       real(rk)                 :: a
       !
       a = rpab(3)
       !
       do i=0,ncheck ! <i|
         do j=0,nmax  ! |j>
           if (i==j-1) then
             dp(i,j) = cmplx(0._rk,sqrt(real(j,  rk)/a)/(2._rk),kind=rk)
           else if (i==j+1) then
             dp(i,j) = cmplx(0._rk,sqrt(real(j+1,rk)/a)/(2._rk),kind=rk)
           else
             dp(i,j) = 0._rk
           end if
         end do
       end do
       !
     end subroutine ghwp_grad_p
     !
     subroutine ghwp_kinetic(nmax,ncheck,mass,rpab,ke)
       integer(ik), intent(in)  :: nmax    ! Highest order of the explicit 
                                           ! basis functions
       integer(ik), intent(in)  :: ncheck  ! Highest order of the implicit 
                                           ! basis functions
       real(rk), intent(in)     :: mass    ! Particle mass
       real(rk), intent(in)     :: rpab(4)
       complex(rk), intent(out) :: ke(0:ncheck,0:nmax) ! Kinetic energy matrix
       !
       integer(ik)              :: i, j
       real(rk)                 :: a, b, p, scale
       complex(rk)              :: v
       !
       a = rpab(3) ; b = rpab(4) ; p = rpab(2)
       scale = - 0.5_rk / mass
       !
       do i=0,ncheck ! <i|
         do j=0,nmax  ! |j>
           if (i==j) then
             v = -p**2 -(2._rk*j+1._rk)*(b**2+a**2)/a
           else if (i==j+1) then 
             v = 2*p*sqrt(real(j+1,rk))*cmplx(-b,-a,kind=rk)/sqrt(a)
           else if (i==j-1) then
             v = 2*p*sqrt(real(j  ,rk))*cmplx(-b, a,kind=rk)/sqrt(a)
           else if (i==j+2) then
             v = sqrt((j+1)*real(j+2,rk))*cmplx((a**2-b**2)/a,-2*b,kind=rk)
           else if (i==j-2) then
             v = sqrt( j   *real(j-1,rk))*cmplx((a**2-b**2)/a, 2*b,kind=rk)
           else
             v = 0._rk
           end if
           ke(i,j) = ke(i,j) + scale * v
         end do
       end do
       !
     end subroutine ghwp_kinetic
     !
     subroutine ghwp_coefficient_normalization(nmax,c)
       integer(ik), intent(in)    :: nmax
       complex(rk), intent(inout) :: c(0:nmax)
       !
       integer(ik)                :: i
       real(rk)                   :: coeff_sum
       real(rk)                   :: norm
       !
       coeff_sum = 0._rk
       check_normalization: do i=0,nmax
         coeff_sum = coeff_sum + conjg(c(i))*c(i)
       end do check_normalization
       !
       norm = 1._rk/sqrt(coeff_sum)
       c(:) =  cmplx(norm*real(c(:),rk),norm*aimag(c(:)),kind=rk)
       !
     end subroutine ghwp_coefficient_normalization
     !
     subroutine ghwp_satisfy_constraints(nmax,ncheck,frpab,h,c,drpab)
                                           ! Highest order of the..
       integer(ik), intent(in) :: nmax     ! ..explicit basis functions
       integer(ik), intent(in) :: ncheck   ! ..implicit basis functions
       complex(rk), intent(in) :: frpab(0:ncheck,0:nmax,4)                                        ! The last index correspond to matrix elements of the operators :
                                           ! 1: \frac{\partial}{\partial r}
                                           ! 2: \frac{\partial}{\partial p}
                                           ! 3: \frac{\partial}{\partial a}
                                           ! 4: \frac{\partial}{\partial b}
       complex(rk), intent(in) :: h    (0:ncheck,0:nmax)   ! Hamiltonian matrix
       complex(rk), intent(in) :: c(0:nmax)                ! wavefunction coeff.
       real(rk), intent(out)   :: drpab(4)                 ! Time derivatives of the                                                               ! w.f. parameters 
       !
       integer(ik)             :: i
       complex(rk)             :: frpablc(0:ncheck,4)
       complex(rk)             :: hlc    (0:ncheck)
       real(rk)                :: weight
       !
       ! Variables passed through dgelss subroutine (linear equation solver)
       !
       real(rk)                :: A(1:2*(ncheck+1),4)
       real(rk)                :: B(1:2*(ncheck+1),1)
       !
       do i=1,4
         frpablc(:,i) = matmul(frpab(:,:,i),c(:))
       end do
       !
       hlc(:) = -(0._rk,1._rk)*matmul(h(:,:),c(:))
       !
       !
       !  DGELSS computes the minimum norm solution to a real linear least
       !  squares problem:
       !                  Minimize 2-norm(| B - A*x |).
       !  using the singular value decomposition (SVD) of A. A is an M-by-N
       !  matrix which may be rank-deficient.
       !
       ! Construction matrices A & B:
       !
       separate_reim: do i=0,ncheck
         weight       = ghwp_weighting_function(i,nmax)
         A(2*i+1,1:4) = weight*real (frpablc(i,1:4),rk)
         A(2*i+2,1:4) = weight*aimag(frpablc(i,1:4))
         B(2*i+1,1)   = weight*real (hlc(i),rk)
         B(2*i+2,1)   = weight*aimag(hlc(i))
       end do separate_reim
       !
       call lapack_gelss(A,B)
       !
       ! The solution vector, x, is returned in the B array.
       !
       drpab(1:4) = B(1:4,1)
       !
     end subroutine ghwp_satisfy_constraints
     !
     function ghwp_weighting_function(n,nmax) result(w)
       integer(ik), intent(in) :: n    ! Basis function order
       integer(ik), intent(in) :: nmax ! Max. explicit order
       real(rk)                :: w    ! The weight; DGELSS will square it later
       !
       real(rk) :: arg
       !
       select case (ghwp1D_weight)
         case default
           stop 'ghwp_eighting_function - unimplemented'
         case ('fermi')
           !
           ! Fermi function :
           !
           arg = (n-nmax - ghwp1D_implicit_fermi_shift)/ghwp1D_implicit_fermi_width ;
           w   = 1.0_rk / sqrt( 1.0_rk + exp(-arg) )
         case ('step')
           !
           ! Heaviside step function :
           !
           w = 1._rk
           if (n<=nmax) w = 0._rk
         end select
       !
       !
     end function ghwp_weighting_function
     ! 
     subroutine ghwp_classical(rpab,mass,v,drpab)
       !
       real(rk), intent(in)  :: rpab(4)
       real(rk), intent(in)  :: mass
       real(rk), external    :: v
       real(rk), intent(out) :: drpab(4)
       !
       real(rk)              :: r, p, a, b
       complex(rk)           :: w, dw
       real(rk)              :: dv
       real(rk)              :: ddv
       !
       r = rpab(1); p = rpab(2); a = rpab(3); b = rpab(4)
       w = cmplx(b,a,kind=rk)
       !
       call potential_derivatives(r,v,dv,ddv)
       !
 !   write(*,*) r,v(r),dv,ddv
 !   read(*,*)
       !
       drpab(1) = p / mass
       drpab(2) = - dv
       dw = - 2._rk*w*w/mass - 0.5_rk*ddv
       drpab(3) = aimag(dw)
       drpab(4) = real(dw,rk)
       !
     end subroutine ghwp_classical
     !
     subroutine potential_derivatives(r,v,dv,ddv)
       real(rk), intent(in)  :: r
       real(rk), external    :: v
       real(rk), intent(out) :: dv
       real(rk), intent(out) :: ddv
       !
       real(rk), parameter   :: limit = 1E-5
       !
       dv = ( v(r+limit)-v(r-limit)) / (2._rk*limit)
       !
       ddv = ( v(r+limit)+v(r-limit)-2._rk*v(r)) / (limit*limit)
       !
     end subroutine potential_derivatives
     !
 end module ghwp1D
