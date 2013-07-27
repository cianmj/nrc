 !
 !  Expectation values for GHWP wavefunction, and calculation of w.f. density.
 !
 !  Last modified: August 17, 2004
 !
 !  Written by: S.P and C.M.-J.
 !
 module ghwp1D_properties
   use ghwp1D_options
   use ghwp_math
   implicit none
   !
   private
   public ghwp_expectation_values, wavefunction_density
   !
   contains
     !
     subroutine ghwp_expectation_values(nmax,rpab,c,rexp,pexp)
       integer(ik), intent(in) :: nmax
       real(rk), intent(in)    :: rpab(4)    ! wavefunction parameters(r,p,a,b)
       complex(rk), intent(in) :: c(0:nmax)  ! basis function coefficients
       real(rk), intent(out)   :: rexp       ! expectation value of position (r)
       real(rk), intent(out)   :: pexp       !      "        "   of momentum (p)
       !
       integer(ik)             :: i
       real(rk)                :: r, p ,a, b
       real(rk)                :: rsum, psum1, psum2
       !
       r = rpab(1)
       p = rpab(2)
       a = rpab(3)
       b = rpab(4)
       !
       rsum = 0._rk; psum1 = 0._rk; psum2 = 0._rk
       do i=0,nmax-1
         rsum  = rsum  + sqrt(real(i,rk)+1)*real (conjg(c(i)  )*c(i+1),rk)
         psum1 = psum1 + sqrt(real(i,rk)+1)*real (conjg(c(i+1))*c(i  ),rk)
         psum2 = psum2 + sqrt(real(i,rk)+1)*aimag(conjg(c(i)  )*c(i+1))
       end do
       !
       rexp = r + rsum/sqrt(a)
       pexp = p + 2._rk*b*psum1/sqrt(a) + 2._rk*sqrt(a)*psum2
       !
     end subroutine ghwp_expectation_values
     !
     subroutine wavefunction_density(x,nmax,rpab,c,density)
       real(rk)                :: x
       integer(ik), intent(in) :: nmax
       real(rk), intent(in)    :: rpab(4)
       complex(rk), intent(in) :: c   (0:nmax)
       real(rk), intent(out)   :: density
       !
       integer(ik)             :: k
       real(rk)                :: r, p, a, b
       real(rk)                :: factrl  (0:nmax)
       real(rk)                :: dbfactrl(0:nmax)
       real(rk)                :: rx  (0:nmax)
       real(rk)                :: Hrx (0:nmax,0:nmax)
       complex(rk)             :: psi (0:nmax)
       complex(rk)             :: wfn
       !
       r = rpab(1); p = rpab(2); a = rpab(3); b = rpab(4)
       !
       call factorials(nmax,factrl,dbfactrl)
       !
       rx(:) = sqrt(2.*a)*(x-r)
       call evaluate_hermites(nmax,rx,Hrx)
       !
       wfn = (0._rk,0._rk)
       do k=0,nmax
         psi(k) = 1._rk / sqrt( (2._rk**k)*factrl(k) ) * (2.*a/Pi)**(0.25_rk) &
                 * Hrx(k,0) * exponential(x)
         wfn = wfn + c(k)*psi(k)
       end do
       !
       density = wfn*conjg(wfn)
       !
       contains
         !
         function exponential(x) result(e)
           real(rk), intent(in) :: x
           complex(rk)          :: e
           !
           e = exp(cmplx(0.,p*(x-r),kind=rk))*exp((cmplx(-a,b,kind=rk))*(x-r)**2)
           !
         end function exponential
         !
     end subroutine wavefunction_density
     !
 end module ghwp1D_properties
