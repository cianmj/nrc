 !
 !  Matrix elements of the potential over Gauss-Hermite basis
 !
 !  Last modified: July 27, 2004
 !
 !  Written by: S.P and C.M.-J.
 !
 ! Will solve the following potential integral using Gauss-Hermite integration:
 !                                                                       2
 !  /         *                  /      1       \1/2 /                 -z
 !  | \phi_m(x) v(x) \phi_n(x) = | ------------ |    | H(z) v(x) H(z) e  dx
 !  /                            |  m+n         |    /  m         n
 !             1/2               \ 2    m! n! Pi/
 !    z = (2 A)   (x-R)                                   - integral -
 !        
 !               max                /    -1/2      \
 !    integral ~ Sum  W(z ) H (z ) v|(2 A)   z  + R| H (z )
 !               i=0     i   m  i   \         i    /  n  i
 !
 ! Find the roots 'z(:)' of 'integration_order' order Hermite polynomial
 ! and evaluate all H. polynomials up to that order at those roots. Then 
 ! using the weights 'W' at those roots approximate the integral.
 !
 !
 module ghwp1D_pe
   use ghwp1D_options
   use ghwp_math
   implicit none
   !
   contains
     !
     subroutine ghwp_potential(nmax,ncheck,v,rpab,pe)
       integer(ik), intent(in)  :: nmax      ! Highest order of the 
                                             ! basis functions
       integer(ik), intent(in)  :: ncheck    ! Highest order of the 
                                             ! illicit basis functions
       real(rk), external       :: v                   ! External potential
       real(rk), intent(in)     :: rpab(4)             ! W.f. parameters
       complex(rk), intent(out) :: pe(0:ncheck,0:nmax) ! Potential energy matrix
       !
       integer(ik)              :: i, j, k, l          ! summation indices
       integer(ik)              :: integration_order
       real(rk)                 :: r                   ! Centre of the basis set
       real(rk)                 :: a                   ! Width parameter
       real(rk)                 :: factrl  (0:ncheck)
       real(rk)                 :: dbfactrl(0:ncheck)
       real(rk),allocatable     :: z      (:)          ! roots of H. polyn.
       real(rk),allocatable     :: wght   (:)          ! weights of roots
       real(rk),allocatable     :: hermite(:,:)        ! H.p. evaluated at roots
       real(rk),allocatable     :: v_val  (:)
       real(rk)                 :: vij
       !
       integration_order = nmax + ncheck + ghwp1D_potential_order + 1
       !
       allocate(z(0:integration_order),wght(0:integration_order), &
               hermite(0:integration_order,0:integration_order), &
               v_val(0:integration_order))
       !
       r = rpab(1) ; a = rpab(3)
       !
       call factorials(ncheck,factrl,dbfactrl)
       !
       call hermite_polynomial(integration_order,z,hermite,wght)
       !
       vals: do i=0,integration_order-1
         v_val(i) = wght(i)*v(sqrt(0.5_rk/a)*z(i)+r)
       end do vals
       !
       do i=0,ncheck
         do j=0,nmax
           vij = 0
           do k=0,integration_order-1
             vij = vij + hermite(i,k)*hermite(j,k)*v_val(k)
           end do
           pe(i,j) = pe(i,j) + vij*sqrt( 1._rk/(Pi*2._rk**(i+j)*factrl(i)*factrl(j)) )
         end do
       end do
       !
       deallocate(z, wght, hermite, v_val)
       !
     end subroutine ghwp_potential
     !
 end module ghwp1D_pe
