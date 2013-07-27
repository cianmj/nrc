 !
 !  Matrix elements of the potential over Gauss-Hermite basis
 !
 !  Last modified: July 27, 2004
 !
 !  Written by: S.P and C.M.-J.
 !
 module ghwp_math
   use ghwp1D_options
   implicit none
   private
   public factorials, hermite_polynomial, evaluate_hermites
   public Pi
   !
   integer(ik), parameter  :: max_hermite_order = 100 
   logical                 :: h_good_order(0:max_hermite_order) = .false.
   real(rk)                :: h_roots     (0:max_hermite_order, &
                                           0:max_hermite_order)
   real(rk)                :: h_weights   (0:max_hermite_order, &
                                           0:max_hermite_order)
   real(rk)                :: h_values    (0:max_hermite_order, &
                                           0:max_hermite_order, &
                                           0:max_hermite_order)
   real(rk)                :: Pi
   !
   contains
     !
     subroutine factorials(num,factrl,dbfactrl)
       integer(ik), intent(in) :: num
       real(rk), intent(out)   :: factrl(0:num)    ! factorial & 
       real(rk), intent(out)   :: dbfactrl(0:num)  ! db. factorial 
                                                   ! values stored
       !
       integer(ik)             :: i, j
       !
       factrl(0) = 1._rk
       do j=1,num
         factrl(j) = factrl(j-1)*j
       end do
       !
       dbfactrl(0) = 1._rk
       do j=2,num
         dbfactrl(1) = 1._rk
         dbfactrl(j) = dbfactrl(j-2)*j
       end do
       !
     end subroutine factorials
     !
     ! The algorithm for calculating the roots of the hermite polynomials
     ! was taken from 'Calculation of Gauss Quadrature Rules*', by Gene H. Golub
     ! and John H. Welsh 
     !
     subroutine hermite_polynomial(n,z,E,weights)
       integer, intent(in)   :: n
       real(rk), intent(out) :: z(0:n)       ! roots of Hermite poly. of degree n
       real(rk), intent(out) :: E(0:n,0:n)   ! Hermite polynomial upto degree n-1
                                             ! evaluated at 'z'
       real(rk), intent(out) :: weights(0:n) ! weights associated with 'z'
       !
       integer             :: i, j
       integer             :: info
       real(rk)            :: a(1000)
       real(rk)            :: b(1000)
       !
       real(rk)            :: factrln
       !
       if ( h_good_order(n) ) then
         z      (0:n)     = h_roots  (0:n,n)
         weights(0:n)     = h_weights(0:n,n)
         E      (0:n,0:n) = h_values (0:n,0:n,n)
         return
       end if
       !
       Pi = 4._rk*atan2(1._rk,1._rk)
       !
       a(1:n) = 0._rk
       do i=1,n-1
         b(i) = sqrt(0.5_rk*i)
       end do
       !
       call dsterf(n,a,b,info)
       if (info/=0) stop 'dsterf failed'
       !
       z(0:n-1) = a(1:n)
       !
       call evaluate_hermites(n,z,E)
       !
       factrln = 1._rk
       do i=1,n
        factrln = factrln * real(i,rk)
       end do
       !
       do j=0,n-1
         weights(j)= 2._rk**(n-1)*factrln*sqrt(Pi)/((real(n,rk)*E(n-1,j))**2)
         if (E(n-1,j)==0) then
           write(*,*) 'Error finding weighting function'
           stop
         end if
       end do
       !
       h_good_order(n) = .true.
       h_roots  (0:n,n)     = z      (0:n)
       h_weights(0:n,n)     = weights(0:n)
       h_values (0:n,0:n,n) = E      (0:n,0:n)
       !
     end subroutine hermite_polynomial
     !
     subroutine evaluate_hermites(n,x,E)
       integer(ik), intent(in) :: n
       real(rk), intent(in)    :: x(0:n)
       real(rk), intent(out)   :: E(0:n,0:n)
       !
       integer(ik)             :: i, k, l
       real(rk)                :: H(0:n,0:n)
       !
       H(0,0)=1._rk
       if (n>=1) then
         H(1,1)=2._rk
         H(1,0)=0._rk
       end if
       !
       do k=0,n
         do l=k+1,n
           H(k,l)=0.0_rk
         end do
       end do
       !
       do k=2,n
         H(k,0)= -2._rk*real(k-1,rk)*H(k-2,0)
       end do
       !
       do k=2,n
         do l=1,k
           H(k,l)=2._rk*H(k-1,l-1) - 2._rk*real(k-1,rk)*H(k-2,l)
         end do
       end do
       !
       E(0,0:n) = 1
       !
       if (n>=1) then
         E(1,0:n) = 2*x(0:n)
       end if
       !
       do l=0,n
         do k=2,n
           E(k,l) = H(k,k)
           do i=1,k           ! evaluate Hermite, order k, at x(l)
             E(k,l) = x(l)*E(k,l) + H(k,k-i)
           end do
         end do
       end do
       !
     end subroutine evaluate_hermites
     !
 end module ghwp_math
