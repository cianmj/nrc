
       !Calculation of an nth order Hermite Polynomial
       !at a point z(i)

      Subroutine hermite(n,z,E)
        Implicit none
        Integer k,l,i
        Integer, intent(in) :: n
        Double Precision :: H(0:n,0:n),Eval(0:n,0:n)
        Double Precision, intent(in) :: z(0:n)
        Double Precision, intent(out) :: E(0:n,0:n)

      !Construction of the coeffient matrix H

        H(0,0)=1.0d0
        H(1,1)=2.0d0
        H(1,0)=0.0d0

        do k=0,n
          do l=k+1,n
            H(k,l)=0.0d0
          enddo
        enddo

        do k=2,n
          H(k,0)= -2.0d0*DBLE(k-1)*H(k-2,0)
        enddo

        do k=2,n
          do l=1,k
            H(k,l)=2.0d0*H(k-1,l-1) - 2.0d0*DBLE(k-1)*H(k-2,l)
          enddo
        enddo

        do l=0,n
          E(0,l)=1.0d0
        enddo

        do l=0,n
            E(1,l)=2.0d0*z(l)
        enddo

        do l=0,n
            do k=2,n 
              Eval(k,l) = H(k,k)
              do i=1,k           ! evaluate Hermite, order k, at z(l)
                 E(k,l) = z(l)*Eval(k,l) + H(k,k-i)
                 Eval(k,l)=E(k,l)
              enddo
            enddo
        enddo

!       do k=0,n
!        do l=0,n-1
!         write(*,*) k,l,E(k,l)
!        enddo
!       enddo

      End Subroutine

