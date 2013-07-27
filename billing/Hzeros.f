       Subroutine Hzeros(n,z,E)

       Implicit none

       Interface
        SUBROUTINE hermite(n,z,E)
         IMPLICIT NONE
         INTEGER, intent(in) :: n
         DOUBLE PRECISION, intent(in) :: z(0:n)
         DOUBLE PRECISION, intent(out) :: E(0:n,0:n)
        END SUBROUTINE
       End Interface

       Integer k,i,j,m,l
       Integer, intent(in) :: n
       Double Precision :: PI,z1(0:n)
       Double Precision, intent(out) :: z(0:n),E(0:n,0:n)
       Double Precision :: H(0:n,0:n)
       Double Precision :: Eval(0:n,0:n),znew(0:n)

       PI=4d0*atan2(1d0,1d0)

!       write(*,*) 'Order of Hermite polynomial= ',n

       if (n.gt.21) then
         write(*,*) 'Order of polynomial too high'
         stop
       endif

!!! Approximate the location of the roots

      IF (Modulo(n,2) .eq. 0) THEN
        do i=0,MIN(n,4)
          if (modulo(i,2) .ne. 0) then
           z(i)= (i*PI)/(2.0d0*sqrt(2.0d0*DBLE(n)+1.0d0))
          else
           z(i) = 0.0d0
          endif
        enddo
        
       if (n.gt.4) then
        do j=5,MIN(n-1,10)
          if (modulo(i,2) .ne. 0) then
            z(j)= z(j-2) + 1.1d0*(z(j-2)-z(j-4))
          else
           z(i) = 0.0d0
          endif
        enddo
       endif

      ELSE
        do i=0,MIN(n,5)
          if (modulo(i,2) .eq. 0) then
            z(i)= (DBLE(i)*PI)/(2.0d0*sqrt(2.0d0*DBLE(n)+1.0d0))
          else
           z(i) = 0.0d0
          endif
        enddo

       if (n.gt.5) then
        do j=6,MIN(n-1,10)
          if (modulo(i,2) .eq. 0) then
            z(j)= z(j-2) + 1.1d0*(z(j-2)-z(j-4))
          else
           z(i) = 0.0d0
          endif
        enddo
       endif

      ENDIF

!        write(*,*) 'Array z'
!        write(*,*) (z(i),i=0,n-1)


      Call hermite(n,z,E)


! Newton-Raphson iteration for each root z(i) to obtain more
! accurate values... works between n=1-14 ????

      if (Modulo(n,2) .eq. 0) then
        if (n.eq.0) then
          write(*,*) "There are no roots"
        else
           i=1
10         do while(.true.)
              znew(i)=z(i)-(E(n,i)/(2.0d0*DBLE(n)*E(n-1,i)))

              if (abs(znew(i)-z(i)).lt. 1E-12) then
!                write(*,*) i,z(i)
               i=i+2
                if (i.gt.n) goto20
                if (i.gt.10) then
                 z(i-1) = 0.0d0
                 z(i)= z(i-2) + (sqrt(6.0d0/(DBLE(n)+1.0d0)) + &
 & sqrt(6.0d0/DBLE(n))/2.0d0)
                 Call hermite(n,z,E)
                endif
               goto10
              
              else
                z(i)=znew(i)
                 Call hermite(n,z,E)

              endif
            enddo
        endif
      else
           i=0
15          do while(.true.)
              znew(i)=z(i)-((E(n,i)/(2.0*n*E(n-1,i))))

               if (abs(znew(i)-z(i)).lt. 1E-12) then
!                write(*,*) z(i)
                i=i+2
                 if (i.gt.n) goto 20
                 if (i.gt.10) then
                  z(i-1) = 0.0d0
                  z(i)= z(i-2) + (sqrt(6.0d0/(DBLE(n)+1.0d0)) + &
 & sqrt(6.0d0/DBLE(n))/2.0d0)
                  Call hermite(n,z,E)
                 endif
                goto 15

               else

                 z(i)=znew(i)
                 Call hermite(n,z,E)

               endif
            enddo
       endif


!!! End of iteration (organize roots in increamenting order)

20    if (Modulo(n,2) .eq. 0) then
       do i=1,n-1,2
        j=(i-1)/2
        z1(j)=z(i)
       enddo
      else
       do i=0,n-1,2
        j=i/2
        z1(j)=z(i)
       enddo
      endif

      do j=INT((n+1)/2),n
       z1(j)=0
      enddo

      do i=0,INT((n-1)/2)
        z(i) = - z1(INT((n-1)/2)-i)
      enddo

      IF (modulo(n,2).eq.0) THEN
        do i=INT((n+1)/2),n-1
         z(i)= z1(-INT((n+1)/2)+i)
        enddo
      ELSE
        do i=INT((n+1)/2),n-1
         z(i)= z1(-INT((n-1)/2)+i)
        enddo
      ENDIF
        
!       write(*,*) "Correct roots:"
!       write(*,*) (z(k),k=0,n-1)


      CALL hermite(n,z,E)

!       do k=0,n
!        do l=0,n-1
!         write(*,*) k,l,E(k,l)
!        enddo
!       enddo

      End Subroutine
