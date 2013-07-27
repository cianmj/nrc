PROGRAM TESTP

USE FNVAR
USE VARSAVE

IMPLICIT NONE

INTEGER :: err
! Input parameters
DOUBLE PRECISION :: hbar, Pi, epsilon, step
DOUBLE PRECISION :: p,v,b
!
DOUBLE PRECISION, ALLOCATABLE :: potvar(:)
!
! External functions
DOUBLE PRECISION, EXTERNAL :: factrl,dbfactrl
DOUBLE PRECISION, EXTERNAL :: DPoten, Potent
!
! Variables used passed through subroutines
DOUBLE PRECISION, ALLOCATABLE :: rvec(:),Er(:,:), Vr(:)

!!! Universal constants
Pi = 4d0*atan2(1d0,1d0)

!! Parameter specification section

step = 1.0E-2
p = mass*v
lim1 = lim + 1

hm = lim + 1!         ! hm must stay below 21 (-roots of Hermite max.)


!! Set variable defaults
loop = 1

!! Allocation section

Allocate(Vr(0:lim1),rvec(0:lim1),Er(0:lim1,0:lim1), &
 & Vrsave(0:lim1),Potvar(0:lim1), STAT = err)

IF (err.ne.0) THEN
 WRITE(*,*) 'Allocation error - #',err
 STOP
ENDIF


write(*,*) 'r = ',r
write(*,*) 'Order/ Approximate    /    Exact'
write(*,*) '0',DPoten(0),Potent(0,r)
write(*,*) '1',DPoten(1),Potent(1,r)
write(*,*) '2',DPoten(2),Potent(2,r)
write(*,*) '3',DPoten(3),Potent(3,r)
write(*,*) '4',DPoten(4),Potent(4,r)
write(*,*) '5',DPoten(5),Potent(5,r)
write(*,*) '6',DPoten(6),Potent(6,r)

write(*,*) 'DONE'

END PROGRAM
