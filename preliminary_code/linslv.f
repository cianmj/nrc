SUBROUTINE LINSLV(loop,dn,dn1,DrLC,DpLC,DaLC,DbLC,HlC,dr,dp,da,db)

IMPLICIT NONE

INTEGER :: k,l,num
INTEGER, SAVE :: SIZE
INTEGER, INTENT(in) :: loop, dn,dn1
DOUBLE COMPLEX, INTENT(in), DIMENSION(0:dn) :: DrLC,DpLC,DaLC,DbLC,HlC
DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
!
DOUBLE PRECISION, INTRINSIC :: DBLE,AIMAG


!! Variables passed through dgelss
INTEGER :: INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
DOUBLE PRECISION :: RCOND,pct,cdnum2
DOUBLE PRECISION :: AA(1:2*dn1+2,1:4), BB(1:2*dn1+2,1), S(4)
DOUBLE PRECISION, ALLOCATABLE :: WORK(:)


pct = 0.85d0 !	-was set to 0.15 ! weighting for each preceding row

IF (loop.eq.1) THEN
 SIZE = 1000
 LWORK = 1000
 ALLOCATE (WORK(SIZE))
ELSE
 LWORK = SIZE
 ALLOCATE (WORK(SIZE))
ENDIF

M = 2*dn1+2 !		! # rows of matrix A
N = 4 !			! # columns   "   "
NRHS = 1 !		! # columns of the matrices B and X
RCOND = 1E-10 !		! treat values less than this as zeros

LDA = M
LDB = M


!!! Construct matrices (break into real and imaginary parts)

num = 0
DO k=0,dn
 num = num + 1
 AA(k+num,1:4)=(pct**(k+num-1))*(/DBLE(DrLC(k)),DBLE(DpLC(k)), &
  & DBLE(DaLC(k)),DBLE(DbLC(k))/)
 AA(k+num+1,1:4)=(pct**(k+num))*(/AIMAG(DrLC(k)),AIMAG(DpLC(k)), &
  & AIMAG(DaLC(k)),AIMAG(DbLC(k))/)

 BB(k+num,1) = (pct**(k+num-1))*DBLE(HlC(k))
 BB(k+num+1,1) = (pct**(k+num))*AIMAG(HlC(k))

ENDDO

!write(*,*) 'linslv'
!DO k=1,2*dn1
! DO l=1,4
!  WRITE(*,*) k,l,AA(k,l)
! ENDDO
! WRITE(*,*) 'b',k,BB(k,1)
!ENDDO
!read(*,*)

CALL DGELSS(M,N,NRHS,AA,LDA,BB,LDB,S,RCOND,RANK,WORK,LWORK,INFO)

IF (INFO.ne.0) THEN
 WRITE(*,*) 'Error solving the equations- (dgelss)'
 WRITE(*,*) 'INFO',INFO
ENDIF


!IF (RANK.ne.(2*dn1)) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WRITE(*,*) 'Rank of AA: ', RANK
! CONTINUE
!ENDIF


cdnum2 = S(1)/S(4)

IF ((cdnum2.gt.200).AND.(MODULO(loop,101).eq.0)) THEN
 WRITE(UNIT=21,FMT='(A)') 'Condition number2:'
 WRITE(UNIT=21,FMT='(I6,A2,F12.4)') loop,'  ',cdnum2
ENDIF

IF (cdnum2.gt.2000) THEN
 WRITE(*,*) 'Program stopped at loop #', loop
 WRITE(*,*) 'Condition number too larger ~ ', cdnum2
 READ(*,*)
ENDIF


dr = BB(1,1)
dp = BB(2,1)
da = BB(3,1)
db = BB(4,1)


SIZE = WORK(1)

END SUBROUTINE
