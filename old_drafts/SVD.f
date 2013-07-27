SUBROUTINE SVD(loop,dn,dn1,DrLC,DpLC,DaLC,DbLC,HlC,dr,dp,da,db)

IMPLICIT NONE

INTEGER :: i,k,l,num
INTEGER, SAVE :: save
INTEGER, INTENT(in) :: loop, dn,dn1
DOUBLE COMPLEX, INTENT(in), DIMENSION(0:dn) :: DrLC,DpLC,DaLC,DbLC,HlC
DOUBLE PRECISION, INTENT(out) :: dr,dp,da,db
DOUBLE PRECISION :: V(4,4),Sinv(4,4),cutoff, Utrans(4,1:2*dn+2)
DOUBLE PRECISION :: UBB(4,1),VSinv(4,4),XX(4,1)
!
DOUBLE PRECISION, INTRINSIC :: DBLE,AIMAG,MATMUL
!
!! Variables passed through DGESVD
CHARACTER ::JOBU, JOBVT
INTEGER :: INFO, LDA, LDU, LDVT, LWORK, M, N
DOUBLE PRECISION :: AA(1:2*dn+2,1:4), BB(1:2*dn+2,1)
DOUBLE PRECISION :: S(4), U(2*dn+2,4), VT(4,4)
DOUBLE PRECISION, ALLOCATABLE :: WORK(:)

cutoff = 1E-12

IF (loop.eq.1) THEN
 LWORK = 5000
 ALLOCATE (WORK(5000))
ELSE
 LWORK = SAVE
 ALLOCATE (WORK(SAVE))
ENDIF

JOBU = 'A'
JOBVT = 'A'

M = 2*dn+2 !		! # rows of matrix A
N = 4 !			! # columns   "   "
LDA = M
LDU = M
LDVT = N


!!! Construct matrices (break into real and imaginary parts)

num = 0
DO k=0,dn
 num = num + 1
 AA(k+num,1:4)=(/DBLE(DrLC(k)),DBLE(DpLC(k)), &
  & DBLE(DaLC(k)),DBLE(DbLC(k))/)
 AA(k+num+1,1:4)=(/AIMAG(DrLC(k)),AIMAG(DpLC(k)), &
  & AIMAG(DaLC(k)),AIMAG(DbLC(k))/)

 BB(k+num,1) = DBLE(HlC(k))
 BB(k+num+1,1) = AIMAG(HlC(k))
ENDDO

CALL DGESVD(JOBU,JOBVT,M,N,AA,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)

IF (INFO.ne.0) THEN
 WRITE(*,*) 'Error with the subroutine DGESVD'
 WRITE(*,*) 'INFO',INFO
ENDIF


!!! Calculation of the solutions to Ax = B

V(:,:) = TRANSPOSE(VT(:,:))

DO i = 1,4
 DO k = 1,4
  IF (i.eq.k) THEN
   IF (S(i).lt.cutoff) THEN
    Sinv(i,k) = 0
!    WRITE(*,*) 'Clipping occurred..', i, S(i)
!    WRITE(UNIT=21,FMT='(A)') 'Clipping occurred'
!    WRITE(UNIT=21,FMT='(F12.6)') S(i) == 0
   ELSE    
    Sinv(i,k) = 1 / S(i)
   ENDIF
  ELSE
   Sinv(i,k) = 0
  ENDIF
 ENDDO
ENDDO

Utrans(:,:) = TRANSPOSE(U(:,:))

UBB(:,:) = MATMUL(Utrans,BB)
VSinv(:,:) = MATMUL(V,Sinv)
XX(:,:) = MATMUL(Vsinv,UBB)

dr = XX(1,1)
dp = XX(2,1)
da = XX(3,1)
db = XX(4,1)

SAVE = WORK(1)
DEALLOCATE(WORK)

RETURN

END SUBROUTINE
