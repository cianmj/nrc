MODULE FNVAR
 IMPLICIT NONE
 SAVE
 INTEGER :: loop = 0
 INTEGER, PARAMETER :: n = 4		! # of basis fns -1
 INTEGER, PARAMETER :: lim = 8		! order of potential exp.
 DOUBLE PRECISION :: dt = 0.02d0, runtime = 4000
 DOUBLE PRECISION, PARAMETER :: epsilon = 1E-12
 DOUBLE PRECISION, PARAMETER :: hbar = 1.0d0
 DOUBLE PRECISION, PARAMETER :: mass = 2.0d0*(1822.88853d0)
 DOUBLE PRECISION :: r = 1.5d0 ! 1.8d0 ! 1.398397231814d0 !
 DOUBLE PRECISION :: v = 0.0005d0
 DOUBLE PRECISION :: a = 25.0d0 ! 2.086759d0 ! 28.0d0
 DOUBLE PRECISION :: b = 0.0d0
 INTEGER :: lim1,hm
 DOUBLE PRECISION :: p, Pi
END MODULE FNVAR

MODULE VARSAVE
 IMPLICIT NONE
 SAVE
 DOUBLE PRECISION :: rsave = 1232.0d0 ! arbitrary !
 DOUBLE PRECISION, ALLOCATABLE :: Vrsave(:)
END MODULE VARSAVE

MODULE ROOT
 IMPLICIT NONE
 SAVE
 INTEGER :: count = 0
 DOUBLE PRECISION :: factvar(0:30)
 DOUBLE PRECISION, ALLOCATABLE :: z(:),E(:,:),WT(:)
END MODULE ROOT
