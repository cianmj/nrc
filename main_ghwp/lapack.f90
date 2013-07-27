module lapack
!
!  Simplistic type-agnostic LAPACK interface
!
!  use accuracy
  implicit none
  !
  interface lapack_gelss
    module procedure lapack_cgelss
    module procedure lapack_zgelss
    module procedure lapack_sgelss
    module procedure lapack_dgelss
  end interface ! lapack_gelss
  !
  contains
  !
  subroutine lapack_cgelss(a,b)
    complex, intent(inout) :: a(:,:)
    complex, intent(inout) :: b(:,:)
    !
    external cgelss
    real                   :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex                :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real                   :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call cgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), rwork, info)
    !
    if (info/=0) then
      write (*,"(' cgelss returned ',i)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_cgelss
  !
  subroutine lapack_zgelss(a,b)
    double complex, intent(inout) :: a(:,:)
    double complex, intent(inout) :: b(:,:)
    !
    external zgelss
    double precision       :: s    (   min(size(a,dim=1),size(a,dim=2)))
    double complex         :: work (50*max(size(a,dim=1),size(a,dim=2)))
    double precision       :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), rwork, info)
    !
    if (info/=0) then
      write (*,"(' cgelss returned ',i)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_zgelss
  !
  subroutine lapack_sgelss(a,b)
    real, intent(inout) :: a(:,:)
    real, intent(inout) :: b(:,:)
    !
    external sgelss
    real                :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real                :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call sgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), info)
    !
    if (info/=0) then
      write (*,"(' sgelss returned ',i)") info
      stop 'lapack_sgelss - sgelss failed'
    end if
  end subroutine lapack_sgelss
  !
  subroutine lapack_dgelss(a,b)
    double precision, intent(inout) :: a(:,:)
    double precision, intent(inout) :: b(:,:)
    !
    external dgelss
    double precision    :: s    (   min(size(a,dim=1),size(a,dim=2)))
    double precision    :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), info)
    !
    if (info/=0) then
      write (*,"(' dgelss returned ',i)") info
      stop 'lapack_dgelss - dgelss failed'
    end if
  end subroutine lapack_dgelss
  !
end module lapack
