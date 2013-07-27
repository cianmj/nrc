!
!  Common processing part of the radial Schroedinger equation code.
!  All problem-specific values go to the main program, this is
!  completely problem-independent!
!
   real(rk), parameter :: Bohr2Angstrom   = 0.529177249_rk 
   real(rk), parameter :: Hartree2KcalMol = 627.51_rk
   real(rk), parameter :: Emass2Kg        = 9.1093897e-31_rk
   real(rk), parameter :: AMU2Kg          = 1.6605402e-27_rk
   real(rk), parameter :: AMU2Emass       = AMU2Kg/Emass2Kg
   integer, parameter  :: ofile           = 10
   logical, parameter  :: verbose         = .false.

   type solution
      logical  :: valid
      integer  :: nroots
      real(rk) :: score
      real(rk) :: energy
      real(rk) :: vector(2,npoints)
   end type solution

   type scan
      logical  :: valid
      integer  :: nroots
      real(rk) :: el, er
      real(rk) :: score
      real(rk) :: energy
   end type scan

   type(solution)      :: roots (0:maxeig, 0:maxl)         ! Solutions
   type(scan)          :: scores(intervals,0:maxl)         ! Results of the solutions scan
   integer             :: ne              (0:maxl)

   real(rk), external  :: potential
   real(rk)            :: estep, el, er
   integer             :: lval, interval, nx, nr
   type(solution)      :: result

   lengthScale = Bohr2Angstrom
   energyScale = Hartree2KcalMol / (2 * mass * AMU2Emass )
   estep       = (emax-emin)/intervals
   call header
!
   write(6,"('Running for ',a)") trim(system)
   roots (:,:)%valid = .false.
   scores(:,:)%valid = .false.
   do lval=minl,maxl
      do interval=1,intervals
         er = emin + estep*interval
         el = er - estep
         scores(interval,lval)%el = el
         scores(interval,lval)%er = er
!
!       numerovSolve works in special units, where hbar**2/(2*m) = 1
!       we'll need to make sure the units we use are consistent!
!
         el = el/energyScale
         er = er/energyScale
         call numerovSolve(lval,rcage/lengthScale,npoints,potential,el,er, &
                           result%energy,result%vector,result%score,result%nroots)
         result%energy      = energyScale * result%energy
         result%vector(1,:) = lengthScale * result%vector(1,:)
!
!       Warning: We don't fix up normalization of the wavefunction here,
!                so it's going to be off by sqrt(lengthScale).
!
!       Store results for the scan output
!
         scores(interval,lval)%valid  = .true.
         scores(interval,lval)%nroots = result%nroots
         scores(interval,lval)%score  = result%score
!
!       Is this the best solution this L and number of roots?
!
!        if (result%score>0.1_rk) cycle
         nr = result%nroots
         if (nr<0 .or. nr>maxeig) cycle
         if (roots(nr,lval)%valid .and. roots(nr,lval)%score<result%score ) cycle
         if (.not.roots(nr,lval)%valid) ne(lval) = ne(lval) + 1
         roots(nr,lval) = result
         roots(nr,lval)%valid = .true.
!
      end do ! interval
      write(6,"(' For L = ',i3,' found ',i4,' roots')") lval, ne(lval)
   end do ! lval
   write(6,"(' Finished, writing out the results')")
!
!  Done scanning for solutions
!
   call print_energies
   call print_functions
   call print_scores

   contains

     subroutine openfile(suffix)
        character(len=*)   :: suffix
        character(len=150) :: buffer
        integer            :: in,out

        write(buffer,"(a,'-',a,'.out')") trim(prefix), trim(suffix)
        out = 1
        do in=1,len(buffer)
           if (buffer(in:in)/=' ') then
              buffer(out:out) = buffer(in:in)
              out             = out + 1
           end if
        end do
        buffer(out:) = ' '
        open(ofile,file=trim(buffer),form='formatted',status='unknown')
     end subroutine openfile

     subroutine header
        call openfile('header')
        write(ofile,"(' name                            = ',a)") trim(system)
        write(ofile,"(' grid points                     = ',i)") npoints
        write(ofile,"(' energy bins                     = ',i)") intervals
        write(ofile,"(' minimum L                       = ',i)") minl
        write(ofile,"(' maximum L                       = ',i)") maxl
        write(ofile,"(' maximum number of eigenvalues/L = ',i)") maxeig
        write(ofile,"(' cage radius (Angstrom)          = ',f12.5)") rcage
        write(ofile,"(' particle mass (AMU)             = ',f12.5)") mass
        write(ofile,"(' minimum energy (kcal/mol)       = ',f12.5)") emin
        write(ofile,"(' maximum energy (kcal/mol)       = ',f12.5)") emax
        write(ofile,"(' energy scale (kcal/mol)         = ',f12.5)") energyScale
        write(ofile,"(' length scale (Angstrom)         = ',f12.5)") lengthScale
        close(ofile)
     end subroutine header
!
     subroutine print_energies
        integer :: lval, ie
        logical :: prev_missing

        call openfile('eigenvalues')
        write(ofile,"('# ',a)") trim(system)
        write(ofile,"(a2,1x,a4,1x,a4,1x,a4,2(1x,a25))") &
              '# ',' L ','Root','Nzero',' Energy (kcal/mol)', 'Score'
        do lval=minl,maxl
           prev_missing = .false.
           do ie=0,maxeig
              if (.not.roots(ie,lval)%valid) then
                 prev_missing = .true.
                 cycle
              end if
              if (prev_missing) then
                 prev_missing = .false.
                 write(6,"(' One or more levels missing for L = ',i4,' before level ',i4)") lval, ie
              end if
              write(ofile,"(2x,3(1x,i4),1x,f25.15,5x,g11.5)") &
                    lval, ie, roots(ie,lval)%nroots, roots(ie,lval)%energy, roots(ie,lval)%score
           end do ! ie
        end do ! lval
        close(ofile)
     end subroutine print_energies

     subroutine print_functions
        character(len=120) :: buffer
        integer            :: lval, ie, ipt

        do lval=minl,maxl
           do ie=0,maxeig
              if (.not.roots(ie,lval)%valid) cycle
              write(buffer,"('func-',i4,'-',i4)") lval, ie
              call openfile(buffer)
              write(ofile,"('# ',a)") trim(system)
              write(ofile,"('# L = ',i4,' Root = ',i4,' Zeros = ',i4,' Eigenvalue = ',f25.15)") &
                    lval, ie, roots(ie,lval)%nroots, roots(ie,lval)%energy
              write(ofile,"('# Multiply psi by ',f12.6,' to get proper normalization')") &
                    sqrt(1._rk/lengthScale)
              write(ofile,"('#',a10,1x,a15)") 'R, Angstrom', 'Psi'
              do ipt=1,npoints
                 write(ofile,"(1x,f10.5,1x,f15.5)") roots(ie,lval)%vector(:,ipt)
              end do
              close(ofile)
           end do ! ie
        end do ! lval

     end subroutine print_functions

     subroutine print_scores
        character(len=120) :: buffer
        integer            :: lval, interval
        real(rk)           :: er,el

        do lval=minl,maxl
           write(buffer,"('score-',i4)") lval
           call openfile(buffer)
           write(ofile,"('# ',a)") trim(system)
           write(ofile,"('# L = ',i4)") lval
           write(ofile,"('#',a15,1x,a15,1x,a15,1x,a15,1x,a6)") &
                 'Emin(kcal)', 'Emax(kcal)', 'Score', 'OptE', 'Nzeros'
           do interval=1,intervals
              el = scores(interval,lval)%el
              er = scores(interval,lval)%er
              write(ofile,"(1x,f15.10,1x,f15.10,1x,g15.5,1x,f15.9,1x,i6)")       &
                    scores(interval,lval)%el,    scores(interval,lval)%er,     &
                    scores(interval,lval)%score, scores(interval,lval)%energy, &
                    scores(interval,lval)%nroots
           end do
           close(ofile)
        end do ! lval

     end subroutine print_scores
