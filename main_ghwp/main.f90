 program main
   use ghwp1D_options
   use ghwp1D
   use ghwp1D_properties 
   implicit none
   character(len=10)        :: optimal
   character(len=10)        :: classical
   integer(ik)              :: nmax = 4                 ! == basis size - 1
   real(rk)                 :: mass = 2.*1822.88853_rk  ! particle mass
   real(rk), external       :: v                        ! external potential
   real(rk)                 :: tstep = 0.02_rk          ! propagation time step
   real(rk)                 :: rpab (4,0:1)             ! wavepacket parameters
   real(rk)                 :: rpab2 (4,0:1)
   complex(rk), allocatable :: c (:,:)                  ! wavepacket coefficients
   complex(rk), allocatable :: c2 (:,:)
   type(ghwp1D_property)    :: prop                     ! w.p. properties
   type(ghwp1D_property)    :: prop2
   integer(ik)              :: loop = 0
   real(rk)                 :: x
   real(rk)                 :: energy, energy2
   real(rk)                 :: density, density2
   !
   allocate( c(0:nmax,0:1),c2(0:nmax,0:1) )
   !
   open (10,file='optimal.10'  ,status='replace',action='write')
   open (11,file='classical.11',status='replace',action='write')
   open (12,file='potential.12',status='replace',action='write')
   open (13,file='density.13',  status='replace',action='write')
   !
   write(*,*) 'Beginning propagation of wavepacket in ',ghwp1D_potential,'potential   &
            &  expanded to',int(ghwp1D_potential_order),'th order.','  Using ',nmax+1,&
            &  ' basis functions',' and a propagation time step of:',real(tstep),'fs',& 
            &  '  The',ghwp1D_weight,' weighting function is used for the propagation & 
            &  of the',ghwp1D_implicit_functions,' implicit basis functions.'
   !   
   ! Initial starting parameter:
   !
   optimal   = 'optimal'        ! use our optimum propagation method
   rpab(1,:) = 1.5_rk           ! center position of wavepacket(r)
   rpab(2,:) = 0.0005_rk*mass   ! momentum        "     "      (p)
   rpab(3,:) = 25._rk           ! wavepacket width (a)
   rpab(4,:) = 0._rk            !     "      phase (b)
   !
   c(:,0) = cmplx(0.,0.,kind=rk)
   c(:,0) = cmplx(1.,0.,kind=rk)
   !
   call ghwp_coefficient_normalization(nmax,c(:,0))
   c(:,1) = c(:,0)
   !
   ! implement classical propagation of 'r' and 'p'
   classical  = 'classical'
   rpab2(:,:) = rpab(:,:)
   c2(:,:)    = c(:,:)
   !
   !
   x = -5.01_rk
   do
     x = x + 0.01_rk
     write(12,fmt='(F6.2,F14.8)') x, v(x)
     if (x>10._rk) exit
   end do 
   !
   do while(.true.)
     loop = loop + 1
     call ghwp1D_propagate(optimal,nmax,mass,v,tstep,rpab,c,prop)
     call ghwp1D_propagate(classical,nmax,mass,v,tstep,rpab2,c2,prop2)
     !
     energy = prop%epot+prop%ekin; energy2 = prop2%epot+prop2%ekin
     !
     if (modulo(loop,100)==0) then
       write(10,fmt='(F8.2,F12.8,F12.8,F12.8)') loop*tstep,prop%r,prop%p,energy 
       write(11,fmt='(F8.2,F12.8,F12.8,F12.8)') loop*tstep,prop2%r,prop2%p,energy2
       if (modulo(loop,200)==0) then
         write(13,fmt='(A,I8)') '#', loop
         x = 1.51_rk
         do
           x = x + 0.01_rk
           call wavefunction_density(x,nmax,rpab (:,1),c (:,1),density)
           call wavefunction_density(x,nmax,rpab2(:,1),c2(:,1),density2)
           write(13,fmt='(F8.2,F12.8,F12.8)') x,0.3*density+energy,0.3*density2+energy2
           if (x>5._rk) exit
         end do
         write(13,fmt='(A1)')
         write(13,fmt='(A1)')
       end if 
     end if
     !
     !
     if (modulo(loop,100000)==0) write(*,*) 'elapsed time = ',int(loop*tstep),'fs'
     if (loop*tstep > 4000) exit
   end do
   !
   !
   close(10)
   close(11)
   !
   write(*,*) 'completed'
   !
 end program main
