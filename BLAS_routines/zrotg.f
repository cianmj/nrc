      subroutine zrotg(ca,cb,c,s)
      double complex ca,cb,s
      double precision c
      double precision norm,scale
      double complex alpha
c     Linux hack:
c        Note that linux fortran doesn't support cdabs, so we wrote one;
c        hopefully, this won't slow things too much!  -TAA 95/08/15
      double precision cdabs
c
      if (cdabs(ca) .ne. 0.0d0) go to 10
         c = 0.0d0
         s = (1.0d0,0.0d0)
         ca = cb
         go to 20
   10 continue
         scale = cdabs(ca) + cdabs(cb)
         norm = scale*dsqrt((cdabs(ca/dcmplx(scale,0.0d0)))**2 +
     *                      (cdabs(cb/dcmplx(scale,0.0d0)))**2)
         alpha = ca /cdabs(ca)
         c = cdabs(ca) / norm
         s = alpha * dconjg(cb) / norm
         ca = alpha * norm
   20 continue
      return
      end


c     Linux hack:
c        Note that linux fortran doesn't support cdabs, so we wrote one;
c        hopefully, this won't slow things too much!  -TAA 95/08/15
      double precision function cdabs(z)
      double complex z
      cdabs=dsqrt(dble(z*dconjg(z)))
      return
      end
c
