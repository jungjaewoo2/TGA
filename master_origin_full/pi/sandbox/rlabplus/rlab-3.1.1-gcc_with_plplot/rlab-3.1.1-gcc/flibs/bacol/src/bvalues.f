      subroutine bac1val(kcol, xsol, nint, x, npde, npts,
     &                  usol, y, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the solution u and the first nderv
c       derivatives of u at the npts points xsol. Then return them in 
c       the array usol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Dec. 23, 2002.
c Modified by Marijan Kostrun, project rlabplus, 3-28-2006
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector.
c
        double precision        xsol(npts)
c                               xsol is an arbitrary set of spatial
c                               points at which the solution and the
c                               first nderv derivative values are
c                               to be calculated.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a,x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        double precision        y(npde*(nint*kcol+nconti))
c                               y is the vector of bspline 
c                               coefficients at the final time step.
c
c       output:
        double precision        usol(npts,npde)
c                               usol is the solution and the spatial
c                               derivatives up to the nderiv-th
c                               derivative at the given points and at
c                               the final time step.
c
c       Work Storage:
        double precision        work((kcol+nconti)
     *                                 +kcol*(nint+1)+2*nconti)
c                               work is a floating point work storage
c                               array.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 mflag
c                               mflag is required by subroutine
c                               interv.
c
c       Pointers into the floating point work array:
        integer                 ixbs
c                               work(ixbs) contains the breakpoint
c                               sequence.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 ii
        integer                 mj
        integer                 mm
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bbsplvd
c                               interv
c
c-----------------------------------------------------------------------
    
c     set up the value for ileft, mflag and ncpts.
      ileft = 0
      mflag = 0
      ncpts = nint * kcol + nconti

c     set the pointer into the floating point work array
      ixbs  = (kcol+nconti) + 1

c     Store the piecewise polynomial space breakpoint sequence in 
c     work(ixbs).
c
      do 10 i = 1, kcol + nconti
         work(ixbs-1+i) = x(1)
         work(ixbs-1+i+ncpts) = x(nint+1)
c          write (*,*) 'x(1) = ', work(ixbs-1+i)
c          write (*,*) 'x(',nint+1,') = ', work(ixbs-1+i+ncpts)
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
c             write (*,*) 'x(',i,') = ', work(ixbs-1+ii+k)
   20 continue
   30 continue

      do 70 i = 1, npts

c                 write (*,*) 'i1 = ', i
c                 write (*,*) 'xsol(',i,') = ', xsol(i)
c
c     interv is called to compute ileft. bbsplvd is called to compute
c     the values of the basis function at the required point.
         call binterv(work(ixbs), ncpts, xsol(i), ileft, mflag)
c              write (*,*) 'i2 = ', i
         call bbsplvd(work(ixbs),kcol+nconti,xsol(i),ileft,work,1)
c             write (*,*) 'i3 = ', i
         ii = ileft - kcol - nconti
         do 50 k = 1, npde
           usol(i,k) = zero
             do 40 m = 1, kcol + nconti
               mm = (m + ii - 1) * npde
               usol(i,k) = usol(i,k) + y(mm+k) * work(m)
   40          continue
   50      continue
   70 continue
      return
      end
