      subroutine revalu(kcol, xsol, nint, x, npde, npts, nstep, usol, 
     &                  y, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the solution u at the npts points xsol 
c       and at the current and previous nstep-1 time step. Then return
c       them in the array usol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Jan 7, 2002.
c
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
        integer                 nstep
c                               nstep-1 is the number of previous steps.
c                               When user wants to calculate solution
c                               at tout, let nstep = 1.
c
        double precision        y(npde*(nint*kcol+nconti), nstep)
c                               y is the vector of bspline 
c                               coefficients at the current time step
c                               and previous nstep-1 steps.
c
c       output:
        double precision        usol(npde, npts, nstep)
c                               usol is the solution at the given 
c                               points and at the current time step 
c                               and previous nstep-1 steps.
c
c       Work Storage:
        double precision        work((kcol+nconti)+kcol*(nint+1)
     *                               +2*nconti)
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
        integer                 k
        integer                 m
        integer                 n
        integer                 ii
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
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
   20 continue
   30 continue

      do 70 n = 1, nstep
         do 60 i = 1, npts
c
c     interv is called to compute ileft. bbsplvd is called to compute
c     the values of the basis function at the required point.
c               write (*,*) '--1'
            call binterv(work(ixbs), ncpts, xsol(i), ileft, mflag)
c                   write (*,*) '--2'
            call bbsplvd(work(ixbs),kcol+nconti,xsol(i),ileft,work,1)
c                   write (*,*) '--3'
            ii = ileft - kcol - nconti
            do 50 k = 1, npde
               usol(k,i,n) = zero
               do 40 m = 1, kcol + nconti
                  mm = (m + ii - 1) * npde  
                  usol(k,i,n) = usol(k,i,n) + y(mm+k,n) * work(m)
   40          continue
   50       continue
   60    continue
   70 continue
      return
      end
