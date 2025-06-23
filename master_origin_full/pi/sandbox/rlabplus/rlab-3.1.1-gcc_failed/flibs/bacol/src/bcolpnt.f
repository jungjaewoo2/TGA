      subroutine bcolpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the piecewise polynomial space breakpoint
c       sequence, and calculates the collocation point sequence.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 3, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points. 
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval. 
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x. 
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
c       Work Storage:
        double precision        work(kcol*kcol)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [a,b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i), 
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
c-----------------------------------------------------------------------
c Local Variables:  
        double precision        rho(mxkcol+1)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+1)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Generate the piecewise polynomial space breakpoint sequence.
      do 10 i = 1, kcol + nconti
         xbs(i) = x(1)
         xbs(i + ncpts) = x(nint + 1)
   10 continue
      do 30 i = 2, nint
         ii = (i - 2) * kcol + kcol + nconti
         do 20 j = 1, kcol
            xbs(ii + j) = x(i)
   20    continue
   30 continue

c-----------------------------------------------------------------------
c     Compute the Gaussian points.
      call gauleg(kcol, kcol*kcol, rho, wts, work, 2)

c     Define the collocation point sequence.
      xcol(1) = x(1)
      do 50 i = 1, nint
         ii = (i - 1) * kcol + 1
         do 40 j = 1, kcol
            xcol(ii + j) = x(i) + h(i) * rho(j)
   40    continue
   50 continue
      xcol(ncpts) = x(nint + 1)

      return
      end
