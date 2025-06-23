      subroutine meshsq(kcol, nint, x, work, h, excol, ewts)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mesh size sequence, then generates
c       the collocation points and Gaussian weights for error
c       estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 5, 2001.
c
c-----------------------------------------------------------------------
c Constants:
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
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
c       Work Storage:
        double precision        work((kcol+3)*(kcol+3))
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        double precision        excol(nint*(kcol+3))
c                               excol is the collocation point sequence
c                               which is used for error estimate.
c
        double precision        ewts(nint*(kcol+3))
c                               ewts is the gaussian weight sequence
c                               which is used for error estimate.
c
c-----------------------------------------------------------------------
c Local Variables:  
        double precision        rho(mxkcol+3)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+3)
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
c     Calculate the mesh step size sequence.
      do 10 i = 1, nint
         h(i) = x(i+1)-x(i)
   10 continue

c     Compute the Gaussian points and Gaussian weights.
      call gauleg(kcol+3, (kcol+3)*(kcol+3), rho, wts,
     &            work, 4)

c     Define the collocation point sequence.
      do 30 i = 1, nint
         ii = (i - 1) * (kcol + 3)
         do 20 j = 1, kcol+3
            excol(ii + j) = x(i) + h(i) * rho(j)
            ewts(ii + j) = h(i) * wts(j)
   20    continue
   30 continue

      return
      end
