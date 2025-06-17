      subroutine eval(npde,kcol,ileft,icpt,ncpts,uval,uxval,uxxval,
     &                fbasis,y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine evaluates u(k), ux(k), and uxx(k), k=1 to npde,
c       at the icpt-th collocation point.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Feb. 11, 2001.
c
c-----------------------------------------------------------------------
c Constants
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 icpt
c                               the index of the collocation point.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        fbasis((kcol+nconti)*3)
c                               Basis function values at the icpt-th
c                               collocation point.
c                               fbasis(k+(j-1)*(kcol+nconti)) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti).
c
        double precision        y(ncpts*npde)
c                               y is the vector of bspline coefficients.
c
c       Output:
        double precision        uval(npde)
c                               uval gives the approximation to
c                               u(t,x).
c
        double precision        uxval(npde)
c                               uxval gives the approximation to
c                               the first spatial derivative of u(t,x).
c
        double precision        uxxval(npde)
c                               uxxval gives the approximation to
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 j
        integer                 m
c-----------------------------------------------------------------------
      do 10 j = 1, npde
         uval(j)   = zero
         uxval(j)  = zero
         uxxval(j) = zero
   10 continue
      if (icpt .ne. 1 .and. icpt .ne. ncpts) then
         do 30 j = 1, npde
            do 20 m = 1, kcol + nconti
               uval(j)   = uval(j) + fbasis(m) *
     &                     y((ileft-kcol-3+m) * npde + j)
               uxval(j)  = uxval(j) + fbasis(m+kcol+nconti) *
     &                     y((ileft-kcol-3+m) * npde + j)
               uxxval(j) = uxxval(j) + fbasis(m+2*(kcol+nconti)) *
     &                     y((ileft-kcol-3+m) * npde + j)
   20          continue
   30       continue
      else
         if (icpt .eq. 1) then
            do 40 j = 1, npde
               uval(j)   = uval(j)   + fbasis(1) * y(j)
               uxval(j)  = uxval(j)  + fbasis(1+kcol+nconti) * y(j)
     &                     + fbasis(2+kcol+nconti) * y(npde + j)
               uxxval(j) = uxxval(j) + fbasis(1+2*(kcol+nconti)) * y(j)
     &                     + fbasis(2+2*(kcol+nconti)) * y(npde + j)
     &                     + fbasis(3+2*(kcol+nconti)) * y(2*npde + j)
   40       continue
         else
            do 50 j = 1, npde
               uval(j)   = uval(j)   + fbasis(kcol+nconti)
     &                     * y((ncpts - 1) * npde + j)
               uxval(j)  = uxval(j)  + fbasis((kcol+nconti)*2)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*2-1)
     &                     * y((ncpts - 2) * npde + j)
               uxxval(j) = uxxval(j) + fbasis((kcol+nconti)*3)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*3-1)
     &                     * y((ncpts - 2) * npde + j)
     &                     + fbasis((kcol+nconti)*3-2)
     &                     * y((ncpts - 3) * npde + j)
   50       continue
         endif
      endif
      return
      end
