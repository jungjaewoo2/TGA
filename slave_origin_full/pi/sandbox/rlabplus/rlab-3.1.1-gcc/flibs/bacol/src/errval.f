      subroutine berrval(kcol, nint, npde, neq, nptse, istart, icount,
     &                  xbs, xsol, y, errbas, usol)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the values of the (kcol+nconti) nonzero
c       bspline basis function at each Gaussian point of xsol.
c       Then determine the solution usol, which is used for error
c       estimate, at xsol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 29, 2001.
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
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(kcol*nint+2) is the number of
c                               bspline coefficients.
c
        integer                 nptse
c                               nptse is the number of Gaussian points
c                               in each subinterval for the error
c                               estimate.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xbs((kcol+1)*nint+nconti+nconti)
c                               The breakpoint sequence.
c
        double precision        xsol(nptse*nint)
c                               xsol is a set of spatial points at which
c                               the solution are to be calculated for
c                               error estimate.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients.
c
c       output:
        double precision        errbas(kcol+nconti, nptse*nint)
c                               errbas is the values of the nonzero
c                               basis functions at xsol.
c
        double precision        usol(npde, nptse*nint)
c                               uval is the solution at xsol.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bbsplvd
c
c-----------------------------------------------------------------------

c     check whether errbas is necessary to be calculated.
      if ((istart .eq. 1) .and. (icount .eq. 0)) goto 30

c     calculate errbas.
      do 20 i = 1, nint
         ileft = kcol + nconti + (i - 1) * kcol
         do 10 j = 1, nptse
            jj = (i - 1) * nptse + j
            call bbsplvd(xbs, kcol+nconti, xsol(jj), ileft,
     &                   errbas(1,jj), 1)
   10    continue
   20 continue

   30 continue

c     compute the values of usol at xsol.
      do 70 i = 1, nint
         do 60 j = 1, nptse
            jj = (i - 1) * nptse + j
            do 50 k = 1, npde
               usol(k,jj) = zero
               do 40 m = 1, kcol + nconti
                  mm = npde * (m + (i - 1) * kcol - 1) + k
                  usol(k,jj) = usol(k,jj) + y(mm) * errbas(m,jj)
   40          continue
   50       continue
   60    continue
   70 continue

      return
      end
