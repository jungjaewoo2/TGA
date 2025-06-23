      subroutine berrest(kcol, nint, npde, neq1, neq2, npts, icount,
     &                  xsol, wts, xbs1, xbs2, y1, y2, istart, mflag2,
     &                  atol, rtol, lenwk, work, errba1, errba2,
     &                  errrat, errint, errcom, ieflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the error estimate at each subinterval
c       and for each component of PDEs, and decides whether a remeshing
c       is necessary or not.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 nintsm
        parameter              (nintsm = 15)
c                               when the current step is the first step
c                               after remeshing, we require
c                               if nint <= nintsm
c                                    errrat < saffa2
c                               else
c                                    saffa1 < errrat < saffa2.
c                               endif
c
        double precision        zero
        parameter              (zero = 0.0d0)
c
        double precision        one
        parameter              (one = 1.0d0)
c
        double precision        two
        parameter              (two = 2.0d0)
c
        double precision        saffa1
        parameter              (saffa1 = 0.1d0)
c
        double precision        saffa2
        parameter              (saffa2 = 0.4d0)
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
        integer                 neq1
c                               neq1=npde*(nint*kcol+nconti) is the
c                               number of bspline coefficients (or
c                               DAEs) when using dassl_kcol.
c
        integer                 neq2
c                               neq2=neq1+npde*nint is the number of
c                               bspline coefficients (or DAEs) when
c                               using dassl_kcol+1.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector, which is equal to
c                               nint*(kcol+3).
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xsol(npts)
c                               xsol is the npts Gauss-Legend
c                               points at which the solution are
c                               to be calculated.
c
        double precision        wts(npts)
c                               wts is the npts Gauss-Legend
c                               weights at the corresponding xsol.
c
        double precision        xbs1((kcol+1)*nint+nconti+nconti)
c                               xbs1 is the breakpoint sequence when
c                               using dassl_kcol.
c
        double precision        xbs2((kcol+2)*nint+nconti+nconti)
c                               xbs2 is the breakpoint sequence when
c                               using dassl_kcol+1.
c
        double precision        y1(neq1)
c                               y1 is the vector of bspline
c                               coefficients when using dassl_kcol.
c
        double precision        y2(neq2)
c                               y2 is the vector of bspline
c                               coefficients when using dassl_kcol+1.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 mflag2
c                               mflag2 = 0, scalar atol and rtol.;
c                               mflag2 = 1, vector atol and rtol.
c
        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        integer                 lenwk
c                               lenwk is the size of the work storage
c                               array and must satisfy:
c                               lenwk >= 2*npde*nint*(kcol+3)
c                                        +npde*nint
c
c       Work Storage:
        double precision        work(lenwk)
c                               work is a floating point work storage
c                               array of size lenwk.
c
c       output:
        double precision        errba1((kcol+nconti)*npts)
c                               errba1 is the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol.
c
        double precision        errba2((kcol+1+nconti)*npts)
c                               errba2 is the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol+1.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of errcom.
c
        double precision        errint(nint)
c                               errint is the error estimate at each
c                               subinterval.
c
        double precision        errcom(npde)
c                               errcom is the error estimate for
c                               each component of pdes at the whole
c                               range, i.e. from x_a to x_b.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        errsum
c                               errsum is the sum of errint.
c
        double precision        errmax
c                               errmax is the maximum value of
c                               errint(i), i = 1, nint.
c
        double precision        aerr
c                               aerr is the average value of errint(i),
c                               i = 1, nint.
c
        double precision        disind
c                               disind is equal to errmax/aerr, and it
c                               indicates the error distribution over
c                               the mesh.
c
c       Pointers into the floating point work array:
        integer                 iusol1
c                               work(iusol1) stores the values at the
c                               npts points when using dassl_kcol.
c
        integer                 iusol2
c                               work(iusol2) stores the values at the
c                               npts points when using dassl_kcol+1.
c
        integer                 ierrci
c                               work(ierrci) stores the error estimate
c                               at each subinterval for each component.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j, m, ij, im, mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               berrval
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      iusol1 = 1
      iusol2 = iusol1 + npde * nint * (kcol + 3)
      ierrci = iusol2 + npde * nint * (kcol + 3)

c-----------------------------------------------------------------------
c     Generate the different values at the npts points xsol, and save
c     in work(iusol1) and work(iusol2).
c          write (*,*) 'errest: 1'
      call berrval(kcol, nint, npde, neq1, kcol+3, istart, icount,
     &            xbs1, xsol, y1, errba1, work(iusol1))
c          write (*,*) 'errest: 2'
      call berrval(kcol+1, nint, npde, neq2, kcol+3, istart, icount,
     &            xbs2, xsol, y2, errba2, work(iusol2))
c          write (*,*) 'errest: 3'

c-----------------------------------------------------------------------
c     Initialization task.
      do 10 i = 1, nint
         errint(i) = zero
   10 continue

      do 20 i = 1, npde
         errcom(i) = zero
   20 continue

      do 30 i = 1, npde * nint
         work(ierrci - 1 + i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Calculate the error estimate at each subinterval for each
c     component of PDEs.
      if (mflag2 .eq. 0) then
         do 60 m = 1, npde
            do 50 i = 1, nint
               do 40 j = 1, kcol + 3
                  ij = (i - 1) * (kcol + 3) + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im) + ((work(iusol1-1+mm)
     &                       - work(iusol2-1+mm)) /
     &                       (atol(1) + rtol(1)*abs(work(iusol1-1+mm))))
     &                       **2 * wts(ij)
   40          continue
   50       continue
   60    continue
      else
         do 90 m = 1, npde
            do 80 i = 1, nint
               do 70 j = 1, kcol + 3
                  ij = (i - 1) * (kcol + 3) + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im) + ((work(iusol1-1+mm)
     &                       - work(iusol2-1+mm)) /
     &                       (atol(m) + rtol(m)*abs(work(iusol1-1+mm))))
     &                       **2 * wts(ij)
   70          continue
   80       continue
   90    continue
      endif

c-----------------------------------------------------------------------
c     Calculate errint and errcom.
      do 110 j = 1, npde
         do 100 i = 1, nint
            ij = ierrci - 1 + (j - 1) * nint + i
            errint(i) = errint(i) + work(ij)
            errcom(j) = errcom(j) + work(ij)
  100    continue
  110 continue

c     Take the square root and update errint and errcom.
      do 120 i = 1, nint
         errint(i) = sqrt(errint(i))
         errint(i) = errint(i) ** (one/dble((kcol+2)))
c        errint(i) = errint(i) ** (one/dble(2*(kcol+2)))
  120 continue

      do 130 i = 1, npde
         errcom(i) = sqrt(errcom(i))
  130 continue

c-----------------------------------------------------------------------
c     Decide whether remeshing is needed.
      ieflag = 0

c     update errrat.
      errrat = zero
      do 140 i = 1, npde
         if (errcom(i) .gt. errrat) then
            errrat = errcom(i)
         endif
  140 continue

c     Calculate errsum to be the sum of the errint. Find the maximum
c     errint(i) and save it in errmax.
      errsum = errint(1)
      errmax = errint(1)
      do 150 i = 2, nint
         if (errmax .lt. errint(i)) errmax = errint(i)
         errsum = errint(i) + errsum
  150 continue

c     Let aerr be the mean value of errint(i).
      aerr = errsum/dble(nint)

c     Calculate disind.
      disind = errmax/aerr
c     print *, 'the rate is ', disind, ',   errrat = ', errrat

      if (disind .gt. two) then
         ieflag = 1
      else
         if ((istart .ne. 1) .or. (icount .ne. 0)) then
            if (nint .gt. nintsm) then
               if ((errrat .ge. saffa2) .or. (errrat .le. saffa1))
     &            ieflag = 1
            else
               if (errrat .ge. saffa2) ieflag = 1
            endif
         else
            if (errrat .ge. one) ieflag = 1
         endif
      endif

      return
      end
