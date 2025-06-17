      subroutine remesh(istart, icount, nintmx, ninpre, ninold, errrat,
     &                  errint, irshfg, xold, nint, kcol, x, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates a new mesh by equidistributing the error
c       in each subinterval.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 22, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        double precision        point5
        parameter              (point5 = 0.5d0)
c
        double precision        one
        parameter              (one    = 1.0d0)
c
        double precision        two
        parameter              (two    = 2.0d0)
c
        double precision        saffac
        parameter              (saffac = 0.2d0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 nintmx
c                               the maximal number of subintervals that
c                               the user requires.
c
        integer                 ninpre
c                               ninpre is the number of subintervals
c                               when icount = 0 before remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        errint(ninold)
c                               errint is the error estimate at
c                               each subintervals.
c
c       Output:
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial call or continuation
c                                           calls;
c                                      = 1, remesh with a hot start.
c                                      = 2, remesh with a cold start.
c
        double precision        xold(ninpre+1)
c                               xold is the spatial mesh when icount = 0
c                               before remeshing.
c
c       In-output: 
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               ninmx >= nint >= 1.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh. As input, it is
c                               the value before remeshing; as output,
c                               it is the value after remeshing. 
c
c       Work storage:
        double precision        work(2*ninold+1)
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        aerr
        double precision        berr
c
c       Pointers into the floating point work array:
        integer                 ierror
c                               work(ierror-1+i) is the L2-norm error
c                               estimate at the first i subintervals.
c
        integer                 ixold
c                               work(ixold) contains a copy of mesh
c                               points before remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
c
c-----------------------------------------------------------------------
c Functions used:
c                               dble
c                               int
c
c-----------------------------------------------------------------------
      
c     Set the pointers into the floating point work array.
      ierror = 1
      ixold = ierror + ninold

c-----------------------------------------------------------------------
c     Update xold.
      if (icount .eq. 0) then
         do 10 i = 1, ninpre + 1
            xold(i) = x(i)
   10    continue
      endif

c-----------------------------------------------------------------------
c     Update icount, irshfg and nint.
      icount = icount + 1

c     If this is the first remesh at the current step which is not the
c     initial step.
      if ((icount .eq. 1) .and. (istart .eq. 1)) then
         irshfg = 1
         goto 20
      endif

c     If after four hot start the code still can not satisfy the error
c     requirement, a cold start will take place.
      If ((icount .eq. 5) .and. (istart .eq. 1)) then
         irshfg = 2
         nint = ninpre
         goto 20
      endif
 
c     Update errrat.
      errrat = (errrat/saffac) ** (one/dble(kcol+2))

c     Set the upper bound and lower bound of the ratio of nint over
c     ninold.
      if (errrat .gt. two) then
         errrat = two
      else
         if (errrat .lt. point5) then
            errrat = point5
         endif
      endif

      nint = int(ninold * errrat)

c     The code does not allow nint = ninold.
      if (nint .eq. ninold) then
         nint = nint + 1
      endif

   20 continue

c-----------------------------------------------------------------------
c     Update work(ixold) to be the mesh before remeshing.
      do 30 i = 1, ninold + 1
         work(ixold-1+i) = x(i)
   30 continue

c-----------------------------------------------------------------------
c     Store work(i) to be the sum of the error at the first i
c     subintervals.
      work(ierror) = errint(1)
      do 40 i = ierror-1+2, ninold
         work(i) = errint(i) + work(i-1)
   40 continue

c     Let aerr to be the mean value of errint(i).
      aerr = work(ninold)/dble(nint)

c     Equidistribute the mesh points.
      berr = aerr
      j = 1

      do 60 i = 2, nint
   50    continue
         if (berr .gt. work(j)) then
            j = j + 1
            goto 50
         else
            if (j .eq. 1) then
               x(i) = work(ixold) + (work(ixold-1+2) - work(ixold))
     &                * berr/work(1)
            else
               x(i) = work(ixold-1+j) + (work(ixold-1+j+1) -
     &                work(ixold-1+j)) * (berr - work(j-1))/errint(j)
            endif
         endif
         berr = berr + aerr
   60 continue

      x(1) = work(ixold)
      x(nint+1) = work(ixold-1+ninold+1)

      return
      end
