      subroutine divdif(neq, nstep, psi, work, y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the divided difference, which is required
c       by DASSL for a hot start, after calculating the bspline 
c       coefficients at the last nstep steps.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 4, 2001.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of bspline 
c                               coefficients after remeshing.
c
        integer                 nstep
c                               nstep is the number of time steps
c                               on which the remeshing is needed.
c 
        double precision        psi(6)
c                               psi is the stepsize vector of the
c                               previous 6 time steps.
c
c       Work Storage:
        double precision        work(6)
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        y(nstep*neq)
c                               y is the vector of bspline coefficients
c                               at the last nstep time steps after
c                               remeshing. 
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 m
c
c-----------------------------------------------------------------------

c     Update y to be the divide difference.

      do 40 i = 1, nstep - 1
         work(1) = psi(i)
         if (nstep .gt. 2) then
            do 10 j = 2, nstep-i
               work(j) = psi(j+i-1) - psi(j-1)
   10       continue
         endif
         do 30 j = nstep, i+1, -1
            do 20 m = 1, neq
               y((j-1)*neq+m) = (y((j-2)*neq+m) - y((j-1)*neq+m))
     *                          / work(j-i)
   20       continue
   30    continue
   40 continue

      work(1) = psi(1)  
      do 50 i = 2, nstep - 1
         work(i) = work(i-1) * psi(i)
   50 continue
      do 70 i = 2, nstep 
         do 60 m = 1, neq
            y((i-1)*neq+m) = y((i-1)*neq+m) * work(i-1)
   60    continue
   70 continue

c-----------------------------------------------------------------------
      return
      end
