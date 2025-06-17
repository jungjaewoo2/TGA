      subroutine sucstp(istep, nstep, icount, neq2, icdas, cdasr, cypre,
     &                  cyprer, idas, dasr, ypre)

c-----------------------------------------------------------------------
c Purpose:
c       This routine stores the necessary information after each
c       accepted time step (i.e. no need remeshing). These information
c       is needed if a remeshing is required next step.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 17, 2001.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input:
        integer                 istep
c                               istep is the number of time steps that
c                               DASSL has taken when using dassl.
c
        integer                 nstep
c                               nstep is the number of previous steps
c                               necessary.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 neq2
c                               neq2 is the number of bspline
c                               coefficients when using dassl_kcol+1.
c
        integer                 icdas(20)
c                               icdas stores the first 20 elements of
c                               the integer work array in dassl.
c
        double precision        cdasr(40)
c                               cdasr stores the first 40 elements of
c                               the floating point work array in dassl.
c
        double precision        cypre(neq2)
c                               cypre is the vector of bspline
c                               coefficients at the current step.
c
        double precision        cyprer(nstep*neq2)
c                               cyprer is the vector of bspline
c                               coefficients at the previous steps.
c
c       output:
        integer                 idas(20)
c                               idas is a copy of icdas.
c
        double precision        dasr(40)
c                               dasr is a copy of cdasr.
c
        double precision        ypre(6*neq2)
c                               ypre1 stores the bspline coefficients
c                               at the past 6 steps of dassl_kcol+1.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 imod
c
c-----------------------------------------------------------------------
c Loop Indices:
        integer                 i
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c Fortran Functions used:
c                               mod
c
c-----------------------------------------------------------------------

      do 10 i = 1, 20
         idas(i) = icdas(i)
   10 continue

      call dcopy(40, cdasr, 1, dasr, 1)

      imod = mod(istep, 6)
      if (imod .eq. 0) then
         imod = 6
      endif
      imod = 7 - imod

      call dcopy(neq2, cypre, 1, ypre((imod-1)*neq2+1), 1)

      if (icount .eq. 0) goto 99

c-----------------------------------------------------------------------
      if (nstep .eq. 6) then
         nstep = 5
      endif

      if (imod .eq. 6) then
         call dcopy(nstep*neq2, cyprer, 1, ypre, 1)
      else
         if ((imod+nstep) .le. 6) then
            call dcopy(nstep*neq2, cyprer, 1, ypre(imod*neq2+1), 1)
         else
            call dcopy((6-imod)*neq2, cyprer, 1, ypre(imod*neq2+1), 1)
            call dcopy((nstep+imod-6)*neq2, cyprer((6-imod)*neq2+1), 1,
     &                 ypre, 1)
         endif
      endif

c-----------------------------------------------------------------------
   99 continue

      return
      end
