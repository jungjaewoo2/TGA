      subroutine reinit(npde, kcol, kold, nint, ninold, ncpts, neq,
     &                  neqold, icount, istep, nstep, x, xold, yold,
     &                  iflag, work, lw, ipivot, h, xbs, xcol,
     &                  fbasis, y, abdblk, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks after remeshing: 
c
c               calculating the mesh step size sequence,
c               generating the piecewise polynomial space breakpoint
c               sequence,
c               calculating the collocation point sequence,
c               calculating the B-spline basis functions,
c               constructing abdblk of the collocation matrices and
c               calculating the bspline coefficients at the last nstep
c               steps which is needed for a warm start.
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
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval after
c                               remeshing. 
c
        integer                 kold
c                               kold is the number of collocation points
c                               to be used in each subinterval before
c                               remeshing. 
c
        integer                 nint
c                               nint is the number of subintervals after
c                               remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals 
c                               before remeshing.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bspline 
c                               coefficients after remeshing.
c
        integer                 neqold
c                               neqold=npde*(kold*ninold+nconti) is
c                               the number of bspline 
c                               coefficients before remeshing.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 istep
c                               istep is the number of time steps that
c                               DASSL has taken when kcol collocation
c                               points are used in each subinterval.
c
        integer                 nstep
c                               nstep is the number of time steps
c                               on which the remeshing is needed.
c 
        double precision        x(nint+1)
c                               x is the spatial mesh after remeshing.
c
        double precision        xold(ninold+1)
c                               xold is the spatial mesh before
c                               remeshing.
c
        double precision        yold(6*neqold)
c                               yold is the vector of bspline
c                               coefficients at the last nstep time
c                               steps before remeshing. 
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        integer                 iflag
c                               iflag is a flag.
c                               iflag = 0, initialization is done for
c                                          dassl_kcol.
c                                       1, initialization is done for
c                                          dassl_kcol+1.
c 
        integer                 lw
c                               lw is the size of the work storage 
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +(kold+nconti)+kold*(ninold+1)
c                                     +2*nconti
c                               Since nint >= ninold/2 and kcol >=
c                               kold+1, it implies that lw >= 3*neqold.
c
c       Work Storage:
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
        integer                 ipivot(neq)
c                               pivoting information from the 
c                               factorization of the temporary matrix.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. 
c
        double precision        y(nstep*neq)
c                               y is the vector of bspline coefficients
c                               at the last nstep time steps after
c                               remeshing. 
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by bcrdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:  
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one 
c                               collocation block of work.
c 
        integer                 imod
c
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since bcrdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since bcrdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the 
c                               bottom block which is required since 
c                               bcrdcmp overwrites the input collocation
c                               matrix.
c
        integer                 ivwork
c                               work(ivwork) is the work storage
c                               required by values.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 jj
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bbsplvd
c                               bcolpnt
c                               bcrdcmp
c                               bcrslve
c                               revalu
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------
c Fortran Functions Used:
c                               mod
c
c-----------------------------------------------------------------------
c      write (*,*) '-1'
c     Generate the piecewise polynomial space breakpoint sequence,
c     and calculates the collocation point sequences.
      call bcolpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c     Update yold.
      if ((iflag .ne. 0) .or. (icount .ne. 1)) goto 5

      imod = mod(istep, 6)
      if (imod .eq. 0) then
         imod = 6
      endif
      imod = 7 - imod
      if (imod .eq. 1) goto 5

      if ((imod+nstep-1) .gt. 6) then
         if (imod .le. 3) then
            call dcopy((imod+nstep-1-6)*neqold, yold, 1, work, 1)
            call dcopy((6-imod+1)*neqold, yold((imod-1)*neqold+1), 1,
     &                 yold, 1)
            call dcopy((imod+nstep-1-6)*neqold, work, 1,
     &                 yold((6-imod+1)*neqold+1), 1)
         else
            call dcopy((6-imod+1)*neqold, yold((imod-1)*neqold+1), 1,
     &                 work, 1)
            call dcopy((imod+nstep-1-6)*neqold, yold, -1,
     &                 yold((6-imod+1)*neqold+1), -1)
            call dcopy((6-imod+1)*neqold, work, 1, yold, 1)
         endif
      else
         call dcopy(nstep*neqold, yold((imod-1)*neqold+1), 1, yold, 1)
      endif

    5 continue

c      write (*,*) '-2'

c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      ivwork = iabdbt + npde*npde*nconti

c-----------------------------------------------------------------------
c     Call revalu to calculate the values at xcol and at the last nstep
c     time step. Then save in y.
      call revalu(kold, xcol, ninold, xold, npde, ncpts, nstep,
     &            y, yold, work(ivwork))
c      write (*,*) '-3'
c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 10 i = 1, npde * npde * nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   10 continue
      do 20 i = 1, nint*nels
         abdblk(i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j) 
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bbsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c             write (*,*) '-4'

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 30 i = 1, npde 
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   30 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 70 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is the position in the y vector where the values for the 
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
            jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bbsplvd(xbs,kcol+nconti,xcol(ii),ileft,
     &                   fbasis(1,1,ii),3)

            do 50 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 40 m = 1, npde
                  mm = ll + (m-1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   40          continue
   50       continue
   60    continue
   70 continue
c                      write (*,*) '-5'
c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bbsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      do 80 i = 1, npde
         ii = ((i-1)+npde)*npde + i  
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   80 continue

c             write (*,*) '-6'

c-----------------------------------------------------------------------   
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nels*nint,abdblk,1,work(iabdbk),1)

c-----------------------------------------------------------------------
c     Generate the vector y.

c     LU decompose the matrix.
c             write (*,*) '-7'
      call bcrdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) go to 999

c     Solve the linear system. This gives the basis function 
c     coefficients for the initial conditions, i.e. y(t0).
      do 90 i = 1, nstep
         ii = (i - 1) * neq + 1
         call bcrslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &               (kcol+nconti)*npde,nint,work(iabdbt),npde,
     &               ipivot,y(ii),0)
         if (icflag .ne. 0) go to 999
   90 continue

c             write (*,*) '-8'

  999 return
      end
