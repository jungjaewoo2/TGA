      subroutine iniyp(t0, npde, kcol, nint, neq, ncpts, xcol,
     &                 abdtop, abdblk, abdbot, fbasis, y, yprime,
     &                 ipivot, work, lw, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks required by 
c       bacol including:
c
c               determining yprime(t0).
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
        double precision        t0
c                               t0 is the initial time.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval. 
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x. 
c                               nint >= 1.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bsplines 
c                               coefficients (or DAEs).
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a,x_b].
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. fbasis(k,j,i) contains the
c                               values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        y(neq)
c                               y = y(t0) is the initial vector of
c                               bspline coefficients. 
c
c       Output:
        double precision        abdtop(npde*npde*nconti)
c                               The first block of the matrix A.
c
        double precision        abdbot(npde*npde*nconti)
c                               The last block of the matrix A.
c
        double precision        yprime(neq)
c                               yprime = yprime(t0) is the initial 
c                               vector of bspline coefficients
c                               for the first temporal derivative. 
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by bcrdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        integer                 lw
c                               lw is the size of the work storage 
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint+
c                                     npde*3+2*npde*npde+npde.
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
c-----------------------------------------------------------------------
c Local Variables:  
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of abdtop
c                               which is required since bcrdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since bcrdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of abdbot
c                               which is required since bcrdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde, 
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde), That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bcrdcmp
c                               bcrslve
c                               eval
c                               bac1dbxa
c                               bac1dbxb
c                               pde1f
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + npde*npde*kcol*(kcol+nconti)*nint
      iu     = iabdbt + npde*npde*nconti
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde

c-----------------------------------------------------------------------
c     Initialize abdtop, abdbot and abdblk to zero.
      do 20 i = 1, npde * npde * nconti
         abdtop(i) = zero
         abdbot(i) = zero
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue

c-----------------------------------------------------------------------

c     Update the values at the left boundary.
      call eval(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y) 
      call bac1dbxa(t0, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde) 

c     Store -work(idbdt), which is the right side of the left boundary 
c     conditions, into yprime.
      do 100 i = 1, npde
         yprime(i) = - work(idbdt-1+i)
  100 continue

c     Set up the top block and save in abdtop.
      do 120 j = 1, npde
         do 110 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            abdtop(ii) = work(idbdu-1+mm) - abdtop(jj)
  110    continue
  120 continue

c-----------------------------------------------------------------------
c     Generate the right side of ODEs at the collocation points
c     and save in yprime(i), i = npde + 1, neq - npde.
      do 140 i = 1, nint

c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol

         do 130 j = 1, kcol

c           jj is the index of the current collocation point.
            jj = (i - 1) * kcol + j + 1

c           mm is the pointer of yprime.
            mm = (jj - 1) * npde + 1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1,1,jj),y)

c           Evaluate the function f defining the PDE at the current
c           collocation point, storing the result in yprime.
            call pde1f(t0,xcol(jj),work(iu),work(iux),work(iuxx),
     &                 yprime(mm),npde)

  130    continue
  140 continue

c----------------------------------------------------------------------- 
c     Update the values at the right boundary.
      call eval(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y) 
      call bac1dbxb(t0,work(iu),work(iux),work(idbdu),
     &            work(idbdux),work(idbdt),npde)

c     Store -work(idbdt), which is the right side of the right boundary 
c     conditions, into yprime.
      do 150 i = 1, npde
         ii = neq - npde + i
         yprime(ii) = - work(idbdt-1+i)
  150 continue

c     Set up the bottom block and save in abdbot.
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            abdbot(jj) = work(idbdu-1+mm) - abdbot(ii)
  160    continue
  170 continue

c-----------------------------------------------------------------------
c     Copy the collocation matrix into temporary storage.
      call dcopy(npde*npde*nconti, abdtop, 1, work(iabdtp), 1)
c
      call dcopy(npde*npde*kcol*(kcol+nconti)*nint, abdblk, 1,
     &           work(iabdbk), 1)
c
      call dcopy(npde*npde*nconti, abdbot, 1, work(iabdbt), 1)

c-----------------------------------------------------------------------
c     Generate the initial vector yp(t0).

c     LU decompose the matrix.

      call bcrdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)


      if (icflag .ne. 0) go to 999
c
c     Solve the linear system. This gives yprime(t0) 
      call bcrslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            yprime,0)

  999 continue
      return
      end
