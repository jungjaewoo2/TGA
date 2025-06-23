      subroutine calres(npde, kcol, nint, ncpts, neq, xcol, abdblk,
     &                  fbasis, t, y, yprime, work, delta)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by res. It provides a lower-level
c       interface to generate the residue at the current time t.
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
        double precision        negone
        parameter              (negone = -1.0D0)
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
        integer                 neq
c                               neq=npde*ncpts is the number of bsplines 
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               T is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(neq)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
c       Work storage:
        double precision        work(4*npde+2*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+2*npde*npde.
c
c       Output:
        double precision        delta(neq)
c                               delta is the residual of the DAE system.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 m
        integer                 k
c
c       Indices:
        integer                 ii
        integer                 jj
        integer                 mm
        integer                 kk
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
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
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bac1bxa
c                               bac1bxb
c                               f
c                               eval
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c BLAS Subroutines:
c       double precision:
c                               dscal
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde

c-----------------------------------------------------------------------
c     Initialize the residual to the zero vector.
      do 10 i = 1, neq 
         delta(i) = zero
   10 continue

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations.

      do 30 i = 1, nint

c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol

         do 20 j = 1, kcol

c           jj is the pointer of collocation point.
            jj = (i - 1) * kcol + j + 1

c           mm is the pointer of delta.
            mm = (jj - 1) * npde + 1

c           kk is the pointer of the basis function values at
c           the current collocation point.
            kk =(jj-1)*(kcol+nconti)*3+1

c           Generate the approximate solution and its spatial 
c           derivatives at the current collocation point.
            call eval(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(kk),y)
c                write (*,*) 'calres: work(iu) = ', work(iu)
c           Evaluate the function f defining the PDE at the current 
c           collocation point, storing the result in delta.
            call pde1f(t, xcol(jj), work(iu), work(iux),
     &              work(iuxx), delta(mm), npde)


   20    continue
   30 continue

c     Scale (delta(i), i=npde+1,npde*(ncpts-1)) with negative one.    
      call dscal(npde*kcol*nint, negone, delta(npde+1), 1)

c-----------------------------------------------------------------------
c     Calculate the portion of the residual vector which depends on the
c     collocation matrix. delta := delta + A*Yprime.

c     Calculate (delta(i), i=1, npde), which depend on the left
c     boundary point.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call bac1bxa(t, work(iu), work(iux), delta(1), npde)

c-----------------------------------------------------------------------
c     Calculate (delta(i), i=1, npde), which depend on the right
c     boundary point.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call bac1bxb(t, work(iu), work(iux), delta(neq-npde+1), npde)

c-----------------------------------------------------------------------
c     Calculate (delta(i), i = npde+1, (ncpts-1)*npde), which depend
c     on the nint blocks in the middle of the collocation matrix A.
      do 70 i = 1, nint
         do 60 j = 1, kcol + nconti
            do 50 k = 1, kcol
               kk = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(k-1)*npde
               do 40 m = 1, npde
                  ii = npde+(i-1)*npde*kcol+(k-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  delta(ii) = delta(ii) + abdblk(kk) * yprime(mm)
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
      return
      end
