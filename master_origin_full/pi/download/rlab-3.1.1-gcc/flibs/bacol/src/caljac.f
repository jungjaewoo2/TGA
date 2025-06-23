      subroutine caljac(npde, kcol, nint, ncpts, neq, xcol, abdtop, 
     &                  abdblk, abdbot, fbasis, t, y, yprime, cj, 
     &                  work, pd)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by jac. It provides a lower-level
c       interface to generate the iteration matrix.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 13, 2001.
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
c                               kcolis the number of collocation points
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
c                               neq=npde*ncpts is the number of bspline 
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdtop(npde*npde*nconti)
c                               abdtop stores the top block of the ABD
c                               matrices.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
        double precision        abdbot(npde*npde*nconti)
c                               abdbot stores the bottom block of the 
c                               ABD matrices.
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
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(neq)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        double precision        cj
c                               cj is a scalar chosen by DASSL to
c                               accelerate convergence of the modified
c                               Newton iteration used to solve the
c                               implicit equations resulting from the
c                               BDF methods.
c
c       Work storage:
        double precision        work(4*npde+5*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+5*npde*npde.
c
c       Output:
        double precision        pd(npde*npde*(2*nconti
     *                             +nint*kcol*(kcol+nconti)))
c                               pd is the ABD Jacobian (iteration)
c                               matrix of the residual of the DAE 
c                               system defined by RES.
c
c-----------------------------------------------------------------------
c Local Variables:  
        integer                 ipdtop
c                               ipdtop is the pointer into pd where the
c                               top block of the ABD Jacobian is stored.
c
        integer                 ipdblk
c                               ipdblk is the pointer into pd where the
c                               nint blocks in the middle of the ABD 
c                               Jacobian are stored.
c
        integer                 ipdbot
c                               ipdbot is the pointer into pd where the
c                               bottom block of the ABD Jacobian is 
c                               stored.
c
        integer                 nsiztb
c                               nsiztb is the size of the top block
c                               as same as the bottom block of the ABD
c                               Jacobian.
c
        integer                 nsizbk
c                               nsizbk is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 n
c
        integer                 ii
        integer                 ij
        integer                 jj
        integer                 kk
        integer                 nn
        integer                 mm
        integer                 jk
        integer                 jk2
        integer                 jk3
        integer                 mn
        integer                 mn2
        integer                 mn3
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
        integer                 idfdu
c                               work(idfdu) stores the Jacobian of f
c                               with respect to u.
c
        integer                 idfdux
c                               work(idfdux) stores the Jacobian of f
c                               with respect to u_x.
c
        integer                 idfuxx
c                               work(idfuxx) stores the Jacobian of f
c                               with respect to u_xx.
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
c Subroutines Called:
c                               pde1fjac
c                               bac1dbxa
c                               bac1dbxb
c                               eval
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               daxpy
c
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde
      idfdu  = iuxx   + npde
      idfdux = idfdu  + npde * npde
      idfuxx = idfdux + npde * npde
      idbdu  = idfuxx + npde * npde
      idbdux = idbdu  + npde * npde
      idbdt  = idbdux + npde * npde

c     Set the indices into pd which define the ABD Jacobian.
      ipdtop = 1
      ipdblk = ipdtop + nconti * npde * npde
      ipdbot = ipdblk + nint * npde * npde * kcol * (kcol + nconti)

c-----------------------------------------------------------------------
c     Calculate the size of top (or bottom) block and the size of a
c     subblock in the middle.
      nsiztb = npde * npde * nconti
      nsizbk = npde * npde * kcol * (kcol + nconti)

c     Initialize pdtop, pdblk and pdbot to zero.
      do 10 i = 1, nsiztb
         pd(ipdtop-1+i) = zero
         pd(ipdbot-1+i) = zero
   10 continue
      do 20 i = 1, nint * nsizbk
         pd(ipdblk-1+i) = zero
   20 continue
c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations and
c     caluculate the portion of dG/dY which depends on them.

      do 70 i = 1, nint

c        ii+1 is the pointer to the first element at the i-th subblock
c        of the jacobian matrix, i = 1, nint.
         ii = ipdblk - 1 + (i - 1) * nsizbk

c        ij is the value of ileft for the current collocation point.
         ij = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c           jj+1 is the pointer to the first element corresponding to
c           the j-th collocation point in the i-th interval.
            jj = ii + (j - 1) * npde

c           mm is the index of the current collocation point.
            mm = (i - 1) * kcol + j + 1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ij,mm,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1+(mm-1)*(kcol+nconti)*3),y)

c           Generate dfdu, dfdux, and dfdux at the current
c           collocation point (the j-th point of the i-th
c           subinterval).
            call pde1fjac(t, xcol(1+(i-1)*kcol+j), work(iu),
     &                  work(iux), work(iuxx), work(idfdu), 
     &                  work(idfdux), work(idfuxx), npde)

            do 50 k = 1, kcol + nconti

c              kk+1 is the pointer to the first element of a npde by
c              npde submatrix, which is corresponding to the j-th 
c              collocation point in the i-th interval, and the k-th 
c              nonzero basis function.
               kk = jj + (k-1) * npde * npde * kcol

c              jk is the pointer to the k-th nonzero function at the
c              mm-th collocation point in the basis function,
c              fbasis(1).
               jk = (mm - 1) * (kcol + nconti) * 3 + k

c              jk2 is the pointer to the first derivative for the 
c              above basis function.
               jk2 = jk + kcol + nconti

c              jk3 is the pointer to the second derivative for the 
c              above basis function.
               jk3 = jk2 + kcol + nconti

               do 40 m = 1, npde
                  do 30 n = 1, npde

c                    nn is the pointer to the (n, m) element of the
c                    npde by npde submatrix.
                     nn = kk + (m-1)*npde*kcol + n

c                    mn is the pointer to the (n, m) element of dfdu.
                     mn = idfdu - 1 + (m - 1) * npde + n

c                    mn2 is the pointer to the (n, m) element of dfdux.
                     mn2 = mn + npde * npde

c                    mn3 is the pointer to the (n, m) element of dfduxx.
                     mn3 = mn2 + npde * npde

c                    now set up the value in pd at the place nn.
                     pd(nn) = - work(mn) * fbasis(jk) 
     &                        - work(mn2) * fbasis(jk2)
     &                        - work(mn3) * fbasis(jk3)

   30             continue
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Update the values at the left boundary.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call bac1dbxa(t, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)

c     Update the top block of the collocation matrix dG/dY'.
      do 90 j = 1, npde
         do 80 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) =
     &            fbasis(2+kcol+nconti) * work(idbdux-1+mm)
            abdtop(ii) =
     &            work(idbdu-1+mm) - abdtop(jj)
   80    continue
   90 continue

c-----------------------------------------------------------------------
c     Update the values at the right boundary.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call bac1dbxb(t, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)

c     Update the bottom block of the collocation matrix.
      do 110 j = 1, npde
         do 100 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) =
     &            fbasis(1+kcol+kcol+nconti+(ncpts-1)*(kcol+nconti)*3)
     &            * work(idbdux-1+mm)
            abdbot(jj) =
     &            work(idbdu-1+mm) - abdbot(ii)
  100    continue
  110 continue

c-----------------------------------------------------------------------
c     Add cj * A to pdmat. This is cj * dG/dY'.

      call daxpy(nsiztb, cj, abdtop, 1, pd(ipdtop), 1)

      call daxpy(nint*nsizbk, cj, abdblk, 1, pd(ipdblk), 1)

      call daxpy(nsiztb, cj, abdbot, 1, pd(ipdbot), 1)

c-----------------------------------------------------------------------
      return
      end
