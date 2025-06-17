      subroutine iniy(t0, npde, kcol, nint, neq, ncpts, ifglin, xcol,
     &                xbs, abdblk, fbasis, y, ipivot, work, lw, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks required by 
c       inital including:
c
c               calculating the Bspline basis functions,
c               constructing abdblk of the collocation matrices and
c               determining y(t0).
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
        integer                 ifglin
c                               ifglin is a flag for the boundary
c                               conditions.
c                               ifglin = 1, indicate both derichlet
c                                           boundary conditions;
c                                      = 0, else.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i), 
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
        integer                 lw
c                               lw is the size of the work storage 
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +2*neq+2*npde+2*npde*npde
c
c       Work Storage:
        integer                 ipivot(neq)
c                               pivoting information from the 
c                               factorization of the temporary matrix.
c
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
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
        integer                 idelta
c                               work(idelta) contains the residual which
c                               indicates how well y satisfies to the
c                               boundary condition and the initial
c                               condition at the internal collocation
c                               points.
c
        integer                 ivcol
c                               work(ivcol) contains the values of u
c                               at the internal collocation points.
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
c                               contains dbdux(npde,npde). That is,
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
        integer                 l
        integer                 m
        integer                 ii
        integer                 jj
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bac1bxa
c                               bac1bxb
c                               bbsplvd
c                               bac1dbxa
c                               bac1dbxb
c                               eval
c                               pde1uint
c                               bcrdcmp
c                               bcrslve
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c                               dscal
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      idelta = iabdbt + npde*npde*nconti
      ivcol  = idelta + neq
      iu     = ivcol  + neq-2*npde
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 20 i = 1, npde*npde*nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue
      do 30 i = 1, nint*nels
         abdblk(i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j) 
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bbsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Uinit is called to evaluate the first npde components at the
c     left boundary point, and save in y.
      call pde1uint(xcol(1), y(1), npde)

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 40 i = 1, npde 
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   40 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 80 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 70 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is the position in the y vector where the values for the 
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
            jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bbsplvd(xbs,kcol+nconti,xcol(ii),ileft,
     &                  fbasis(1,1,ii),3)
            call pde1uint(xcol(ii), y(jj), npde)

            do 60 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 50 m = 1, npde
                  mm = ll + (m -1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   50          continue
   60       continue
   70    continue
   80 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bbsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      ii = neq - npde + 1
      call pde1uint(xcol(ncpts), y(ii), npde)
      do 90 i = 1, npde
         ii = ((i-1)+npde)*npde + i  
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   90 continue

c-----------------------------------------------------------------------   
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)

c     Check whether both boundary conditions are derichlet boundary
c     conditions. If no, copy the values at the internal collocation
c     points to work(ivcol), which will be used for newton iterations.
      if (ifglin .eq. 0) then
         call dcopy(neq-2*npde,y(npde+1),1,work(ivcol),1)
         call dscal(neq-2*npde,negone,work(ivcol),1)
      endif

c-----------------------------------------------------------------------
c     Generate the initial vector y(t0).
c-----------------------------------------------------------------------

c     LU decompose the matrix.
      call bcrdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999
 
c     Solve the linear system. If derichlet boundary conditions are
c     given, this gives the basis function coefficients for the initial
c     conditions, i.e. y(t0). If not, this gives the predictor of y(t0).
      call bcrslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,y,0)
      if (icflag .ne. 0) goto 999
 
c     Check whether both boundary conditions are derichlet boundary
c     conditions.
      if (ifglin .eq. 1) goto 999

c-----------------------------------------------------------------------
c     Newton iteration loop.

c     Calculate (work(idelta-1+i), i = npde+1, neq-npde), which depends
c     on the nint blocks in the middle of the collocation matrix A.
      call dcopy(neq-2*npde,work(ivcol),1,work(idelta+npde),1)
      do 130 i = 1, nint
         do 120 j = 1, kcol + nconti
            do 110 l = 1, kcol
               ll = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(l-1)*npde
               do 100 m = 1, npde
                  ii = idelta-1+npde+(i-1)*npde*kcol+(l-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  work(ii) = work(ii) + abdblk(ll) * y(mm)
  100          continue
  110       continue
  120    continue
  130 continue
                                                                                
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)
 
c     Update the values at the left boundary.
      call eval(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y)
      call bac1bxa(t0, work(iu), work(iux), work(idelta), npde)
      call bac1dbxa(t0, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)
 
c     Set up the top block and save in work(iabdtp).
      do 150 j = 1, npde
         do 140 i = 1, npde
            ii = iabdtp - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            work(ii) = work(idbdu-1+mm) - work(jj)             
  140    continue
  150 continue  

c     Update the values at the right boundary.
      call eval(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y)
      call bac1bxb(t0, work(iu), work(iux), work(idelta+neq-npde),
     &             npde)
      call bac1dbxb(t0,work(iu),work(iux),work(idbdu),
     &            work(idbdux),work(idbdt),npde)
 
c     Set up the bottom block and save in work(iabdbt).
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = iabdbt - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i  
            work(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            work(jj) = work(idbdu-1+mm) - work(ii)
  160    continue
  170 continue 

c     LU decompose the matrix.
      call bcrdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

c                write (*,*) '(bcrdcmp) icflag = ', icflag

      if (icflag .ne. 0) goto 999
 
c     Solve the corrector equation.
      call bcrslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            work(idelta),0)

c                write (*,*) '(bcrslve) icflag = ', icflag

      if (icflag .ne. 0) goto 999
 
c     Now generate the corrector of y(t0).
      do 180 i = 1, neq
         y(i) = y(i) - work(idelta-1+i)
  180 continue

c-----------------------------------------------------------------------

  999 return
      end
