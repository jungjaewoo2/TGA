c-----------------------------------------------------------------------
c
c  complex eigenvalue problem
c
c-----------------------------------------------------------------------
      subroutine zndrv1(v, ldv, workl, lworkl, workev,
     &                  workd, d, resid, ax, iselect, iparam, ipntr,
     &                  a, n, nev, ncv,
     &                  which, sigma, info, ierr, ipiv, rwork, rd)
c
      integer           iparam(11), ipntr(14)
      logical           iselect(ncv)
      Complex*16 
     &                  ax(n), d(n), a(n,n),
     &                  v(ldv,ncv), workd(3*n), 
     &                  workev(3*ncv), resid(n), 
     &                  workl(lworkl)
      Double precision  
     &                  rwork(ncv), rd(ncv,3)
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info,
     &                  ierr, ipiv(n)
      Complex*16 
     &                  sigma
      Double precision 
     &                  tol
      logical           rvec
c
      Double precision 
     &                  dznrm2 , dlapy2
      external          dznrm2 , zaxpy , dlapy2

      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)


      bmat  = 'I'

      tol    = zero
      ido    = 0
      info   = 0

c
c     calculate [sigma!=0] A = (A-sigma*I)^{-1}
c
      if (zabs(sigma) .ne. zero) then
c         a -> a-sigma*I
          do 9 i=1,n
            a(i,i) = a(i,i) - sigma
   9      continue
c         a -> P L U
          call zgetrf(n, n, a, n, ipiv, info)
          if (info .ne. 0 ) goto 9000
c         a -> a^{-1}
          call zgetri(n, a, n, ipiv, workd, 3*n, info)
          if (info .ne. 0 ) goto 9000
      end if

c
 10   continue
c
         call znaupd  ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork,info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
            call zarpax (a, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c
      if ( info .lt. 0 ) goto 9000
c
      rvec = .true.
c
      call zneupd  (rvec, 'A', iselect, d, v, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, ierr)
c
      if ( ierr .ne. 0) goto 9000
c
 9000 continue
c
      end

c-----------------------------------------------------------------------
c
c  complex eigenvalue problem  eigs(a,b,..)
c
c-----------------------------------------------------------------------
      subroutine zndrv2(v, workl, lworkl, workev,
     &                  workd, d, resid, iselect, iparam, ipntr,
     &                  a, b, n, nev, ncv,
     &                  which, sigma, info, ierr, ipiv, rwork, rd)
c
      integer           iparam(11), ipntr(14)
      logical           iselect(ncv)
      Complex*16 
     &                  d(ncv), a(n,n), b(n,n),
     &                  v(n,ncv), workd(3*n), 
     &                  workev(2*ncv), resid(n), 
     &                  workl(lworkl)
      Double precision  
     &                  rwork(n), rd(ncv,3)
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,
     &                  ierr, ipiv(n), ldv
      Complex*16 
     &                  sigma
      Double precision 
     &                  tol
      logical           rvec
c
      external          zarpax, zcopy, zgetrf, zgetri,
     &                  znaupd, zneupd
c
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
c
c
      bmat  = 'G'
c
      ldv   = n
c
      tol    = zero
      ido    = 0
      info   = 0
      ierr   = 0
c
c
c     calculate [sigma!=0] A = (A-sigma*I)^{-1}
c
      if (zabs(sigma) .ne. zero) then
c         a -> a-sigma*M
          do 9 i=1,n
            do 8 j=1,n
              a(i,j) = a(i,j) - sigma * b(i,j)
   8        continue
   9      continue
      endif
c
c     a -> P L U
      call zgetrf(n, n, a, n, ipiv, info)
      if (info .ne. 0 ) goto 9000
c     a -> a^{-1}
      call zgetri(n, a, n, ipiv, workd, 3*n, info)
      if (info .ne. 0 ) goto 9000
c
c
c
 10   continue
c
         call znaupd  ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, info )
c
         if (ido .eq. -1 ) then
c           (2) <- b . (1)
            call zarpax (b, n, workd(ipntr(1)), workd(ipntr(2)))
c           (2) -> (1)
            call zcopy( n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
c           (2) <- a . (1)
            call zarpax (a, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c
         if (ido .eq. 1 ) then
c           (2) <- a . (3)
            call zarpax (a, n, workd(ipntr(3)), workd(ipntr(2)))
            go to 10
         end if
c 
         if (ido .eq. 2 ) then
c           (2) <- b . (1)
            call zarpax (b, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c
      if ( info .lt. 0 ) goto 9000
c
      rvec = .true.
c
c
      call zneupd  (rvec, 'A', iselect, d, v, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv, v,
     &        ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, ierr)
c
c
 9000 continue
c
      
      end

c-----------------------------------------------------------------------
c
c  symmetric simple eigenvalue problem
c
c-----------------------------------------------------------------------

      subroutine dsdrv1(v, ldv, workl, lworkl,
     &                  workd, d, resid, ax, iselect, iparam, ipntr,
     &                  a, n, nev, ncv,
     &                  which, sigma, info, ierr, ipiv)
c-----------------------------------------------------------------------
c
      integer          ldv
      Double precision
     &                 v(ldv,ncv), workl(lworkl),
     &                 workd(3*n), d(ncv,2), resid(n),
     &                 ax(n), a(n,n)

      logical          iselect(ncv)

      integer          iparam(11), ipntr(11), ipiv(n)

      character        bmat*1, which*2

      integer          ido, nev, ncv, lworkl, info, ierr

      logical          rvec

      Double precision      
     &                 tol, sigma

      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
c  
      Double precision
     &                 dnrm2
      external         dnrm2, daxpy, dsaupd, dseupd, darpax,
     &                 dgetrf, dgetri

      intrinsic        abs

      bmat = 'I'

      tol = zero
      info = 0
      ido = 0

c
c     calculate [sigma!=0] A = (A-sigma*I)^{-1}
c
      if (sigma .ne. zero) then
c         a -> a-sigma*I
          do 9 i=1,n
            a(i,i) = a(i,i) - sigma
   9      continue
c         a -> P L U
          call dgetrf(n, n, a, n, ipiv, info)
          if (info .ne. 0 ) goto 9000
c         a -> a^{-1}
          call dgetri(n, a, n, ipiv, workd, 3*n, info)
          if (info .ne. 0 ) goto 9000
      end if
 
 10   continue

         call dsaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
         if (ido .eq. -1 .or. ido .eq. 1) then
            call darpax(a, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if

c
      if ( info .lt. 0 ) goto 9000

      rvec = .true.
c
      call dseupd ( rvec, 'All', iselect, d, v, ldv, sigma,
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
c
      if ( ierr .ne. 0) goto 9000
c
 9000 continue
      return
      end



c-----------------------------------------------------------------------
c
c  nonsymmetric simple eigenvalue problem  eigs(a,..)
c
c-----------------------------------------------------------------------

      subroutine dndrv1(v, ldv, workl, lworkl, workev,
     &                  workd, d, resid, ax, iselect, iparam, ipntr,
     &                  a, n, nev, ncv,
     &                  which, sigmar, info, ierr, ipiv)
c-----------------------------------------------------------------------
c
      integer           ldv
      integer           iparam(11), ipntr(14)
      logical           iselect(ncv)
c
      Double precision
     &                  a(n,n), ax(n), d(ncv,3), resid(n),
     &                  v(ldv,ncv), workd(3*n),
     &                  workev(3*ncv),
     &                  workl(3*ncv*ncv+6*ncv)
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info,
     &                  ierr
      Double precision
     &                  tol, sigmar, sigmai
      logical           rvec
c
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
c
c
      Double precision
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy, darpax
c
      intrinsic         abs
c
      bmat  = 'I'
c
c
      lworkl  = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
      sigmai = 0

c
c     calculate [sigma!=0] A = (A-sigma*I)^{-1}
c
      if (sigmar .ne. zero) then
c         a -> a-sigma*I
          do 9 i=1,n
            a(i,i) = a(i,i) - sigmar
   9      continue
c         a -> P L U
          call dgetrf(n, n, a, n, ipiv, info)
          if (info .ne. 0 ) goto 9000
c         a -> a^{-1}
          call dgetri(n, a, n, ipiv, workd, 3*n, info)
          if (info .ne. 0 ) goto 9000
      end if
 

 10   continue
         call dnaupd ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        info )
         if (ido .eq. -1 .or. ido .eq. 1) then
            call darpax (a, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c 
      if ( info .lt. 0 ) goto 9000
c
      rvec = .true.
c
      call dneupd ( rvec, 'A', iselect, d, d(1,2), v, ldv,
     &      sigmar, sigmai, workev, bmat, n, which, nev, tol,
     &      resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &      lworkl, ierr )
c
      if ( ierr .ne. 0) goto 9000
c
c
 9000 continue
c
      end


c-----------------------------------------------------------------------
c
c  nonsymmetric general eigenvalue problem  eigs(a,b,..)
c
c-----------------------------------------------------------------------

      subroutine dndrv2(v, ldv, workl, lworkl, workev,
     &                  workd, d, resid, iselect, iparam, ipntr,
     &                  a, b, n, nev, ncv,
     &                  which, sigma, info, ierr, ipiv)
c-----------------------------------------------------------------------
c
      integer           ldv
      integer           iparam(11), ipntr(14)
      logical           iselect(ncv)
c
      Double precision
     &                  a(n,n), b(n,n), d(ncv,3), resid(n),
     &                  v(ldv,ncv), workd(3*n),
     &                  workev(3*ncv),
     &                  workl(3*ncv*ncv+6*ncv)
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,
     &                  ierr
      Double precision
     &                  tol, sigma
      logical           rvec
c
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
c
c
      Double precision
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy, darpax
c
      intrinsic         abs
c
      bmat  = 'G'
      iparam(7) = 3
c
c
      lworkl  = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
c
c     calculate: [sigma!=0]&& a <- (a-sigma*b)^{-1}
c                [sigma==0]&& a <- (a)^{-1}
c
      if (sigma .ne. zero) then
c         a -> a-sigma*b
          do 9 i=1,n
            do 8 j=1,n
                  a(i,j) = a(i,j) - sigma * b(i,j)
   8        continue 
   9      continue
      end if
c
c     a -> P L U
      call dgetrf(n, n, a, n, ipiv, info)
      if (info .ne. 0 ) goto 9000
c     a -> a^{-1}
      call dgetri(n, a, n, ipiv, workd, 3*n, info)
      if (info .ne. 0 ) goto 9000
 

 10   continue
         call dnaupd ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        info )
c
         if (ido .eq. -1 ) then
c           (2) <- b . (1)
            call darpax (b, n, workd(ipntr(1)), workd(ipntr(2)))
c           (2) -> (1)
            call dcopy( n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
c           (2) <- a . (1)
            call darpax (a, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c
         if (ido .eq. 1 ) then
c           (2) <- a . (3)
            call darpax (a, n, workd(ipntr(3)), workd(ipntr(2)))
            go to 10
         end if
c 
         if (ido .eq. 2 ) then
c           (2) <- a . (1)
            call darpax (b, n, workd(ipntr(1)), workd(ipntr(2)))
            go to 10
         end if
c 
      if ( info .lt. 0 ) goto 9000
c
      rvec = .true.
c
      call dneupd ( rvec, 'A', iselect, d, d(1,2), v, ldv,
     &      sigma, zero, workev, bmat, n, which, nev, tol,
     &      resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &      lworkl, ierr )
c
      if ( ierr .ne. 0) goto 9000
c
c
 9000 continue
c
      end
c
c ------------------------------------------------------------------
c     matrix vector subroutine
c     Computes ax <--- a*x, 
c
      subroutine darpax(a, n, x, ax)
      integer           n, i, j
      Double precision
     &                  a(n,n), x(n), ax(n),
     &                  xsum
c
c      print *, 'x(1) =', x(1)
      do 1000 i=1,n
        xsum = 0.0;
c        print *, 'x(', i, ') =', x(i)
        do 1001 j=1,n
c          print *, 'a(', i, ',', j ,') =', a(i,j)
          xsum = xsum + a(i,j)*x(j)
1001      continue
        ax(i) = xsum
c        print *, 'a.x(', i, ') =', ax(i)
1000    continue

      return
      end
c

c
c ------------------------------------------------------------------
c     matrix vector subroutine
c     Computes ax <--- a*x, 
c
      subroutine zarpax(a, n, x, ax)
      integer           n, i, j
      Complex*16
     &                  a(n,n), x(n), ax(n),
     &                  xsum
c
c      print *, 'x(1) =', dble(x(1)), dimag(x(1))
      do 1000 i=1,n
        xsum = 0.0;
c        print *, 'x(', i, ') =', dble(x(i)), dimag(x(i))
        do 1001 j=1,n
c          print *, 'a(', i, ',', j ,') =', dble(a(i,j)), dimag(a(i,j))
          xsum = xsum + a(i,j)*x(j)
1001      continue
        ax(i) = xsum
c        print *, 'a.x(', i, ') =', dble(ax(i)), dimag(ax(i))
1000    continue

      return
      end
c
