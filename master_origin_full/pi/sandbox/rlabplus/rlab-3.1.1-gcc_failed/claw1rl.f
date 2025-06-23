c     ===============================================================
      subroutine CLAW1RL(meqn,mwaves,mbc,maux,mwork,mthlim,
     &                   q, work, aux,
     &                   xlower, dx, maxmx,
     %                   t1, t2,
     &                   method, cflv, nv,
     &                   mthbc,
     &                   info, dtv)
c     ===============================================================
c
c     A clawpack driver for RLaB. Based on claw1ez.f by
c       Randall J. LeVeque, distributed with CLAWPACK Version 4.1
c       Version of August, 2002. Original source code is
c       freely downloadable from
c       http://www.amath.washington.edu/~claw/
c
c     This version by Marijan Kostrun, project rlabplus, VIII-2005

      implicit double precision (a-h,o-z)
      external bc1,rp1,src1,b4step1

      dimension    q(1-mbc:maxmx+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, maux)
      dimension work(mwork)
      dimension mthlim(mwaves)
c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(2)

c
      integer info
c
c      dtv(1) = dtmin
c      dtv(2) = dtmax

      info = 0
c
      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) ' require mthbc(1) and mthbc(2) BOTH be set to 2'
         return
         endif

c
c     # grid spacing
c     dx = (xupper - xlower) / float(mx)
c

c
c     # set aux array:
c
c      if (maux .gt. 0)  then
c         call setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
c         endif

c
c     ----------
c     Main loop:
c     ----------
c
         call claw1(maxmx,meqn,mwaves,mbc,maxmx,
     &           q,aux,maux,xlower,dx,t1,t2,dtv,cflv,nv,
     &           method,mthlim,
     &           mthbc,work,mwork,info,bc1,rp1,src1,b4step1)
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) '*** ERROR in claw1 ***  info =',info
            if (info.eq.1) then
               write(6,*) '***   either mx > maxmx or mbc < 2'
               endif
            if (info.eq.2) then
               write(6,*) '***   dt does not divide (tend - tstart)'
               write(6,*) '***   and dt is fixed since method(1)=0'
               endif
            if (info.eq.3) then
               write(6,*) '***   method(1)=1 and cflv(2) > cflv(1)'
               endif
            if (info.eq.4) then
               write(6,*) '***   mwork is too small'
               endif
            if (info.eq.11) then
               write(6,*) '***   Too many time steps, n > nv(1)'
               endif
            if (info.eq.12) then
               write(6,*)
     &          '***   The Courant number is greater than cflv(1)'
               write(6,*) '***   and dt is fixed since method(1)=0'
               endif

            go to 999
            endif

c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c
  999 continue
c
      return
      end


c
c
c     =====================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for 1d problem
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
      dimension amdq(1-mbc:maxmx+mbc, meqn)

      common /comxt/ dtcom,dxcom,tcom

      external claw1r

c     the riemann solver is a c-function which uses a user function within
c     RLaB to solve the riemann problem or part of it.
      call claw1r(tcom,ql,qr,wave,s,amdq,apdq)

      return
      end

