c sparskit.f
c
c Wrapper for sparskit.a library
c
c by Marijan Kostrun, (c) VI 2005

c
c  openfile/closfile for INOUT modules
c
      subroutine OPENFILE(fname,lenfn,idn)
      Integer   lenfn, idn
      Character*256 fname
      If(lenfname .gt. 1) Open (UNIT=idn,FILE=fname,STATUS='OLD')
      Return
      End

      subroutine CLOSFILE(idn)
      Integer   idn
      Close (UNIT=idn)
      Return
      End

c
c  required by iterative solvers
c
      function DISTDOT(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
      end

c
c  runrc: modified to be quiter, and run a solver depending on an parameter isolver
c
      subroutine RUNRC(n,rhs,sol,ipar,fpar,wk,a,ja,ia,
     +     au,jau,ju,isolver)
      implicit none
      integer n,ipar(16),ia(n+1),ja(*),ju(*),jau(*)
      integer isolver
      Double Precision fpar(16),rhs(n),sol(n), wk(*), a(*), au(*)
      external gmres
c-----------------------------------------------------------------------
c     the actual tester. It starts the iterative linear system solvers
c     with a initial guess suppied by the user.
c
c     The structure {au, jau, ju} is assumed to have the output from
c     the ILU* routines in ilut.f.
c
c-----------------------------------------------------------------------
c     local variables
c
      integer i, its
      Double precision res, dnrm2
c     real dtime, dt(2), time
c     external dtime
      external dnrm2
      save its,res
c
c     ipar(2) can be 0, 1, 2, please don't use 3
c
      if (ipar(2).gt.2) ipar(2)=2
c
c     normal execution
c
      its = 0
      res = 0.0D0
c
      ipar(1) = 0
c
10    continue 

      goto(101,102,103,104,105,106,107,108,109) isolver+100

101   call cg(n,rhs,sol,ipar,fpar,wk)
      goto 200

102   call cgnr(n,rhs,sol,ipar,fpar,wk)
      goto 200

103   call bcg(n,rhs,sol,ipar,fpar,wk)
      goto 200

104   call bcgstab(n,rhs,sol,ipar,fpar,wk)
      goto 200

105   call tfqmr(n,rhs,sol,ipar,fpar,wk)
      goto 200

106   call gmres(n,rhs,sol,ipar,fpar,wk)
      goto 200

107   call fgmres(n,rhs,sol,ipar,fpar,wk)
      goto 200

108   call dqgmres(n,rhs,sol,ipar,fpar,wk)
      goto 200

109   call dbcg(n,rhs,sol,ipar,fpar,wk)

200   continue
c
c     output the residuals
c
      res = fpar(5)
c
      if (ipar(1).eq.1) then
         call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.2) then
         call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
         call splusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      endif
c
      return
      end
c-----end-of-runrc
