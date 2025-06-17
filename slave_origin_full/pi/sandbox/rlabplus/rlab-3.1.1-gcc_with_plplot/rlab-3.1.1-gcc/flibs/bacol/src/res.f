      subroutine res(t, y, yprime, delta, ires, rpar, ipar)

c-----------------------------------------------------------------------
c Purpose:
c       This is the subroutine which defines the differential/algebraic
c       system to be solved by DASSL. It returns the residual (delta) of
c       the DAE system:
c                          delta := G(t, Y, Y')
c       The DAE G(t, Y, Y') = 0 arises from applying the method-of-lines
c       and bspline collocation to the system of NPDE PDES of
c       the form:
c                       u_t = f(t, x, u, u_x, u_xx)
c       In the discretized form this yields:
c                      G(t, Y, Y') = A*Y' - F~
c       The abd matrix A contains the collocation equations and some
c       boundary condition information, the vector F~ contains the rhs
c       of the collocation equations and the corresponding boundary
c       conditions.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t
c                               T is the current time.
c
        double precision        y(*)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(*)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        double precision        rpar(*)
c                               rpar is the BACOL floating point work 
c                               array.
c
        integer                 ipar(*)
c                               rpar is the BACOL integer work array.
c
c       Output:
        integer                 ires
c                               ires is a user flag set to alert DASSL
c                               of an illegal y vector. ires is not 
c                               used in this subroutine since there is
c                               no restriction on the vector of 
c                               bspline coefficients.
c
        double precision        delta(*)
c                               delta is the residual of the DAE system.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
c                               ipar(inpde) = npde
c
        integer                 ikcol
c                               ipar(ikcol) = kcol.
c
        integer                 inint
c                               ipar(inint) = nint.
c
        integer                 incpt1
c                               ipar(incpt1) = ncpts1.
c
        integer                 ineq1
c                               ipar(ineq1) = neq1.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ixcol1
c                               rpar(ipar(ixcol1)) stores the 
c                               collocation points when using
c                               dassl_kcol.
c
        integer                 ixcol2
c                               rpar(ipar(ixcol2)) stores the 
c                               collocation points when using
c                               dassl_kcol+1.
c
        integer                 iabbk1
c                               rpar(ipar(iabbk1)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_kcol.
c
        integer                 iabbk2
c                               rpar(ipar(iabbk2)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_kcol+1.
c
        integer                 iwkrj
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi1
c                               rpar(ipar(ibasi1)) stores the basis
c                               function values at the collocation
c                               points when using dassl_kcol.
c                               rpar(ipar(ibasi1)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 ibasi2
c                               rpar(ipar(ibasi2)) stores the basis
c                               function values at the collocation
c                               points when using dassl_kcol+1.
c                               rpar(ipar(ibasi2)) contains
c                               a three dimensional array A of size
c                               (kcol+1+nconti,3,ncpts). A(k,j,i)
c                               contains the values of the (j-1)st
c                               derivative (j=1,2,3) of the k-th
c                               non-zero basis function (k=1,...,
c                               kcol+1+nconti) at the i-th collocation
c                               point.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts1
        integer                 ncpts2
        integer                 neq1
        integer                 neq2
c
c-----------------------------------------------------------------------
c       Direct IPAR indices:
        parameter              (inpde  =  1)
        parameter              (ikcol  =  2)
        parameter              (inint  =  3)
        parameter              (incpt1 =  4)
        parameter              (ineq1  =  5)
c
c-----------------------------------------------------------------------
c       IPAR indices which serve as an indirect pointer into RPAR:
        parameter              (ixcol1 = 22)
        parameter              (iabbk1 = 27)
        parameter              (iwkrj  = 30)
        parameter              (ibasi1 = 31)
c
        parameter              (ixcol2 = 41)
        parameter              (iabbk2 = 46)
        parameter              (ibasi2 = 48)
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                              calres 
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts1 = ipar(incpt1)
      neq1   = ipar(ineq1)
      ncpts2 = ncpts1 + nint
      neq2   = neq1 + nint * npde

c     Calculate residual for dassl_kcol.
      call calres(npde, kcol, nint, ncpts1, neq1, rpar(ipar(ixcol1)),
     &            rpar(ipar(iabbk1)), rpar(ipar(ibasi1)), t, y,
     &            yprime, rpar(ipar(iwkrj)), delta)

c     Calculate residual for dassl_kcol+1.
      call calres(npde, kcol+1, nint, ncpts2, neq2, rpar(ipar(ixcol2)),
     &            rpar(ipar(iabbk2)), rpar(ipar(ibasi2)), t, y(neq1+1),
     &            yprime(neq1+1), rpar(ipar(iwkrj)), delta(neq1+1))

      return
      end
