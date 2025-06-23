      subroutine bacol(t0, tout, atol, rtol, npde, kcol, nintmx, nint,
     &                 x, mflag, rpar, lrp, ipar, lip, y, idid,
     &                 lstdout, stdout
     &                 )

c-----------------------------------------------------------------------
c Purpose:
c       The purpose of BACOL is to solve NPDE dimensional systems of
c       second order parabolic partial differential equations (PDEs)
c       in one space variable of the form:
c
c            dU            /                                \
c            -- (t,x) = f | t, x, U(t,x), U (t,x), U   (t,x) | ,
c            dt            \               x         xx     /
c
c       where x_a < x < x_b and t > t0, with initial conditions at
c       time t = t0 are given by:
c
c                          u(t0,x) = u_0(x),
c
c       for x_a <= x <= x_b, subject to separated boundary conditions
c       given by:
c
c                         /                      \
c                   b    | t, U(t,x_a), U (t,x_a) | = 0,
c                    x_a  \              x       /
c
c                         /                      \
c                   b    | t, U(t,x_b), U (t,x_b) | = 0,
c                    x_b  \              x       /
c
c       for t > t0 and x = x_a, x = x_b, respectively.
c
c       Guide to the above notation:
c          dU
c          -- (t,x) = denotes the first partial derivative of U(t,x)
c          dt         with respect to the time variable t.
c
c          U (t,x) = denotes the first partial derivative of U(t,x)
c           x        with respect to space variable x.
c
c          U  (t,x) = denotes the second partial derivative of U(t,x)
c           xx       with respect to space variable x.
c
c       Further, the above functions NPDE dimensional vector functions.
c
c       BACOL is a method of lines algorithm which uses bspline
c       collocation to discretize the spatial domain [x_a,x_b].
c       The output is a vector of bspline coefficients which
c       can be used to calculate the approximate solution U(t,x) and
c       it's spatial derivatives at (tout,x) where x_a <= x <= x_b
c       and t0 < tout.
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Setup of BACOL:
c       BACOL requires that the user specifies the system of PDEs and
c       the related initial and boundary conditions as well as setting
c       input parameters (which define the bspline space and the
c       requested error tolerances) and allocating work storage.
c
c       The calling sequence of BACOL is:
c
c       call bacol(t0, tout, atol, rtol, npde, kcol, nint, nintmx, x,
c    &             mflag, rpar, lrp, ipar, lip, y, idid)
c
c       which will generate the vector y = Y(tout) upon successful
c       completion. Generally, the call to BACOL will be followed by a
c       call to VALUES to calculate the solution at a set of points:
c
c       call bac1val(kcol, xsol, nint, x, npde, npts, nderiv,
c     &             usol, y, work)
c
c       The details of the parameters to VALUES are documented within
c       the source code for that routine. The input parameters for
c       BACOL are dealt with in detail below, but a quick summary is:
c
c       [t0, tout] is the time domain of the problem.
c       NPDE is the number of compenents in the PDE system.
c       atol is the absolute error tolerance.
c       rtol is the relative error tolerance.
c       kcol, nint, and x define the bspline space.
c       nintmx is the maximum number of subintervals allowed.
c       mflag(1:6) is used to control the operation of BACOL.
c       rpar(lrp) is a floating point work array.
c       ipar(lip) is an integer work array.
c
c       The user must check idid to determine what further action needs
c       to be taken.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               On input, t0 < tout is the initial
c                               time. On output, t0 is the current time,
c                               t0 <= tout.
c
        double precision        tout
c                               tout is the desired final output time.
c
c                               After a successful return from BACOL,
c                               the time stepping may be resumed by
c                               changing tout so that t0 < tout and
c                               setting mflag(1) = 1 to indicate a
c                               continuation of the previous problem.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
c
c                               If PDE components vary in importance,
c                               then vector error tolerances may be
c                               used by setting mflag(2) = 1. In this
c                               case, the dimension of atol must be
c                               npde. The user will define atol(1),
c                               atol(2), ..., atol(npde) appropriately.
c                               Note that a change from scalar to vector
c                               tolerances (or vice versa) constitutes
c                               a new problem, and BACOL will have to
c                               be reinitialized.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
c
c                               If PDE components vary in importance,
c                               then vector error tolerances may be
c                               used by setting mflag(2) = 1. In this
c                               case, the dimension of rtol must be
c                               npde. The user will define rtol(1),
c                               rtol(2), ..., rtol(npde) appropriately.
c                               Note that a change from scalar to vector
c                               tolerances (or vice versa) constitutes
c                               a new problem, and BACOL will have to
c                               be reinitialized.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                               1 < kcol <= mxkcol.
c
c                               The order of the bsplines used
c                               will be (kcol+2).
c
        integer                 nint
c                               at input, nint is the number of
c                               subintervals defined by the spatial
c                               mesh x at the initial time t0.
c                               at output, nint is the number of
c                               subintervals at tout.
c                               nint >= 1.
c
        integer                 nintmx
c                               the maximum number of subintervals that
c                               the user requires.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a,x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c                               at input, x(1:nint+1) stores the mesh
c                               points at the initial time t0.
c                               at output, x(1:nint+1) stores the mesh
c                               points at tout.
c
        integer                 mflag(6)
c                               This vector determines the interaction
c                               of BACOL with DASSL.
c
c                                         How to set mflag(1):
c                               On the initial call to BACOL with a
c                               new problem, set mflag(1) = 0, which
c                               indicates that BACOL and DASSL should
c                               perform the initialization steps that
c                               are required by each code, respectively.
c
c                               In order to continue time stepping in
c                               the current problem after a successful
c                               return from BACOL, set mflag(1) = 1,
c                               idid = 1, and ensure that t0 < tout.
c
c                                         How to set mflag(2):
c                               If scalar absolute and relative error
c                               tolerances (atol and rtol) are desired,
c                               then set mflag(2) = 0.
c
c                               For vector absolute and relative error
c                               tolerances, set mflag(2) = 1, define
c                               atol(1), ..., atol(npde), and
c                               rtol(1), ..., rtol(npde), as described
c                               above, ensuring that the dimension of
c                               each of atol and rtol is at least npde.
c
c                                         How to set mflag(3):
c                               If there are no restrictions on t0,
c                               then set mflag(3) = 0.
c                               Since DASSL may actually "step over"
c                               tout and then interpolate, there is the
c                               option to enforce t0 <= tstop <= tout.
c                               If this is desirable, set mflag(3) = 1,
c                               and define rpar(1) = tstop.
c
c                                         How to set mflag(4):
c                               If the user wishes, BACOL will return
c                               the computed solution and derivative
c                               after a certain number of accepted time
c                               steps or at TOUT, whichever comes first.
c                               This is a good way to proceed if the
c                               user wants to see the behavior of the
c                               solution.
c                               If the user only wants the solution at
c                               TOUT, set mflag(4) = 0;
c                               else, set mflag(4) = 1, and assign a
c                               positive integer for ipar(8) that
c                               defines the number of time steps before
c                               BACOL is stopped.
c
c                                         How to set mflag(5):
c                               If both boundary conditions are
c                               dirichlet, set mflag(5) = 1;
c                                    else, set mflag(5) = 0.
c
c                                         How to set mflag(6):
c                               If the user wants to specify an initial
c                               stepsize, set mflag(6) = 1,
c                                             and define rpar(2) = the
c                                             initial stepsize;
c                                   else, set mflag(6) = 0;
c
        integer                 lrp
c                               lrp is the size of the rpar storage
c                               array and must satisfy:
c                               lrp>=134+nintmx*(35+35*kcol+31*npde
c    +                               +38*npde*kcol+8*kcol*kcol)+14*kcol
c    +                               +79*npde+npde*npde*(21
c    +                               +4*nintmx*kcol*kcol+12*nintmx*kcol
c    +                               +6*nintmx)
c
        integer                 lip
c                               lip is the size of the ipar integer
c                               work array and must satisfy:
c                               lip>=115+npde*(nintmx*(2*kcol+1)+4)
c
c       Work Storage:
        double precision        rpar(lrp)
c                               rpar is a floating point work array
c                               of size lrp.
c
        integer                 ipar(lip)
c                               ipar is an integer work array
c                               of size lip.
c
c       Output:
        double precision        y(npde*(kcol*nintmx+2))
c                               On successful return from BACOL,
c                               y(1:npde*(kcol*nint+2)) is
c                               the vector of bspline coefficients at
c                               the current time.
c
        integer                 idid
c                               idid is the BACOL exit status flag
c                               which is based on the exit status from
c                               DASSL plus some additional status codes
c                               based on error checking performed by
c                               BACOL on initialization. Positive
c                               values of idid indicate a successful
c                               return. Negative values of idid indicate
c                               an error which may or may not be fatal.
c                               The exact descriptions of idid return
c                               values will be discussed below.
c
c                               For calls other than the first call
c                               (mflag(1) = 1), the user must check
c                               idid, set idid = 1 (if necessary), and
c                               take other actions which are necessary
c                               such as defining a new value of tout.
c
c                               An excerpt from the DASSL source code
c                               documentation is included to define
c                               the IDID return codes and clarify
c                               the operation of DASSL within BACOL.
c
c-----------------------------------------------------------------------
c  The following is an excerpt from the DASSL source code documentation:
c
C  -------- OUTPUT -- AFTER ANY RETURN FROM DASSL ---------------------
C
C  The principal aim of the code is to return a computed solution at
C  TOUT, although it is also possible to obtain intermediate results
C  along the way. To find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the IDID parameter.
C
C  IDID -- Reports what the code did.
C
C                     *** Task completed ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- A step was successfully taken in the
C                   intermediate-output mode. The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- The integration to TSTOP was successfully
C                   completed (T=TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- The integration to TOUT was successfully
C                   completed (T=TOUT) by stepping past TOUT.
C                   Y(*) is obtained by interpolation.
C                   YPRIME(*) is obtained by interpolation.
C
C                    *** Task interrupted ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- A large amount of work has been expended.
C                   (About 500 steps)
C
C           IDID = -2 -- The error tolerances are too stringent.
C
C           IDID = -3 -- The local error test cannot be satisfied
C                   because you specified a zero component in ATOL
C                   and the corresponding computed solution
C                   component is zero. Thus, a pure relative error
C                   test is impossible for this component.
C
C           IDID = -6 -- DASSL had repeated error test
C                   failures on the last attempted step.
C
C           IDID = -7 -- The corrector could not converge.
C
C           IDID = -8 -- The matrix of partial derivatives
C                   is singular.
C
C           IDID = -9 -- The corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           IDID =-10 -- The corrector could not converge
C                   because IRES was equal to minus one.
C
C           IDID =-11 -- IRES equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           IDID =-12 -- DASSL failed to compute the initial
C                   YPRIME.
C
C                    *** Task terminated ***
C                Reported by the value of IDID=-33
C
C           IDID = -33 -- The code has encountered trouble from which
C                   it cannot recover. A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. For example, this occurs
C                   when invalid input is detected.
C
C  RTOL, ATOL -- These quantities remain unchanged except when
C               IDID = -2. In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration. However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
C                    (CALLS AFTER THE FIRST)
C
C  This code is organized so that subsequent calls to continue the
C  integration involve little (if any) additional effort on your
C  part. You must monitor the IDID parameter in order to determine
C  what to do next.
C
C  Recalling that the principal task of the code is to integrate
C  from T to TOUT (the interval mode), usually all you will need
C  to do is specify a new TOUT upon reaching the current TOUT.
C
C  Do not alter any quantity not specifically permitted below,
C  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
C  or the differential equation in subroutine RES. Any such
C  alteration constitutes a new problem and must be treated as such,
C  i.e., you must start afresh.
C
C  You cannot change from vector to scalar error control or vice
C  versa (mflag(2)), but you can change the size of the entries of
C  RTOL, ATOL. Increasing a tolerance makes the equation easier
C  to integrate. Decreasing a tolerance will make the equation
C  harder to integrate and should generally be avoided.
C
C  If it has been necessary to prevent the integration from going
C  past a point TSTOP (mflag(3), rpar(itstop)), keep in mind that
C  the code will not integrate to any TOUT beyond the currently
C  specified TSTOP. Once TSTOP has been reached you must change
C  the value of TSTOP or set mflag(3)=0. You may change mflag(3)
C  or TSTOP at any time but you must supply the value of TSTOP in
C  rpar(itstop) whenever you set mflag(3)=1.
C
C  -------- ERROR MESSAGES ---------------------------------------------
C
C      The SLATEC error print routine XERMSG is called in the event of
C   unsuccessful completion of a task.  Most of these are treated as
C   "recoverable errors", which means that (unless the user has directed
C   otherwise) control will be returned to the calling program for
C   possible action after the message has been printed.
C
C   In the event of a negative value of IDID other than -33, an appro-
C   priate message is printed and the "error number" printed by XERMSG
C   is the value of IDID.  There are quite a number of illegal input
C   errors that can lead to a returned value IDID=-33.  The conditions
C   and their printed "error numbers" are as follows:
C
C   Error number       Condition
C
C        1       Some element of INFO vector is not zero or one.
C        2       NEQ .le. 0
C        3       MAXORD not in range.
C        4       LRW is less than the required length for RWORK.
C        5       LIW is less than the required length for IWORK.
C        6       Some element of RTOL is .lt. 0
C        7       Some element of ATOL is .lt. 0
C        8       All elements of RTOL and ATOL are zero.
C        9       INFO(4)=1 and TSTOP is behind TOUT.
C       10       HMAX .lt. 0.0
C       11       TOUT is behind T.
C       12       INFO(8)=1 and H0=0.0
C       13       Some element of WT is .le. 0.0
C       14       TOUT is too close to T to start integration.
C       15       INFO(4)=1 and TSTOP is behind T.
C       16       --( Not used in this version )--
C       17       ML illegal.  Either .lt. 0 or .gt. NEQ
C       18       MU illegal.  Either .lt. 0 or .gt. NEQ
C       19       TOUT = T.
c----------------------------------------------------------------------
C
C   If DASSL is called again without any action taken to remove the
C   cause of an unsuccessful return, XERMSG will be called with a fatal
C   error flag, which will cause unconditional termination of the
C   program.  There are two such fatal errors:
C
C   Error number -998:  The last step was terminated with a negative
C       value of IDID other than -33, and no appropriate action was
C       taken.
C
C   Error number -999:  The previous call was terminated because of
C       illegal input (IDID=-33) and there is illegal input in the
C       present call, as well.  (Suspect infinite loop.)
C

c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 MAXORD
        parameter              (MAXORD = 5)
c                               MAXORD is the maximum order of the
c                               backward differentiation formula (BDF)
c                               methods used by DASSL. MAXORD = 5 is
c                               the default used by DASSL.
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
        integer                 maxrsh
        parameter              (maxrsh = 20)
c                               maxrsh is the maximum number of
c                               remesh times at one time step,
c                               i.e., icount must less than or equal
c                               to maxrsh
c
        double precision        one
        parameter              (one = 1.0D0)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Local variables:
c
        integer                 neq1
c                               neq1=npde*ncpts1 is the number of
c                               bspline coefficients (or DAEs) when
c                               using dassl_{kcol}.
c
        integer                 neq2
c                               neq2=neq1+npde*nint is the number of
c                               bspline coefficients (or DAEs) when
c                               using dassl_{kcol+1}.
c
        integer                 neq
c                               neq = neq1 + neq2.
c
        integer                 leniw
c                               leniw = 20 + neq is the length of the
c                               integer work array required by dassl.
c
        integer                 lenpd1
c                               lenpd1 is the size of the Almost Block
c                               Diagonal (ABD) Jacobian required by
c                               dassl_{kcol}.
c                               lenpd1=npde*npde*(2*nconti
c                                      +kcol*(kcol+nconti)*nint)
c
        integer                 lenpd2
c                               lenpd2 is the size of the Almost Block
c                               Diagonal (ABD) Jacobian required by
c                               dassl_{kcol+1}.
c                               lenpd2=lenpd1+npde*npde*nint
c                                      *(2*kcol+nconti+1)
c
        integer                 lenpd
c                               lenpd = lenpd1 + lenpd2 .
c
        integer                 lenrw
c                               lenrw = 40+(MAXORD+4)*neq+lenpd
c                               is the total size of the floating point
c                               work array required by dassl_{kcol}.
c
        integer                 lenin1
c                               lenin1 is the size of the floating
c                               point work array used by INIY and INIYP
c                               when using dassl_{kcol}.
c                               lenin1>=lenpd1+2*neq1+npde*2+2*npde*npde
c
        integer                 lenin2
c                               lenin2 is the size of the floating
c                               point work array used by INIY and INIYP
c                               when using dassl_{kcol+1}.
c                               lenin2>=lenpd2+2*neq2+npde*2+2*npde*npde
c
        integer                 lenri1
c                               lenri1 is the size of the floating
c                               point work array used by REINIT when
c                               using dassl_{kcol}.
c
        integer                 lenri2
c                               lenri2 is the size of the floating
c                               point work array used by REINIT when
c                               using dassl_{kcol+1}.
c
        integer                 lenrj
c                               lenrj is the size of the floating
c                               point work array used by RES and JAC.
c                               lenrj>=4*npde+5*npde*npde.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c                               lenerr>=2*npde*necpts+npde*nint.
c
        integer                 ncpts1
c                               ncpts1=(kcol*nint+nconti) is the total
c                               number of collocation points when using
c                               dassl_{kcol}.
c
        integer                 ncpts2
c                               ncpts2=ncpts1+nint is the total number
c                               of collocation points when using
c                               dassl_{kcol+1}.
c
        integer                 necpts
c                               necpts=(kcol+3)*nint is the total number
c                               of collocation points used for
c                               error estimate.
c
        integer                 icflag
c                               This is the status flag from the almost
c                               block diagnonal factorization routine,
c                               CRDCMP.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        torign
c                               torign is the initial time, i.e. = t0
c                               at the beginning.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 isstep
c                               isstep is the number of accepted time
c                               steps since we restart BACOL in the
c                               case that mflag(4) = 1.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before the current remeshing.
c
        integer                 ninpre
c                               ninpre is the number of subintervals
c                               when icount = 0 before remeshing.
c
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial step or any step
c                                           without needing remesh;
c                                      = 1, remesh with a hot start.
c                                      = 2, remesh with a cold start.
c
        integer                 neqpre
c                               neqpre is the number of bspline
c                               coefficients when icount = 0 before
c                               remeshing when using dassl_{kcol+1}.
c
        integer                 irold
c                               irold is the value of ipar(ixold) before
c                               remeshing.
c
        integer                 nstep
c                               nstep is the number of steps which are
c                               necessary for a remesh.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, ii, jj, kk
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 itstop
        parameter              (itstop =  1)
c                               rpar(itstop) = tstop as defined when
c                               mflag(3) = 1.
c
        integer                 iiniss
        parameter              (iiniss =  2)
c                               rpar(iiniss) = the initial stepsize when
c                               mflag(6) = 1.
c
        integer                 irpstr
        parameter              (irpstr = 11)
c                               rpar(1:irpstr-1) are reserved to store
c                               floating point scalar quantities.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt1
        parameter              (incpt1 =  4)
c                               ipar(incpt1) = ncpts1.
c
        integer                 ineq1
        parameter              (ineq1  =  5)
c                               ipar(ineq1) = neq1.
c
        integer                 iipstp
        parameter              (iipstp =  6)
c                               ipar(iipstp) = the minimum size of ipar.
c
        integer                 irpstp
        parameter              (irpstp =  7)
c                               ipar(irpstp) = the minimum size of rpar.
c
        integer                 iinstp
        parameter              (iinstp =  8)
c                               As input, ipar(iinstp) is the number of
c                               intermediate time steps before the user
c                               asks BACOL to return the computed
c                               solution.
c                               As output, ipar(iinstp) will keep the
c                               same value unless TOUT is reached before
c                               ipar(iinstp) time steps. If TOUT is
c                               reached, ipar(iinstp) will be the number
c                               of time steps before TOUT is reached.
c                               See mflag(4) for more details.
c
        integer                 irshin
        parameter              (irshin =  9)
c                               ipar(irshin) is the number of remeshing
c                               times at the initial step.
c
        integer                 isteps
        parameter              (isteps = 10)
c                               ipar(isteps) is the number of time steps
c                               on the current problem.
c
        integer                 irmesh
        parameter              (irmesh = 11)
c                               ipar(irmesh) is the number of remeshing
c                               times after BACOL starts the initial
c                               step.
c
        integer                 istalr
        parameter              (istalr = 12)
c                               ipar(istalr) is the number of accepted
c                               steps after the last successful
c                               remeshing.
c
        integer                 istblc
        parameter              (istblc = 13)
c                               ipar(istblc) is the number of steps
c                               BACOL has taken before the latest cold
c                               start.
c
        integer                 icolds
        parameter              (icolds = 14)
c                               ipar(icolds) is the number of the times
c                               when BACOL performs a cold start after
c                               remeshing.
c
        integer                 idasi
        parameter              (idasi  = 61)
c                               ipar(idasi) stores, before remeshing,
c                               the first 20 elements of the integer
c                               point work array in dassl.
c
        integer                 iinfo
        parameter              (iinfo  = 81)
c                               ipar(iinfo) is an integer array
c                               required by dassl_{kcol}.
c                               ipar(iinfo)   = mflag(1),
c                               ipar(iinfo+1) = mflag(2),
c                               ipar(iinfo+3) = mflag(3),
c                               ipar(iinfo+7) = mflag(6).
c                               (See the documentation of DDASSL for
c                               details.)
c
        integer                 iiwork
        parameter              (iiwork = 96)
c                               ipar(iiwork) is the integer work array
c                               for dassl.
c
        integer                 ipivot
        parameter              (ipivot = 116)
c                               ipar(ipivot-1+i), i = 1, neq, contains
c                               the pivoting information from the
c                               factorization of the temporary matrix
c                               for dassl.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ih
        parameter              (ih     = 21)
c                               rpar(ipar(ih)) stores the mesh step
c                               size sequence.
c
        integer                 ixcol1
        parameter              (ixcol1 = 22)
c                               rpar(ipar(ixcol1)) stores the
c                               collocation points when using
c                               dassl_{kcol}.
c
        integer                 ixbs1
        parameter              (ixbs1  = 23)
c                               rpar(ipar(ixbs1)) stores the breakpoint
c                               sequence when using dassl_{kcol}.
c
        integer                 iy1
        parameter              (iy1    = 24)
c                               rpar(ipar(iy1)) stores the vector of
c                               solution components to the DAE system
c                               when using dassl_{kcol}.
c
        integer                 iyp1
        parameter              (iyp1   = 25)
c                               rpar(ipar(iyp1)) stores the vector of
c                               solution component derivatives of the
c                               DAE system when using dassl_{kcol}.
c
        integer                 iabtp1
        parameter              (iabtp1 = 26)
c                               rpar(ipar(iabtp1)) stores the top block
c                               of the ABD collocation matrices when
c                               using dassl_{kcol}.
c
        integer                 iabbk1
        parameter              (iabbk1 = 27)
c                               rpar(ipar(iabbk1)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_{kcol}.
c
        integer                 iabbt1
        parameter              (iabbt1 = 28)
c                               rpar(ipar(iabbt1)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using dassl_{kcol}.
c
        integer                 irwork
        parameter              (irwork = 29)
c                               rpar(ipar(irwork)) stores the floating
c                               point work array for DASSL. And it
c                               is also used to be a work storage for
c                               the subroutine INIY and INIYP to get
c                               the initial condition.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by RES and JAC.
c
        integer                 ibasi1
        parameter              (ibasi1 = 31)
c                               rpar(ipar(ibasi1)) stores the basis
c                               function values at the collocation
c                               points when using dassl_{kcol}.
c                               rpar(ipar(ibasi1)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts1). A(k,j,i) stores
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 itatol
        parameter              (itatol = 32)
c                               rpar(ipar(itatol)) is the absolute error
c                               tolerance request on the integration
c                               error.
c
        integer                 itrtol
        parameter              (itrtol = 33)
c                               rpar(ipar(itrtol)) is the relative error
c                               tolerance request on the integration
c                               error.
c
        integer                 iexcol
        parameter              (iexcol = 34)
c                               rpar(ipar(iexcol)) stores the
c                               collocation points which are used for
c                               error estimate.
c
        integer                 iewts
        parameter              (iewts  = 35)
c                               rpar(ipar(iewts)) stores the gaussian
c                               weights which are used for error
c                               estimate.
c
        integer                 iebas1
        parameter              (iebas1 = 36)
c                               rpar(ipar(iebas1)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol}.
c
        integer                 iebas2
        parameter              (iebas2 = 37)
c                               rpar(ipar(iebas2)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol+1}.
c
        integer                 iercom
        parameter              (iercom = 38)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 39)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 40)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 ixcol2
        parameter              (ixcol2 = 41)
c                               rpar(ipar(ixcol2)) stores the
c                               collocation points when using
c                               dassl_{kcol+1}.
c
        integer                 ixbs2
        parameter              (ixbs2  = 42)
c                               rpar(ipar(ixbs2)) stores the breakpoint
c                               sequence when using dassl_{kcol+1}.
c
        integer                 iy2
        parameter              (iy2    = 43)
c                               rpar(ipar(iy2)) stores the vector of
c                               solution components to the DAE system
c                               when using dassl_{kcol+1}.
c
        integer                 iyp2
        parameter              (iyp2   = 44)
c                               rpar(ipar(iyp2)) stores the vector of
c                               solution component derivatives of the
c                               DAE system when using dassl_{kcol+1}.
c
        integer                 iabtp2
        parameter              (iabtp2 = 45)
c                               rpar(ipar(iabtp2)) stores the top block
c                               of the ABD collocation matrices when
c                               using dassl_{kcol+1}.
c
        integer                 iabbk2
        parameter              (iabbk2 = 46)
c                               rpar(ipar(iabbk2)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_{kcol+1}.
c
        integer                 iabbt2
        parameter              (iabbt2 = 47)
c                               rpar(ipar(iabbt2)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using dassl_{kcol+1}.
c
        integer                 ibasi2
        parameter              (ibasi2 = 48)
c                               rpar(ipar(ibasi2)) stores the basis
c                               function values at the collocation
c                               points when using dassl_{kcol+1}.
c                               rpar(ipar(ibasi2)) contains
c                               a three dimensional array A of size
c                               (kcol+1+nconti,3,ncpts2). A(k,j,i)
c                               stores the values of the (j-1)st
c                               derivative (j=1,2,3) of the k-th
c                               non-zero basis function (k=1,...,
c                               kcol+1+nconti) at the i-th collocation
c                               point.
c
        integer                 iwkdnm
        parameter              (iwkdnm = 49)
c                               rpar(ipar(iwkdnm)) is the work storage
c                               for the modification version of the
c                               subroutine DDANRM.
c
        integer                 ixold
        parameter              (ixold  = 51)
c                               rpar(ipar(ixold)) stores the mesh point
c                               sequence when icount = 0 before
c                               remeshing.
c
        integer                 idasr
        parameter              (idasr  = 52)
c                               rpar(ipar(idasr)) stores, before
c                               remeshing, the first 40 elements of the
c                               floating point work array in dassl.
c
        integer                 iypre
        parameter              (iypre  = 53)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy2)) at the previous 6 steps.
c                               It is required for a hot restart after
c                               remeshing.
c
        integer                 iyprer
        parameter              (iyprer = 54)
c                               rpar(ipar(iyprer)) stores the
c                               information of at the previous steps
c                               after remeshing.
c
        integer                 iey1
        parameter              (iey1   = 55)
c                               rpar(ipar(iey1)) stores the bspline
c                               coefficients for dassl_{kcol} at the
c                               farthest point integration has reached.
c
        integer                 iey2
        parameter              (iey2   = 56)
c                               rpar(ipar(iey2)) stores the bspline
c                               coefficients for dassl_{kcol+1} at the
c                               farthest point integration has reached.

c Marijan Kostrun, project rlabplus
      save

      integer       lstdout
      character*256 stdout
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bcolpnt
c                               ddassl
c                               divdif
c                               errest
c                               iniy
c                               iniyp
c                               meshsq
c                               reinit
c                               remesh
c                               sucstp
      external                jac
      external                res

      if (lstdout .gt. 1) then
        lout = 6
        open (lout, file=stdout, status='old' )
      else
        lout = 0
      endif


c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, January 2, 2001.
c
c-----------------------------------------------------------------------

c     Check validity of the mflag vector.
      do 1 i = 1, 6
         if ((mflag(i) .lt. 0) .or. (mflag(i) .gt. 1)) goto 710
   1  continue

      if (mflag(4) .eq. 1) then
         isstep = 0
         if (ipar(iinstp) .le. 0) goto 715
      endif

      irshfg = 0

c     Check for continuation of a previous problem.
      if (mflag(1) .eq. 1) then
         istart = 1
         goto 200
      else

c        Check if the user specifies an initial stepsize
         if (mflag(6) .eq. 1) then
            if (((tout-t0)*rpar(iiniss)) .lt. zero) goto 720
            if (rpar(iiniss) .eq. zero) goto 725
         endif

         ipar(irmesh) = 0
         ipar(irshin) = 0
         ipar(istalr) = 0
         ipar(istblc) = 0
         ipar(icolds) = 0
         istart = 0
         irshfg = 0
         torign = t0
      endif

c-----------------------------------------------------------------------
c     On the initial call or after remeshing, check for valid input and
c     initialize the workspace.
c-----------------------------------------------------------------------

  100 continue

c     Check validity of npde, kcol, and nint.
      if (npde .le. 0) goto 730
      if ((kcol .le. 1) .or. (kcol .gt. mxkcol)) goto 740
      if ((nint .le. 0) .or. (nint .gt. nintmx)) goto 750

c     Check for a monotone mesh.
      do 110 i = 1, nint
         if (x(i) .ge. x(i+1)) goto 760
  110 continue

c-----------------------------------------------------------------------
c     Calculate the extra storage requirements of res and jac.
      lenrj = (4 + 5 * npde) * npde

c     Calculate the number of collocation points when using
c     dassl_{kcol}.
      ncpts1 = nint * kcol + nconti

c     Calculate the number of DAEs when using dassl_{kcol}.
      neq1 = npde * ncpts1

c     Size of the ABD iteration matrix when using dassl_{kcol}.
      lenpd1 = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)

c     Calculate the extra storage requirements of iniy and iniyp when
c     using dassl_{kcol}.
      lenin1 = lenpd1 + 2 * neq1 + 2 * npde * (1 + npde)

c-----------------------------------------------------------------------
c     Calculate the number of collocation points when using
c     dassl_{kcol+1}.
      ncpts2 = ncpts1 + nint

c     Calculate the number of DAEs when using dassl_{kcol+1}.
      neq2 = neq1 + npde * nint

c     Size of the ABD iteration matrix when using dassl_{kcol+1}.
      lenpd2 = lenpd1 + npde * npde * nint * (kcol + kcol + nconti + 1)

c     Calculate the extra storage requirements of iniy and iniyp when
c     using dassl_{kcol+1}.
      lenin2 = lenpd2 + 2 * neq2 + 2 * npde * (1 + npde)

c-----------------------------------------------------------------------
c     Calculate the total number of variables given to dassl.
      neq = neq1 + neq2

c     Size of the total iteration matrix in dassl.
      lenpd = lenpd1 + lenpd2

c     Total size of the DASSL floating point work array.
      lenrw = 40 + (MAXORD + 4) * neq + lenpd

c     Total size of the DASSL integer work array.
      leniw = 20 + neq

c-----------------------------------------------------------------------
c     Calculate the number of collocation point used for error estimate.
      necpts = (kcol + 3) * nint

c     Calculate the extra storage requirements of errest.
      lenerr = (2 * necpts + nint) * npde

c-----------------------------------------------------------------------
c     Save the input parameters in the ipar integer communication
c     storage array.
      ipar(inpde)  = npde
      ipar(ikcol)  = kcol
      ipar(inint)  = nint
      ipar(incpt1) = ncpts1
      ipar(ineq1)  = neq1

c-----------------------------------------------------------------------
c     Calculate the offsets into the rpar floating point storage array.
c-----------------------------------------------------------------------
      ipar(itatol) = irpstr
      ipar(itrtol) = ipar(itatol) + npde
      ipar(ih)     = ipar(itrtol) + npde

      ipar(iy1)    = ipar(ih)     + nint
      ipar(iy2)    = ipar(iy1)    + neq1
      ipar(iyp1)   = ipar(iy2)    + neq2
      ipar(iyp2)   = ipar(iyp1)   + neq1

      ipar(ixcol1) = ipar(iyp2)   + neq2
      ipar(ixbs1)  = ipar(ixcol1) + ncpts1
      ipar(iabtp1) = ipar(ixbs1)  + ncpts1 + kcol + nconti
      ipar(iabbk1) = ipar(iabtp1) + npde * npde * nconti
      ipar(iabbt1) = ipar(iabbk1) + npde * npde * nint * kcol
     &                              * (kcol + nconti)
      ipar(ibasi1) = ipar(iabbt1) + npde * npde * nconti

      ipar(irwork) = ipar(ibasi1) + (kcol + nconti) * 3 * ncpts1
      ipar(iwkrj)  = ipar(irwork) + lenrw

      ipar(iexcol) = ipar(iwkrj)  + lenrj
      ipar(iewts)  = ipar(iexcol) + necpts
      ipar(ierint) = ipar(iewts)  + necpts
      ipar(iercom) = ipar(ierint) + nint
      ipar(iebas1) = ipar(iercom) + npde
      ipar(iebas2) = ipar(iebas1) + (kcol + nconti) * necpts
      ipar(iework) = ipar(iebas2) + (kcol + 1 + nconti) * necpts

      ipar(ixcol2) = ipar(iework) + lenerr
      ipar(ixbs2)  = ipar(ixcol2) + ncpts2
      ipar(iabtp2) = ipar(ixbs2)  + ncpts2 + kcol + 1 + nconti
      ipar(iabbk2) = ipar(iabtp2) + npde * npde * nconti
      ipar(iabbt2) = ipar(iabbk2) + npde * npde * nint * (kcol+1)
     &                              * (kcol + 1 + nconti)
      ipar(ibasi2) = ipar(iabbt2) + npde * npde * nconti

      ipar(iwkdnm) = ipar(ibasi2) + (kcol + 1 + nconti) * 3 * ncpts2

      ipar(iyprer) = ipar(iwkdnm) + neq
      ipar(ixold)  = ipar(iyprer) + 6 * neq2
      ipar(idasr)  = ipar(ixold)  + nintmx + 1
      ipar(iypre)  = ipar(idasr)  + 40

      ipar(iey1)   = ipar(irwork) + 40 + 3 * neq
      ipar(iey2)   = ipar(iey1) + neq1

c     The offset is different between the initial call and remeshing.
      if ((irshfg .ne. 0) .and. (istart .eq. 1)) then
         ipar(irpstp) = ipar(iypre) + 6 * neqpre - 1
      else
         ipar(irpstp) = ipar(iypre) + 6 * neq2 - 1
      endif

c     Check for a sufficiently large rpar floating point work array.
      if (lrp .lt. ipar(irpstp)) goto 770

c     Calculate the offsets into the integer storage array.
      ipar(iipstp) = ipivot + neq - 1

c     Check for a sufficiently large ipar integer work array.
      if (lip .lt. ipar(iipstp)) goto 780

c     Check whether it is initial call or for remeshing.
      if ((irshfg .ne. 0) .and. (istart .ne. 0)) goto 300

c-----------------------------------------------------------------------
c     Perform initializations for using dassl_{kcol}.
c-----------------------------------------------------------------------
c     Check whether it is initial call or for remeshing.
      if (irshfg .eq. 0) then

c        Set the info vector required by DASSL.
         ipar(iinfo)   = mflag(1)
         ipar(iinfo+1) = mflag(2)
         ipar(iinfo+2) = 1
         ipar(iinfo+3) = mflag(3)

         do 120 i = 6, 14
            ipar(iinfo+i-1) = 0
  120    continue

c        Indicate whether an user-supplied initial stepsize is used.
         ipar(iinfo+7) = mflag(6)

c        Indicate an analytic user supplied Jacobian (iteration) matrix.
         ipar(iinfo+4) = 1

c        Indicate an ABD Jacobian (Iteration) matrix.
         ipar(iinfo+14) = 1

      else
         ipar(irshin) = ipar(irshin) + 1
      endif

!                 write (*,*) '1'
      call meshsq(kcol, nint, x, rpar(ipar(irwork)), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))
!                 write (*,*) '2'

      call bcolpnt(kcol, nint, ncpts1, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol1)),
     &            rpar(ipar(ixbs1)))
!                 write (*,*) '3'

      icflag = 0

      call iniy(t0, npde, kcol, nint, neq1, ncpts1, mflag(5),
     &          rpar(ipar(ixcol1)), rpar(ipar(ixbs1)),
     &          rpar(ipar(iabbk1)), rpar(ipar(ibasi1)),
     &          rpar(ipar(iy1)), ipar(ipivot),
     &          rpar(ipar(irwork)), lenin1, icflag)
!                 write (*,*) '4'

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

!                 write (*,*) '5'
      call iniyp(t0, npde, kcol, nint, neq1, ncpts1,
     &           rpar(ipar(ixcol1)), rpar(ipar(iabtp1)),
     &           rpar(ipar(iabbk1)), rpar(ipar(iabbt1)),
     &           rpar(ipar(ibasi1)), rpar(ipar(iy1)),
     &           rpar(ipar(iyp1)), ipar(ipivot),
     &           rpar(ipar(irwork)), lenin1, icflag)
!                 write (*,*) '6'
      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

      irshfg = 0

c     Set the parameters defining the ABD Jacobian (iteration) matrix.
c     This is inside the integer work array of dassl. Those parameter
c     include npde, kcol, nint.
      ipar(iiwork-1+17) = npde
      ipar(iiwork-1+18) = kcol
      ipar(iiwork-1+19) = nint

c     Set initial idid to be zero.
      idid = 0

c-----------------------------------------------------------------------
c     Perform initializations for using dassl_{kcol+1}.
c-----------------------------------------------------------------------

c     Copy rtol and atol to be the relative and absolute error request
c     for the integration error.
      if (irshfg .eq. 0) then
         do 130 i = 1, npde
            rpar(ipar(itatol)-1+i) = atol(i)
  130    continue
         do 140 i = 1, npde
            rpar(ipar(itrtol)-1+i) = rtol(i)
  140    continue
      endif
!                 write (*,*) '7'
      call bcolpnt(kcol+1, nint, ncpts2, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol2)),
     &            rpar(ipar(ixbs2)))

      icflag = 0
!                 write (*,*) '8'
      call iniy(t0, npde, kcol+1, nint, neq2, ncpts2, mflag(5),
     &          rpar(ipar(ixcol2)), rpar(ipar(ixbs2)),
     &          rpar(ipar(iabbk2)), rpar(ipar(ibasi2)),
     &          rpar(ipar(iy2)), ipar(ipivot),
     &          rpar(ipar(irwork)), lenin2, icflag)

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif
!                 write (*,*) '9'
      call iniyp(t0, npde, kcol+1, nint, neq2, ncpts2,
     &           rpar(ipar(ixcol2)), rpar(ipar(iabtp2)),
     &           rpar(ipar(iabbk2)), rpar(ipar(iabbt2)),
     &           rpar(ipar(ibasi2)), rpar(ipar(iy2)),
     &           rpar(ipar(iyp2)), ipar(ipivot),
     &           rpar(ipar(irwork)), lenin2, icflag)

!                 write (*,*) '(bacol) icflag = ', icflag

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif
!                 write (*,*) '10'
c     Copy rpar(ipar(iy2)) to rpar(ipar(iypre)).
      call dcopy(neq2, rpar(ipar(iy2)), 1, rpar(ipar(iypre)), 1)

      if (mflag(3) .eq. 1) then
         rpar(ipar(irwork)) = rpar(itstop)
      endif

c     Set the initial stepsize if applicable.
      if (mflag(6) .eq. 1) rpar(ipar(irwork)+2) = rpar(iiniss)

      goto 400

c-----------------------------------------------------------------------
c     When an adaptive mesh is used, this is not the first call for
c     the problem, and integration is to continue.
c-----------------------------------------------------------------------
  200 continue

c     Examine idid to determine if DASSL can be called again.
      if (idid .ne. 1) goto 790

c     If t0 must not go beyond tstop in the time stepping, then update
c     the DASSL floating point work array with the value of tstop.
      if (mflag(3) .eq. 1) then
         rpar(ipar(irwork)) = rpar(itstop)
         rpar(ipar(idasr))  = rpar(itstop)
      endif

c     Reset t0. If necessary, update rpar(ipar(iy1)), rpar(ipar(iyp1)),
c     rpar(ipar(iy2)) and rpar(ipar(iyp2)).
      if (rpar(ipar(irwork)-1+4) .lt. tout) then
         t0 = rpar(ipar(irwork)-1+4)
         goto 400
      else
         t0 = tout
         call ddatrp(rpar(ipar(irwork)-1+4), tout, rpar(ipar(iy1)),
     &               rpar(ipar(iyp1)), neq, ipar(iiwork-1+8),
     &               rpar(ipar(irwork)+40+3*neq),
     &               rpar(ipar(irwork)-1+29))
         idid = 3
         goto 500
      endif

c-----------------------------------------------------------------------
c     Initialization after remeshing.
c-----------------------------------------------------------------------
  300 continue

      ipar(irmesh) = ipar(irmesh) + 1

      do 310 i = 1, 20
         ipar(iiwork-1+i) = ipar(idasi-1+i)
  310 continue

c     Tell DASSL to calculate the new iteration matrix and reset the
c     size of the iteration matrix.
      ipar(iiwork-1+5) = -1
      ipar(iiwork-1+16) = lenpd

c     Reset the ABD information in DASSL.
      ipar(iiwork-1+19) = nint

      if (nint .lt. ninold) then
         call dcopy(nintmx+1+40+6*neqpre, rpar(irold),
     &              1, rpar(ipar(ixold)), 1)
      else
         if (nint .gt. ninold) then
            call dcopy(nintmx+1+40+6*neqpre, rpar(irold),
     &                 -1, rpar(ipar(ixold)), -1)
         endif
      endif

c     Reset the first 40 elements of the floating point work array in
c     DASSL.
      call dcopy(40, rpar(ipar(idasr)), 1, rpar(ipar(irwork)), 1)

      lenri1 = lenpd1 + kcol + 1 + nconti + (kcol + 1) * (ninpre + 1)
     &         + 2 * nconti
      lenri2 = lenpd2 + kcol + 1 + nconti + (kcol + 1) * (ninpre + 1)
     &         + 2 * nconti

      call meshsq(kcol, nint, x, rpar(ipar(irwork)+40), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

      if (irshfg .ne. 2) then
         nstep = ipar(idasi-1+8) + 1
      else
         nstep = 1
      endif

!              write (*,*) 'before reinit1: icount = ', icount

      call reinit(npde, kcol, kcol+1, nint, ninpre, ncpts1, neq1,
     &            neqpre, icount, ipar(idasi-1+11), nstep, x,
     &            rpar(ipar(ixold)), rpar(ipar(iypre)), 0,
     &            rpar(ipar(irwork)+40+6*neq1), lenri1, ipar(ipivot),
     &            rpar(ipar(ih)), rpar(ipar(ixbs1)),
     &            rpar(ipar(ixcol1)),
     &            rpar(ipar(ibasi1)), rpar(ipar(irwork)+40),
     &            rpar(ipar(iabbk1)), icflag)

!           write (*,*) 'after reinit1: icount,icflag = ', icount, icflag

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

c     cold start.
      if (irshfg .eq. 2) then
         call dcopy(neq1, rpar(ipar(irwork)+40), 1, rpar(ipar(iy1)), 1)
         call iniyp(t0, npde, kcol, nint, neq1, ncpts1,
     &              rpar(ipar(ixcol1)), rpar(ipar(iabtp1)),
     &              rpar(ipar(iabbk1)), rpar(ipar(iabbt1)),
     &              rpar(ipar(ibasi1)), rpar(ipar(iy1)),
     &              rpar(ipar(iyp1)), ipar(ipivot),
     &              rpar(ipar(irwork)), lenin1, icflag)
         goto 320
      endif

      call divdif(neq1, ipar(idasi-1+8)+1, rpar(ipar(idasr)-1+29),
     &            rpar(ipar(iework)), rpar(ipar(irwork)+40))

  320 continue

!            write (*,*) 'before reinit2: icount = ', icount

      call reinit(npde, kcol+1, kcol+1, nint, ninpre, ncpts2,
     &            neq2, neqpre, icount, ipar(idasi-1+11),
     &            nstep, x, rpar(ipar(ixold)), rpar(ipar(iypre)), 1,
     &            rpar(ipar(irwork)+40+6*neq), lenri2,
     &            ipar(ipivot), rpar(ipar(ih)), rpar(ipar(ixbs2)),
     &            rpar(ipar(ixcol2)), rpar(ipar(ibasi2)),
     &            rpar(ipar(irwork)+40+6*neq1), rpar(ipar(iabbk2)),
     &            icflag)

!           write (*,*) 'after reinit12 icount,icflag = ', icount, icflag

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

c     cold start.
      if (irshfg .eq. 2) then
         call dcopy(neq2, rpar(ipar(irwork)+40+6*neq1), 1,
     &              rpar(ipar(iy2)), 1)
         call iniyp(t0, npde, kcol+1, nint, neq2, ncpts2,
     &              rpar(ipar(ixcol2)), rpar(ipar(iabtp2)),
     &              rpar(ipar(iabbk2)), rpar(ipar(iabbt2)),
     &              rpar(ipar(ibasi2)), rpar(ipar(iy2)),
     &              rpar(ipar(iyp2)), ipar(ipivot),
     &              rpar(ipar(irwork)), lenin2, icflag)
         goto 400
      endif

      call dcopy(ipar(idasi-1+8)*neq2, rpar(ipar(irwork)+40+6*neq1),
     &           1, rpar(ipar(iyprer)), 1)

      call divdif(neq2, ipar(idasi-1+8)+1, rpar(ipar(idasr)-1+29),
     &            rpar(ipar(iework)), rpar(ipar(irwork)+40+6*neq1))

      ii = ipar(irwork) + 40
      jj = ii + 3 * neq
      kk = ii + 6 * neq1

      do 330 i = ipar(iiwork-1+8)+1, 1, -1
         call dcopy(neq2, rpar(kk+(i-1)*neq2), -1,
     &              rpar(jj+(i-1)*neq+neq1), -1)
         call dcopy(neq1, rpar(ii+(i-1)*neq1), -1, rpar(jj+(i-1)*neq),
     &              -1)
  330 continue

c-----------------------------------------------------------------------
c     Time integration loop for DASSL.
c-----------------------------------------------------------------------

  400 continue

      call ddassl(res, neq, t0, rpar(ipar(iy1)), rpar(ipar(iyp1)),
     &            tout, ipar(iinfo), rpar(ipar(itrtol)),
     &            rpar(ipar(itatol)), idid, rpar(ipar(irwork)), lenrw,
     &            ipar(iiwork), leniw, rpar, ipar, jac)

c-----------------------------------------------------------------------
c     Check for a successful time step and decide whether to continue
c     integration or to perform a remeshing.
c-----------------------------------------------------------------------
!              write (*,*) '(bacol) after ddassl: idid = ', idid
      if (idid .le. 0) goto 600
!             write (*,*) '(bacol) before errest'

      call berrest(kcol, nint, npde, neq1, neq2, necpts, icount,
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)),
     &            rpar(ipar(ixbs1)), rpar(ipar(ixbs2)),
     &            rpar(ipar(iey1)), rpar(ipar(iey2)), istart, mflag(2),
     &            atol, rtol, lenerr, rpar(ipar(iework)),
     &            rpar(ipar(iebas1)), rpar(ipar(iebas2)),
     &            errrat, rpar(ipar(ierint)), rpar(ipar(iercom)),
     &            ieflag)

!             write (*,*) '(bacol) after errest: ieflag,icount = ',
!      &        ieflag, icount

      if (ieflag .eq. 0) then

c        The current step is accepted.
         if (icount .ne. 0) then

            ipar(istalr) = 1

            ipar(irpstp) = ipar(iypre) + 6 * neq2 - 1
c           Check for a sufficiently large rpar floating point
c           work array.
            if (lrp .lt. ipar(irpstp)) goto 770

         else
            ipar(istalr) = ipar(istalr) + 1
         endif

c        Update the backup information.
         call sucstp(ipar(iiwork-1+11), ipar(iiwork-1+8)+1, icount,
     &               neq2, ipar(iiwork), rpar(ipar(irwork)),
     &               rpar(ipar(iey2)), rpar(ipar(iyprer)), ipar(idasi),
     &               rpar(ipar(idasr)), rpar(ipar(iypre)))
         icount = 0
         istart = 1
         irshfg = 0

c        Check whether the integration is done or not.
         if (mflag(4) .eq. 0) then
            if (t0 .lt. tout) then
               goto 400
            else
               goto 500
            endif
         else
            isstep = isstep + 1
            if ((t0 .lt. tout) .and. (isstep .lt. ipar(iinstp))) then
               goto 400
            else
               goto 500
            endif
         endif

      else
!             write (*,*) '(bacol) before remesh: icount = ', icount
c        The current step is rejected.
         if (icount .eq. maxrsh) goto 610

c        For the first remeshing at the current step, save nintpre and
c        neqpre at the last successful step.
         if (icount .eq. 0) then
            ninpre = nint
            neqpre = neq2
         endif

         ninold = nint
         irold = ipar(ixold)

c         print *, 'before remesh, nint = ', nint

         call remesh(istart, icount, nintmx, ninpre, ninold,
     &               errrat, rpar(ipar(ierint)), irshfg,
     &               rpar(ipar(ixold)), nint, kcol, x,
     &               rpar(ipar(iework)))
c         print *, 'remesh here', ', nint now is ', nint

         if (istart .eq. 1) then

c           This is not the initial step.
            t0 = rpar(ipar(idasr)-1+4)

c           In the first step after a remeshing, we do not allow DASSL
c           to increase the step size.
            if (rpar(ipar(idasr)-1+3) .gt. rpar(ipar(idasr)-1+7)) then
               rpar(ipar(idasr)-1+3) = rpar(ipar(idasr)-1+7)
            endif

c           In the first step after a remeshing, we do not allow DASSL
c           to increase the order of BDF method.
            if (ipar(idasi-1+7) .gt. ipar(idasi-1+8)) then
               ipar(idasi-1+7) = ipar(idasi-1+8)
            endif

            if (irshfg .eq. 2) then
c              This is a cold start
               ipar(istblc) = ipar(istblc) + ipar(idasi-1+11)
               ipar(icolds) = ipar(icolds) + 1
               ipar(iinfo) = 0
            endif

         else

c           This is the initial step.
            t0 = torign
            ipar(iinfo) = 0

         endif
         goto 100

      endif

c-----------------------------------------------------------------------
c     Successful return section.
c-----------------------------------------------------------------------
  500 continue

c     Retrieve the value of mflag(1).
      mflag(1) = ipar(iinfo)

c     Retrieve the output vector y from the rpar communication array.
      do 510 i = 1, neq1
         y(i) = rpar(ipar(iy1)-1+i)
  510 continue

c     Retrieve information on the time stepping from the ipar array.
      ipar(isteps) = ipar(istblc) + ipar(iiwork-1+11)

c     Retrieve the value of ipar(iinstp) when mflag(4) = 1.
      if (mflag(4) .eq. 1) ipar(iinstp) = isstep

      return

c-----------------------------------------------------------------------
c     Unsuccessful return section.
c-----------------------------------------------------------------------
  600 continue
!        write(6,9999) 'ERROR: BACOL runtime error in time stepping.'
!        write(6,9999) '       An error code and message should have'
!        write(6,9999) '       been issued by DASSL.'
      return
  610 continue
!        if (istart .eq. 1) then
!           write(6,9998) 'ERROR: BACOL has remeshed ', maxrsh,
!      &    ' times at',' t0 =', rpar(ipar(idasr)-1+4)
!        else
!           write(6,9998) 'ERROR: BACOL has remeshed ', maxrsh,
!      &    ' times at',' t0 =', torign
!        endif
      return

c-----------------------------------------------------------------------
c     The following section is the return point for invalid input.
c-----------------------------------------------------------------------

  710 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require:  0 <= mflag(i) <= 1, i = 1, 2, ..., 6.'
      endif
      idid = -51
      return

  715 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999)
     &  'Require:  if mflag(4) = 1, ipar(8) must be set to'
        write(6,9999) 'be a positive integer.'
      endif
      idid = -52
      return

  720 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require:  if mflag(6) = 1, tout must be in front'
        write(6,9999) 'of t0.'
      endif
      idid = -53
      return

  725 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require:  if mflag(6) = 1, rpar(2) must be the'
        write(6,9999) 'initial stepsize, thus nonzero.'
      endif
      idid = -54
      return

  730 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: npde > 0.'
      endif
      idid = -55
      return

  740 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: 1 < kcol <=', mxkcol, '.'
      endif
      idid = -56
      return

  750 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: 0 < nint <=', nintmx, '.'
      endif
      idid = -57
      return

  760 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: x(1) < x(2) < ... < x(nint+1).'
      endif
      idid = -58
      return

  770 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: lrp >= ', ipar(irpstp), '.'
      endif
      idid = -59
      return

  780 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'Require: lip >= ', ipar(iipstp), '.'
      endif
      idid = -60
      return

  790 continue
      if (lout .gt. 1) then
        write(6,9999) 'ERROR: BACOL input violation.'
        write(6,9999) 'IDID .ne. 1, on a continuation call of BACOL'
        write(6,9999) 'If IDID > 1, set idid = 1 and tout (t0 < tout)'
        write(6,9999)
     &  'If IDID < -1, the code cannot be continued due to'
        write(6,9999) '              a previous error.'
      endif
      idid = -61
      return

c-----------------------------------------------------------------------
 9998 format(a,i4,a,a,e12.5)
 9999 format(a,i4,a,i4,a,i4,a,i4,a,i4)
c-----------------------------------------------------------------------
      end
