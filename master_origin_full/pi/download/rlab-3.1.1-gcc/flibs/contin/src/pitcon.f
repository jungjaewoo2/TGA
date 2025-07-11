C  PCON61.FOR  27 February 1991
C
      SUBROUTINE PITCON(DF, FPAR, FX, IERROR, IPAR, IWORK, LIW,
     *                  NVAR, RWORK, LRW, XR, SLNAME,
     *                  STDOUT, LSTDOUT
     *                  )
C
C***********************************************************************
C
C  This is the interface program between the user and the PITCON package.
C
C
C  A. Introduction
C
C
C  PCON61 is version 6.1 of PITCON, the University of Pittsburgh continuation
C  program.
C
C  PITCON was written by
C
C    Professor Werner C Rheinboldt and John Burkardt,
C    Department of Mathematics and Statistics
C    University of Pittsburgh,
C    Pittsburgh, Pennsylvania, 15260, USA.
C
C    E-Mail: wcrhein@vms.cis.pitt.edu
C            burkardt@psc.edu
C
C  The original work on this package was partially supported by the National
C  Science Foundation under grants MCS-78-05299 and MCS-83-09926.
C
C  PITCON computes a sequence of solution points along a one dimensional
C  manifold of a system of nonlinear equations F(X)=0 involving NVAR-1
C  equations and an NVAR dimensional unknown vector X.
C
C  The operation of PITCON is somewhat analogous to that of an initial value
C  ODE solver.  In particular, the user must begin the computation by
C  specifying an approximate initial solution, and subsequent points returned
C  by PITCON lie on the curve which passes through this initial point and is
C  implicitly defined by F(X)=0.  The extra degree of freedom in the system is
C  analogous to the role of the independent variable in a differential
C  equations.
C
C  However, PITCON does not try to solve the algebraic problem by turning it
C  into a differential equation system.  Unlike differential equations, the
C  solution curve may bend and switch back in any direction, and there may be
C  many solutions for a fixed value of one of the variables.  Accordingly,
C  PITCON is not required to parametrize the implicitly defined curve with a
C  fixed parameter.  Instead, at each step, PITCON selects a suitable variable
C  as the current parameter and then determines the other variables as
C  functions of it.  This allows PITCON to go around relatively sharp bends.
C  Moreover, if the equations were actually differentiated - that is, replaced
C  by some linearization - this would introduce an inevitable "drift" away from
C  the true solution curve.  Errors at previous steps would be compounded in a
C  way that would make later solution points much less reliable than earlier
C  ones.  Instead, PITCON solves the algebraic equations explicitly and each
C  solution has to pass an acceptance test in an iterative solution process
C  with tolerances provided by the user.
C
C  PITCON is only designed for systems with one degree of freedom.  However,
C  it may be used on systems with more degrees of freedom if the user reduces
C  the available degrees of freedom by the introduction of suitable constraints
C  that are added to the set of nonlinear equations.  In this sense, PITCON may
C  be used to investigate the equilibrium behavior of physical systems with
C  several degrees of freedom.
C
C  Program options include the ability to search for solutions for which a
C  given component has a specified value.  Another option is a search for a
C  limit or turning point with respect to a given component; that is, of a
C  point where this particular solution component has a local extremum.
C
C  Another feature of the program is the use of two work arrays, IWORK and
C  RWORK.  All information required for continuing any interrupted computation
C  is saved in these two arrays.
C
C
C  B. PITCON Calling Sequence
C
C
C  SUBROUTINE PITCON(DF,FPAR,FX,IERROR,IPAR,IWORK,LIW,NVAR,RWORK,LRW,XR,SLVNAM)
C
C  On the first call, PITCON expects a point XR and a routine FX defining a
C  nonlinear function F.  As noted earlier, XR and FX specify a curve which
C  PITCON is to trace. On the first call, PITCON simply verifies that F(XR)=0,
C  and if this is not the case, the program attempts to correct XR to a new
C  value satisfying the equation.
C
C  On subsequent calls, PITCON assumes that the input vector XR contains the
C  point which had been computed on the previous call.  It also assumes that
C  the work arrays IWORK and RWORK contain the results of the prior
C  calculations.  PITCON estimates an appropriate stepsize, computes the
C  tangent direction to the curve at the given input point, and calculates a
C  predicted new point on the curve.  A form of Newton's method is used to
C  correct this point so that it lies on the curve.  If the iteration is
C  successful, the code returns with a new point XR.  Otherwise, the stepsize
C  may be reduced, and the calculation retried.
C
C  Aside from its ability to produce successive points on the solution curve,
C  PITCON may be asked to search for "target points" or "limit points".
C  Target points are solution vectors for which a certain component has a
C  specified value.  One might ask for all solutions for which XR(17)=4.0, for
C  instance.  Limit points occur when the curve turns back in a given
C  direction, and have the property that the corresponding component of the
C  tangent vector vanishes there.
C
C  If the user has asked for the computation of target or limit points, then
C  PITCON will usually proceed as described earlier, producing a new
C  continuation point on each return.  But if a target or limit point is
C  found, such a point is returned as the value of XR, temporarily interrupting
C  the usual form of the computation.
C
C
C  C. Overview of PITCON parameters:
C
C
C  DF     Input,        EXTERNAL DF, routine for evaluating the Jacobian of F.
C  FPAR   Input/output, REAL FPAR(*), user defined parameter array.
C  FX     Input,        EXTERNAL FX, routine for evaluating the function F.
C  IERROR       Output, INTEGER IERROR, error return flag.
C  IPAR   Input/output, INTEGER IPAR(*), user defined parameter array.
C  IWORK  Input/output, INTEGER IWORK(LIW).  Communication and work array.
C  LIW    Input,        INTEGER LIW, the dimension of IWORK.
C  NVAR   Input,        INTEGER NVAR, number of variables, the dimension of X.
C  RWORK  Input/output, REAL RWORK(LRW), Communication and work space array.
C  LRW    Input,        INTEGER LRW, the dimension of RWORK.
C  XR     Input/output, REAL XR(NVAR), the current solution point.
C  SLVNAM Input,        EXTERNAL SLVNAM, the solver to use on linear systems.
C
C
C  D. Details of PITCON parameters:
C
C
C  DF     Input, EXTERNAL DF, the name of the Jacobian evaluation routine.
C         This name must be declared EXTERNAL in the calling program.
C
C         DF is not needed if the finite difference option is used
C         (IWORK(9)=1 or 2). In that case, only a dummy name is needed for DF.
C
C         Otherwise, the user must write a routine which evaluates the
C         Jacobian matrix of the function FX at a given point X and stores it
C         in the FJAC array in accordance with the format used by the solver
C         specified in SLVNAM.
C
C         In the simplest case, when the full matrix solverDENSLV solver
C         provided with the package is used, DF must store  D F(I)/D X(J) into
C         FJAC(I,J).
C
C         The array to contain the Jacobian will be zeroed out before DF is
C         called, so that only nonzero elements need to be stored.  DF must
C         have the form:
C
C           SUBROUTINE DF(NVAR,FPAR,IPAR,X,FJAC,IERROR)
C
C           NVAR   Input, INTEGER NVAR, number of variables.
C
C           FPAR   Input, REAL FPAR(*), vector for passing real parameters.
C                  This vector is not used by the program, and is only provided
C                  for the user's convenience.
C
C           IPAR   Input, INTEGER IPAR(*), vector for passing integer
C                  parameters.  This vector is not used by the program, and is
C                  only provided for the user's convenience.
C
C           X      Input, REAL X(NVAR), the point at which the Jacobian is
C                  desired.
C
C           FJAC   Output, REAL FJAC(*), the array containing the Jacobian.
C
C                  If DENSLV is the solver:  FJAC must be dimensioned
C                  FJAC(NVAR,NVAR) as shown above, and DF sets
C                  FJAC(I,J)=D F(I)/DX(J).
C
C                  If BANSLV is the solver:  the main portion of the Jacobian,
C                  rows and columns 1 through NVAR-1, is assumed to be a banded
C                  matrix in the standard LINPACK form with lower bandwidth ML
C                  and upper bandwidth MU.  However, the final column of the
C                  Jacobian is allowed to be full.
C
C                  BANSLV will pass to DF the beginning of the storage for
C                  FJAC, but it is probably best not to doubly dimension FJAC
C                  inside of DF, since it is a "hybrid" object.  The first
C                  portion of it is a (2*ML+MU+1, NEQN) array, followed by a
C                  single column of length NEQN (the last column of the
C                  Jacobian).  Thus the simplest approach is to declare FJAC to
C                  be a vector, and then then to store values as follows:
C
C                    If J is less than NVAR, then
C                      if I-J .LE. ML and J-I .LE. MU,
C                        set K=(2*ML+MU)*J + I - ML
C                        set FJAC(K)=D F(I)/DX(J).
C                      else
C                        do nothing, index is outside the band
C                      endif
C                    Else if J equals NVAR, then
C                      set K=(2*ML+MU+1)*(NVAR-1)+I,
C                      set FJAC(K)=D F(I)/DX(J).
C                    endif.
C
C           IERROR Output, INTEGER IERROR, error return flag.  DF should set
C                  this to 0 for normal return, nonzero if trouble.
C
C  FPAR   Input/output, REAL FPAR(*), a user defined parameter array.
C
C         FPAR is not used in any way by PITCON.  It is provided for the user's
C         convenience.  It is passed to DF and FX, and hence may be used to
C         transmit information between the user calling program and these user
C         subprograms. The dimension of FPAR and its contents are up to the
C         user.  Internally, the program declares DIMENSION FPAR(*) but never
C         references its contents.
C
C  FX     Input, EXTERNAL FX, the name of the routine which evaluates the
C         function.
C
C         FX computes the value of the nonlinear function.  This name must be
C         declared EXTERNAL in the calling program.  FX should evaluate the
C         NVAR-1 function components at the input point X, and store the result
C         in the vector FVEC.  An augmenting equation will be stored in entry
C         NVAR of FVEC by the PITCON program.
C
C         FX should have the following form:
C
C           SUBROUTINE FX(NVAR,FPAR,IPAR,X,FVEC,IERROR)
C
C           NVAR   Input, INTEGER NVAR, number of variables.
C
C           FPAR   Input/output, REAL FPAR(*), array of user parameters.
C
C           IPAR   Input/output, INTEGER IPAR(*), array of user parameters.
C
C           X      Input, REAL X(NVAR), the point at which function evaluation
C                  is required.
C
C           FVEC   Output, REAL FVEC(NVAR), the value of the function at point
C                  X.  Only the first NVAR-1 entries of FVEC are to be set by
C                  the routine.  PITCON sets the final value itself.
C
C           IERROR Output, INTEGER IERROR, error flag.  FX should set this to 0
C                  if there are no problems, or to 2 if there is a problem.
C
C  IERROR Output, INTEGER IERROR, error return flag.
C
C         On return from PITCON, a nonzero value of IERROR is a warning of some
C         problem which may be minor, serious, or fatal.
C
C         0, No errors occurred.
C
C         1, Insufficient storage provided in RWORK and IWORK, or NVAR is less
C            than 2.  This is a fatal error, which occurs on the first call to
C            PITCON.
C
C         2, A user defined error condition occurred in the FX or DF
C            subroutines.  PITCON treats this as a fatal error.
C
C         3, A numerically singular matrix was encountered.  Continuation
C            cannot proceed without some redefinition of the problem.  This is
C            a fatal error.
C
C         4, Unsuccessful corrector iteration.  Loosening the tolerances
C            RWORK(1) and RWORK(2), or decreasing the maximum stepsize RWORK(4)
C            might help.  This is a fatal error.
C
C         5, Too many corrector steps.  The corrector iteration was proceeding
C            properly, but too slowly.  Increase number of Newton steps
C            IWORK(17), increase the error tolerances RWORK(1) or RWORK(2), or
C            decrease RWORK(4).  This is a fatal error.
C
C         6, Null tangent vector.  A serious error which indicates a data
C            problem or singularity in the nonlinear system.  This is a fatal
C            error.
C
C         7, Root finder failed while searching for a limit point.
C            This is a warning.  It means that the target or limit point
C            computation has failed, but the main computation (computing the
C            continuation curve itself) may continue.
C
C         8, Limit point iteration took too many steps.  This is a warning
C            error.  It means that the limit point computation has failed, but
C            the main computation (computing the continuation curve itself) may
C            continue.
C
C         9, Undiagnosed error condition.  This is a fatal error.
C
C  IPAR   Input/output, INTEGER IPAR(*), user defined parameter array.
C
C         IPAR is not used in any way by PITCON.  It is provided for the user's
C         convenience in transmitting parameters between the calling program
C         and the user routines FX and DF.  IPAR is declared in the PITCON
C         program and passed through it to FX and DF, but otherwise ignored.
C         Note, however, that if BANSLV is used for the solver routine, then
C         IPAR(1) must contain the lower bandwidth, and IPAR(2) the upper
C         bandwidth of the Jacobian matrix.
C
C  IWORK  Input/output, INTEGER IWORK(LIW).  Communication and workspace array.
C
C         The specific allocation of IWORK is described in the section devoted
C         to the work arrays.  Some elements of IWORK MUST be set by the user,
C         others may be set to change the way PITCON works.
C
C  LIW    Input, INTEGER LIW, the dimension of IWORK.
C
C         The minimum acceptable value of LIW depends on the solver chosen,
C         but for either DENSLV or BANSLV, setting LIW=29+NVAR is sufficient.
C
C  NVAR   Input, INTEGER NVAR, the number of variables, the dimension of X.
C
C         This is, of course, one greater than the number of equations or
C         functions.  NVAR must be at least 2.
C
C  RWORK  Input/output, REAL RWORK(LRW), work array.
C
C         The specific allocation of RWORK is described in the section
C         devoted to the work arrays.
C
C  LRW    Input, INTEGER LRW, the dimension of RWORK.
C
C         The minimum acceptable value depends heavily on the solver options.
C         There is storage required for scalars, vectors, and the Jacobian
C         array.  The minimum acceptable value of LRW is the sum of three
C         corresponding numbers.
C
C         For DENSLV with user-supplied Jacobian,
C
C           LRW=29 + 4*NVAR + NVAR*NVAR.
C
C         For DENSLV with internally approximated Jacobian,
C
C           LRW=29 + 6*NVAR + NVAR*NVAR.
C
C         For BANSLV, with a Jacobian matrix with upper bandwidth MU and lower
C         bandwidth ML, and NBAND=2*ML+MU+1, with user supplied Jacobian,
C
C           LRW=29 + 6*NVAR + (NVAR-1)*NBAND.
C
C         For BANSLV with internally approximated Jacobian,
C
C           LRW=29 + 9*NVAR + (NVAR-1)*NBAND.
C
C  XR     Input/output, REAL XR(NVAR), the current solution point.
C
C         On the first call, the user should set XR to a starting point which
C         at least approximately satisfies F(XR)=0.  The user need never
C         update XR again.
C
C         Thereafter, on each return from the program with IERROR=0, XR will
C         hold the most recently computed point, whether a continuation, target
C         or limit point.
C
C  SLVNAM Input, EXTERNAL SLVNAM, the name of the solver to use on linear
C         systems.
C
C         The linear systems have the form A*x=b, where A is the augmented
C         Jacobian matrix.  A will be square, and presumably nonsingular.
C         The routine must return the value of the solution x.
C
C         Two possible choices for SLVNAM are "DENSLV" and "BANSLV", which are
C         the names of routines provided with the package.  DENSLV is
C         appropriate for a full storage jacobian, and BANSLV for a jacobian
C         which is banded except for the last column.
C
C         The advanced user may study the source code for these two routines
C         and write an equivalent solver more suited to a given problem.
C
C
C  E. The Integer Work Array IWORK
C
C
C  Input to the program includes the setting of some of the entries in IWORK.
C  Some of this input is optional.  The user input section of IWORK involves
C  entries 1 through 8, and, possibly also 17 and 29.
C
C  IWORK(1) must be set by the user.  All other entries have default values.
C
C
C  IWORK(1)        On first call only, the user must set IWORK(1)=0.
C                  Thereafter, the program sets IWORK(1) before return to
C                  explain what kind of point is being returned.  This return
C                  code is:
C
C                      1 return of corrected starting point.
C                      2 return of continuation point.
C                      3 return of target point.
C                      4 return of limit point.
C
C                  NOTE:  At any time, PITCON may be called with a negative
C                  value of IWORK(1). This requests a check of the
C                  jacobian routine against a finite difference approximation.
C                  The program will call the jacobian routine, and
C                  then estimate the jacobian.  If IWORK(1)=-1, then it will
C                  print out the value of the maximum difference, and the row
C                  and column of the jacobian in which it appears.  Otherwise,
C                  the program will print out the entire matrix
C                  FP(I,J)-DEL(J)F(I), where DEL(J)F(I) represents the finite
C                  difference approximation.
C
C                  Before a call with negative IWORK(1), the current value of
C                  IWORK(1) should be saved, and then restored to the previous
C                  value after the call, in order to resume calculation.
C
C                  IWORK(1) does not have a default value.  The user MUST set
C                  it.
C
C  IWORK(2)        The component of the current continuation point XR which is
C                  to be used as the continuation parameter.  On first call,
C                  the program is willing to use the index NVAR as a default,
C                  but the user should set this value if better information is
C                  available.
C
C                  After the first call, the program sets this value for each
C                  step automatically unless the user prevents this by setting
C                  the parameterization option IWORK(3) to a non-zero valus.
C                  Note that a poor choice of IWORK(2) may cause the algorithm
C                  to fail.  IWORK(2) defaults to NVAR on the first step.
C
C  IWORK(3)        Parameterization option.  The program would prefer to be
C                  free to choose a new local parameter from step to step.
C                  The value of IWORK(3) allows or prohibits this action.
C                  IWORK(3)=0 allows the program to vary the parameter,
C                  IWORK(3)=1 forces the program to use whatever the contents
C                  of IWORK(2) are, which will not be changed from the user's
C                  input or the default.  The default is IWORK(3)=0.
C
C  IWORK(4)        Newton iteration Jacobian update option.
C                  0, the Jacobian is reevaluated at every step of the
C                     Newton iteration.  This is costly, but may result in
C                     fewer Newton steps and fewer Newton iteration rejections.
C                  1, the Jacobian is evaluated only on the first and
C                     IWORK(17)-th steps of the Newton process.
C                  2, the Jacobian is evaluated only when absolutely
C                     necessary, namely, at the first step, and when the
C                     process fails. This option is most suitable for problems
C                     with mild nonlinearities.
C
C                  The default is IWORK(4)=0.
C
C  IWORK(5)        Target point index.  If IWORK(5) is not zero, it is presumed
C                  to be the component index between 1 and NVAR for which
C                  target points are sought.  In this case, the value of
C                  RWORK(7) is assumed to be the target value.  The program
C                  will monitor every new continuation point, and if it finds
C                  that a target point may lie between the new point and the
C                  previous point, it will compute this target point and
C                  return.  This target point is defined by the property that
C                  its component with the index prescribed in IWORK(5) will
C                  have the value given in RWORK(7).  For a given problem there
C                  may be zero, one, or many target points requested.
C                  The default of IWORK(5) is 0.
C
C  IWORK(6)        Limit point index.  If IWORK(6) is nonzero, then the program
C                  will search for limit points with respect to the component
C                  with index IWORK(6); that is, of points for which the
C                  IWORK(6)-th variable has a local extremum, or equivalently
C                  where the IWORK(6)-th component of the tangent vector is
C                  zero.  The default of IWORK(6) is zero.
C
C  IWORK(7)        Control of the amount of intermediate output produced by the
C                  program. IWORK(7) may have a value between 0 and 3.
C                  For IWORK(7) = 0 there is no intermediate output while for
C                  IWORK(7) = 3 the most intermediate output is produced.
C                  The default is 1.
C
C  IWORK(8)        FORTRAN unit to which output is to be written.  The
C                  default value is site dependent but should be the standard
C                  output device.
C                  The default is 6 on the Cray, Vax and PC, or 9 on the
C                  Macintosh.
C
C  IWORK(9)        Control of the Jacobian option specifying whether the user
C                  has supplied a Jacobian routine, or wants the program
C                  to approximate the Jacobian.
C                  0, the user has supplied the Jacobian.
C                  1, program is to use forward difference approximation.
C                  2, program is to use central difference approximation.
C                  IWORK(9) defaults to 0.
C
C  IWORK(10)       State indicator of the progress of the program.
C                  The values are:
C                  0, start up with unchecked starting point.
C                  1, first step.  Corrected starting point available.
C                  2, two successive continuation points available, as well
C                     as the tangent vector at the oldest of them.
C                  3, two successive continuation points available, as well
C                     as the tangent vector at the newest of them.
C
C  IWORK(11)       Index of the last computed target point. This is used to
C                  avoid repeated computation of a target point.  If a target
C                  point has been found, then the target index IWORK(5) is
C                  copied into IWORK(11).
C
C  IWORK(12)       Second best choice for the local parameterization index.
C                  This index may be tried if the first choice causes poor
C                  performance in the Newton corrector.
C
C  IWORK(13)       Beginning location in IWORK of unused integer work space
C                  available for use by the solver.
C
C  IWORK(14)       LIW, the user declared dimension of the array IWORK.
C
C  IWORK(15)       Beginning location in RWORK of unused real work space
C                  available for use by the solver.
C
C  IWORK(16)       LRW, the user declared dimension of RWORK.
C
C  IWORK(17)       Maximum number of corrector steps allowed during one run
C                  of the Newton process in which the Jacobian is updated at
C                  every step.  If the Jacobian is only evaluated at
C                  the beginning of the Newton iteration then 2*IWORK(17) steps
C                  are allowed.
C                  IWORK(17) must be greater than 0.  It defaults to 10.
C
C  IWORK(18)       Number of stepsize reductions that were needed for
C                  producing the last continuation point.
C
C  IWORK(19)       Total number of calls to the user Jacobian routine DF.
C
C  IWORK(20)       Total number of calls to the matrix factorization routine.
C                  If DENSLV is the chose solver then factorization is done by
C                  the LINPACK routine SGEFA.  If BANSLV is the solver, the
C                  LINPACK routine SGBFA will be used.
C
C  IWORK(21)       Total number of calls to the back-substitution routine.
C                  If DENSLV is the chosen solver, then back substitution is
C                  done by the LINPACK routine SGESL.  If BANSLV is used, then
C                  the LINPACK routine SGBSL will be used.
C
C  IWORK(22)       Total number of calls to the user function routine FX.
C
C  IWORK(23)       Total number of steps taken in limit point iterations.
C                  Each step involves determining an approximate limit point
C                  and applying a Newton iteration to correct it.
C
C  IWORK(24)       Total number of Newton corrector steps used during the
C                  computation of target points.
C
C  IWORK(25)       Total number of Newton steps taken during the correction
C                  of a starting point or the continuation points.
C
C  IWORK(26)       Total number of predictor stepsize-reductions needed
C                  since the start of the continuation procesds.
C
C  IWORK(27)       Total number of calls to the program.  This also
C                  corresponds to the number of points computed.
C
C  IWORK(28)       Total number of Newton steps taken during current iteration.
C
C  IWORK(30)       and on are reserved for use by the linear equation solver,
C                  and typically are used for pivoting.
C
C
C  F. The Real Work Array RWORK
C
C
C  Input to the program includes the setting of some of the entries in RWORK.
C  Some of this input is optional.  The user input section of RWORK involves
C  entries 1 through 7 and possibly 20.  All entries of RWORK have default
C  values.
C
C
C  RWORK(1)        Absolute error tolerance.   This value is used mainly during
C                  the Newton iteration.  RWORK(1) defaults to SQRT(EPMACH)
C                  where EPMACH is the machine relative precision stored in
C                  RWORK(8).
C
C  RWORK(2)        Relative error tolerance.  This value is used mainly during
C                  the Newton iteration.  RWORK(2) defaults to SQRT(EPMACH)
C                  where EPMACH is the machine relative precision stored in
C                  RWORK(8).
C
C  RWORK(3)        Minimum allowable predictor stepsize.  If failures of
C                  the Newton correction force the stepsize down to this level,
C                  then the program will give up.  The default value is
C                  SQRT(EPMACH).
C
C  RWORK(4)        Maximum allowable predictor step.  Too generous a value
C                  may cause erratic behavior of the program.  The default
C                  value is SQRT(NVAR).
C
C  RWORK(5)        Predictor stepsize.  On first call, it should be set by
C                  the user.  Thereafter it is set by the program.
C                  RWORK(5) should be positive.  In order to travel in the
C                  negative direction, see RWORK(6).
C                  The default initial value equals 0.5*(RWORK(3)+RWORK(4)).
C
C  RWORK(6)        The local continuation direction, which is either +1.0
C                  or -1.0 .  This asserts that the program is moving in the
C                  direction of increasing or decreasing values of the local
C                  continuation variable, whose index is in IWORK(2).  On first
C                  call, the user must choose IWORK(2).  Therefore, by setting
C                  RWORK(6), the user may also specify whether the program is
C                  to move initially to increase or decrease the variable whose
C                  index is IWORK(2).
C                  RWORK(6) defaults to +1.
C
C  RWORK(7)        A target value.  It is only used if a target index
C                  has been specified through IWORK(5).  In that case, solution
C                  points with the IWORK(5) component equal to RWORK(7) are
C                  to be computed. The code will return each time it finds such
C                  a point.  RWORK(7) does not have a default value.  The
C                  program does not set it, and it is not referenced unless
C                  IWORK(5) has been set.
C
C  RWORK(8)        EPMACH, the value of the machine precision.  The computer
C                  can distinguish 1.0+EPMACH from 1.0, but it cannot
C                  distinguish 1.0+(EPMACH/2) from 1.0. This number is used
C                  when estimating a reasonable accuracy request on a given
C                  computer.  PITCON computes a value for EPMACH internally.
C
C  RWORK(9)        Not currently used.
C
C  RWORK(10)       A minimum angle used in the steplength computation,
C                  equal to 2.0*ARCCOS(1-EPMACH).
C
C  RWORK(11)       Estimate of the angle between the tangent vectors at the
C                  last two continuation points.
C
C  RWORK(12)       The pseudo-arclength coordinate of the previous continuation
C                  pointl; that is, the sum of the Euclidean distances between
C                  all computed continuation points beginning with the start
C                  point.  Thus each new point should have a larger coordinate,
C                  except for target and limit points which lie between the two
C                  most recent continuation points.
C
C  RWORK(13)       Estimate of the pseudo-arclength coordinate of the current
C                  continuation point.
C
C  RWORK(14)       Estimate of the pseudoarclength coordinate of the current
C                  limit or target point, if any.
C
C  RWORK(15)       Size of the correction of the most recent continuation
C                  point; that is, the maximum norm of the distance between the
C                  predicted point and the accepted corrected point.
C
C  RWORK(16)       Estimate of the curvature between the last two
C                  continuation points.
C
C  RWORK(17)       Sign of the determinant of the augmented matrix at the
C                  last continuation point whose tangent vector has been
C                  calculated.
C
C  RWORK(18)       Not currently used.
C
C  RWORK(19)       Not currently used.
C
C  RWORK(20)       Maximum growth factor for the predictor stepsize based
C                  on the previous secant stepsize.  The stepsize algorithm
C                  will produce a suggested step that is no less that the
C                  previous secant step divided by this factor, and no greater
C                  than the previous secant step multiplied by that factor.
C                  RWORK(20) defaults to 3.
C
C  RWORK(21)       The (Euclidean) secant distance between the last two
C                  computed continuation points.
C
C  RWORK(22)       The previous value of RWORK(21).
C
C  RWORK(23)       A number judging the quality of the Newton corrector
C                  convergence at the last continuation point.
C
C  RWORK(24)       Value of the component of the current tangent vector
C                  corresponding to the current continuation index.
C
C  RWORK(25)       Value of the component of the previous tangent vector
C                  corresponding to the current continuation index.
C
C  RWORK(26)       Value of the component of the current tangent vector
C                  corresponding to the limit index in IWORK(6).
C
C  RWORK(27)       Value of the component of the previous tangent vector
C                  corresponding to the limit index in IWORK(6).
C
C  RWORK(28)       Value of RWORK(7) when the last target point was
C                  computed.
C
C  RWORK(29)       Sign of the determinant of the augmented matrix at the
C                  previous continuation point whose tangent vector has been
C                  calculated.
C
C  RWORK(30)       through RWORK(30+4*NVAR-1) are used by the program to hold
C                  an old and new continuation point, a tangent vector and a
C                  work vector.  Subsequent entries of RWORK are used by the
C                  linear solver.
C
C
C  G. Programming Notes
C
C
C  The minimal input and user routines required to apply the program are as
C  follows:
C    Write a function routine FX of the form described above.
C    Use DENSLV as the linear equation solver by setting SLVNAM to DENSLV.
C    Skip writing a Jacobian routine by using the finite difference option.
C    Pass the name of FX as the Jacobian name as well.
C    Declare the name of the function FX as EXTERNAL.
C    Set NVAR in accordance with your problem.
C
C  Then:
C
C    Dimension the vector IWORK to the size LIW=29+NVAR.
C    Dimension the vector RWORK to the size LRW=29+NVAR*(NVAR+6).
C    Dimension IPAR(1) and FPAR(1) as required in the function routine.
C    Dimension XR(NVAR) and set it to an approximate solution of F(XR)=0.
C
C  Set the work arrays as follows:
C
C    Initialize IWORK to 0 and RWORK to 0.0.
C
C    Set IWORK(1)=0 (Problem startup)
C    Set IWORK(7)=3 (Maximum internally generated output)
C    Set IWORK(9)=1 (Forward difference Jacobian)
C
C  Now call the program repeatedly, and never change any of its arguments.
C  Check IERROR to decide whether the code is working satisfactorily.
C  Print out the vector XR to see the current solution point.
C
C  The most obvious input to try to set appropriately after some practice
C  would be the error tolerances RWORK(1) and RWORK(2), the minimum, maximum
C  and initial stepsizes RWORK(3), RWORK(4) and RWORK(5), and the initial
C  continuation index IWORK(2).
C
C  For speed and efficiency, a Jacobian routine should be written. It can be
C  checked by comparing its results with the finite difference Jacobian.
C
C  For a particular problem, the target and limit point input can be very
C  useful.  For instance, in the case of a discretized boundary value problem
C  with a real parameter it may be desirable to compare the computed solutions
C  for different discretization-dimensions and equal values of the parameter.
C  For this the target option can be used with the prescribed values of the
C  parameter. Limit points usually are of importance in connection with
C  stability considerations.
C
C
C  H. References
C
C
C  1.
C  Werner Rheinboldt,
C  Solution Field of Nonlinear Equations and Continuation Methods,
C  SIAM Journal of Numerical Analysis,
C  Volume 17, 1980, pages 221-237.
C
C  2.
C  Cor den Heijer and Werner Rheinboldt,
C  On Steplength Algorithms for a Class of Continuation Methods,
C  SIAM Journal of Numerical Analysis,
C  Volume 18, 1981, pages 925-947.
C
C  3.
C  Werner Rheinboldt,
C  Numerical Analysis of Parametrized Nonlinear Equations
C  John Wiley and Sons, New York, 1986
C
C  4.
C  Werner Rheinboldt and John Burkardt,
C  A Locally Parameterized Continuation Process,
C  ACM Transactions on Mathematical Software,
C  Volume 9, Number 2, June 1983, pages 215-235.
C
C  5.
C  Werner Rheinboldt and John Burkardt,
C  Algorithm 596, A Program for a Locally Parameterized Continuation Process,
C  ACM Transactions on Mathematical Software,
C  Volume 9, Number 2, June 1983, Pages 236-241.
C
C  6.
C  J J Dongarra, J R Bunch, C B Moler and G W Stewart,
C  LINPACK User's Guide,
C  Society for Industrial and Applied Mathematics,
C  Philadelphia, 1979.
C
C  7.
C  Richard Brent,
C  Algorithms for Minimization without Derivatives,
C  Prentice Hall, 1973.
C
C  8.
C  Tony Chan,
C  Deflated Decomposition of Solutions of Nearly Singular Systems,
C  Technical Report 225,
C  Computer Science Department,
C  Yale University,
C  New Haven, Connecticut, 06520,
C  1982.
C
C
C  I.  Converting between REAL and DOUBLE PRECISION versions
C
C
C  PITCON was written to use single precision floating point arithmetic, and
C  for most problems, this has sufficed, even on machines for which the
C  default REAL variable is stored in 32 bits.  However, if PITCON cannot
C  solve a problem, the failure is sometimes caused by poor scaling of the
C  functions, a near-singular jacobian, or similar conditions that are
C  sometimes helped by using higher precision.
C
C  The current version of the code is fairly easy to convert.  The type of
C  each variable is declared in every routine.  Generic FORTRAN functions
C  are used.  All constants are declared in PARAMETER statements at the
C  beginning of each routine.
C
C  To create a DOUBLE PRECISION version from a REAL version, a user should
C  make the following global substitutions:
C
C    Replace         by
C
C    REAL            DOUBLE PRECISION  or  REAL*8
C
C    ISAMAX          IDAMAX
C    SAXPY           DAXPY
C    SCOPY           DCOPY
C    SDOT            DDOT
C    SGBDI           DGBDI
C    SGBFA           DGBFA
C    SGBSL           DGBSL
C    SGEDI           DGEDI
C    SGEFA           DGEFA
C    SGESL           DGESL
C    SNRM2           DNRM2
C    SSCAL           DSCAL
C    SSWAP           DSWAP
C
C  and to convert from DOUBLE PRECISION to REAL, the changes should be
C  reversed.
C
C  The same changes can be made to the sample programs.  In most cases, this
C  should be sufficient to create proper double precision programs.
C
C
C  J.  A sample calling program
C
C
C  The following is a sample program that calls PITCON to handle the
C  Freudenstein-Roth function.
C
C
C  C  PCPRB1.FOR  22 February 1991
C  C
C        PROGRAM PCPRB1
C  C
C  C  The Freudenstein-Roth function.
C  C
C  C  For a more complicated version of this example, see PCPRB8.
C  C
C  C  Reference
C  C
C  C  F Freudenstein, B Roth,
C  C  Numerical Solutions of Nonlinear Equations,
C  C  Journal of the Association for Computing Machinery,
C  C  Volume 10, 1963, Pages 550-556.
C  C
C  C  The function F(X) is of the form
C  C
C  C    FX(1) = X1 - X2**3 + 5*X2**2 -  2*X2 - 13 + 34*(X3-1)
C  C    FX(2) = X1 + X2**3 +   X2**2 - 14*X2 - 29 + 10*(X3-1)
C  C
C  C  Starting from the point (15,-2,0), the program is required to produce
C  C  solution points along the curve until it reaches a solution point
C  C  (*,*,1).  It also may be requested to look for limit points in the
C  C  first or third components.
C  C
C  C  The correct value of the solution at X3=1 is (5,4,1).
C  C
C  C  Limit points in the first variable occur at:
C  C
C  C    (14.28309, -1.741377,  0.2585779)
C  C    (61.66936,  1.983801, -0.6638797)
C  C
C  C  Limit points for the third variable occur at:
C  C
C  C    (20.48586, -0.8968053, 0.5875873)
C  C    (61.02031,  2.230139, -0.6863528)
C  C
C  C
C        INTEGER   LIW
C        INTEGER   LRW
C        INTEGER   NVAR
C  C
C        PARAMETER (NVAR=3)
C        PARAMETER (LIW=NVAR+29)
C        PARAMETER (LRW=29+(6+NVAR)*NVAR)
C  C
C        EXTERNAL  FXROTH
C        EXTERNAL  DFROTH
C        EXTERNAL  DENSLV
C        EXTERNAL  PITCON
C  C
C        REAL      FPAR(1)
C        INTEGER   I
C        INTEGER   IERROR
C        INTEGER   IPAR(1)
C        INTEGER   IWORK(LIW)
C        INTEGER   J
C        CHARACTER NAME*12
C        REAL      RWORK(LRW)
C        REAL      XR(NVAR)
C  C
C  C  Set work arrays to zero:
C  C
C        DO 10 I=1,LIW
C          IWORK(I)=0
C  10      CONTINUE
C        DO 20 I=1,LRW
C          RWORK(I)=0.0
C  20      CONTINUE
C  C
C  C  Set some entries of work arrays.
C  C
C  C  IWORK(1)=0 ; This is a startup
C  C  IWORK(2)=2 ; Use X(2) for initial parameter
C  C  IWORK(3)=0 ; Program may choose parameter index
C  C  IWORK(4)=0 ; Update jacobian every newton step
C  C  IWORK(5)=3 ; Seek target values for X(3)
C  C  IWORK(6)=1 ; Seek limit points in X(1)
C  C  IWORK(7)=1 ; Control amount of output.
C  C  IWORK(8)=6 ; Output unit
C  C  IWORK(9)=2 ; Jacobian choice.
C  C
C        IWORK(1)=0
C        IWORK(2)=2
C        IWORK(3)=0
C        IWORK(4)=0
C        IWORK(5)=3
C        IWORK(6)=1
C        IWORK(7)=1
C        IWORK(8)=6
C        IWORK(9)=0
C  C
C  C  RWORK(1)=0.00001; Absolute error tolerance
C  C  RWORK(2)=0.00001; Relative error tolerance
C  C  RWORK(3)=0.01   ; Minimum stepsize
C  C  RWORK(4)=20.0   ; Maximum stepsize
C  C  RWORK(5)=0.3    ; Starting stepsize
C  C  RWORK(6)=1.0    ; Starting direction
C  C  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
C  C
C        RWORK(1)=0.00001
C        RWORK(2)=0.00001
C        RWORK(3)=0.01
C        RWORK(4)=20.0
C        RWORK(5)=0.3
C        RWORK(6)=1.0
C        RWORK(7)=1.0
C  C
C  C  Set starting point.
C  C
C        XR(1)=15.0
C        XR(2)=-2.0
C        XR(3)=0.0
C  C
C        WRITE(6,*)' '
C        WRITE(6,*)'PCPRB1:'
C        WRITE(6,*)' '
C        WRITE(6,*)'PITCON sample program.'
C        WRITE(6,*)'Freudenstein-Roth function'
C        WRITE(6,*)'2 equations, 3 variables.'
C        WRITE(6,*)' '
C        WRITE(6,*)'Step  Type of Point     '//
C       *'X(1)         X(2)         X(3)'
C        WRITE(6,*)' '
C        I=0
C        NAME='Start point  '
C        WRITE(6,'(1X,I3,2X,A12,2X,3G14.6)')I,NAME,(XR(J),J=1,NVAR)
C  C
C        DO 30 I=1,30
C  C
C          CALL PITCON(DFROTH,FPAR,FXROTH,IERROR,IPAR,IWORK,LIW,
C       *  NVAR,RWORK,LRW,XR,DENSLV)
C  C
C          IF(IWORK(1).EQ.1)THEN
C            NAME='Corrected    '
C          ELSEIF(IWORK(1).EQ.2)THEN
C            NAME='Continuation '
C          ELSEIF(IWORK(1).EQ.3)THEN
C            NAME='Target point '
C          ELSEIF(IWORK(1).EQ.4)THEN
C            NAME='Limit point  '
C            ENDIF
C  C
C          WRITE(6,'(1X,I3,2X,A12,2X,3G14.6)')I,NAME,(XR(J),J=1,NVAR)
C  C
C          IF(IWORK(1).EQ.3)THEN
C            WRITE(6,*)' '
C            WRITE(6,*)'We have reached the point we wanted.'
C            WRITE(6,*)'The code may stop now.'
C            STOP
C            ENDIF
C  C
C          IF(IERROR.NE.0)THEN
C            WRITE(6,*)' '
C            WRITE(6,*)'An error occurred.'
C            WRITE(6,*)'We will terminate the computation now.'
C            STOP
C            ENDIF
C  C
C  30      CONTINUE
C  C
C        WRITE(6,*)' '
C        WRITE(6,*)'PITCON did not reach the point of interest.'
C        STOP
C        END
C        SUBROUTINE FXROTH(NVAR,FPAR,IPAR,X,F,IERROR)
C  C
C  C  Evaluate the function at X.
C  C
C  C  ( X1 - ((X2-5.0)*X2+2.0)*X2 - 13.0 + 34.0*(X3-1.0)  )
C  C  ( X1 + ((X2+1.0)*X2-14.0)*X2 - 29.0 + 10.0*(X3-1.0) )
C  C
C        REAL      F(*)
C        REAL      FPAR(*)
C        INTEGER   IERROR
C        INTEGER   IPAR(*)
C        REAL      X(NVAR)
C  C
C        F(1)=X(1)
C       *     -((X(2)-5.0)*X(2)+2.0)*X(2)
C       *     -13.0
C       *     +34.0*(X(3)-1.0)
C  C
C        F(2)=X(1)
C       *    +((X(2)+1.0)*X(2)-14.0)*X(2)
C       *    -29.0
C       *    +10.0*(X(3)-1.0)
C  C
C        RETURN
C        END
C        SUBROUTINE DFROTH(NVAR,FPAR,IPAR,X,FJAC,IERROR)
C  C
C  C  Evaluate the Jacobian:
C  C
C  C  ( 1.0   (-3.0*X(2)+10.0)*X(2)-2.0      34.0   )
C  C  ( 1.0   (3.0*X(2)+2.0)*X(2)-14.0       10.0   )
C  C
C        REAL      FJAC(NVAR,NVAR)
C        REAL      FPAR(*)
C        INTEGER   IERROR
C        INTEGER   IPAR(*)
C        REAL      X(NVAR)
C  C
C        FJAC(1,1)=1.0
C        FJAC(1,2)=(-3.0*X(2)+10.0)*X(2)-2.0
C        FJAC(1,3)=34.0
C  C
C        FJAC(2,1)=1.0
C        FJAC(2,2)=(3.0*X(2)+2.0)*X(2)-14.0
C        FJAC(2,3)=10.0
C  C
C        RETURN
C        END
C
C
C  Here is the output from the run of this sample problem
C
C
C  PCPRB1:
C
C  PITCON sample program.
C  Freudenstein-Roth function
C  2 equations, 3 variables.
C
C  Step  Type of Point     X(1)         X(2)         X(3)
C
C    0  Start point      15.0000      -2.00000      0.000000E+00
C
C  PITCON 6.1
C  University of Pittsburgh Continuation Code
C
C  Last modification on 25 February 1991
C
C    1  Corrected        15.0000      -2.00000      0.000000E+00
C    2  Continuation     14.7105      -1.94205      0.653814E-01
C    3  Continuation     14.2846      -1.72915      0.268742
C    4  Limit point      14.2831      -1.74138      0.258578
C    5  Continuation     16.9061      -1.20941      0.546845
C    6  Continuation     24.9179     -0.599064      0.555136
C    7  Continuation     44.8783      0.487670      0.595261E-01
C    8  Continuation     60.0889       1.57585     -0.542365
C    9  Continuation    -11.1752       4.23516       1.55667
C   10  Target point     5.00000       4.00000       1.00000
C
C  We have reached the point we wanted.
C  The code may stop now.
C  FORTRAN STOP
C
      DOUBLE PRECISION EIGHT
      DOUBLE PRECISION HALF
      DOUBLE PRECISION HUNDRD
      DOUBLE PRECISION ONE
      DOUBLE PRECISION THREE
      DOUBLE PRECISION TWO
      DOUBLE PRECISION ZERO
C
      PARAMETER (EIGHT=8.0)
      PARAMETER (HALF=0.5)
      PARAMETER (HUNDRD=100.0)
      PARAMETER (ONE=1.0)
      PARAMETER (THREE=3.0)
      PARAMETER (TWO=2.0)
      PARAMETER (ZERO=0.0)
C
      EXTERNAL  COQUAL
      EXTERNAL  CORECT
      EXTERNAL  DF
      EXTERNAL  FX
      EXTERNAL  IDAMAX
      EXTERNAL  REPS
      EXTERNAL  ROOT
      EXTERNAL  DDOT
      EXTERNAL  SETSTP
      EXTERNAL  SLNAME
      EXTERNAL  DNRM2
      EXTERNAL  DSCAL
      EXTERNAL  START
      EXTERNAL  TANGNT
C
      INTRINSIC ABS
      INTRINSIC ATAN2
      INTRINSIC MAX
      INTRINSIC SIGN
      INTRINSIC SQRT
C
      INTEGER   LIW
      INTEGER   LRW
      INTEGER   NVAR
C
      DOUBLE PRECISION A
      DOUBLE PRECISION ATCIPC
      DOUBLE PRECISION ATCJPC
      DOUBLE PRECISION B
      DOUBLE PRECISION DETS
      DOUBLE PRECISION DIRLPC
      DOUBLE PRECISION FA
      DOUBLE PRECISION FB
      DOUBLE PRECISION FPAR(*)
      DOUBLE PRECISION HTAN
      INTEGER   I
      INTEGER   ICALL
      INTEGER   ICRIT
      INTEGER   IERROR
      INTEGER   IFLAG
      INTEGER   IHOLD
      INTEGER   IMITL
      INTEGER   IPAR(*)
      INTEGER   IPC
      INTEGER   IPL
      INTEGER   IDAMAX
      INTEGER   IT
      INTEGER   ITMP
      INTEGER   IWORK(LIW)
      INTEGER   IWRITE
      INTEGER   JOB
      INTEGER   JPC
      INTEGER   LIM
      INTEGER   LOUNIT
      INTEGER   LPC
      INTEGER   LTC
      INTEGER   LTC1
      INTEGER   LWK
      INTEGER   LWK1
      INTEGER   LXC
      INTEGER   LXC1
      INTEGER   LXF
      INTEGER   LXF1
      INTEGER   MLSTEP
      INTEGER   MODSAV
c      DOUBLE PRECISION REPS
      DOUBLE PRECISION RWORK(LRW)
      DOUBLE PRECISION SCIPL
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION SKALE
      DOUBLE PRECISION SN
      DOUBLE PRECISION SNL
      DOUBLE PRECISION DNRM2
      DOUBLE PRECISION STEPX
      DOUBLE PRECISION TCIPC
      DOUBLE PRECISION TCIPL
      DOUBLE PRECISION TCOS
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION TLIPC
      DOUBLE PRECISION TMAX
      DOUBLE PRECISION TMAX2
      DOUBLE PRECISION TMP
      DOUBLE PRECISION TMP1
      DOUBLE PRECISION TSIN
      DOUBLE PRECISION TSN
      DOUBLE PRECISION XABS
      DOUBLE PRECISION XDIF
      DOUBLE PRECISION XLOW
      DOUBLE PRECISION XR(NVAR)
      DOUBLE PRECISION XUP
C
      SAVE      ICALL
C
      DATA      ICALL /0/
C
C     Definition of stream control variables 256 appears in gsl_findroot.c
      INTEGER       LSTDOUT
      Character*256 STDOUT

      IF (LSTDOUT .gt. 1) THEN
          OPEN ( IWORK(8), FILE=STDOUT, STATUS='OLD' )
      ELSE
          IWORK(8) = 0
      ENDIF
C     End of it. KMK


      ICALL=ICALL+1
C
C**********************************************************************
C
C  1.  Preparations.
C
C  Check entries of IWORK and RWORK, compute some constants.
C  For negative values of IWORK(1), compare user and finite difference
C    jacobian and return.
C  For IWORK(1)=0, set up starting data, check starting point and return.
C  For positive IWORK(1), prepare to compute next point.
C
C**********************************************************************
C
      IERROR=0
      IF(IWORK(8).LT.1.OR.IWORK(8).GT.99)IWORK(8)=6
      LOUNIT=IWORK(8)
      IF(IWORK(7).LT.0.OR.IWORK(7).GT.3)IWORK(7)=1
      IWRITE=IWORK(7)
C
C  Initialization for first call
C
      IF(IWORK(1).EQ.0.OR.ICALL.EQ.1)THEN
        IF(IWRITE.GE.1)THEN
          WRITE(LOUNIT,*)' '
          WRITE(LOUNIT,*)'PITCON 6.1'
          WRITE(LOUNIT,*)'University of Pittsburgh Continuation Code'
          WRITE(LOUNIT,*)' '
          WRITE(LOUNIT,*)'Last modification on 27 February 1991'
          WRITE(LOUNIT,*)' '
          ENDIF
        DO 10 I=11,17
          RWORK(I)=ZERO
10        CONTINUE
        DO 20 I=21,29
          RWORK(I)=ZERO
20        CONTINUE
        IWORK(10)=0
        IWORK(11)=0
        IWORK(12)=0
        DO 30 I=18,28
          IWORK(I)=0
30        CONTINUE
C
C  This is the only place REPS is called.  It might be replaced by
C  a call to the SLATEC machine parameter function R1MACH:
C
        RWORK(8)=D1MACH(4)
C
C  or by simply setting RWORK(8) to the "machine epsilon", usually
C  taken as the smallest power of 2, EPMACH, such that (1.0+EPMACH) is
C  greater than 1.0, but (1.0+(EPMACH/2.0)) equals 1.0.
C
c        RWORK(8)=REPS()
        IF(IWRITE.GE.2)WRITE(LOUNIT,919)RWORK(8)
        TCOS=SQRT(ONE-RWORK(8))
        TSIN=SQRT(RWORK(8))
        RWORK(10)=TWO*ATAN2(TSIN,TCOS)
        IF(NVAR.LE.1)THEN
          IERROR=1
          IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *    'PITCON - Number of variables must be at least 2.'
          GO TO 900
          ENDIF
        ENDIF
C
      LXC=29
      LXC1=LXC+1
      LXF=LXC+NVAR
      LXF1=LXF+1
      LTC=LXF+NVAR
      LTC1=LTC+1
      LWK=LTC+NVAR
      LWK1=LWK+1
      IWORK(13)=30
      IWORK(14)=LIW
      IWORK(15)=LWK+NVAR+1
      IWORK(16)=LRW
      IF(LIW.LT.IWORK(13))THEN
        WRITE(LOUNIT,*)'PITCON - Input value of LIW=',LIW
        WRITE(LOUNIT,*)'         Minimal acceptable value=',IWORK(13)
        IERROR=1
        GO TO 900
      ELSEIF(LRW.LT.IWORK(15))THEN
        WRITE(LOUNIT,*)'PITCON - Input value of LRW=',LRW
        WRITE(LOUNIT,*)'         Minimal acceptable value=',IWORK(15)
        IERROR=1
        GO TO 900
        ENDIF
C
C  Check entries of IWORK
C
      IF(IWORK(2).LT.1.OR.IWORK(2).GT.NVAR)IWORK(2)=NVAR
      IPC=IWORK(2)
      IF(IWORK(3).NE.1)IWORK(3)=0
      IF(IWORK(4).LT.0.OR.IWORK(4).GT.2)IWORK(4)=0
      IF(IWORK(5).LT.1.OR.IWORK(5).GT.NVAR)IWORK(5)=0
      IT=IWORK(5)
      IF(IWORK(6).LT.1.OR.IWORK(6).GT.NVAR)IWORK(6)=0
      LIM=IWORK(6)
      IF(IWORK(9).LT.0.OR.IWORK(9).GT.2)IWORK(9)=0
      IF(IWORK(17).LT.1)IWORK(17)=10
C
C  Check entries of RWORK
C
      IF(RWORK(1).LE.ZERO)RWORK(1)=SQRT(RWORK(8))
      IF(RWORK(2).LE.ZERO)RWORK(2)=SQRT(RWORK(8))
      IF(RWORK(3).LT.SQRT(RWORK(8)))RWORK(3)=SQRT(RWORK(8))
      IF(RWORK(4).LT.RWORK(3))THEN
        TEMP=NVAR
        RWORK(4)=MAX(RWORK(3),SQRT(TEMP))
        ENDIF
      IF(RWORK(5).LT.ZERO)THEN
        RWORK(5)=-RWORK(5)
        RWORK(6)=-RWORK(6)
        ENDIF
      IF(RWORK(5).LT.RWORK(3).OR.RWORK(5).GT.RWORK(4))THEN
        RWORK(5)=HALF*(RWORK(3)+RWORK(4))
        ENDIF
      HTAN=RWORK(5)
      IF(RWORK(6).NE.(-ONE))RWORK(6)=ONE
      IF(RWORK(20).LT.ONE.OR.RWORK(20).GT.HUNDRD)RWORK(20)=THREE
C
C  Check user Jacobian routine versus finite difference calculation.
C
      IF(IWORK(1).LT.0)THEN
        JOB=3
        CALL SLNAME(DETS,FX,DF,FPAR,IERROR,IPC,IPAR,IWORK,LIW,
     1  JOB,NVAR,RWORK,LRW,XR,RWORK(LWK1))
        GO TO 900
        ENDIF
C
C  Starting point check
C
      IF(IWORK(1).EQ.0)THEN
        CALL START(DF,FPAR,FX,IERROR,IPAR,IPC,IWORK,IWRITE,LIW,
     *  LOUNIT,LRW,NVAR,RWORK,RWORK(LTC1),RWORK(LWK1),RWORK(LXC1),
     *  RWORK(LXF1),XR,SLNAME)
        GO TO 900
        ENDIF
C
C***********************************************************************
C
C  2.  Target point
C
C  If (IT.NE.0) target points are sought.  Check to see if target component
C  IT has value XIT lying between XC(IT) and XF(IT).  If so, get linearly
C  interpolated starting point, and use Newton's method to get target point.
C
C***********************************************************************
C
      IF(IT.EQ.0.OR.IWORK(10).LE.1)GO TO 300
      IF( (IWORK(1).EQ.3) .AND. (RWORK(7).EQ.RWORK(28)) .AND.
     *    (IT.EQ.IWORK(11)) ) GO TO 300
      MODSAV=IWORK(4)
210   CONTINUE
      XLOW=RWORK(LXC+IT)
      XUP=RWORK(LXF+IT)
      IF(RWORK(7).LT.XLOW.AND.RWORK(7).LT.XUP)GO TO 300
      IF(RWORK(7).GT.XLOW.AND.RWORK(7).GT.XUP)GO TO 300
      IF(IWRITE.GE.2)WRITE(LOUNIT,*)
     *'PITCON - Attempt correction of approximate target point.'
C
C  Approximate the target point using the bracketing solutions.
C
      IF(XLOW.NE.XUP)THEN
        SKALE=(RWORK(7)-XLOW)/(XUP-XLOW)
      ELSE
        SKALE=ONE
        ENDIF
      CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
      CALL DSCAL(NVAR,SKALE,XR,1)
      SKALE=ONE-SKALE
      CALL DAXPY(NVAR,SKALE,RWORK(LXC1),1,XR,1)
      XR(IT)=RWORK(7)
C
C  Call CORECT to compute the exact target point, holding index IT fixed.
C
      ICRIT=0
      CALL CORECT(DF,FPAR,FX,IERROR,IT,IPAR,IWORK,
     1 NVAR,RWORK,TMP,RWORK(LWK1),XR,LRW,LIW,ICRIT,SLNAME)
      IWORK(24)=IWORK(24)+IWORK(28)
      IF(IERROR.NE.0.AND.IWORK(4).GT.0)THEN
        IERROR=0
        IWORK(4)=IWORK(4)-1
        IF(IWRITE.GE.1)WRITE(LOUNIT,1080)IWORK(4)
        GO TO 210
        ENDIF
      IWORK(4)=MODSAV
      IF(IERROR.NE.0)THEN
        WRITE(LOUNIT,*)'PITCON - Target point calculation failed.'
        GO TO 900
        ENDIF
C
C  Record the values of IT and XIT, and compute the arclength to the target
C  point.
C
      IWORK(1)=3
      IWORK(11)=IT
      RWORK(28)=RWORK(7)
      SKALE=-ONE
      CALL DCOPY(NVAR,XR,1,RWORK(LWK1),1)
      CALL DAXPY(NVAR,SKALE,RWORK(LXC1),1,RWORK(LWK1),1)
      RWORK(14)=RWORK(12)+DNRM2(NVAR,RWORK(LWK1),1)
      IWORK(27)=IWORK(27)+1
      GO TO 900
C
C***********************************************************************
C
C  3.  Tangent and local continuation parameter calculation.
C
C  Unless the tangent and limit point calculations were already performed,
C  (because the loop was interrupted for a limit point calculation),
C  set up and solve the equation for the tangent vector.
C
C  Force the tangent vector to be of unit length, and try to preserve
C  the "sense" or "direction" of the curve by forcing the IPL-th component
C  of the tangent vector to agree in sign with the IPL-th component of the
C  previous secant vector.  (On the first step, however, we have to use
C  the user's input direction to choose the sign).
C
C  Set the local continuation parameter IPC as the index of the
C  component of greatest magnitude in the tangent vector, unless a
C  limit point appears to be coming in that direction, and another
C  choice is available.
C
C***********************************************************************
C
  300 CONTINUE
      IF(IWORK(10).EQ.3) GO TO 600
C
C  STORE OLD TANGENT IN WORK ARRAY RWORK(LWK1)-RWORK(LWK+NVAR), COMPUTE
C  NEW TANGENT FOR CURRENT POINT XF IN RWORK(LXF1)-RWORK(LXF+NVAR).
C  THE TANGENT IS CALLED TC AND STORED IN RWORK(LTC1)-RWORK(LTC+NVAR).
C  NOTE THAT FOR A NORMAL RETURN THE CURRENT POINT XF WILL BE STORED
C  AS XC IN RWORK(LXC1)-RWORK(LXC+NVAR) AND HENCE THAT TC WILL THEN
C  CORRESPOND TO THE POINT XC
C
      IPL=IPC
      CALL DCOPY(NVAR,RWORK(LTC1),1,RWORK(LWK1),1)
      RWORK(29)=RWORK(17)
      CALL TANGNT(RWORK(17),FX,DF,FPAR,IERROR,IPC,IPAR,IWORK,NVAR,
     1 RWORK,RWORK(LTC1),RWORK(LXF1),LIW,LRW,SLNAME)
      IF(IERROR.NE.0)THEN
        IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *  'PITCON - Tangent calculation failed.'
        GO TO 900
        ENDIF
C
C  FIND LARGEST AND SECOND LARGEST COMPONENTS OF TANGENT
C
      IPC=IDAMAX(NVAR,RWORK(LTC1),1)
      TMAX=RWORK(LTC+IPC)
      RWORK(LTC+IPC)=ZERO
      JPC=IDAMAX(NVAR,RWORK(LTC1),1)
      RWORK(LTC+IPC)=TMAX
      TMAX=ABS(TMAX)
      TMAX2=ABS(RWORK(LTC+JPC))
      IF(JPC.EQ.0)JPC=IPC
      IF(IWORK(3).EQ.1)JPC=IPC
C
C  ADJUST SIGN OF THE TANGENT VECTOR.
C  COMPARE THE SIGN OF THE COMPONENT TC(IPL) WITH THE
C  SIGN OF XF(IPL)-XC(IPL). IF A TARGET OR LIMIT POINT
C  HAS BEEN COMPUTED BETWEEN XC AND XF THEN USE XF(IPL)-XR(IPL)
C  INSTEAD. AT THE FIRST STEP COMPARE WITH THE USER INPUT DIRECTION.
C  IF THESE SIGNS DIFFER, CHANGE THE SIGN OF THE TANGENT VECTOR
C  AND THE SIGN OF THE DETERMINANT.
C
      SCIPL=RWORK(6)
      TCIPL=RWORK(LTC+IPL)
      IF(IWORK(10).GT.1)THEN
        SCIPL=RWORK(LXF+IPL)
        TMP=RWORK(LXC+IPL)
        IF(IWORK(1).EQ.3.OR.IWORK(1).EQ.4)TMP=XR(IPL)
        SCIPL=SCIPL-TMP
        ENDIF
      IF(SIGN(ONE,TCIPL).NE.SIGN(ONE,SCIPL))THEN
        SKALE=-ONE
        CALL DSCAL(NVAR,SKALE,RWORK(LTC1),1)
        RWORK(17)=-RWORK(17)
        ENDIF
C
C  Unless we are computing a starting point, record the new state.
C
      IF(IWORK(10).GT.1)IWORK(10)=3
C
C  THE LOCATION OF THE LARGEST COMPONENT OF THE TANGENT VECTOR
C  WILL BE USED FOR THE LOCAL PARAMETERIZATION OF THE CURVE UNLESS
C  A LIMIT POINT IN IPC APPEARS TO BE COMING.
C  TO CHECK THIS, WE COMPARE TCIPC=TC(IPC) AND THE
C  SECOND LARGEST COMPONENT TCJPC=TC(JPC).  IF TCJPC IS NO LESS
C  THAN .1 OF TCIPC, AND TC(JPC) IS LARGER THAN TL(JPC),
C  WHEREAS TC(IPC) IS LESS THAN TL(IPC), WE WILL RESET THE
C  LOCAL PARAMETER IPC TO JPC.
C  BUT IF NOT ALLOWED TO SWITCH PARAMETERS, IGNORE THE CHECK
C
      IF(IWORK(3).NE.1)THEN
        ATCIPC=ABS(RWORK(LTC+IPC))
        ATCJPC=ABS(RWORK(LTC+JPC))
        IF(JPC.NE.IPC.AND.IWORK(10).GT.1)THEN
          TLIPC=RWORK(LWK+IPC)
          TCIPC=RWORK(LTC+IPC)
          TMP=ABS(RWORK(LWK+JPC))
          IF((SIGN(ONE,TCIPC).EQ.SIGN(ONE,TLIPC))
     1      .AND.(ATCIPC.LT.ABS(TLIPC))
     2      .AND.(ATCJPC.GE.MAX(0.1*ATCIPC,TMP)))THEN
            IF(IWRITE.GE.3)WRITE(LOUNIT,350)IPC
            ITMP=IPC
            IPC=JPC
            JPC=ITMP
            ENDIF
          ENDIF
        ENDIF
      IF(IWRITE.GE.3.AND.IWORK(3).EQ.0)THEN
        WRITE(LOUNIT,360)IPC
        WRITE(LOUNIT,370)JPC
        ENDIF
C
C  RECORD VALUES OF THE COMPONENT TCIPC=TC(IPC) OF
C  NEW TANGENT VECTOR, AS WELL AS OF THE CORRESPONDING COMPONENT
C  TLIPC=WK(IPC) OF THE OLD TANGENT VECTOR.
C  SET SIGN OF THE DETERMINANT.
C  FOR USE IN THE LIMIT POINT COMPUTATION RECORD THE
C  COMPONENT TC(LIM) OF THE NEW TANGENT VECTOR
C
      IWORK(2)=IPC
      IWORK(12)=JPC
      RWORK(24)=RWORK(LTC+IPC)
      RWORK(25)=RWORK(LWK+IPC)
      RWORK(6)=RWORK(17)
      IF(LIM.GT.0)THEN
        RWORK(27)=RWORK(26)
        RWORK(26)=RWORK(LTC+LIM)
        IF(IWRITE.GE.2)WRITE(LOUNIT,1020)RWORK(26)
        ENDIF
C
C  COMPUTE ALPHLC, THE ANGLE BETWEEN OLD TANGENT AND TANGENT AT XC
C
      IF(IWORK(10).LE.1) GO TO 600
      TCOS=DDOT(NVAR,RWORK(LWK1),1,RWORK(LTC1),1)
      IF(TCOS.GT.ONE)TCOS=ONE
      IF(TCOS.LT.(-ONE))TCOS=-ONE
      TSIN=SQRT(ONE-TCOS*TCOS)
      RWORK(11)=ATAN2(TSIN,TCOS)
      IF((IWORK(10).LE.1).OR.(RWORK(17).EQ.RWORK(29)))GO TO 400
      IF(IWRITE.GE.2)WRITE(LOUNIT,*)
     *'PITCON - Possible bifurcation point detected.'
C
C***********************************************************************
C
C  4.  Limit point check.
C
C  Skip this section if LIM=0.
C
C  Otherwise, user has requested a search for limit points in a given
C  index by setting LIM to a nonzero value.
C
C  Compare LIM-th components of previous and current tangent vectors.
C  If a sign difference occurs, we assume a limit point has been passed.
C  Attempt to compute a point XR between the previous and current points,
C  for which the LIM-th component of the tangent is zero.
C
C  This search will be guided by a rootfinder.  The scalar function
C  to be zeroed out is the LIM-th tangent vector component.
C
C***********************************************************************
C
  400 CONTINUE
      IF((LIM.LE.0).OR.
     *   (IWORK(10).LE.1).OR.
     *   (SIGN(ONE,RWORK(26)).EQ.SIGN(ONE,RWORK(27)))) GO TO 600
      IF(IWRITE.GE.2)WRITE(LOUNIT,*)
     *'PITCON - Attempt correction of approximate limit point.'
      MODSAV=IWORK(4)
410   CONTINUE
C
C  SEE IF EITHER ENDPOINT CAN BE ACCEPTED AS LIMIT POINT.
C  IF XC IS LIMIT POINT, WORK ALREADY CONTAINS TANGENT AT XC.
C
      IF(ABS(RWORK(27)).LE.HALF*RWORK(1)) THEN
        CALL DCOPY(NVAR,RWORK(LXC1),1,XR,1)
        GO TO 490
        ENDIF
      IF(ABS(RWORK(26)).LE.HALF*RWORK(1)) THEN
        CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
        CALL DCOPY(NVAR,RWORK(LTC1),1,RWORK(LWK1),1)
        GO TO 490
        ENDIF
C
C  If interval is extremely small, simply assign one
C  endpoint of the interval as the answer.
C
      XDIF=ABS(RWORK(LXC+LIM)-RWORK(LXF+LIM))
      XABS=MAX(ABS(RWORK(LXC+LIM)),ABS(RWORK(LXF+LIM)))
      IF(XDIF.LE.EIGHT*RWORK(8)*(ONE+XABS)) THEN
        IF(ABS(RWORK(27)).GT.ABS(RWORK(26))) THEN
          CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
          CALL DCOPY(NVAR,RWORK(LTC1),1,RWORK(LWK1),1)
        ELSE
          CALL DCOPY(NVAR,RWORK(LXC1),1,XR,1)
          ENDIF
        GO TO 490
        ENDIF
C
C  BEGIN ROOT-FINDING ITERATION WITH INTERVAL (0,1) AND
C  FUNCTION VALUES TLLIM, TCLIM.
C
      A=ZERO
      FA=RWORK(27)
      B=ONE
      FB=RWORK(26)
C
C  FIND LPC, INDEX OF MAXIMUM ENTRY OF SECANT (EXCEPT THAT
C  LPC MUST NOT EQUAL LIM) AND SAVE SIGN IN DIRLPC
C  SO THAT NEW TANGENTS MAY BE PROPERLY SIGNED.
C
      TEMP=RWORK(LXC+LIM)
      RWORK(LXC+LIM)=RWORK(LXF+LIM)
      CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
      SKALE=-ONE
      CALL DAXPY(NVAR,SKALE,RWORK(LXC1),1,XR,1)
      LPC=IDAMAX(NVAR,XR,1)
      RWORK(LXC+LIM)=TEMP
      DIRLPC=SIGN(ONE,RWORK(LXF+LPC)-RWORK(LXC+LPC))
C
C  SET FIRST APPROXIMATION TO ROOT TO WHICHEVER ENDPOINT
C  HAS SMALLEST LIM-TH COMPONENT OF TANGENT
C
      IF(ABS(RWORK(26)).GE.ABS(RWORK(27))) THEN
        SN=ZERO
        TSN=RWORK(27)
        CALL DCOPY(NVAR,RWORK(LXC1),1,XR,1)
      ELSE
        SN=ONE
        TSN=RWORK(26)
        CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
        ENDIF
C
C  CALL ROOTFINDER FOR APPROXIMATE ROOT SN.  USE LINEAR COMBINATION
C  OF PREVIOUS ROOT SNL, AND ONE OF 0.0 AND 1.0 TO GET A STARTING
C  POINT FOR CORRECTION.  RETURN TO CURVE ON LINE X(LPC)=CONSTANT,
C  COMPUTE TANGENT THERE, AND SET FUNCTION VALUE TO TANGENT(LIM)
C
      MLSTEP=25
      IMITL=0
      IF(IWRITE.GE.2)THEN
        TMP=ZERO
        WRITE(LOUNIT,1100)TMP,RWORK(27)
        TMP=ONE
        WRITE(LOUNIT,1100)TMP,RWORK(26)
        ENDIF
441   CONTINUE
      SNL=SN
      CALL ROOT(A,FA,B,FB,SN,TSN,IMITL,IFLAG,IERROR,RWORK(8))
      IWORK(23)=IWORK(23)+1
      IF(IERROR.NE.0)THEN
        IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *  'PITCON - Bad interval for rootfinder'
        GO TO 900
        ENDIF
      IF(IFLAG.EQ.(-1).OR.IFLAG.EQ.0)GO TO 490
C
C  FIND WHETHER SN LIES IN (0.0,SNL) OR (SNL,1.0).
C  USE APPROPRIATE LINEAR COMBINATION TO GET STARTING POINT.
C
      IF(SN.LE.SNL) THEN
C
C  SET X(SN)=(SNL-SN)/(SNL-0.0)*X(0.0)+(SN-0.0)/(SNL-0.0)*X(SNL)
C
        IF(SNL.LE.ZERO)THEN
          SKALE=ZERO
        ELSE
          SKALE=SN/SNL
          ENDIF
        IF(SKALE.LE.ZERO)SKALE=ZERO
        IF(SKALE.GE.ONE)SKALE=ONE
        CALL DSCAL(NVAR,SKALE,XR,1)
        SKALE=ONE-SKALE
        CALL DAXPY(NVAR,SKALE,RWORK(LXC1),1,XR,1)
      ELSE
C
C  SET X(SN)=(SN-SNL)/(1.0-SNL)*X(1.0)+(1.0-SN)/(1.0-SNL)*X(SNL)
C
        IF(SNL.GE.ONE)THEN
          SKALE=ZERO
        ELSE
          SKALE=(ONE-SN)/(ONE-SNL)
          ENDIF
        IF(SKALE.LE.ZERO)SKALE=ZERO
        IF(SKALE.GE.ONE)SKALE=ONE
        CALL DSCAL(NVAR,SKALE,XR,1)
        SKALE=ONE-SKALE
        CALL DAXPY(NVAR,SKALE,RWORK(LXF1),1,XR,1)
        ENDIF
C
582   CONTINUE
      ICRIT=0
      CALL CORECT(DF,FPAR,FX,IERROR,LPC,IPAR,IWORK,
     1 NVAR,RWORK,TMP,RWORK(LWK1),XR,LRW,LIW,ICRIT,SLNAME)
      IF(IERROR.NE.0.AND.IWORK(4).GT.0)THEN
        IWORK(4)=IWORK(4)-1
        CALL DCOPY(NVAR,RWORK(LXC1),1,XR,1)
        SKALE=ONE-SN
        CALL DSCAL(NVAR,SKALE,XR,1)
        SKALE=SN
        CALL DAXPY(NVAR,SKALE,RWORK(LXF1),1,XR,1)
        IF(IWRITE.GE.1)WRITE(LOUNIT,1090)IWORK(4)
        GO TO 582
        ENDIF
      IWORK(4)=MODSAV
      IF(IERROR.NE.0)THEN
        WRITE(LOUNIT,*)
     *  'PITCON - Corrector failed on limit point.'
        GO TO 900
        ENDIF
      CALL TANGNT(TEMP,FX,DF,FPAR,IERROR,LPC,IPAR,IWORK,NVAR,
     1 RWORK,RWORK(LWK1),XR,LIW,LRW,SLNAME)
      IF(IERROR.NE.0)THEN
        IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *  'PITCON - Tangent error in limit point calculation.'
        GO TO 900
        ENDIF
C
C  Adjust the sign of the tangent vector so the LPC-th component
C  has the same sign as the LPC-th component of the secant.
C
      IF(DIRLPC.NE.SIGN(ONE,RWORK(LWK+LPC))) THEN
        SKALE=-ONE
        CALL DSCAL(NVAR,SKALE,RWORK(LWK1),1)
        ENDIF
C
C  SEE IF WE CAN ACCEPT THE NEW POINT BECAUSE TANGENT(LIM) IS SMALL
C  OR MUST STOP BECAUSE 25 ITERATIONS TAKEN,
C  OR IF WE MUST GO ON.
C
      TSN=RWORK(LWK+LIM)
      IF(IWRITE.GE.2)WRITE(LOUNIT,1100)SN,TSN
      IF(ABS(TSN).LE.RWORK(1))GO TO 490
      IF(IMITL.LT.MLSTEP)GO TO 441
C
C  Limit point iteration has failed.  We set an error flag,
C  but we return the partial results of the abortive computation
C  anyway.
C
      IERROR=8
      IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *'PITCON - Too many steps in limit point calculation.'
C
C  Limit point iteration is over.
C  Compute and store information.
C
  490 CONTINUE
      CALL DCOPY(NVAR,XR,1,RWORK(LWK1),1)
      SKALE=-ONE
      CALL DAXPY(NVAR,SKALE,RWORK(LXC1),1,RWORK(LWK1),1)
      RWORK(14)=RWORK(12)+DNRM2(NVAR,RWORK(LWK1),1)
      IWORK(27)=IWORK(27)+1
      IWORK(1)=4
      GO TO 900
C
C***********************************************************************
C
C  5.  Compute next predictor step length, HTAN.
C
C***********************************************************************
C
  600 CONTINUE
      IF(IWORK(10).GT.1)CALL SETSTP(HTAN,IWORK,LIW,LRW,RWORK)
C
C***********************************************************************
C
C  6.  Continuation step
C
C  Our current data is the current point XC, its tangent vector TC, and
C  a steplength HTAN.  We predict the location of the next point on the
C  curve using the Euler approximation XR=XC+HTAN*TC.
C
C  Newton iteration is applied to this point, to force it to lie on the
C  curve.  In order to make the system square, an augmenting equation
C  is added to the system, specifying that XR(IPC)=XC(IPC)+HTAN*TC(IPC).
C  (The right hand side is a constant.)
C
C  If the Newton correction process fails, the stepsize is reduced and
C  prediction and correction retried.  Failure will most likely be
C  signaled by repeated step reductions, until the minimum allowable
C  stepsize is reached.  If this occurs, PITCON has failed, and cannot
C  proceed along the curve any more.
C
C***********************************************************************
C
  700 CONTINUE
      IWORK(18)=0
      IPC=IWORK(2)
      JPC=IWORK(12)
      IHOLD=IPC
      MODSAV=IWORK(4)
710   CONTINUE
      IF(IWRITE.GE.2)WRITE(LOUNIT,715)HTAN
      CALL DCOPY(NVAR,RWORK(LXF1),1,XR,1)
      CALL DAXPY(NVAR,HTAN,RWORK(LTC1),1,XR,1)
  740 CONTINUE
      IF(IWRITE.GE.2.AND.IWORK(3).EQ.0)WRITE(LOUNIT,725)IHOLD
      ICRIT=0
      CALL CORECT(DF,FPAR,FX,IERROR,IHOLD,IPAR,IWORK,
     1 NVAR,RWORK,STEPX,RWORK(LWK1),XR,LRW,LIW,ICRIT,SLNAME)
      IWORK(25)=IWORK(25)+IWORK(28)
      IF(IERROR.EQ.0)GO TO 800
C
C  Only VERY fatal errors should abort the correction process.
C  Abort this process only after kicking and screaming.
C
      IF(IERROR.EQ.2)THEN
        WRITE(LOUNIT,*)'PITCON - Fatal error during Newton correction'
        WRITE(LOUNIT,*)'         of a continutation point!'
        GO TO 900
        ENDIF
C
C  FOR CASES WHERE CORRECTOR FAILED TO CONVERGE,
C  AND THE VALUE OF IHOLD WAS NOT JPC,
C  AND THE USER HAS NOT FIXED THE CONTINUATION PARAMETER
C  KEEP THE SAME STEPSIZE BUT USE JPC INSTEAD OF IPC.
C
      IF(IWORK(3).NE.1.AND.IHOLD.NE.JPC.AND.JPC.NE.0)THEN
        IHOLD=JPC
        IF(IERROR.EQ.5)GO TO 740
        GO TO 710
        ENDIF
C
C  IF JPC FAILS AS IHOLD VALUE, RETURN TO IPC.
C
      IHOLD=IPC
C
C  IF WE ARE USING SOME FORM OF MODIFIED NEWTON, BUT
C  WE HAVE REACHED MINIMUM STEPSIZE OR
C  HAVE HAD TWO FAILURES IN A ROW ON THIS STEP,
C  USE A BETTER NEWTON METHOD
C
      IF(IWORK(4).GT.0.AND.
     *  (HTAN.LT.RWORK(20)*RWORK(3).OR.IWORK(18).GE.2) )THEN
        IWORK(4)=IWORK(4)-1
        WRITE(LOUNIT,1130)IWORK(4)
        IF(IERROR.NE.5)GO TO 710
        GO TO 740
        ENDIF
C
C  NO CONVERGENCE, TRY A SMALLER STEPSIZE
C
      IF(HTAN.LT.RWORK(20)*RWORK(3))THEN
        IERROR=4
        IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *  'PITCON - Stepsize fell below minimum.'
        GO TO 900
        ENDIF
      HTAN=HTAN/RWORK(20)
      IWORK(18)=IWORK(18)+1
      IF(IERROR.NE.5)GO TO 710
      SKALE=ONE/RWORK(20)
      CALL DSCAL(NVAR,SKALE,XR,1)
      SKALE=ONE-SKALE
      CALL DAXPY(NVAR,SKALE,RWORK(LXF1),1,XR,1)
      GO TO 740
C
C***********************************************************************
C
C  7.  Successful step.  Update information.
C
C***********************************************************************
C
  800 CONTINUE
      IWORK(1)=2
      IF(IWORK(3).NE.1)IWORK(2)=IHOLD
      IWORK(4)=MODSAV
      IWORK(10)=2
      IWORK(26)=IWORK(26)+IWORK(18)
      IWORK(27)=IWORK(27)+1
      RWORK(5)=HTAN
      RWORK(22)=RWORK(21)
      IF(IWRITE.GE.1.AND.IWORK(18).GT.0)WRITE(LOUNIT,755)IWORK(18)
      CALL DCOPY(NVAR,XR,1,RWORK(LWK1),1)
      SKALE=-ONE
      CALL DAXPY(NVAR,SKALE,RWORK(LXF1),1,RWORK(LWK1),1)
      RWORK(21)=DNRM2(NVAR,RWORK(LWK1),1)
      RWORK(12)=RWORK(13)
      RWORK(13)=RWORK(12)+RWORK(21)
      RWORK(14)=RWORK(13)
C
C  UPDATE VALUES OF COMPUTED POINTS.  STORE OLD POINT IN
C  RWORK(LXC1)-RWORK(LXC+NVAR) AND NEW POINT IN RWORK(LXF1)-
C  RWORK(LXF+NVAR). COMPUTE AND STORE CORDIS, THE MAXIMUM NORM
C  OF THE CORRECTOR DISTANCE.
C
      RWORK(15)=ZERO
      IF(IWORK(28).NE.0)THEN
        DO 810 I=1,NVAR
          TMP1=XR(I)-(RWORK(LXF+I)+HTAN*RWORK(LTC+I))
          IF(ABS(TMP1).GT.RWORK(15))RWORK(15)=ABS(TMP1)
 810      CONTINUE
        ENDIF
      CALL DCOPY(NVAR,RWORK(LXF1),1,RWORK(LXC1),1)
      CALL DCOPY(NVAR,XR,1,RWORK(LXF1),1)
C
C  COMPUTE CONVERGENCE QUALITY
C
      CALL COQUAL(STEPX,IWORK,LIW,RWORK,LRW)
C
C**********************************************************************
C
C  8.  All returns to calling program should exit here.
C
C**********************************************************************
C
  900 CONTINUE
      IF(IERROR.EQ.0)THEN
      ELSEIF(IERROR.EQ.1)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Storage or data error.'
      ELSEIF(IERROR.EQ.2)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)
     *  '         User-defined function or jacobian error.'
      ELSEIF(IERROR.EQ.3)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Solver failed.'
      ELSEIF(IERROR.EQ.4)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Iteration failed.'
      ELSEIF(IERROR.EQ.5)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Too many corrector steps.'
      ELSEIF(IERROR.EQ.6)THEN
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Zero tangent vector.'
      ELSEIF(IERROR.EQ.7)THEN
        WRITE(LOUNIT,*)'PITCON - Nonfatal error.'
        WRITE(LOUNIT,*)'         Root finder failed on limit point.'
      ELSEIF(IERROR.EQ.8)THEN
        WRITE(LOUNIT,*)'PITCON - Nonfatal error.'
        WRITE(LOUNIT,*)'         Too many limit point steps.'
      ELSE
        IERROR=9
        WRITE(LOUNIT,*)'PITCON - FATAL ERROR.'
        WRITE(LOUNIT,*)'         Unidentified error.'
      ENDIF
      IF (LSTDOUT .gt. 1)
     *  CLOSE ( IWORK(8) )
      RETURN
350   FORMAT(' PITCON - Expect limit point in index ',I5)
360   FORMAT(' PITCON - Continuation index: First choice',I5)
370   FORMAT(' PITCON -                    Second choice',I5)
715   FORMAT(' PITCON - Predictor using stepsize ',G14.6)
725   FORMAT(' PITCON - Corrector fixing index',I5)
755   FORMAT(' PITCON - Predictor stepsize was reduced ',I6,' times.')
919   FORMAT(' PITCON - Machine epsilon=',G14.6)
1020  FORMAT(' PITCON - Tangent vector has limit component =',G14.6)
1080  FORMAT(' PITCON - Retry target computation with IWORK(4)=',I6)
1090  FORMAT(' PITCON - Retry limit computation with IWORK(4)=',I6)
1100  FORMAT(' PITCON - For S=',G14.6,' TAN(X(S))(LIM)=',G14.6)
1130  FORMAT(' PITCON - Retrying step with IWORK(4)=',I6)
      END
