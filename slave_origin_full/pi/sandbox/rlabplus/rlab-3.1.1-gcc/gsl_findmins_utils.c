//
//
// findmins, without derivatives, utility functions
//
//
static double mins_gslrlab_f (const gsl_vector * x, void *params)
{
  Ent *rent = 0;
  double rval;

  MDPTR(xmdr) = (void *) x->data;

  if (pent)
    rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_f_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}

//
// jacobian function for the function being minimized
//
static void mins_gslrlab_df (const gsl_vector * x, void *params, gsl_vector * df)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(xmdr) = (void *) x->data;

  if (pent)
    rent = ent_call_rlab_script_2args(cmx_df_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_df_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != cmx_dimx)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (i=0; i<cmx_dimx; i++)
    gsl_vector_set (df, i, MdrV0 (retm, i));

  ent_Clean (rent);
  return;
}

static void func_fdf (const gsl_vector * x, void *params, double *f, gsl_vector * df)
{
  *f = mins_gslrlab_f (x, params);
  mins_gslrlab_df (x, params, df);
  return;
}

//
// using GSL SIMPLEX solver
//
//
MDR * mdr_findmins_simplex(MDR * ystart, MDR * h, char * outs, double abserr, double target, int maxiter, int istand)
{
  int j, i, cmx_dimx=SIZE(ystart), status, irep=0;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_vector *ss = gsl_vector_alloc (cmx_dimx);
  if(SIZE(h) > 0)
  {
    for (i = 0; i < cmx_dimx; i++)
      gsl_vector_set (ss, i, mdrV0_safe(h,i));
  }
  else
  {
    for (i = 0; i < cmx_dimx; i++)
      gsl_vector_set (ss, i, 1.0);
  }

  gsl_vector *x = gsl_vector_copy_mdr( ystart );
  gsl_multimin_function f =
  {
    &mins_gslrlab_f,
    cmx_dimx,
    NULL
  };
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, cmx_dimx);
  gsl_multimin_fminimizer_set (s, &f, x, ss);

  int louts = isvalidstring( outs );
  FILE *fptr=0;
  time_t t1, t2;
  double prev_val=-1, curr_val=-1, timer;

  gsl_set_error_handler_off();

  //
  // do not write any run-time messages
  //
  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr, "RLaB: using gsl '%s' minimizer in multidimensions.\n",
             gsl_multimin_fminimizer_name (s));
    fprintf (fptr, "RLaB: Stream %s: Messages from the solver follow.\n", outs);
  }
  else
    louts = 0;

  i = 0;
  do
  {
    i++;

    // do an iteration
    status = gsl_multimin_fminimizer_iterate (s);
    if (status)
      break;
    status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s), abserr);

    // update
    prev_val = curr_val;
    curr_val = s->fval;
    if (curr_val == prev_val)
      irep++;
    else
      irep=0;

    if (fptr)
    {
            // write text messages
      fprintf (fptr, "iteration = %i", i);
      for (j = 0; j < cmx_dimx; j++)
        fprintf (fptr, " : x[%i] = %g", j + 1, gsl_vector_get (s->x, j));
      fprintf (fptr, " : f = %g", s->fval);
      fprintf (fptr, " : size = %g", gsl_multimin_fminimizer_size (s));
      fprintf (fptr, "\n");
    }

          // did we reach the target?
    if (!isnand(target))
      if (target == s->fval && status == GSL_CONTINUE)
        status = GSL_SUCCESS;
  }
  while ((status == GSL_CONTINUE) && (i < maxiter) && (irep < istand));

  if (fptr)
  {
    fprintf (fptr,
             "RLaB: gsl multidimension minimizer '%s' reports",
             gsl_multimin_fminimizer_name (s));
    fprintf (fptr, " '%s' .\n", gsl_strerror (status) );
          // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr,
             "RLaB: minimization in multidimensions lasted %g sec.\n",
             timer
            );
    fflush (fptr);
    fclose (fptr);
  }

  MDR * w = copy_mdr_gsl_vector(s->x);

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);

  return w;
}

//
// LIBMJDPOWELL fortran wrapper
//
static int newuoa_f (int * idummy, double * xval, double * f)
{
  Ent *rent = 0;

  // x
  MDPTR(xmdr) = xval;

  if (pent)
    rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_f_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  *f = class_double (rent);
  ent_Clean (rent);
  return 1;
}

MDR * mdr_findmins_powell(MDR * ystart, MDR * a, MDR *b, char * outs, int npt,
                          double rhobeg, double rhoend, int maxiter)
{
  int i, j, k, cmx_dimx=SIZE(ystart), n_con=0;
  int louts = isvalidstring( outs );
  FILE *fptr=0;
  time_t t1, t2;
  double timer;

  MDR *bounds=0, *linop=0, *w=mdr_Float_BF(ystart);
  double * wk = GC_MALLOC(((npt+13)*(npt+cmx_dimx)+3*cmx_dimx*(cmx_dimx+3)/2+1)*sizeof(double));

  if (SIZE(a)>0)
  {
    if (MNC(a) == cmx_dimx)
    {
      if (MNC(b)==1)
      {
        bounds = mdr_Float_BF(b);
        linop = mdr_Float_BF(a);
      }
      else if (MNC(b)==2)
      {
        // count constraints. ignore infs
        for(i=0; i<MNR(b); i++)
          if ((mdr0(b,i,0)!= -create_inf()) && (mdr0_safe(b,i,0) < mdr0_safe(b,i,1)))
            n_con++;
        for(i=0; i<MNR(b); i++)
          if (mdr0(b,i,1)!= create_inf() && (mdr0_safe(b,i,0) < mdr0_safe(b,i,1)))
            n_con++;

        if (n_con>0)
        {
          linop = mdr_Create(n_con, SIZE(ystart));
          bounds = mdr_Create(n_con,1);
          mdr_Zero(bounds);
          mdr_Zero(linop);
          j=0;
          for(i=0; i<MNR(b); i++)
          {
            if ((mdr0(b,i,1)!= create_inf()) && (mdr0_safe(b,i,0) < mdr0_safe(b,i,1)))
            {
              MdrV0(bounds,j) = mdr0_safe(b,i,1);
              for (k=0; k<cmx_dimx; k++)
                Mdr0(linop,j,k)=mdr0(a,i,k);
              j++;
            }
          }
          for(i=0; i<MNR(b); i++)
          {
            if ((mdr0(b,i,0)!= -create_inf()) && (mdr0_safe(b,i,0) < mdr0_safe(b,i,1)))
            {
              MdrV0(bounds,j) = -mdr0_safe(b,i,0);
              for (k=0; k<cmx_dimx; k++)
                Mdr0(linop,j,k) = -mdr0(a,i,k);
              j++;
            }
          }
        }
      }
      mdr_Transpose_inplace(linop); // powell likes his transposed
    }
  }
  else if (SIZE(b)>0)
  {
    if (MNC(b) == 2)
    {
      if (MNR(b) == cmx_dimx)
        bounds = mdr_Float_BF(b);
      else
      {
        bounds = mdr_Create(cmx_dimx,2);
        for (i=0; i<cmx_dimx; i++)
        {
          Mdr0(bounds,i,0) = mdr0_safe(b,i,0);
          Mdr0(bounds,i,1) = mdr0_safe(b,i,1);
        }
      }
    }
  }

  //
  // run-time messages
  //
  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
    t1 = clock ();
  else
    louts = 0;

  if (n_con>0)
  {
    if (fptr)
    {
      fprintf (fptr, "RLaB: using LINCOA minimizer in multidimensions.\n");
      fprintf (fptr, "RLaB: Stream %s: Messages from the solver follow.\n", outs);
    }

    // findmin with constraint (b is one-column vector)
    //    linop * w <= b
    // or
    //    b[;1] <= linop * w <= b[;2] (b is two-column matrix)
    LINCOA (&cmx_dimx, &npt, &n_con, MDRPTR(linop), &cmx_dimx,
             MDRPTR(bounds), MDRPTR(w), &rhobeg, &rhoend, &louts, &maxiter, wk, newuoa_f, outs);
  }
  else if (SIZE(b)>0)
  {
    if (fptr)
    {
      fprintf (fptr, "RLaB: using BOBYQA minimizer in multidimensions.\n");
      fprintf (fptr, "RLaB: Stream %s: Messages from the solver follow.\n", outs);
    }

    // findmin with constraint
    //    b[;1] <= w <= b[;2]
    BOBYQA (&cmx_dimx, &npt, MDRPTR(w), &Mdr0(bounds,0,0), &Mdr0(bounds,0,1),
             &rhobeg, &rhoend, &louts, &maxiter, wk, newuoa_f, outs);
  }
  else
  {
    if (fptr)
    {
      fprintf (fptr, "RLaB: using NEWUOA minimizer in multidimensions.\n");
      fprintf (fptr, "RLaB: Stream %s: Messages from the solver follow.\n", outs);
    }

    // unconstrained optimization
    NEWUOA (&cmx_dimx, &npt, MDRPTR(w), &rhobeg, &rhoend, &louts, &maxiter, wk, newuoa_f, outs);
  }

  if (fptr)
  {
    // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: minimization in multidimensions lasted %g sec.\n", timer);
    fflush (fptr);
    fclose (fptr);
  }

  if (bounds)
    mdr_Destroy (bounds);
  if (linop)
    mdr_Destroy (linop);

  GC_FREE (wk);

  return w;
}

//
// LIBMJDPOWELL fortran wrapper
//
static int constr_f (int * idummy1, int * idummy2,
                     double * xval, double * f, double * con_f)
{
  Ent *rent = 0;
  int i;

  // x
  MDPTR(xmdr) = xval;

  // function value
  if (pent)
  {
    rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
  }
  else
  {
    rent = ent_call_rlab_script_1arg (cmx_f_name, xent);
  }
  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);
  if (f)
    *f = class_double (rent);
  ent_Clean (rent);

  // constraint value
  if (pent)
    rent = ent_call_rlab_script_2args(cmx_g_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_g_name, xent);
  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  MDR * c = ent_data(rent);
  if (cmx_dimg > 0)
  {
    if (cmx_dimg != SIZE(c))
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    for (i=0; i<cmx_dimg; i++)
      con_f[i] = mdrV0(c,i);
  }
  else
  {
    cmx_dimg = SIZE(c);
  }
  ent_Clean (rent);
  return 1;
}

MDR * mdr_findmins_powell_cobyla(MDR * ystart, char * outs, double rhobeg, double rhoend, int maxiter)
{
  int cmx_dimx=SIZE(ystart);
  int louts = isvalidstring( outs );
  FILE *fptr=0;
  time_t t1, t2;
  double timer;

  MDR *w=mdr_Float_BF(ystart);

  // call constr_f once to figure out the dimension of the constraint function 'g'
  // after the first call, 'cmx_dimg' contains dimension of 'g'
  if (cmx_dimg == -1)
  {
    constr_f (NULL, NULL, MDRPTR(w), NULL, NULL);
  }
  double * wk = GC_MALLOC((cmx_dimx*(3*cmx_dimx+2*cmx_dimg+11)+4*cmx_dimg+6)*sizeof(double));
  int    * ia = GC_MALLOC((cmx_dimg+1)*sizeof(int));

  //
  // run-time messages
  //
  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr, "RLaB: using COBYLA minimizer in multidimensions.\n");
    fprintf (fptr, "RLaB: Stream %s: Messages from the solver follow.\n", outs);
  }
  else
    louts = 0;

  // optimization with nonlinear bounds
  COBYLA (&cmx_dimx, &cmx_dimg, MDRPTR(w), &rhobeg, &rhoend, &louts, &maxiter, wk, ia,
           constr_f, outs);

  if (fptr)
  {
    // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: minimization in multidimensions lasted %g sec.\n", timer);
    fflush (fptr);
    fclose (fptr);
  }

  GC_FREE (wk);
  GC_FREE (ia);

  return w;
}


//
//
// findmins, with derivatives, utility functions
//
//
MDR * mdr_findmins_f_df_gsl(int imethod, MDR * ystart, double h, char * outs, double abserr, double tol,
                            double target, int maxiter)
{
  MDR *w=0;
  FILE *fptr=0;
  int i=0, j, louts=isvalidstring(outs), status;
  time_t t1, t2;
  double timer;

  gsl_set_error_handler_off();

  if (imethod<0 || imethod>3)
    return w;

  //
  // unconstrained minimization with the GSL
  //
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x = gsl_vector_copy_mdr( ystart );

  gsl_multimin_function_fdf f =
  {
    &mins_gslrlab_f,
    &mins_gslrlab_df,
    &func_fdf,
    cmx_dimx,
    NULL
  };

  switch (imethod)
  {
    case 0:
      T = gsl_multimin_fdfminimizer_conjugate_fr;
      break;
    case 1:
      T = gsl_multimin_fdfminimizer_conjugate_pr;
      break;
    case 2:
      T = gsl_multimin_fdfminimizer_vector_bfgs;
      break;
    case 3:
      T = gsl_multimin_fdfminimizer_steepest_descent;
      break;
    default:
      T = gsl_multimin_fdfminimizer_conjugate_fr;
  }
  s = gsl_multimin_fdfminimizer_alloc (T, cmx_dimx);
  gsl_multimin_fdfminimizer_set (s, &f, x, h, tol);

  // write run-time messages if any
  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr, "RLaB: using gsl '%s' minimizer in multidimensions.\n",
             gsl_multimin_fdfminimizer_name (s));
    fprintf (fptr,
             "RLaB: Stream %s: Messages from the solver follow.\n", outs
            );
  }
  else
    louts = 0;


  i = 0;
  do
  {
    i++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    if (status != GSL_SUCCESS)
      break;

    status = gsl_multimin_test_gradient (s->gradient, abserr);

    if (fptr)
    {
        // write run-time messages if any
      fprintf (fptr, "iteration = %i", i);
      for (j = 0; j < cmx_dimx; j++)
        fprintf (fptr, " : x[%i] = %g", j + 1, gsl_vector_get (s->x, j));
      fprintf (fptr, " : f = %g", s->f);
      fprintf (fptr, "\n");
    }

    if (!isnand(target))
      if (target == s->f && status == GSL_CONTINUE)
        status = GSL_SUCCESS;
  }
  while ((status == GSL_CONTINUE) && (i < maxiter));

  // write run-time messages if any
  if (fptr)
  {
    fprintf (fptr, "RLaB: gsl multidimension minimizer '%s' reports",
             gsl_multimin_fdfminimizer_name (s));
    fprintf (fptr, " '%s' .\n", gsl_strerror (status));

      // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: minimization in multidimensions lasted %g sec.\n", timer);
    fclose (fptr);
  }

  w = mdr_Float_BF ( ystart );
  for (i = 0; i < cmx_dimx; i++)
    MdrV0 (w, i) = gsl_vector_get (s->x, i);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return w;
}

//
// unconstrained optimization with CONMAX
//
MDR * mdr_findmins_f_df_conmax(int cmx_dimf, MDR * ystart, char * outs,
                               double encsm, int limsm, int nrkstep, double tolcon, int istep,
                              int maxiter)
{
  MDR *w=0;
  FILE *fptr=0;
  int i=0, louts=isvalidstring(outs);
  time_t t1, t2;
  double timer;

  int ioptn = 0, iptb, indb, liwrk, lwrk, iter, nparm;
  MDR *fun=0, *pttbl=0, *iwork=0, *work=0, *err=0;

  cmx_icntyp = mdi_Create(cmx_numgr, 1);
  for (i=0; i < cmx_dimf; i++)
    MdiV0(cmx_icntyp,i) = 1;              // default optimization goal: f(x) <= w

  //
  // fun: optimization goal
  //
  fun = mdr_Create(cmx_dimf, 1);
  mdr_Zero (fun);

  //
  // x:
  //
  nparm = cmx_dimx;
  iptb = cmx_dimg;
  indb = cmx_dimx;
  pttbl = mdr_Create(iptb,indb);

  liwrk = 7 * cmx_numgr + 7 * nparm + 3;
  iwork = mdi_Create(liwrk,1);

  lwrk = 2 * nparm * nparm + 4 * cmx_numgr * nparm + 11 * cmx_numgr +
      27 * nparm + 13;
  work  = mdr_Create(lwrk,1);

  // copy the initial point as column vector
  w = mdr_Float_BF(ystart);
  err = mdr_Create(cmx_numgr + 3, 1);

  //
  // the other arrays needed by the code
  //
  ioptn = 10000; // call  CMXFNSET  once per each 'x'
  // 10
  if (limsm > 0 && encsm >0)
  {
    ioptn = ioptn + 10;
    MdrV0 (work, 0) = encsm;
    MdiV0 (iwork,0) = limsm;
  }
  // 100
  if (nrkstep > 1)
  {
    MdiV0 (iwork,1) = nrkstep;
    ioptn = ioptn + 100;
  }
  if (tolcon > 0)
  {
    MdrV0 (work, 1) = tolcon;
    ioptn = ioptn + 200;
  }
  // 1000
  ioptn = ioptn + 1000 * istep;

  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr,
             "RLaB: using solver CONMAX for unconstrained"
                 " optimization in multidimensions.\n"
            );
    fprintf (fptr,
             "RLaB: Stream %s: Messages from the solver follow.\n", outs
            );
  }
  else
    louts = 0;

  //
  // call CONMAX now
  //
  CONMAX2(&ioptn, &nparm, &cmx_numgr, &maxiter, MDRPTR(fun), &cmx_dimf,
           MDRPTR(pttbl), &iptb, &indb, MDIPTR(iwork), &liwrk, MDRPTR(work), &lwrk,
                  &iter, MDRPTR(w), MDRPTR(err), outs, &louts, CMXFNSET);

  if (fptr)
  {
    fprintf (fptr, "RLaB: solver CONMAX finished the optimization search.\n");
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: optimization lasted %g sec.\n", timer);
    fclose (fptr);
  }

    // clean-up allocated memory arrays
  mdr_Destroy (err);
  mdr_Destroy (work);
  mdr_Destroy (iwork);
  mdr_Destroy (pttbl);
  mdr_Destroy (cmx_icntyp);
  mdr_Destroy (fun);

  return w;
}

//
// unconstrained optimization with Proximal Bundle: 811.f
//
MDR * mdr_findmins_f_df_pbun(int imethod, int cmx_dimf, MDR * ystart, char * outs,
                             int pb_met, int pb_mes, int pb_mtesx, int pb_mtesf, int maxiter,
                             double pb_tolx, double pb_tolf, double pb_tolb,
                             double pb_tolg, double pb_told, double pb_tols,
                             double pb_tolp, double pb_eta,  double pb_xmax
                             )
{
  MDR *w=0;
  FILE *fptr=0;
  int louts=isvalidstring(outs);
  time_t t1, t2;
  double timer;
  int na, nia, nra, nf, ipar[7], iterm , ihes=0;
  double rpar[9], fp = 0, gmax;

  MDR *ia=0, *ra=0;

  nf  = cmx_dimx;
  na  = nf + 3;

  nia = nf + na + 1;
  ia = mdi_Create(nia, 1);

  if (imethod == 5)
  {
    // pbunu
    nra = nf*(nf+1)/2+nf*(na+5)+5*na+4;
    ipar[1 -1] = pb_met;
  }
  else
  {
    // pbnew
    nra = nf*(nf+1)*(na+3)/2+nf*(na+6)+5*na+4;
    ipar[1 -1] = 1;   // mos - distance measure exponent, 1 or 2
  }
  ra  = mdr_Create(nra,1);

  // integer parameters
  ipar[2 -1] = pb_mes;
  ipar[3 -1] = pb_mtesx;
  ipar[4 -1] = pb_mtesf;
  ipar[5 -1] = maxiter;
  ipar[6 -1] = maxiter;
  if (louts > 1)
    ipar[7 -1] = -2;  // maximum verbosity

    // double parameter
  rpar[1 -1] = pb_tolx;
  rpar[2 -1] = pb_tolf;
  rpar[3 -1] = pb_tolb;
  rpar[4 -1] = pb_tolg;
  rpar[5 -1] = pb_told;
  rpar[6 -1] = pb_tols;
  rpar[7 -1] = pb_tolp;
  rpar[8 -1] = pb_eta;
  rpar[9 -1] = pb_xmax;

  // copy the initial point as column vector
  w = mdr_Float_BF(ystart);

  //
  // call PBUNU now
  //
  if (louts > 1)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr,
             "RLaB: using solver PBUN/PNEW for unconstrained"
                 " optimization in multidimensions.\n"
            );
    fprintf (fptr,
             "RLaB: Stream %s: Messages from the solver follow.\n", outs
            );
  }
  else
    louts = 0;


  if (imethod == 5)
    PBUNU (&nf, &na, MDRPTR(w), MDIPTR(ia), MDRPTR(ra),
            ipar, rpar, &fp, &gmax, &iterm,
            outs, &louts);
  else
    PNEWU (&nf, &na, MDRPTR(w), MDIPTR(ia), MDRPTR(ra),
            ipar, rpar, &fp, &gmax, &ihes, &iterm,
            outs, &louts);

  if (fptr)
  {
    fprintf (fptr, "RLaB: solver PBUN/PNEW finished the optimization search.\n");
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: optimization lasted %g sec.\n", timer);
    fclose (fptr);
  }

  // clean-up allocated memory arrays
  mdr_Destroy (ia);
  mdr_Destroy (ra);

  return w;
}









