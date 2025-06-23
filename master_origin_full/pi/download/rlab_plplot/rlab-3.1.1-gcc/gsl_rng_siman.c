// rlabplus (C) 2003-2005 Marijan Kostrun
//
// GSL Science Library - simulated annealing, basic algorithm
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING
// **********************************************************************

// other variables
static Ent *xe_ent, *xn_ent;
static MDR *xe, *xn;
static Ent *energy=0, *nextstep=0;

// functions in this library
static void mdr_copy_one_to_two(MDR *x1, MDR *x2);
static void gsl_siman_nextstep (MDR * x);
static double gsl_siman_energyfn (MDR * x);


MDR *
rlab_gsl_siman_solve (const gsl_rng * r, MDR * x0_p, gsl_siman_params_t params, FILE *fptr)
{
  double E, new_E, best_E;
  int i, j, done;
  double T, timer;
  int n_evals = 1, n_iter = 0, n_accepts, n_rejects, n_eless, xrow = 0, xcol =
    0;

  time_t t1;

  MDR *xp     = mdr_Float_BF (x0_p);
  MDR *new_x  = mdr_Float_BF (x0_p);
  MDR *best_x = mdr_Float_BF (x0_p);

  E = gsl_siman_energyfn (xp);

  best_E = E;
  T = params.t_initial;
  done = 0;

  if (fptr)
  {
    t1 = clock ();
    fprintf (fptr, "RLaB: using GSL simulated annealing solver 'siman'.\n");
    fprintf (fptr, "RLaB: Using the GSL random number generator %s.\n",
             gsl_rng_name (r));
    xrow = MNR (xp);
    xcol = MNC (xp);
  }

  while (!done)
  {
    n_accepts = 0;
    n_rejects = 0;
    n_eless = 0;
    // print configuration if requested
    if (fptr)
    {
      for (i = 0; i < xrow; i++)
      {
        fprintf (fptr, "siman:");
        for (j = 0; j < xcol; j++)
          fprintf (fptr, " x[%i,%i] = %g ,", i + 1, j + 1, Mdr0 (best_x, i, j));
        fprintf (fptr, " energy = %g", best_E);
        fprintf (fptr, "\n");
      }
    }
    for (i = 0; i < params.iters_fixed_T; ++i)
    {
      mdr_copy_one_to_two(xp, new_x);
      gsl_siman_nextstep(new_x);
      new_E = gsl_siman_energyfn (new_x);
      if (new_E <= best_E)
      {
        mdr_copy_one_to_two(new_x, best_x);
        best_E = new_E;
      }
      ++n_evals;
      if (new_E < E)
      {
        // yay! take a step
        mdr_copy_one_to_two(new_x, xp);
        E = new_E;
        ++n_eless;
      }
      else if (gsl_rng_uniform (r) < exp (-(new_E - E) / (params.k * T)))
      {
        // yay! take a step
        mdr_copy_one_to_two(new_x, xp);
        E = new_E;
        ++n_accepts;
      }
      else
      {
        ++n_rejects;
      }
    }
    // apply the cooling schedule to the temperature
    T /= params.mu_t;
    ++n_iter;
    if (T < params.t_min)
      done = 1;
  }
  // at the end, copy the result onto the initial point, so we pass it
  // back to the caller:
  if (fptr)
  {
    fprintf (fptr,
             "RLaB: the GSL simulated annealing solver reports 'annealing done' .\n");
    timer = (clock () - t1) / 1e6;
    fprintf (fptr, "RLaB: annealing lasted %g sec.\n", timer);
  }

  mdr_Destroy (xp);
  mdr_Destroy (new_x);

  return best_x;
}

#undef THIS_SOLVER
#define THIS_SOLVER "siman"
Ent *
ent_gsl_siman (int nargs, Datum args[])
{
  // call parameters:
  //  d1 - energy function  f(x)
  //      d2 - new step function s(x)
  //  e3 - output stream
  //  e4 - [n(tries), n(iterT)]
  //  e5 - [k, Tinit, mu, Tfinal]

  Ent *e3=0, *eo=0, *rent;
  MDR *x0=0;
  char *outs=0;
  int niter = 10, ntry = 200, idummy;
  double k = 1, tinit = 0.002, tfinal = 2e-6, mu = 1.005, ddummy;
  Btree *bw = btree_Create ();
  ent_type (bw) = BTREE;
  ListNode *node;

  FILE *fptr=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //

  // SIMAN uses its own IRNG:
  int irng = RLAB_IRNG_SIMAN - 1;

  if (rlab_gsl_irand[irng].index < 1)
  {
    rlab_setup_default_gsl_irng (irng);
  }

  if (nargs < 4)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Basic GSL simulated annealing solver.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   y=siman(energy,newstep,x0 /,options/)\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where, energy=function(x), given the configuration x,\n");
    fprintf (rlab_stderr, THIS_SOLVER ": newstep=function(x), determines a new configuration xnew,\n");
    fprintf (rlab_stderr, THIS_SOLVER ": 'x0' is the initial configuration.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": options=<<stdout;Ti;Tf;mu;nT>> controls the annealing process,\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where 'stdout' is the name of the output stream, 'Ti' is the\n");
    fprintf (rlab_stderr, THIS_SOLVER ": inital and 'Tf' the final temperature (Ti>Tf) and 'mu' is the\n");
    fprintf (rlab_stderr, THIS_SOLVER ": cooling rate (cr>1), while 'nT' is the number of attempts at a\n");
    fprintf (rlab_stderr, THIS_SOLVER ": particular temperature.\n");
    rerror (THIS_SOLVER ": requires at least 4 arguments !");
  }
  //
  // Get function ptr for energy
  //
  energy = bltin_get_ent(args[0]);
  if (!isfuncent(energy))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get function ptr for new step
  //
  nextstep = bltin_get_ent(args[1]);
  if (!isfuncent(nextstep))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // Get  initial configuration
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("siman: 'x0' has to be a real matrix!");
  x0 = class_matrix_real (e3);

  //
  // is there an output stream to which the run-time errors should be sent
  //
  //
  // options for the solver
  //
  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs  = class_char_pointer(var_ent (node));
        }
      }
      // initial temperature
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_SIMAN_TI);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0)
          tinit = ddummy;
      }
      // final temperature
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_SIMAN_TF);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0 && ddummy > tinit)
          tfinal = ddummy;
        else
          tfinal = tinit/10.0;
      }
      // cooling rate
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_SIMAN_MU);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 1.0)
          mu = ddummy;
      }
      // number of iterations at fixed temperature np
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_SIMAN_NT);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          ntry = idummy;
      }
    }
  }

  gsl_siman_params_t params = { ntry, niter, 1e100, k, tinit, mu, tfinal };

  //
  // Set up ENTITIES for user-function.
  //
  xe = mdr_CreateEmpty(MNR(x0),MNC(x0));
  xe_ent = ent_Assign_Rlab_MDR (xe);
  ent_IncRef (xe_ent);

  xn = mdr_CreateEmpty(MNR(x0),MNC(x0));
  xn_ent = ent_Assign_Rlab_MDR (xn);
  ent_IncRef (xn_ent);

  //
  // return the result as a list
  //
  if(outs)
    fptr = fopen(outs,"a");

  rent = ent_Assign_Rlab_MDR ( rlab_gsl_siman_solve (rlab_gsl_irand[irng].generator, x0, params, fptr) );

  if (fptr)
    fclose(fptr);

  // Clean Up
  ent_Clean (eo);
  ent_Clean (e3);
  ent_Clean (energy);
  ent_Clean (nextstep);

  MDPTR(xe)=0;
  ent_DecRef (xe_ent);
  ent_Destroy (xe_ent);

  MDPTR(xn) = 0;
  ent_DecRef (xn_ent);
  ent_Destroy (xn_ent);

  return rent;
}

static double
gsl_siman_energyfn (MDR * x)
{
  Ent *rent = 0;
  double rval;

  MDPTR(xe) = MDPTR(x);

  rent = ent_call_rlab_script_1arg (energy, xe_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}

static void
mdr_copy_one_to_two(MDR *x1, MDR *x2)
{
  int i, n = SIZE(x1);

  if (n != SIZE(x2))
    rerror(RLAB_ERROR_TERRIBLE_INTERNAL_ERROR ": " THIS_SOLVER ": x1 and x2 mismatch!");

  for (i=0; i<n; i++)
    MdrV0(x2,i) = MdrV0(x1,i);

  return;
}

static void
gsl_siman_nextstep (MDR * x)
{
  Ent *rent = 0;
  MDR *retm;
  int i, n;

  MDPTR(xn) = MDPTR(x);

  rent = ent_call_rlab_script_1arg (nextstep, xn_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);
  n = SIZE(retm);

  if (n != SIZE(x))
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i=0; i<n; i++)
    MdrV0(x,i) = MdrV0(retm,i);

  return;
}
