// rlabplus (C) 2003-2008 Marijan Kostrun
//
// GSL Science Library - numerical integration using monte carlo
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

/////////////////////////////////////////////////////////
//                                                     //
//                                                     //
//    M O N T E   C A R L O   I N T E G R A T I O N    //
//                                                     //
//                                                     //
/////////////////////////////////////////////////////////

static double
mcnint_gslrlab_f (double *x, size_t dim, void *dummy);

static MDR *xmdr;
static Ent *pent, *xent;
static Ent *fname_ent;

// **************************************************************
// * Builtin interface to gsl monte carlo numerical integration *
// **************************************************************
//
// montecarlo: all
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mcnint"
Ent * ent_nintegrate_mc (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *eo=0;
  MDR *x3, *xlo, *xhi;

  //
  // vegas parameters
  //
  int vegas_mode = 2, vegas_verbose = -1;
  double eabs = 1e-4, chicomp = 0.5, alpha = 0;
  char *vegas_outstream=0;

  //
  // miser parameters
  //
  double efrac = 0.1, dither = 0.1;

  int xr, idummy, imethod = 0;
  double result, ddummy;
  size_t calls = 100000, therm_calls = 0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  ListNode *node;

  // MC uses its own IRNG:
  int irng = RLAB_IRNG_MC - 1;

  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);

  // parameters and initial message
  if (nargs < 3 || nargs > 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Monte Carlo integration of a real function\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   mcnint(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where f = function(x/,p/), 'p' is its parameter array,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": X=[xlo,xhi] are the finite bounds of integration in\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": multidimensions. The solver depends on options=<<imethod;..>>\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where, imethod=0 for plain, 1 for miser and 2 for vegas.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Check the manual for other options.\n");
    rerror (THIS_SOLVER "requires at least 4 arguments");
  }
  //
  // Get function ptr
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get function parameters
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy(e2);

  //
  // x(lo), lower end of integration interval for each dimension
  //
  e3  = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_REAL)
    rerror ("mcnint: A matrix of integration bounds X=[xlo,xhi] is not a real matrix!");
  x3 = class_matrix_real (e3);
  if (MNC(x3)!=2)
    rerror ("mcnint: X=[xlo,xhi] has to be a two-column matrix!");
  xlo = mdr_PartitionCol(x3,1);
  xhi = mdr_PartitionCol(x3,2);

  xr = MNR (xlo);

  gsl_monte_function F = { &mcnint_gslrlab_f, xr, 0 };

  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // imethod
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_IMETHOD);
      if (node != 0)
      {
        idummy = class_double ( var_ent (node) );
        if (idummy >= 0 && idummy <= 2)
        {
          imethod = idummy;
          if (imethod == 1) alpha = 2.0;
          if (imethod == 2) alpha = 1.5;
        }
      }
      // ncalls: plain, miser, vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_NCALL);
      if (node != 0)
      {
        idummy = class_double ( var_ent (node) );
        if (idummy > 0)
          calls = idummy;
      }
      // ntherm: vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_NTH);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy > 0)
          therm_calls = idummy;
      }
      // chicomp: vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_CHI);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 2)
          chicomp = ddummy;
      }
      // mode: vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_MODE);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        switch (idummy)
        {
          case 0:
            vegas_mode = GSL_VEGAS_MODE_IMPORTANCE;
            break;
          case 1:
            vegas_mode = GSL_VEGAS_MODE_STRATIFIED;
            break;
          default:
            vegas_mode = GSL_VEGAS_MODE_IMPORTANCE_ONLY;
            break;
        }
      }
      // alpha: miser, vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_ALPHA);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 10.0)
          alpha = ddummy;
      }
      // efrac: miser
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_EF);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 2.0)
          efrac = ddummy;
      }
      // dither: miser
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_DITH);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 10.0)
          dither = ddummy;
      }
      // eabs: plain, miser, vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MC_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 1.0)
          eabs = ddummy;
      }
      // stdout: vegas
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          vegas_outstream = class_char_pointer(var_ent (node));
          if (strlen(vegas_outstream) > 1)
            vegas_verbose = 2;
        }
      }
    }
    else if (ent_type(eo) == MATRIX_DENSE_REAL)
    {
      imethod = 0;
      // ncalls
      idummy = class_double ( eo );
      if (idummy > 0) calls = idummy;
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  //
  // Set up ENTITIES for rlab function being integrated
  // x
  xmdr = mdr_CreateEmpty (xr, 1);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  if (imethod == 0)
  {
    //
    // monte carlo plain
    //
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (xr);
    gsl_monte_plain_integrate (&F, MDRPTR(xlo), MDRPTR(xhi), xr, calls,
                                rlab_gsl_rng_r[irng], s, &result, &eabs);
    gsl_monte_plain_free (s);
  }
  else if (imethod == 1)
  {
    //
    // monte carlo miser
    //
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (xr);
    // user supplied parameters for MISER read the manual
    s->estimate_frac = efrac;     // default value = 0.1
    s->alpha = alpha;             // default value = 2
    s->dither = dither;           // default value = 0.1
    gsl_monte_miser_integrate (&F, MDRPTR(xlo), MDRPTR(xhi), xr, calls,
                                rlab_gsl_rng_r[irng], s, &result, &eabs);
    gsl_monte_miser_free (s);
  }
  else
  {
    //
    // monte carlo vegas
    //
    if (therm_calls == 0) therm_calls = calls / 10;
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (xr);
    // first check if the outstream works, if not reset vegas_verbose=-1
    FILE *fptr = NULL;
    if (vegas_outstream)
      fptr = fopen (vegas_outstream, "a");
    if (!fptr)
      vegas_verbose = -1;
    else
    {
      s->verbose = vegas_verbose;
      s->ostream = fptr;
    }
    // put the rest of the parameters to
    s->alpha = alpha;
    s->mode = vegas_mode;
    s->stage = 0; // use thermalization to tune the integrator
    // vegas: thermalization
    gsl_monte_vegas_integrate (&F, MDRPTR(xlo), MDRPTR(xhi), xr, therm_calls,
                                rlab_gsl_rng_r[irng], s, &result, &eabs);
    // vegas: integrate and check convergence using the original gsl criterion
    s->stage = 1; // use the tuned mesh created by the thermalization
    do
    {
      gsl_monte_vegas_integrate (&F, MDRPTR(xlo), MDRPTR(xhi), xr, calls,
                                  rlab_gsl_rng_r[irng], s, &result, &eabs);
    }
    while (fabs (s->chisq - 1.0) > chicomp);
    gsl_monte_vegas_free (s);

    if (fptr)
      fclose (fptr);
  }

  //
  // clean-up
  //
  ent_Clean (pent);
  ent_Clean (e2);

  mdr_Destroy (xlo);
  mdr_Destroy (xhi);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean(fname_ent);

  ent_Clean (e3);

  return ent_Create_Rlab_Double(result);
}



//
// The interface to the user-specified function.
//  MDR params is passed by default, no need for it to be the argument of the function
//
static double
mcnint_gslrlab_f (double *x, size_t dim, void *dummy)
{
  Ent *rent=0;
  double rval=0;

  MDPTR(xmdr) = (void *) x;

  if (pent)
    rent = ent_call_rlab_script_2args(fname_ent, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname_ent, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean(rent);
  return rval;
}
