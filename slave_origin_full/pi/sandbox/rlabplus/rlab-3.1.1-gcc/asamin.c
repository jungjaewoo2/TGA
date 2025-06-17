// rlabplus (C) 2003-2008 Marijan Kostrun
//
// ASA - simulated annealing, advanced algorithm
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

// rlab headers, located in variable $RLAB_SDK
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mds.h"
#include "mdc.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

// gsl headers
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_siman.h>

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// asa library, use garbage collector for memory allocations
#include "gsl_rng.h"
#include "asa_usr_asa.h"
#define  calloc(a,b)  GC_malloc((a)*(b))
#define  free(a)      GC_free((a))
#include "asa.c"

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
//
// arguments for the user functions:
//
//

// cost function
// cost = function(x)
static Ent  *cost_name=0;
static Ent  *cost_ent_x=0;
static MDR  *cost_mdr_x=0;

static double
user_cost_function(
    double *parameter, double *parameter_minimum,
    double *parameter_maximum,
    double *tangents,
    double *curvature,
    int *number_parameters,
    int *parameter_type,
    int *valid_state_generated_flag,
    int *exit_status,
    USER_DEFINES *OPTIONS);

// new generated state
// newst = function(rnd,x[i],T[i],i)
static Ent *newst_name=0;
static Ent  *newst_ent_rnd=0, *newst_ent_x=0, *newst_ent_t=0, *newst_ent_i=0;
static MDR  *newst_mdr_rnd=0, *newst_mdr_x=0, *newst_mdr_t=0, *newst_mdr_i=0;

static double
user_generating_distrib (LONG_INT * seed,
                         ALLOC_INT * parameter_dimension,
                         ALLOC_INT index_v,
                         double temperature_v,
                         double init_param_temp_v,
                         double temp_scale_params_v,
                         double parameter_v,
                         double parameter_range_v,
                         double *last_saved_parameter,
                         USER_DEFINES * USER_OPTIONS);

// acceptance test
// acct = function(rnd, current_cost, previous_cost, cost_T)
static Ent *acct_name=0;
static Ent  *acct_ent_rnd=0, *acct_ent_ccost=0, *acct_ent_pcost=0, *acct_ent_costt=0 ;
static MDR  *acct_mdr_rnd=0, *acct_mdr_ccost=0, *acct_mdr_pcost=0, *acct_mdr_costt=0;
static void
user_acceptance_test (double current_cost,
                      double *parameter_lower_bound,
                      double *parameter_upper_bound,
                      ALLOC_INT * parameter_dimension,
                      USER_DEFINES * USER_OPTIONS);

// cost scheduling
// csched = function(cost_T);
static Ent  *csched_name=0;
static Ent  *csched_ent_costt=0;
static MDR  *csched_mdr_costt=0;
static double
user_cost_schedule (double test_temperature, USER_DEFINES * USER_OPTIONS);

// index of the RNG used by asa:
static int idrng=0;

//
// asa solver
//
#undef  THIS_SOLVER
#define THIS_SOLVER "asamin"
Ent *
ent_asa_min (int nargs, Datum args[])
{
  // call parameters:
  //  d1 - cost function  f(x)
  //  e2 - initial point, x0
  //  e3 - range of x  [xmin,xmax]

  Ent *e2=0, *e3=0, *eo=0;
  MDR *x0=0, *rx=0, *ptype=0, *w=0;

  double ddummy;

  char *out_name=0;
  char *out_devnull = "/dev/null";
  FILE *out_file=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  time_t timer;
  int i;

  //Btree *bw = btree_Create ();
  //ent_type (bw) = BTREE;
  ListNode *node;

  // asa variables
  double *parameter=0, *parameter_minimum=0, *parameter_maximum=0;
  ALLOC_INT number_parameters=0;
  int *parameter_type=0;
  int valid_state_generated_flag=0, exit_status=0;
  double *tangents, *curvature;
  double *current_user_parameter_temp=0, cost = 0;
  LONG_INT fake_seed=0;

  USER_DEFINES * ASA_USER_OPTIONS;
  ASA_USER_OPTIONS = (USER_DEFINES *) GC_malloc (sizeof (USER_DEFINES));

  ASA_USER_OPTIONS->Limit_Acceptances = 1000;
  ASA_USER_OPTIONS->Limit_Generated = 99999;
  ASA_USER_OPTIONS->Limit_Invalid_Generated_States = 1000;
  ASA_USER_OPTIONS->Accepted_To_Generated_Ratio = 1.0E-4;

  ASA_USER_OPTIONS->Cost_Precision = 1.0E-18;
  ASA_USER_OPTIONS->Maximum_Cost_Repeat = 5;
  ASA_USER_OPTIONS->Number_Cost_Samples = 5;
  ASA_USER_OPTIONS->Temperature_Ratio_Scale = 1.0E-5;
  ASA_USER_OPTIONS->Cost_Parameter_Scale_Ratio = 1.0;
  ASA_USER_OPTIONS->Temperature_Anneal_Scale = 100.0;

  ASA_USER_OPTIONS->Include_Integer_Parameters = FALSE;
  ASA_USER_OPTIONS->User_Initial_Parameters = FALSE;
  ASA_USER_OPTIONS->Sequential_Parameters = -1;
  ASA_USER_OPTIONS->Initial_Parameter_Temperature = 1.0;

  ASA_USER_OPTIONS->Acceptance_Frequency_Modulus = 100;
  ASA_USER_OPTIONS->Generated_Frequency_Modulus = 10000;
  ASA_USER_OPTIONS->Reanneal_Cost = 1;
  ASA_USER_OPTIONS->Reanneal_Parameters = TRUE;

  ASA_USER_OPTIONS->Delta_X = 0.0001;
  ASA_USER_OPTIONS->User_Tangents = FALSE;
  ASA_USER_OPTIONS->Curvature_0 = TRUE;

  // user function: if user does not specify them then ASA defaults
  // are used - yes I did read the original code:
  ASA_USER_OPTIONS->Cost_Schedule      = user_cost_schedule;
  ASA_USER_OPTIONS->Acceptance_Test    = user_acceptance_test;
  ASA_USER_OPTIONS->Generating_Distrib = user_generating_distrib;

  MDR *mdr_curvature=0, *mdr_tangents=0;

  // ASA uses rlabplus rng facility and a dedicated system irng.
  // user can choose its type, but cannot sample from it.
  int irng = RLAB_IRNG_SIMAN - 1;
  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);
  idrng = RLAB_RNG_ASA_UNIFORM - 1;
  if (rlab_bltin_rng[idrng].idist<1)
  {
    // this is for book-keeping purpose only
    int    idist  = RLAB_RNG_UNIFORM_IDX;
    int    nparam = RLAB_RNG_UNIFORM_NPARAM;
    double param[RLAB_RNG_UNIFORM_NPARAM]= RLAB_RNG_UNIFORM_PARAMS;
    rlab_setup_rng (&idrng, &idist, &nparam, param, &irng);
  }

  if (nargs < 4)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Adaptive simulated annealing solver.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   y=asa(cost,newstep,x0,b /,options/)\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where, cost=function(x), given the configuration x,\n");
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
  // cost function: function pointer
  //
  cost_name = bltin_get_ent(args[0]);
  if (!isfuncent(cost_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // newstate generating function: function pointer
  //
  newst_name = bltin_get_ent(args[1]);
  if (!isfuncent(newst_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // initial configuration: x0
  //
  e2 = bltin_get_ent (args[2]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
  x0 = class_matrix_real (e2);
  if (!EQVECT(x0))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");

  w  = mdr_Float_BF (x0);
  parameter = MDRPTR(w);
  number_parameters = SIZE(w);
  exit_status = 0;

  //
  // parameter range: rangeX
  //
  e3 = bltin_get_ent (args[3]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_2COLMATRIX "\n");
  rx = class_matrix_real (e3);
  if (SIZE(rx)<2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_2COLMATRIX "\n");
  if (MNC(rx)!=2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_2COLMATRIX "\n");

  parameter_minimum = GC_MALLOC(SIZE(x0) * sizeof(double));
  parameter_maximum = GC_MALLOC(SIZE(x0) * sizeof(double));

  // did user specify the valid state?
  valid_state_generated_flag = 1;
  for (i = 0; i < number_parameters; i++)
  {
    parameter_minimum[i] = mdr0(rx,MIN(i,MNR(rx)-1),0);
    parameter_maximum[i] = mdr0(rx,MIN(i,MNR(rx)-1),1);
    valid_state_generated_flag &= (parameter_minimum[i] <= mdrV0(x0,i)         );
    valid_state_generated_flag &= (         mdrV0(x0,i) <= parameter_maximum[i]);
  }

  //
  // acceptance test: function pointer
  //
  acct_name = 0;
  if (nargs > 4)
  {
    acct_name = bltin_get_ent (args[4]);
    if (!isfuncent(acct_name))
    {
      ent_Clean(acct_name);
      acct_name = 0;
    }
  }

  //
  // cost scheduling: function pointer
  //
  csched_name=0;
  if (nargs > 5)
  {
    csched_name = bltin_get_ent (args[5]);
    if (!isfuncent(csched_name))
    {
      ent_Clean(csched_name);
      csched_name = 0;
    }
  }

  //
  // options: last argument
  //
  eo = bltin_get_ent (args[nargs-1]);
  if (eo!=0)
    if (ent_type (eo) == BTREE)
    {
      // parameter_type
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_PT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ptype = ent_data (var_ent (node));
          if (ptype->nrow * ptype->ncol == number_parameters)
          {
            parameter_type = (int *) GC_malloc(number_parameters * sizeof(int));
            for (i=0; i<number_parameters; i++)
            {
              parameter_type[i] = mdrV0(ptype, i);
              if (parameter_type[i] != -2 &&
                  parameter_type[i] != -1 &&
                  parameter_type[i] !=  1 &&
                  parameter_type[i] !=  2)
                parameter_type[i] = -1; // default: real type
            }
          }
        }
      }
      // initial parameter temperature
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_INITPTEMP);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0)
            ASA_USER_OPTIONS->Initial_Parameter_Temperature = ddummy;
        }
      }
      // curvature
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_CURV);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy == 1)
            ASA_USER_OPTIONS->Curvature_0 = FALSE;
        }
      }
      // dx
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_DX);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0)
            ASA_USER_OPTIONS->Delta_X = ddummy;
        }
      }
      // sequential "cooling' of the parameters
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_SEQ);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy == 1)
            ASA_USER_OPTIONS->Sequential_Parameters = 1;
        }
      }
      // ASA_USER_OPTIONS->Acceptance_Frequency_Modulus
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_AFM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Acceptance_Frequency_Modulus =
                (int) ddummy;
        }
      }
      // ASA_USER_OPTIONS->Generated_Frequency_Modulus
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_GFM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Generated_Frequency_Modulus =
                (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Cost_Precision = 1.0E-18;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_COSTP);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0)
            ASA_USER_OPTIONS->Cost_Precision = ddummy;
        }
      }
      //ASA_USER_OPTIONS->Maximum_Cost_Repeat = 5;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_MAXCOSTR);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Maximum_Cost_Repeat = (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Number_Cost_Samples = 5;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_NOCOSTSM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Number_Cost_Samples = (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Temperature_Ratio_Scale = 1.0E-5;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_TEMPRATS);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0)
            ASA_USER_OPTIONS->Temperature_Ratio_Scale = ddummy;
        }
      }
      //ASA_USER_OPTIONS->Cost_Parameter_Scale_Ratio = 1.0;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_COSTPSRAT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Cost_Parameter_Scale_Ratio = ddummy;
        }
      }
      //ASA_USER_OPTIONS->Temperature_Anneal_Scale = 100.0;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_COSTPSRAT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Temperature_Anneal_Scale = ddummy;
        }
      }
      //ASA_USER_OPTIONS->Limit_Acceptances = 1000;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_LIMITACC);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Limit_Acceptances = (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Limit_Generated = 99999;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_LIMITGEN);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Limit_Generated = (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Limit_Invalid_Generated_States = 1000;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_LIMIGS);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy >= 1)
            ASA_USER_OPTIONS->Limit_Invalid_Generated_States = (int) ddummy;
        }
      }
      //ASA_USER_OPTIONS->Accepted_To_Generated_Ratio = 1.0E-4;
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ASA_AGR);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0)
            ASA_USER_OPTIONS->Accepted_To_Generated_Ratio = ddummy;
        }
      }
      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          out_name = class_char_pointer (var_ent (node));
          if (!isvalidstring(out_name))
            out_name=0;
        }
      }
    }

  // check that parameter_type is given
  if (!parameter_type)
  {
    parameter_type = (int *) GC_malloc(number_parameters * sizeof(int));
    for (i=0; i<number_parameters; i++)
      parameter_type[i] = -1; // default: real type
  }

  //
  // cost = function(x)
  //
  cost_mdr_x = mdr_CreateEmpty(w->nrow,w->ncol);
  cost_ent_x = ent_Assign_Rlab_MDR(cost_mdr_x);
  ent_IncRef (cost_ent_x);

  if (csched_name)
  {
    csched_mdr_costt = mdr_CreateEmpty(1,1);
    csched_ent_costt = ent_Assign_Rlab_MDR(csched_mdr_costt);
    ent_IncRef (csched_ent_costt);
  }

  //
  // newstate = function (rnd, x_i, T_i, i)
  //
  if (newst_name)
  {
    // rnd:
    newst_mdr_rnd = mdr_CreateEmpty(1,1);
    newst_ent_rnd = ent_Assign_Rlab_MDR (newst_mdr_rnd);
    ent_IncRef  (newst_ent_rnd);

    // x:
    newst_mdr_x = mdr_CreateEmpty(1,1);
    newst_ent_x = ent_Assign_Rlab_MDR (newst_mdr_x);
    ent_IncRef  (newst_ent_x);
    // T:
    newst_mdr_t = mdr_CreateEmpty(1,1);
    newst_ent_t = ent_Assign_Rlab_MDR (newst_mdr_t);
    ent_IncRef  (newst_ent_t);
    // i:
    newst_mdr_i = mdr_CreateEmpty(1,1);
    newst_ent_i = ent_Assign_Rlab_MDR (newst_mdr_i);
    ent_IncRef  (newst_ent_i);
  }

  //
  // acceptance_test = function (rnd, curr_cost, last_cost, cost_temp)
  //
  if (acct_name)
  {
    // rnd
    acct_mdr_rnd = mdr_CreateEmpty(1,1);
    acct_ent_rnd = ent_Assign_Rlab_MDR (acct_mdr_rnd);
    ent_IncRef  (acct_ent_rnd);
    // current cost
    acct_mdr_ccost = mdr_CreateEmpty(1,1);
    acct_ent_ccost = ent_Assign_Rlab_MDR (acct_mdr_ccost);
    ent_IncRef  (acct_ent_ccost);
    // previous cost
    acct_mdr_pcost = mdr_CreateEmpty(1,1);
    acct_ent_pcost = ent_Assign_Rlab_MDR (acct_mdr_pcost);
    ent_IncRef  (acct_ent_pcost);
    // cost_temp
    acct_mdr_costt = mdr_CreateEmpty(1,1);
    acct_ent_costt = ent_Assign_Rlab_MDR (acct_mdr_costt);
    ent_IncRef  (acct_ent_costt);
  }

  if (!out_name)
    out_name = out_devnull;
  else
    out_file = fopen(out_name,"w");

  ASA_USER_OPTIONS->Asa_Out_File = out_name;

  // write greetings
  if (out_file)
  {
    // start the timer
    timer = clock ();
    // write down boring info
    fprintf (out_file, "RLaB: Minimization using simulated annealing ASA code.\n");
  }

  //
  // calculation of tangents and curvatures
  //
  mdr_tangents = mdr_Create(1,number_parameters);
  tangents = MDRPTR(mdr_tangents);

  if (ASA_USER_OPTIONS->Curvature_0)
  {
    curvature = (double *) NULL;
  }
  else
  {
    mdr_curvature = mdr_Create(number_parameters,number_parameters);
    curvature = MDRPTR(mdr_curvature);
  }

  //
  // if user does not specify initial parameter temperature we have to:
  //
  if (!current_user_parameter_temp)
  {
    current_user_parameter_temp = (double *) GC_malloc(number_parameters * sizeof(double));
    for (i=0; i<number_parameters; i++)
      current_user_parameter_temp[i]= ASA_USER_OPTIONS->Initial_Parameter_Temperature;
  }
  ASA_USER_OPTIONS->User_Parameter_Temperature = current_user_parameter_temp;

  //
  // asa do my bidding!
  //
  cost = asa (user_cost_function,
              0, &fake_seed,
              parameter, parameter_minimum, parameter_maximum,
              tangents, curvature,
             &number_parameters, parameter_type,
             &valid_state_generated_flag, &exit_status,
              ASA_USER_OPTIONS);

  // tell me how you did
  if (out_file)
  {
    fprintf (out_file, "RLaB: built-in solver ASA reports status %i (", exit_status);

    // report what asa did
    switch (exit_status)
    {
      case 0: // normal exit
        fprintf (out_file, "success)\n");
        break;

      case 1:
        fprintf (out_file, "parameter temperature too small)\n");
        break;

      case 2:
        fprintf (out_file, "cost temperature too small)\n");
        break;

      case 3:
        fprintf (out_file, "cost started to repeat)\n");
        break;

      case 4:
        fprintf (out_file, "too many invalid states generated)\n");
        break;

      case 5:
        fprintf (out_file, "user function requested exit)\n");
        break;

      case 7:
        fprintf (out_file, "invalid input by the user)\n");
        break;

      case 8:
        fprintf (out_file, "invalid cost function)\n");
        break;

      case 9:
        fprintf (out_file, "invalid derivatives of the cost function)\n");
        break;

      case -1:
        fprintf (out_file, "memory allocation error)\n");
        break;
    }

    timer -= clock ();
    fprintf (out_file, "RLaB: ASA calculation lasted %g sec.\n", -timer / 1e6);

    fclose (out_file);
  }

  //
  // clean-up
  //

  // cost function
  MDPTR(cost_mdr_x) = 0;
  ent_DecRef (cost_ent_x);
  ent_Destroy (cost_ent_x);

  // new state
  if (newst_name)
  {
    // rnd:
    MDPTR(newst_mdr_rnd) = 0;
    ent_DecRef (newst_ent_rnd);
    ent_Destroy (newst_ent_rnd);
    // x[i]:
    MDPTR(newst_mdr_x) = 0;
    ent_DecRef (newst_ent_x);
    ent_Destroy (newst_ent_x);
    // T[i]:
    MDPTR(newst_mdr_t) = 0;
    ent_DecRef (newst_ent_t);
    ent_Destroy (newst_ent_t);
    // i:
    MDPTR(newst_mdr_i) = 0;
    ent_DecRef (newst_ent_i);
    ent_Destroy (newst_ent_i);
  }
  
  // acceptance test
  if (acct_name)
  {
    // rnd
    MDPTR(acct_mdr_rnd) = 0;
    ent_DecRef  (acct_ent_rnd);
    ent_Destroy (acct_ent_rnd);
    // current cost
    MDPTR(acct_mdr_ccost) = 0;
    ent_DecRef  (acct_ent_ccost);
    ent_Destroy (acct_ent_ccost);
    // previous cost
    MDPTR(acct_mdr_pcost) = 0;
    ent_DecRef  (acct_ent_pcost);
    ent_Destroy (acct_ent_pcost);
    // cost_temp
    MDPTR(acct_mdr_costt) = 0;
    ent_DecRef  (acct_ent_costt);
    ent_Destroy (acct_ent_costt);
  }

  // cost scheduling
  if (csched_name)
  {
    // cost_T:
    MDPTR(csched_mdr_costt) = 0;
    ent_DecRef  (csched_ent_costt);
    ent_Destroy (csched_ent_costt);
  }

  GC_free(parameter_type);
  GC_free(current_user_parameter_temp);
  GC_free(ASA_USER_OPTIONS);

  ent_Clean (cost_name);
  ent_Clean (newst_name);
  ent_Clean (acct_name);
  ent_Clean (csched_name);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (eo);
  GC_FREE (parameter_maximum);
  GC_FREE (parameter_minimum);

  //
  // return the result: cost_ent_x contains the last value of parameters
  //
  Btree * rtree = btree_Create ();
  install (rtree, "status", ent_Create_Rlab_Double (exit_status) );
  install (rtree, "cost", ent_Create_Rlab_Double(cost) );
  install (rtree, "coef", ent_Assign_Rlab_MDR(w) );
  install (rtree, "tangent", ent_Assign_Rlab_MDR(mdr_tangents) );
  if (mdr_curvature)
    install (rtree, "curvature", ent_Assign_Rlab_MDR (mdr_curvature) );
  return ent_Assign_Rlab_BTREE(rtree);
}


//
// cost = function(x)
//
static double user_cost_function(
    double *parameter, double *parameter_minimum,
    double *parameter_maximum,
    double *tangents,
    double *curvature,
    int *number_parameters,
    int *parameter_type,
    int *valid_state_generated_flag,
    int *exit_status,
    USER_DEFINES *OPTIONS)
{
  Ent *rent = 0;
  double ret;

  //
  // pass parameter as an entity for user-scripted cost function
  // and execute its script
  //
  MDPTR(cost_mdr_x) = (void *) parameter;
  rent = ent_call_rlab_script_1arg (cost_name, cost_ent_x);
  ret  = class_double(rent);
  ent_Clean (rent);

  // if cost function returned 'nan' declare
  // invalid state
  if (isnand(ret))
    *valid_state_generated_flag = 0;
  else
    *valid_state_generated_flag = 1;

  return ret;
}

//
// newstate = function(rnd, x[i], T[i], i)
//
static double
user_generating_distrib (LONG_INT * seed,
                         ALLOC_INT * parameter_dimension,
                         ALLOC_INT index_v,
                         double temperature_v,
                         double init_param_temp_v,
                         double temp_scale_params_v,
                         double parameter_v,
                         double parameter_range_v,
                         double *last_saved_parameter,
                         USER_DEFINES * USER_OPTIONS)
{
  double ret;
  int ione=1;

  double x;
  gsl_random_number_generator (&ione, &ione, &x, &idrng);

  if (newst_name)
  {
    // user specified function for generating new set of parameters
    Ent *rent = 0;
    double idx= (double)index_v + 1.0; // index_v starts from 0

    MDPTR(newst_mdr_rnd) = (void *) &x;
    MDPTR(newst_mdr_x)   = (void *) &parameter_v;
    MDPTR(newst_mdr_t)   = (void *) &temperature_v;
    MDPTR(newst_mdr_i)   = (void *) &idx;

    rent = ent_call_rlab_script_4args (newst_name, newst_ent_rnd, newst_ent_x, newst_ent_t, newst_ent_i);
    ret  = class_double(rent);
    ent_Clean (rent);
  }
  else
  {
    // user did not specify function to generate a new state of parameters
    // we merely copy generate_asa_state
    double y, z;
    y = x < 0.5 ? -1 : 1;
    z = y * temperature_v * (pow((1.0 + 1.0/(temperature_v)), ABS (2.0*x - 1.0)) - 1.0);
    ret = parameter_v + z * parameter_range_v;
  }

  return ret;
}

//
// user_acceptance_test = fn(rnd, old cost, new cost, cost_temp)
//
static void
user_acceptance_test (double current_cost,
                      double *parameter_lower_bound,
                      double *parameter_upper_bound,
                      ALLOC_INT * parameter_dimension,
                      USER_DEFINES * USER_OPTIONS)
{
  double uniform_test;
  double last_cost=*(USER_OPTIONS->Last_Cost);
  double curr_cost_temp;

  int ione=1;

  //
  // did user specify her own cost_schedule?
  //
  if (csched_name)
  {
    curr_cost_temp =
        user_cost_schedule (USER_OPTIONS->Cost_Temp_Curr,USER_OPTIONS)
        + (double) EPS_DOUBLE;
  }
  else
  {
    curr_cost_temp = USER_OPTIONS->Cost_Temp_Curr;
  }

  // generate a random number U[0,1]: use dedicated RNG for asa/siman
  gsl_random_number_generator (&ione, &ione, &uniform_test, &idrng);

  if (acct_name)
  {
    //
    // call user function, to determine if the new cost passes the test
    //
    Ent *rent = 0;
    double ret;

    // rnd
    MDPTR(acct_mdr_rnd)   = (void *) &uniform_test;
    // current cost
    MDPTR(acct_mdr_ccost) = (void *) &current_cost;
    // previous (accepted) cost
    MDPTR(acct_mdr_pcost) = (void *) &last_cost;
    // cost temperature
    MDPTR(acct_mdr_costt) = (void *) &curr_cost_temp;

    rent = ent_call_rlab_script_4args (acct_name, acct_ent_rnd, acct_ent_ccost, acct_ent_pcost, acct_ent_costt);
    ret  = class_double(rent);
    ent_Clean (rent);

    //user reports 1 (accept) or 0 (reject)
    if (ret)
      USER_OPTIONS->User_Acceptance_Flag = TRUE;
    else
      USER_OPTIONS->User_Acceptance_Flag = FALSE;
  }
  else
  {
    // use Boltzmann's test
    double delta_cost = (current_cost - last_cost) / (curr_cost_temp + (double) EPS_DOUBLE);
    double x = MIN ( 1.0, exp(-delta_cost) );       /* Boltzmann test */

    USER_OPTIONS->Prob_Bias = x;
    if (x >= uniform_test)
      USER_OPTIONS->User_Acceptance_Flag = TRUE;
    else
      USER_OPTIONS->User_Acceptance_Flag = FALSE;
  }
  return;
}



//
// user_cost_schedule = function (cost_T)
//
static double
user_cost_schedule (double test_temperature, USER_DEFINES * USER_OPTIONS)
{
  double ret = test_temperature;
  Ent *rent=0;

  if (csched_name)
  {
    MDPTR(csched_mdr_costt) = (void *) &test_temperature;
    rent = ent_call_rlab_script_1arg (csched_name, csched_ent_costt);
    ret = class_double(rent);
    ent_Clean (rent);
  }

  return ret;
}
