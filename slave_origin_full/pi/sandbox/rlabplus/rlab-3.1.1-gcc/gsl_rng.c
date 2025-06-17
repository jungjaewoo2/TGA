// rlabplus (C) 2003-2005 Marijan Kostrun
//
// GSL Science Library - random number generators and dependent functions
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
#include "mdr_mdc.h"


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
#include <gsl/gsl_histogram.h>
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

// ***********************************************************************
//
// Integer random number generator facility in rlab+rlabplus
//
// ***********************************************************************
#define RLAB_IRNG_MAXNO     64
#define RLAB_IRNG_MAXUSERNO 32

#define RLAB_IRNG_STAT    (RLAB_IRNG_MAXNO-8)
#define RLAB_IRNG_SIMAN   (RLAB_IRNG_MAXNO-7)
#define RLAB_IRNG_MC      (RLAB_IRNG_MAXNO-6)
#define RLAB_IRNG_DRNG    (RLAB_IRNG_MAXNO-5)
#define RLAB_IRNG_HRNG    (RLAB_IRNG_MAXNO-4)
#define RLAB_IRNG_SAMPLE  (RLAB_IRNG_MAXNO-3)
#define RLAB_IRNG_SHUFFLE (RLAB_IRNG_MAXNO-2)
#define RLAB_IRNG_NORMAL  (RLAB_IRNG_MAXNO-1)
#define RLAB_IRNG_UNIFORM (RLAB_IRNG_MAXNO)

#define RLAB_GSL_IRNG_DEFAULT gsl_rng_gfsr4

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

int RLAB_IRNG_DEFAULT = 0;
gsl_rng *rlab_gsl_rng_r[RLAB_IRNG_MAXNO] = {0};

struct _int_ran_gen
{
  int  index;
  char name[32];
};

static struct _int_ran_gen
bltin_iran[] =
{
  // index, name, would have put the pointer to the generator but
  // cannot do it
  {-1, "none" },
  { 1, "borosh13" },
  { 2, "cmrg" },
  { 3, "coveyou" },
  { 4, "fishman18" },
  { 5, "fishman20" },
  { 6, "fishman2x" },
  { 7, "gfsr4" },
  { 8, "knuthran" },
  { 9, "knuthran2" },
  {10, "lecuyer21" },
  {11, "minstd" },
  {12, "mrg" },
  {13, "mt19937" },
  {14, "mt19937_1999" },
  {15, "mt19937_1998" },
  {16, "r250" },
  {17, "ran0" },
  {18, "ran1" },
  {19, "ran2" },
  {20, "ran3" },
  {21, "rand" },
  {22, "rand48" },
  {23, "random128_bsd" },
  {24, "random128_glibc2" },
  {25, "random128_libc5" },
  {26, "random256_bsd" },
  {27, "random256_glibc2" },
  {28, "random256_libc5" },
  {29, "random32_bsd" },
  {30, "random32_glibc2" },
  {31, "random32_libc5" },
  {32, "random64_bsd" },
  {33, "random64_glibc2" },
  {34, "random64_libc5" },
  {35, "random8_bsd" },
  {36, "random8_glibc2" },
  {37, "random8_libc5" },
  {38, "random_bsd" },
  {39, "random_glibc2" },
  {40, "random_libc5" },
  {41, "randu" },
  {42, "ranf" },
  {43, "ranlux" },
  {44, "ranlux389" },
  {45, "ranlxd1" },
  {46, "ranlxd2" },
  {47, "ranlxs0" },
  {48, "ranlxs1" },
  {49, "ranlxs2" },
  {50, "ranmar" },
  {51, "slatec" },
  {52, "taus" },
  {53, "taus2" },
  {54, "taus113" },
  {55, "transputer" },
  {56, "tt800" },
  {57, "uni" },
  {58, "uni32" },
  {59, "vax" },
  {60, "waterman14" },
  {61, "zuf" },
  {62, "default"},
  { 0, ""}
};

static int rlab_iran = 7; // gfsr4

//
// get an integer from local source of entropy
//
unsigned long int
iseed_uran (void)
{
  unsigned long int seed = 0;
  FILE *fh = NULL;
  // getting the random number seed from /dev/urandom
  fh = fopen ("/dev/urandom", "r");
  fread (&seed, sizeof (seed), 1, fh);
  fclose (fh);
  return seed;
}

// int
// rlab_setup_default_gsl_irng (void)
// {
//   gsl_rng_default = gsl_rng_gfsr4;
//   gsl_rng_default_seed = iseed_uran ();
//   gsl_rng_r = gsl_rng_alloc (gsl_rng_default);
//   return 0;
// }

int
rlab_setup_default_gsl_irng (int idx)
{
  if (!rlab_gsl_rng_r[idx])
  {
    //
    // setting up new generator
    //
    rlab_gsl_rng_r[idx]
        = gsl_rng_alloc (RLAB_GSL_IRNG_DEFAULT);
    gsl_rng_set(rlab_gsl_rng_r[idx], iseed_uran());

    return 0;
  }

  // user messed up - she tried to setup a generator that
  // has already been setup
  return 1;
}

//
// set the integer random number generator
//
Ent *
ent_gsl_irng (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int i, n, ix;

  unsigned long int nseed=0L;

  //
  // no arguments: write the information about the avaialable generators
  //
  if (nargs==0)
  {
    fprintf(stdout, "Currently defined integer random number generators :\n");
    for (i=0; i<RLAB_IRNG_MAXUSERNO; i++)
    {
      if (rlab_gsl_rng_r[i])
      {
        if (RLAB_IRNG_DEFAULT == i)
          fprintf(stdout, " user * ");
        else
          fprintf(stdout, " user   ");

        fprintf(stdout, "%3i: %s, size = %6i\n", i+1,
                gsl_rng_name (rlab_gsl_rng_r[i]),
                (int) gsl_rng_size (rlab_gsl_rng_r[i])
               );
      }
    }
    for (i=RLAB_IRNG_MAXUSERNO; i<RLAB_IRNG_MAXNO; i++)
    {
      if (rlab_gsl_rng_r[i])
      {
        n = gsl_rng_size (rlab_gsl_rng_r[i]);
        fprintf(stdout, " system %3i: %s, size = %6i\n",
                i+1, gsl_rng_name(rlab_gsl_rng_r[i]),
                n);
      }
    }
    fprintf(stdout, "(*) is a user-default generator.\n");

    return ent_Create_Rlab_Success();
  }


  //
  // irng(ix): set ix as a default IRNG
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("irng: missing index of the generator");

  ix = (int) class_double( e1 );
  if (ix<1 || ix>RLAB_IRNG_MAXNO)
    rerror("irng: index out of range");

  ix--;
  if (ix < RLAB_IRNG_MAXUSERNO)
  {
    // default generator can be only in range
    // 0 .. RLAB_IRNG_MAXUSERNO-1
      RLAB_IRNG_DEFAULT = ix;
  }

  if (nargs==1)
  {
    //
    // irng(ix): set ix as a default IRNG
    //
    if (!rlab_gsl_rng_r[ix])
    {
      // initialize the integer random number generator
      rlab_setup_default_gsl_irng(ix);
    }

    // cleanup
    ent_Clean (e1);

    // go back
    return ent_Create_Rlab_Double(RLAB_IRNG_DEFAULT+1.0);
  }

  //
  // two arguments  : irng(I,genname), or
  // three arguments: irng(I,genname,seed)
  //

  // get the name of the generator
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
    rerror("irng: missing name of the generator");
  char *s1 = class_char_pointer(e2);
  if (!s1)
    rerror ("irng: missing name of the generator");
  if (strlen(s1) == 0)
    rerror ("irng: missing name of the generator");

  if (!strcmp (s1, "help") || !strcmp (s1, "?"))
  {
    fprintf (stdout, "irng: available integer random generators"
        " (second argument):\nirng: ");
    for (i = 1; bltin_iran[i].index; i++)
    {
      if (!(i%5))
        fprintf (stdout, "\nirng: ");
      else
        fprintf (stdout, "'%s',", bltin_iran[i].name);
    }
    fprintf (stdout, "\nirng: please check the manual.\n");
    fflush (stdout);
    rerror ("improper first argument");
  }

  //
  // check the name and if is OK set it as pdfname
  //
  for (i = 0; bltin_iran[i].index; i++)
  {
    if (!strcmp (s1, bltin_iran[i].name))
    {
      rlab_iran = bltin_iran[i].index;
      break;
    }
  }

  // did the user provide the seed
  if (nargs==3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      nseed = (unsigned long int) class_double (e3);
    else
      nseed = iseed_uran ();
  }
  else
    nseed = iseed_uran ();

  //
  // decide on random number generator
  //
  // clean previous instantiation
  if (rlab_gsl_rng_r[ix])
    gsl_rng_free(rlab_gsl_rng_r[ix]);
  rlab_gsl_rng_r[ix] = 0;

  // new instantiation
  switch (rlab_iran)
  {
    case  0:
      // not really a generator: use it to clear all memory
      // associated with the irng
      break;

    case  1:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_borosh13);
      break;

    case  2:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_cmrg);
      break;

    case  3:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_coveyou);
      break;

    case  4:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman18);
      break;

    case  5:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman20);
      break;

    case  6:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman2x);
      break;

    case 62: case  7:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_gfsr4);
      break;

    case  8:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_knuthran);
      break;

    case  9:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_knuthran2);
      break;

    case 10:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_lecuyer21);
      break;

    case 11:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_minstd);
      break;

    case 12:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mrg);
      break;

    case 13:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937);
      break;

    case 14:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937_1999);
      break;

    case 15:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937_1998);
      break;

    case 16:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_r250);
      break;

    case 17:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran0);
      break;

    case 18:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran1);
      break;

    case 19:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran2);
      break;

    case 20:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran3);
      break;

    case 21:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_rand);
      break;

    case 22:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_rand48);
      break;

    case 23:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_bsd);
      break;

    case 24:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_glibc2);
      break;

    case 25:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_libc5);
      break;

    case 26:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_bsd);
      break;

    case 27:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_glibc2);
      break;

    case 28:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_libc5);
      break;

    case 29:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_bsd);
      break;

    case 30:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_glibc2);
      break;

    case 31:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_libc5);
      break;

    case 32:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_bsd);
      break;

    case 33:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_glibc2);
      break;

    case 34:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_libc5);
      break;

    case 35:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_bsd);
      break;

    case 36:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_glibc2);
      break;

    case 37:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_libc5);
      break;

    case 38:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_bsd);
      break;

    case 39:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_glibc2);
      break;

    case 40:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_libc5);
      break;

    case 41:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_randu);
      break;

    case 42:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranf);
      break;

    case 43:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlux);
      break;

    case 44:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlux389);
      break;

    case 45:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxd1);
      break;

    case 46:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxd2);
      break;

    case 47:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs0);
      break;

    case 48:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs1);
      break;

    case 49:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs2);
      break;

    case 50:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranmar);
      break;

    case 51:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_slatec);
      break;

    case 52:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus);
      break;

    case 53:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus2);
      break;

    case 54:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus113);
      break;

    case 55:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_transputer);
      break;

    case 56:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_tt800);
      break;

    case 57:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_uni);
      break;

    case 58:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_uni32);
      break;

    case 59:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_vax);
      break;

    case 60:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_waterman14);
      break;

    case 61:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_zuf);
      break;
  }

  //
  // adjust the default generator if it got erased with "none"
  // start from the top; assign 0 if nothing found
  //
  if (!rlab_gsl_rng_r[RLAB_IRNG_DEFAULT])
  {
    for (i = RLAB_IRNG_MAXUSERNO; i>0; i--)
    {
      RLAB_IRNG_DEFAULT = i-1;
      if (rlab_gsl_rng_r[RLAB_IRNG_DEFAULT])
        break;
    }
  }

  // set the seed
  if (rlab_gsl_rng_r[ix])
    gsl_rng_set(rlab_gsl_rng_r[ix], nseed);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

//
// get/set the state of an integer random number generator
//
Ent *
ent_gsl_irng_state (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i, n, ix, irng=0;

  char *s1=0;

  MDR *istate=0;
  unsigned char *state=0;

  Btree *bw=0;

  if (nargs!=1 && nargs !=2)
    rerror ("irng_state: one or two arguments required");

  if (nargs == 1)
  {
    // obtain the state of the generator
    //

    // irng(ix): set ix as a default IRNG
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror("irng_state: missing index of the generator");

    ix = (int) class_double( e1 );
    if (ix<1 || ix>RLAB_IRNG_MAXNO)
      rerror("irng_state: index out of range");

    ix--;
    if (!rlab_gsl_rng_r[ix])
      rerror("irng_state: index out of range");

    s1     = (char*) gsl_rng_name (rlab_gsl_rng_r[ix]);
    n      = gsl_rng_size (rlab_gsl_rng_r[ix]);
    state  = gsl_rng_state(rlab_gsl_rng_r[ix]);

    istate = mdi_Create(1,n);

    for (i=0;i<n;i++)
      MdiV0(istate, i) = state[i];

    ent_Clean (e1);

    bw = btree_Create ();
    install (bw, "name", ent_Create_Rlab_String (s1));
    install (bw, "state", ent_Assign_Rlab_MDR(istate));

    return ent_Assign_Rlab_BTREE (bw);
  }

  ListNode *node=0;

  //
  // irng_state(i,x):
  //   set IRNG 'i' to data from 'x'
  //

  // first argument: i, index of the IRNG
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("irng_state: missing index of the generator");

  ix = (int) class_double( e1 );
  if (ix<1 || ix>RLAB_IRNG_MAXNO)
    rerror("irng_state: index out of range");

  ix--;
  if (rlab_gsl_rng_r[ix])
  {
    // free previous instantiation of the generator
    gsl_rng_free(rlab_gsl_rng_r[ix]);
    rlab_gsl_rng_r[ix] = 0;
  }

  // second argument: x, state of the IRNG
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror("irng_state: missing state of the generator");

  node = btree_FindNode (ent_data (e2), RLAB_NAME_RNG_NAME);
  if (!node)
    rerror("irng_state: missing state of the generator");
  s1 = class_char_pointer (var_ent (node));

  node = btree_FindNode (ent_data (e2), RLAB_NAME_RNG_STATE);
  if (!node)
    rerror("irng_state: missing state of the generator");
  istate = ent_data (var_ent (node));

  //
  // from the name find the generator
  //
  for (i = 0; bltin_iran[i].index; i++)
  {
    if (!strcmp (s1, bltin_iran[i].name))
    {
      irng = bltin_iran[i].index;
      break;
    }
  }
  if (!bltin_iran[i].index)
    rerror("irng_state: no such generator");


  // new instantiation
  switch (irng)
  {
    case -1:
      // not really a generator: use it to clear all memory
      // associated with the irng
      break;

    case  1:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_borosh13);
      break;

    case  2:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_cmrg);
      break;

    case  3:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_coveyou);
      break;

    case  4:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman18);
      break;

    case  5:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman20);
      break;

    case  6:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_fishman2x);
      break;

    case 62: case  7:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_gfsr4);
      break;

    case  8:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_knuthran);
      break;

    case  9:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_knuthran2);
      break;

    case 10:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_lecuyer21);
      break;

    case 11:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_minstd);
      break;

    case 12:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mrg);
      break;

    case 13:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937);
      break;

    case 14:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937_1999);
      break;

    case 15:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_mt19937_1998);
      break;

    case 16:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_r250);
      break;

    case 17:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran0);
      break;

    case 18:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran1);
      break;

    case 19:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran2);
      break;

    case 20:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ran3);
      break;

    case 21:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_rand);
      break;

    case 22:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_rand48);
      break;

    case 23:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_bsd);
      break;

    case 24:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_glibc2);
      break;

    case 25:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random128_libc5);
      break;

    case 26:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_bsd);
      break;

    case 27:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_glibc2);
      break;

    case 28:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random256_libc5);
      break;

    case 29:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_bsd);
      break;

    case 30:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_glibc2);
      break;

    case 31:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random32_libc5);
      break;

    case 32:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_bsd);
      break;

    case 33:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_glibc2);
      break;

    case 34:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random64_libc5);
      break;

    case 35:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_bsd);
      break;

    case 36:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_glibc2);
      break;

    case 37:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random8_libc5);
      break;

    case 38:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_bsd);
      break;

    case 39:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_glibc2);
      break;

    case 40:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_random_libc5);
      break;

    case 41:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_randu);
      break;

    case 42:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranf);
      break;

    case 43:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlux);
      break;

    case 44:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlux389);
      break;

    case 45:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxd1);
      break;

    case 46:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxd2);
      break;

    case 47:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs0);
      break;

    case 48:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs1);
      break;

    case 49:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranlxs2);
      break;

    case 50:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_ranmar);
      break;

    case 51:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_slatec);
      break;

    case 52:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus);
      break;

    case 53:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus2);
      break;

    case 54:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_taus113);
      break;

    case 55:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_transputer);
      break;

    case 56:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_tt800);
      break;

    case 57:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_uni);
      break;

    case 58:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_uni32);
      break;

    case 59:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_vax);
      break;

    case 60:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_waterman14);
      break;

    case 61:
      rlab_gsl_rng_r[ix]
          = gsl_rng_alloc(gsl_rng_zuf);
      break;
  }

  n      = gsl_rng_size (rlab_gsl_rng_r[ix]);
  state  = gsl_rng_state(rlab_gsl_rng_r[ix]);

  if (n != MNR(istate)*MNC(istate) )
    rerror ("irng_state: mismatch in state size of the generator");

  if (istate->type == RLAB_TYPE_INT32)
  {
    for (i=0;i<n;i++)
      state[i] = (unsigned char) MdiV0(istate, i);
  }
  else
  {
    for (i=0;i<n;i++)
      state[i] = (unsigned char) MdrV0(istate, i);
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

//
// sample the integer random number generator
//
Ent *
ent_gsl_irand (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr = 1, nc = 1, i, n;
  MDR *w;

  //
  // one argument: it is a matrix, read its dimensions and sample rng
  //
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else
      rerror ("irand: improper first argument");
  }

  //
  // two arguments: both reals - nrow and ncol of sample matrix
  //
  if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) == MATRIX_DENSE_REAL
        && ent_type (e2) == MATRIX_DENSE_REAL)
    {
      nr = class_double (e1);
      nc = class_double (e2);
    }
  }

  //
  // sample the generator now: only default generator can be sampled
  //
  if (!rlab_gsl_rng_r[RLAB_IRNG_DEFAULT])
    rlab_setup_default_gsl_irng (RLAB_IRNG_DEFAULT);

  w = mdi_Create (nr, nc);
  n = nr * nc;
  for (i = 0; i < n; i++)
    MdiV0 (w, i) = (int) gsl_rng_get (rlab_gsl_rng_r[RLAB_IRNG_DEFAULT]);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


// *******************************************************************
//
// rand()
//
// *******************************************************************

//
// list of distributions: NOT THE GENERATORS
//
struct _ran_pdf
{
  int  index;       // index (ordinal number)
  char name[32];    // name of the distribution
  int  nparam;                            // no. of parameters
};

static struct _ran_pdf
bltin_rng[] =
{
  // index, name of pdf, number of parameters supplied by user
  {-1, "none", 0},
  { 1, "normal", 2},
  { 2, "normaltail", 3},
  { 3, "exp", 1},
  { 4, "laplace", 2},
  { 5, "exppow", 3},
  { 6, "cauchy", 2},
  { 7, "rayleigh", 1},
  { 8, "rayleightail", 2},
  { 9, "levy", 3},
  {10, "levyskew", 4},
  {11, "gamma", 2},
  {12, "uniform", 2},
  {13, "default", 2},
  {14, "lognormal", 2},
  {15, "chisq", 1},
  {16, "F", 2},
  {17, "t", 1},
  {18, "beta", 2},
  {19, "logistic", 2},
  {20, "pareto", 2},
  {21, "weibull", 2},
  {22, "gumbel1", 2},
  {23, "gumbel2", 2},
  {24, "poisson", 1},
  {25, "binomial", 2},
  {26, "negbinomial", 2},
  {27, "geometric", 1},
  {28, "hypergeom", 3},
  {29, "log", 1},
  { 0, {0}, 0}
};

//
// list of RNG
//
// size of list
#define RLAB_RNG_MAXNO     64
#define RLAB_RNG_MAXUSERNO 32
#define RLAB_RNG_NORMAL      (RLAB_RNG_MAXNO-1)
#define RLAB_RNG_UNIFORM     (RLAB_RNG_MAXNO)

// default generator: uniform in [0,1]
#define RLAB_RNG_UNIFORM_IDX    12
#define RLAB_RNG_UNIFORM_NPARAM 2
#define RLAB_RNG_UNIFORM_PARAMS {0,1}

// additional generator: normal(0,1)
#define RLAB_RNG_NORMAL_IDX    1
#define RLAB_RNG_NORMAL_NPARAM 2
#define RLAB_RNG_NORMAL_PARAMS {0,1}

int RLAB_RNG_DEFAULT = 0;

struct _rlab_ran_pdf
{
  int     idist;    // index of the distribution
  int     nparam;   // no. of parameters
  double *param;    // the parameter array
  int     irng;     // index of the underlying integer rng
};

struct _rlab_ran_pdf
    rlab_bltin_rng[RLAB_RNG_MAXNO] = {{0, 0, 0, 0}};

//
// setup a rng using all necessary information
//
int
rlab_setup_rng (int *idrng, int *idist, int *nparam, double *param,
                int *irng)
{
  int i;

  // pass the fixed parameters
  rlab_bltin_rng[*idrng].idist  = *idist;
  rlab_bltin_rng[*idrng].nparam = *nparam;
  rlab_bltin_rng[*idrng].irng   = *irng;
  if(!rlab_gsl_rng_r[*irng])
    rlab_setup_default_gsl_irng (*irng);

  // reallocate memory for param arrays
  if (rlab_bltin_rng[*idrng].param)
    GC_free(rlab_bltin_rng[*idrng].param);
  rlab_bltin_rng[*idrng].param
      = (double *) GC_malloc(*nparam * sizeof(double));

  // copy param array
  for (i=0; i<*nparam;i++)
    rlab_bltin_rng[*idrng].param[i] = param[i];

  return 1;
}


void
gsl_random_number_generator (int *nr, int *nc, double *w, int *idx)
{
  int i, n = (*nr) * (*nc);
  int idrng, irng;

  if(idx)
    idrng = *idx;
  else
    idrng = RLAB_RNG_DEFAULT;

  //
  // setup the default rng: uniform on [0,1)
  //
  int     idist;
  int     nparam;
  double *param=0;

  // .idist==-1 by default, meaning 'none'
  if (rlab_bltin_rng[idrng].idist<1)
  {
    // setup default rng
    idist   = RLAB_RNG_UNIFORM_IDX;
    nparam  = RLAB_RNG_UNIFORM_NPARAM;
    double dparam[RLAB_RNG_UNIFORM_NPARAM] = RLAB_RNG_UNIFORM_PARAMS;
    rlab_setup_rng (&idrng, &idist,
                    &nparam, dparam, &RLAB_IRNG_DEFAULT);
  }

  // get the data for the generator
  idist = rlab_bltin_rng[idrng].idist;
  param = rlab_bltin_rng[idrng].param;
  irng  = rlab_bltin_rng[idrng].irng;

  switch (idist)
  {
    case -1:
      for (i = 0; i < n; i++)
        w[i] = 0;
      break;

    case  1:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_gaussian (rlab_gsl_rng_r[irng], param[1]);
      break;

    case  2:
      for (i = 0; i < n; i++)
        w[i] =
            param[0]
            + gsl_ran_gaussian_tail (rlab_gsl_rng_r[irng],
                                     param[2] - param[0],
                                     param[1]);
      break;

    case  3:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_exponential (rlab_gsl_rng_r[irng], param[0]);
      break;

    case  4:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_laplace (rlab_gsl_rng_r[irng], param[1]);
      break;

    case  5:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_exppow (rlab_gsl_rng_r[irng], param[1], param[2]);
      break;

    case  6:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_cauchy (rlab_gsl_rng_r[irng], param[1]);
      break;

    case  7:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_rayleigh (rlab_gsl_rng_r[irng], param[0]);
      break;

    case  8:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_rayleigh_tail (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case  9:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_levy (rlab_gsl_rng_r[irng], param[1], param[2]);
      break;

    case 10:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_levy_skew (rlab_gsl_rng_r[irng],
                                 param[1], param[2], param[3]);
      break;

    case 11:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_gamma (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 12: case 13:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_flat (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 14:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_lognormal (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 15:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_chisq (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 16:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_fdist (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 17:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_tdist (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 18:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_beta (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 19:
      for (i = 0; i < n; i++)
        w[i] = param[0]
            + gsl_ran_logistic (rlab_gsl_rng_r[irng], param[1]);
      break;

    case 20:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_pareto (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 21:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_weibull (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 22:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_gumbel1 (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 23:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_gumbel2 (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 24:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_poisson (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 25:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_binomial (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 26:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_negative_binomial (rlab_gsl_rng_r[irng], param[0], param[1]);
      break;

    case 27:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_geometric (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 28:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_hypergeometric (rlab_gsl_rng_r[irng],
                                       param[0], param[1], param[2]);
      break;

    case 29:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_logarithmic (rlab_gsl_rng_r[irng], param[0]);
  }
  return;
}

//
// fill-in function: submitted an MDI
//
void
gsl_random_number_generator_int (int *nr, int *nc, unsigned int *w, int *idx)
{
  int i, n = (*nr) * (*nc);
  int idrng;

  if(!idx)
    idrng = *idx;
  else
    idrng = RLAB_RNG_DEFAULT;

  //
  // setup the default rng: uniform on [0,1)
  //
  int     idist                   = RLAB_RNG_UNIFORM_IDX;
  int     nparam                  = RLAB_RNG_UNIFORM_NPARAM;
  double  dparam[RLAB_RNG_UNIFORM_NPARAM]
        = RLAB_RNG_UNIFORM_PARAMS;
  double *param=0;
  int     irng;

  // .idist==0 by default, meaning 'none'
  if (!rlab_bltin_rng[idrng].idist)
  {
    // setup default rng
    rlab_setup_rng (&idrng, &idist,
                    &nparam, dparam, &RLAB_IRNG_DEFAULT);
  }

  // get the data for the generator
  idist  = rlab_bltin_rng[idrng].idist;
  param  = rlab_bltin_rng[idrng].param;
  irng   = rlab_bltin_rng[idrng].irng;

  switch (idist)
  {
    case 24:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_poisson (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 25:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_binomial (rlab_gsl_rng_r[irng],
                                 param[0], param[1]);
      break;

    case 26:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_negative_binomial (rlab_gsl_rng_r[irng],
                                          param[0], param[1]);
      break;

    case 27:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_geometric (rlab_gsl_rng_r[irng], param[0]);
      break;

    case 28:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_hypergeometric (rlab_gsl_rng_r[irng],
                                       param[0], param[1], param[2]);
      break;

    case 29:
      for (i = 0; i < n; i++)
        w[i] = gsl_ran_logarithmic (rlab_gsl_rng_r[irng], param[0]);
      break;

    default:
      printf("error: requested integers from RNG where the distribution\n"
             "error: is continuous. Zero matrix returned!\n");
      for (i = 0; i < n; i++)
        w[i] = 0;
  }
  return;
}

//
// sample normal(0,1) distribution: user can configure IRNG only
// subroutine
//

// normal scalar
double
gslrngnf_ ( void )
{
  int irng = RLAB_IRNG_NORMAL - 1;
  int idx  = RLAB_RNG_NORMAL - 1;

  if (rlab_bltin_rng[idx].idist != RLAB_RNG_NORMAL_IDX)
  {
    //
    // setup NORMAL RNG and its IRNG
    //
    int     idist                   = RLAB_RNG_NORMAL_IDX;
    int     nparam                  = RLAB_RNG_NORMAL_NPARAM;
    double  dparam[RLAB_RNG_NORMAL_NPARAM]
    = RLAB_RNG_NORMAL_PARAMS;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }
  else
    irng = rlab_bltin_rng[idx].irng;

  return gsl_rng_uniform (rlab_gsl_rng_r[irng]);
}

// normal matrix
int
gslrngns_ (int *nr, int *nc, double *w)
{
  int irng = RLAB_IRNG_NORMAL - 1;
  int idx  = RLAB_RNG_NORMAL - 1;
  int i, j;

  if (rlab_bltin_rng[idx].idist != RLAB_RNG_NORMAL_IDX)
  {
    //
    // setup NORMAL RNG and its IRNG
    //
    int     idist                   = RLAB_RNG_NORMAL_IDX;
    int     nparam                  = RLAB_RNG_NORMAL_NPARAM;
    double  dparam[RLAB_RNG_NORMAL_NPARAM]
          = RLAB_RNG_NORMAL_PARAMS;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }

  if (w)
  {
    for (i=0; i<*nr; i++)
      for (j=0; j<*nc; j++)
        w[i +(*nr)*j]= gsl_ran_gaussian (rlab_gsl_rng_r[irng], 1.0);
    return 0;
  }

  return 1;
}

//
// sample uniform distribution on [0,1]: user can configure IRNG only
// subroutine
//

// uniform scalar
double gslrnguf_ ( void )
{
  int idx  = RLAB_RNG_UNIFORM  - 1;
  int irng = RLAB_IRNG_UNIFORM - 1;

  if (rlab_bltin_rng[idx].idist != RLAB_RNG_NORMAL_IDX)
  {
    //
    // setup UNIFORM RNG and its IRNG
    //
    int     idist                   = RLAB_RNG_UNIFORM_IDX;
    int     nparam                  = RLAB_RNG_UNIFORM_NPARAM;
    double  dparam[RLAB_RNG_UNIFORM_NPARAM]
          = RLAB_RNG_UNIFORM_PARAMS;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }
  else
    irng = rlab_bltin_rng[idx].irng;


  return gsl_rng_uniform (rlab_gsl_rng_r[irng]);
}


// uniform matrix
int
gslrngus_ (int *nr, int *nc, double *w)
{
  int irng = RLAB_IRNG_UNIFORM - 1;
  int idx  = RLAB_RNG_UNIFORM - 1;
  int i, j;

  if (rlab_bltin_rng[idx].idist != RLAB_RNG_UNIFORM_IDX)
  {
    //
    // setup UNIFORM RNG and its IRNG
    //
    int     idist                   = RLAB_RNG_UNIFORM_IDX;
    int     nparam                  = RLAB_RNG_UNIFORM_NPARAM;
    double  dparam[RLAB_RNG_UNIFORM_NPARAM]
          = RLAB_RNG_UNIFORM_PARAMS;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }
  else
    irng = rlab_bltin_rng[idx].irng;

  if (w)
  {
    for (i=0; i<*nr; i++)
      for (j=0; j<*nc; j++)
        w[i +(*nr)*j]= gsl_rng_uniform (rlab_gsl_rng_r[irng]);

    return 0;
  }

  return 1;
}


//
// sample default distribution
// subroutine (fortran)
//
int gslrngs_ (int *nr, int *nc, double *w)
{
  int idx = RLAB_RNG_DEFAULT;

  if (rlab_bltin_rng[idx].idist<1)
  {
    int     idist                   = RLAB_RNG_UNIFORM_IDX;
    int     nparam                  = RLAB_RNG_UNIFORM_NPARAM;
    double  dparam[RLAB_RNG_UNIFORM_NPARAM]
          = RLAB_RNG_UNIFORM_PARAMS;
    int     irng                    = RLAB_IRNG_DEFAULT;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }

  gsl_random_number_generator (nr, nc, w, &idx);
  return 0;
}

//
// sample default distribution
// function
//
double gslrngf_ (void)
{
  int idx = RLAB_RNG_DEFAULT;
  int nr=1, nc=1;
  double rval;

  if (rlab_bltin_rng[idx].idist<1)
  {
    int     idist                   = RLAB_RNG_UNIFORM_IDX;
    int     nparam                  = RLAB_RNG_UNIFORM_NPARAM;
    double  dparam[RLAB_RNG_UNIFORM_NPARAM]
          = RLAB_RNG_UNIFORM_PARAMS;
    int     irng                    = RLAB_IRNG_DEFAULT;
    rlab_setup_rng (&idx, &idist, &nparam, dparam, &irng);
  }

  gsl_random_number_generator (&nr, &nc, &rval, &idx);
  return rval;
}

//
// binomial distribution
//
//
int gslrngb_ ( int n, double p )
{
  rlab_setup_default_gsl_irng(RLAB_IRNG_STAT);
  return gsl_ran_binomial (rlab_gsl_rng_r[RLAB_IRNG_STAT], p, n);
}

Ent *
ent_gsl_rng (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *p;
  int i, j, irng, idist, npar, idx;
  char *s1;

  //
  // no arguments: return a single rng
  //
  if (nargs != 4 && nargs!=1 && nargs!=0 && nargs!=3)
    rerror("rng: none, one or four arguments required");

  if (nargs==0)
  {
    fprintf(stdout, "Currently defined random number generators :\n");
    for (i=0; i<RLAB_RNG_MAXNO; i++)
    {
      if (rlab_bltin_rng[i].idist>0)
      {
        if (i < RLAB_RNG_MAXUSERNO)
        {
          if (RLAB_RNG_DEFAULT == i)
            fprintf(stdout, " user * ");
          else
            fprintf(stdout, " user   ");
        }
        else
          fprintf(stdout, " system ");

        fprintf(stdout, "%3i: %s, with irng no. %i\n", i+1,
                bltin_rng[rlab_bltin_rng[i].idist].name,
                rlab_bltin_rng[i].irng+1
               );
        fprintf(stdout, "  nparam = %1i\n", rlab_bltin_rng[i].nparam);
        for (j=0; j<rlab_bltin_rng[i].nparam;j++)
          fprintf(stdout, "  p[%i] = %g\n", j+1, rlab_bltin_rng[i].param[j]);
      }
    }
    printf("(*) is a user-default generator.\n");

    return ent_Create_Rlab_Double(RLAB_RNG_DEFAULT+1.0);
  }

  //
  // first argument: the index of the RNG
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    rerror("rng: improper first argument");
  idx = (int) class_double (e1) - 1;
  if (idx < 0 || idx > RLAB_RNG_MAXNO-1)
    rerror("rng: improper first argument");

  if (nargs==1)
  {
    // only one argument: index of default rng
    if (idx < 0 || idx > RLAB_RNG_MAXUSERNO-1)
      rerror("rng: improper first argument");

    if (rlab_bltin_rng[idx].idist>0)
      RLAB_RNG_DEFAULT = idx;
    else
      rerror("rng: non-existing RNG cannot be default");

    ent_Clean (e1);
    return ent_Create_Rlab_Double(RLAB_RNG_DEFAULT+1.0);
  }


  //
  // second argument: name of the distribution
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
    rerror ("rng: improper second argument");

  s1 = class_char_pointer(e2);
  if (!s1)
    rerror ("rand: improper first argument 'name'");
  if (strlen(s1) == 0)
    rerror ("rand: improper first argument 'name'");

  // if given name is 'help' or '?' print available
  // distributions
  if (!strcmp (s1, "help") || !strcmp (s1, "?"))
  {
    fprintf
        (stdout, "rand: acceptable values for distribution (first argument):");
    for (i = 0; bltin_rng[i].index; i++)
    {
      if (!(i%5))
        fprintf (stdout, "\nrand: ");
      else
        fprintf (stdout, "'%s', ", bltin_rng[i].name);
    }
    fprintf (stdout, "please check the manual.\n");
    fflush (stdout);
    rerror ("improper first argument");
  }

  // check the name against list of names
  idist = 0;
  npar  = 0;
  for (i = 0; bltin_rng[i].index; i++)
  {
    if (!strcmp (s1, bltin_rng[i].name))
    {
      idist = bltin_rng[i].index;
      npar  = bltin_rng[i].nparam;
      break;
    }
  }

  //
  // third argument: parameters of the distribution
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("rng: improper third argument");
  p = class_matrix_real (e3);

  if (MNR(p)*MNC(p) < npar)
    rerror ("rng: improper third argument");

  //
  // fourth argument: index of the IRNG to be used
  //
  if (nargs==4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("rng: improper fourth argument");
    irng = (int) class_double (e4) - 1;
  }
  if (irng < 0 || irng > RLAB_IRNG_MAXNO-1)
    irng = RLAB_IRNG_DEFAULT;

  // does user try to setup the default rng?
  if (idx < RLAB_RNG_MAXUSERNO)
    RLAB_RNG_DEFAULT = idx;

  if (idx != RLAB_RNG_UNIFORM-1)
    rlab_setup_rng (&idx, &idist, &npar, MDRPTR(p), &irng);

  // cleanup
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}


Ent *
ent_gsl_rand (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1=0, *e2=0;
  MDR *w;
  int nr = 1, nc = 1;

  //
  // no arguments: return a single rng
  //
  if (nargs == 0)
  {
    nr = 1;
    nc = 1;
  }
  else if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = MNR (x1);
      nc = MNC (x1);
     }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = MNR (x1);
      nc = MNC (x1);
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = MNR (x1);
      nc = MNC (x1);
    }
    else
      rerror ("improper first argument");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) == MATRIX_DENSE_REAL
        && ent_type (e2) == MATRIX_DENSE_REAL)
    {
      nr = class_double (e1);
      nc = class_double (e2);
    }
  }

  nr = nr < 1 ? 1 : nr;
  nc = nc < 1 ? 1 : nc;

  //
  // call the default RNG to do thy bidding
  //
  w = mdr_Create (nr, nc);
  gsl_random_number_generator (&nr, &nc, MDRPTR(w), &RLAB_RNG_DEFAULT);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


Ent *
ent_gsl_urand (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1 = 0, *e2 = 0;
  MDR *w;
  int nr = 1, nc = 1;

  //
  // no arguments: return a single rng
  //
  if (nargs == 0)
  {
    nr = 1;
    nc = 1;
  }
  else if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else
      rerror ("improper first argument");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) == MATRIX_DENSE_REAL
        && ent_type (e2) == MATRIX_DENSE_REAL)
    {
      nr = class_double (e1);
      nc = class_double (e2);
    }
  }

  nr = nr < 1 ? 1 : nr;
  nc = nc < 1 ? 1 : nc;

  //
  // call the system uniform RNG to do thy bidding
  //
  w = mdr_Create (nr, nc);
  gslrngus_ (&nr, &nc, MDRPTR(w));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_gauss (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1=0, *e2=0;
  MDR *w;
  int nr = 1, nc = 1;

  //
  // no arguments: return a single rng
  //
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else
      rerror ("improper first argument");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) == MATRIX_DENSE_REAL
      && ent_type (e2) == MATRIX_DENSE_REAL)
    {
      nr = class_double (e1);
      nc = class_double (e2);
    }
  }

  nr = nr < 1 ? 1 : nr;
  nc = nc < 1 ? 1 : nc;

  //
  // call the system uniform RNG to do thy bidding
  //
  w = mdr_Create (nr, nc);
  gslrngns_ (&nr, &nc, MDRPTR(w));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

// *******************************************************************
//
// randomize (fill-in)
//
// *******************************************************************

Ent *
ent_gsl_randomize (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0, *eo=0;
  ListNode *node=0;
  ListNode *var=0;

  // Long live FORTRAN!
  extern int dlatmr_ ( int *,  int *, char *,  int *, char *, double *,  int *, double *, double *,
                       char *, char *, double *,  int *, double *, double *,  int *,
                       double *, char *,  int *,  int *,  int *, double *, double *,
                       char *, double *,  int *,  int *,  int * );
  extern int dlatms_ ( int *, int *, char *, int *, char *, double *, int *, double *, double *,
                       int *, int *, char *, double *, int *, double *, int * );
  extern int zlatmr_ ( int *, int *, char *, int *, char *, Complex *, int *, double *, Complex *,
                       char *, char *, Complex *, int *, double *, Complex *, int *,
                       double *, char *, int *, int *, int *, double *, double *,
                       char *, Complex *, int *, int *, int * );
  extern int zlatms_ ( int *, int *, char *, int *, char *, Complex *, int *, double *, Complex *,
                       int *, int *, char *, Complex *, int *, Complex *, int * );
  // Ah brings back the memories.

  if (nargs < 1 || nargs > 3)
    rerror("randomize: one, two or three arguments required");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  var = (ListNode *) (args[0].u.ptr);
  if (args[0].type != VAR)
    rerror("\nrandomize: first argument has to be variable and not expression!");
  e1 = var_ent (var);
  if (e1->refc > 1)
  {
    ent_DecRef (e1);
    e1 = ent_Duplicate (e1);
    listNode_AttachEnt (var, e1);
  }

  if (ent_type (e1) == MATRIX_DENSE_REAL && nargs == 1)
  {
    //
    // use default RNG to fill-in the values
    //
    MDR * x1 = ent_data (e1);
    int nr = x1->nrow;
    int nc = x1->ncol;

    if (x1->type == RLAB_TYPE_INT32)
      gsl_random_number_generator_int ( &nr, &nc, (unsigned int *) MDIPTR(x1), &RLAB_RNG_DEFAULT );
    else
      gsl_random_number_generator     ( &nr, &nc, MDRPTR(x1), &RLAB_RNG_DEFAULT );

  }
  else if (ent_type (e1) == MATRIX_DENSE_REAL && nargs > 1)
  {
    //
    // d is given: assume user wants to create a random matrix of desired
    // properties (d is either eigen- or singular values)
    //
    MDR * x1 = ent_data (e1);
    MDR * d  = 0;
    int idummy;
    double ddummy;
    int nr = x1->nrow;
    int nc = x1->ncol;
    int iseed[4];
    char sym  = 'N';
    double sparsity = 0.0;
    int kl = nr - 1;
    int ku = nc - 1;
    int mode = 6;

    int isv = 0; //

    if (nargs == 3 || nargs == 2)
    {
      eo = bltin_get_ent (args[nargs-1]);
      if (ent_type (eo) == BTREE)
      {
        // bandwidth_upper
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_BWUP);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy >= 1 && idummy <= ku)
            ku = idummy;
        }
        // bandwidth_lower
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_BWLO);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy >= 1 && idummy <= kl)
            kl = idummy;
        }
        // sparsity
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SPAR);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0 && ddummy <= 1)
            sparsity = ddummy;
        }
        // symmetric
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SYMM);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            sym = 's';
        }
        // hermitean
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_HERM);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            sym = 's';
        }
        // singular_values in d
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SVAL);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            isv = 1;
        }
      }
      if (eo)
        if (ent_type(eo)!=UNDEF)
          ent_Clean (eo);
    }

    if (nargs >= 2)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == MATRIX_DENSE_REAL)
      {
        d = class_matrix_real (e2);
        if (MNR(d) * MNC(d) == MIN(nr, nc))
          mode = 0;
      }
    }

    if (mode != 0 && isv == 0)
    {
      // no singular_values requested
      d = mdr_Create(1, MIN(nr, nc));
    }
    else if (mode != 0 && isv == 1)
    {
      // singular_values requested but not given: complain
      fprintf (stdout, "randomize: singular_values = 1 but not given");
      d = mdr_Create(1, MIN(nr, nc));
    }

    double cond = 0;
    double dmax = 0;
    char rsign = 'F';
    char grade = 'N';
    char pivtng = 'N';
    char dist = 'S';

    double anorm = -1;
    char pack = 'N';
    int info;

    if (! isv)
    {
      // find general matrix with (given) diagonal elements
      dlatmr_ ( &nr, &nc, &dist, iseed, &sym, MDRPTR(d), &mode, &cond,
                 &dmax, &rsign, &grade, NULL, NULL, NULL, NULL, NULL,
                 NULL, &pivtng, NULL, &kl, &ku, &sparsity, &anorm,
                 &pack, MDRPTR(x1), &nr, NULL, &info );
    }
    else
    {
      MDR * work = mdr_Create(1, 3 * MAX(nr, nc));

      // find general matrix with given singular (eigen) values
      dlatms_ ( &nr, &nc, &dist, iseed, &sym, MDRPTR(d), &mode, &cond,
                 &dmax, &kl, &ku, &pack, MDRPTR(x1), &nr, MDRPTR(work), &info );

      mdr_Destroy (work);
    }

    if (mode != 0)
      mdr_Destroy (d);
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC * x1 = ent_data (e1);
    MDC * dc = 0;
    MDR * d  = 0;
    int idummy;
    double ddummy;
    int nr = x1->nrow;
    int nc = x1->ncol;
    char dist = 'n';
    int iseed[4];
    char sym  = 'n';
    double sparsity = 0.0;
    int kl = nr - 1;
    int ku = nc - 1;
    int mode = 6;

    int isv = 0; //

    if (nargs == 3 || nargs == 2)
    {
      eo = bltin_get_ent (args[nargs-1]);
      if (ent_type (eo) == BTREE)
      {
        // bandwidth_upper
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_BWUP);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy >= 1 && idummy <= ku)
            ku = idummy;
        }
        // bandwidth_lower
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_BWLO);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy >= 1 && idummy <= kl)
            kl = idummy;
        }
        // sparsity
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SPAR);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0 && ddummy <= 1)
            sparsity = ddummy;
        }
        // symmetric
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SYMM);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            sym = 's';
        }
        // hermitean
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_HERM);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            sym = 'h';
        }
        // singular_values in d
        node = btree_FindNode (ent_data (eo), RLAB_NAME_RAND_SVAL);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 1)
            isv = 1;
        }
      }
      if (eo)
        if (ent_type(eo)!=UNDEF)
          ent_Clean (eo);
    }

    if (nargs >= 2)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == MATRIX_DENSE_REAL)
      {
        d = class_matrix_real (e2);
        if (MNR(d) * MNC(d) == MIN(nr, nc))
        {
          mode = 0;
          dc = (MDC *) mdr_coerce_mdc (d);
        }
      }
      else if (ent_type (e2) == MATRIX_DENSE_COMPLEX)
      {
        dc = mdc_Copy (ent_data (e2));
        if (MNR(dc) * MNC(dc) == MIN(nr, nc))
          mode = 0;
      }
    }

    if (mode != 0 && isv == 0)
    {
      // no singular_values requested
      dc = mdc_Create(1, MIN(nr, nc));
    }
    else if (mode != 0 && isv == 1)
    {
      // singular_values requested but not given: complain
      fprintf (stdout, "randomize: singular_values = 1 but not given");
      dc = mdc_Create(1, MIN(nr, nc));
    }

    double cond = 0;
    Complex dmax;
    char rsign = 'F';
    char grade = 'N';
    char pivtng = 'N';

    double anorm = -1;
    char pack = 'N';
    int info;

    if (! isv)
    {
      // find general matrix with (given) diagonal elements
      zlatmr_ ( &nr, &nc, &dist, iseed, &sym, MDCPTR(dc), &mode, &cond,
                 &dmax, &rsign, &grade, NULL, NULL, NULL, NULL, NULL,
                 NULL, &pivtng, NULL, &kl, &ku, &sparsity, &anorm,
                 &pack, MDCPTR(x1), &nr, NULL, &info );
    }
    else
    {
      MDC * work = mdc_Create(1, 3 * MAX(nr, nc));

      // find general matrix with given singular (eigen) values
      zlatms_ ( &nr, &nc, &dist, iseed, &sym, MDCPTR(dc), &mode, &cond,
                 &dmax, &kl, &ku, &pack, MDCPTR(x1), &nr, MDCPTR(work), &info );

      mdc_Destroy (work);
    }

    mdc_Destroy (dc);
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}


// *******************************************************************
//
// shuffle and sample a row-wise data matrix
//
// *******************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "shuffle"
Ent *
ent_gsl_shuffle (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1=0, *e2=0, *rent=0;
  MD  *x=0;
  MDR *wr=0;
  MDS *ws=0;
  MDC *wc=0;
  int i, j, nr, nc, ns=0, *idx=0;
  size_t n=0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // shuffling uses its own IRNG:
  int irng = RLAB_IRNG_SHUFFLE - 1;

  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);

  if (nargs!=1 && nargs!=2)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Shuffle vector, or rows of data matrix.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   shuffle(x/,n/) , with 'x' any string or numeric matrix.\n");
    rerror ("No parameters given!");
  }

  gsl_set_error_handler_off ();

  //
  // get data matrix x:
  //
  e1 = bltin_get_ent (args[0]);
  if (!isdensematrix(e1))
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MD);
    goto _exit_shuffle;
  }
  x  = ent_data(e1);
  nr = MNR(x);
  nc = MNC(x);

  // special cases:
  //  zero size matrix or scalar
  if (nr*nc==1 || nr*nc==0)
  {
    rent = ent_Copy(e1);
    goto _exit_shuffle;
  }

  if (EQVECT(x))
    n = nr * nc;
  else
    n = nr;

  idx = (int *) GC_malloc(n * sizeof(int));
  for (i = 0; i < n; i++)
    idx[i] = i;

  // get number of samples
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_REAL)
      ns = (int) class_double(e2);
  }
  if (ns<1 || ns>n)
    ns = n;

  // shuffle the index matrix
  gsl_ran_shuffle (rlab_gsl_rng_r[irng], idx, n, sizeof (int));


  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    //
    // real or integer matrix
    //
    if (nr==1 || nc==1)
    {
      switch (x->type)
      {
        case RLAB_TYPE_INT32:
          wr = mdi_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+(nc==1));
          for (i=0; i<ns; i++)
            MdiV0 (wr, i) = MdiV0 (x, idx[i]);
          break;

        case RLAB_TYPE_DOUBLE:
          wr = mdr_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+ (nc==1));
          for (i=0; i < ns; i++)
            MdrV0 (wr, i) = MdrV0 (x, idx[i]);
          break;

        default:
          wr = mdr_Create(0,0);
          break;
      }
    }
    else
    {
      switch (x->type)
      {
        case RLAB_TYPE_INT32:
          wr = mdi_Create (ns, nc);
          for (i = 0; i < ns; i++)
            for (j = 0; j < nc; j++)
              Mdi0 (wr, i, j) = Mdi0 (x, idx[i], j);
            break;

        case RLAB_TYPE_DOUBLE:
          wr = mdr_Create (ns, nc);
          for (i = 0; i < ns; i++)
            for (j = 0; j < nc; j++)
              Mdr0 (wr, i, j) = Mdr0 (x, idx[i], j);
            break;

        default:
          wr = mdr_Create(0,0);
          break;
      }
    }

    rent = ent_Assign_Rlab_MDR(wr);
    goto _exit_shuffle;
  }
  else if (ent_type (e1) == MATRIX_DENSE_STRING)
  {
    //
    // string matrix
    //
    if (nr==1 || nc==1)
    {
      ws = mds_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+ (nc==1));
      for (i=0; i<ns; i++)
        MdsV0 (ws, i) = cpstr(MdsV0 (x, idx[i]));
    }
    else
    {
      ws = mds_Create (ns, nc);
      for (i = 0; i < ns; i++)
        for (j = 0; j < nc; j++)
          Mds0 (ws, i, j) = cpstr (Mds0 (x, idx[i], j));
    }

    rent = ent_Assign_Rlab_MDS(ws);
    goto _exit_shuffle;
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    //
    // complex matrix
    //
    if (nr==1 || nc==1)
    {
      wc = mdc_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+ (nc==1));
      for (i=0; i < ns; i++)
        MdcV0 (wc, i) = MdcV0 (x, idx[i]);
    }
    else
    {
      wc = mdc_Create (ns, nc);
      for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
          Mdc0 (wc, i, j) = Mdc0 (x, idx[i], j);
    }
    rent = ent_Assign_Rlab_MDC (wc);
  }


_exit_shuffle:

  ent_Clean (e1);
  ent_Clean (e2);
  if (idx)
    GC_free (idx);

  return rent;
}


Ent *
ent_gsl_sample (int nargs, Datum args[])
{
  //  e1  - row-wise data matrix
  //  e2  - size of sample
  Ent *e1=0, *e2=0;
  MDR *xr=0, *wr=0;
  MDS *xs=0, *ws=0;
  MDC *xc=0, *wc=0;
  int i, j, nr, nc, ns = 0, *idx = 0, *sidx = 0, n;

  // sampling uses its own IRNG:
  int irng = RLAB_IRNG_SAMPLE - 1;

  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);

  if (nargs != 2)
  {
    printf ("sample: Sample with replacement vector, or rows of data matrix.\n");
    printf ("sample: Format:\n");
    printf ("sample:   sample(x,n) , where 'n' is the size of sample, and\n");
    printf ("sample: 'x' is a matrix of row-data.\n");
    rerror ("No parameters given!");
  }

  //
  // size of the sample
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("sample: 'n' has to be an integer.");
  ns = class_double (e2);

  if (ns<1)
    rerror ("sample: 'n' has to be positive integer.");

  e1 = bltin_get_ent (args[0]);
  if (    (ent_type(e1) != MATRIX_DENSE_REAL) && (ent_type(e1) != MATRIX_DENSE_STRING)
      &&  (ent_type(e1) != MATRIX_DENSE_COMPLEX)    )
    rerror ("sample: parameter does not exist.");

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    xr = ent_data (e1);
    nr = xr->nrow;
    nc = xr->ncol;

    sidx = (int *) GC_malloc(ns * sizeof(int));;

    if (nr==1 || nc==1)
      n = nr * nc;
    else
      n = nr;

    idx  = (int *) GC_malloc(n * sizeof(int));
    for (i=0; i<n; i++)
        idx[i] = i;

    gsl_ran_sample (rlab_gsl_rng_r[irng], sidx, ns, idx, n, sizeof (int));

    if (nr==1 || nc==1)
    {
      switch (xr->type)
      {
        case RLAB_TYPE_INT32:
          wr = mdi_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+(nc==1));
          for (i=0; i<ns; i++)
            MdiV0 (wr, i) = MdiV0 (xr, sidx[i]);
          break;

        case RLAB_TYPE_DOUBLE:
          wr = mdr_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+ (nc==1));
          for (i=0; i < ns; i++)
            MdrV0 (wr, i) = MdrV0 (xr, sidx[i]);
          break;

        default:
          wr = mdr_Create(0,0);
          break;
      }
    }
    else
    {
      switch (xr->type)
      {
        case RLAB_TYPE_INT32:
          wr = mdi_Create (ns, nc);
          for (i = 0; i < ns; i++)
            for (j = 0; j < nc; j++)
              Mdi0 (wr, i, j) = Mdi0 (xr, sidx[i], j);
          break;

        case RLAB_TYPE_DOUBLE:
          wr = mdr_Create (ns, nc);
          for (i = 0; i < ns; i++)
            for (j = 0; j < nc; j++)
              Mdr0 (wr, i, j) = Mdr0 (xr, sidx[i], j);
          break;

        default:
          wr = mdr_Create(0,0);
          break;
      }
    }

    ent_Clean (e1);
    ent_Clean (e2);
    GC_free (idx);
    GC_free (sidx);

    return ent_Assign_Rlab_MDR (wr);
  }
  else if (ent_type (e1) == MATRIX_DENSE_STRING)
  {
    xs = ent_data (e1);
    nr = xs->nrow;
    nc = xs->ncol;

    sidx = (int *) GC_malloc(ns * sizeof(int));;

    if (nr==1 || nc==1)
      n = nr * nc;
    else
      n = nr;

    idx  = (int *) GC_malloc(n * sizeof(int));
    for (i=0; i<n; i++)
      idx[i] = i;

    gsl_ran_sample (rlab_gsl_rng_r[irng], sidx, ns, idx, n, sizeof (int));

    if (nr==1 || nc==1)
    {
      ws = mds_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+(nc==1));
      for (i=0; i<ns; i++)
        MdsV0 (ws, i) = cpstr(MdsV0 (xs, sidx[i]));
    }
    else
    {
      ws = mds_Create (ns, nc);
      for (i = 0; i < ns; i++)
        for (j = 0; j < nc; j++)
          Mds0 (ws, i, j) = cpstr (Mds0 (xs, sidx[i], j));
    }

    ent_Clean (e1);
    ent_Clean (e2);
    GC_free (idx);
    GC_free (sidx);

    return ent_Assign_Rlab_MDS(ws);
  }

  //
  // shuffle rows of complex matrix
  //
  xc = ent_data (e1);
  nr = xc->nrow;
  nc = xc->ncol;

  sidx = (int *) GC_malloc(ns * sizeof(int));;

  if (nr==1 || nc==1)
    n = nr * nc;
  else
    n = nr;

  idx  = (int *) GC_malloc(n * sizeof(int));
  for (i=0; i<n; i++)
    idx[i] = i;

  gsl_ran_sample (rlab_gsl_rng_r[irng], sidx, ns, idx, n, sizeof (int));

  if (nr==1 || nc==1)
  {
    wc = mdc_Create ((nr==1)+ns*(nc==1), ns*(nr==1)+(nc==1));
    for (i=0; i<ns; i++)
    {
      MdcV0r (wc, i) = MdcV0r (xc, sidx[i]);
      MdcV0i (wc, i) = MdcV0i (xc, sidx[i]);
    }
  }
  else
  {
    wc = mdc_Create (ns, nc);
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        Mdc0r (wc, i, j) = Mdc0r (xc, sidx[i], j);
        Mdc0i (wc, i, j) = Mdc0i (xc, sidx[i], j);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);
  GC_free (idx);
  GC_free (sidx);

  return ent_Assign_Rlab_MDC(wc);
}


// *******************************************************
//
// discrete random number generators
//
// *******************************************************

#define RLAB_GSL_DRNG_SIZE 32

static MDS *drng_string  [RLAB_GSL_DRNG_SIZE] = {0};
static MDR *drng_real    [RLAB_GSL_DRNG_SIZE] = {0};
static MDR *drng_weights [RLAB_GSL_DRNG_SIZE] = {0};
static MDC *drng_complex [RLAB_GSL_DRNG_SIZE] = {0};
static int drng_nr       [RLAB_GSL_DRNG_SIZE] = {0};
static int drng_idx      [RLAB_GSL_DRNG_SIZE] = {0};
static int drng_irng_idx [RLAB_GSL_DRNG_SIZE] = {0};
static int last_gen = 0;
static gsl_ran_discrete_t *drng_g[RLAB_GSL_DRNG_SIZE] = {0};

Ent *
ent_gsl_setdrng (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int nw = 0, igen, i;
  if (nargs !=0 && nargs !=1 && nargs != 3 && nargs != 4)
    rerror ("drng: three or four parameters required");

  if (nargs == 0)
  {
    // list all DRNGs currently in use
    fprintf(stdout, "Currently defined discrete random number generators :\n");
    for (i=0; i<RLAB_GSL_DRNG_SIZE; i++)
    {
      if (drng_g[i])
      {
        fprintf(stdout, "user");
        if (i==last_gen)
          fprintf (stdout, " * ");
        else
          fprintf (stdout, "   ");
        fprintf(stdout, "%3i: size= %i, with irng no. %3i,",
                i+1,
                drng_nr[i],
                drng_irng_idx[i]+1
               );
        if (drng_real[i])
          fprintf(stdout," type real\n");
        else if (drng_string[i])
          fprintf(stdout," type string\n");
        else if (drng_complex[i])
          fprintf(stdout," type complex\n");
      }
    }
    fprintf(stdout, "(*) is a default generator.\n");

    return ent_Create_Rlab_Success();
  }
  else if (nargs == 1)
  {
    // only set the default DRNG
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      igen = (int) class_double (e1);
      if (igen >= 1 && igen <= RLAB_GSL_DRNG_SIZE)
        last_gen = igen - 1;
      else
        rerror ("drng: DRNG index out of range");
    }

    ent_Clean (e1);

    return ent_Create_Rlab_Double(last_gen+1);
  }

  //
  // get the index of the generator and clean up its arrays
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    igen = (int) class_double (e1);
    if (igen >= 1 && igen <= RLAB_GSL_DRNG_SIZE)
      last_gen = igen - 1;
    else
      rerror ("drng: DRNG index out of range");
  }

  // previous initiation
  if (drng_weights[last_gen])
  {
    mdr_Destroy (drng_weights[last_gen]);
    drng_weights[last_gen] = 0;
  }
  if (drng_g[last_gen])
  {
    gsl_ran_discrete_free(drng_g[last_gen]);
    drng_g[last_gen] = 0;
  }
  if (drng_real[last_gen])
  {
    mdr_Destroy (drng_real[last_gen]);
    drng_real[last_gen] = 0;
  }
  if (drng_complex[last_gen])
  {
    mdc_Destroy (drng_complex[last_gen]);
    drng_complex[last_gen] = 0;
  }
  if (drng_string[last_gen])
  {
    mds_Destroy (drng_string[last_gen]);
    drng_string[last_gen] = 0;
  }
  drng_idx[last_gen] = 0;

  //
  // get the weights
  //
  e2 = bltin_get_ent (args[2]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("setdrng: third argument 'wgt' is real vector");

  // set weights
  drng_weights[last_gen] = mdr_Float_BF ( ent_data (e2) );
  nw = MNR (drng_weights[last_gen]) * MNC (drng_weights[last_gen]);
  if (MNC (drng_weights[last_gen]) != 1 && MNR (drng_weights[last_gen]) != 1)
    rerror ("setdrng: third argument 'wgt' is real vector");

  //
  // get the values
  //
  e3 = bltin_get_ent (args[1]);
  if (ent_type (e3) == MATRIX_DENSE_REAL)
  {
    drng_real[last_gen] = mdr_Copy ( ent_data(e3) );
    drng_nr[last_gen]   = MNR (drng_real[last_gen]) * MNC (drng_real[last_gen]);
    if (drng_nr[last_gen] != nw)
    {
      // failure: clean up everything
      mdr_Destroy (drng_weights[last_gen]);
      mdr_Destroy (drng_real[last_gen]);
      drng_real[last_gen]    = 0;
      drng_weights[last_gen] = 0;
      drng_idx[last_gen]     = 0;
      rerror ("setdrng: Sizes of values and weights should not differ!");
    }
    drng_idx[last_gen] = 2;               // for reals
  }
  else if (ent_type (e3) == MATRIX_DENSE_COMPLEX)
  {
    drng_complex[last_gen] = mdc_Copy (ent_data (e3));
    drng_nr[last_gen] = MNR (drng_complex[last_gen]) * MNC (drng_complex[last_gen]);
    if (drng_nr[last_gen] != nw)
    {
      // failure: clean up everything
      mdr_Destroy (drng_weights[last_gen]);
      mdc_Destroy (drng_complex[last_gen]);
      drng_complex[last_gen] = 0;
      drng_weights[last_gen] = 0;
      drng_idx[last_gen]     = 0;
      rerror ("setdrng: Sizes of values and weights should not differ!");
    }
    drng_idx[last_gen] = 3;               // for complex
  }
  else if (ent_type (e3) == MATRIX_DENSE_STRING)
  {
    drng_string[last_gen] = mds_Copy (ent_data (e3));
    drng_nr[last_gen] = MNR (drng_string[last_gen]) * MNC (drng_string[last_gen]);
    if (drng_nr[last_gen] != nw)
    {
      // failure: clean up everything
      mdr_Destroy (drng_weights[last_gen]);
      mds_Destroy (drng_string[last_gen]);
      drng_string[last_gen]  = 0;
      drng_weights[last_gen] = 0;
      drng_idx[last_gen]     = 0;
      rerror ("setdrng: Sizes of values and weights should not differ!");
    }
    drng_idx[last_gen] = 1;               // for strings
  }
  else
    rerror ("setdrng: parameter does not exist.");

  //
  // index of IRNG to be used with respective DRNG
  //
  drng_irng_idx[last_gen] = RLAB_IRNG_DRNG-1;
  if (nargs==4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
    {
      int irng = (int) class_double (e4);
      if (irng > 0 && irng <= RLAB_IRNG_MAXNO)
        drng_irng_idx[last_gen] = irng-1;
    }
  }

  drng_g[last_gen] = gsl_ran_discrete_preproc
      (drng_nr[last_gen], MDRPTR(drng_weights[last_gen]));

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}


Ent *
ent_gsl_drand (int nargs, Datum args[])
{
  // call parameters:
  //  e1 - nr
  //  e2 - nc
  Ent *e1=0, *e2=0;
  MDR *wr;
  int i, j, nr=1, nc=1, i1, igen = last_gen;

  //
  // each DRNG uses its own IRNG:
  //
  int irng = drng_irng_idx[last_gen];
  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);

  //
  // figure out nr and nc
  //
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else
      rerror ("drand: improper second argument 'x'\n");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("drand: improper second argument 'nr'\n");
    nr = class_double (e1);


    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("drng: improper third argument 'nc'\n");
    nc = class_double (e2);
  }

  ent_Clean (e1);
  ent_Clean (e2);

  nr = nr < 1 ? 1 : nr;
  nc = nc < 1 ? 1 : nc;

  if (drng_idx[igen] == 1 && drng_string[igen])
  {
    // string values
    MDS *ws = mds_Create (nr, nc);
    for (i = 0; i < nr; i++)
      for (j = 0; j < nc; j++)
        Mds0 (ws, i, j) =
            cpstr (MdsV0 (drng_string[igen],
                   gsl_ran_discrete (rlab_gsl_rng_r[irng], drng_g[igen])));

    return ent_Assign_Rlab_MDS(ws);
  }
  if (drng_idx[igen] == 2 && drng_real[igen])
  {
    // real values
    if (drng_real[igen]->type == RLAB_TYPE_INT32)
    {
      // integers
      wr = mdi_Create (nr, nc);
      for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
          Mdi0 (wr, i, j) =
              MdiV0 (drng_real[igen], gsl_ran_discrete (rlab_gsl_rng_r[irng], drng_g[igen]));
    }
    else
    {
      // real
      wr = mdr_Create (nr, nc);
      for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
          Mdr0 (wr, i, j) =
              MdrV0 (drng_real[igen], gsl_ran_discrete (rlab_gsl_rng_r[irng],
                     drng_g[igen]));
    }
    return ent_Assign_Rlab_MDR(wr);
  }
  if (drng_idx[igen] == 3 && drng_complex[igen])
  {
    // complex values
    MDC *wc = mdc_Create (nr, nc);
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        i1 = gsl_ran_discrete (rlab_gsl_rng_r[irng], drng_g[igen]);
        Mdc1r (wc, i, j) = Mdc1r (drng_complex[igen], i1 + 1, 1);
        Mdc1i (wc, i, j) = Mdc1i (drng_complex[igen], i1 + 1, 1);
      }

    return ent_Assign_Rlab_MDC(wc);
  }

  return ent_Assign_Rlab_MDR( NULL);
}


// *******************************************************
//
// histogrammatic random number generators
//
// *******************************************************

#define RLAB_GSL_HRNG_SIZE 32

static int hrng_irng_idx [RLAB_GSL_HRNG_SIZE] = {0};
static gsl_histogram_pdf * hrng_g [RLAB_GSL_HRNG_SIZE] = {0};
static int last_hrng = 0;

Ent *
ent_gsl_sethrng (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int igen, i;

  int nhist, nbins;
  MDR *bin=0, *range=0;

  gsl_histogram *hist=0;

  // count parameters
  if (nargs !=0 && nargs !=1 &&nargs != 3 && nargs != 4)
    rerror ("hrng: one, three or four parameters required");

  if (nargs == 0)
  {
    // list all HRNGs currently in use
    fprintf(stdout, "Currently defined histogrammic random number generators :\n");
    for (i=0; i<RLAB_GSL_HRNG_SIZE; i++)
    {
      if (hrng_g[i])
      {
        fprintf(stdout, "user");
        if (i==last_hrng)
          fprintf (stdout, " * ");
        else
          fprintf (stdout, "   ");
        fprintf(stdout, "%3i: with irng no. %3i,", i+1,
                hrng_irng_idx[i]+1);
      }
    }
    fprintf(stdout, "(*) is a default generator.\n");

    return ent_Create_Rlab_Success();
  }
  else if (nargs == 1)
  {
    // only set the default HRNG
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      igen = (int) class_double (e1);
      if (igen >= 1 && igen <= RLAB_GSL_HRNG_SIZE)
        last_hrng = igen - 1;
      else
        rerror ("hrng: HRNG index out of range");
    }

    ent_Clean (e1);
    return ent_Create_Rlab_Double(last_hrng+1);
  }

  //
  // 1: index of the generator
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    igen = (int) class_double (e1);
    if (igen >= 1 && igen <= RLAB_GSL_HRNG_SIZE)
      last_hrng = igen - 1;
    else
      rerror ("hrng: HRNG index out of range");
  }

  // clean up previous initiation
  if (hrng_g[last_hrng])
  {
    gsl_histogram_pdf_free(hrng_g[last_hrng]);
    hrng_g[last_hrng] = 0;
  }

  //
  // 2: range of the histogram
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("hrng: second argument is not 'range' of histogram");
  range = class_matrix_real (e2);
  // get the size
  nbins = MNR(range)*MNC(range) - 1;

  //
  // 3: bins the histogram
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("hrng: third argument is not 'bin' of histogram");
  // get the values
  bin = class_matrix_real (e3);

  // get the numbers of histograms
  if (MNR(bin)==1 || MNC(bin)==1)
    nhist = 1;
  else
    nhist = MNC(bin);

  // some basic checks
  if (nhist > 1)
    rerror ("hrng: third argument is not single histogram");
  if (nbins != MNR(bin)*MNC(bin))
    rerror ("hrng: 'bin' and 'range' are not proper for histogram");

  // set up a histogram
  hist = gsl_histogram_alloc(nbins);
  hist->range = MDRPTR(range);
  hist->bin   = MDRPTR(bin);

  // allocate a hrng and set it up using histogram
  hrng_g[last_hrng] = gsl_histogram_pdf_alloc(nbins);
  if(gsl_histogram_pdf_init(hrng_g[last_hrng], hist) == GSL_EDOM)
  {
    // error: possibly negative bins
    hist->range = 0;
    hist->bin   = 0;
    gsl_histogram_free(hist);
    rerror ("hrng: 'bin' and 'range' are not proper for histogram");
  }

  //
  // index of IRNG to be used with respective HRNG
  //
  hrng_irng_idx[last_hrng] = RLAB_IRNG_HRNG-1;
  if (nargs==4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
    {
      int irng = (int) class_double (e4);
      if (irng > 0 && irng <= RLAB_IRNG_MAXNO)
        hrng_irng_idx[last_hrng] = irng-1;
    }
  }

  // Cleanup
  hist->range = 0;
  hist->bin   = 0;

  gsl_histogram_free(hist);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}


Ent *
ent_gsl_hrand (int nargs, Datum args[])
{
  // call parameters:
  //  e1 - nr
  //  e2 - nc
  Ent *e1 = 0, *e2 = 0;
  MDR *w=0;

  int i, j, nr=1, nc=1;
  double r;

  //
  // each DRNG uses its own IRNG:
  //
  int irng = drng_irng_idx[last_hrng];
  if (!rlab_gsl_rng_r[irng])
    rlab_setup_default_gsl_irng (irng);

  //
  // figure out nr and nc
  //
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      nr = x1->nrow;
      nc = x1->ncol;
    }
    else
      rerror ("hrand: improper second argument 'x'\n");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("hrand: improper second argument 'nr'\n");
    nr = class_double (e1);

    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("hrand: improper third argument 'nc'\n");
    nc = class_double (e2);


  }

  nr = nr < 1 ? 1 : nr;
  nc = nc < 1 ? 1 : nc;

  w = mdr_Create(nr, nc);

  if (hrng_g[last_hrng])
  {
    // hrng exists
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
    {
      r = gsl_rng_uniform(rlab_gsl_rng_r[irng]);
      Mdr0 (w, i, j) = gsl_histogram_pdf_sample(hrng_g[last_hrng], r);
    }
  }
  else
    rerror("hrand: sampling nonexisting histogrammic RNG");

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}



//
// load other rng-dependent files
//
#include "gsl_rng_siman.c"
#include "gsl_rng_nintegrate_mc.c"
#include "gsl_rng_deoptim.c"

