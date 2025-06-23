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
#include <gsl/gsl_rng.h>

//
// list of IRNG
//
#define RLAB_IRNG_MAXNO     64
#define RLAB_IRNG_MAXUSERNO 32

#define RLAB_IRNG_SIMAN   (RLAB_IRNG_MAXNO-7)
#define RLAB_IRNG_MC      (RLAB_IRNG_MAXNO-6)
#define RLAB_IRNG_DRNG    (RLAB_IRNG_MAXNO-5)
#define RLAB_IRNG_HRNG    (RLAB_IRNG_MAXNO-4)
#define RLAB_IRNG_SAMPLE  (RLAB_IRNG_MAXNO-3)
#define RLAB_IRNG_SHUFFLE (RLAB_IRNG_MAXNO-2)
#define RLAB_IRNG_NORMAL  (RLAB_IRNG_MAXNO-1)
#define RLAB_IRNG_UNIFORM (RLAB_IRNG_MAXNO)

extern int RLAB_IRNG_DEFAULT;
extern gsl_rng *rlab_gsl_rng_r[];
extern int rlab_setup_default_gsl_irng (int idx);

struct _rlab_ran_pdf
{
  int     idist;    // index of the distribution
  int     nparam;   // no. of parameters
  double *param;    // the parameter array
  int     irng;     // index of the underlying integer rng
};
extern struct _rlab_ran_pdf rlab_bltin_rng[];

//
// list of RNG
//
// size of list
#define RLAB_RNG_MAXNO     64
#define RLAB_RNG_MAXUSERNO 32
#define RLAB_RNG_ASA_UNIFORM (RLAB_RNG_MAXNO-1)
// default generator: uniform in [0,1]
#define RLAB_RNG_UNIFORM_IDX    12
#define RLAB_RNG_UNIFORM_NPARAM 2
#define RLAB_RNG_UNIFORM_PARAMS {0,1}

// normal (0,1) generator
#define RLAB_RNG_NORMAL_IDX    1
#define RLAB_RNG_NORMAL_NPARAM 2
#define RLAB_RNG_NORMAL_PARAMS {0,1}

// default generator: uniform in [0,1]
extern int  RLAB_RNG_DEFAULT;
extern int  rlab_setup_rng (int *idrng, int *idist, int *nparam, double *param, int *irng);

// fill-in
extern void gsl_random_number_generator (int *nr, int *nc, double *w, int *idx);

// mostly fortran
extern int    gslrngns_ (int *nr, int *nc, double *w);    // sample system U[0,1] distribution: matrix, or
extern double gslrngnf_ ( void );                         //   single value
extern int    gslrngus_ (int *nr, int *nc, double *w);    // sample system U[0,1] distribution: matrix, or
extern double gslrnguf_ ( void );                         //   single value
extern int    gslrngs_  (int *nr, int *nc, double *w);    // sample user-default distribution: matrix, or
extern double gslrngf_  (void);                           //   single value

