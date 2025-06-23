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


// 
// the rng's for the rest of rlab+rlabplus
// 
extern int RLAB_IRNG_DEFAULT;
extern gsl_rng *rlab_gsl_rng_r[];

// sample dedicated uniform [0,1) RNG
// subroutine:
extern int    gslrngus_ ( int *nr, int *nc, double *w );
// function:
extern double gslrnguf_ ( void );

// sample default RNG
// subroutine
extern double gslrngs_( int *nr, int *nc, double *w );
// function
extern double gslrngf_( void );





