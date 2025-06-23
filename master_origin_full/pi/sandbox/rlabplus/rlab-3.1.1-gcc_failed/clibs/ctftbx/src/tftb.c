/* tftb.c: time-frequency toolbox */

/*  This file is a part of rlabplus
   Copyright (C) 2008  Marijan Kostrun

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   See the file ./COPYING
   ********************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "gc.h"


/*   local functions   */
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#define sgn(A)          ((A) > 0.0 ? 1.0 : -1.0)
#define ABS(a)          (((a) >= (0)) ? (a) : (-a))
#define SWAP(a,b)       {temp = (a); (a)=(b); (b)=temp;}
#define ROUND1(x)       (((((x)-(int)(x))>=0) &&(((x)-(int)(x))<0.5)) ? ((int)(x)) : ((int)(x+1)))
#define ROUND(x)        ((int)(sgn((x))*ROUND1((x)*sgn((x)))))
#define ISODD(x)        ((x/2.0)== ((int)(x/2)) ? 0 : 1)


/* local constants */
#define pi               3.141592653589793
#define EPS              0.0000000001

#define TRUE             1
#define FALSE            0

extern int ctftbx_debug;

/* definition of the distance identifiers */
enum ctftbx_dist
{
  CTFTBX_UNKNOWN_DIST,
  CTFTBX_LQ,
  CTFTBX_QUADRATIC,
  CTFTBX_CORRELATION,
  CTFTBX_KOLMOGOROV,
  CTFTBX_KULLBACK,
  CTFTBX_CHERNOFF,
  CTFTBX_MATUSITA,
  CTFTBX_NLQ,
  CTFTBX_LSD,
  CTFTBX_JENSEN
};

/* definition of the kernel shapes */
enum ctftbx_kernel_shape
{
  CTFTBX_UNKNOWN_KERNEL,
  CTFTBX_MTEK,
  CTFTBX_RGK,
  CTFTBX_GMCWK,
  CTFTBX_WIGNER,
  CTFTBX_SPECTRO
};

/* parametres for the MTEK */
#define NB_PARAM_MTEK    7
#define ALPHA            parameters[0]
#define BETA             parameters[1]
#define GAMMA            parameters[2]
#define R                parameters[3]
#define TAU_0            parameters[4]
#define NU_0             parameters[5]
#define LAMBDA           parameters[6]

/*-------------------------------------------*/
/*  definition of the structures and types   */
/*-------------------------------------------*/

#include "../../../complex.h"
#include "../../../fftp.h"

// struct _rlab_complex
// {
//   double r;     /* Real part */
//   double i;     /* Imaginary part */
// };
// typedef struct _rlab_complex Complex;

/* Signal structure */
typedef struct SIG
{
    int            length;          // Length of the signal in points
    double         sample_freq;     // Sample frequency of the signal
    double        *time_instants;   // instants of sample for the signal
    unsigned char  is_complex;      // TRUE if there exists an imag part
    Complex       *signal;          // real part of the signal

}
type_signal, kernel_type_signal;

typedef struct Time_freq_rep
{
    int            N_freq;        // number of freq bins in the TFR matrix
    int            N_time;        // number of time_bins in the TFR matrix
    double        *freq_bins;     // fqs for each line of the matrix
    double        *time_instants; // instant for each column of the TFR
    unsigned char  is_complex;    // TRUE if there exists an imag part
    Complex       *tfr;           // tfr
}
type_TFR;

typedef struct Ambi_func
{
    int            N_doppler;     // number of doppler bins in the AF
    int            N_delay;	      // number of delay bins in the AF matrix
    double        *doppler_bins;  // doppler bin for each line of the AF: x-coordinate
    double        *delay_bins;    // delay bin for each column of the AF: y-coordinate
    unsigned char  is_complex;    // TRUE if there exists an imag part
    Complex       *af;            // AF
}
type_AF, kernel_type_AF;

//
// Local function
//
#include "divers.c"

//
// load the rest of the library
//
#include "af.c"
#include "af2tfr.c"
#include "bj.c"
#include "bud.c"
#include "cw.c"
#include "distance.c"
#include "grd.c"
#include "hough.c"
#include "kernel.c"
#include "mh.c"
#include "mhs.c"
#include "mmce.c"
#include "page.c"
#include "pmh.c"
#include "ppage.c"
#include "pwv.c"
#include "ri.c"
#include "ridb.c"
#include "ridbn.c"
#include "ridh.c"
#include "ridt.c"
#include "spwv.c"
#include "wv.c"
#include "zam.c"





