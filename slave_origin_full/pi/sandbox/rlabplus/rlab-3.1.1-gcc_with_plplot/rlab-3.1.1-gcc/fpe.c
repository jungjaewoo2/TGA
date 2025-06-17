/* fpe.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1994,  Ian R. Searle

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

/*
 * Set up floating point handling. What we want:
 *
 * Ignore everything, BUT overflow.
 *
 * Set rounding mode, round to nearest representable number, tie -> even 
 *
 * Now the trick is to figure out how to do this on each type of
 * hardware/OS combination.
 *
 */

#include "config.h"

#ifdef SETUP_FPE

#ifdef HAVE_FPSETMASK

/*
 * SVR3.2 of SVR4.x type floating point exception handling.
 * HP-UX for PA-RISC maybe.
 */

#ifdef HAVE_IEEEFP_H

#include <ieeefp.h>
static fp_except fpmask, fpmask_orig;
static fp_rnd fprnd, fprnd_orig;

#else

/* Act of desparation! (this works on HP-UX) */
#include <math.h>
static fp_except fpmask, fpmask_orig;
static fp_rnd fprnd, fprnd_orig;

#endif /* HAVE_IEEEFP_H */

void
setup_fpe_handling ()
{
  fpmask = fpgetmask ();

#ifdef FP_X_INV
  fpmask = fpmask | !FP_X_INV;
#endif

#ifdef FP_X_DNML
  fpmask = fpmask | !FP_X_DNML;
#endif

#ifdef FP_X_DZ
  fpmask = fpmask | !FP_X_DZ;
#endif

#ifdef FP_X_OFL
  fpmask = fpmask | !FP_X_OFL;
#endif

#ifdef FP_X_UFL
  fpmask = fpmask | !FP_X_UFL;
#endif

#ifdef FP_X_IMP
  fpmask = fpmask | !FP_X_IMP;
#endif

  fpmask_orig = fpsetmask (fpmask);

#ifdef FP_RN
  fprnd_orig = fpsetround (FP_RN);
#endif
}
#else

/*****************************************************************************/
#if defined(HAVE_IEEE_HANDLER)

/*
 * Berkeley/SunOS-4 type floating point exception handling.
 */

/*
 * This code block is NOT FINISHED YET.
 * I have not had a chance to test it, or
 * to figure out how to restore the old 
 * handlers.
 */

#ifdef HAVE_FLOATINGPOINT_H
#include <floatingpoint.h>
static sigfpe_handler_type old_handler1, old_handler2;
static sigfpe_handler_type hdl;
extern int ieee_handler ();
#include "util.h"
#endif

void
setup_fpe_handling ()
{
  ieee_handler ("get", "overflow", old_handler1);
  ieee_handler ("get", "invalid", old_handler2);

  hdl = (sigfpe_handler_type) fpecatch;
  ieee_handler ("set", "overflow", hdl);
  ieee_handler ("set", "inexact", SIGFPE_IGNORE);
  ieee_handler ("set", "division", SIGFPE_IGNORE);
  ieee_handler ("set", "underflow", SIGFPE_IGNORE);
  ieee_handler ("set", "invalid", SIGFPE_IGNORE);
}

/*****************************************************************************/
#elif defined(HAVE_IEEE_SET_FP_CONTROL)

#ifdef HAVE_MACHINE_FPU_H
#include <machine/fpu.h>
#endif

void
setup_fpe_handling ()
{
  unsigned long csr = ieee_get_fp_control ();
  csr |= IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF;
  ieee_set_fp_control (csr);
}

/*****************************************************************************/
#elif defined(HAVE___SETFPUCW)

/*
 * Linux systems.
 */

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

void
setup_fpe_handling ()
{
  __setfpucw (_FPU_IEEE);
}

/*****************************************************************************/
#elif defined(win32)

/*
 * win32 (Windows 95/NT) floating point setup
 */

#include <float.h>

void
setup_fpe_handling ()
{
  /* Link with with 80bit precision library. */
  /* 64 bit precision. */
  /* _controlfp (_PC_64, _MCW_PC); */

  /* Infinity */
  _controlfp (_IC_AFFINE, _MCW_IC);
}

/*****************************************************************************/
#else /* HAVE_IEEE_HANDLER */

void
setup_fpe_handling ()
{
  ;				/* Do nothing */
}

#endif
#endif /* HAVE_FPSETMASK */
#else /* SETUP_FPE */

void
setup_fpe_handling ()
{
  ;				/* Do nothing. */
}

#endif /* SETUP_FPE */
