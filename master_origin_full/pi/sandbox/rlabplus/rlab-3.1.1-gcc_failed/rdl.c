/* rdl.c */
/* rdl.c: RLaB Dynamic linking */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

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

#include "rlab.h"
#include "symbol.h"
#include "bltin.h"
#include "class.h"
#include "ent.h"
#include "util.h"
#include "mem.h"

#include <stdio.h>

#ifdef HAVE_SO			/* Don't do anything if Shared Objects not supported. */

typedef Ent *(*EFPTR) ();

/* **************************************************************
 * Completely separate the different DLopen functions for each 
 * platform, since there will be differences for most platforms.
 * ************************************************************** */

#ifdef HAVE_DLOPEN

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#ifdef HAVE_DL_H
#include <dl.h>
#endif

Ent *
DLopen (int nargs, Datum args[])
{
  char *fn, *fun_name;
  const char *errmsg;
  Ent *e1, *e2, *rent;
  Bltin *newb;
#ifdef HAVE_DL_H
  shl_t handle;
#else
  void *handle;
#endif
  void *fptr = 0;

  if (nargs != 2)
    rerror ("dlopen: requires 2 arguments");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  fn = class_char_pointer (e1);
  fun_name = class_char_pointer (e2);

#ifdef RTLD_NOW
  handle = dlopen (fn, RTLD_NOW);
#else
#ifdef DL_HP
  handle = dlopen (fn, BIND_DEFERRED, 0L);
#else
  handle = dlopen (fn, RTLD_LAZY);
#endif
#endif

  if (handle == 0)
  {
    ent_Clean (e1);
    ent_Clean (e2);
#ifndef DL_HP
    fprintf (stderr, "%s\n", dlerror ());
#else
    fprintf (stderr, "%d\n", errno);
#endif
    rerror ("dlopen: dlopen error while opening file");
  }

#ifndef DL_HP
  fptr = dlsym (handle, fun_name);
  if ((errmsg = dlerror ()) != 0)
#else
  if (dlsym (&handle, fun_name, TYPE_UNDEFINED, &fptr))
#endif
  {
    ent_Clean (e1);
    ent_Clean (e2);
#ifndef DL_HP
    fprintf (stderr, "%s\n", errmsg);
#else
    fprintf (stderr, "%d\n", errno);
#endif
    rerror ("dlopen: dlsym error while obtaining handle");
  }

  /*
   * Now create and return the built-in function.
   */

  newb = (Bltin *) GC_MALLOC (sizeof (Bltin));
  newb->type = BLTIN;
  newb->name = 0;
  newb->func = (EFPTR) fptr;

  rent = ent_Create ();
  ent_SetType (rent, BLTIN);
  ent_data (rent) = newb;

  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}

#endif /* HAVE_DLOPEN */

/* **************************************************************
 * OS2 Version of dynamic linking...
 * ************************************************************** */

#ifdef OS2

#define INCL_DOSERRORS
#define INCL_DOSMODULEMGR
#include <os2.h>

void
DLopen (return_ptr, n_args, d_arg)
     VPTR *return_ptr;
     int n_args;
     Datum *d_arg;
{
  char *errmsg, *fn, *fun_name;
  Bltin *newb;
  ListNode *FN, *FUNC;
  HMODULE handle;
  PFN fptr;

  if (n_args != 2)
    error_1 ("dlopen: requires 2 arguments", 0);

  FN = bltin_get_string ("dlopen", d_arg, 1);
  FUNC = bltin_get_string ("dlopen", d_arg, 2);

  fn = string_GetString (e_data (FN));
  fun_name = string_GetString (e_data (FUNC));

  {
    ULONG ordinal = 0;
    UCHAR errstr[128];
    UCHAR errmsg[256];
    APIRET rc;

    rc = DosLoadModule (errstr, sizeof (errstr), fn, &handle);

    if (rc != NO_ERROR)
    {
      remove_tmp_destroy (FN);
      remove_tmp_destroy (FUNC);
      sprintf (errmsg, "OS2 API return code = %ld (%s)", rc, errstr);
      error_1 ("DosLoadModule", errmsg);
    }

    rc = DosQueryProcAddr (handle, ordinal, fun_name, &fptr);

    if (rc != NO_ERROR)
    {
      remove_tmp_destroy (FN);
      remove_tmp_destroy (FUNC);
      sprintf (errmsg, "OS2 API return code = %ld", rc);
      error_1 ("DosQueryProcAddr", errmsg);
    }

  }

  /*
   * Now create and return the built-in function.
   */

  newb = (Bltin *) MALLOC (sizeof (Bltin));
  newb->type = BLTIN;
  newb->name = 0;
  newb->func = (EFPTR) fptr;

  *return_ptr = (VPTR) newb;

  remove_tmp_destroy (FN);
  remove_tmp_destroy (FUNC);
}

#endif /* OS2 */
#endif /* HAVE_SO */
