/* rlab_macros.h */

/*  This file is a part of RLaB ("Our"-LaB) and rlabplus
   Copyright (C) 2014-2017 M.Kostrun

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
#ifndef RLAB_MACROS_L
#define RLAB_MACROS_L

#define RLAB_MDR_CALL_2ARGS_FUNC(init,name,descr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j;\
  double a, x; \
  MDR *x1=0, *x2=0, *w=0; \
  if (nargs!=1 && nargs!=2) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" init "(x1,x2), or " init "( [x1,x2] )\n"); \
    rerror ("One or two arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror (init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=2) \
      rerror ( init ": single argument has to be two-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = Mdr0 (x1, i, 0); \
      a = Mdr0 (x1, i, 1); \
      MdrV0 (w, i) = name (x, a); \
    } \
  } \
  else if (nargs == 2) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2) < 1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        Mdr0 (w, i, j) = name (x, a); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}



#define RLAB_MDR_CALL_3ARGS_FUNC(init,name,descr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0; \
  int nr, nc, nr2, nc2, i, j;\
  double a, b, x; \
  MDR *x1=0, *x2=0, *x3=0, *w=0; \
  if (nargs!=1 && nargs!=2 && nargs!=3) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" init "([x1,x2,x3]), or " init "( x1,[x2,x3] ), or " init "(x1,x2,x3)\n"); \
    rerror ("One, two or three arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=3) \
      rerror ( init ": single argument has to be 3-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      MdrV0 (w, i) = name (x, a, b); \
    } \
  } \
  else if (nargs == 2) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    nr2 = MNR(x2); \
    nc2 = MNC(x2); \
    if (nc2!=2) \
      rerror ( init ": second argument has to be two-column matrix"); \
    nr = MAX(nr, nr2); \
    nc = MAX(nc, nc2); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, 0); \
        b = mdr0_safe (x2, i, 1); \
        Mdr0 (w, i, j) = name (x, a, b); \
      } \
    } \
  } \
  else if (nargs == 3) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        Mdr0 (w, i, j) = name (x, a, b); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  return ent_Assign_Rlab_MDR(w); \
}



#define RLAB_MDR_CALL_4ARGS_FUNC(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0, *e4=0; \
  int nr, nc, i, j;\
  double a, b, c, x; \
  MDR *x1=0, *x2=0, *x3=0, *x4=0, *w=0; \
  if (nargs!=1 && nargs!=4) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or four arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=4) \
      rerror ( init ": single argument has to be 4-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      c = mdr0 (x1, i, 3); \
      MdrV0 (w, i) = name (x, a, b, c); \
    } \
  } \
  else if (nargs == 4) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    e4 = bltin_get_ent (args[3]); \
    if (ent_type (e4) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x4 = class_matrix_real (e4); \
    if (SIZE(x4)<1) \
      rerror ( init ": fourth argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x4)); \
    nc = MAX(nc, MNC(x4)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        c = mdr0_safe (x4, i, j); \
        Mdr0 (w, i, j) = name (x, a, b, c); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  ent_Clean (e4); \
  return ent_Assign_Rlab_MDR(w); \
} 


#define RLABENT_CALL_D_FUNC_D(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0; \
  int i; \
  MDR *x1=0, *w; \
  if (nargs != 1) \
  { \
    printf ( init ": " descr); \
    printf ( init ": Format:\n"); \
    printf ( init ": \t" calldescr); \
    rerror ( init ": One parameter required\n" ); \
  }\
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = ent_data (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  w = mdr_Create_SameSize (x1); \
  for (i=0; i<SIZE(x1); i++) \
  { \
    MdrV0 (w, i) = name (mdrV0(x1,i)); \
  } \
  ent_Clean (e1); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_C(init,name,descr,calldescr,constant) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0; \
  int i; \
  MDR *x1=0, *w; \
  if (nargs != 1) \
  { \
    printf ( init ": " descr); \
    printf ( init ": Format:\n"); \
    printf ( init ": \t" calldescr); \
    rerror ( init ": One parameter required\n" ); \
  }\
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = ent_data (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  w = mdr_Create_SameSize (x1); \
  for (i=0; i<SIZE(x1); i++) \
  { \
    MdrV0 (w, i) = name (mdrV0(x1,i),constant); \
  } \
  ent_Clean (e1); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_I_D(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j;\
  double x; \
  int n; \
  MDR *x1=0, *x2=0, *w=0; \
  if (nargs!=1 && nargs!=2) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One or two arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror (init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=2) \
      rerror ( init ": single argument has to be two-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      n = mdi0 (x1, i, 0); \
      x = mdr0 (x1, i, 1); \
      MdrV0 (w, i) = name (n,x); \
    } \
  } \
  else if (nargs == 2) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2) < 1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        n = mdi0_safe (x1, i, j); \
        x = mdr0_safe (x2, i, j); \
        Mdr0 (w, i, j) = name (n, x); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j;\
  double x,y; \
  MDR *x1=0, *x2=0, *w=0; \
  if (nargs!=1 && nargs!=2) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One or two arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror (init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=2) \
      rerror ( init ": single argument has to be two-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdi0 (x1, i, 0); \
      y = mdr0 (x1, i, 1); \
      MdrV0 (w, i) = name (x,y); \
    } \
  } \
  else if (nargs == 2) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2) < 1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        y = mdr0_safe (x2, i, j); \
        Mdr0 (w, i, j) = name (x, y); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_I_DA_UA(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j; \
  MDR *x1=0, *x2=0, *w; \
  if (nargs != 2) \
  { \
    printf ( init ": " descr); \
    printf ( init ": Format:\n"); \
    printf ( init ": \t" calldescr); \
    rerror ( init ": Two parameters required\n" ); \
  }\
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = ent_data (e1); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  e2 = bltin_get_ent (args[1]); \
  if (ent_type (e2) != MATRIX_DENSE_REAL) \
    rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
  x2 = ent_data (e2); \
  if (SIZE(x2) < 1) \
    rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
  if (nc != MNC(x2)) \
    rerror ( init ": arguments has to have same number of columns"); \
  double *p = GC_malloc(nc * sizeof(double)); \
  unsigned int *n = GC_malloc(nc * sizeof(unsigned int)); \
  nr = MAX(nr, MNR(x2)); \
  w = mdr_Create (nr, 1); \
  for (i=0; i<nr; i++) \
  { \
    for (j = 0; j < nc; j++) \
    { \
      p[j] = mdr0_safe (x1, i, j); \
      n[j] = (unsigned int) mdi0_safe (x2, i, j); \
    } \
    MdrV0 (w, i) = name (nc, p, n); \
  } \
  GC_free(n); \
  GC_free(p); \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}



#define RLABENT_CALL_D_FUNC_I_DA_DA(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j; \
  MDR *x1=0, *x2=0, *w; \
  if (nargs != 2) \
  { \
    printf ( init ": " descr); \
    printf ( init ": Format:\n"); \
    printf ( init ": \t" calldescr); \
    rerror ( init ": Two parameters required\n" ); \
  }\
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = ent_data (e1); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  e2 = bltin_get_ent (args[1]); \
  if (ent_type (e2) != MATRIX_DENSE_REAL) \
    rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
  x2 = ent_data (e2); \
  if (SIZE(x2) < 1) \
    rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
  if (nc != MNC(x2)) \
    rerror ( init ": arguments has to have same number of columns"); \
  double *p = GC_malloc(nc * sizeof(double)); \
  double *t = GC_malloc(nc * sizeof(double)); \
  nr = MAX(nr, MNR(x2)); \
  w = mdr_Create (nr, 1); \
  for (i=0; i<nr; i++) \
  { \
    for (j = 0; j < nc; j++) \
    { \
      p[j] = mdr0_safe (x1, i, j); \
      t[j] = mdr0_safe (x2, i, j); \
    } \
    MdrV0 (w, i) = name (nc, p, t); \
  } \
  GC_free(t); \
  GC_free(p); \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D_C(init,name,descr,calldescr,const) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0; \
  int nr, nc, i, j;\
  double x, y; \
  MDR *x1=0, *x2=0, *w=0; \
  if (nargs!=1 && nargs!=2) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or two arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=2) \
      rerror ( init ": single argument has to be 2-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      y = mdr0 (x1, i, 1); \
      MdrV0 (w, i) = name (x, y, const); \
    } \
  } \
  else if (nargs == 2) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be two-column matrix"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        y = mdr0_safe (x2, i, j); \
        Mdr0 (w, i, j) = name (x, y, const); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D_D(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0; \
  int nr, nc, i, j;\
  double a, b, x; \
  MDR *x1=0, *x2=0, *x3=0, *w=0; \
  if (nargs!=1 && nargs!=3) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or three arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=3) \
      rerror ( init ": single argument has to be 3-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      MdrV0 (w, i) = name (x, a, b); \
    } \
  } \
  else if (nargs == 3) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        Mdr0 (w, i, j) = name (x, a, b); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D_D_C(init,name,descr,calldescr,const) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0; \
  int nr, nc, i, j;\
  double a, b, x; \
  MDR *x1=0, *x2=0, *x3=0, *w=0; \
  if (nargs!=1 && nargs!=3) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or three arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=3) \
      rerror ( init ": single argument has to be 3-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      MdrV0 (w, i) = name (x, a, b,const); \
    } \
  } \
  else if (nargs == 3) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        Mdr0 (w, i, j) = name (x, a, b,const); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D_D_D(init,name,descr,calldescr) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0, *e4=0; \
  int nr, nc, i, j;\
  double a, b, c, x; \
  MDR *x1=0, *x2=0, *x3=0, *x4=0, *w=0; \
  if (nargs!=1 && nargs!=4) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or four arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=4) \
      rerror ( init ": single argument has to be 4-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      c = mdr0 (x1, i, 3); \
      MdrV0 (w, i) = name (x, a, b, c); \
    } \
  } \
  else if (nargs == 4) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    e4 = bltin_get_ent (args[3]); \
    if (ent_type (e4) != MATRIX_DENSE_REAL) \
      rerror ( init ": fourth argument has to be MATRIX-DENSE-REAL"); \
    x4 = class_matrix_real (e4); \
    if (SIZE(x4)<1) \
      rerror ( init ": fourth argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x4)); \
    nc = MAX(nc, MNC(x4)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        c = mdr0_safe (x4, i, j); \
        Mdr0 (w, i, j) = name (x, a, b, c); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  ent_Clean (e4); \
  return ent_Assign_Rlab_MDR(w); \
}


#define RLABENT_CALL_D_FUNC_D_D_D_D_C(init,name,descr,calldescr,const) \
Ent * ent_##name (int nargs, Datum args[]) \
{ \
  Ent *e1=0, *e2=0, *e3=0, *e4=0; \
  int nr, nc, i, j;\
  double a, b, c, x; \
  MDR *x1=0, *x2=0, *x3=0, *x4=0, *w=0; \
  if (nargs!=1 && nargs!=4) \
  { \
    printf (init ": " descr); \
    printf (init ": Format:\n"); \
    printf (init ": \t" calldescr); \
    rerror ("One, or four arguments required"); \
  } \
  e1 = bltin_get_ent (args[0]); \
  if (ent_type (e1) != MATRIX_DENSE_REAL) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  x1 = class_matrix_real (e1); \
  if (SIZE(x1) < 1) \
    rerror ( init ": first argument has to be MATRIX-DENSE-REAL"); \
  nr = MNR(x1); \
  nc = MNC(x1); \
  if (nargs == 1) \
  { \
    if (nc!=4) \
      rerror ( init ": single argument has to be 4-column MATRIX-DENSE-REAL"); \
    w = mdr_Create (nr, 1); \
    for (i = 0; i < nr; i++) \
    { \
      x = mdr0 (x1, i, 0); \
      a = mdr0 (x1, i, 1); \
      b = mdr0 (x1, i, 2); \
      c = mdr0 (x1, i, 3); \
      MdrV0 (w, i) = name (x, a, b, c, const); \
    } \
  } \
  else if (nargs == 4) \
  { \
    e2 = bltin_get_ent (args[1]); \
    if (ent_type (e2) != MATRIX_DENSE_REAL) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    x2 = class_matrix_real (e2); \
    if (SIZE(x2)<1) \
      rerror ( init ": second argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x2)); \
    nc = MAX(nc, MNC(x2)); \
    e3 = bltin_get_ent (args[2]); \
    if (ent_type (e3) != MATRIX_DENSE_REAL) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    x3 = class_matrix_real (e3); \
    if (SIZE(x3)<1) \
      rerror ( init ": third argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x3)); \
    nc = MAX(nc, MNC(x3)); \
    e4 = bltin_get_ent (args[3]); \
    if (ent_type (e4) != MATRIX_DENSE_REAL) \
      rerror ( init ": fourth argument has to be MATRIX-DENSE-REAL"); \
    x4 = class_matrix_real (e4); \
    if (SIZE(x4)<1) \
      rerror ( init ": fourth argument has to be MATRIX-DENSE-REAL"); \
    nr = MAX(nr, MNR(x4)); \
    nc = MAX(nc, MNC(x4)); \
    w = mdr_Create (nr, nc); \
    for (i = 0; i < nr; i++) \
    { \
      for (j = 0; j < nc; j++) \
      { \
        x = mdr0_safe (x1, i, j); \
        a = mdr0_safe (x2, i, j); \
        b = mdr0_safe (x3, i, j); \
        c = mdr0_safe (x4, i, j); \
        Mdr0 (w, i, j) = name (x, a, b, c, const); \
      } \
    } \
  } \
  ent_Clean (e1); \
  ent_Clean (e2); \
  ent_Clean (e3); \
  ent_Clean (e4); \
  return ent_Assign_Rlab_MDR(w); \
}

#endif