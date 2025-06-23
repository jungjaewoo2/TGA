/* ent.h */

/* This file is a part of RLaB ("Our"-LaB)
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

/*
 * An Rlab entity. This structure isolates the data, from the
 * variable. Since there can be many variables pointing to a
 * single entity, this struct contains a reference count. Ents
 * also contain the type info since, in the future, rlab will pick
 * the acutal type values to avoid collisions.
 */

#ifndef RLAB_ENT_H
#define RLAB_ENT_H

#include "mds.h"
#include "mdr.h"
#include "complex.h"
#include "mdc.h"

extern int  not_ent_double_vector(Ent *e);
extern int  ismde(Ent *ent);
extern int  isfuncent(Ent *ent);
extern int  isdensematrix(Ent *x);
extern Ent *ent_Create_Rlab_Success (void);
extern Ent *ent_Create_Rlab_Failure (void);
extern Ent *ent_Create_Rlab_Error   (void);
extern Ent *ent_Create_Rlab_Double (double);
extern Ent *ent_Create_Rlab_Complex (double,double);
extern Ent *ent_Create_Rlab_Int (int);
extern Ent *ent_Create_Rlab_String (char *);
extern Ent *ent_Assign_Rlab_String (char *);
extern Ent *ent_Assign_Rlab_MDR (MDR *x);
extern Ent *ent_Assign_Rlab_MDC (MDC *x);
extern Ent *ent_Assign_Rlab_MDS (MDS *x);
extern Ent *ent_Assign_Rlab_MDE (MDE *x);
extern Ent *ent_Assign_Rlab_BTREE (Btree *b);
extern Ent *ent_Assign_Rlab_Rtype (void *, int);

extern Ent  *ent_Create (void);
extern void *ent_Destroy (Ent *ent);
extern int   ent_Clean (Ent *ent);
extern void  ent_double_Destroy (Ent *ent);
extern void  ent_undef_Destroy (Ent *ent);
extern void  ent_undef_Print (Ent *ent, FILE *fptr);

extern void ent_SetType (Ent *ent, int type);

extern Ent *ent_Copy (Ent *ent);
extern Ent *ent_Duplicate (Ent *ent);

extern Ent *ent_Inc (Ent *ent);
extern Ent *ent_Dec (Ent *ent);

extern void ent_Double (Ent *ent, double d);
extern char *ent_Double_CharPointer (Ent *e);
extern char *double_Class (void *ptr);
extern char *undef_Class (void *ptr);
extern MDS *double_Type_BF (void *ptr);
extern size_t double_Sizeof (void *ptr);
extern MDS *undef_Type_BF (void *ptr);
extern void *undef_Copy (void *u);

#define ent_IncRef(e)  (((Ent *) (e))->refc++)
#define ent_DecRef(e)  (((Ent *) (e))->refc--)

#define ent_Ref(e)     (((Ent *) (e))->refc)

#define ent_type(a)    (((Ent *) (a))->type)
#define ent_double(a)  (((Ent *) (a))->d)
#define ent_data(a)    (((Ent *) (a))->data)

#endif  /* RLAB_ENT_H */
