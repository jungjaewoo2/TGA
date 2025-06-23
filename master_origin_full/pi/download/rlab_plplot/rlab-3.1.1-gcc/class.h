/* class.h */

/*
 * Provide the framework for class and class-method functionality.
 */

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

#ifndef RLAB_CLASS_H
#define RLAB_CLASS_H

#include "rlab.h"
#include "ent.h"

extern char *etd (Ent * e);
extern void class_init (void);
extern void ent_data_class_copy(Ent *old, Ent *new_content);
extern Ent *class_copy (Ent * e);
extern void class_destroy (Ent * e);
extern void class_print (Ent * e);

extern Ent *class_negate (Ent * e);
extern Ent *class_add (Ent * e1, Ent * e2);
extern Ent *class_addto (Ent * e1, Ent * e2);
extern Ent *class_subtract (Ent * e1, Ent * e2);
extern Ent *class_subtractfrom (Ent * e1, Ent * e2);
extern Ent *class_multiply (Ent * e1, Ent * e2);
extern Ent *class_el_multiply (Ent * e1, Ent * e2);
extern Ent *class_elmultiplyby (Ent * e1, Ent * e2);
extern Ent *class_rdivide (Ent * e1, Ent * e2);
extern Ent *class_el_rdivide (Ent * e1, Ent * e2);
extern Ent *class_elrdivideby (Ent * e1, Ent * e2);
extern Ent *class_ldivide (Ent * e1, Ent * e2);
extern Ent *class_el_ldivide (Ent * e1, Ent * e2);
extern Ent *class_power (Ent * e1, Ent * e2);
extern Ent *class_el_power (Ent * e1, Ent * e2);

extern int class_logical_scalar (Ent * e);
extern size_t class_sizeof (Ent * e);
extern char * class_char_pointer (Ent * e);
extern double class_double (Ent * e);
extern int class_int (Ent * e);
extern Complex class_complex (Ent * e);
extern MDR *class_matrix_real (Ent * e);
extern MDR *class_matrix_double (Ent * e);
extern MDS *class_matrix_string (Ent * e);

extern int class_size (Ent * e);
extern Ent *class_append (Ent * e1, Ent * e2);
extern Ent *class_stack (Ent * e1, Ent * e2);

extern Ent *class_matrix_sub_1 (Ent * var, Ent * i, Ent * j);
extern Ent *class_matrix_sub_2 (Ent * var, Ent * i);
extern Ent *class_matrix_sub_3 (Ent * var, Ent * j);

extern Ent *class_matrix_assign_1 (Ent * var, Ent * i, Ent * j, Ent * rhs);
extern Ent *class_matrix_assign_2 (Ent * var, Ent * i, Ent * rhs);
extern Ent *class_matrix_assign_3 (Ent * var, Ent * j, Ent * rhs);

extern Ent *class_vector_sub (Ent * var, Ent * ind);
extern Ent *class_matrix_vector_assign (Ent * var, Ent * i, Ent * rhs);

extern Ent *class_forloop_value (Ent * e, int i);
extern Ent *class_empty (void);

extern void *class_member_ref (Ent * e, char *name, int *rtype);

extern Ent *class_eq (Ent * e1, Ent * e2);
extern Ent *class_ne (Ent * e1, Ent * e2);
extern Ent *class_lt (Ent * e1, Ent * e2);
extern Ent *class_le (Ent * e1, Ent * e2);
extern Ent *class_gt (Ent * e1, Ent * e2);
extern Ent *class_ge (Ent * e1, Ent * e2);
extern Ent *class_and (Ent * e1, Ent * e2);
extern Ent *class_or (Ent * e1, Ent * e2);
extern Ent *class_not (Ent * e);

extern Ent *class_transpose (Ent * e);
extern Ent *class_nc_transpose (Ent * e);
extern Ent *class_reshape_col (Ent * e);

extern Ent *class_increment (Ent * e);
extern Ent *class_decrement (Ent * e);

extern int class_attribute (Ent *e, char *name);

#endif /* RLAB_CLASS_H */
