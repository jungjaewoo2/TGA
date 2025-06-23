#ifndef GSL_RLAB_H
#define GSL_RLAB_H

// rlab
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdc.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

//gsl:
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

gsl_vector_int * alloc_gsl_vector_int_mdr( MDR * x );
void free_gsl_vector_int_mdr( gsl_vector_int * v );

gsl_vector * alloc_gsl_vector_mdr( MDR * x );
void free_gsl_vector_mdr( gsl_vector * v );

gsl_matrix * alloc_gsl_matrix_mdr( MDR * x );
void free_gsl_matrix_mdr( gsl_matrix * v );

gsl_vector * gsl_vector_copy_mdr(MDR * x);  // use GSL free
gsl_matrix * gsl_matrix_copy_mdr(MDR * x);  // use GSL free

MDR * alloc_mdr_gsl_vector( gsl_vector * x );
void  free_mdr_gsl_vector( MDR *v );

MDR * copy_mdr_gsl_vector( gsl_vector * x );
MDR * copy_mdr_gsl_matrix( gsl_matrix * x );

#endif /* RLAB_MDR_H */
