//
//
//
#include "gsl_rlab.h"

//
// submit MDR as a gsl vector
//
gsl_vector * alloc_gsl_vector_mdr( MDR * x )
{
  gsl_vector * v;

  if (x->type == RLAB_TYPE_INT32)
    rerror("gsl: cannot allocate an integer vector as a double");

  v = (gsl_vector *) GC_MALLOC (sizeof (gsl_vector));

  if (!v)
    rerror("gsl: failed to allocate space for vector struct");

  v->data = MDRPTR(x);
  v->size = (x->nrow) * (x->ncol);
  v->stride = 1;
  v->block = 0;
  v->owner = 0;

  return v;
}

void free_gsl_vector_mdr( gsl_vector * v )
{
  if (!v)
    return;

  v->data = NULL;
  v->size = 0;
  v->stride = 1;
  v->block = 0;
  v->owner = 0;

  GC_free (v);
  v=0;
  return;
}


//
// copy MDR to gsl vector: use GSL free to release
//
gsl_vector * gsl_vector_copy_mdr(MDR * x)
{
  int i;
  size_t n = SIZE(x);
  gsl_vector * v=0;

  if (n<0)
  {
    fprintf(stderr, "gsl_vector_copy_mdr: null-pointer sent instead of valid MDR!\n");
    return v;
  }

  v = gsl_vector_alloc(n);

  for (i=0; i<n; i++)
    gsl_vector_set(v, i, mdrV0(x,i));

  return v;
}

//
// attach MDR to a gsl_matrix
//
gsl_matrix * alloc_gsl_matrix_mdr( MDR * x )
{
  gsl_matrix * m=0;
  if (x->type == RLAB_TYPE_INT32)
    rerror("alloc_gsl_matrix_mdr: cannot allocate an integer matrix as a double");

  m = (gsl_matrix *) GC_MALLOC (sizeof (gsl_matrix));

  if (!m)
    rerror("alloc_gsl_matrix_mdr: failed to allocate space for matrix struct");

  m->data  = MDRPTR(x);
  m->size1 = x->nrow;
  m->size2 = x->ncol;
  m->tda   = x->nrow;
  m->block = 0;
  m->owner = 0;

  return m;
}

//
// free so created gsl_matrix
//
void free_gsl_matrix_mdr( gsl_matrix * m )
{
  if (!m)
    return;

  m->data  = 0;
  m->size1 = 0;
  m->size2 = 0;
  m->tda   = 0;
  m->block = 0;
  m->owner = 0;

  GC_free (m);
  m=0;

  return;
}

//
// copy MDR to gsl vector: use GSL free to release
//
gsl_matrix * gsl_matrix_copy_mdr(MDR * x)
{
  int i,j;
  size_t n = SIZE(x);
  gsl_matrix * m=0;

  if (n<1)
  {
    fprintf(stderr, "gsl_vector_copy_mdr: null-pointer sent instead of valid MDR!\n");
    return m;
  }

  m = gsl_matrix_alloc(MNR(x), MNC(x));

  for (i=0; i<MNR(x); i++)
    for (j=0; j<MNC(x); j++)
      gsl_matrix_set(m, i, j, mdr0(x,i,j));

  return m;
}

MDR * alloc_mdr_gsl_vector( gsl_vector * x )
{
  MDR *v = 0;

  if (!x)
  {
    fprintf(stderr, "alloc_mdr_gsl_vector: null-pointer sent instead of valid gsl_vector!\n");
    return v;
  }

  v = mdr_CreateEmpty(1, x->size);
  MDPTR(v) =  x->data;
  return v;
}

void free_mdr_gsl_vector( MDR * v )
{
  if (!v)
    return;

  MNR(v) = 0;
  MNC(v) = 0;
  MDPTR(v) = 0;
  GC_free (v);
  v=0;
  return;
}

MDR * copy_mdr_gsl_vector( gsl_vector * x )
{
  MDR *v = 0;

  if (!x)
  {
    fprintf(stderr, "copy_mdr_gsl_vector: null-pointer sent instead of valid gsl_vector!\n");
    return v;
  }

  int i,n = x->size;

  if (n>0)
  {
    v = mdr_Create(1, n);
    for(i=0; i<n; i++)
      MdrV0(v,i) = gsl_vector_get(x, i);
  }

  return v;
}

MDR * copy_mdr_gsl_matrix( gsl_matrix * x )
{
  size_t i,j;
  MDR *v = 0;
  if (!x)
  {
    fprintf(stderr, "copy_mdr_gsl_matrix: null-pointer sent instead of valid gsl_vector!\n");
    return v;
  }

  v = mdr_Create(x->size1,x->size2);
  for(i=0; i<x->size1;i++)
    for(j=0; j<x->size2;j++)
      Mdr0(v,i,j) = gsl_matrix_get(x, i, j);
  return v;
}




