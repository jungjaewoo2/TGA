// mmio.c

/* This file is a part of RLaB+rlaplus
   Copyright (C) 2007-2008  Marijan Kostrun

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
 * Contains parts of the code from:
 *   Matrix Market I/O library for ANSI C
 *   See http://math.nist.gov/MatrixMarket for details.
 * modifications:
 *   malloc -> GC_malloc, as GC is official garbage collector of rlabplus
 */
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"
#include "listnode.h"
#include "mdc.h"
#include "msr.h"
#include "msc.h"


#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "mmio.h"

extern MSR *mdr_real_spconvert (MDR * m);
extern MSC *mdr_complex_spconvert (MDR * m);

// static int
// mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
//                            double **val_, int **I_, int **J_)
// {
//   FILE *f;
//     MM_typecode matcode;
//     int M, N, nz;
//     int i;
//     double *val=0;
//     int *I=0, *J=0;
//
//     if ((f = fopen(fname, "r")) == NULL)
//             return -1;
//
//
//     if (mm_read_banner(f, &matcode) != 0)
//     {
//         printf("mm_read_unsymetric: Could not process Matrix Market banner ");
//         printf(" in file [%s]\n", fname);
//         return -1;
//     }
//
//
//
//     if ( !(mm_is_real(matcode) && mm_is_matrix(matcode) &&
//             mm_is_sparse(matcode)))
//     {
//         fprintf(stderr, "Sorry, this application does not support ");
//         fprintf(stderr, "Market Market type: [%s]\n",
//                 mm_typecode_to_str(matcode));
//         return -1;
//     }
//
//     /* find out size of sparse matrix: M, N, nz .... */
//
//     if (mm_read_mtx_crd_size(f, &M, &N, &nz) !=0)
//     {
//         fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
//         return -1;
//     }
//
//     *M_ = M;
//     *N_ = N;
//     *nz_ = nz;
//
//     /* reseve memory for matrices */
//
//     val = (double *) GC_malloc(nz * sizeof(double));
//     I = (int *) GC_malloc(nz * sizeof(int));
//     J = (int *) GC_malloc(nz * sizeof(int));
//
//     *val_ = val;
//     *I_ = I;
//     *J_ = J;
//
//     /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
//     /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
//     /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
//
//     for (i=0; i<nz; i++)
//     {
//         fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
//         I[i]--;  /* adjust from 1-based to 0-based */
//         J[i]--;
//     }
//     fclose(f);
//
//     return 0;
// }

// static int
// mm_is_valid(MM_typecode matcode)
// {
//     if (!mm_is_matrix(matcode)) return 0;
//     if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
//     if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
//     if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
//                 mm_is_skew(matcode))) return 0;
//     return 1;
// }

static int
mm_read_banner(FILE *f, MM_typecode *matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;


    mm_clear_typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
        storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
            storgae) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
            mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    return 0;
}

static int
mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
    if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

static int
mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;

    else
    do
    {
        num_items_read = fscanf(f, "%d %d %d", M, N, nz);
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);

    return 0;
}


static int
mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    /* set return null parameter values, in case we exit with errors */
    *M = *N = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;

    else /* we have a blank line */
    do
    {
        num_items_read = fscanf(f, "%d %d", M, N);
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);

    return 0;
}

static int
mm_write_mtx_array_size(FILE *f, int M, int N)
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

// static int
// mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
//                      double val[], MM_typecode matcode)
// {
//     int i;
//     if (mm_is_complex(matcode))
//     {
//         for (i=0; i<nz; i++)
//             if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
//                 != 4) return MM_PREMATURE_EOF;
//     }
//     else if (mm_is_real(matcode))
//     {
//         for (i=0; i<nz; i++)
//         {
//             if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
//                 != 3) return MM_PREMATURE_EOF;
//
//         }
//     }
//
//     else if (mm_is_pattern(matcode))
//     {
//         for (i=0; i<nz; i++)
//             if (fscanf(f, "%d %d", &I[i], &J[i])
//                 != 2) return MM_PREMATURE_EOF;
//     }
//     else
//         return MM_UNSUPPORTED_TYPE;
//
//     return 0;
//
// }

// static int
// mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
//                       double *real, double *imag, MM_typecode matcode)
// {
//     if (mm_is_complex(matcode))
//     {
//             if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
//                 != 4) return MM_PREMATURE_EOF;
//     }
//     else if (mm_is_real(matcode))
//     {
//             if (fscanf(f, "%d %d %lg\n", I, J, real)
//                 != 3) return MM_PREMATURE_EOF;
//
//     }
//
//     else if (mm_is_pattern(matcode))
//     {
//             if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
//     }
//     else
//         return MM_UNSUPPORTED_TYPE;
//
//     return 0;
//
// }


/************************************************************************
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
************************************************************************/

// static int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I,
//                            int **J, double **val, MM_typecode *matcode)
// {
//     int ret_code;
//     FILE *f;
//
//     if (strcmp(fname, "stdin") == 0) f=stdin;
//     else
//     if ((f = fopen(fname, "r")) == NULL)
//         return MM_COULD_NOT_READ_FILE;
//
//
//     if ((ret_code = mm_read_banner(f, matcode)) != 0)
//         return ret_code;
//
//     if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
//             mm_is_matrix(*matcode)))
//         return MM_UNSUPPORTED_TYPE;
//
//     if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
//         return ret_code;
//
//
//     *I = (int *)  GC_malloc(*nz * sizeof(int));
//     *J = (int *)  GC_malloc(*nz * sizeof(int));
//     *val = NULL;
//
//     if (mm_is_complex(*matcode))
//     {
//         *val = (double *) GC_malloc(*nz * 2 * sizeof(double));
//         ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
//                 *matcode);
//         if (ret_code != 0) return ret_code;
//     }
//     else if (mm_is_real(*matcode))
//     {
//         *val = (double *) GC_malloc(*nz * sizeof(double));
//         ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
//                 *matcode);
//         if (ret_code != 0) return ret_code;
//     }
//
//     else if (mm_is_pattern(*matcode))
//     {
//         ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
//                 *matcode);
//         if (ret_code != 0) return ret_code;
//     }
//
//     if (f != stdin) fclose(f);
//     return 0;
// }

static int
mm_write_banner(FILE *f, MM_typecode matcode)
{
    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    GC_free(str);
    if (ret_code !=2 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

// static int
// mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
//                  double val[], MM_typecode matcode)
// {
//     FILE *f;
//     int i;
//
//     if (strcmp(fname, "stdout") == 0)
//         f = stdout;
//     else
//     if ((f = fopen(fname, "w")) == NULL)
//         return MM_COULD_NOT_WRITE_FILE;
//
//     /* print banner followed by typecode */
//     fprintf(f, "%s ", MatrixMarketBanner);
//     fprintf(f, "%s\n", mm_typecode_to_str(matcode));
//
//     /* print matrix sizes and nonzeros */
//     fprintf(f, "%d %d %d\n", M, N, nz);
//
//     /* print values */
//     if (mm_is_pattern(matcode))
//         for (i=0; i<nz; i++)
//             fprintf(f, "%d %d\n", I[i], J[i]);
//     else
//     if (mm_is_real(matcode))
//         for (i=0; i<nz; i++)
//             fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
//     else
//     if (mm_is_complex(matcode))
//         for (i=0; i<nz; i++)
//             fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i],
//                         val[2*i+1]);
//     else
//     {
//         if (f != stdout) fclose(f);
//         return MM_UNSUPPORTED_TYPE;
//     }
//
//     if (f !=stdout) fclose(f);
//
//     return 0;
// }


/**
*  Create a new copy of a string s.  strdup() is a common routine, but
*  not part of ANSI C, so it is included here.  Used by mm_typecode_to_str().
*
*/
static char *mystrdup(const char *s)
{
	int len = strlen(s);
	char *s2 = (char *) GC_malloc((len+1)*sizeof(char));
	return strcpy(s2, s);
}

static char
* mm_typecode_to_str(MM_typecode matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
	  //char *mystrdup(const char *);

    /* check for MTX type */
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return mystrdup(buffer);

}

//
// readmm: read matrix in the matrix market format (nist)
//
Ent *
ent_readmm (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  char *fname;
  int i, j, k, p, nr=0, nc=0;
  double dd, di;

  FILE *f=0;

  int iold=0;
  MM_typecode matcode;
  int m, n, nz;
  long int id;

  //
  // matrices to be read in
  //
  MDR *drm=0;
  MDC *dcm=0;
  MSR *srm=0;
  MSC *scm=0;

  ListNode *var;

  // Check n_args
  if (nargs != 1 && nargs != 2)
    rerror ("readmm: one or two arguments required!");

  // filename
  e1 = bltin_get_ent (args[0]);
  fname = class_char_pointer (e1);
  if (!fname)
    rerror ("readmm: empty filename");
  if (!strlen(fname))
    rerror ("readmm: empty filename");

  // use the memory storage of VAR2 to store
  // what ever you are reading in
  if (nargs == 2)
  {
    // check reference count for e2. if more then 1 duplicate it
    // and reassign it to VAR2
    var = (ListNode *) (args[1].u.ptr);
    e2 = var_ent (var);
    if (e2->refc > 1)
    {
      ent_DecRef (e2);
      e2 = ent_Duplicate (e2);
      listNode_AttachEnt (var, e2);
    }
    // check the type of e2:
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      drm = (MDR *) ent_data (e2);
      nc  = drm->ncol;
      nr  = drm->nrow;
      iold = 1;
    }
    else if (ent_type (e2) == MATRIX_DENSE_COMPLEX)
    {
      dcm = (MDC *) ent_data (e2);
      nr  = dcm->nrow;
      nc  = dcm->ncol;
      iold = 1;
    }
    else if (ent_type (e2) == MATRIX_SPARSE_REAL)
    {
      srm = (MSR *) ent_data (e2);
      nr  = srm->nr;
      nc  = srm->nc;
      iold = 1;
    }
    if (ent_type (e2) == MATRIX_SPARSE_COMPLEX)
    {
      scm = (MSC *) ent_data (e2);
      nr  = scm->nr;
      nc  = scm->nc;
      iold = 1;
    }
  }

  // open the file. do not use rlab's built-in file system
  if ((f = fopen(fname, "r")) == NULL)
    rerror ("readmm: cannot open file!");

  // banner
  if (mm_read_banner(f, &matcode) != 0)
    rerror ("readmm: file not in matrix market format (no banner)!\n");

  // type of matrix
  if (!mm_is_matrix(matcode))
    rerror ("readmm: file not in matrix market format (not a matrix)!\n");

  // output entity
  rent = ent_Create ();

  if (mm_is_dense(matcode))
  {
    //
    // dense matrix: real, integer or complex
    //

    // find the size of the dense matrix: m=nrows, n=ncols
    if (mm_read_mtx_array_size(f, &m, &n) == MM_PREMATURE_EOF)
    {
      fclose (f);
      rerror ("readmm: file not in matrix market format (not a matrix)!\n");
    }

    if (mm_is_real(matcode))
    {
      // real:

      // create storage if it does not exists
      if (!drm)
        drm = mdr_Create(m,n);
      else if (m!=nr || n!=nc)
        rerror ("readmm: mismatch in size of supplied variable \n");

      // read it in: general, symmetric or skew
      if (mm_is_general(matcode))
      {
        for (j=0; j<n; j++)
          for (i=0; i<m; i++)
        {
           fscanf(f, "%lg\n", &dd);
           Mdr0(drm,i,j) = dd;
        }
      }
      else if (mm_is_symmetric(matcode))
      {
        for (j=0; j<n; j++)
          for (i=j; i<m; i++)
        {
          fscanf(f, "%lg\n", &dd);
          Mdr0(drm,i,j) = dd;
          Mdr0(drm,j,i) = dd;
        }
      }
      else if (mm_is_skew(matcode))
      {
        for (j=0; j<n; j++)
        {
          Mdr0(drm,j,j) = 0;
          for (i=j+1; i<m; i++)
          {
            fscanf(f, "%lg\n", &dd);
            Mdr0(drm,i,j) =  dd;
            Mdr0(drm,j,i) = -dd;
          }
        }
      }

      if (iold)
        ent_data (rent) = mdr_CreateScalar(1.0);
      else
        ent_data (rent) = drm;

      ent_SetType (rent, MATRIX_DENSE_REAL);
    }
    else if (mm_is_integer(matcode))
    {
      // integer:

      // create storage if it does not exists
      if (!drm)
        drm = mdi_Create(m,n);
      else if (m!=nr || n!=nc )
        rerror ("readmm: mismatch in size of supplied variable \n");

      // read it in
      if (mm_is_general(matcode))
      {
        if (drm->type == RLAB_TYPE_INT32)
        {
          for (j=0; j<n; j++)
            for (i=0; i<m; i++)
          {
            fscanf(f, "%li\n", &id);
            Mdi0(drm,i,j) = id;
          }
        }
        else
        {
          for (j=0; j<n; j++)
            for (i=0; i<m; i++)
          {
            fscanf(f, "%li\n", &id);
            Mdr0(drm,i,j) = (double) id;
          }
        }
      }
      else if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
      {
        if (drm->type == RLAB_TYPE_INT32)
        {
          for (j=0; j<n; j++)
            for (i=j; i<m; i++)
          {
            fscanf(f, "%li\n", &id);
            Mdi0(drm,i,j) = id;
            Mdi0(drm,j,i) = id;
          }
        }
        else
        {
          for (j=0; j<n; j++)
            for (i=j; i<m; i++)
          {
            fscanf(f, "%li\n", &id);
            Mdr0(drm,i,j) = (double) id;
            Mdr0(drm,j,i) = (double) id;
          }
        }
      }
      else if (mm_is_skew(matcode))
      {
        if (drm->type == RLAB_TYPE_INT32)
        {
          for (j=0; j<n; j++)
          {
            Mdi0(drm,j,j) = 0;
            for (i=j+1; i<m; i++)
            {
              fscanf(f, "%li\n", &id);
              Mdi0(drm,i,j) =  id;
              Mdi0(drm,j,i) = -id;
            }
          }
        }
        else
        {
          for (j=0; j<n; j++)
          {
            Mdr0(drm,j,j) = 0;
            for (i=j+1; i<m; i++)
            {
              fscanf(f, "%li\n", &id);
              Mdr0(drm,i,j) = (double)  id;
              Mdr0(drm,j,i) = (double) -id;
            }
          }
        }
      }


      ent_SetType (rent, MATRIX_DENSE_REAL);
      if (iold)
        ent_data (rent) = mdr_CreateScalar(1.0);
      else
        ent_data (rent) = drm;
    }
    else if (mm_is_complex(matcode))
    {
      // complex

      // create storage if it does not exists
      if (!dcm)
        dcm = (MDC*) mdc_Create(m,n);
      else if (m!=nr || n!=nc )
        rerror ("readmm: mismatch in size of supplied variable \n");

      if (mm_is_general(matcode))
      {
        for (j=0; j<n; j++)
         for (i=0; i<m; i++)
        {
           fscanf(f, "%lg %lg\n", &dd, &di);
           Mdc0r(dcm,i,j) = dd;
           Mdc0i(dcm,i,j) = di;
        }
      }
      else if (mm_is_hermitian(matcode))
      {
        for (j=0; j<n; j++)
          for (i=j; i<m; i++)
        {
          fscanf(f, "%lg %lg\n", &dd, &di);
          // [i;j]
          Mdc0r(dcm,i,j) = dd;
          Mdc0i(dcm,i,j) = di;
          // [j;i]
          Mdc0r(dcm,j,i) =  dd;
          Mdc0i(dcm,j,i) = -di;
        }
      }
      if (iold)
      {
        // old storage used: report  a ok
        ent_SetType (rent, MATRIX_DENSE_REAL);
        ent_data (rent) = mdr_CreateScalar(1.0);
      }
      else
      {
        // new storage created: return it
        ent_data (rent) = dcm;
        ent_SetType (rent, MATRIX_DENSE_COMPLEX);
      }

    }
  }
  else if (mm_is_sparse(matcode))
  {
    //
    // sparse matrix: real, integer or complex
    //

    // find the size of the dense matrix: m=nrows, n=ncols,
    // nz=no of nonzero elements in the file (that may be different from
    // the no of  nonzero elements in the matrix for symmetric and skew-symmetric
    // matrices)
    if (mm_read_mtx_crd_size(f, &m, &n, &nz) == MM_PREMATURE_EOF)
    {
      fclose (f);
      rerror ("readmm: file not in matrix market format (not a matrix)!\n");
    }

    MDR *xd=0;

    if (mm_is_real(matcode))
    {
      // general format
      if (mm_is_general(matcode))
      {
        // coordinate storage
        xd = mdr_Create(nz,3);

        // read in coordinate formatted data
        for(k=0;k<nz;k++)
        {
          fscanf(f, "%d %d %lg\n", &i, &j, &dd);
          Mdr0(xd,k,0) = i;
          Mdr0(xd,k,1) = j;
          Mdr0(xd,k,2) = dd;
        }

        // convert coordinate to MSR
        // define the matrix
        if (!srm)
        {
          srm = (MSR*) mdr_real_spconvert( xd );
          mdr_Destroy (xd);
        }

        // new storage created: return it
        ent_data (rent) = srm;
        ent_SetType (rent, MATRIX_SPARSE_REAL);

      }
      else if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
      {
        // coordinate storage
        xd = mdr_Create(2*nz,3);
        mdr_Zero(xd);

        // read in coordinate formatted data
        p = 0;
        for(k=0;k<nz;k++)
        {
          fscanf(f, "%d %d %lg\n", &i, &j, &dd);
          Mdr0(xd,p,0) = i;
          Mdr0(xd,p,1) = j;
          Mdr0(xd,p,2) = dd;
          p++;
          if (i!=j)
          {
            Mdr0(xd,p,0) = j;
            Mdr0(xd,p,1) = i;
            Mdr0(xd,p,2) = dd;
            p++;
          }
        }

        // convert coordinate to MSR
        // define the matrix
        if (!srm)
        {
          srm = (MSR*) mdr_real_spconvert( xd );
          mdr_Destroy (xd);
        }

        // new storage created: return it
        ent_data (rent) = srm;
        ent_SetType (rent, MATRIX_SPARSE_REAL);
      }

    }
    else if (mm_is_complex(matcode))
    {
      // general format
      if (mm_is_general(matcode))
      {
        // coordinate storage
        xd = mdr_Create(nz,4);

        // read in coordinate formatted data
        for(k=0;k<nz;k++)
        {
          fscanf(f, "%d %d %lg %lg\n", &i, &j, &dd, &di);
          Mdr0(xd,k,0) = i;
          Mdr0(xd,k,1) = j;
          Mdr0(xd,k,2) = dd;
          Mdr0(xd,k,3) = di;
        }

        // convert coordinate to MSR
        // define the matrix
        if (!srm)
        {
          scm = (MSC*) mdr_complex_spconvert( xd );
          mdr_Destroy (xd);
        }

        ent_data (rent) = scm;
        ent_SetType (rent, MATRIX_SPARSE_COMPLEX);
      }
      else if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
      {
        // coordinate storage
        xd = mdr_Create(2*nz,4);
        mdr_Zero(xd);

        // read in coordinate formatted data
        p = 0;
        for(k=0;k<nz;k++)
        {
          fscanf(f, "%d %d %lg %lg\n", &i, &j, &dd, &di);
          Mdr0(xd,p,0) = i;
          Mdr0(xd,p,1) = j;
          Mdr0(xd,p,2) = dd;
          Mdr0(xd,p,3) = di;
          p++;
          if (i!=j)
          {
            Mdr0(xd,p,0) = j;
            Mdr0(xd,p,1) = i;
            Mdr0(xd,p,2) = dd;
            Mdr0(xd,p,3) = di;
            p++;
          }
        }

        // convert coordinate to MSR
        // define the matrix: it ignores 0,0-entries
        if (!srm)
        {
          scm = (MSC*) mdr_complex_spconvert( xd );
          mdr_Destroy (xd);
        }

        // new storage created: return it
        ent_data (rent) = scm;
        ent_SetType (rent, MATRIX_SPARSE_COMPLEX);
      }
    }
  }

  if (f !=stdin)
    fclose(f);

  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}


//
// writemm: write matrix in the matrix market format (nist)
//
Ent *
ent_writemm (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  char *fname;
  int ia, i, j, k;

  FILE *f=0;

  MM_typecode matcode;
  int n;

  //
  // matrices to be read in
  //
  MDR *drm=0;
  MDC *dcm=0;
  MSR *srm=0;
  MSC *scm=0;

  // Check n_args
  if (nargs != 2)
    rerror ("writemm: two arguments required!");

  // filename
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
    rerror ("writemm: missing filename");
  fname = class_char_pointer (e1);
  if (!fname)
    rerror ("writemm: empty filename");
  if (!strlen(fname))
    rerror ("writemm: empty filename");

  // filename
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL    &&
      ent_type(e2)!=MATRIX_DENSE_COMPLEX &&
      ent_type(e2)!=MATRIX_SPARSE_REAL   &&
      ent_type(e2)!=MATRIX_SPARSE_COMPLEX)
    rerror ("writemm: trying to write non-matrix");

  // open the file. do not use rlab's built-in file system
  if ((f = fopen(fname, "w")) == NULL)
    rerror ("writemm: cannot open file for writing");

  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_general(&matcode);

  if (ent_type(e2)==MATRIX_DENSE_REAL)
  {
    drm = (MDR *) ent_data (e2);

    mm_set_dense(&matcode);

    if (drm->type == RLAB_TYPE_INT32)
    {
      mm_set_integer(&matcode);

      mm_write_banner(f, matcode);
      mm_write_mtx_array_size(f, drm->nrow, drm->ncol);

      for (j=0; j<drm->ncol; j++)
        for (i=0; i<drm->nrow; i++)
      {
        fprintf(f, "%i\n", Mdi0(drm,i,j));
      }
    }
    else
    {
      mm_set_real(&matcode);

      mm_write_banner(f, matcode);
      mm_write_mtx_array_size(f, drm->nrow, drm->ncol);

      for (j=0; j<drm->ncol; j++)
        for (i=0; i<drm->nrow; i++)
          fprintf(f, "%lg\n", Mdr0(drm,i,j));
    }
  }
  else if (ent_type(e2)==MATRIX_DENSE_COMPLEX)
  {
    dcm = (MDC *) ent_data (e2);

    mm_set_dense(&matcode);
    mm_set_complex(&matcode);

    mm_write_banner(f, matcode);
    mm_write_mtx_array_size(f, dcm->nrow, dcm->ncol);

    for (j=0; j<dcm->ncol; j++)
      for (i=0; i<dcm->nrow; i++)
        fprintf(f, "%lg %lg\n", Mdc0r(dcm,i,j), Mdc0i(dcm,i,j));
  }
  else if (ent_type(e2)==MATRIX_SPARSE_REAL)
  {
    srm = (MSR *) ent_data(e2);

    mm_set_sparse(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    mm_write_banner(f, matcode);
    mm_write_mtx_crd_size(f, srm->nr, srm->nc, srm->nnz);

    k = 1;
    for (i = 0; i < srm->nr; i++)
    {
      n  = srm->ia[i + 1] - srm->ia[i];
      ia = srm->ia[i];
      for (j = 0; j < n; j++)
      {
        fprintf (f, "%i %i %g\n",
                 k, srm->ja[ia + j - 1], srm->d[ia + j - 1]);
      }
      k++;
    }
  }
  else if (ent_type(e2)==MATRIX_SPARSE_COMPLEX)
  {
    scm = (MSC *) ent_data (e2);

    mm_set_sparse(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_complex(&matcode);

    mm_write_banner(f, matcode);
    mm_write_mtx_crd_size(f, scm->nr, scm->nc, scm->nnz);

    k = 1;
    for (i = 0; i < scm->nr; i++)
    {
      n  = scm->ia[i + 1] - scm->ia[i];
      ia = scm->ia[i];
      for (j = 0; j < n; j++)
      {
        fprintf (f, "%i %i %lg %lg\n",
                 k, scm->ja[ia + j - 1],
                 RE(scm->c[ia + j - 1]), IM(scm->c[ia + j - 1]));
      }
      k++;
    }
  }

  if (f !=stdin)
    fclose(f);

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);

  if (e2)
    if (ent_type(e2)!=UNDEF)
      ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar(1.0);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

