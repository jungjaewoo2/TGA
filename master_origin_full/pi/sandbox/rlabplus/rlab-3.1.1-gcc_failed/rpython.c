// Copyright (C) 2003-2010 Marijan Kostrun
//   part of rlabplus for linux project on rlabplus.sourceforge.net
//
// linux-gpib wrapper for rlabplus
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

//
// rlab headers
//
#include "rlab.h"
#include "mdr.h"
#include "mdc.h"
#include "mdr_mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "complex.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"

//
// standard headers
//
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef HAVE_PYTHON

#include "python2.7/Python.h"

#define _THIS_LIB "rlabplus_python.c"

static int
    lib_debug=0;

static char *def_module = "__main__";

enum PYTHON_RLAB_TYPES
{
  PYTHON_NONE,
  PYTHON_INT,
  PYTHON_LONG,
  PYTHON_FLOAT,
  PYTHON_COMPLEX,
  PYTHON_STRING
};

enum PYTHON_RLAB_TYPES
    _py_object_is_scalar (PyObject * py_v)
{
  if (PyInt_Check(py_v))
    return PYTHON_INT;
  else if(PyLong_Check(py_v))
    return PYTHON_LONG;
  else if(PyFloat_Check(py_v))
    return PYTHON_FLOAT;
  else if (PyComplex_Check(py_v))
    return PYTHON_COMPLEX;
  else if (PyString_Check(py_v))
    return PYTHON_STRING;
  else
    return PYTHON_NONE;
}

enum PYTHON_RLAB_TYPES
    _py_object_is_vector (PyObject * py_v, int *nr, int *nc)
{
  enum PYTHON_RLAB_TYPES t1, t2;
  PyObject *f1, *f2;

  *nr = 0;
  *nc = 0;

  int i, isame;

  if ( PyList_Check(py_v) )
  {
    *nc = PyList_Size( py_v );

    if (!(*nc))
      return PYTHON_NONE;

    *nr = 1;

    // get info about the first element in the list
    f1 = PyList_GetItem (py_v,0);
    t1 = _py_object_is_scalar(f1);

    // does it have only one element
    if ((*nc) == 1)
      return t1;

    // it has many elements, are they all of the same type as the first one?
    isame = 1;
    for (i=1; i<(*nc) && isame; i++)
    {
      f2 = PyList_GetItem (py_v,i);
      t2 = _py_object_is_scalar(f2);
      if (t1 != t2)
        isame = 0;
    }

    if (isame)
      return t1;
    else
    {
      *nc = 0;
      *nr = 0;
      return PYTHON_NONE;
    }
  }
  else
    return PYTHON_NONE;
}

enum PYTHON_RLAB_TYPES
    _py_object_is_matrix (PyObject * py_v, int *nr, int *nc)
{
  enum PYTHON_RLAB_TYPES t1, t2;
  PyObject *f1, *f2;

  int nc1=0, nr1=0;
  int nc2=0, nr2=0;

  int i, isame;

  if ( PyList_Check(py_v) )
  {
    *nr = PyList_Size( py_v );

    if (!(*nr))
    {
      *nc = 0;
      return PYTHON_NONE;
    }

    // get info about the first element in the list
    f1 = PyList_GetItem (py_v,0);
    t1 = _py_object_is_vector (f1, &nr1, &nc1);

    *nc = nr1*nc1;

    // does it have only one element
    if ((*nr) == 1)
      return t1;

    // it has many elements, are they all of the same type as the first one?
    isame = 1;
    for (i=1; i<(*nr) && isame; i++)
    {
      f2 = PyList_GetItem (py_v,i);
      t2 = _py_object_is_vector(f2, &nr2, &nc2);
      if ((t1 != t2) || ((*nc) != nr2*nc2))
        isame = 0;
    }

    if (isame)
      return t1;
  }

  *nc = 0;
  *nr = 0;
  return PYTHON_NONE;
}


Ent *
ent_Py_Initialize (int nargs, Datum args[])
{
//   char *_this_function = "ent_Py_Initialize";

  Ent *e1=0, *e2=0;
  char *fn=0;
  int initsigs=0;

  //
  // first argument: python executable
  //
  if (nargs >= 1)
  {
    e1 = bltin_get_ent(args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      fn = class_char_pointer (e1);
  }

  //
  // second argument: signal handlers?
  //
  if (nargs >= 2)
  {
    e2 = bltin_get_ent(args[1]);
    if (ent_type(e2) == MATRIX_DENSE_REAL)
      initsigs = (int) class_double ( e2 );
  }

  if (fn)
    Py_SetProgramName (fn);

  if ( !Py_IsInitialized () )
    Py_InitializeEx (initsigs);

  ent_Clean (e1);
  ent_Clean (e2);

  Ent *rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (!Py_IsInitialized());
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

Ent *
ent_Py_GetVersion (int nargs, Datum args[])
{
//   char *_this_function = "ent_Py_GetVersion";

  Ent *rent = ent_Create();

  if ( !Py_IsInitialized() )
    rerror (": embedded python has not been initialized!" );
  else
  {
    ent_data(rent) = mds_CreateScalar ( ( (char *) Py_GetVersion () ) );
    ent_type(rent) = MATRIX_DENSE_STRING;
  }

  return rent;
}

Ent *
ent_Py_Finalize (int nargs, Datum args[])
{
  if (Py_IsInitialized())
    Py_Finalize ();

  Ent *rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (Py_IsInitialized());
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

Ent *
ent_PyRun_SimpleString (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;
  MDS *s1;

  int istatus=0, iline=1, n, i;

  //
  if (nargs != 1)
    rerror("one argument required!");

  if (!Py_IsInitialized())
    rerror("initialize python interpreter first!");

  //
  // commands as a string vector
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("first argument has to be string vector!");
  s1 = ent_data ( e1 );

  if ((s1->nrow)!=1 && (s1->ncol)!=1)
    rerror("first argument has to be string vector!");

  n  = (s1->nrow) * (s1->ncol);
  for (i=0; i<n; i++)
  {
    istatus += PyRun_SimpleString ( MdsV0(s1,i) );
    if (istatus)
      break;
    iline = i+1;
  }

  ent_Clean(e1);

  rent = ent_Create();
  if (istatus)
    ent_data(rent) = mdr_CreateScalar ( -iline );
  else
    ent_data(rent) = mdr_CreateScalar ( 0 );
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

Ent *
ent_PyRun_File (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;
  char *f=0;
  FILE *F=0;

  int istatus=0;

  //
  if (nargs != 1)
    rerror("one argument required!");

  if (!Py_IsInitialized())
    rerror("initialize python interpreter first!");

  //
  // commands as a string vector
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("first argument has to be string vector!");
  f = class_char_pointer ( e1 );

  if (!f)
    rerror("first argument has to be string scalar");

  F = fopen(f, "r");
  if (!F)
    rerror("valid filename required as first argument");

  istatus = PyRun_SimpleFileEx(F, f, 1);

  ent_Clean(e1);

  rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar ( istatus );
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

Ent *
ent_PyRun_EvaluateExpression (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  char *v=0, *g=0, *l=0;

  //
  if (nargs < 1 || nargs>3)
    rerror("one, two or three arguments required!");

  if (!Py_IsInitialized())
    rerror("initialize python interpreter first!");

  //
  // commands as a string vector
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("first argument has to be string vector!");
  v = class_char_pointer ( e1 );

  rent = ent_Create();
  if (!v)
  {
    ent_data(rent) = mdr_Create (0,0);
    ent_type(rent) = MATRIX_DENSE_REAL;
    goto return_here;
  }
  MDR *wr=0;
  MDC *wc=0;
  MDS *ws=0;
  int nr, nc, i, j;

  //
  // global dictionary
  //
  if (nargs >= 2)
  {
    e2 = bltin_get_ent(args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
      g = class_char_pointer ( e2 );
  }
  if (!g)
    g = def_module;

  // get the object pointing to the module first
  PyObject *py_g = PyImport_AddModule( g );
  PyObject *py_dict_g = PyModule_GetDict( py_g );
  if (!py_dict_g)
  {
    PyErr_Clear();
    rerror("Terrible Internal Error: Failed to get module dictionary!");
  }

  //
  // local dictionary
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent(args[2]);
    if (ent_type(e3) == MATRIX_DENSE_STRING)
      l = class_char_pointer ( e3 );
  }
  PyObject *py_l=0;
//   PyObject *py_dict_l=0;
  if (l)
  {
    py_l = PyImport_AddModule( l );
//     py_dict_l = PyModule_GetDict( py_l );
  }

  // evaluate string
  PyObject *py_v = PyRun_String ( v, Py_eval_input, py_dict_g, py_l );

  if ( !py_v )
  {
    // default definition
    if (lib_debug)
      fprintf(stdout, "expression '%s' cannot be evaluated in dictionaries "
          "%s and %s!\n", v, g, l);

    // clear error
    PyErr_Clear();

    wr = mdr_Create(0,0);
    ent_data(rent) = wr;
    ent_type(rent) = MATRIX_DENSE_REAL;
    goto return_here;
  }
  Py_INCREF(py_v);


  // figure out if it is a scalar, or a list
  switch ( _py_object_is_scalar (py_v) )
  {
    case PYTHON_INT:

      // integer scalar
      wr = mdi_CreateScalar((int) PyInt_AsLong(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long scalar
      wr = mdi_CreateScalar((int) PyLong_AsLong(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float
      wr = mdr_CreateScalar( PyFloat_AsDouble(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(1,1);
      MdcV0r(wc,0) = PyComplex_RealAsDouble(py_v);
      MdcV0i(wc,0) = PyComplex_ImagAsDouble(py_v);
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ws = mds_Create(1,1);
      MdsV0(ws,0) = cpstr( PyString_AsString(py_v) );
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not scalar
      break;
  }

  switch ( _py_object_is_vector (py_v, &nr, &nc) )
  {
    case PYTHON_INT:
      // integer vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdiV0(wr, i)    = (int) PyInt_AsLong(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdiV0(wr, i)    = (int) PyLong_AsLong(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float vector
      wr = mdr_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdrV0(wr, i)    = PyFloat_AsDouble(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(nr,nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdcV0r(wc,i) = PyComplex_RealAsDouble(py_vi);
        MdcV0i(wc,i) = PyComplex_ImagAsDouble(py_vi);
      }
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ws = mds_Create(nr,nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdsV0(ws,i) = cpstr( PyString_AsString(py_vi) );
      }
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not vector either
      break;
  }

  // is it matrix?
  switch ( _py_object_is_matrix (py_v, &nr, &nc) )
  {
    case PYTHON_INT:
      // integer vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdi0(wr,i,j)    = (int) PyInt_AsLong(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdi0(wr,i,j)    = (int) PyLong_AsLong(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float vector
      wr = mdr_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdr0(wr,i,j)     = PyFloat_AsDouble(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(nr,nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdc0r(wc,i,j) = PyComplex_RealAsDouble(py_vij);
          Mdc0i(wc,i,j) = PyComplex_ImagAsDouble(py_vij);
        }
      }
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ws = mds_Create(nr,nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mds0(ws,i,j) = cpstr( PyString_AsString(py_vij) );        }
      }
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not vector either so we apply the default definition
      wr = mdr_Create(0,0);
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      break;
  }

  // is this necessary?
  Py_DECREF(py_v);

return_here:
  //clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  return rent;
}

Ent *
ent_Py_ManipulateObject (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  MDR *r=0;
  MDC *c=0;
  MDS *s=0;

  char *m=0, *v=0;

  int rval=0, i, j, n;

  PyObject *py_m=0, *py_v=0;
  Py_complex *cpy_c=0;

  //
  if (nargs != 3 && nargs != 2)
    rerror("two or three arguments required!");

  if (!Py_IsInitialized())
    rerror("initialize python interpreter first!");

  //
  // module name: make it tard-proof
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) == MATRIX_DENSE_STRING)
    m = class_char_pointer ( e1 );
  if (!m)
    m = def_module;

  // get the object pointing to the module first
  py_m = PyImport_AddModule( m );
  if (!py_m)
  {
    PyErr_Clear();
    rerror("Module Import Failed!");
  }

  //
  // variable name
  //
  e2 = bltin_get_ent(args[1]);
  if (ent_type(e2) != MATRIX_DENSE_STRING)
    rerror("second argument has to be string scalar!");
  v = class_char_pointer ( e2 );
  if (!v)
    rerror("second argument has to be string scalar!");


  //
  // do we retrieve the value of the argument or set its value?
  //
  if (nargs == 3)
  {
    //
    // variable content
    //
    e3 = bltin_get_ent(args[2]);
    if (ent_type(e3) == UNDEF)
      rerror("third argument cannot be undefined!");

    // delete previous instance of the object if it exists
    if ( PyObject_HasAttrString(py_m, v) )
      PyObject_DelAttrString(py_m, v);

    // instantiante the new object
    switch(ent_type(e3))
    {
      case MATRIX_DENSE_REAL:

        r = ent_data(e3);

        if ((r->nrow)*(r->ncol)>0)
        {
          if (r->nrow == 1 && r->ncol==1)
          {
            // scalar
            if(r->type == RLAB_TYPE_INT32)
              py_v = PyInt_FromLong( (long) MdiV0 (r,0) );
            else
              py_v = PyFloat_FromDouble( MdrV0 (r,0) );

          }
          else if (r->nrow == 1 || r->ncol==1)
          {
            // vector as a python list
            n = (r->nrow) * (r->ncol);
            py_v = PyList_New ( n );
            if(r->type == RLAB_TYPE_INT32)
              for (i = 0; i<n; i++)
                PyList_SetItem(py_v, i, PyInt_FromLong    ( (long) MdiV0 (r,i) ) );
            else
              for (i = 0; i<n; i++)
                PyList_SetItem(py_v, i, PyFloat_FromDouble(        MdrV0 (r,i) ) );
          }
          else
          {
            // matrix as a stacked list of rows
            py_v = PyList_New ( r->nrow );
            for (i=0; i<r->nrow; i++)
            {
              // create a list that contains row of matrix
              PyObject *py_vi = PyList_New ( r->ncol );
              if(r->type == RLAB_TYPE_INT32)
                for (j=0; j<r->ncol; j++)
                  PyList_SetItem(py_vi, j, PyInt_FromLong    ( (long) Mdi0 (r,i,j) ) );
              else
                for (j=0; j<r->ncol; j++)
                  PyList_SetItem(py_vi, j, PyFloat_FromDouble(        Mdr0 (r,i,j) ) );

              // stack a row to the python object
              PyList_SetItem(py_v, i, py_vi);
            }
          }
        }
        break;

      case MATRIX_DENSE_COMPLEX:

        c = ent_data(e3);

        if ((c->nrow)*(c->ncol)>0)
        {
          if (c->nrow == 1 && c->ncol==1)
          {
            // scalar
            cpy_c = (Py_complex *) &MdcV0 (c,0);
            py_v = PyComplex_FromCComplex( *cpy_c );

          }
          else if (c->nrow == 1 || c->ncol==1)
          {
            // vector as a python list
            n = (c->nrow) * (c->ncol);
            py_v = PyList_New ( n );
            for (i=0; i<n; i++)
            {
              cpy_c = (Py_complex *) &MdcV0 (c,i);
              PyList_SetItem(py_v, i, PyComplex_FromCComplex( *cpy_c ) );
            }
          }
          else
          {
            // matrix as a stacked list of rows
            py_v = PyList_New ( c->nrow );
            for (i=0; i<c->nrow; i++)
            {
              // create a list that contains row of matrix
              PyObject *py_vi = PyList_New ( c->ncol );
              for (j=0; j<c->ncol; j++)
              {
                cpy_c = (Py_complex *) &Mdc0 (c,i,j);
                PyList_SetItem(py_vi, j, PyComplex_FromCComplex( *cpy_c ) );
              }

              // stack a row to the python object
              PyList_SetItem(py_v, i, py_vi);
            }
          }
        }
        break;

      case MATRIX_DENSE_STRING:

        s = ent_data(e3);

        if ((s->nrow)*(s->ncol)>0)
        {
          if (s->nrow == 1 && s->ncol==1)
          {
            // scalar
            py_v = PyString_FromString( (MdsV0(s,0)) );

          }
          else if (s->nrow == 1 || s->ncol==1)
          {
            // vector as a python list
            n = (s->nrow) * (s->ncol);
            py_v = PyList_New ( n );
            for (i=0; i<n; i++)
            {
              PyList_SetItem(py_v, i, PyString_FromString( (MdsV0(s,i)) ) );
            }
          }
          else
          {
            // matrix as a stacked list of rows
            py_v = PyList_New ( s->nrow );
            for (i=0; i<s->nrow; i++)
            {
              // create a list that contains row of matrix
              PyObject *py_vi = PyList_New ( s->ncol );
              for (j=0; j<s->ncol; j++)
                PyList_SetItem(py_vi, j, PyString_FromString( (Mds0(s,i, j)) ) );

              // stack a row to the python object
              PyList_SetItem(py_v, i, py_vi);
            }
          }
        }
        break;

      default:

        rerror("The function is not defined for this type of argument!");
        break;

    }

    // at this point we should have PyObject py_v or we are in trouble!
    if (py_v)
    {
      Py_INCREF(py_v);
      PyObject_SetAttrString(py_m, v, py_v);
    }
    else
      rerror ("Terrible internal error: SetAttr failed!");


    // clean up
    ent_Clean(e1);
    ent_Clean(e2);
    ent_Clean(e3);

    // return to caller
    rent = ent_Create();
    ent_data(rent) = mdr_CreateScalar ( rval );
    ent_type(rent) = MATRIX_DENSE_REAL;
    return rent;
  }

  //
  // retrieve content of python variable
  //
  // object exists
  rent = ent_Create();

  MDR *wr=0;
  MDC *wc=0;
  MDS *ws=0;

  int nr, nc;

  // does the object exist? if not go back
  py_v = PyObject_GetAttrString (py_m, v);
  if ( !py_v )
  {
    // default definition
    if (lib_debug)
      fprintf(stdout, "object %s does not exist in module"
          "%s!\n", v, m);

    // clear error
    PyErr_Clear();

    wr = mdr_Create(0,0);
    ent_data(rent) = wr;
    ent_type(rent) = MATRIX_DENSE_REAL;
    goto return_here;
  }
  // Is this necessary?
  Py_INCREF(py_v);

  // figure out if it is a scalar, or a list
  switch ( _py_object_is_scalar (py_v) )
  {
    case PYTHON_INT:

      // integer scalar
      wr = mdi_CreateScalar((int) PyInt_AsLong(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long scalar
      wr = mdi_CreateScalar((int) PyLong_AsLong(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float
      wr = mdr_CreateScalar( PyFloat_AsDouble(py_v) );
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(1,1);
      MdcV0r(wc,0) = PyComplex_RealAsDouble(py_v);
      MdcV0i(wc,0) = PyComplex_ImagAsDouble(py_v);
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ent_data(rent) = mds_CreateScalar( PyString_AsString(py_v) );
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not scalar
      break;
  }

  switch ( _py_object_is_vector (py_v, &nr, &nc) )
  {
    case PYTHON_INT:
      // integer vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdiV0(wr, i)    = (int) PyInt_AsLong(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdiV0(wr, i)    = (int) PyLong_AsLong(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float vector
      wr = mdr_Create (nr, nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdrV0(wr, i)    = PyFloat_AsDouble(py_vi);
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(nr,nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdcV0r(wc,i) = PyComplex_RealAsDouble(py_vi);
        MdcV0i(wc,i) = PyComplex_ImagAsDouble(py_vi);
      }
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ws = mds_Create(nr,nc);
      for (i=0; i<nr*nc; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        MdsV0(ws,i) = cpstr( PyString_AsString(py_vi) );
      }
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not vector either
      break;
  }

  // is it matrix?
  switch ( _py_object_is_matrix (py_v, &nr, &nc) )
  {
    case PYTHON_INT:
      // integer vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdi0(wr,i,j)    = (int) PyInt_AsLong(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_LONG:
      // long vector
      wr = mdi_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdi0(wr,i,j)    = (int) PyLong_AsLong(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_FLOAT:
      // float vector
      wr = mdr_Create (nr, nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdr0(wr,i,j)     = PyFloat_AsDouble(py_vij);
        }
      }
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      goto return_here;

    case PYTHON_COMPLEX:
      // complex
      wc = mdc_Create(nr,nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mdc0r(wc,i,j) = PyComplex_RealAsDouble(py_vij);
          Mdc0i(wc,i,j) = PyComplex_ImagAsDouble(py_vij);
        }
      }
      ent_data(rent) = wc;
      ent_type(rent) = MATRIX_DENSE_COMPLEX;
      goto return_here;

    case PYTHON_STRING:
      // string
      ws = mds_Create(nr,nc);
      for (i=0; i<nr; i++)
      {
        PyObject *py_vi = PyList_GetItem (py_v, i);
        for (j=0; j<nc; j++)
        {
          PyObject *py_vij = PyList_GetItem (py_vi, j);
          Mds0(ws,i,j) = cpstr( PyString_AsString(py_vij) );        }
      }
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
      goto return_here;

    case PYTHON_NONE:
      // it is not vector either so we apply the default definition
      wr = mdr_Create(0,0);
      ent_data(rent) = wr;
      ent_type(rent) = MATRIX_DENSE_REAL;
      break;
  }

  // is this necessary?
  Py_DECREF(py_v);


return_here:

  // clean up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // go back now
  return rent;
}

#endif














