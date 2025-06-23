/* rlab_macros_code.h */

/*  This file is a part of RLaB ("Our"-LaB) and rlabplus
   Copyright (C) 2017 M.Kostrun

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
#ifndef RLAB_MACROS_CODE_L
#define RLAB_MACROS_CODE_L

/**
  * return value:
 *  -1,0   - x is not proper
 *   1     - x is proper, w is not proper
 *   2     - x, w exist are proper
 */
static inline int isent_stat_list(Ent *e, MDR **x, MDR **w)
{
  ListNode *node=0;
  Btree *b=0;
  Ent *ex=0, *ew=0;

  *x = *w = NULL;

  if (!e)
    return -2;

  if (ent_type(e)!=BTREE)
    return -1;

  b = ent_data(e);
  if (!b)
    return -1;

  // check for 'val'
  node = btree_FindNode(b, RLAB_NAME_GEN_VALUE);
  if (!node)
    return 0;
  ex = var_ent (node);
  if (ent_type(ex)!=MATRIX_DENSE_REAL)
    return 0;
  *x = ent_data(ex);
  if (SIZE(*x)<1)
  {
    *x = NULL;
    return 0;
  }

  // check for 'wgt'
  node = btree_FindNode(b, RLAB_NAME_GEN_WEIGHT);
  if (!node)
  {
    return 1;
  }

  ew = var_ent (node);

  if (ent_type(ew)!=MATRIX_DENSE_REAL)
  {
    *w = NULL;
    return 1;
  }

  *w = ent_data(ew);
  if (SIZE(*w)<1)
  {
    *w = NULL;
    return 1;
  }

  // they both have to be same size and double
  if (!EQSIZE(*x,*w) || !MD_TYPE_DOUBLE(*x) || !MD_TYPE_DOUBLE(*w))
  {
    *w = NULL;
    return 1;
  }

  return 2;
}



static inline int md_are_same_size_or_one_is_scalar (MD *x1, MD *x2)
{
  if (SIZE(x1)>0)
    if (SIZE(x2)>0)
    {
      if (EQSIZE(x1,x2))
        return (1);

      if (EQSCAL(x1) || EQSCAL(x2))
        return (1);
    }

  return (0);
}

static inline int md_are_same_size(MD *x1, MD *x2)
{
  if (SIZE(x1)>0)
    if (SIZE(x2)>0)
      if (EQSIZE(x1,x2))
        return (1);

  return (0);
}

static inline int md_is_vector(MD *x1)
{
  if (SIZE(x1)>0)
    if (EQVECT(x1))
      return (1);

  return (0);
}

#define RLABCODE_GETARG_VAR_OR_ENTITY(init,args_idx_p1,ent,node) \
{\
  switch (args[args_idx_p1-1].type) \
  { \
    case VAR: \
      var = (ListNode *) (args[args_idx_p1-1].u.ptr); \
      ent = var_ent (var); \
      if (ent->refc > 1) \
      { \
        ent_DecRef (ent); \
        ent = ent_Duplicate (ent); \
        listNode_AttachEnt (node, ent); \
      } \
      break; \
    \
    case ENTITY: \
      ent = (Ent *) args[args_idx_p1-1].u.ptr; \
      if (ent->refc > 1) \
      { \
        ent_DecRef (ent); \
        ent = ent_Duplicate (ent); \
      } \
      break; \
    \
    default: \
      rerror(init ": " #args_idx_p1 " " RLAB_ERROR_ENTITY_OR_VAR "\n"); \
  }\
}


#define RLABCODE_INSTALL_ENTITY_IN_BTREE(b,ent,etyp,edat,enode) \
{\
  ent = ent_Create(); \
  ent_type(ent) = etyp; \
  ent_data(ent) = edat; \
  install  (b, enode, ent); \
}

#define RLABCODE_PROCESS_ARG1_S(init,arg_idx,ent,s,minlen) \
{\
  ent = bltin_get_ent( args[(arg_idx)] ); \
  if (ent_type(ent) != MATRIX_DENSE_STRING) \
    rerror (init ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n"); \
  s = class_char_pointer(ent); \
  if (isvalidstring(s) < (minlen)) \
    rerror (init ": " RLAB_ERROR_ARG1_MDS_SCALAR " Its length must be greater then " #minlen "\n"); \
}

#define RLABCODE_PROCESS_ARG1_S_NOERRORS(init,arg_idx,ent,s,minlen) \
{\
  s = NULL; \
  ent = bltin_get_ent( args[(arg_idx)] ); \
  if (ent_type(ent) == MATRIX_DENSE_STRING) \
  { \
    s = class_char_pointer(ent); \
    if (isvalidstring(s) < (minlen)) \
      s = NULL; \
  } \
}


#define RLABCODE_PROCESS_ARG_S(init,arg_idx,ent,s,minlen,errdescr) \
{\
  ent = bltin_get_ent( args[(arg_idx)] ); \
  if (ent_type(ent) != MATRIX_DENSE_STRING) \
    rerror (init ": " errdescr "\n"); \
  s = class_char_pointer(ent); \
  if (isvalidstring(s) < (minlen)) \
    rerror (init ": " errdescr " Its length must be greater then " #minlen "\n"); \
}

#define RLABCODE_PROCESS_ARG_MD_IGNORE_ERROR(init,arg_idx,ent,etype,d,entfunc,dfunc,op,dvalerr,derr) \
{\
  ent = bltin_get_ent( args[(arg_idx)] ); \
  if (ent_type(ent) == etype) \
  {\
    d = entfunc(ent); \
    if (dfunc(d) op (dvalerr)) \
      d = derr; \
  }\
}


#define RLABCODE_PROCESS_ARG_F_S(init,arg_idx,ent,etype,d,funcd,op,dval,errdescr) \
{\
  ent = bltin_get_ent( args[(arg_idx)] ); \
  if (ent_type(ent) != etype) \
    rerror (init ": " errdescr "\n"); \
  d = ent_data( ent ); \
  if (funcd(d) op dval) \
  { \
    rerror (init ": " errdescr "\n" ); \
  } \
}


#define RLABCODE_PROCESS_BTREE_ENTRY_S(ent,node,name,s,minlen,dfltval) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    s = (char *) class_char_pointer( var_ent(node) ); \
    if (isvalidstring(s) < (minlen)) \
      s = dfltval; \
  } \
}

#define RLABCODE_PROCESS_BTREE_ENTRY_MD(ent,node,name,s,stype,f,fval,dfltval) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    s = (stype *) ent_data( var_ent(node) ); \
    if (f(s) < (fval)) \
      s = dfltval; \
  } \
}

#define RLABCODE_PROCESS_BTREE_ENTRY_F_S(ent,node,name,s,funcent,sval,dfltval) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    s = funcent( var_ent(node) ); \
    if (s < (sval)) \
      s = dfltval; \
  } \
}

#define RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(ent,node,name,d) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    d  = class_double ( var_ent(node) ); \
  } \
}
#define RLABCODE_PROCESS_BTREE_ENTRY_D(ent,node,name,d,funcent,op,dval) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    double _dummy= funcent ( var_ent(node) ); \
    if (_dummy op (dval)) \
      d = _dummy; \
  } \
}

#define RLABCODE_PROCESS_BTREE_ENTRY_D_INRANGE(ent,node,name,d,funcent,dmin,dmax) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    double _dummy= funcent( var_ent(node) ); \
    if ((_dummy >= (dmin)) && (_dummy <= (dmax))) \
      d = _dummy; \
  } \
}

#define RLABCODE_PROCESS_BTREE_ENTRY_BOOL(ent,node,name,d) \
{\
  node = btree_FindNode( ent_data(ent) , name ); \
  if (node) \
  { \
    int _dummy = class_int( var_ent(node) ); \
    if (_dummy) \
      d = 1; \
    else \
      d = 0; \
  } \
}


#define RLABCODE_CHECK_IF_FILENAME_EXISTS(rf,fname,close_after_rw) \
{\
  rf = rfile_find(fname); \
  if (!rf) \
  { \
    if (*fname == '|') \
    { \
      rf = get_rfile_ds (fname, RFILE_FILE, "r", 0); \
    } \
    else if (    !strncmp(fname, "h5://", 5) \
             ||  !strncmp(fname, "hdf5://", 7) \
             ||  (strstr (fname, ".h5")) \
       ) \
    { \
      rf = get_rfile_ds (fname, RFILE_H5, "r" ); \
    } \
    else if (!strncmp(fname, "tcp://", 6)) \
    { \
      rf = get_rfile_ds (fname, RFILE_SOCKET, 1); \
    } \
    else if (    (!strncmp(fname, "http://", 7)) \
                ||  (!strncmp(fname, "https://", 8)) \
                ||  (!strncmp(fname, "ftp://", 6)) \
         ) \
    { \
      rf = get_rfile_ds(fname, RFILE_CURL); \
    } \
    else \
    { \
      rf = get_rfile_ds (fname, RFILE_FILE, "r", 0); \
    } \
    close_after_rw = 1; \
  } \
}


#endif
