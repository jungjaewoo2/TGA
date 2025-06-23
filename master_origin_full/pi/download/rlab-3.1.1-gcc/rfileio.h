/* fileio.h */

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

#ifndef RLAB_RFILEIO_H
#define RLAB_RFILEIO_H

#include "rlab.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"
#include "btree.h"

#include "stdio.h"

extern void chomp_string(char *);
extern void chomp2_string(char *);
extern FILE * RLAB_STDERR_DS;
extern FILE * get_file_ds      (char *name, char *mode, int buffsize);
extern int    get_int_file_ds    (char *name);
extern char * get_eol_file_name  (char *name);
extern char * get_file_ds_name (FILE *fn);
extern char * get_eol_file_ds_name (FILE * fn);
extern char * get_csp_file_ds_name (FILE * fn);
extern MDS  * get_fmt_file_ds_name (FILE * fn);
extern int    get_file_ds_intfd(char *name, int ispeed, int idatb, int iparity,
                                int istop,  int iflow, int iraw, int ibuff );
#ifdef HAVE_LIBCURL
# include <curl/curl.h>
extern CURL * get_file_ds_curl(char *name);
#endif

extern int close_file_ds (char *name);

extern void class_io_init (void);

extern Ent *ReadM (int nargs, Datum args[]);
extern Ent *WriteM (int nargs, Datum args[]);

extern Ent *ReadB (int nargs, Datum args[]);
extern Ent *WriteB (int nargs, Datum args[]);

extern Ent *ReadASCII (int nargs, Datum args[]);
extern Ent *WriteASCII (int nargs, Datum args[]);

extern void *matrix_ReadB (FILE * fn, int T, int P, char **name, int *rtype);
extern Btree *btree_ReadB (FILE * fn, char **name);

extern void mdr_WriteB (MDR * m, FILE * fn, char *name);
extern void mdc_WriteB (MDC * m, FILE * fn, char *name);
extern void mds_WriteB (MDS * m, FILE * fn, char *name);
extern void btree_WriteB (Btree * btree, FILE * fn, char *name);

extern Ent *Fread (int nargs, Datum args[]);
extern Ent *Fwrite (int nargs, Datum args[]);
extern Ent *Fseek (int nargs, Datum args[]);

extern Ent *ListOpenFiles (int nargs, Datum args[]);
extern Ent *ReadS (int nargs, Datum args[]);

#ifdef HAVE_HDF5_SO
// HDF5 specific functions
extern Ent *ent_hdf5_isfile (int nargs, Datum args[]);
extern Ent *ent_h5ls (int nargs, Datum args[]);
extern Ent *ent_h5mv (int nargs, Datum args[]);
extern Ent *ent_h5cp (int nargs, Datum args[]);
extern Ent *ent_h5ln (int nargs, Datum args[]);
#endif

#endif /* RLAB_RFILEIO__H */
