// Copyright (C) 2003-2006 Marijan Kostrun
//   part of rlabplus for linux project on rlabplus.sourceforge.net
//
// gnuplot through pipe for rlabplus
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

#include "rlab.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"
#include "mdr_mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "complex.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <stdarg.h>

#define  GNUPLOT_MAX_NUM_WINDOWS 32
#define  GNUPLOT_MAX_NAME_SIZE   512
#define  GNUPLOT_SET_TERM        "set term x11"
#define  GNUPLOT_MAX_CMD_SIZE    2048

/** Maximal size of a name in the PATH */
#define PATH_MAXNAMESZ       4096

/** Define P_tmpdir if not defined (this is normally a POSIX symbol) */
#ifndef P_tmpdir
#define P_tmpdir "."
#endif


/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.h
  @author   N. Devillard
  @date     Sep 1998
  @version  $Revision: 1.11 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.
  
 */
/*--------------------------------------------------------------------------*/

/*
    $Id: gnuplot_i.h,v 1.11 2003/01/27 08:58:04 ndevilla Exp $
    $Author: ndevilla $
    $Date: 2003/01/27 08:58:04 $
    $Revision: 1.11 $
 */

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_
#endif

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

/** Maximal size of a temporary file name/pipe */


/*---------------------------------------------------------------------------
                                New Types
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.c
  @author   N. Devillard
  @date Sep 1998
  @version  $Revision: 2.10 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

 */


/*-------------------------------------------------------------------------*/
/**
  @brief    Find out where a command lives in your PATH.
  @param    pname Name of the program to look for.
  @return   pointer to statically allocated character string.

  This is the C equivalent to the 'which' command in Unix. It parses
  out your PATH environment variable to find out where a command
  lives. The returned character string is statically allocated within
  this function, i.e. there is no need to free it. Beware that the
  contents of this string will change from one call to the next,
  though (as all static variables in a function).

  The input character string must be the name of a command without
  prefixing path of any kind, i.e. only the command name. The returned
  string is the path in which a command matching the same name was
  found.

  Examples (assuming there is a prog named 'hello' in the cwd):

  @verbatim
  gnuplot_get_program_path("hello") returns "."
  gnuplot_get_program_path("ls") returns "/bin"
  gnuplot_get_program_path("csh") returns "/usr/bin"
  gnuplot_get_program_path("/bin/ls") returns NULL
  @endverbatim
  
 */
/*-------------------------------------------------------------------------*/

char * gnuplot_get_program_path(char * pname)
{
    int         i, j, lg;
    char    *   path;
    static char buf[PATH_MAXNAMESZ];

    /* Trivial case: try in CWD */
    sprintf(buf, "./%s", pname) ;
    if (access(buf, X_OK)==0) {
        sprintf(buf, ".");
        return buf ;
    }
    /* Try out in all paths given in the PATH variable */
    buf[0] = 0;
    path = getenv("PATH") ;
    if (path!=NULL) {
        for (i=0; path[i]; ) {
            for (j=i ; (path[j]) && (path[j]!=':') ; j++);
            lg = j - i;
            strncpy(buf, path + i, lg);
            if (lg == 0) buf[lg++] = '.';
            buf[lg++] = '/';
            strcpy(buf + lg, pname);
            if (access(buf, X_OK) == 0) {
                /* Found it! */
                break ;
            }
            buf[0] = 0;
            i = j;
            if (path[i] == ':') i++ ;
        }
    } else {
        fprintf(stderr, "PATH variable not set\n");
    }
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0) return NULL ;
    /* Otherwise truncate the command name to yield path only */
    lg = strlen(buf) - 1 ;
    while (buf[lg]!='/') {
        buf[lg]=0 ;
        lg -- ;
    }
    buf[lg] = 0;
    return buf ;
}


static FILE *
gnuplot_pipe_init( char * errname )
{
  FILE *  handle ;
  char pipcmd[128] = {'\0'};

  if (getenv("DISPLAY") == NULL)
  {
    fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
  }
  if (gnuplot_get_program_path("gnuplot")==NULL)
  {
    fprintf(stderr, "cannot find gnuplot in your PATH");
    return NULL ;
  }

  if (errname)
  {
    if (strlen(errname)>0)
      sprintf(pipcmd, "gnuplot 1>%s 2>&1", errname);
    else
      sprintf(pipcmd, "gnuplot 1>/dev/null 2>&1");
  }
  else
    sprintf(pipcmd, "gnuplot 1>/dev/null 2>&1");

  handle = popen(pipcmd, "w") ;

  if (handle == NULL)
  { fprintf(stderr, "error starting gnuplot\n"); }

  return handle;
}

static FILE *
gnuplot_file_init( char * filename )
{
  FILE *  handle ;
  char pipcmd[128] = {'\0'};

  if (filename)
  {
    if (strlen(filename)>0)
      handle = fopen(filename, "w") ;
  }

  if (handle == NULL)
  { fprintf(stderr, "error writing to file\n"); }

  return handle;
}


static void
gnuplot_cmd(FILE * handle, char *  cmd, ...)
{
  va_list ap ;
  char    local_cmd[GNUPLOT_MAX_CMD_SIZE];

  va_start(ap, cmd);
  vsprintf(local_cmd, cmd, ap);
  va_end(ap);

  strcat(local_cmd, "\n");

  fputs(local_cmd, handle) ;
  fflush(handle) ;
  return ;
}

static int GNUPLOT_NUM_OPEN_WINDOWS = 0;
static int GNUPLOT_DEFAULT_WINDOW      = 0;
static FILE * gnuplot_h     [ GNUPLOT_MAX_NUM_WINDOWS ] = {0};
static int    gnuplot_ispipe[ GNUPLOT_MAX_NUM_WINDOWS ] = {0};
static char * gnuplot_fnames[ GNUPLOT_MAX_NUM_WINDOWS ] = {0};

//
// read or set the default gnuplot window
//  I = gnuwin();   // read the default window
//  gnuwin(I);      // set I as a default window
//
Ent * ent_gnuplot_win(int nargs, Datum args[])
{
  Ent *rent=0, *e1=0;
  MDR *w;
  int i;

  if (nargs > 1)
    rerror("gnuwin: no or one argument required");

  if (nargs == 0)
  {
    if (GNUPLOT_NUM_OPEN_WINDOWS == 0)
    { w = mdr_Create(0,0); }
    else
    { w = mdr_CreateScalar(GNUPLOT_DEFAULT_WINDOW); }
  }
  else if (nargs == 1)
  {
    if (GNUPLOT_NUM_OPEN_WINDOWS == 0)
    { w = mdr_Create(0,0); }
    else
    {
      e1 = bltin_get_ent (args[0]);
      i  = (int) class_double (e1);
      i  = i > 0 ? i : 1;
      i  = i < GNUPLOT_MAX_NUM_WINDOWS ? i : GNUPLOT_MAX_NUM_WINDOWS;
      if (gnuplot_h[i-1] != NULL)
        GNUPLOT_DEFAULT_WINDOW = i;
       w = mdr_CreateScalar(GNUPLOT_DEFAULT_WINDOW);
    }
  }

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);

  rent = ent_Create();
  ent_data(rent) = w;
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

Ent * ent_gnuplot_allwin(int nargs, Datum args[])
{
  Ent   *rent=0, *R=0;
  MDR   *wd, *wa;
  MDS   *ws;
  int    i, j=0;
  Btree *bw;

  //
  // check arguments:
  //    gnuwin()    - returns info about available gnuplot devices
  //    gnuwin(...) - sets the devices
  //
  if (nargs == 0)
  {
    if (GNUPLOT_NUM_OPEN_WINDOWS == 0)
    {
      wa = wd = mdr_Create(0,0);
      ws = mds_Create(0,0);
    }
    else
    {
      wd = mdr_CreateScalar(GNUPLOT_DEFAULT_WINDOW);
      wa = mdr_Create(1,GNUPLOT_NUM_OPEN_WINDOWS);
      ws = mds_Create(1,GNUPLOT_NUM_OPEN_WINDOWS);
      i = 0;
      for (j = 0; j < GNUPLOT_MAX_NUM_WINDOWS; j++)
      {
        if (gnuplot_h[j] != NULL)
        {
          MdrV0(wa,i) = j + 1;
          MdsV0(ws,i) = cpstr(gnuplot_fnames[i]);
          ++i;
        }
      }
    }

    bw = btree_Create ();

    R  = ent_Create ();
    ent_data (R) = wa;
    ent_type (R) = MATRIX_DENSE_REAL;
    install (bw, "available", R);
    ent_Clean (R);

    R  = ent_Create ();
    ent_data (R) = wd;
    ent_type (R) = MATRIX_DENSE_REAL;
    install (bw, "default", R);
    ent_Clean (R);

    R  = ent_Create ();
    ent_data (R) = ws;
    ent_type (R) = MATRIX_DENSE_STRING;
    install (bw, "device", R);
    ent_Clean (R);

    rent = ent_Create ();
    ent_data (rent) = btree_Copy(bw);
    ent_type (rent) = BTREE;
    return rent;
  }
  else if (nargs>=1 && nargs <=3)
  {
  }
}

//
// write to file:
//  I = gnustart(filename)
// write to pipe:
//  I = gnustart()
//  I = gnustart(,term)
//  I = gnustart(,term, stderr)
//
Ent * ent_gnuplot_init(int nargs, Datum args[])
{
  Ent *rent=0, *e1=0, *e2=0, *e3=0;
  int i;
  char * fn = 0, * se = 0;

  //
  if (nargs > 3)
    rerror("gnustart: requires no, one or two arguments");

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    fn = class_char_pointer (e1);
    if (strlen(fn)==0)
      fn = 0;
  }
  // third argument is the stderr for messages from gnuplot
  else if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    se = class_char_pointer (e3);
  }

  //
  // gnuplot
  //
  if (GNUPLOT_NUM_OPEN_WINDOWS < GNUPLOT_MAX_NUM_WINDOWS)
  {
    if (GNUPLOT_NUM_OPEN_WINDOWS == GNUPLOT_DEFAULT_WINDOW)
    {
      if (fn)
      {
        gnuplot_h     [ GNUPLOT_DEFAULT_WINDOW ] = gnuplot_file_init (fn);
        gnuplot_fnames[ GNUPLOT_DEFAULT_WINDOW ] = cpstr(fn);
        gnuplot_ispipe[ GNUPLOT_DEFAULT_WINDOW ] = 0;
      }
      else
      {
        gnuplot_h     [ GNUPLOT_DEFAULT_WINDOW ] = gnuplot_pipe_init (se);
        gnuplot_fnames[ GNUPLOT_DEFAULT_WINDOW ] = cpstr( "|gnuplot" );
        gnuplot_ispipe[ GNUPLOT_DEFAULT_WINDOW ] = 1;
      }
      GNUPLOT_DEFAULT_WINDOW ++;
      GNUPLOT_NUM_OPEN_WINDOWS ++;
    }
    else
    {
      for (i=0; i < GNUPLOT_DEFAULT_WINDOW; i++)
      {
        if (gnuplot_h[i] == NULL)
        {
          if (fn)
          {
            gnuplot_h     [ i ] = gnuplot_file_init (fn);
            gnuplot_fnames[ i ] = cpstr( fn );
            gnuplot_ispipe[ i ] = 0;
          }
          else
          {
            gnuplot_h     [ i ] = gnuplot_pipe_init (se);
            gnuplot_fnames[ i ] = cpstr( "|gnuplot" );
            gnuplot_ispipe[ i ] = 1;
          }

          GNUPLOT_NUM_OPEN_WINDOWS ++;
          GNUPLOT_DEFAULT_WINDOW = i+1;
          break;
        }
      }
    }
    if (!fn)
    {
      gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ], "set mouse" );
      gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ],
                 GNUPLOT_SET_TERM " title 'Gnuplot %i: rlabplus'",
                 GNUPLOT_DEFAULT_WINDOW
               );
    }

  }
  else
    rerror("gnustart: maximum number of gnuplot windows reached!");

  //
  // send commands to the pipe
  //
  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      char * s1 = class_char_pointer (e2);
      if (!s1)
        if (strlen(s1)>0)
          gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ],
                       "%s", s1 );
    }
  }
  else if (nargs == 0)
  {
    gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ],
                 GNUPLOT_SET_TERM " title 'Gnuplot %i: rlabplus'",
                 GNUPLOT_DEFAULT_WINDOW
               );
  }

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);
  if (e2)
    if (ent_type(e2)!=UNDEF)
      ent_Clean (e2);
  if (e3)
    if (ent_type(e3)!=UNDEF)
      ent_Clean (e3);


  rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (GNUPLOT_DEFAULT_WINDOW);
  ent_type(rent) = MATRIX_DENSE_REAL;

  return rent;
}

//
// gnuclose (I)
//
Ent * ent_gnuplot_close (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;
  MDR *w;

  int i, j;

  if (nargs != 1)
    rerror("gnuclose: requires single argument!");

  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror("gnuclose: requires single argument!");
  i = (int) class_double ( e1 );

  if (i>0 && i<GNUPLOT_MAX_NUM_WINDOWS)
  {
    if (gnuplot_h[ i-1 ])
    {
      if (gnuplot_ispipe[i-1])
      {
        gnuplot_cmd(gnuplot_h[ i-1 ], "quit;\n"); // tell Gnuplot to turn itself off
        pclose(gnuplot_h[ i-1 ]);
      }
      else
        fclose(gnuplot_h[ i-1 ]);

      gnuplot_h     [ i-1 ] = NULL;
      gnuplot_ispipe[ i-1 ] = 0;
      GC_free(gnuplot_fnames[ i-1 ]);
      gnuplot_fnames[ i-1 ] = 0;
      if (GNUPLOT_DEFAULT_WINDOW == i)
      {
        for (j = GNUPLOT_MAX_NUM_WINDOWS; j > 0; j--)
        {
          if (gnuplot_h[ j-1 ])
          {
            GNUPLOT_DEFAULT_WINDOW = j;
            break;
          }
        }
      }
      GNUPLOT_NUM_OPEN_WINDOWS --;
    }

  }

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);

  rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (1.0);
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

//
// gnucmd(cmd)
//
Ent * ent_gnuplot_cmd (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;
  MDR *w;
  MDS *s1;

  int i, j, n;

  if (nargs != 1)
    rerror("gnucmd: single argument required");

  i = GNUPLOT_DEFAULT_WINDOW;

  if (!i || !gnuplot_h[ i-1 ])
    rerror("gnucmd: something is terribly wrong, accessing non-window!");

  //
  // s, commands for gnuplot
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("gnucmd: string vector for first argument required");
  s1 = ent_data ( e1 );

  if (!s1)
    rerror("gnucmd: string vector for first argument required");

  if (MNR(s1)!=1 && MNC(s1)!=1)
    rerror("gnucmd: string vector for first argument required");

  n = MNR (s1) * MNC (s1);
  for (j=0; j<n; j++)
  {
    if (strlen(MdsV0(s1,j))>0)
    {
      gnuplot_cmd( gnuplot_h[ i-1 ], MdsV0(s1,j) );
    }
  }

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);

  rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (1.0);
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

//
// gnuprint (x, y, z), or
// gnuprint (m)
//
Ent * ent_gnuplot_print (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  MDR *w ;

  int i, j, k, n;

  if (nargs != 1 && nargs != 3)
    rerror("gnuprint: one or three arguments required");

  //
  // i, index of the gnuplot window
  //
  i = GNUPLOT_DEFAULT_WINDOW;

  if (!i || !gnuplot_h[ i-1 ])
    rerror("gnucmd: something is terribly wrong, accessing non-window!");

  if (nargs == 1)
  {
    e1 = bltin_get_ent(args[0]);
    if (ent_type(e1) != MATRIX_DENSE_REAL)
      rerror("gnuprint: second argument has to be a matrix!");
    MDR * x = ent_data ( e1 );
    int nr = MNR (x);
    int nc = MNC (x);
    for (j = 0; j < nr; j++)
    {
      for (k = 0; k < nc ; k++)
      {
        if (x->isint)
          fprintf(gnuplot_h[ i-1 ], "  %i", Mdi0(x,j,k)) ;
        else
          fprintf(gnuplot_h[ i-1 ], "  %g", Mdr0(x,j,k)) ;
      }
      fprintf(gnuplot_h[ i-1 ], "\n") ;
    }
    fprintf(gnuplot_h[ i-1 ], "e\n") ;
    fflush (gnuplot_h[ i-1 ]) ;
  }
  else if (nargs == 3)
  {
    e1 = bltin_get_ent(args[0]);
    if (ent_type(e1) != MATRIX_DENSE_REAL)
      rerror("gnuprint: second argument has to be a vector!");
    MDR *x = ent_data ( e1 );
    if (MNR(x)!=1 && MNC(x)!=1)
      rerror("gnuprint: second argument has to be a vector!");
    int nr = MNR(x) * MNC(x);

    e2 = bltin_get_ent(args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror("gnuprint: third argument has to be a vector!");
    MDR *y = ent_data ( e2 );
    if (MNR(y)!=1 && MNC(y)!=1)
      rerror("gnuprint: third argument has to be a vector!");
    int nc = MNR(y) * MNC(y);
    
    e3 = bltin_get_ent(args[2]);
    if (ent_type(e3) != MATRIX_DENSE_REAL)
      rerror("gnuprint: fourth argument has to be a matrix!");
    MDR *z = ent_data ( e3 );
    if (MNR(z)!=nr || MNC(z)!=nc)
      rerror("gnuprint: the size of 'z' does not match that of 'x' and 'y'!");

    for (j = 0; j < nr; j++)
    {
      for (k = 0; k < nc; k++)
      {
        // x
        if (x->isint)
          fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(x,j)) ;
        else
          fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(x,j)) ;
        // y
        if (y->isint)
          fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(y,k)) ;
        else
          fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(y,k)) ;
        // z
        if (z->isint)
          fprintf(gnuplot_h[ i-1 ], "  %i\n", Mdi0(z,j,k)) ;
        else
          fprintf(gnuplot_h[ i-1 ], "  %g\n", Mdr0(z,j,k)) ;
      }
      fprintf(gnuplot_h[ i-1 ], "\n") ;
    }
    fprintf(gnuplot_h[ i-1 ], "e\n") ;
    fflush (gnuplot_h[ i-1 ]) ;

    if (e2)
      if (ent_type(e2)!=UNDEF)
        ent_Clean (e2);
    if (e3)
      if (ent_type(e3)!=UNDEF)
        ent_Clean (e3);
  }

  if (e1)
    if (ent_type(e1)!=UNDEF)
      ent_Clean (e1);

  rent = ent_Create();
  ent_data(rent) = mdr_CreateScalar (1.0);
  ent_type(rent) = MATRIX_DENSE_REAL;
  return rent;
}

