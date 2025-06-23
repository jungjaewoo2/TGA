// Copyright (C) 2003-2008 Marijan Kostrun
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
#include "mdrf1.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "complex.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "symbol.h"
#include "list.h"
#include "btree.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>

#define  GNUPLOT_MAX_NUM_WINDOWS 32
#define  GNUPLOT_MAX_NAME_SIZE   512
#define  GNUPLOT_SET_TERM        "set term x11 enhanced title 'Gnuplot %i: rlabplus'"
#define  GNUPLOT_DEFAULT_PFLAGS  "-noraise"
#define  GNUPLOT_MAX_CMD_SIZE    32768

/** Maximal size of a name in the PATH */
#define PATH_MAXNAMESZ       4096

/** Define P_tmpdir if not defined (this is normally a POSIX symbol) */
#ifndef P_tmpdir
#define P_tmpdir "."
#endif

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_
#endif

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

static int rlab_gnuplot_idebug=0;


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
      sprintf(pipcmd, "gnuplot " GNUPLOT_DEFAULT_PFLAGS " 1>%s 2>&1", errname);
    else
      sprintf(pipcmd, "gnuplot " GNUPLOT_DEFAULT_PFLAGS " 1>/dev/null 2>&1");
  }
  else
    sprintf(pipcmd, "gnuplot " GNUPLOT_DEFAULT_PFLAGS " 1>/dev/null 2>&1");

  handle = popen(pipcmd, "w") ;

  if (handle == NULL)
  { fprintf(stderr, "error starting gnuplot\n"); }

  return handle;
}

static FILE *
gnuplot_file_init( char * filename )
{
  FILE *  handle=0;

  if (filename)
  {
    if (strlen(filename)>0)
      handle = fopen(filename, "w") ;
  }

  if (handle == NULL)
  { fprintf(stderr, "error writing to file\n"); }

  return handle;
}

// static void
// gnuplot_cmd(FILE * handle, char *  cmd, ...)
// {
//   va_list ap ;
//   char    local_cmd[GNUPLOT_MAX_CMD_SIZE];
//
//   va_start(ap, cmd);
//   vsprintf(local_cmd, cmd, ap);
//   va_end(ap);
//
//   strcat(local_cmd, "\n");
//
//   fputs(local_cmd, handle) ;
//   fflush(handle) ;
//   return ;
// }

static int GNUPLOT_NUM_OPEN_WINDOWS = 0;
static int GNUPLOT_DEFAULT_WINDOW   = 0;
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
  Ent *e1=0;
  MDR *w=0;
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
      {
        GNUPLOT_DEFAULT_WINDOW = i;
        fputs ("reset\n", gnuplot_h[ i-1 ]);
        fflush(gnuplot_h[ i-1 ]);
      }
      w = mdr_CreateScalar(GNUPLOT_DEFAULT_WINDOW);
    }
  }

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

Ent * ent_gnuplot_allwin(int nargs, Datum args[])
{
  MDR   *wd=0, *wa=0;
  MDS   *ws=0;
  int    i, j=0;
  Btree *bw;

  //
  // check arguments:
  //    gnuwins()    - returns info about available gnuplot devices
  //
  if (GNUPLOT_NUM_OPEN_WINDOWS > 0)
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
        MdsV0(ws,i) = cpstr(gnuplot_fnames[j]);
        ++i;
      }
    }
  }

  //
  //
  //
  bw = btree_Create ();
  install (bw, RLAB_NAME_GNUPLOT_GNUWINS_WINDOWS, ent_Assign_Rlab_MDR(wa));
  install (bw, RLAB_NAME_GNUPLOT_GNUWINS_ACTIVE, ent_Assign_Rlab_MDR(wd));
  install (bw, RLAB_NAME_GNUPLOT_GNUWINS_DEVICE, ent_Assign_Rlab_MDS(ws));
  return ent_Assign_Rlab_BTREE(bw);
}

//
// write to file:
//  I = gnustart(filename)
// write to pipe:
//  I = gnustart(fn)
//  I = gnustart();
//  I = gnustart(  ,term)
//  I = gnustart(  ,term, stderr)
//
Ent * ent_gnuplot_init(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int i;
  char *fn=0, *se=0, *s1=0;

  char gnuwin_label[128];

  //
  if (nargs > 3)
    rerror("gnustart: requires none, one, two or three arguments");

  //
  // first argument: name of the file to which all commands in a particular
  // session will be written to
  //
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
    {
      fn = class_char_pointer (e1);
      if (!fn)
        if (strlen(fn)==0)
          fn = 0;
    }
  }

  //
  // second argument: commands to be sent to a gnuplot session that is
  // about to be open
  //
  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      s1 = class_char_pointer (e2);
      if (!s1)
        if (strlen(s1)==0)
          s1 = 0;
    }
  }

  // third argument is the stderr for messages from gnuplot
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) == MATRIX_DENSE_STRING)
    {
      se = class_char_pointer (e3);
      if (se)
        if (strlen(se)==0)
          se = 0;
    }
  }
  if (!fn && !s1 && !se)
    rerror("gnustart: inconsistent input - cannot set all arguments !");

  //
  // start a gnuplot session
  //
  if (GNUPLOT_NUM_OPEN_WINDOWS > GNUPLOT_MAX_NUM_WINDOWS)
    rerror("gnustart: maximum number of gnuplot windows reached!");

  if (GNUPLOT_NUM_OPEN_WINDOWS == GNUPLOT_DEFAULT_WINDOW)
  {
    //
    // all gnuplot sessions 1:GNUPLOT_NUM_OPEN_WINDOWS are being used: open a new
    // gnuplot session and make it default
    //
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
    //
    // in the range 1:GNUPLOT_NUM_OPEN_WINDOWS some sessions were previously
    // closed. go over all indices and check if their session is still on.
    // if not, use that index to open new session.
    //
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


  //
  // send commands to the pipe with two possibilities:
  // (1) user did not provide any parameters (e2) for the terminal - assume X11 and
  // initiate it appropriately,
  // (2) user provided arguments for the terimanal as e2, send them to the default
  // gnuplot session
  //

  // user did not provide any parameters: this is an gnuplot session with x11 output
  if (nargs==0 || (!s1))
  {
//     gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ], "set mouse" );
//     gnuplot_cmd( gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ],
//                  GNUPLOT_SET_TERM,
//                  GNUPLOT_DEFAULT_WINDOW
//     );
    fputs("set mouse", gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    fputs("\n", gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    // use GNUPLOT_DEFAULT_WINDOW to index the open gnuwin plots in the title bar
    sprintf(gnuwin_label, GNUPLOT_SET_TERM, GNUPLOT_DEFAULT_WINDOW);
    fputs(gnuwin_label, gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    //
    fputs("\n", gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    fflush(gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
  }
  else if (nargs>=2 && s1)
  {
    // allow user to index open gnuplot sessions through '%i' in their specification
    sprintf(gnuwin_label, s1, GNUPLOT_DEFAULT_WINDOW);
    fputs(gnuwin_label, gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
//   fputs(s1, gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    fputs("\n", gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
    fflush(gnuplot_h[ GNUPLOT_DEFAULT_WINDOW - 1 ]);
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Int(GNUPLOT_DEFAULT_WINDOW);
}

//
// gnuclose (I)
//
Ent * ent_gnuplot_close (int nargs, Datum args[])
{
  Ent *e1=0;

  MDR *x=0;

  int i, n, j, k;

  if (nargs != 1)
    rerror("gnuclose: requires single argument!");

  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
  {
    ent_Clean (e1);
    return ent_Create_Rlab_Int(2);
  }

  x = ent_data(e1);
  n = x->nrow * x->ncol;

  for (k=0; k<n; k++)
  {
    switch (x->type)
    {
      case RLAB_TYPE_INT32:
        i = MdiV0(x,k);
        break;

      case RLAB_TYPE_DOUBLE:
        i = (int) MdrV0(x,k);
        break;

      default:
        i = -1;
        break;
    }

    if (i<0 || i>=GNUPLOT_MAX_NUM_WINDOWS)
      continue;

    if (!gnuplot_h[ i-1 ])
      continue;

    if (gnuplot_ispipe[i-1])
    {
      // tell Gnuplot to off itself
      fputs("quit;\n", gnuplot_h[ i-1 ]);
      fflush(gnuplot_h[ i-1 ]);
      pclose(gnuplot_h[ i-1 ]);
    }
    else
    {
      fflush(gnuplot_h[ i-1 ]);
      fclose(gnuplot_h[ i-1 ]);
    }

    gnuplot_h     [ i-1 ] = NULL;
    gnuplot_ispipe[ i-1 ] = 0;
    GC_free (gnuplot_fnames[ i-1 ]);
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

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

//
// gnucmd(cmd)
//
#undef  THIS_SOLVER
#define THIS_SOLVER "gnucmd"
Ent * ent_gnuplot_cmd (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDS *s1=0;
  int i, j, n;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  char * src=0;

  if (nargs != 1 && nargs !=2)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");

  //
  // s, commands for gnuplot
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) == MATRIX_DENSE_STRING)
  {
    s1 = ent_data ( e1 );
    n = SIZE(s1);
    if (n<1)
      rerror("gnucmd: string vector for first argument required");
    if (!EQVECT(s1))
      rerror("gnucmd: string vector for first argument required");
  }
  else if (ent_type(e1) == MATRIX_DENSE_REAL)
  {
    int idummy = (int) class_double(e1);
    if (idummy)
      rlab_gnuplot_idebug = 1;
    else
      rlab_gnuplot_idebug = 0;
    goto _exit_gnucmd;
  }
  else
    rerror("gnucmd: string vector for first argument required");

  if (nargs>1)
  {
    e2 = bltin_get_ent(args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      src = class_char_pointer(e2);
    }
  }
  if (isvalidstring(src)<1)
    src="gnucmd";
  
  i = GNUPLOT_DEFAULT_WINDOW;
  
  if (!i || !gnuplot_h[ i-1 ])
    rerror(THIS_SOLVER ": Horrible Internal Error: Trying to access non-window!");
  

  for (j=0; j<n; j++)
  {
    if (!MdsV0(s1,j))
      continue;
    if (!strlen(MdsV0(s1,j)))
      continue;
    // print command to gnuplot device
    fputs(MdsV0(s1,j), gnuplot_h[ i-1 ]);
    fputs("\n", gnuplot_h[ i-1 ]);
    fflush( gnuplot_h[ i-1 ]);

    // for debugging purposes print 
    if (!rlab_gnuplot_idebug)
      continue;
    fprintf(rlab_stderr, "gnuwin(%i): (%s) %s\n", i, src, MdsV0(s1,j));
  }

_exit_gnucmd:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

//
// gnuprint (x, y, z), or
// gnuprint (m)
//
Ent * ent_gnuplot_print (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int i, j, k;

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

    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      MDR * x = ent_data ( e1 );
      int nr = x->nrow;
      int nc = x->ncol;
      for (j = 0; j < nr; j++)
      {
        for (k = 0; k < nc ; k++)
        {
          if (x->type == RLAB_TYPE_INT32)
            fprintf(gnuplot_h[ i-1 ], "  %i", Mdi0(x,j,k)) ;
          else
          {
            if (isnand(Mdr0(x,j,k)))
              fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
            else
              fprintf(gnuplot_h[ i-1 ], "  %g", Mdr0(x,j,k)) ;
          }
        }
        fprintf(gnuplot_h[ i-1 ], "\n") ;
      }
    }
    else if (ent_type(e1) == MATRIX_DENSE_STRING)
    {
      MDS *s = ent_data ( e1 );
      int nr = s->nrow;
      int nc = s->ncol;
      if (nr == 1 || nc == 1)
      {
        for (j = 0; j < nr*nc; j++)
            fprintf(gnuplot_h[ i-1 ], "%s\n", MdsV0(s,j)) ;
      }
      else
      {
        for (j = 0; j < nr; j++)
        {
          for (k = 0; k < nc ; k++)
            fprintf(gnuplot_h[ i-1 ], " %s", Mds0(s,j,k)) ;

          fprintf(gnuplot_h[ i-1 ], "\n") ;
        }
      }
    }
    else
      rerror("gnuprint: second argument has to be a matrix!");
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
      rerror("gnuprint: third argument has to be a matrix!");
    MDR *z = ent_data ( e3 );

    if (MNR(z)==nr && MNC(z)==nc)
    {
      // data on the grid
      // x[i],y[j],z[i,j]=z(x[i],y[j])
      for (j = 0; j < nr; j++)
      {
        for (k = 0; k < nc; k++)
        {
          // x
          if (x->type == RLAB_TYPE_INT32)
            fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(x,j)) ;
          else
          {
            if (isnand(MdrV0(x,j)))
              fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
            else
              fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(x,j)) ;
          }
          // y
          if (y->type == RLAB_TYPE_INT32)
            fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(y,k)) ;
          else
          {
            if (isnand(MdrV0(y,k)))
              fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
            else
              fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(y,k)) ;
          }
          // z
          if (z->type == RLAB_TYPE_INT32)
            fprintf(gnuplot_h[ i-1 ], "  %i\n", Mdi0(z,j,k)) ;
          else
          {
            if (isnand(Mdr0(z,j,k)))
              fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN);
            else
              fprintf(gnuplot_h[ i-1 ], "  %g\n", Mdr0(z,j,k)) ;          }
        }
        fprintf(gnuplot_h[ i-1 ], "\n") ;
      }
      fprintf(gnuplot_h[ i-1 ], "e\n") ;
      fflush (gnuplot_h[ i-1 ]) ;
    }
    else if (MNR(z)*MNC(z)==nr && MNR(z)*MNC(z)==nc)
    {
      // scatter data set
      // x[i],y[i],z[i]=z(x[i],y[i])
      for (j = 0; j < nr; j++)
      {
        // x
        if (x->type == RLAB_TYPE_INT32)
          fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(x,j)) ;
        else
        {
          if (isnand(MdrV0(x,j)))
            fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
          else
            fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(x,j)) ;
        }
        // y
        if (y->type == RLAB_TYPE_INT32)
          fprintf(gnuplot_h[ i-1 ], "  %i", MdiV0(y,j)) ;
        else
        {
          if (isnand(MdrV0(y,j)))
            fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
          else
            fprintf(gnuplot_h[ i-1 ], "  %g", MdrV0(y,j)) ;
        }
        // z
        if (z->type == RLAB_TYPE_INT32)
          fprintf(gnuplot_h[ i-1 ], "  %i\n", MdiV0(z,j)) ;
        else
        {
          if (isnand(MdrV0(z,j)))
            fprintf(gnuplot_h[ i-1 ], " %s", RLAB_SPRINTF_NAN) ;
          else
            fprintf(gnuplot_h[ i-1 ], "  %g\n", MdrV0(z,j)) ;
        }
      }
      fprintf(gnuplot_h[ i-1 ], "e\n") ;
      fflush (gnuplot_h[ i-1 ]) ;
    }
    else
      rerror("gnuprint: the size of 'z' does not match that of 'x' and 'y'!");

  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

