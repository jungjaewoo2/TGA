/* rfile.c
 * Handle RLaB rfiles */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994  Ian R. Searle

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

#include "rlab.h"
#include "symbol.h"
#include "mem.h"
#include "mds.h"
#include "util.h"
#include "rfile.h"
#include "rlab_solver_parameters_names.h"

#ifdef __STDC__
#include <stdlib.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <sys/types.h>
#ifdef __riscos
#include <swis.h>
#endif
#include <stdio.h>
#include <errno.h>
#ifdef HAVE_DIRENT_H
#include <dirent.h>
#else
#ifndef __riscos
#include <sys/dir.h>
#else
#include "riscos_ut.h"
#endif /* __riscos */
#endif

#ifdef THINK_C
#include <myconsole16.h>
extern int mac_rows;		/* no. of rows of console (defined in RLAB_ROWS env.) */
static int line_counter;
static int show = 1;		/* show = 0 if q (quit) is pressed */
#endif

/* The RLAB_SEARCH_PATH separator character */
#define RLAB_PATH_FILENAME_LEN  (256)

#ifdef unix
static char PATH_SEP = ':';
static char PATH_DELIM[] = "/";
#define cat_help_file cat_help_file_unix
#endif
#ifdef OS2
static char PATH_SEP = ';';
static char PATH_DELIM[] = "\\";
#define cat_help_file cat_help_file_unix
#endif
#ifdef DOS
static char PATH_SEP = ':';
static char PATH_DELIM[] = "\\";
#define cat_help_file cat_DOS
#endif
#ifdef THINK_C
static char PATH_SEP = ';';
static char PATH_DELIM[] = ":";
#define cat_help_file cat_help_file_unix
#endif
#ifdef __riscos
static char PATH_SEP = ';';
static char PATH_DELIM[] = ".";
#define cat_help_file cat_help_file_riscos
#define cat_help_file_riscos cat_help_file_unix
#endif

extern int rpclose _PROTO ((FILE * fp));
extern int run_program _PROTO ((char *name));

/*
 * A structure for holding the parsed up directory names.
 */

struct _path
{
  int ndir;
  char **dira;
};

typedef struct _path Path;

void print_rfiles _PROTO ((char *dirname, FILE * fptr));
int find_open _PROTO ((char *filenm, char *dirname));
Path *path_split _PROTO ((char *path));
Path *name_split _PROTO ((char *path));
void path_Destroy _PROTO ((Path * path));

static void cswap _PROTO ((char *v[], int i, int j));
void nsort _PROTO ((char *v[], int right, int left));
static int nlength _PROTO ((char *carray[], int n));
static void cat_help_file_unix _PROTO ((char *name));
#ifdef DOS
static void cat_DOS _PROTO ((char *name));
#endif

static char *find_full_filename _PROTO ((char *rfilenm, char *dirname));

static char *l_search_path;
static char *get_rlab_search_path _PROTO ((void));

#ifdef __riscos
void open_arch_dir _PROTO ((char *));
void set_dir_opener _PROTO ((char *));
#endif

/* Make this static so we can close after catching a SIGPIPE */
static FILE *pfile;

/* **************************************************************
 * Some interface functions to avoid global variables.
 * ************************************************************** */

static char *search_path;

void
set_search_path (value)
     char *value;
{
  search_path = cpstr (value);
}

static char *help_dir;

void
set_help_dir (value)
     char *value;
{
  help_dir = cpstr (value);
}

static char *lib_dir;

void
set_lib_dir (value)
     char *value;
{
  lib_dir = cpstr (value);
}

static char *pager;
static char *help_pager;

void
set_pager (value)
     char *value;
{
  pager = cpstr (value);
}

void
set_help_pager (value)
     char *value;
{
  help_pager = cpstr (value);
}

#ifdef __riscos
static char *dir_opener;

void
set_dir_opener (value)
     char *value;
{
  dir_opener = cpstr (value);
}
#endif

/* **************************************************************
 * Rfile interface to parser
 * ************************************************************** */
#if ((!defined THINKC) && (!defined __riscos))
void rfile ()
{
  int i;
  Path *path;

  if (pfile)
    rpclose (pfile);

  if ((pfile = popen (pager, "w")) == 0)
  {
    fprintf (stderr, "ERROR, could not open %s pager for write\n", pager);
    return;
  }

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  for (i = 0; i < path->ndir; i++)
  {
    fprintf (pfile, "%s :\n", path->dira[i]);
    print_rfiles (path->dira[i], pfile);
  }

  fflush (pfile);
  path_Destroy (path);
  rpclose (pfile);
  pfile = 0;
}
#else

#ifdef __riscos
void rfile ()
{
  int i;
  Path *path;


  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;

  path = path_split (l_search_path);

  for (i = 0; i < path->ndir; i++)
  {
    open_arch_dir (path->dira[i]);
  }

  path_Destroy (path);
}

#else /* THINK_C */

void rfile (void)
{
  int i;
  Path *path=0;

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  for (i = 0; i < path->ndir; i++)
  {
    show = more (stdout, "");
    if (show)
      fprintf (stdout, "%s :\n", path->dira[i]);

    print_rfiles (path->dira[i], stdout);
  }
  fflush (stdout);
  path_Destroy (path);
}

#endif /* __riscos */
#endif /* ! THINKC && !__riscos */

/*
 * Add the `.r' extension to a name if it is not
 * there. Do nothing if it is not.
 */
char * rfile_ext (char *name)
{
  char *nm, *np;
  int len=isvalidstring(name), lext=isvalidstring(RLAB_FILENAME_EXTENSION);;

  /* On the off chance someone hands a zero length name. */
  if (len < 1)
    return (0);

  /* Check for the RLAB_FILENAME_EXTENSION extension. */
#ifndef __riscos
  np = name;
  if (!strcmp(&np[len - lext], RLAB_FILENAME_EXTENSION))
  {
    /* name already has extension, do nothing. */
    return (nm = cpstr (name));
  }
  else
  {
    /* Add the `.r' to name. */
    nm = GC_MALLOC (sizeof (char) * (len + lext + 1));
    if (!nm)
      rerror ("out of memory");
    strcpy (nm, name);
    strcat (nm, RLAB_FILENAME_EXTENSION);
    nm[len + lext] = '\0';
  }
  return (nm);
#else /* __riscos */
  return (nm = cpstr (name));
#endif
}

/*
 * Interpret a rfile (load it).
 * Argument `name' is destroyed by rfile_load.
 * name can be a single name-string or white-space separated
 * string of names. Each string can have the postfix .r or not.
 */

void rfile_load (char *name)
{
  char *nm=0;
  int fval, i, j;
  Path *names=0, *path=0;

  /* Initialize */
  fval = 0;

  /* Split up name if necessary. */
  names = name_split (name);
  if (names == 0)
  {
    fprintf (stderr, "WARNING, invalid name supplied to rfile\n");
    return;
  }

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  /* Load each file in the command. */
  for (j = 0; j < names->ndir; j++)
  {
    nm = rfile_ext (names->dira[j]);
    if (nm == 0)
      continue;

    /* Check each directory in the path. */
    for (i = 0; i < path->ndir; i++)
    {
      if ((fval = find_open (nm, path->dira[i])) == 1)
      {
        break;      /* We found/loaded it. */
      }
    }

    if (fval == 0)
      fprintf (stderr, "WARNING, could not find/open: %s\n", nm);

    GC_FREE (nm);
  }

  GC_FREE (name);
  path_Destroy (names);
  path_Destroy (path);
}

void rfile_load_one (char *name)
{
  char *nm;
  int fval, i;
  Path *path;

  /* Initialize */
  fval = 0;

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  nm = rfile_ext (name);

  /* Check each directory in the path. */
  for (i = 0; i < path->ndir; i++)
  {
    if ((fval = find_open (nm, path->dira[i])) == 1)
    {
      break;			/* We found/loaded it. */
    }
  }

  if (fval == 0)
    fprintf (stderr, "WARNING, could not find/open: %s\n", nm);

  GC_FREE (nm);
  path_Destroy (path);
}

/*
 * Split a colon separated path into distinct character strings.
 * EX: .:/u1/ian/misc/rlab/examples:/usr/local/rlab
 * Return a pointer to an array of character pointers.
 */

Path *
path_split (path)
     char *path;
{
  int i, n, ndir, start;
  char **dira;
  Path *new = (Path *) GC_MALLOC (sizeof (Path));
  if (new == 0)
    rerror ("out of memory");

  ndir = 1;
  i = 0;
  /* Count the number of directories */
  while (path[i] != '\0')
  {
    if (path[i] == PATH_SEP)
      ndir++;
    i++;
  }
  /* Get space for char pts */
  dira = (char **) GC_MALLOC (ndir * sizeof (char *));
  if (dira == 0)
    rerror ("out of memory");

  /* Now get directory strings */
  start = 0;
  for (i = 0; i < ndir; i++)
  {
    /* Count length of string */
    n = 0;
    while (path[start + n] != PATH_SEP && path[start + n] != '\0')
      n++;

    /* Now copy string into dira */
    dira[i] = (char *) GC_MALLOC ((n + 1) * sizeof (char));
    if (dira[i] == 0)
      rerror ("out of memory");

    strncpy (dira[i], &path[start], n);
    dira[i][n] = '\0';
    start += (n + 1);
  }
  new->ndir = ndir;
  new->dira = dira;
  return (new);
}

/*
 * Split a whitespace separated namelist into distinct
 * character strings.
 * EX: `poly roots bode'
 * Return a pointer to an array of character pointers.
 */

Path *
name_split (name)
     char *name;
{
  int i, n, ndir, start;
  char **dira;
  Path *new;

  if (name == 0 && strlen (name) == 0)
    return (0);

  new = (Path *) GC_MALLOC (sizeof (Path));
  if (new == 0)
    rerror ("out of memory");

  ndir = 1;
  i = 0;
  /* Count the number of directories */
  while (name[i] != '\0')
  {
    if (name[i] == ' ')
      ndir++;
    i++;
  }
  /* Get space for char pts */
  dira = (char **) GC_MALLOC (ndir * sizeof (char *));
  if (dira == 0)
    rerror ("out of memory");

  /* Now get directory strings */
  start = 0;
  for (i = 0; i < ndir; i++)
  {
    /* Count length of string */
    n = 0;
    while (name[start + n] != ' ' && name[start + n] != '\0')
      n++;

    /* Now copy string into dira */
    dira[i] = (char *) GC_MALLOC ((n + 1) * sizeof (char));
    if (dira[i] == 0)
      rerror ("out of memory");

    strncpy (dira[i], &name[start], n);
    dira[i][n] = '\0';
    start += (n + 1);
  }
  new->ndir = ndir;
  new->dira = dira;
  return (new);
}

#ifdef THINK_C
static int
more (FILE * fptr, char *msg)
{
  int x, y, c;

  if (!show)
    return 0;

  fprintf (fptr, "More %s ('q' to quit) ... ", msg);
  cgetxy (&x, &y, fptr);	/* get cursor position */
  c = getchar ();

  if (c == '\r' || c == '\n')
    cgotoxy (1, y - 1, fptr);	/* goto cursor position col 1  row y-1 */
  else
    cgotoxy (1, y, fptr);
  fprintf (fptr, "                                    ");	/* erase more .. */

  if (c == '\r' || c == '\n')
    cgotoxy (1, y - 1, fptr);
  else
    cgotoxy (1, y, fptr);

  return c == 'q' ? 0 : 1;	/* return 0 if 'q' (quit) is pressed */
}

static int
page_break (FILE * fptr, int line_counter)
{
  if (line_counter % (mac_rows - 3) == 0)
    return more (fptr, "");
  else
    return show;
}

/*
 * my_pager
 * build a simple pager here
 * t.s.yang  Nov/93  ucb
 */

void
my_pager (char *name)
{
  FILE *fp;
  char line[101], msg[64];
  int counter;
  long csize, fsize;

  show = 1;
  if ((fp = fopen (name, "r")) == NULL)
  {
    fprintf (stderr, "No help available.\n");
    return;
  }
  fseek (fp, 0L, SEEK_END);
  fsize = ftell (fp);
  rewind (fp);
  csize = 0;
  counter = 0;
  while (fgets (line, 100, fp) != NULL)
  {
    fputs (line, stdout);
    counter++;
    csize += strlen (line);
    if ((counter % (mac_rows - 3)) == 0)
    {
      sprintf (msg, "(%d%%)", csize * 100 / fsize);
      if (!more (stdout, msg))
	break;
    }
  }
  fclose (fp);
}
#endif /* THINK_C */

/*
 * Print all of the rfiles in a given directory.
 */

void print_rfiles (char *dirname, FILE *fptr)
{
  char **names, temps[200];
  int i, length, ncol, nfile, swidth;
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif

  nfile = 0;

  /* Create array of rfiles names in this directory */
  if ((dirp = opendir (dirname)) != 0)
  {
    /* Count rfiles */
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (rfile_test (direntp->d_name))
        nfile++;
    }
    rewinddir (dirp);

    /* Create array of rfile names */
    if (!nfile)
    {
      fprintf (fptr, "\n");
      return;
    }
    names = (char **) GC_MALLOC (nfile * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");

    i = 0;

    while ((direntp = readdir (dirp)) != NULL)
    {
      if (rfile_test (direntp->d_name))
      {
        names[i] = cpstr (direntp->d_name);

        /* remove `.r' */
        names[i][strlen (names[i]) - 2] = '\0';
        i++;
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", dirname);
    return;
  }

  /* Sort rfile names */
  nsort (names, 0, nfile - 1);

  length = nlength (names, nfile);

#ifdef HAVE_GETENV
  sprintf (temps, "%s", getenv ("COLUMNS"));
  swidth = strtod (temps, NULL) - 8;
  if (swidth <= 0)
    swidth = TERM_WIDTH;
#else
  swidth = TERM_WIDTH;
#endif

  //ncol = TERM_WIDTH / (length + 2);
  ncol = swidth / (length + 2);

  /* Now print out rfile names */
  for (i = 1; i <= nfile; i++)
  {

#ifndef THINK_C
    fprintf (fptr, "%*s ", length + 2, names[i - 1]);
#else
    if (show)
      fprintf (fptr, "%*s ", length + 2, names[i - 1]);
#endif /* ! THINK_C */

    if ((i % ncol) == 0)
    {

#ifndef THINK_C
      fprintf (fptr, " \n");
      fflush (fptr);
#else
      if (show)
        fprintf (fptr, " \n");
      fflush (fptr);
      line_counter++;
      show = page_break (fptr, line_counter);
#endif /* ! THINK_C */

    }
    GC_FREE (names[i - 1]);
  }

#ifndef THINK_C
  fprintf (fptr, "\n\n");
#else
  if (show)
    fprintf (fptr, "\n\n");
#endif /* ! THINK_C */

  fflush (fptr);

  GC_FREE (names);
}

/*
 * Given a potential file name, and a directory name,
 * search the directory for the file. If it exists,
 * call new_file().
 *
 * Return TRUE if successful, FALSE if not.
 */

int find_open (char *filenm, char *dirname)
{
  /* FIX */
  char tmp[RLAB_PATH_FILENAME_LEN];
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif

  if ((dirp = opendir (dirname)) != 0)
  {
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (!strcmp (direntp->d_name, filenm))
      {
        strcpy (tmp, dirname);
        strcat (tmp, PATH_DELIM);
        strcat (tmp, direntp->d_name);
        closedir (dirp);  /* just in case we error */
        run_program (tmp);
        return (1);
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", dirname);
  }

  return (0);     /* File not found */
}

/*
 * Test if filename is a valid rfile.
 * Must end in a `.r'
 */

int rfile_test (char *filename)
{
  int len = isvalidstring (filename);
  int lex = isvalidstring (RLAB_FILENAME_EXTENSION);

  if (len<lex)
    return (0);

  /* Check the last two characters */
  if (!strcmp(&filename[len - lex], RLAB_FILENAME_EXTENSION))
    return (1);

  return (0);
}

/*
 * Destroy a Path struct.
 */

void
path_Destroy (path)
     Path *path;
{
  int i;

  for (i = 0; i < path->ndir; i++)
    GC_FREE (path->dira[i]);

  GC_FREE (path->dira);
  path->ndir = 0;
  GC_FREE (path);
}

/*
 * Simple character qsort.
 */

void
nsort (v, left, right)
     char *v[];
     int left, right;
{
  int i, last;

  if (left >= right)
    return;
  cswap (v, left, (left + right) / 2);
  last = left;
  for (i = left + 1; i <= right; i++)
    if (strcmp (v[i], v[left]) < 0)
      cswap (v, ++last, i);
  cswap (v, left, last);
  nsort (v, left, last - 1);
  nsort (v, last + 1, right);
}

/*
 * Interchange v[i] and v[j]
 */

static void
    cswap (v, i, j)
     char *v[];
     int i, j;
{
  char *temp;

  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}

/*
 * Find the length of the longest string in the array.
 */

static int
nlength (carray, n)
     char *carray[];
     int n;
{
  int i, tmp, length;

  length = strlen (carray[0]);
  for (i = 1; i < n; i++)
  {
    tmp = strlen (carray[i]);
    if (tmp > length)
      length = tmp;
  }
  return (length);
}

/* **************************************************************
 * Help related functions
 * ************************************************************** */

/*
 * List contents of help directory,, and directories in
 * RLAB_SEARCH_PATH.
 */

#if ((! defined THINK_C) && (!defined __riscos))
void
help ()
{
  char **names, temps[200];
  int i, length, ncol, nfile, swidth;
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif
  Path *path;

  nfile = 0;

  if (pfile)
    rpclose (pfile);

  if ((pfile = popen (help_pager, "w")) == 0)
  {
    fprintf (stderr, "ERROR, could not open %s pager for write\n", pager);
    return;
  }

  /*
   * Print out the RLAB_HELP_DIR files 1st
   */

  if ((dirp = opendir (help_dir)) != 0)
  {
    /* Count ALL files */
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (direntp->d_name[0] != '.')
	nfile++;
    }
    rewinddir (dirp);

    /* Create array of rfile names */
    names = (char **) GC_MALLOC (nfile * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");

    i = 0;

    while ((direntp = readdir (dirp)) != NULL)
    {
      if (direntp->d_name[0] != '.')
      {
	names[i] = cpstr (direntp->d_name);
	i++;
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", help_dir);
    return;
  }

  /* Sort rfile names */
  nsort (names, 0, nfile - 1);

  length = nlength (names, nfile);

#ifdef HAVE_GETENV
  sprintf (temps, "%s", getenv ("COLUMNS"));
  swidth = strtod (temps, NULL) - 8;
  if (swidth <= 0)
    swidth = TERM_WIDTH;
#else
    swidth = TERM_WIDTH;
#endif

  ncol = swidth / (length + 2);

  /* Now print out rfile names */

  for (i = 1; i <= nfile; i++)
  {
    fprintf (pfile, "%*s ", length + 2, names[i - 1]);
    if ((i % ncol) == 0)
      fprintf (pfile, "\n");
    GC_FREE (names[i - 1]);
  }

  fprintf (pfile, "\n");
  fflush (pfile);

  GC_FREE (names);

  /*
   * Now do the rfiles in RLAB_LIB_DIR.
   */

  fprintf (pfile, "\n%s :\n", lib_dir);
  print_rfiles (lib_dir, pfile);

  /*
   * Now do the same for the directories in RLAB_SEARCH_PATH
   */

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  for (i = 0; i < path->ndir; i++)
  {
    fprintf (pfile, "%s :\n", path->dira[i]);
    print_rfiles (path->dira[i], pfile);
  }
  fflush (pfile);
  path_Destroy (path);

  rpclose (pfile);
  pfile = 0;
}

#else

#ifdef __riscos

void
help ()
{
  int i;
  Path *path;
  path = path_split (help_dir);
  for (i = 0; i < path->ndir; i++)
  {
    open_arch_dir (path->dira[i]);
  }
  path = path_split (lib_dir);
  for (i = 0; i < path->ndir; i++)
  {
    open_arch_dir (path->dira[i]);
  }
}

#else /* THINK_C */

void
help ()
{
  char **names, temps[200];
  int i, length, ncol, nfile, swidth;
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif
  Path *path;

  nfile = 0;
  pfile = stdout;

  /*
   * Print out the RLAB_HELP_DIR files 1st
   */

  if ((dirp = opendir (help_dir)) != 0)
  {
    /* Count ALL files */
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (direntp->d_name[0] != '.')
	nfile++;
    }
    rewinddir (dirp);

    /* Create array of rfile names */
    names = (char **) GC_MALLOC (nfile * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");

    i = 0;

    while ((direntp = readdir (dirp)) != NULL)
    {
      if (direntp->d_name[0] != '.')
      {
	names[i] = cpstr (direntp->d_name);
	i++;
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", help_dir);
    return;
  }

  /* Sort rfile names */
  nsort (names, 0, nfile - 1);

  length = nlength (names, nfile);

#ifdef HAVE_GETENV
  sprintf (temps, "%s", getenv ("COLUMNS"));
  swidth = strtod (temps, NULL) - 8;
  if (swidth <= 0)
    swidth = TERM_WIDTH;
#else
  swidth = TERM_WIDTH;
#endif
  ncol = swidth / (length + 2);

  /* Now print out rfile names */
  line_counter = 0;
  show = 1;
  for (i = 1; i <= nfile; i++)
  {
    if (show)
      fprintf (stdout, "%*s ", length + 2, names[i - 1]);
    if ((i % ncol) == 0)
    {
      if (show)
	fprintf (stdout, "\n");
      line_counter++;
      show = page_break (stdout, line_counter);
    }
    GC_FREE (names[i - 1]);
  }
  if (show)
    fprintf (stdout, "\n");
  fflush (stdout);
  GC_FREE (names);

  /*
   * Now do the rfiles in RLAB_LIB_DIR.
   */

  show = more (stdout, "");
  if (show)
    fprintf (stdout, "\n%s :\n", lib_dir);
  line_counter = 3;
  print_rfiles (lib_dir, stdout);

  /*
   * Now do the same for the directories in RLAB_SEARCH_PATH
   */

  /* Re-split the PATH, in case it has changed */
  path = path_split (search_path);
  for (i = 0; i < path->ndir; i++)
  {
    show = more (stdout, "");
    if (show)
      fprintf (stdout, "\n%s :\n", path->dira[i]);
    line_counter = 3;
    print_rfiles (path->dira[i], stdout);
  }
  fflush (stdout);
  path_Destroy (path);
}

#endif /* __riscos */
#endif /* ((! defined THINK_C) && (!defined __riscos)) */

/*
 * Search the help directory, and then RLAB_SEARCH_PATH
 * for a file identified by `name'. Print out the
 * "help contents".
 */

void
help_name (name)
     char *name;
{
  char *tmp;
  int i, length;
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif
  Path *path;

  /*
   * Search the default help directory 1st
   */

  /* Create array of rfiles names in this directory */
#ifndef __riscos
  if ((dirp = opendir (help_dir)) != 0)
  {
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (!strcmp (name, direntp->d_name))
      {
	length = strlen (help_dir) + strlen (direntp->d_name) + 2;
	tmp = (char *) GC_MALLOC (length * sizeof (char));
	if (tmp == 0)
	  rerror ("out of memory");

	strcpy (tmp, help_dir);
	strcat (tmp, PATH_DELIM);
	strcat (tmp, direntp->d_name);
	cat_help_file (tmp);
	closedir (dirp);
	GC_FREE (tmp);
	GC_FREE (name);
	return;
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", help_dir);
  }

  /*
   * Next, search the RLAB_LIB_DIR.
   */

  if ((tmp = find_full_filename (name, lib_dir)) != 0)
  {
    cat_help_file (tmp);
    GC_FREE (name);
    GC_FREE (tmp);
    return;
  }

  /*
   * Search RLAB_SEARCH_PATH
   */

  /* Re-split the PATH, in case it has changed */
  if ((l_search_path = get_rlab_search_path ()) == 0)
    l_search_path = search_path;
  path = path_split (l_search_path);

  for (i = 0; i < path->ndir; i++)
  {
    if ((tmp = find_full_filename (name, path->dira[i])) != 0)
    {
      cat_help_file (tmp);
      GC_FREE (name);
      GC_FREE (tmp);
      path_Destroy (path);
      return;
    }
  }
  fprintf (stderr, "ERROR, could not find/open: %s\n", name);
  GC_FREE (name);
  path_Destroy (path);

#else /*__riscos*/
  path = path_split (help_dir);
  for (i = 0; i < path->ndir; i++)
  {
    length = strlen (path->dira[i]) + strlen (name) + 2;
    tmp = (char *) GC_MALLOC (length * sizeof (char));
    strcpy (tmp, path->dira[i]);
    strcat (tmp, PATH_DELIM);
    strcat (tmp, name);
    cat_help_file (tmp);
    GC_FREE (tmp);
  }
  GC_FREE (name);
#endif /* ! _riscos */
}

/*
 * Cat a help file to stdout.
 */

static void
cat_help_file_unix (name)
     char *name;
{
#ifdef THINK_C
  my_pager (name);
#else
  char *sys_string;
  int length;

  /* Build system string */
  length = strlen (help_pager) + strlen (name) + 2;
  sys_string = (char *) GC_MALLOC (length * sizeof (char));
  if (sys_string == 0)
    rerror ("out of memory");

  sprintf (sys_string, "%s %s", help_pager, name);
  system (sys_string);
  GC_FREE (sys_string);
#endif /* THINK_C */
}

/*
 * A simple builtin cat for OSes that
 * cannot handle more than one process.
 * See The ANSI-C Programming Language.
 */

#ifdef DOS
static void filecat _PROTO ((FILE * fin));
static void
cat_DOS (char *name)
{
  FILE *fp;
  if ((fp = fopen (name, "r")) == 0)
  {
    fprintf (stderr, "cannot open file %s for write\n", name);
    rerror ("help file error");
  }
  else
  {
    filecat (fp);
    fclose (fp);
  }
}

static void
filecat (fn)
     FILE *fn;
{
  int c;
  while ((c = getc (fn)) != EOF)
  {
    putc (c, stdout);
  }
}
#endif

/*
 * Given a partial rfile name, return the full
 * pathname to that file.
 */

static char *
find_full_filename (rfilenm, dirname)
     char *rfilenm, *dirname;
{
  /* FIX */
  char tmp[RLAB_PATH_FILENAME_LEN], *fullname;
  int length;
  DIR *dirp;
#ifdef HAVE_DIRENT_H
  struct dirent *direntp;
#else
  struct direct *direntp;
#endif

  strcpy (tmp, rfilenm);
  strcat (tmp, ".r");

  if ((dirp = opendir (dirname)) != 0)
  {
    while ((direntp = readdir (dirp)) != NULL)
    {
      if (!strcmp (direntp->d_name, tmp))
      {
	length = strlen (dirname) + strlen (direntp->d_name) + 2;
	fullname = (char *) GC_MALLOC (length * sizeof (char));
	if (fullname == 0)
	  rerror ("out of memory");

	strcpy (fullname, dirname);
	strcat (fullname, PATH_DELIM);
	strcat (fullname, direntp->d_name);
	closedir (dirp);
	return (fullname);
      }
    }
    closedir (dirp);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: %s\n", dirname);
  }
  return (0);
}

/*
 * Check the global variable _rlab_search_path
 * for a new path string.
 */

static char *
get_rlab_search_path ()
{
  Ent *ent;
  ListNode *lpath;
  char *path;

  if ((lpath = btree_FindNode (get_symtab_ptr (), "_rlab_search_path")) == 0)
    return (0);

  ent = var_ent (lpath);
  if (ent_type (ent) == MATRIX_DENSE_STRING)
  {
    path = cpstr (mds_GetString (ent_data (ent)));
    return (path);
  }

  return (0);
}

/*
 * Make sure there is no .r
 */

char *
no_ext (name)
     char *name;
{
  char *nm, *np;
  int len = strlen (name);

  /* On the off chance someone hands a zero length name. */
  if (len == 0)
    return (0);

  /* Check for the `.r' extension. */
  np = name;
  if (np[len - 1] == 'r' && np[len - 2] == '.')
  {
    /* name already has extension, strip it. */
    nm = GC_MALLOC (sizeof (char) * len - 1);
    if (nm == 0)
      rerror ("out of memory");

    strncpy (nm, name, len - 2);
    return (nm);
  }
  else
  {
    /* Do nothing. */
    return (nm = cpstr (name));
  }
  return (nm);
}

void
require (name)
     char *name;
{
  char *func_name;
  int i;
  ListNode *ent;
  Path *names;

  names = name_split (name);
  if (names == 0)
  {
    fprintf (stderr, "WARNING, invalid name supplied to require\n");
    return;
  }

  for (i = 0; i < names->ndir; i++)
  {
    /* Get the function name. */
    func_name = no_ext (names->dira[i]);

    /* Now look for it in the symbol table. */
    if ((ent = lookup (0, func_name)) != 0)
    {
      /* Now we must check to see if there is something there */
      if (ent_type (var_ent (ent)) == UNDEF)
      {
	/* Load the file. */
	rfile_load_one (names->dira[i]);
      }
      /* Else, do nothing. */
    }
    else
    {
      /* Load the file. */
      rfile_load_one (names->dira[i]);
    }
    GC_FREE (func_name);
  }

  GC_FREE (name);
  path_Destroy (names);
}

#ifdef __riscos
void
open_arch_dir (name)
     char *name;
{
  char *sys_string;
  int length;

  /* Build system string */
  length = strlen (dir_opener) + strlen (name) + 2;
  sys_string = (char *) GC_MALLOC (length * sizeof (char));

  sprintf (sys_string, "%s %s", dir_opener, name);
  system (sys_string);

  GC_FREE (sys_string);
}
#endif

/* **************************************************************
 * Run a whole directory of rfiles through the machine.
 * ************************************************************** */

static DIR *dirp;
#ifdef HAVE_DIRENT_H
static struct dirent *direntp;
#else
static struct direct *direntp;
#endif

int rfile_dir (char *p)
{
  /* FIX */
  char tmp[RLAB_PATH_FILENAME_LEN], *dirname=p;
  int ldir = isvalidstring(dirname), found_colon=0;

  if (ldir<1)
    return (0);

  while(dirname)
  {
    // temporarily terminate string at ':'
    char *c = strchr(dirname, PATH_SEP);
    if (c)
    {
      found_colon=1;
      *c = '\0';
    }

    // check directory name for '::'
    if (isvalidstring(dirname) > 0)
    {
      if ((dirp = opendir (dirname)) != 0)
      {
        while ((direntp = readdir (dirp)) != NULL)
        {
          if (rfile_test (direntp->d_name))
          {
            strcpy (tmp, dirname);
            strcat (tmp, PATH_DELIM);
            strcat (tmp, direntp->d_name);
            run_program (tmp);
          }
        }
        closedir (dirp);
      }
    }

    // recover string before ':' was removed and move one place right from it
    if(found_colon)
    {
      found_colon=0;
      *c = PATH_SEP;
      c++;
    }
    dirname = c;
  }

  return (0);
}
