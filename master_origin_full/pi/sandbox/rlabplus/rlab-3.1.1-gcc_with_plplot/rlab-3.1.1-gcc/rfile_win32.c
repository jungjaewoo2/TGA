/* rfile_win32.c
 * Handle RLaB rfiles */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1997  Ian R. Searle

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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <io.h>

/* The RLAB_SEARCH_PATH separator character */
static char PATH_SEP = ';';
static char PATH_DELIM[] = "\\";

extern int rpclose (FILE * fp);
extern int run_program (char *name);

/*
 * A structure for holding the parsed up directory names.
 */

struct _path
{
  int ndir;
  char **dira;
};

typedef struct _path Path;

int rfile_test (char *filename);
void print_rfiles (char *dirname, FILE * fp);
int find_open (char *filenm, char *dirname);
Path *path_split (char *path);
Path *name_split (char *path);
void path_Destroy (Path * path);

static void cat_help_file (char *fname, FILE * fp);
static void swap (char *v[], int i, int j);
void nsort (char *v[], int right, int left);
static int nlength (char *carray[], int n);
static char *find_full_filename (char *rfilenm, char *dirname);

static char *l_search_path;
static char *get_rlab_search_path (void);

/* Make this static so we can close after catching a SIGPIPE */
static FILE *pfile;

/* **************************************************************
 * Some interface functions to avoid global variables.
 * ************************************************************** */

static char *search_path;

void
set_search_path (char *value)
{
  search_path = cpstr (value);
}

static char *help_dir;

void
set_help_dir (char *value)
{
  help_dir = cpstr (value);
}

static char *lib_dir;

void
set_lib_dir (char *value)
{
  lib_dir = cpstr (value);
}

static char *pager;
static char *help_pager;

void
set_pager (char *value)
{
  pager = cpstr (value);
}

void
set_help_pager (char *value)
{
  help_pager = cpstr (value);
}

/* **************************************************************
 * Rfile interface to parser
 * ************************************************************** */

void
rfile (void)
{
  int i;
  Path *path;

  /* Open a pipe to more < */
  if (pfile)
    rpclose (pfile);

  if ((pfile = popen (pager, "w")) == 0)
  {
    fprintf (stderr, "ERROR, could not open RLAB2_PAGER pager for write\n");
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

  path_Destroy (path);
  fflush (pfile);
  rpclose (pfile);
  pfile = 0;
}

/*
 * Add the `.r' extension to a name if it is not
 * there. Do nothing if it is not.
 */

char *
rfile_ext (name)
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
    /* name already has extension, do nothing. */
    return (nm = cpstr (name));
  }
  else
  {
    /* Add the `.r' to name. */
    nm = GC_MALLOC (sizeof (char) * len + 3);
    if (nm == 0)
      rerror ("out of memory");
    strcpy (nm, name);
    strcat (nm, ".r");
    nm[len + 2] = '\0';
  }
  return (nm);
}

/*
 * Interpret a rfile (load it).
 * Argument `name' is destroyed by rfile_load.
 * name can be a single name-string or white-space separated
 * string of names. Each string can have the postfix .r or not.
 */

void
rfile_load (name)
     char *name;
{
  char *nm;
  int fval, i, j;
  Path *names, *path;

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
	break;			/* We found/loaded it. */
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

void
rfile_load_one (name)
     char *name;
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

/*
 * Print all of the rfiles in a given directory.
 */

void
print_rfiles (char *dirname, FILE * fp)
{
  char *fpathname, **names;
  int i, length, ncol, nfile, size;
  struct _finddata_t fileinfo;
  long int lhandle;
  int handle;

  nfile = 0;

  /* Create the full pathname. */
  size = strlen (dirname) + 5;
  fpathname = (char *) GC_MALLOC (size * sizeof (char));
  strcpy (fpathname, dirname);

  /* Get the rest in the hard way. */
  fpathname[size - 5] = '\\';
  fpathname[size - 4] = '*';
  fpathname[size - 3] = '.';
  fpathname[size - 2] = 'r';
  fpathname[size - 1] = '\0';

  /* Create array of rfiles names in this directory */
  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    nfile++;

    /* Count rfiles */
    while ((handle = _findnext (lhandle, &fileinfo)) != -1)
    {
      nfile++;
    }
    _findclose (lhandle);

    /* Create array of rfile names */
    if (!nfile)
    {
      fprintf (fp, "\n");
      return;
    }
    names = (char **) GC_MALLOC (nfile * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");

    i = 0;
    if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
    {
      names[0] = cpstr (fileinfo.name);

      /* remove `.r' */
      names[0][strlen (names[0]) - 2] = '\0';
      i++;

      while ((handle = _findnext (lhandle, &fileinfo)) != -1)
      {
	names[i] = cpstr (fileinfo.name);

	/* remove `.r' */
	names[i][strlen (names[i]) - 2] = '\0';
	i++;
      }
    }
    _findclose (lhandle);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: .\n");
    return;
  }

  /* Sort rfile names */
  nsort (names, 0, nfile - 1);

  length = nlength (names, nfile);
  ncol = TERM_WIDTH / (length + 2);

  /* Now print out rfile names */
  for (i = 1; i <= nfile; i++)
  {
    fprintf (fp, "%*s ", length + 2, names[i - 1]);

    if ((i % ncol) == 0)
    {
      fprintf (fp, " \n");
      fflush (fp);
    }
    GC_FREE (names[i - 1]);
  }

  fprintf (fp, "\n\n");
  fflush (fp);

  GC_FREE (fpathname);
  GC_FREE (names);
}

/*
 * Given a potential file name, and a directory name,
 * search the directory for the file. If it exists,
 * call new_file().
 *
 * Return TRUE if successful, FALSE if not.
 */

int
find_open (char *filenm, char *dirname)
{
  /* FIX */
  char *fpathname;
  int size_dir, size_fn;
  struct _finddata_t fileinfo;
  long int lhandle;

  /* Fix up file specification. */
  size_dir = strlen (dirname);
  size_fn = strlen (filenm);
  fpathname = (char *) GC_MALLOC ((size_dir + size_fn + 2) * sizeof (char));

  strcpy (fpathname, dirname);
  strcat (fpathname, PATH_DELIM);
  strcat (fpathname, filenm);
  fpathname[size_dir + size_fn + 1] = '\0';

  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    run_program (fpathname);
    _findclose (lhandle);
    GC_FREE (fpathname);
    return (1);
  }
  else
  {
    return (0);			/* File not found */
  }
  return (0);			/* Shut up the compiler. */
}

/* **************************************************************
 * Run a whole directory of rfiles through the machine.
 * ************************************************************** */

int
rfile_dir (char *dirname)
{
  /* FIX */
  struct _finddata_t fileinfo;
  long int lhandle;
  char *fpathname, *fname;
  int handle;
  int size;

  /* Create the full pathname. */
  size = strlen (dirname) + 5;
  fpathname = (char *) GC_MALLOC (size * sizeof (char));
  strcpy (fpathname, dirname);

  /* Get the rest in the hard way. */
  fpathname[size - 5] = '\\';
  fpathname[size - 4] = '*';
  fpathname[size - 3] = '.';
  fpathname[size - 2] = 'r';
  fpathname[size - 1] = '\0';

  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    size = strlen (dirname) + strlen (fileinfo.name);
    fname = (char *) GC_MALLOC ((size + 2) * sizeof (char));
    strcpy (fname, dirname);
    strcat (fname, PATH_DELIM);
    strcat (fname, fileinfo.name);
    run_program (fname);
    GC_FREE (fname);

    while ((handle = _findnext (lhandle, &fileinfo)) != -1)
    {
      size = strlen (dirname) + strlen (fileinfo.name);
      fname = (char *) GC_MALLOC ((size + 2) * sizeof (char));
      strcpy (fname, dirname);
      strcat (fname, PATH_DELIM);
      strcat (fname, fileinfo.name);
      run_program (fname);
      GC_FREE (fname);
    }
    _findclose (lhandle);
    return (1);
  }
  else
  {
    fprintf (stderr, "\tERROR opening directory: %s\n", dirname);
    return (0);			/* File not found */
  }
  return (0);			/* Shut up the compiler. */
}

/*
 * Test if filename is a valid rfile.
 * Must end in a `.r'
 */

int
rfile_test (char *filename)
{
  int len;
  len = strlen (filename);

  /* Check the last two characters */
  if (filename[len - 1] == 'r' && filename[len - 2] == '.')
    return (1);

  return (0);
}

/*
 * Destroy a Path struct.
 */

void
path_Destroy (Path * path)
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
nsort (char *v[], int left, int right)
{
  int i, last;

  if (left >= right)
    return;
  swap (v, left, (left + right) / 2);
  last = left;
  for (i = left + 1; i <= right; i++)
  {
    if (strcmp (v[i], v[left]) < 0)
      swap (v, ++last, i);
  }

  swap (v, left, last);
  nsort (v, left, last - 1);
  nsort (v, last + 1, right);
}

/*
 * Interchange v[i] and v[j]
 */

static void
swap (char *v[], int i, int j)
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

void
help (void)
{
  char *fpathname, **names;
  int i, length, ncol, nfile, size;
  struct _finddata_t fileinfo;
  long int lhandle;
  int handle;
  Path *path;

  nfile = 0;

  /* Open a pipe to more < */
  if (pfile)
    rpclose (pfile);

  if ((pfile = popen (help_pager, "w")) == 0)
  {
    fprintf (stderr, "ERROR, could not open RLAB2_HELP_PAGER for write\n");
    return;
  }

  /* Create the full pathname. */
  size = strlen (help_dir) + 5;
  fpathname = (char *) GC_MALLOC (size * sizeof (char));
  strcpy (fpathname, help_dir);

  /* Get the rest in the hard way. */
  fpathname[size - 5] = '\\';
  fpathname[size - 4] = '*';
  fpathname[size - 3] = '.';
  fpathname[size - 2] = '*';
  fpathname[size - 1] = '\0';

  /* Create array of rfiles names in this directory */
  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    nfile++;

    /* Count help files */
    while ((handle = _findnext (lhandle, &fileinfo)) != -1)
    {
      nfile++;
    }
    _findclose (lhandle);

    /* Create array of rfile names */
    if (!nfile)
    {
      fprintf (pfile, "\n");
      return;
    }
    names = (char **) GC_MALLOC (nfile * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");

    i = 0;
    if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
    {
      names[0] = cpstr (fileinfo.name);
      i++;

      while ((handle = _findnext (lhandle, &fileinfo)) != -1)
      {
	names[i] = cpstr (fileinfo.name);
	i++;
      }
    }
    _findclose (lhandle);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: .\n");
    return;
  }

  /* Sort rfile names */
  nsort (names, 0, nfile - 1);

  length = nlength (names, nfile);
  ncol = TERM_WIDTH / (length + 2);

  /* Now print out rfile names */
  for (i = 1; i <= nfile; i++)
  {
    fprintf (pfile, "%*s ", length + 2, names[i - 1]);

    if ((i % ncol) == 0)
    {
      fprintf (pfile, " \n");
      fflush (pfile);
    }
    GC_FREE (names[i - 1]);
  }

  fprintf (pfile, "\n");
  fflush (pfile);

  GC_FREE (fpathname);
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
  rpclose (pfile);
  path_Destroy (path);
}

/*
 * Search the help directory, and then RLAB_SEARCH_PATH
 * for a file identified by `name'. Print out the
 * "help contents".
 */

void
help_name (char *name)
{
  char *fpathname, *tmp;
  int i, length, nfile, size;
  struct _finddata_t fileinfo;
  long int lhandle;
  int handle;
  Path *path;

  nfile = 0;

  /* Open a pipe to MORE < */
  if (pfile)
    rpclose (pfile);

  if ((pfile = popen (help_pager, "w")) == 0)
  {
    fprintf (stderr, "ERROR, could not open RLAB2_HELP_PAGER for write\n");
    return;
  }

  /* Create the full pathname. */
  size = strlen (help_dir) + 5;
  fpathname = (char *) GC_MALLOC (size * sizeof (char));
  strcpy (fpathname, help_dir);

  /* Get the rest in the hard way. */
  fpathname[size - 5] = '\\';
  fpathname[size - 4] = '*';
  fpathname[size - 3] = '.';
  fpathname[size - 2] = '*';
  fpathname[size - 1] = '\0';

  /* Search for the named helpfile. */
  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    /* Check 1st file (what are the odds). */
    if (!strcmp (name, fileinfo.name))
    {
      length = strlen (help_dir) + strlen (fileinfo.name) + 2;
      tmp = (char *) GC_MALLOC (length * sizeof (char));
      if (tmp == 0)
	rerror ("out of memory");

      strcpy (tmp, help_dir);
      strcat (tmp, PATH_DELIM);
      strcat (tmp, fileinfo.name);
      cat_help_file (tmp, pfile);
      rpclose (pfile);
      _findclose (lhandle);
      GC_FREE (tmp);
      GC_FREE (name);
      return;
    }

    /* Now, check the rest of the files in the directory. */
    while ((handle = _findnext (lhandle, &fileinfo)) != -1)
    {
      if (!strcmp (name, fileinfo.name))
      {
	length = strlen (help_dir) + strlen (fileinfo.name) + 2;
	tmp = (char *) GC_MALLOC (length * sizeof (char));
	if (tmp == 0)
	  rerror ("out of memory");

	strcpy (tmp, help_dir);
	strcat (tmp, PATH_DELIM);
	strcat (tmp, fileinfo.name);
	cat_help_file (tmp, pfile);
	rpclose (pfile);
	_findclose (lhandle);
	GC_FREE (tmp);
	GC_FREE (name);
	return;
      }
    }
    _findclose (lhandle);
  }
  else
  {
    fprintf (stderr, "ERROR opening directory: .\n");
    return;
  }

  /*
   * Next, search the RLAB_LIB_DIR.
   */

  if ((tmp = find_full_filename (name, lib_dir)) != 0)
  {
    cat_help_file (tmp, pfile);
    rpclose (pfile);
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
      cat_help_file (tmp, pfile);
      rpclose (pfile);
      GC_FREE (name);
      GC_FREE (tmp);
      path_Destroy (path);
      return;
    }
  }
  fprintf (stderr, "ERROR, could not find/open: %s\n", name);
  rpclose (pfile);
  GC_FREE (name);
  path_Destroy (path);
}

/*
 * Given a partial rfile name, return the full
 * pathname to that file.
 */

static char *
find_full_filename (char *rfilenm, char *dirname)
{
  /* FIX */
  char *fpathname, tmp[256];
  int size_dir, size_fn;
  struct _finddata_t fileinfo;
  long int lhandle;
  int handle;

  /* Fix up file specification. */
  strcpy (tmp, rfilenm);
  strcat (tmp, ".r");

  size_dir = strlen (dirname);
  size_fn = strlen (tmp);
  fpathname = (char *) GC_MALLOC ((size_dir + size_fn + 2) * sizeof (char));

  strcpy (fpathname, dirname);
  strcat (fpathname, PATH_DELIM);
  strcat (fpathname, tmp);
  fpathname[size_dir + size_fn + 1] = '\0';

  if ((lhandle = _findfirst (fpathname, &fileinfo)) != -1)
  {
    if (!strcmp (tmp, fileinfo.name))
    {
      _findclose (lhandle);
      return (fpathname);
    }
    while ((handle = _findnext (lhandle, &fileinfo)) != -1)
    {
      if (!strcmp (tmp, fileinfo.name))
      {
	_findclose (lhandle);
	return (fpathname);
      }
    }
  }
  else
  {
    return (0);			/* File not found */
  }
  return (0);			/* Shut up the compiler. */
}

/*
 * Check the global variable _rlab_search_path 
 * for a new path string.
 */

static char *
get_rlab_search_path (void)
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
no_ext (char *name)
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
require (char *name)
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

/*
 * Print a file to a FILE descriptor, character at a time.
 */

static void
cat_help_file (char *fname, FILE * fp)
{
  int ichar;
  FILE *fnamep;

  /* Open the file for read. */
  if ((fnamep = fopen (fname, "r")) == 0)
  {
    fprintf (stderr, "ERROR, could not open %s for read\n", fname);
    return;
  }

  /* Read the file, and put characters to FILE *fp */
  while ((ichar = fgetc (fnamep)) != EOF)
  {
    fputc (ichar, fp);
  }
  fputc ('\n', fp);
  fflush (fp);

  /* Close the help file. */
  fclose (fnamep);
}
