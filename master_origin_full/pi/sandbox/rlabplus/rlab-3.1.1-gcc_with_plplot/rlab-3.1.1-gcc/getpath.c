#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */

#ifdef WIN32
#include <windows.h>
#endif

#include "util.h"

#define SEP '\\'		/* Normal path separator. */
#define ALTSEP '/'		/* Alternate (windows can handle it from open(). */
#define RLABPATH  ".;.\\rlib;.\\toolbox;.\\controls-toolbox"
#define LANDMARK  "bin\\rlab.exe"

#define MAXPATHLEN  256
#define DELIM ';'

static char prefix[MAXPATHLEN + 1];
static char progpath[MAXPATHLEN + 1];

/*
 * Return TRUE if ch is a separator, FALSE if not.
 */
static int
is_sep (char ch)
{
  return ch == SEP || ch == ALTSEP;
}

/*
 * Reduce / Remove the last directory from the path argument.
 * Walk the string backward until we find a path separator.
 * When we find it, mark the location as string end.
 */
static void
reduce (char *dir)
{
  int i = strlen (dir);
  while (i > 0 && !is_sep (dir[i]))
    --i;
  dir[i] = '\0';
}

/*
 * simplistic use of stat() to determine if a file exists.
 * returns 1 if it does exist
 * returns 0 if it does not exist.
 */
static int
exists (char *filename)
{
  struct stat buf;
  return stat (filename, &buf) == 0;
}

/*
 * Join two character buffers together "buffer" + "stuff"
 * buffer must be lart enough to hold both!
 * We are joining two strings to form a pathname, so insert the
 * directory separator if necessary.
 */
static void
join (char *buffer, char *stuff)
{
  int n;			/* length of buffer */
  int k;			/* length of stuff */
  if (is_sep (stuff[0]))
    n = 0;
  else
  {
    n = strlen (buffer);
    if (n > 0 && !is_sep (buffer[n - 1]) && n < MAXPATHLEN)
      buffer[n++] = SEP;
  }
  k = strlen (stuff);
  if (n + k > MAXPATHLEN)
    k = MAXPATHLEN - n;
  strncpy (buffer + n, stuff, k);
  buffer[n + k] = '\0';
}

/*
 * search argv0_path for landmark.
 * If found return: 1
 * else return:     0
 */
static int
search_for_prefix (char *argv0_path, char *landmark)
{
  int n;

  /* Search from argv0_path, until root is found */
  strcpy (prefix, argv0_path);
  do
  {
    n = strlen (prefix);
    join (prefix, landmark);
    if (exists (prefix))
    {
      prefix[n] = '\0';
      return 1;
    }
    prefix[n] = '\0';
    reduce (prefix);
  }
  while (prefix[0]);
  return 0;
}

/*
 * Get the program's path....
 * It is either specified in argv[0], or we have to find
 * it by joining elements of PATH and progname.
 */
static void
get_progpath ()
{
  char *path = getenv ("PATH");
  char *prog = get_progname ();

  if (prog == NULL || *prog == '\0' || (strlen (prog) == 3))
    prog = "rlab.exe";

  /* If there is no slash in the argv0 path, then we have to
   * assume rlab is in the user's $PATH, since there's no
   * other way to find a directory to start the search from.  If
   * $PATH isn't exported, you lose.
   */
  if (strchr (prog, SEP) || strchr (prog, ALTSEP))
  {
    /* There was a slash in the filename, use it, and return. */
    strcpy (progpath, prog);
  }
  else if (path)
  {
    while (1)			/* step through the PATH options, trying one till we find it. */
    {
      char *delim = strchr (path, DELIM);

      if (delim)
      {
	int len = delim - path;
	strncpy (progpath, path, len);
	*(progpath + len) = '\0';
      }
      else
	strcpy (progpath, path);

      join (progpath, prog);
      if (exists (progpath))
      {
	break;			/* Found it, we can leave now. */
      }
      if (!delim)
      {
	progpath[0] = '\0';
	break;
      }
      path = delim + 1;
    }
  }
  else
  {
    progpath[0] = '\0';
  }
}

/*
 * Calculate:  prefix
 *             progpath
 */
static void
calculate_path ()
{
  char argv0_path[MAXPATHLEN + 1];
  char *rlab2_home = getenv ("RLAB2_HOME");

  /* Get the path to this program. */
  get_progpath ();
  strcpy (argv0_path, progpath);
  reduce (argv0_path);
  if (rlab2_home == NULL || *rlab2_home == '\0')
  {
    if (search_for_prefix (argv0_path, LANDMARK))
      rlab2_home = prefix;
    else
      rlab2_home = NULL;
  }
  else
    strcpy (prefix, rlab2_home);
}

/* External interface. */
char *
calculate_home ()
{
  calculate_path ();
  return prefix;
}
