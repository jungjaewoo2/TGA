/* misc.c */

#include "rlab.h"

#ifndef HAVE_RINDEX

/* **************************************************************
 * A BSD index, and rindex look-alike.
 * ************************************************************** */
/*
   #include <string.h>

   char *
   index (s, c)
   char *s, c;
   {
   return (strchr (s, (int) c));
   }

   char *
   rindex (s, c)
   char *s, c;
   {
   return (strrchr (s, (int) c));
   }
 */
#endif /* HAVE_RINDEX */
