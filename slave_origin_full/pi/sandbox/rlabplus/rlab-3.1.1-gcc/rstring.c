/*
 * These functions are pirated from the GNU libc. We have to do
 * this because of strtoks pointer behavior.
 */

/* Copyright (C) 1991 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   
   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
   
   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.
   
   You should have received a copy of the GNU Library General Public
   License along with the GNU C Library; see the file COPYING.LIB.  If
   not, write to the Free Software Foundation, Inc., 675 Mass Ave,
   Cambridge, MA 02139, USA.  */

#include <errno.h>
#include <string.h>

#ifdef __riscos
#define EINVAL		20
#endif

static char *olds = NULL;
static size_t Rstrspn (const char *s, const char *accept);
static char *Rstrpbrk (const char *s, const char *accept);

/* Parse S into tokens separated by characters in DELIM.
   If S is NULL, the last string strtok() was called with is
   used.  For example:
        char s[] = "-abc=-def";
	x = strtok(s, "-");		// x = "abc"
	x = strtok(NULL, "=-");		// x = "def"
	x = strtok(NULL, "=");		// x = NULL
		// s = "abc\0-def\0"
*/

char *
Rstrtok (char *s, const char *delim)
{
  char *token;

  if (s == NULL)
  {
    if (olds == NULL)
    {
      errno = EINVAL;
      return NULL;
    }
    else
      s = olds;
  }

  /* Scan leading delimiters.  */
  s += Rstrspn (s, delim);
  if (*s == '\0')
  {
    olds = NULL;
    return NULL;
  }

  /* Find the end of the token.  */
  token = s;
  s = Rstrpbrk (token, delim);
  if (s == NULL)
    /* This token finishes the string.  */
    olds = NULL;
  else
  {
    /* Terminate the token and make OLDS point past it.  */
    *s = '\0';
    olds = s + 1;
  }
  return token;
}

/* Find the first ocurrence in S of any character in ACCEPT.  */

char *
Rstrpbrk (const char *s, const char *accept)
{
  while (*s != '\0')
    if (strchr (accept, *s) == NULL)
      ++s;
    else
      return (char *) s;

  return NULL;
}

/* Return the length of the maximum initial segment
   of S which contains only characters in ACCEPT.  */

static size_t
Rstrspn (const char *s, const char *accept)
{
  register const char *p;
  register const char *a;
  register size_t count = 0;

  for (p = s; *p != '\0'; ++p)
  {
    for (a = accept; *a != '\0'; ++a)
      if (*p == *a)
        break;
    if (*a == '\0')
      return count;
    else
      ++count;
  }

  return count;
}
