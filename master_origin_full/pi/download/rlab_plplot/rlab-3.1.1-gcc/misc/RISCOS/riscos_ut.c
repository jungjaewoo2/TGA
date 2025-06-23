/* riscos_ut.c */

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
   RISC OS note:
   what is this file for??
   Well, it is a collection of unrelated functions neede here and there
   in the RISC OS version of RLAB
   It is quite messy.
   A few functions come from Niklas Rojemo (throwback)
   Others (unix like stuff, sys/dir.h) Gustav (pmoore@cix.compulink.co.uk)
   Some is my own stuff

   ********************************************************************** */

#ifdef __riscos

#include <string.h>
#include <setjmp.h>
#include <stddef.h>
#include <swis.h>
#include <stdlib.h>

#include "rlab.h"
#include "util.h"
#include "riscos_ut.h"


/* throwback handling - begin */
static int     ThrowbackStarted = 0;
static char    *Filename;
static char    *ErrorFile;
/* throwback handling - end */

/* pipe from pipe file - begin */
static int     pipe_on = 0;
static char    *rlab_pipe_in;
/* pipe from pipe file - end */

/* filetype checking - begin*/
static char    *rlab_filetype;
/* filetype checking - end*/


void
throwback_init(int throwback, char *np)
{
        if (throwback)
                Filename = np;
        else
                Filename = 0;
#ifdef DEBUG
        printf("filename = %s\n", np);
#endif
}

void
throwback_finish(void)
{
        if (ThrowbackStarted > 0) {
                ThrowbackStarted = 0;
                ThrowbackEnd();
        }
}

void
throwback(int level, int lineno, char *error)
{
        os_error       *err;
        if (!Filename)
                return;
        if (!ThrowbackStarted) {
                err = ThrowbackStart();
                if (err) {
                        fprintf(stderr, "ThrowbackStart %s\n", err->errmess);
                        ThrowbackStarted = 0;
                        return;
/*                        exit(-1);*/
                }
                err = ThrowbackSendStart(Filename);
                if (err) {
                        fprintf(stderr, "ThrowbackSendStart %s", err->errmess);
                        ThrowbackStarted = -1;
                } else
                        ThrowbackStarted = 1;
        }
        if (ThrowbackStarted > 0) {
                ThrowbackSendError(level, lineno, error);
        }
}

os_error *ThrowbackStart(void)
{
          _kernel_swi_regs regs;
          return (os_error *) _kernel_swi(os_X | DDEUtils_ThrowbackStart, &regs, &regs);
}


os_error *ThrowbackSendStart(char *filename)
{
        _kernel_swi_regs regs;
        ErrorFile = filename;

        regs.r[0] = Throwback_ReasonProcessing;
        regs.r[1] = 0;
        regs.r[2] = (int) filename;
        return (os_error *) _kernel_swi(os_X |  DDEUtils_ThrowbackSend,  &regs, &regs);
}

os_error *ThrowbackSendError(int level, int lineno, char *error)
{
        _kernel_swi_regs regs;
        regs.r[1] = 0;
        regs.r[2] = (int) ErrorFile;
        regs.r[3] = lineno;
        regs.r[5] = (int) error;
        if (level == ThrowbackInfo) {
          regs.r[0] = Throwback_ReasonInfoDetails;
          regs.r[4] = 0;
        } else {
          regs.r[0] = Throwback_ReasonErrorDetails;
          regs.r[4] = level;
        }
        return (os_error *) _kernel_swi(os_X |  DDEUtils_ThrowbackSend,  &regs, &regs);

}

os_error *ThrowbackEnd(void)
{
        _kernel_swi_regs regs;

        return (os_error *) _kernel_swi(os_X | DDEUtils_ThrowbackEnd, &regs, &regs);

}
/* throwback handling - end */

/* readline lex_yy handling - begin */
char *kgets(char *s, int n)
{
        _kernel_swi_regs        r;

        if (n > 0)
        {
                r.r[0] = (int) s;
                r.r[1] = n-1;
                r.r[2] = 32;
                r.r[3] = 255;
                _kernel_swi(OS_ReadLine, &r, &r);
                s[r.r[1]] = '\n';
                s[r.r[1]+sizeof(char)] = '\0';
        }
        return s;
}
/* readline lex_yy handling - end */


/* filetype checking - begin*/
void
set_rlab_filetype (value)
     char *value;
{
  rlab_filetype = cpstr (value);
}

int
ro_check_type (char *name)
{
  int ret_var;
  _kernel_swi_regs regs;


  regs.r[0] = 5; /*check file existance*/
  regs.r[1] = (int) name;
  _kernel_swi(OS_File, &regs, &regs);
  if (regs.r[0] == 1) {                     /*file exists*/
    regs.r[0] = 31;   			/*convert filetype text to number*/
    regs.r[1] = (int) rlab_filetype;
    _kernel_swi(OS_FSControl, &regs, &regs);
    if (filetype(name) == regs.r[2]) {
        ret_var = CORRECT_TYPE;
    } else {
        ret_var = WRONG_TYPE;
    }
  } else {
        ret_var = NOT_EXISTING;
  }
  return(ret_var);
}

int filetype (const char *file)
{
	int type;
	_kernel_osfile_block osfile;

	type = _kernel_osfile(ReadCat, file, &osfile);

	if (type < 0)
		return F_ERROR;

	else if (type == 0)
		return F_NONE;

	else if (type == 2)
		return F_DIR;

	else if ((osfile.load & 0xFFF00000) != 0xFFF00000)
		return F_UNSTAMPED;

	type = (osfile.load >> 8) & 0xFFF;
	return type;
}

/* filetype checking - end*/


void
ro_output_warn(int interactive, char* name, int type_eval)
{
    if (interactive != 0) {
    /* interactive mode is ON */
        if (type_eval == WRONG_TYPE) {
  	   warning_1 ("has wrong filetype, cannot load");
        } else {
          if (type_eval == NOT_EXISTING) {
             warning_1 ("does not exist, cannot load");
          } else {
             /* this should never happen */
             warning_1 ("cannot load, unknown cause");
          }
        }
    } else {
    /* interactive mode is OFF */
        if (type_eval == WRONG_TYPE) {
  	   printf ("%s has wrong filetype, cannot load\n", name);
        } else {
          if (type_eval == NOT_EXISTING) {
             printf ("%s does not exist, cannot load\n", name);
          } else {
             /* this should never happen */
             printf ("cannot load %s, unknown cause\n", name);
          }
        }
    }
}

/* pipe from pipe file - begin */
void
set_rlab_pipe_in (value)
     char *value;
{
  rlab_pipe_in = cpstr (value);
}

char *ro_gets(char *s, int n)
{
  FILE *pipe;

  if (pipe_on) {
    do {
      pipe=fopen(rlab_pipe_in, "r");
    } while (pipe == NULL);

    fgets(s, n, pipe);
    fclose (pipe);
    if (!strcmp(s,"pipe_off\n")) {
       pipe_on = 0;
       strcpy(s,"\n");
    }
  }
  else {
    kgets(s, n);
    if (!strcmp(s,"pipe_on\n")) {
       pipe_on = 1;
       strcpy(s,"\n");
    }
  }
  return s;
}

/* pipe from pipe file - end */
/* chdir - begin */

int chdir (const char *dir)
{
        _kernel_swi_regs regs;

        regs.r[0] = SetDir;
        regs.r[1] = (int)dir;

        if ( _kernel_swi(OS_FSControl,&regs,&regs) )
                return -1;
        else
                return 0;
}

/* chdir - end */
/* C.Dir: Directory handling */

DIR *opendir (char *name)
{
	register char *cache;
	register DIR *dirp;
	_kernel_swi_regs regs;

	dirp = malloc(sizeof(DIR));
	if (dirp == 0)
		return NULL;

	cache = malloc(CACHE_SIZE);
	if (cache == NULL)
	{
		free(dirp);
		return NULL;
	}

	regs.r[0] = ReadDir;
	regs.r[1] = (int) name;
	regs.r[2] = (int) cache;
	regs.r[3] = CACHE_SIZE;
	regs.r[4] = 0;
	regs.r[5] = CACHE_SIZE;
	regs.r[6] = (int) "*";

	if (_kernel_swi(OS_GBPB, &regs, &regs))
	{
		free (cache);
		free (dirp);
		return NULL;
	}

	dirp->dd_pos = 0;
	dirp->dd_num = regs.r[3];
	dirp->dd_loc = cache;
	dirp->dd_cache = cache;

	return dirp;
}

struct direct *readdir (DIR *dirp)
{
	register char *from;
	register char *to;
	static struct direct dir;

	/* If no more entries, return NULL */
	if (dirp->dd_pos >= dirp->dd_num)
		return NULL;

	/* Copy the name into the struct direct */
	to = dir.d_name;
	from = dirp->dd_loc;
	while ((*to++ = *from++) != 0)
		;

	/* Update the DIR structure */
	++dirp->dd_pos;
	dirp->dd_loc = from;

	/* Set up the rest of the struct direct */
	dir.d_ino = 0;
	dir.d_reclen = dir.d_namlen = strlen(dir.d_name);

	return &dir;
}

void closedir (DIR *dirp)
{
	if (dirp != NULL)
	{
		if (dirp->dd_cache != NULL)
			free(dirp->dd_cache);

		free(dirp);
	}
}

int seekdir (DIR *dirp, int pos)
{
	register int old_pos = dirp->dd_pos;
	register int i = 0;
	register char *p = dirp->dd_cache;

	while (i < pos)
	{
		while (*p++ != 0)
			;
		++i;
	}

	dirp->dd_pos = pos;
	dirp->dd_loc = p;

	return old_pos;
}

#endif
