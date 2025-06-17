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
#ifndef RLAB_RISCOS_UT_H
#define RLAB_RISCOS_UT_H

#ifdef __riscos


/* throwback handling - begin */
#define ThrowbackInfo        -1
#define ThrowbackWarning      0
#define ThrowbackError        1
#define ThrowbackSeriousError 2
#define DDEUtils_ThrowbackRegister   0x42585
#define DDEUtils_ThrowbackUnRegister 0x42586
#define DDEUtils_ThrowbackStart      0x42587
#define DDEUtils_ThrowbackSend       0x42588
#define DDEUtils_ThrowbackEnd        0x42580
#define Throwback_ReasonProcessing     0
#define Throwback_ReasonErrorDetails   1
#define Throwback_ReasonInfoDetails    2
#define os_X (0x00020000)

typedef struct os_error os_error;

struct os_error
   {  unsigned int errnum;
      char errmess [252];
   };

os_error *ThrowbackStart _PROTO ((void));
os_error *ThrowbackSendStart _PROTO ((char *filename));
os_error *ThrowbackSendError _PROTO ((int level,int lineno,char *error));
os_error *ThrowbackEnd _PROTO ((void));
void throwback _PROTO ((int, int, char*));
void throwback_finish _PROTO ((void));
void throwback_init _PROTO ((int, char *));
/* throwback handling - end */

/* readline lex_yy handling - begin */
char *kgets _PROTO((char *s, int n));
int keyb_check _PROTO((int code));
/* readline lex_yy handling - end */

/* filetype checking - begin*/
#define RLAB_FILE_TYPE 0xFFF
#define CORRECT_TYPE  1
#define WRONG_TYPE 0
#define NOT_EXISTING -1
#define ReadCat	5
#define F_NONE		(-1)
#define F_DIR		(-2)
#define F_UNSTAMPED	(-3)
#define F_ERROR		(-4)

int ro_check_type   _PROTO((char *));
void set_rlab_filetype _PROTO((char *));
void ro_output_warn _PROTO((int, char*, int));
int filetype _PROTO((const char *file));
/* filetype checking - end*/

/* pipe from pipe file - begin */
void set_rlab_pipe_in _PROTO((char *));
char *ro_gets _PROTO((char *s, int n));
/* pipe from pipe file - end */

/* H.Dir: Directory handling */

#define SetDir		0
#define ReadDir 	9
#define MAX_DIR		255
#define CACHE_SIZE	(MAX_DIR * (MAXNAMELEN + 1))
#define MAXNAMELEN      10      /* Name must be no longer than this */

/*these should never happen with Acorn's C - they are here
  just not to introduce another ifdef in ian's code */
#define ENOENT		1
/* No process matches the specified process ID.  */
#define EACCES		11
/* Bad address. An invalid pointer was detected.  */
#define ENOTDIR 	18
/* File was a directory when something else was wanted.  */


#define MAXNAMELEN      10      /* Name must be no longer than this */

struct direct
{
        long    d_ino;                  /* inode number of entry */
        short   d_reclen;               /* length of this record */
        short   d_namlen;               /* length of d_name string */
        char    d_name[MAXNAMELEN + 1]; /* directory name */
};

#define DIRSIZ(dp) \
                ((sizeof (struct direct) - (MAXNAMELEN+1)) \
                + (((dp)->d_namlen+1 + 3) & ~3))

typedef struct
{
	int	dd_pos;			/* current directory entry */
	int	dd_num;			/* number of directory entries */
        char	*dd_loc;		/* current position in the cache */
        char	*dd_cache;		/* cache of directory names */
}
DIR;

int chdir _PROTO((const char *));
DIR *opendir _PROTO((char *name));
struct direct *readdir _PROTO((DIR *dirp));
void closedir _PROTO((DIR *dirp));
int seekdir _PROTO((DIR *dirp, int pos));


#define telldir(dirp) \
	((dirp)->dd_pos)

#define rewinddir(dirp)	\
	(void)((dirp)->dd_pos = 0, (dirp)->dd_loc = (dirp)->dd_cache)

#endif

#endif /*RLAB_RISCOS_UT_H*/
