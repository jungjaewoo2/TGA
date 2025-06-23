// rfileio.c

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2007-2009  Marijan Kostrun

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
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "bltin.h"
#include "listnode.h"
#include "mem.h"
#include "list.h"
#include "util.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"
#include "msr.h"
#include "msc.h"
#include "mathl.h"
#include "btree.h"
#include "mdrf1.h"

#include "getline.h"

#include "rlabplus_termcap.h"
#include "rlab_solver_parameters_names.h"

FILE *RLAB_STDERR_DS=0;


#include <math.h>
#include <stdio.h>

#ifdef HAVE_HDF5_SO
static int RLAB_HDF5_DEBUG = 0;
#define HDF5_RLAB_TRANSPOSE_NAME "Transpose"
#define HDF5_RLAB_TRANSPOSE_VAL  1
#endif

#ifdef HAVE_LIBCURL
static int CURLOPT_RLABPLUS_DEBUG = 0;
static int CURL_MAX_NUM_BYTES = 32768;
FILE *curl_stderr_fd=0;
#endif


#ifndef FOPEN_MAX
# define FOPEN_MAX 20
#endif

int RLABPLUS_SERIAL_DEBUG;  // 1, for hex; 2, for char; 3, for int. Set in print.c

#define SERIAL_MAX_NUM_BYTES 32768
int RLABPLUS_SERIAL_2BYTETIME_US  = 0;
int RLABPLUS_SERIAL_TIMEOUT_US    = 100000;

static char
    buffer[SERIAL_MAX_NUM_BYTES];

char * RLABPLUS_FILE_STDOUT = "stdout";
char * RLABPLUS_FILE_STDERR = "stderr";
char * RLABPLUS_FILE_STDIN  = "stdin";
char * RLABPLUS_FILE_EOL_UNIX = "\n";
char * RLABPLUS_FILE_CSP = " ";
char * RLABPLUS_FILE_FMT = "?";
char * RLABPLUS_FILE_NAN = "nan";
char * RLABPLUS_FILE_INF_POS = "inf";
char * RLABPLUS_FILE_INF_NEG = "-inf";

static int SOCKET_MAX_NUM_BYTES = 32768;

extern int rpclose (FILE * fp); // close pipe safely, see main.c

FILE *get_file_ds (char *name, char *mode, int buffsize);
int   close_file_ds (char *name);

void *matrix_ReadB (FILE * fn, int T, int P, char **name, int *rtype);
void *msparse_ReadB (FILE * fn, int T, int P, char **name, int *rtype);

void btree_WriteB (Rfile *rf, char *name, Btree * btree);
void mdr_WriteB (Rfile *rf, char *name, MDR * m);
void mdc_WriteB (Rfile *rf, char *name, MDC * m);
void mds_WriteB (Rfile *rf, char *name, MDS * m);
void msr_WriteB (Rfile *rf, char *name, MSR * m);
void msc_WriteB (Rfile *rf, char *name, MSC * m);

// void btree_WriteB (Btree * btree, FILE * fn, char *name);
// void mdr_WriteB (MDR * m, FILE * fn, char *name);
// void mdc_WriteB (MDC * m, FILE * fn, char *name);
// void mds_WriteB (MDS * m, FILE * fn, char *name);
// void msr_WriteB (MSR * m, FILE * fn, char *name);
// void msc_WriteB (MSC * m, FILE * fn, char *name);

Btree *btree_ReadB (FILE * fn, char **name);

static int ferr_check (char *name, FILE * fn);
static int fread_check (FILE * fn);
static void reverse_word (void *word_ptr, int nbytes);

static void mdr_WriteASCII (MDR * m, FILE * fn, char *name);
static void mdc_WriteASCII (MDC * m, FILE * fn, char *name);
static void mds_WriteASCII (MDS * m, FILE * fn, char *name);
static void btree_WriteASCII (Btree * btree, FILE * fn, char *name);

/*
 * Interface arrays for binary read/write operations.
 */

static OpDef writeb_method[NCL];
static OpDef writem_method[NCL];
static OpDef write_method[NCL];

#include "rfileio_write.c"

void class_io_init (void)
{
  /*
   * writeb ()
   */

  writeb_method[BTREE].type = 0;
  writeb_method[BTREE].op = (void *) btree_WriteB;

  writeb_method[MATRIX_DENSE_REAL].type = 0;
  writeb_method[MATRIX_DENSE_REAL].op = (void *) mdr_WriteB;

  writeb_method[MATRIX_DENSE_COMPLEX].type = 0;
  writeb_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_WriteB;

  writeb_method[MATRIX_SPARSE_REAL].type = 0;
  writeb_method[MATRIX_SPARSE_REAL].op = (void *) msr_WriteB;

  writeb_method[MATRIX_SPARSE_COMPLEX].type = 0;
  writeb_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_WriteB;

  writeb_method[MATRIX_DENSE_STRING].type = 0;
  writeb_method[MATRIX_DENSE_STRING].op = (void *) mds_WriteB;

  /*
   * writem ()
   */

  writem_method[MATRIX_DENSE_REAL].type = 0;
  writem_method[MATRIX_DENSE_REAL].op = (void *) mdr_WriteGeneric;

  writem_method[MATRIX_DENSE_COMPLEX].type = 0;
  writem_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_WriteGeneric;

  writem_method[MATRIX_SPARSE_REAL].type = 0;
  writem_method[MATRIX_SPARSE_REAL].op = (void *) msr_WriteGeneric;

  writem_method[MATRIX_SPARSE_COMPLEX].type = 0;
  writem_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_WriteGeneric;

  writem_method[MATRIX_DENSE_STRING].type = 0;
  writem_method[MATRIX_DENSE_STRING].op = (void *) mds_WriteGeneric;

  /*
   * write ()
   */

  write_method[BTREE].type = 0;
  write_method[BTREE].op = (void *) btree_WriteASCII;

  write_method[MATRIX_DENSE_REAL].type = 0;
  write_method[MATRIX_DENSE_REAL].op = (void *) mdr_WriteASCII;

  write_method[MATRIX_DENSE_COMPLEX].type = 0;
  write_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_WriteASCII;

  write_method[MATRIX_DENSE_STRING].type = 0;
  write_method[MATRIX_DENSE_STRING].op = (void *) mds_WriteASCII;

}


static Rfile *rfile_list = 0;	  // Ptr to 1st Rfile
static int nrfile = 0;		      // Number of open files

#ifdef HAVE_HDF5_SO
#include "rfileio_hdf5.c"
#endif


#include "serial.c"
extern int RLABPLUS_SERIAL_DEBUG; // see rfileio.c, writes information about serial port

//
// Create a new Rfile struct, and tack it on
// the front of the list. Return a ptr to the
// new Rfile struct. Make sure and bump the
// list element count.
//
Rfile *
    rfile_Create (void)
{
  Rfile *new = (Rfile *) GC_MALLOC (sizeof (Rfile));
  if (new == 0)
    rerror ("out of memory");

  // zero the file structure
  new->name   = 0;
  new->mode   = 0;
  new->filetype = RFILE_NULL;   // class of file
  new->fileds_f =  0;            // for files
  new->fileds_i = -1;            // for serial ports and sockets
  new->buffer = 0;
  new->eol = 0;
  new->csp = 0;
  new->fmt = 0;
  new->buffer_len = 0;
  new->buffer_pos = 0;

#ifdef HAVE_LIBCURL
  // zero the libcurl objects:
  new->curl = 0;
  new->curl_stderr_fd = 0;
#endif

#ifdef HAVE_HDF5_SO
  // zero the Hierarhical Data Format Ver. 5 objects
  new->h5_file     = -1;
  new->h5_flags    =  0;
  new->h5_fapl_id  =  0;
  new->h5_status   =  0;
#endif

  // zero the socket structure
  memset(&new->remote_server, 0, sizeof(new->remote_server));
  new->socket_timeout = 0;

  // finish
  new->next = rfile_list;
  nrfile++;
  rfile_list = new;
  return (new);
}

int
rfile_Destroy (char *name)
{
  Rfile *rf, *next, *prev;

  rf = rfile_list;
  next = rf->next;
  prev = rf;

  while (next)
  {
//     fprintf(stdout, "rf->name = %s : name = %s\n",rf->name, name );
    if (!strcmp (rf->name, name))
    {
      // Found it!

      // Re-Hook up list
      prev->next = next;

      // Check two special cases
      if (prev == rf)     // 1st element
        rfile_list = next;
      else if (next == 0) // Last element
        prev->next = 0;

      // Now destroy the element
      switch (rf->filetype)
      {
        case RFILE_FILE:
          if (rf->name[0]=='|')
          {
            if (rf->fileds_f)
              rpclose (rf->fileds_f);
          }
          else
          {
            if (rf->fileds_f)
              fclose (rf->fileds_f);
          }
          break;

        case RFILE_COMM:
          serial_Close ( rf );
          break;

#ifdef HAVE_LIBCURL
        case RFILE_CURL:
          curl_easy_cleanup(rf->curl);
          if(rf->curl_stderr_fd)
            fclose(rf->curl_stderr_fd);
          break;
#endif
#ifdef HAVE_HDF5_SO
        case RFILE_H5:
        {
          // flush and close the file
          H5Fflush (rf->h5_file, H5F_SCOPE_LOCAL);
          H5Fclose (rf->h5_file);
          break;
        }
#endif
        case RFILE_SOCKET:
        {
          if (rf->fileds_i > -1)
            close(rf->fileds_i);
          break;
        }

        default:
          break;
      }

      if (rf->name)
      {
        GC_FREE (rf->name);
        rf->name = 0;
      }
      if (rf->mode)
      {
        GC_free (rf->mode);
        rf->mode = 0;
      }
      if (rf->buffer)
      {
        GC_free (rf->buffer);
        rf->buffer = 0;
      }
      if (rf->eol)
      {
        GC_free (rf->eol);
        rf->eol = 0;
      }
      if (rf->nan)
      {
        GC_free (rf->nan);
        rf->nan= 0;
      }
      if (rf->inf_pos)
      {
        GC_free (rf->inf_pos);
        rf->inf_pos = 0;
      }
      if (rf->inf_neg)
      {
        GC_free (rf->inf_neg);
        rf->inf_neg = 0;
      }
      if (rf->csp)
      {
        GC_free (rf->csp);
        rf->csp = 0;
      }
      if (rf->fmt)
      {
        mds_Destroy (rf->fmt);
        rf->fmt = 0;
      }

      rf->filetype = RFILE_NULL;
      rf->next = 0;

      GC_FREE (rf);
      rf = 0;

      nrfile--;
      return (1);   /* Success: an rfile is destroyed */
    }

    prev = rf;
    rf = next;
    next = rf->next;
  }

  return (0);     /* Failure: no rfile with such name was found, so none was destroyed */
}

//
// Initialize the file list
//
void
    init_file_list (void)
{
  rfile_list = rfile_Create ();
}

//
// Walk the list, destroying each node
//
void
    destroy_file_list (void)
{
  Rfile *rf, *next;

  rf = rfile_list;
  next = rf->next;

  while (next)
  {
    rfile_Destroy (rf->name);
    rf = next;
    next = rf->next;
  }
}

//
// Walk the list, look for a particular node
//
Rfile *
rfile_find (char *name)
{
  Rfile *rf, *next;

  rf = rfile_list;
  next = rf->next;

  // if 'name' is a pipe increase it by one - otherwise
  // the pipes are not closed correctly
//   if (name[0]=='|')
//     name++;

  while (next)
  {
    if (!strcmp (rf->name, name))
    { return (rf); }
    rf = next;
    next = rf->next;
  }

  return (0);
}


//
// get rfile pointer with name and check if it is of proper type.
// if it does not exist, or the filetype is RFILE_NULL create a
// new object of desired type
//
#include <stdarg.h>
Rfile *
get_rfile_ds (char * name,  enum RFILE_TYPE rftype, ...)
{
  Rfile *rf=0;
  va_list ap;

  // does rfile with 'name' exists that is of desired type?
  rf = rfile_find (name);
  if (rf)
  {
    if (rf->filetype == rftype)
      return (rf);
    else
    {
      rfile_Destroy (name);
      return 0;
    }
  }

  // Make sure we don't try and open more files than the system will allow
//   fprintf(stderr, "nrfile = %i, %s\n", nrfile, name);

  if (nrfile >= FOPEN_MAX - 3)
    rerror ("get_rfile_ds: Cannot open a new file. Too many files already opened");

  //
  // create new rfile and do all the necessary actions to be of specified type
  //
  rf = rfile_Create();
  rf->name = cpstr(name);
  rf->filetype = rftype;
  rf->buffer = 0;

  if (rftype == RFILE_FILE)
  {
    char *real_name=0;
    //
    // operations with files and pipes:
    //    rf = get_rfile_ds(fname, RFILE_FILE, mode, buffersize)
    //
    va_start(ap, rftype);
    char *mode     = va_arg(ap, char*);   // argument No. 3
    int   buffsize = va_arg(ap, int);     // argument No. 4
    va_end(ap);

    // ignore 'file://' in front of the file name
    if (!strncmp(name, "file://", 7))
      real_name = &name[7];
    else
      real_name = name;

    // is file name there?
    if (!strlen(real_name))
    {
      rfile_Destroy(name);
      fprintf(stderr, "get_rfile_ds: Cannot have zero-length filename\n");
      return 0;
    }

    rf->mode = cpstr (mode);

    if (real_name[0] == '|')
    {
      // '|' indicates pipe
      rf->fileds_f = popen (&real_name[1], rf->mode);
      if (!rf->fileds_f)
      {
        rfile_Destroy (name);
        return 0;
      }
    }
    else
    {
      // we are trying to open a file
      rf->fileds_f = fopen (real_name, rf->mode);
      if (!rf->fileds_f)
      {
        rfile_Destroy (name);
        return 0;
      }
      if (buffsize != 0)
      {
        rf->buffer = (char *) GC_MALLOC (sizeof (char) * buffsize);
        if (rf->buffer == 0)
          rerror ("get_rfile_ds: Cannot create buffer of desired size. Out of memory");
        if (setvbuf (rf->fileds_f, rf->buffer, _IOFBF, buffsize))
        {
          fprintf (stderr, "get_rfile_ds: Cannot create buffer of desired size. Out of memory");
          rfile_Destroy (name);
          return (0);
        }
      }
    }
    return rf;
  } // RFILE_FILE

#ifdef HAVE_HDF5_SO
  if(rftype == RFILE_H5)
  {
    char *real_name=0;
    //
    // operations with hierarhical data format ver. 5 files:
    //    rf = get_rfile_ds (name, RFILE_H5, mode);
    //
    va_start(ap, rftype);
    char *mode = va_arg(ap, char*);   // argument No. 3
    va_end(ap);

    // ignore protocol 'h5://' or 'hdf5://' in front of the file name
    if (!strncmp(name, "h5://", 5))
      real_name = &name[5];
    else if (!strncmp(name, "hdf5://", 7))
      real_name = &name[7];
    else
      real_name = name;

    // is file name there?
    if (!strlen(real_name))
    {
      rfile_Destroy(name);
      fprintf(stderr, "get_rfile_ds: Cannot have zero-length filename!\n");
      return 0;
    }

    rf->mode = cpstr (mode);

    switch (mode[0])
    {
      case 'r':
      case 'R':
        // user want to read from the file. does the file exist?
        if(H5Fis_hdf5(real_name) < 1)
        {
          // file either exists but is not HDF5, or it does not exist
          rfile_Destroy(name);
          fprintf(stderr, "get_rfile_ds: Attempting to read from non-HDF5 file\n");
          return 0;
        }

        // open file for reading
        rf->h5_flags = H5F_ACC_RDONLY;
        rf->h5_file  = H5Fopen(real_name, rf->h5_flags, H5P_DEFAULT);
        break;

      case 'a':
      case 'A':
        // open file for read/write operations
        rf->h5_flags = H5F_ACC_RDWR;
        rf->h5_file  = H5Fopen(real_name, rf->h5_flags, H5P_DEFAULT);
        break;

      case 'w':
      case 'W':
        // create new HDF5 file
        rf->h5_flags = H5F_ACC_TRUNC;
        rf->h5_file = H5Fcreate (real_name, rf->h5_flags, H5P_DEFAULT, H5P_DEFAULT);
        break;
    }

    // get back to user
    return rf;
  }
#endif
  else if (rftype == RFILE_COMM)
  {
    //
    // operations with serial port
    //    rf = get_rfile_ds(name, RFILE_COMM, ispeed, idatb, iparity, istop, iflow, iraw);
    //
    va_start(ap, rftype);
    int ispeed  = va_arg(ap, int);     // argument No. 3
    int idatb   = va_arg(ap, int);     // argument No. 4
    int iparity = va_arg(ap, int);     // argument No. 5
    int istop   = va_arg(ap, int);     // argument No. 6
    int iflow   = va_arg(ap, int);     // argument No. 7
    int iraw    = va_arg(ap, int);     // argument No. 8
    int hupcl   = va_arg(ap, int);     // argument No. 9
    va_end(ap);
//     int i;

    // this is how long in usec it takes to transfer two bytes at the current speed:
    RLABPLUS_SERIAL_2BYTETIME_US = (int) 16e6 / (ispeed + 0.0);

    // print debugging information
    if (RLABPLUS_SERIAL_DEBUG)
    {
      fprintf(stdout, "Configuration for Serial port %s\n", name);
      fprintf(stdout, "  speed       : %i\n", ispeed);

      //
      fprintf(stdout, "  data format : %1i", idatb);
      if (iparity == 0)
        fprintf(stdout, "N");
      else if (iparity == 1)
        fprintf(stdout, "O");
      else if (iparity == 2)
        fprintf(stdout, "E");
      fprintf(stdout, "%1i\n", istop);

      //
      fprintf(stdout, "  flow control: ");
      if (iflow == 2)
        fprintf(stdout, "hardware");
      else if (iflow == 1)
        fprintf(stdout, "xon|xoff");
      else
        fprintf(stdout, "none");
      fprintf(stdout, "\n");

      //
      fprintf(stdout, "  raw         : ");
      if (iraw)
        fprintf(stdout, "yes\n");
      else
        fprintf(stdout, "canonical\n");
      if (hupcl)
        fprintf(stdout, "  hangup on close\n");
    }
    rf->buffer = 0;

    if(serial_OpenDevice (rf, ispeed, idatb, iparity, istop, iflow, iraw, hupcl))
    {
      fprintf (stderr, "get_rfile_ds: Terrible internal error. Serial port failed to open\n");
      rfile_Destroy (name);
      return (0);
    }

    return rf;
  }
#ifdef HAVE_LIBCURL
  else if (rftype == RFILE_CURL)
  {
    // libcurl easy interface magic
    rf->curl = curl_easy_init();
    if(!rf->curl)
    {
      fprintf (stderr, "libcurl: Terrible internal error. 'curl_easy_init' failed.\n");
      rfile_Destroy (name);
      return (0);
    }

    // do a default initialization. user can rewrite by specifying CURLOPT_URL argument
    // later
    curl_easy_setopt(rf->curl, CURLOPT_URL, (name));
    return rf;
  }
#endif
  else if (rftype == RFILE_SOCKET)
  {
    // get_rfile_ds (name, RFILE_SOCKET, int );
    va_start(ap, rftype);
    int opt_proto  = va_arg(ap, int);     // argument No. 3
    int opt_mode   = va_arg(ap, int);     // argument No. 4
    double opt_dto = va_arg(ap, double);  // argument No. 5
    va_end(ap);

    // figure out the timeout
    rf->socket_timeout = opt_dto;

    if (opt_proto == 1)
      rf->fileds_i = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
    else
      rerror ("get_rfile_ds: UDP packets not yet supported");

//     fprintf(stderr, "socket handle = %i\n", rf->fileds_i);

    if (rf->fileds_i < 0)
    {
      rfile_Destroy(name);
      fprintf(stderr, "get_rfile_ds: failed to create socket");
      return 0;
    }

    if (opt_mode==1)
    {
      // connect:
      // process the name:
      //    tcp://{host}:port
      //    udp://{host}:port
      // if 'host' is omitted assume localhost (127.0.0.1)
      // if 'port' for telnet is omitted assume 25
      char *hostport = cpstr( &name[6] );
      char *host = 0;
      char lastchar = '\0';
      char *localhost = "127.0.0.1";
      char *port_string = 0;
      int i=0;

      // did user provide host name?
      if (*hostport == ':')
        host = localhost;
      else
      {
        while(hostport[i] && hostport[i]!=':')
        { i++; }
        lastchar = hostport[i];
        hostport[i] = '\0';
        host = hostport;
      }
      memset(&rf->remote_server, 0, sizeof(rf->remote_server));
      rf->remote_server.sin_family = AF_INET;
      rf->remote_server.sin_addr.s_addr = inet_addr( host );
      if (!rf->remote_server.sin_addr.s_addr)
      {
        struct hostent *myhost = gethostbyname( host );
        struct in_addr h_addr;
        if (!myhost)
        {
          fprintf(stderr, "get_rfile_ds: cannot resolve host name \"%s\"\n", host);
          rfile_Destroy (name);
          if (hostport)
            GC_free (hostport);
          if (myhost)
            free (myhost);
          return 0;
        }

        h_addr.s_addr = *((unsigned long *) myhost->h_addr_list[0]);
        rf->remote_server.sin_addr.s_addr = inet_addr( inet_ntoa(h_addr) );

        if (myhost)
          free (myhost);
      }

      //     fprintf(stderr, "host = %s\n", host);
      if (lastchar)
        hostport[i] = lastchar;

      // did user provide port number(s)?
      if (hostport[i] == ':')
      {
        port_string = &hostport[i+1];
        //       fprintf(stderr, "port = %s\n", port_string);
        if (*port_string)
          rf->remote_server.sin_port = htons( atoi( port_string ) );
        else
        {
          fprintf(stderr, "get_rfile_ds: cannot resolve port range from \"%s\"\n", name);
          rfile_Destroy (name);
          if (hostport)
            GC_free (hostport);
          return 0;
        }
      }
      else
      {
        fprintf(stderr, "get_rfile_ds: cannot resolve port range from \"%s\"\n", name);
        rfile_Destroy (name);
        if (hostport)
          GC_free (hostport);
        return 0;
      }

      if (connect(rf->fileds_i,
        (struct sockaddr *) &rf->remote_server,sizeof(rf->remote_server)) < 0)
      {
        rfile_Destroy (name);
        if (hostport)
          GC_free (hostport);
        fprintf(stderr, "get_rfile_ds: failed to connect to remote server\n");
        return 0;
      }
      if (hostport)
        GC_free (hostport);
    }
    else if (opt_mode == 0)
    {
      // listen

    }

    return rf;
  }

  return (0);
}

//
// Get the file/pipe descriptor asociated with the char string. This function does
// all the work, searches the list, uses popen() instead of fopen when necessary.
// It is used by other file manipulating codes that do not need full access to Rfiles.
//
FILE *
get_file_ds (char *name, char *mode, int buffsize)
{
  int Pipe;
  FILE *fileds_f=0;
  Rfile *rf=0, *new=0;

  //
  // Check for stdout, stderr. We don't want to put them
  // on the file_list, since we don't want to close stdout
  // or stderr ever
  //
  if (!strcmp ("stdout", name))
    return (stdout);
  else if (!strcmp ("stderr", name))
    return (stderr);
  else if (!strcmp ("stdin", name))
    return (stdin);

  /* Check for pipe symbol '|' */
  if (*name == '|')
  {
    name++;
    Pipe = 1;
  }
  else
    Pipe = 0;

  // Check list for previously opened occurence of file
  // in particular return FILE pointer if such exists
  if ((rf = rfile_find (name)))
  {
    if (rf->filetype == RFILE_FILE)
      return (rf->fileds_f);
    else
      return 0;
  }
  else
  {
    /*
     * Make sure we don't try and open more files
     * than the system will allow
     */

    if (nrfile >= FOPEN_MAX - 3)
      rerror ("exceeded system limit for # of open files");

#ifdef HAVE_PIPE
    if (Pipe)
    {
      if ((fileds_f = popen (name, mode)) != 0)
      {
        // Put new fileds_f on list
        new = rfile_Create ();
        new->name     = cpstr (name);
        new->fileds_f = fileds_f;
        new->filetype = RFILE_FILE;
        new->mode     = cpstr (mode);
      }
      else
        return (0);
    }
    else
    {
#endif
      if ((fileds_f = fopen (name, mode)) != 0)
      {
        // Put new fileds_f on list
        new = rfile_Create ();
        new->name = cpstr (name);
        new->fileds_f = fileds_f;
        new->mode = cpstr (mode);
        new->filetype = RFILE_FILE;
        if (buffsize != 0)
        {
          new->buffer = (char *) GC_MALLOC (sizeof (char) * buffsize);
          if (new->buffer == 0)
            rerror ("out of memory");
          if (setvbuf (fileds_f, new->buffer, _IOFBF, buffsize))
          {
            fprintf (stderr, "open: cannot create I/O buffer\n");
            rfile_Destroy (name);
            return (0);
          }
        }
        else
        {
          new->buffer = 0;
        }
      }
      else
        return (0);
#ifdef HAVE_PIPE
  }
#endif
  return (fileds_f);
  }
}

int get_int_file_ds (char *name)
{
  Rfile *rf;

  // Check list for previously opened occurence of file
  // in particular return FILE pointer if such exists
  if ((rf = rfile_find (name)))
  {
    if (rf->filetype == RFILE_COMM)
      return (rf->fileds_i);
  }

  return -1;
}

char * get_eol_file_name (char *name)
{
  Rfile *rf;

  if ((rf = rfile_find (name)))
    return (rf->eol);

  return 0;
}


/* Close the file asociated with the file name */
int
close_file_ds (char *name)
{
  /* Check list */
  if (name == 0)
    return 0;

//   fprintf(stdout,"close_file_ds: closing %s\n", name);
//   if (*name == '|')
//     name++;

  if (rfile_Destroy (name))
  {
    return 1;     /* Success */
  }
  else
  {
    return 0;     /* Failure */
  }
}

/*
 * Give a file descriptor (fn), look up the name of the file.
 * Return NULL if we cannot find it.
 */

char *
get_file_ds_name (FILE * fn)
{
  Rfile *rf, *next;
  /* Check for system file descriptors first. */
  if (fn == stdout)
    return ( RLABPLUS_FILE_STDOUT );
  else if (fn == stderr)
    return ( RLABPLUS_FILE_STDERR );
  else if (fn == stdin)
    return ( RLABPLUS_FILE_STDIN );
  else
  {
    /* Go through the Rfile list, and check each node manually. */
    rf = rfile_list;
    next = rf->next;
    while (next)
    {
      if (fn == rf->fileds_f)
        return (rf->name);  // Found it!
      rf = next;
      next = rf->next;
    }
  }

  /* Hit the end of the list, and we still haven't found it. */
  return NULL;
}

char *
get_eol_file_ds_name (FILE * fn)
{
  Rfile *rf, *next;
  /* Check for system file descriptors first. */
  if (fn == stdout)
    return RLABPLUS_FILE_EOL_UNIX;
  else if (fn == stderr)
    return RLABPLUS_FILE_EOL_UNIX;
  else if (fn == stdin)
    return RLABPLUS_FILE_EOL_UNIX;
  else
  {
    /* Go through the Rfile list, and check each node manually. */
    rf = rfile_list;
    next = rf->next;
    while (next)
    {
      if (fn == rf->fileds_f)
      {
        if (rf->eol)
          return (rf->eol);
        else
          return RLABPLUS_FILE_EOL_UNIX;
      }
      rf = next;
      next = rf->next;
    }
  }

  /* Hit the end of the list, and we still haven't found it. */
  return RLABPLUS_FILE_EOL_UNIX;
}

char *
get_csp_file_ds_name (FILE * fn)
{
  Rfile *rf, *next;
  /* Check for system file descriptors first. */
  if (fn == stdout)
    return RLABPLUS_FILE_EOL_UNIX;
  else if (fn == stderr)
    return RLABPLUS_FILE_EOL_UNIX;
  else if (fn == stdin)
    return RLABPLUS_FILE_EOL_UNIX;
  else
  {
    /* Go through the Rfile list, and check each node manually. */
    rf = rfile_list;
    next = rf->next;
    while (next)
    {
      if (fn == rf->fileds_f)
      {
        if (rf->csp)
          return (rf->csp);
        else
          return RLABPLUS_FILE_CSP;
      }
      rf = next;
      next = rf->next;
    }
  }

  /* Hit the end of the list, and we still haven't found it. */
  return RLABPLUS_FILE_CSP;
}

MDS *
get_fmt_file_ds_name (FILE * fn)
{
  Rfile *rf, *next;

  // Check for system file descriptors first
  if (fn == stdout || fn == stderr || fn == stdin)
    return mds_CreateScalar( RLABPLUS_FILE_FMT );
  else
  {
    // Go through the Rfile list, and check each node manually.
    rf = rfile_list;
    next = rf->next;
    while (next)
    {
      if (fn == rf->fileds_f)
      {
        if (rf->fmt)
          return (rf->fmt);
        else
          return mds_CreateScalar( RLABPLUS_FILE_FMT );
      }
      rf = next;
      next = rf->next;
    }
  }

  // Hit the end of the list, and we still haven't found it
  return mds_CreateScalar( RLABPLUS_FILE_FMT );
}


/* **************************************************************
 * Read a generic matrix from a file. The format is:
 * All lines:    whitespace separated numbers.
 * ************************************************************** */
#ifdef HAVE_LIBCURL

// transfers data from 'buffer' to the rlab variable 'userp'
static size_t write_callback(char *buffer, size_t size, size_t nitems, void *userp)
{
  char *newbuff;
  int   rembuff;
  Rfile * url = (Rfile *) userp;

  size *= nitems;

  if (!size)
    return 0;

  rembuff = url->buffer_len - url->buffer_pos; /* remaining space in buffer */

  if(size > rembuff)
  {
    /* not enough space in buffer */
    newbuff = GC_realloc(url->buffer, url->buffer_len + (size - rembuff));
    if(newbuff==NULL)
    {
      if (CURLOPT_RLABPLUS_DEBUG)
        fprintf(stderr,"libcurl: Terrible internal error: Callback buffer grow failed\n");
      size = rembuff;
    }
    else
    {
      /* realloc suceeded increase buffer size*/
      url->buffer_len += size - rembuff;
      url->buffer = newbuff;

      if (CURLOPT_RLABPLUS_DEBUG)
        fprintf(stderr, "libcurl: Callback buffer grown to %d bytes\n",url->buffer_len);
    }
  }

  memcpy(&url->buffer[url->buffer_pos], buffer, size);
  url->buffer_pos += size;

  if (CURLOPT_RLABPLUS_DEBUG)
    fprintf(stderr, "libcurl: Callback %d size bytes\n", (int) size);

  return size;
}

#endif

#include "r_curl.c"

//
//
//
#define RLAB_READM_DEFAULT_CSP       "'(BLANK|',')"
#define RLAB_READM_DEFAULT_JOINCSP   " "
#undef THIS_SOLVER
#define THIS_SOLVER "readm"
Ent * ReadM (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  char *fname=0, c;
  int block_size = 0;
  int i, j, rv;
  int close_after_rw = 0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  Rfile *rf=0;

//   MDS *s1=0;
  MDR *m=0, *a2=0;

  // Check n_args
  if (nargs != 1 && nargs != 2 && nargs != 3)
    rerror ("readm: one, two or three arguments required!");

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror ("readm: first argument must be a single string");
  fname = class_char_pointer( e1);
  if (!fname)
    rerror ("readm: empty filename");
  if (!strlen(fname))
    rerror ("readm: empty filename");

  // Check if file exists in our list.
  // If the file does not exist pursue default behavior:
  //    User is trying to read an existing file on local file system
  //    However, after reading it, remove it from our list of open files.
  rf = rfile_find(fname);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after reading
    if (*fname == '|')
    {
      // just assume it is a pipe
      rf = get_rfile_ds (fname, RFILE_FILE, "r", 0);
    }
    else if (    !strncmp(fname, "h5://", 5)
             ||  !strncmp(fname, "hdf5://", 7)
             ||  (strstr (fname, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (fname, RFILE_H5, "r" );
    }
    else if (!strncmp(fname, "tcp://", 6))
    {
      // do the socket
      rf = get_rfile_ds (fname, RFILE_SOCKET, 1);
    }
    else
      if (    (!strncmp(fname, "http://", 7))
                ||  (!strncmp(fname, "https://", 8))
                ||  (!strncmp(fname, "ftp://", 6))
         )
    {
      // initialize CURL with 'name' if it already does not exist
      rf = get_rfile_ds(fname, RFILE_CURL);
    }
    else
    {
      // just assume it is a text file
      rf = get_rfile_ds (fname, RFILE_FILE, "r", 0);
    }
    close_after_rw = 1;
  }
  if (!rf)
    rerror ("readm: cannot open file for reading");

  if (rf->filetype == RFILE_COMM)
  {
    //
    // reading from serial port
    //
    MDS *ws = 0;
    MDR *wi = 0;
    int  inb=0, j;

//     int i, ist, icn;
//     int  n=1;
//     int lenq=0;
//     MDS *sq = 0;

    if (nargs > 1)
    {
      // second argument is buffer size
      e2 = bltin_get_ent(args[1]);
      if (ent_type(e2) == MATRIX_DENSE_REAL)
        inb = (int) class_double ( e2 );
    }

    rent = ent_Create();

    if (RLABPLUS_SERIAL_DEBUG)
    {
      fprintf(stderr, "RLABPLUS_SERIAL_TIMEOUT_US=%i, RLABPLUS_SERIAL_2BYTETIME_US=%i\n",
              RLABPLUS_SERIAL_TIMEOUT_US, RLABPLUS_SERIAL_2BYTETIME_US);
    }

    if (inb)
    {
      // read character at the time until 'cr' or 'lf' is encountered
      char *bufptr = buffer;
      int bytes=0, bytes_prev=0, still_arriving=0;

      // sit on the serial port for at most 1 second and wait for output
      for (j=0; j<RLABPLUS_SERIAL_TIMEOUT_US; j+=RLABPLUS_SERIAL_2BYTETIME_US)
      {
        // how many bytes there are at the input buffer
        bytes_prev = bytes;
        ioctl(rf->fileds_i, FIONREAD, &bytes);
        if (bytes_prev != bytes)
          still_arriving = 1;
        else
        {
          if(still_arriving)
            still_arriving = 0;
        }

        if (RLABPLUS_SERIAL_DEBUG)
        {
          fprintf(stderr, "bytes=%i, bytes_prev=%i, still_arriving=%i\n", bytes, bytes_prev, still_arriving);
        }

        // more then zero: go out and read them
        if (bytes > 0 && bytes_prev==bytes && still_arriving==0)
          break;

        // sleep a little to see if more characters are coming?
        usleep(RLABPLUS_SERIAL_2BYTETIME_US);
      }

      // user specified size of buffer - this is how many characters are
      // read and returned to user
      inb = inb > bytes ? bytes : inb;
      wi  = mdi_Create(1, inb);

      // copy data from I/O buffer to local memory
      for(j=0; j<bytes && j<inb; j++)
      {
        // read it first
        read(rf->fileds_i, bufptr, 1);

        // debug: write raw data byte
        if (RLABPLUS_SERIAL_DEBUG == 1)
          fprintf(stdout, " 0x%02x", (*bufptr)&0xff);
        else if (RLABPLUS_SERIAL_DEBUG == 2)
          fprintf(stdout, "%c", *bufptr);
        else if (RLABPLUS_SERIAL_DEBUG == 3)
          fprintf(stdout, " %i", *bufptr);


        MdiV0(wi, j) = ((unsigned int) *bufptr) & 0xff;

        // clean-up
        *bufptr = '\0';

        // move up in the buffer
        bufptr++;
      }

      if (RLABPLUS_SERIAL_DEBUG)
        fprintf(stdout, "\n");

      ent_data(rent) = wi;
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else
    {
      // read character at the time until 'cr' or 'lf' is encountered
      char *bufptr = buffer;
      int bytes=0, bytes_prev=0, still_arriving=0;

      // sit on the serial port for at most 1 second and wait for output
      for (j=0; j<RLABPLUS_SERIAL_TIMEOUT_US; j+=RLABPLUS_SERIAL_2BYTETIME_US)
      {
        // how many bytes there are at the input buffer
        bytes_prev = bytes;
        ioctl(rf->fileds_i, FIONREAD, &bytes);
        if (bytes_prev != bytes)
        {
          still_arriving = 10;
        }
        else
        {
          if(still_arriving)
            still_arriving--;
        }

        if (RLABPLUS_SERIAL_DEBUG)
        {
          fprintf(stderr, "bytes=%i, bytes_prev=%i, still_arriving=%i\n", bytes, bytes_prev, still_arriving);
        }

        // more then zero: go out and read them
        if (bytes > 0 && bytes_prev==bytes && still_arriving==0)
          break;

        // sleep a little to see if more characters are coming?
        usleep(10 * RLABPLUS_SERIAL_2BYTETIME_US);
      }
      // we have reached waiting time without getting any output
      if (j>=RLABPLUS_SERIAL_TIMEOUT_US || bytes==0)
        ws = mds_Create(0,0);
      else
      {
        // copy data from I/O buffer to local memory
        for(j=0; j<bytes; j++)
        {
          // read it first
          read(rf->fileds_i, bufptr, 1);

          // debug: write raw data byte
          if (RLABPLUS_SERIAL_DEBUG == 1)
            fprintf(stdout, " 0x%x", *bufptr);
          else if (RLABPLUS_SERIAL_DEBUG == 2)
            fprintf(stdout, "%c", *bufptr);
          else if (RLABPLUS_SERIAL_DEBUG == 3)
            fprintf(stdout, " %i", *bufptr);

          // move up in the buffer
          bufptr++;
        }
        *bufptr = '\0';

        if (RLABPLUS_SERIAL_DEBUG)
          fprintf(stdout, "\n");

        if (bufptr - buffer > 1)
          ws = mds_CreateScalar( buffer );
        else
          ws = mds_Create(0,0);
      }

      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
    } // if (j>=RLABPLUS_SERIAL_TIMEOUT_US || bytes==0)


    ent_Clean(e1);
    ent_Clean(e2);
    ent_Clean(e3);
    ent_Clean(e4);

    return rent;
  }
  else if (rf->filetype == RFILE_SOCKET)
  {

    // second argument is buffer size
    int inb=0, i, res;
    MDS *ws = 0;
    MDR *w = 0;
    char *rs=0;
    struct timeval timeout;  /* Timeout for select */
    fd_set socks;

    if (nargs > 1)
    {
      e2 = bltin_get_ent(args[1]);
      if (ent_type(e2) == MATRIX_DENSE_REAL)
        inb = (int) class_double ( e2 );
      else
        rerror("readm: improper type of second argument");
    }

    // is there a buffer
    if (!rf->buffer)
    {
      // no buffer. create new
      if (inb < SOCKET_MAX_NUM_BYTES)
      {
        rf->buffer = GC_malloc(SOCKET_MAX_NUM_BYTES * sizeof(unsigned char));
        rf->buffer_len = SOCKET_MAX_NUM_BYTES;
      }
      else
      {
        rf->buffer = GC_malloc(inb * sizeof(unsigned char));
        rf->buffer_len = inb;
      }
    }
    else
    {
      // buffer exists adjust its size only if not bi enough
      if (inb > rf->buffer_len)
      {
        // buffer exists, extend it
        GC_realloc(rf->buffer, inb - rf->buffer_len);
        rf->buffer_len = inb;
      }
    }

    // reading from the socket
    // set timeout to what user has provided us with
    timeout.tv_sec  = (long int) rf->socket_timeout;
    timeout.tv_usec = (long int) (1e6 * (double)(rf->socket_timeout - timeout.tv_sec));

    FD_ZERO(&socks);
    FD_SET (rf->fileds_i, &socks);
//     select(rf->fileds_i+1, &socks, NULL, NULL, &timeout);
    int ir = select(rf->fileds_i+1, &socks, NULL, NULL, &timeout);

    rent = ent_Create();
    if (ir <= 0)
    {
      if (ir < 0)
        fprintf(stderr, "readm: (select failed) with errno %i (%s)\n", errno, strerror(errno));

      ent_Clean(e1);
      ent_Clean(e2);

      w = mdr_Create(0,0);
      ent_data(rent) = w;
      ent_type(rent) = MATRIX_DENSE_REAL;
      return (rent);
    }

    // we need these for reading
    struct sockaddr address;
    socklen_t address_len;
    memset(&address, 0, sizeof(address));

    if (inb > 0)
    {
      // user wants us to read only 'inb' bytes from the socket

      // we peek first so we know how much data is coming
      res = recvfrom(rf->fileds_i, rf->buffer, rf->buffer_len-1, MSG_PEEK,
                     &address, &address_len);

      if (res > 0)
      {
        if (inb < res)
        {
          w = mdi_Create(1,inb);
          res = recvfrom(rf->fileds_i, rf->buffer, inb, 0, &address, &address_len);
          for (i=0; i<inb; i++)
          {
            MdiV0(w,i) = (int) rf->buffer[i];
            rf->buffer[i] = '\0';
          }
        }
        else
        {
          w = mdi_Create(1,res);
          res = recvfrom(rf->fileds_i, rf->buffer, res, 0, &address, &address_len);
          for (i=0; i<res; i++)
          {
            MdiV0(w,i) = (int) rf->buffer[i];
            rf->buffer[i] = '\0';
          }
        }

        // return values:
        ent_data(rent) = w;
        ent_type(rent) = MATRIX_DENSE_REAL;
      }
      else
      {
        w = mdi_Create(0,0);
        ent_type(rent) = MATRIX_DENSE_REAL;
      }
    }
    else
    {
      // inb=0, we read what ever we can and copy it to
      // the output string

      // we peek first so we know how much data is coming
      res = recvfrom(rf->fileds_i, rf->buffer, rf->buffer_len-1, MSG_PEEK,
                       &address, &address_len);
      if (res > 0)
      {
        rs = GC_malloc((res+1) * sizeof(char));
        res = recvfrom(rf->fileds_i, rs, res, 0, &address, &address_len);
        rs[res]='\0';

        // make a string out of it
        ws = mds_Create(1,1);
        MdsV0(ws,0) = rs;
      }
      else
      {
//         fprintf(stdout, "readm: (warning) Data expected but none received!\n");
        ws = mds_Create(0,0);
      }

      // return values:
      ent_data(rent) = ws;
      ent_type(rent) = MATRIX_DENSE_STRING;
    }

    if (close_after_rw)
      rfile_Destroy(fname);

    ent_Clean(e1);
    ent_Clean(e2);

    return (rent);
  }
#ifdef HAVE_HDF5_SO
  else if (rf->filetype == RFILE_H5)
  {
    int nobjs=0, k;
    MDR *coord=0;
    int coord_dim_space=0, coord_npoints=0;

    MDS  *objs=0;
    char *name=0;

    Btree *rtree=btree_Create();
    Btree *btree_prev=0;

    hid_t work_obj_id;
    H5I_type_t work_obj_type = H5I_BADID;
    H5T_class_t t_class = H5T_NO_CLASS;

    if (nargs < 2)
      rerror("readm: (h5) second argument 'dataset name' expected");

    // second argument is list of absolute paths of the objects which
    // are to be retrieved from file
    e2 = bltin_get_ent(args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      objs = ent_data ( e2 );
      nobjs = (objs->nrow) * (objs->ncol);
      if (!nobjs)
        rerror("readm: (h5) zero size of second argument 'dataset name'");
    }
    else
      rerror("readm: (h5) missing second argument 'dataset name'");

    // third argument is used only if a single object is given in as the
    // second argument: then it represents a list of coordinates of the entries
    // from the object that are to be retrieved (e.g., entry in a matrix
    // or a vector)
    if (nargs > 2 && nobjs==1)
    {
      e3 = bltin_get_ent(args[2]);
      if (ent_type(e3) == MATRIX_DENSE_REAL)
      {
        coord = ent_data ( e3 );
        coord_dim_space = coord->ncol;
        coord_npoints   = coord->nrow;
        if (coord_dim_space * coord_npoints == 0)
          rerror("readm: (h5) zero size of third argument 'coordinates'");
      }
    }
    else if (nargs > 2 && nobjs>1)
      rerror("readm: (h5) providing third argument is possible for single 'dataset name' only");

    rent = ent_Create();
    for(k=0; k<nobjs; k++)
    {
      // for each object accessed we need these:
      btree_prev=rtree;

      // iscompound object:
      //    /level{1}/level{2}/ .. /level{N}:compound_level{1}:compound_level{2}:...
      int compound_nlevels;
      MDS *compound_lev=0;

      MDR *w=0;
      MDS *ws=0;
      MDC *wc=0;
      MSR *wsr=0;
      MSC *wsc=0;

      int i=0, j=0, nlevel=1;
      char *rtree_label=0, *name2=0;

      // basic checking of the object name:
      //  it has to start with '/' = absolute path
      //  and it has to be greater in length than 1
      name = MdsV0(objs,k);

      if (name[0]!='/')
      {
        fprintf(stderr, "readm: (h5) only absolute path names allowed for 'dataset name'\n");
        continue;
      }
      if (strlen(name)==1) // we ignore '/' as a object name!
        continue;

      // go through name and obtain the count of levels of compound object
      // we will use that to obtain the name of the entries within the compound
      compound_nlevels = 0;

      // does object name exists in the lists of objects associated with
      // the current file name, and is it of desired type?
      name = MdsV0(objs,k);

      // Turn off error handling
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      work_obj_id = H5Oopen(rf->h5_file, MdsV0(objs,k), H5P_DEFAULT);
      if (work_obj_id > -1)
      {
        work_obj_type = H5Iget_type(work_obj_id);
      }
      else
      {
//         fprintf(stderr, "readm: Object %s not found in file!\n", MdsV0(objs,k));
        continue;
      }

      if (compound_nlevels)
      { *name = ':'; }

      // get the object ID from our list
      if (work_obj_type == H5I_DATASET)
      {
        if (nobjs > 1)
        {
          // there are more then one object being retrieved:
          // read the content of the data set into the local tree:
          // 1. create the tree from object name /level1/level2/..../level{N}
          // to
          //    res=level1.[level2]. .. .[level{N}]
          // determine no. of levels in the object name:
          name = MdsV0(objs,k);
          name++;
          for(i=0; i<strlen(name); i++)
          { if(name[i]=='/') nlevel++; }

          // now go down nlevels-1 and create a series of sub-btree's
          j = 0;
          name2 = name;
          for (i=0; i<nlevel-1;i++)
          {
            j=0;
            // this stops at the next '/' or at the end of the string
            while(name2[j]!='/' && name2[j]) j++;
            if (name2[j]=='/')
            {
              name2[j] = '\0';
              rtree_label = cpstr(name2);
              name2[j] = '/';
              j++;
            }
            else if (name2[j]=='\0')
              rtree_label = cpstr(name2);

            // does the entry exist in the current level of the tree?
            ListNode * node = btree_FindNode (btree_prev, rtree_label);
            if (!node)
            {
              // node does not exist. insert one using the label
              Ent * P = ent_Create ();
              ent_data (P) = btree_Create();
              ent_type (P) = BTREE;
              install (btree_prev, (rtree_label), P);
//               ent_Clean (P);

              // go to it
              node = btree_FindNode (btree_prev, rtree_label);
            }

            // continue the tree from existing node
            btree_prev = (Btree *) ent_data(node->ent);
            name2 = &name2[j];

          } //for (i=0; i<nlevel-1;i++)
        } // if (nobjs > 1)

        // what is the class of object:
        hid_t datatype = H5Dget_type(work_obj_id);
        if (datatype > -1)
        {
          t_class = H5Tget_class(datatype);
          H5Tclose(datatype);
        }
        if (    t_class == H5T_INTEGER
            ||  t_class == H5T_FLOAT
           )
        {
          // retrieve object if it is atomic data
          w  = (MDR *) get_hobject_atomic_data(work_obj_id, coord);
        }
        else if (t_class == H5T_STRING)
        {
          // retrieve object if it is atomic data
          ws = (MDS *) get_hobject_atomic_data(work_obj_id, coord);
        }
        else if (t_class == H5T_COMPOUND)
        {
            // treat compound special complex as an atomic data
          hid_t datatype = ascertain_hobject_is_compound_special_complex (work_obj_id);
          if (datatype > -1)
            wc = (MDC *) get_hobject_compound_special_complex(work_obj_id, datatype, coord);
          else
          {
              // try next rlab compound
            hid_t datatype1 = ascertain_hobject_is_compound_special_msr (work_obj_id);
            if (datatype1 > -1)
            {
              wsr = (MSR *) get_hobject_compound_special_msr (work_obj_id, datatype1);
            }
            else
            {
                // try next rlab compound
              hid_t datatype2 = ascertain_hobject_is_compound_special_msc (work_obj_id);
              if (datatype2 > -1)
              {
                wsc = (MSC *) get_hobject_compound_special_msc (work_obj_id, datatype2);
              }
            }
          }
        }

        // if we had only one object than prepare it for returning it to the
        // user
        if (nobjs == 1)
        {
          btree_Destroy (rtree);
          if (w)
          {
            ent_data(rent) = (w);
            ent_type(rent) = MATRIX_DENSE_REAL;
          }
          else if (wc)
          {
            ent_data(rent) = (wc);
            ent_type(rent) = MATRIX_DENSE_COMPLEX;
          }
          else if (ws)
          {
            ent_data(rent) = (ws);
            ent_type(rent) = MATRIX_DENSE_STRING;
          }
          else if (wsr)
          {
            ent_data(rent) = (wsr);
            ent_type(rent) = MATRIX_SPARSE_REAL;
          }
          else if (wsc)
          {
            ent_data(rent) = (wsc);
            ent_type(rent) = MATRIX_SPARSE_COMPLEX;
          }
          else
          {
            ent_data(rent) = mdr_Create(0,0);
            ent_type(rent) = MATRIX_DENSE_REAL;
          }
        }
        else
        {
          // at this point btree_prev points to the last btree to
          // which the data set will be attached
          Ent * P = ent_Create ();
          if (w)
          {
            ent_type (P) = MATRIX_DENSE_REAL;
            ent_data (P) = (w);
          }
          else if (wc)
          {
            ent_data(P) = (wc);
            ent_type(P) = MATRIX_DENSE_COMPLEX;
          }
          else if (ws)
          {
            ent_type (P) = MATRIX_DENSE_STRING;
            ent_data (P) = (ws);
          }
          else if (wsr)
          {
            ent_data(P) = (wsr);
            ent_type(P) = MATRIX_SPARSE_REAL;
          }
          else if (wsc)
          {
            ent_data(P) = (wsc);
            ent_type(P) = MATRIX_SPARSE_COMPLEX;
          }
          else
          {
            ent_data(P) = mdr_Create(0,0);
            ent_type(P) = MATRIX_DENSE_REAL;
          }
          install (btree_prev, name2, P);
          ent_Clean (P);
        }

      } // if (work_obj->object_type == H5I_DATASET)

      H5Oclose(work_obj_id);

      if (compound_lev)
        mds_Destroy (compound_lev);

    } // next object from user list: for(i=0; i<nobjs; i++)

    if (close_after_rw)
      close_file_ds(fname);

    // we had more then one object. we have thus built a BTREE
    if (nobjs > 1)
    {
      ent_data(rent) = rtree;
      ent_type(rent) = BTREE;
    }
    // rlab stuff:
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);

    return (rent);
  }
#endif
#ifdef HAVE_LIBCURL
  else if (rf->filetype == RFILE_CURL)
  {
    char *fn=0;

      // do we download to file or to the rlab variable?
    if (nargs > 1)
    {
      e2 = bltin_get_ent(args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror("readm: (libcurl) A file name is expected");
      fn = class_char_pointer(e2);
    }
    FILE *fn_fd=0;

    if (fn)
    {
        // the data from http/https is downloaded to the file

        // open file
      fn_fd = fopen(fn, "wb");

        // tell curl to get ready
      curl_easy_setopt(rf->curl, CURLOPT_WRITEFUNCTION, NULL);
      curl_easy_setopt(rf->curl, CURLOPT_WRITEDATA, fn_fd);
    }
    else
    {
        // the data from http/https is downloaded to the rlab variable

        // init local buffer for receiving the data
      rf->buffer = GC_malloc(CURL_MAX_NUM_BYTES * sizeof(char));
      rf->buffer[0] = '\0';
      rf->buffer_len = CURL_MAX_NUM_BYTES;
      rf->buffer_pos = 0;

        // tell curl to get ready:
      curl_easy_setopt(rf->curl, CURLOPT_WRITEDATA, rf);
      curl_easy_setopt(rf->curl, CURLOPT_WRITEFUNCTION, write_callback);
    }

      // i don't know why
    CURLcode res;
    MDS *w=0;
    double val;

    res = curl_easy_perform(rf->curl);

      // how did we do it?
    if (res == CURLE_OK)
    {
      if (CURLOPT_RLABPLUS_DEBUG)
      {
            // download time
        res = curl_easy_getinfo(rf->curl, CURLINFO_SIZE_DOWNLOAD, &val);
        if((CURLE_OK == res) && val)
          fprintf(stdout, "libcurl: Total download time: %0.3f sec.\n", val);
            // download speed
        res = curl_easy_getinfo(rf->curl, CURLINFO_SIZE_DOWNLOAD, &val);
        if((CURLE_OK == res) && val)
          fprintf(stdout, "Average download speed: %0.3f kbyte/sec.\n", val / 1024);
            // bytes downloaded
        res = curl_easy_getinfo(rf->curl, CURLINFO_SIZE_DOWNLOAD, &val);
        if((CURLE_OK == res) && val)
          fprintf(stdout, "Data downloaded: %0.0f bytes.\n", val);
      }

        // figure out the report to the user:
      if (fn)
      {
          // nothing to say
        w = mds_Create(0,0);
      }
      else
      {
          // return the retrieved data as text. this is a trouble because
          // it cannot contain 0's (->end of string).
          // caveat emptor
        if (strlen(rf->buffer)>0)
          w = mds_CreateScalar( rf->buffer );
        else
          w = mds_Create(0,0);
      }
    }
    else
    {
      //
      fprintf(stderr, "libcurl: Error while fetching '%s' : %s\n",
              rf->name, (char *)curl_easy_strerror(res));
      w  = mds_Create(0,0);
    }

    if (!fn)
    {
      GC_free (rf->buffer);
      rf->buffer = 0;
      rf->buffer_len = 0;
      rf->buffer_pos = 0;
    }
    else
    {
      fflush(fn_fd);
      fclose(fn_fd);
    }

    if (close_after_rw)
      rfile_Destroy(fname);

      // rlab stuff:
    ent_Clean (e1);
    ent_Clean (e2);

    rent = ent_Create();
    ent_data(rent) = w;
    ent_type(rent) = MATRIX_DENSE_STRING;
    return (rent);
  }
#endif

  //
  // filename is regular file or a pipe
  //
  MDS *comment=0, *note=0, *lstrip=0, *grep=0;
  MDR *userows=0, *usecols=0;
  char *delim=0, *join_csp=0;
  int min_line_len=1, iskiprows=0, join_rows=1;
  MDS *start=0, *stop=0;

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      iskiprows = (int) class_double (e2);
    }
    else if (ent_type (e2) == BTREE)
    {
      ListNode *
          node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_CSP);

      if (node)
      {
        delim = (char *) class_char_pointer (var_ent (node));
        if (isvalidstring(delim)<1)
          delim = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_MIN_LINE_LEN);
      if (node)
      {
        min_line_len = class_int (var_ent (node));
        if (min_line_len < 1 )
          min_line_len = 1;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_JOINROWS);
      if (node)
      {
        join_rows = class_int (var_ent (node));
        if (join_rows < 1 )
          join_rows = 1;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_JOINCSP);
      if (node)
      {
        join_csp = (char *) class_char_pointer (var_ent (node));
        if (isvalidstring(join_csp)<1)
          join_csp = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_COMMENT);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          comment = ent_data(var_ent (node));
          if (SIZE(comment)<1)
            comment=0;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_NOTE);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          note = ent_data(var_ent (node));
          if (SIZE(note)<1)
            note=0;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_LSTRIP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          lstrip = ent_data(var_ent (node));
          if (SIZE(lstrip)<1)
            lstrip=0;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_START);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          start = ent_data(var_ent (node));
          if (SIZE(start)<1)
            start=0;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_STOP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          stop = ent_data(var_ent (node));
          if (SIZE(stop)<1)
            stop=0;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_SKIPROWS);
      if (node)
      {
        iskiprows = class_int (var_ent (node));
        if (iskiprows <= 0 )
          iskiprows = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_USEROWS);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          userows = ent_data(var_ent (node));
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_USECOLS);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          usecols = ent_data(var_ent (node));
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_GREP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          grep = ent_data(var_ent (node));
        }
      }
    }
  }

  if (nargs >= 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
    {
      //
      // duplicate entity if it is refered to by more than once. we are changing
      // its content which would affect all variables pointing to it.
      //
      ListNode *var = (ListNode *) (args[2].u.ptr);
      e3 = var_ent (var);
      if (e3->refc > 1)
      {
        ent_DecRef (e3);
        e3 = ent_Duplicate (e3);
        listNode_AttachEnt (var, e3);
      }

      // Get block_size from argument list
      a2 = (MDR *) ent_data (e3);
      int nr = (a2->nrow);
      int nc = (a2->ncol);
      if (nr * nc == 1)
        block_size = (int) MdrV0 (a2, 0);
      else if ((a2->nrow) * (a2->ncol) >= 2)
      {
        /* Just do it here. */
        /* Try and open() the file */

        if (!rf->fileds_f)
        {
          fprintf (stdout, "%s, cannot open for read\n", fname);
          rent = ent_Create ();
          ent_data (rent) = mdr_CreateScalar (0.0);
          ent_SetType (rent, MATRIX_DENSE_REAL);
          return (rent);
        }
        m = mdr_CreateScalar (1.0);
        // do the offset
        if (iskiprows)
        {
          c = fgetc (rf->fileds_f);
          while (c!=EOF)
          {
            // skip first iskiprows lines
            if (c == '\n')
              iskiprows--;
            if (iskiprows)
              c = fgetc (rf->fileds_f);
            else
              break;
          }
        }
        for (i = 0; i < nr; i++)
          for (j = 0; j < nc; j++)
        {
          // read the data
          rv = fscanf (rf->fileds_f, "%lf", &Mdr0 (a2, i, j));
          if (rv == EOF || rv == 0)
            Mdr0 (a2, i, j) = create_nan ();
          // read the next character (should work on a single character separators)
          fscanf (rf->fileds_f, "%c", &c);
        }

        ent_Clean (e1);
        ent_Clean (e3);

        if (close_after_rw)
          rfile_Destroy(fname);

        rent = ent_Create ();
        ent_data (rent) = m;
        ent_SetType (rent, MATRIX_DENSE_REAL);
        return (rent);
      }
    }
  }

  /* Try and open() the file */
  if (!rf->fileds_f)
  {
    fprintf (stderr, "%s, cannot open for read\n", fname);
    goto _exit_readm;
  }

  // define default delimiters if user failed to do so
  if (join_rows>1 && !join_csp)
    join_csp = RLAB_READM_DEFAULT_JOINCSP;
  if (!delim)
    delim = RLAB_READM_DEFAULT_CSP;

  if ((m = mdr_ReadGeneric (rf->fileds_f, block_size, iskiprows, delim, min_line_len,
       join_rows, join_csp, comment, note, lstrip, grep, userows,
       usecols, start, stop)) == 0)
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_READM_GENERIC "\n");

_exit_readm:

  if (close_after_rw)
    rfile_Destroy(fname);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(m);
}

#ifdef HAVE_LIBCURL
static size_t
    read_callback(void *ptr, size_t size, size_t nmemb, void *stream)
{
  //
  char *data = (char *) stream;
  int  n;

  n = strlen(data);
  if (n < size * nmemb)
  {
    memcpy(ptr, data, n);
    return 0;
  }
  else
  {
    memcpy(ptr, data, size * nmemb);
    // is this legal? we are offseting position of stream so that if another calls
    // follows, we do not send the same data
//     &stream = &stream[size * nmemb];
    return (n-(size * nmemb));
  }
}
#endif

#undef  THIS_SOLVER
#define THIS_SOLVER "writem"
Ent *
WriteM (int nargs, Datum args[])
{
  ListNode *node=0;
  char *name;
  int type;
  Rfile * rf;
  Ent *e1=0, *e2=0, *e3=0, *rent;
  int close_after_rw=0, i;
  void (*vfptr) ();

  /* Check n_args */
  if (nargs < 2)
    rerror (THIS_SOLVER ": at least 2 arguments required");

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": first argument must be a single string");
  name = class_char_pointer( e1);
  if (!name)
    rerror (THIS_SOLVER ": empty filename");
  if (!strlen(name))
    rerror (THIS_SOLVER ": empty filename");

  // was the file previously accessed?
  rf = rfile_find(name);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(name, "h5://", 5)
        ||  !strncmp(name, "hdf5://", 7)
        ||  (strstr (name, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (name, RFILE_H5, "w" );
    }
    else if (!strncmp(name, "tcp://", 6))
    {
      // do the socket
      rf = get_rfile_ds (name, RFILE_SOCKET, 1);
    }
    else
      if (    (!strncmp(name, "http://", 7))
          ||  (!strncmp(name, "https://", 8))
          ||  (!strncmp(name, "ftp://", 6))
       )
    {
      // initialize CURL with 'name' if it already does not exist
      rf = get_rfile_ds(name, RFILE_CURL);
    }
    else
    {
      // just assume it is a text file
      rf = get_rfile_ds (name, RFILE_FILE, "w", 0);
    }
    close_after_rw = 1;
  }

  if (!rf)
  {
    fprintf (stderr, "writem: Trying to write to URL that is not accessible\n");
    rerror ("Cannot write to URL");
  }

  //
  // now check the second argument
  //
  e2 = bltin_get_ent (args[1]);
  type = ent_type (e2);
  if ((type == UNDEF) && rf->filetype != RFILE_H5)
    rerror("writem: second argument 'data' cannot be undefined");

  if (rf)
  {
    if (rf->filetype == RFILE_COMM)
    {
      // user is trying to write to a previously opened file.
      // check if the file's int fd exists (serial port) and
      // use write() instead of fprintf() to send the data.
      // rf exists:
      if (rf->fileds_i<0)
        rerror("writem: Terrible internal error! Writing to a non-existing serial port");

      if (type == MATRIX_DENSE_STRING)
        mds_intfd_WriteGeneric(name, ent_data(e2));
      else if (type == MATRIX_DENSE_REAL)
        mdi_intfd_WriteGeneric(name, ent_data(e2));
      else if (type == BTREE)
      {
        // DTR
        node = btree_FindNode (ent_data (e2), RLAB_NAME_WRITEM_SETDTR);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            i = (int) class_int(var_ent(node));
            if (i==0 || i==1 || i==-1)
              serial_SetDTR (rf, i);
          }
        }
        // RTS
        node = btree_FindNode (ent_data (e2), RLAB_NAME_WRITEM_SETRTS);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            i = (int) class_int(var_ent(node));
            if (i==0 || i==1 || i==-1)
              serial_SetRTS(rf, i);
          }
        }
      }
      else
        rerror("writem: only strings or integers can be written to serial port");

      ent_Clean (e1);
      ent_Clean (e2);

      rent = ent_Create ();
      ent_data (rent) = mdr_CreateScalar (0.0);
      ent_SetType (rent, MATRIX_DENSE_REAL);
      return (rent);
    }
#ifdef HAVE_HDF5_SO
    else if (rf->filetype == RFILE_H5)
    {
      int ires = 0;
      char *obj_name=0;

      // check if the location for the variable is given, otherwise use its name
      if (nargs == 3)
      {
        e3 = bltin_get_ent (args[2]);
        if (ent_type(e3) == MATRIX_DENSE_STRING)
          obj_name = class_char_pointer(e3);
      }
      else
      {
        if (args[1].type == VAR)
          obj_name = var_key (args[1].u.ptr);
      }

      if (!obj_name)
        rerror ("writem: (h5) empty object name: are you trying to save expression?");
      if (!strlen(obj_name))
        rerror ("writem: (h5) object name of zero length: what are you doing?");

      // delete object if it exists: but only if the file has been opened before
      if (!close_after_rw)
      {
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);
        H5Ldelete(rf->h5_file, obj_name, H5P_DEFAULT);
      }

      switch( type )
      {
        case MATRIX_DENSE_REAL:
          h5_write_rlab_atomic_mdr (rf, obj_name, ent_data(e2));
          break;

        case MATRIX_DENSE_COMPLEX:
          h5_write_rlab_atomic_mdc (rf, obj_name, ent_data(e2));
          break;

        case MATRIX_DENSE_STRING:
          h5_write_rlab_atomic_mds (rf, obj_name, ent_data(e2));
          break;

        case MATRIX_SPARSE_REAL:
          h5_write_rlab_atomic_msr (rf, obj_name, ent_data(e2));
          break;

        case MATRIX_SPARSE_COMPLEX:
          h5_write_rlab_atomic_msc (rf, obj_name, ent_data(e2));
          break;

        case UNDEF:
          // nothing to do: user called WriteM to delete object from HDF file
          break;

        default:
          fprintf(stdout, "writem (HDF5): can save only atomic objects (are you trying to save list?)");
      }

      if (close_after_rw)
        rfile_Destroy(name);

      // rlab stuff:
      ent_Clean (e1);
      ent_Clean (e2);
      ent_Clean (e3);

      rent = ent_Create();
      ent_data(rent) = mdr_CreateScalar(ires);
      ent_type(rent) = MATRIX_DENSE_REAL;
      return (rent);
    }
#endif
#ifdef HAVE_LIBCURL
    else if (rf->filetype == RFILE_CURL)
    {
      // writes a file or a content of a variable to the internet site using
      // http, https or ftp protocols
      char *fn=0;
      FILE *fn_fd=0;
      char  *msg=0, *msg2=0;

      // do we upload a file (given as a second argument) or a content of an
      // rlab variable (third argument) ?
      if (nargs == 2)
      {
/*        e2 = bltin_get_ent(args[1]);*/
        if (ent_type(e2) != MATRIX_DENSE_STRING)
          rerror("readm: (libcurl) File name is expected as a second argument");
        fn = class_char_pointer(e2);
      }
      else if (nargs == 3)
      {
        e3 = bltin_get_ent(args[2]);
        if (ent_type(e3) != MATRIX_DENSE_STRING)
          rerror("readm: (libcurl) String is expected as a third argument");
        msg  = class_char_pointer (e3);
        msg2 = msg;
      }

      if (fn)
      {
        // a file is uploaded to the http/https/ftp site
        // open file
        if ((fn_fd = fopen(fn, "rb")))
        {
          // tell curl to get ready
          curl_easy_setopt(rf->curl, CURLOPT_READFUNCTION, NULL);
          curl_easy_setopt(rf->curl, CURLOPT_READDATA, fn_fd);
        }
      }
      else if (msg)
      {
        // the data is uploaded from the rlab variable to the http/https/ftp site
        curl_easy_setopt(rf->curl, CURLOPT_READDATA, msg);
        curl_easy_setopt(rf->curl, CURLOPT_READFUNCTION, read_callback);
      }
      else
        rerror("readm: (libcurl) Terrible internal error. Sigh.");

      CURLcode res = curl_easy_perform(rf->curl);
      msg = msg2; // we fiddled with 'msg' in the function 'read_callback'
      MDR *w = mdr_CreateScalar(res);

      if (fn_fd)
        fclose(fn_fd);

      // rlab stuff:
      ent_Clean (e1);
      ent_Clean (e2);
      ent_Clean (e3);

      rent = ent_Create();
      ent_data(rent) = w;
      ent_type(rent) = MATRIX_DENSE_REAL;
      return (rent);
    }
#endif
    else if (rf->filetype == RFILE_SOCKET)
    {
      //
      // write second argument to the socket
      //

      if (type == MATRIX_DENSE_STRING)
      {

        //
        // writing a string matrix to a socket
        //

        MDS *s = ent_data(e2);
        int nr = s->nrow;
        int nc = s->ncol;
        int i, j, k;
        char *data;
        int data_len, res=0;

        for (i=0; i<nr; i++)
        {
          for (j=0; j<nc; j++)
          {
            // specify the data
            data = Mds0(s,i,j);
            if (data)
            {
              data_len = strlen(Mds0(s, i, j));
              k = send(rf->fileds_i, data, data_len, 0);
              if (k != data_len)
              {
                fprintf(stderr,"Mismatch in number of sent bytes length=%i, sent=%i\n", data_len, k);
              }
            }
            else
              fprintf(stderr,"Attempting to write null pointer. Is there a problem with string matrix?\n");
          }
        }

        if (close_after_rw)
          rfile_Destroy(name);

        // rlab stuff:
        ent_Clean (e1);
        ent_Clean (e2);

        rent = ent_Create();
        ent_data(rent) = mdr_CreateScalar(res);
        ent_type(rent) = MATRIX_DENSE_REAL;
        return (rent);

      } // MDS
      else if (type == MATRIX_DENSE_REAL)
      {
        //
        // writing a real matrix to a socket
        //
        MDR *r = ent_data(e2);
        int nr = r->nrow;
        int nc = r->ncol;
        int i, j;
        unsigned char *data;
        int data_len, res=0;

        for (i=0; i<nr; i++)
        {
          for (j=0; j<nc; j++)
          {
            // specify the data
            if (r->type == RLAB_TYPE_INT32)
            {
              data = (unsigned char*) &Mdi0(r, i, j);
              data_len = sizeof(int);
            }
            else
            {
              data = (unsigned char*) &Mdr0(r, i, j);
              data_len = sizeof(double);
            }

            if (send(rf->fileds_i, data, data_len, 0) != data_len)
            {
              fprintf(stderr,"Mismatch in number of sent bytes\n");
            }
          }
        }

        if (close_after_rw)
          rfile_Destroy(name);

        // rlab stuff:
        ent_Clean (e1);
        ent_Clean (e2);

        rent = ent_Create();
        ent_data(rent) = mdr_CreateScalar(res);
        ent_type(rent) = MATRIX_DENSE_REAL;
        return (rent);
      }
    }
  }

  //
  // writing to a file or a pipe
  //
  char *eol=0, *csp=0, *nan=0, *inf_pos=0, *inf_neg=0;
  MDS *fmt=0;
  if (nargs == 3)
  {
    e3 = bltin_get_ent(args[2]);
    if (ent_type (e3) == BTREE)
    {
      // eol
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_ENDOFLINE);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          eol = class_char_pointer (var_ent(node));
      }
      // nan
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_NAN);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          nan = class_char_pointer (var_ent(node));
      }
      // inf_pos
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_INF_POS);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          inf_pos = class_char_pointer (var_ent(node));
      }
      // inf_neg
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_INF_NEG);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          inf_neg = class_char_pointer (var_ent(node));
      }
      // csp
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_COLSEP);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          csp = class_char_pointer (var_ent(node));
      }
      // fmt
      node = btree_FindNode (ent_data (e3), RLAB_NAME_WRITEM_FORMAT);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
          fmt = ent_data(var_ent(node));
      }
    }
  }

  if (nan)
  {
    if (rf->nan)
      GC_free (rf->nan);
    rf->nan = cpstr(nan);
  }
  else
    rf->nan = cpstr(RLABPLUS_FILE_NAN);

  if (inf_pos)
  {
    if (rf->inf_pos)
      GC_free (rf->inf_pos);
    rf->inf_pos = cpstr(inf_pos);
  }
  else
    rf->inf_pos = cpstr(RLABPLUS_FILE_INF_POS);

  if (inf_neg)
  {
    if (rf->inf_neg)
      GC_free (rf->inf_neg);
    rf->inf_neg = cpstr(inf_neg);
  }
  else
    rf->inf_neg = cpstr(RLABPLUS_FILE_INF_NEG);

  if (eol)
  {
    if (rf->eol)
      GC_free (rf->eol);
    rf->eol = cpstr(eol);
  }
  else
    rf->eol = cpstr(RLABPLUS_FILE_EOL_UNIX);

  if (csp)
  {
    if (rf->csp)
      GC_free (rf->csp);
    rf->csp = cpstr(csp);
  }
  else
    rf->csp = cpstr(RLABPLUS_FILE_CSP);

  if (fmt)
  {
    if (rf->fmt)
      mds_Destroy (rf->fmt);
    rf->fmt = mds_Copy( fmt );
  }
  else
    rf->fmt = mds_CreateScalar( RLABPLUS_FILE_FMT );

  vfptr = (VVFPTR) writem_method[type].op;
  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e2));
    rerror ("writem() operation not supported");
  }

  rent = ent_Create ();
  if (!rf)
  {
    fprintf (stderr, "%s, cannot open for write\n", name);
    ent_data (rent) = mdr_CreateScalar (RLAB_STATUS_FAILURE);
  }
  else
  {
    (*vfptr) (ent_data (e2), rf);
    ent_data (rent) = mdr_CreateScalar (RLAB_STATUS_SUCCESS);
  }

  if (close_after_rw)
    rfile_Destroy(name);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

/* **************************************************************
 * Binary input / output related functions.
 * ************************************************************** */

#ifdef HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>
#endif

#ifndef HAVE_FREAD_DEC
extern size_t fread ();
#endif

/*
 * Set some fixed sizes.
 * Try to do this in a way that will work on
 * older C compilers.
 */

#if SIZEOF_LONG_INT == 4
typedef long int FOURB;
#else
#if SIZEOF_INT == 4
typedef int FOURB;
#else
#if SIZEOF_SHORT_INT == 4
typedef short int FOURB;
#else
Cannot build proper binary file read functions ... without having a
  four - byte size data type ...
#endif				/* SHORT_INT */
#endif				/* INT */
#endif				/* LONG_INT */
typedef size_t (*FREADER) ();

size_t fread_swap (void *ptr, size_t size, size_t nitems, FILE * stream);

FREADER freadptr;		/* Use this to point one or another fread variants */

#define FREAD(ptr, size, n, fn) \
        if ((*freadptr) (ptr, size, n, fn) != ((size_t) n) ) \
          if (!ferr_check ("readb", fn)) { return (0); }

/*
 * Set a component of type that determines whether this
 * is a big or little endian word machine.
 */

#ifdef WORDS_BIGENDIAN
static FOURB word = 1000;
#else
static FOURB word = 0;
#endif

/*
 * Integers used for entity identification
 */

static FOURB SM = 5;		/* String Matrix */
static FOURB BT = 6;		/* Binary Tree */
static FOURB MT = 0;		/* MaTrix */
static FOURB MI = 20;		/* MaTrix */
static FOURB SP = 7;		/* SParse matrix */

//
// Binary write
// The interface writes out object(s), using the correct
// binary-write method.
//
// Each object it completely responsible for writing/reading itself..
//
// Works with two file formats:
// (1) standard, or the old one
// (2) HDF5
#undef THIS_SOLVER
#define THIS_SOLVER "writeb"
Ent *
WriteB (int nargs, Datum args[])
{
  Ent *e1, *e, *rent;
  ListNode *lnode;
  char *filename, *name;
  int i, type;
  int close_after_rw = 0;
  Rfile *rf=0;

  void (*vfptr) ();

  if (nargs < 2)
    rerror (THIS_SOLVER ": two or more arguments required");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror(THIS_SOLVER ": First argument 'filename' must be string scalar !");
  filename = class_char_pointer (e1);
  if (!filename)
    rerror(THIS_SOLVER ": First argument 'filename' must be string scalar !");

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "w" );
    }
    else
    {
      // must be local/regular file then
      rf = get_rfile_ds(filename, RFILE_FILE, "wb", 0);
    }
    close_after_rw = 1;
  }

  if (!rf)
  {
    fprintf (stderr, THIS_SOLVER ": Cannot open file %s for write", filename);

    ent_Clean (e1);
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_SetType (rent, MATRIX_DENSE_REAL);
    return (rent);
  }

  // Loop over the number of arguments, writing them out...
  for (i = 1; i < nargs; i++)
  {
    if (args[i].type != VAR)
      rerror ("writeb: arguments must be variables");

    lnode = (ListNode *) args[i].u.ptr;
    e = convert_datum_to_ent (args[i]);
    name = var_key (lnode);
    type = ent_type (e);

    // Write the data (finally)
    vfptr = (VVFPTR) writeb_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, "writeb: Entity type: %s\n", etd (e));
      rerror ("writeb: Operation not supported");
    }

    (*vfptr) (rf, name, ent_data (e));

    ent_Clean (e);
  }

  if (close_after_rw)
    rfile_Destroy(filename);

  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

//
// Binary read: regular and HDF5 files
//

// ****************************************************************************
// HDF5 functionality based on h5ex_traverse.c, part of C-API available
// from
//    http://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/api16-c.html
// ****************************************************************************

//   Define operator data structure type for H5Giterate callback.
//   During recursive iteration, these structures will form a
//   linked list that can be searched for duplicate groups,
//   preventing infinite recursion.
struct opdata
{
  Rfile           * rfile;        // Rfile being processed
  Btree           * btree;        // rlab structure to pick up data
  char            * name;
  unsigned          recurs;       // recursion level.  0=root
  struct opdata   * prev;         // pointer to previous opdata
  unsigned long     groupno[2];   // unique group number
};

//  This function recursively searches the linked list of
//  opdata structures for one whose groupno field matches
//  target_groupno.  Returns 1 if a match is found, and 0
//  otherwise.
static int
    group_check (struct opdata *od, unsigned long target_groupno[2])
{
  if (    (od->groupno[0] == target_groupno[0])
      &&  (od->groupno[1] == target_groupno[1])
     )
    return 1; // Group numbers match
  else if (!od->recurs)
    return 0; // Root group reached with no matches
  else
    return group_check (od->prev, target_groupno); // Recursively examine the previous node
}

//
// Operator function to be called by H5Giterate.
//
static herr_t
h5_rlab_op_func (hid_t loc_id, const char *name, void *operator_data)
{
  herr_t          return_val = 0;
  H5G_stat_t      statbuf;
  struct opdata  *od = (struct opdata *) operator_data;

  //
  // Get type of the object and its name
  //
  H5Gget_objinfo (loc_id, name, 0, &statbuf);
  if (statbuf.type == H5G_GROUP)
  {
      // Check group objno against linked list of operator
      // data structures.  Only necessary if there is more
      // than 1 link to the group.
      if (    (statbuf.nlink > 1)
          &&  group_check (od, statbuf.objno)
         )
      {
        fprintf (stderr, "Warning! The HDF5 file contains loops!\n");
      }
      else
      {
        // Initialize new operator data structure and
        // begin recursive iteration on the discovered
        // group.  The new opdata structure is given a
        // pointer to the current one.
        struct opdata nextod;
        od->name = (char *) name;
        nextod.rfile  = od->rfile;
        nextod.btree  = od->btree;
        nextod.recurs = od->recurs + 1;
        nextod.prev = od;
        nextod.groupno[0] = statbuf.objno[0];
        nextod.groupno[1] = statbuf.objno[1];
        return_val = H5Giterate (loc_id, name, NULL, h5_rlab_op_func,
                                 (void *) &nextod);
      }
  }
  else if (statbuf.type == H5G_DATASET)
  {
    // We install only the atomic data, and create necessary groups as
    // the list super-entries

    // But first, obtain the names of the groups up the chain
    // (ah, witness the beauty of recursive functions!)
    Btree *symtab = od->btree;
    int i, ng = od->recurs+1;
    MDS * obj_names = mds_Create(1,ng);
    struct opdata  *x = od;
    MdsV1(obj_names, ng) = cpstr((char*)name);
    i = ng-1;
    while  (x->prev)
    {
      if (x->prev->name)
        MdsV1(obj_names,i) = cpstr(x->prev->name);
      i--;
      x = x->prev;
    }

    // insert the upper levels into the 'symtab'
    ListNode *node;
    for (i=1; i<=ng-1;i++)
    {
      // create the tree in 'symtab' if it already does not exist
      node = btree_FindNode (symtab, MdsV1(obj_names,i));
      if (!node)
      {
        // node does not exist. insert one using the label
        Ent * P = ent_Create ();
        ent_data (P) = btree_Create();
        ent_type (P) = BTREE;
        install (symtab, MdsV1(obj_names,i), P);

        // now, go to it
        node = btree_FindNode (symtab, MdsV1(obj_names,i));
      }
      // at this point 'node' exists. is it a BTREE?
      if (ent_type(node->ent) != BTREE)
      {
        // It is not: Destroy it, and replace it with a list.
        Ent *E = var_ent(node->ent);
        ent_Destroy (E);

        E = ent_Create ();
        ent_data (E) = btree_Create ();
        ent_type (E) = BTREE;
        ent_IncRef (E);
        listNode_AttachEnt (node, E);
      }

      // continue the tree from existing node
      symtab = (Btree *) ent_data(node->ent);
    }
    node = btree_FindNode (symtab, MdsV1(obj_names,ng));

    // now we need to create the true object name for accessing it
    // at the end we should have:
    // obj = /name1/name2.../nameN
    int sl = ng+1;
    for (i=1; i<=ng;i++)
      sl += strlen(MdsV1(obj_names,i));
    char * obj = (char *) GC_malloc(sl*sizeof(char));
    int ic=0;
    for (i=1; i<=ng;i++)
    {
      obj[ic] = '/';
      ic ++;
      memcpy(&obj[ic], MdsV1(obj_names,i), strlen(MdsV1(obj_names,i)));
      ic += strlen(MdsV1(obj_names,i));
    }
    obj[ic] = '\0';

    hid_t new = H5Oopen(od->rfile->h5_file, obj, H5P_DEFAULT);
    if (!new)
    {
      GC_free (obj);
      return return_val;
    }

    //
    // at this point we should install in 'symtab' the content of the atomic data
    //
    MDR *w=0;
    MDC *wc=0;
    MDS *ws=0;
    MSR *wsr=0;
    MSC *wsc=0;

    H5T_class_t t_class=H5T_NO_CLASS;

    // what is the class of object:
    hid_t      datatype = H5Dget_type(new);
    if (datatype > -1)
    {
      t_class = H5Tget_class(datatype);
      H5Tclose(datatype);
    }

    if (    t_class == H5T_INTEGER
        ||  t_class == H5T_FLOAT
       )
    {
      // retrieve object if it is atomic data
      w  = (MDR *) get_hobject_atomic_data(new, NULL);
    }
    else if (t_class == H5T_STRING)
    {
      // retrieve object if it is atomic data
      ws = (MDS *) get_hobject_atomic_data(new, NULL);
    }
    else if (t_class == H5T_COMPOUND)
    {
      // treat compound special complex as an atomic data
      hid_t datatype = ascertain_hobject_is_compound_special_complex (new);
      if (datatype > -1)
        wc = (MDC *) get_hobject_compound_special_complex(new, datatype, NULL);
      else
      {
        // try next rlab compound
        hid_t datatype1 = ascertain_hobject_is_compound_special_msr (new);
        if (datatype1 > -1)
        {
          wsr = (MSR *) get_hobject_compound_special_msr (new, datatype1);
        }
        else
        {
          // try next rlab compound
          hid_t datatype2 = ascertain_hobject_is_compound_special_msc (new);
          if (datatype2 > -1)
            wsc = (MSC *) get_hobject_compound_special_msc (new, datatype2);
        }
      }
    }

    H5Oclose(new);

    if (w)
    {
      if (node)
        ent_Destroy (var_ent (node));
      Ent * Q = ent_Create ();
      ent_data (Q) = w;
      ent_type (Q) = MATRIX_DENSE_REAL;
      if (node)
      {
        ent_IncRef (Q);
        listNode_AttachEnt (node, Q);
      }
      else
        install (symtab, (MdsV1(obj_names,ng)), Q);
    }
    else if (wc)
    {
      if (node)
        ent_Destroy (var_ent (node));
      Ent * Q = ent_Create ();
      ent_data (Q) = wc;
      ent_type (Q) = MATRIX_DENSE_COMPLEX;
      if (node)
      {
        ent_IncRef (Q);
        listNode_AttachEnt (node, Q);
      }
      else
        install (symtab, MdsV1(obj_names,ng), Q);
    }
    else if (ws)
    {
      if (node)
        ent_Destroy (var_ent (node));
      Ent * Q = ent_Create ();
      ent_data (Q) = ws;
      ent_type (Q) = MATRIX_DENSE_STRING;
      if (node)
      {
        ent_IncRef (Q);
        listNode_AttachEnt (node, Q);
      }
      else
        install (symtab, MdsV1(obj_names,ng), Q);
    }
    else if (wsr)
    {
      if (node)
        ent_Destroy (var_ent (node));
      Ent * Q = ent_Create ();
      ent_data (Q) = wsr;
      ent_type (Q) = MATRIX_SPARSE_REAL;
      if (node)
      {
        ent_IncRef (Q);
        listNode_AttachEnt (node, Q);
      }
      else
        install (symtab, (MdsV1(obj_names,ng)), Q);
    }
    else if (wsc)
    {
      if (node)
        ent_Destroy (var_ent (node));
      Ent * Q = ent_Create ();
      ent_data(Q) = (wsc);
      ent_type(Q) = MATRIX_SPARSE_COMPLEX;
      if (node)
      {
        ent_IncRef (Q);
        listNode_AttachEnt (node, Q);
      }
      else
        install (symtab, (MdsV1(obj_names,ng)), Q);
    }

    mds_Destroy(obj_names);
    GC_free (obj);
  }

  return return_val;
}

#include "rfileio_h5ls.c"

#undef THIS_SOLVER
#define THIS_SOLVER "readb"
Ent *
ReadB (int nargs, Datum args[])
{
  Ent *ent, *e1, *e2, *rent;
  int P, T, rtype;
  int close_after_rw = 0;
  char *name, *filename;
  Rfile *rf=0;
  Datum arg2;
  size_t stat;
  Btree *btree, *symtab;
  ListNode *lnode;
  void *m;

  FOURB type;

  /* Check n_args */
  if (nargs > 2 || nargs == 0)
    rerror (THIS_SOLVER ": One or two arguments required");

  /* get file name from argument list, always the 1st argument */
  e1 = bltin_get_ent (args[0]);
  filename = class_char_pointer (e1);
  if (!filename)
    rerror (THIS_SOLVER ": First argument 'filename' must be string scalar !");

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "r" );
    }
    else
    {
      // must be local/regular file then
      rf = get_rfile_ds(filename, RFILE_FILE, "rb", 0);
    }
    close_after_rw = 1;
  }
  if (!rf)
  {
    fprintf (stderr, THIS_SOLVER  ": Cannot open file %s for read", filename);

    ent_Clean (e1);

    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_SetType (rent, MATRIX_DENSE_REAL);
    return (rent);
  }

  if (nargs == 1)
  {
    /* Use the global symbol table */
    symtab = get_symtab_ptr ();
  }
  else
  {
    /* Use a user-supplied list */
    arg2 = args[1];

    /* Check to see if variable exists, etc... */
    if (arg2.type != VAR)
      rerror ("readb: 2nd argument must be a variable");

    // If the specified variable is already a list, then just
    // create (or overwrite) members. If it is not a list, then
    // destroy it, and create a new list.
    e2 = var_ent (arg2.u.ptr);
    if (ent_type (e2) == BTREE)
    {
      /* Create a new member. */
      symtab = ent_data (e2);
    }
    else
    {
      /* Destroy it, and create a list. */
      ent_Destroy (e2);

      e2 = ent_Create ();
      symtab = btree_Create ();
      ent_data (e2) = symtab;
      ent_SetType (e2, BTREE);
      ent_IncRef (e2);
      listNode_AttachEnt (arg2.u.ptr, e2);
    }
  }

  if (rf->filetype == RFILE_H5)
  {

    // read HDF5 file either into the global symbol table (workspace)
    // or into the given variable:
    // iterate over the content of HDF5 file and load all atomic datasets (not compounds)
    // using the conversion:
    //  /name0/.../nameN -> name0.name1.name.....nameN
    int ires;
    H5G_stat_t      statbuf;
    struct opdata   od;

    // Open the root, and find its info
    ires = H5Gget_objinfo (rf->h5_file, "/", 0, &statbuf);

    od.rfile  = rf;
    od.btree  = symtab;
    od.name   = "/";
    od.recurs = 0;
    od.prev = NULL;
    od.groupno[0] = statbuf.objno[0];
    od.groupno[1] = statbuf.objno[1];

    ires = H5Giterate (rf->h5_file, "/", NULL, h5_rlab_op_func, (void *) &od);

    if (close_after_rw)
      rfile_Destroy(filename);

    ent_Clean (e1);

    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (ires);
    ent_SetType (rent, MATRIX_DENSE_REAL);
    return (rent);
  }
  else if (rf->filetype == RFILE_FILE)
  {
    //
    // binary reads starts with:
    // Read the first 4 bytes, and deciper the code...
    //
    while ((stat = fread (&type, sizeof (FOURB), 1, rf->fileds_f)) == 1)
    {
      if (type < 0 || type > 9999)
      {
        /* try flipping word to see if we can fix it? */
        reverse_word (&type, sizeof (FOURB));
      }

      /* Check for BIG or LITTLE Endian. */
      if (type >= 1000 && type < 2000)  /* BIG Endian */
      {
        if (word == 1000)   /* RLaB is running BIG. */
        {
          freadptr = (FREADER) fread;
        }
        else if (word == 0) /* RLaB is running LITTLE. */
        {
          freadptr = (FREADER) fread_swap;
        }
      }
      else if (type < 1000 && type >= 0)  /* LITTE Endian */
      {
        if (word == 1000)   /* RLaB is running BIG. */
        {
          freadptr = (FREADER) fread_swap;
        }
        else if (word == 0) /* RLaB is running LITTLE. */
        {
          freadptr = (FREADER) fread;
        }
      }
      else
      {
        close_file_ds (filename);
        rerror ("read: cannot read VAX or Cray binary formats\n");
      }

      //
      // Now, figure out the matrix type, and precision.
      //
      T = type % 10;    /* Matrix Type. */
      P = ((type - T) % 100) / 10;  /* Matrix Precision. */

      if (T == 0 || T == SM)  /* Numeric or String matrix. */
      {
        if ((m = matrix_ReadB (rf->fileds_f, T, P, &name, &rtype)) == 0)
        {
          close_file_ds (filename);
          rerror ("readb: cannot read matrix");
        }

        if ((lnode = btree_FindNode (symtab, name)) != 0)
        {
          ent_Destroy (var_ent (lnode));
          ent = ent_Create ();
          ent_data (ent) = m;
          ent_SetType (ent, rtype);
          ent_IncRef (ent);
          listNode_AttachEnt (lnode, ent);
        }
        else
        {
          ent = ent_Create ();
          ent_data (ent) = m;
          ent_SetType (ent, rtype);
          install (symtab, name, ent);
        }
      }
      else if (T == SP)   /* Sparse Numeric Matrix. */
      {
        if ((m = msparse_ReadB (rf->fileds_f, T, P, &name, &rtype)) == 0)
        {
          close_file_ds (filename);
          rerror ("readb: cannot read sparse matrix");
        }

        if ((lnode = btree_FindNode (symtab, name)) != 0)
        {
          ent_Destroy (var_ent (lnode));
          ent = ent_Create ();
          ent_data (ent) = m;
          ent_SetType (ent, rtype);
          ent_IncRef (ent);
          listNode_AttachEnt (lnode, ent);
        }
        else
        {
          ent = ent_Create ();
          ent_data (ent) = m;
          ent_SetType (ent, rtype);
          install (symtab, name, ent);
        }
      }
      else if (T == BT)   /* RLaB List (Binary Tree). */
      {
        if ((btree = btree_ReadB (rf->fileds_f, &name)) == 0)
        {
          close_file_ds (filename);
          rerror ("readb: cannot read list");
        }

        if ((lnode = btree_FindNode (symtab, name)) != 0)
        {
        /*
          One can destroy list only if it is not protected
        */
          ent_Destroy (var_ent (lnode));
          ent = ent_Create ();
          ent_data (ent) = btree;
          ent_SetType (ent, BTREE);
          ent_IncRef (ent);
          listNode_AttachEnt (lnode, ent);
        }
        else
        {
          ent = ent_Create ();
          ent_data (ent) = btree;
          ent_SetType (ent, BTREE);
          install (symtab, name, ent);
        }
      }
      else
      {
        close_file_ds (filename);
        rerror ("readb: unknown data type ?");
      }
    }


    rent = ent_Create ();
    if (feof (rf->fileds_f))
      ent_data (rent) = mdr_CreateScalar (0.0);
    else
    {
      fprintf (stderr, "ERROR: abnormal exit, did not encounter EOF");
      ent_data (rent) = mdr_CreateScalar (1.0);
    }

    if (close_after_rw)
      rfile_Destroy(filename);

    ent_Clean (e1);

    ent_SetType (rent, MATRIX_DENSE_REAL);
    return (rent);
  }

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (2);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

/***********************
 * Supporting Routines *
 ***********************/

/**********************
 * Write out Matrices *
 * type:
 * REAL    0000            (MT)
 * INTEGER 0002            (MI) - rlabplus
 * COMPLEX 0000            (MT)
 * STRING  0005            (SM)
 * BTREE   0006            (BT)
 * SPARSE  0007            (SP)
 **********************/

/*
 * RlaB Numeric Matrices are written out with double precision.
 */

//mdr_WriteB (MDR * m, FILE * fn, char *name)

void
    mdr_WriteB (Rfile *rf, char *name, MDR * m)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;

    FOURB imagf, nlen, type;
    FOURB nrow, ncol;

    // Write header info
    // Figure out header info
    if (m->type == RLAB_TYPE_INT32)
      type = word + MI;
    else
      type = word + MT;
    imagf = 0;
    nrow = (FOURB) MNR (m);
    ncol = (FOURB) MNC (m);
    nlen = (FOURB) strlen (name) + 1;

    fwrite (&type, sizeof (FOURB), 1, fn);
    fwrite (&nrow, sizeof (FOURB), 1, fn);
    fwrite (&ncol, sizeof (FOURB), 1, fn);
    fwrite (&imagf, sizeof (FOURB), 1, fn);
    fwrite (&nlen, sizeof (FOURB), 1, fn);
    fwrite (name, sizeof (char), nlen, fn);

    if (m->type == RLAB_TYPE_INT32)
      fwrite ( MDIPTR(m), sizeof (int), MNR (m) * MNC (m), fn);
    else if (m->type == RLAB_TYPE_DOUBLE)
      fwrite ( MDRPTR(m), sizeof (double), MNR (m) * MNC (m), fn);

    if (ferror (fn))
      warning_1 ("writeb: error occurred during write of MDR!");
    if (feof (fn))
      warning_1 ("writeb: EOF occurred during write");
  }
  else if (rf->filetype == RFILE_H5)
  {
    if (h5_write_rlab_atomic_mdr (rf, name, m) < 0)
    {
      fprintf   (stderr, "writeb: (h5) object = %s\n", name);
      warning_1 ("writeb: (h5) error occurred during write");
    }
  }
}

//void mdc_WriteB (MDC * m, FILE * fn, char *name)
void
    mdc_WriteB (Rfile *rf, char *name, MDC * m)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;
    int i, size;
    FOURB imagf, nlen, type;
    FOURB nrow, ncol;

    /* Write header info */

    type = word + MT;
    imagf = 1;

    nrow = (FOURB) MNR (m);
    ncol = (FOURB) MNC (m);
    nlen = (FOURB) strlen (name) + 1;

    fwrite (&type, sizeof (FOURB), 1, fn);
    fwrite (&nrow, sizeof (FOURB), 1, fn);
    fwrite (&ncol, sizeof (FOURB), 1, fn);
    fwrite (&imagf, sizeof (FOURB), 1, fn);
    fwrite (&nlen, sizeof (FOURB), 1, fn);
    fwrite (name, sizeof (char), nlen, fn);

    /* We must write out the real and complex parts separately */
    size = MNR (m) * MNC (m);

    for (i = 0; i < size; i++)
    {
      fwrite (&MdcV0r(m,i), sizeof (double), 1, fn);
    }
    for (i = 0; i < size; i++)
    {
      fwrite (&MdcV0i(m,i), sizeof (double), 1, fn);
    }

    if (ferror (fn))
      warning_1 ("writeb: error occurred during write");
    if (feof (fn))
      warning_1 ("writeb: EOF occurred during write");
  }
  else if (rf->filetype == RFILE_H5)
  {
    if (h5_write_rlab_atomic_mdc (rf, name, m) < 0)
    {
      fprintf   (stderr, "writeb: (h5) object = %s\n", name);
      warning_1 ("writeb: (h5) error occurred during write");
    }
  }
}

//void mds_WriteB (MDS * m, FILE * fn, char *name)
void
    mds_WriteB (Rfile *rf, char *name, MDS * m)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;
    int i, size;
    FOURB imagf, nlen, type;
    FOURB nrow, ncol;

    /* Write header info */

    type = word + SM;
    imagf = 0;

    nrow = (FOURB) MNR (m);
    ncol = (FOURB) MNC (m);
    nlen = (FOURB) strlen (name) + 1;

    fwrite (&type, sizeof (FOURB), 1, fn);
    fwrite (&nrow, sizeof (FOURB), 1, fn);
    fwrite (&ncol, sizeof (FOURB), 1, fn);
    fwrite (&imagf, sizeof (FOURB), 1, fn);
    fwrite (&nlen, sizeof (FOURB), 1, fn);
    fwrite (name, sizeof (char), nlen, fn);

    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      nlen = (FOURB) strlen (MdsV0(m,i)) + 1;
      fwrite (&nlen, sizeof (FOURB), 1, fn);
      fwrite (MdsV0(m,i), sizeof (char), nlen, fn);
    }

    if (ferror (fn))
      warning_1 ("writeb: error occurred during write");
    if (feof (fn))
      warning_1 ("writeb: EOF occurred during write");
  }
  else if (rf->filetype == RFILE_H5)
  {
    if (h5_write_rlab_atomic_mds (rf, name, m) < 0)
    {
      fprintf   (stderr, "writeb: (h5) object = %s\n", name);
      warning_1 ("writeb: (h5) error occurred during write");
    }
  }
}

/*
 * Binary write for matrix-sparse-real.
 */

//void msr_WriteB (MSR * m, FILE * fn, char *name)
void
    msr_WriteB (Rfile *rf, char *name, MSR * m)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;
    FOURB imagf, nlen, type;
    FOURB nrow, ncol, nnz;

    /* Gather information. */
    nrow = (FOURB) m->nr;
    ncol = (FOURB) m->nc;
    nnz = (FOURB) m->nnz;
    nlen = (FOURB) strlen (name) + 1;

    type = word + SP;
    imagf = 0;

    /* Write the header line. */
    fwrite (&type, sizeof (FOURB), 1, fn);
    fwrite (&nrow, sizeof (FOURB), 1, fn);
    fwrite (&ncol, sizeof (FOURB), 1, fn);
    fwrite (&nnz, sizeof (FOURB), 1, fn);
    fwrite (&imagf, sizeof (FOURB), 1, fn);
    fwrite (&nlen, sizeof (FOURB), 1, fn);
    fwrite (name, sizeof (char), nlen, fn);

    /* Write out IA. */
    fwrite (m->ia, sizeof (int), nrow + 1, fn);

    /* Write out JA, and D. */
    fwrite (m->ja, sizeof (int), nnz, fn);
    fwrite (m->d, sizeof (double), nnz, fn);

    if (ferror (fn))
      warning_1 ("writeb: error occurred during write");
    if (feof (fn))
      warning_1 ("writeb: EOF occurred during write");
  }
  else if (rf->filetype == RFILE_H5)
  {
    if (h5_write_rlab_atomic_msr (rf, name, m) < 0)
    {
      fprintf   (stderr, "writeb: (h5) object = %s\n", name);
      warning_1 ("writeb: (h5) error occurred during write");
    }
  }
}

/*
 * Binary write for matrix-sparse-complex.
 */

//void msc_WriteB (MSC * m, FILE * fn, char *name)
void
    msc_WriteB (Rfile *rf, char *name, MSC * m)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;
    FOURB imagf, nlen, type;
    FOURB nrow, ncol, nnz;

    /* Gather information. */
    nrow = (FOURB) m->nr;
    ncol = (FOURB) m->nc;
    nnz = (FOURB) m->nnz;
    nlen = (FOURB) strlen (name) + 1;

    type = word + SP;
    imagf = 1;

    /* Write the header line. */
    fwrite (&type, sizeof (FOURB), 1, fn);
    fwrite (&nrow, sizeof (FOURB), 1, fn);
    fwrite (&ncol, sizeof (FOURB), 1, fn);
    fwrite (&nnz, sizeof (FOURB), 1, fn);
    fwrite (&imagf, sizeof (FOURB), 1, fn);
    fwrite (&nlen, sizeof (FOURB), 1, fn);
    fwrite (name, sizeof (char), nlen, fn);

    /* Write out IA. */
    fwrite (m->ia, sizeof (int), nrow + 1, fn);

    /* Write out JA, and D. */
    fwrite (m->ja, sizeof (int), nnz, fn);
    fwrite (m->c, sizeof (Complex), nnz, fn);

    if (ferror (fn))
      warning_1 ("writeb: error occurred during write");
    if (feof (fn))
      warning_1 ("writeb: EOF occurred during write");
  }
  else if (rf->filetype == RFILE_H5)
  {
    if (h5_write_rlab_atomic_msc (rf, name, m) < 0)
    {
      fprintf   (stderr, "writeb: (h5) object = %s\n", name);
      warning_1 ("writeb: (h5) error occurred during write");
    }
  }
}

/* **************************************************************
 * Read a matrix. Don't bother trying to make this fit the
 * class interface due to Matlab's "gross" binary file format.
 * ************************************************************** */

void *
matrix_ReadB (FILE * fn, int T, int P, char **name, int *rtype)
{
  int i;
  FOURB nlen, imagf, nrow, ncol;
  void *m = 0;

  /* Read rest of header info */
  FREAD (&nrow, sizeof (FOURB), 1, fn);
  FREAD (&ncol, sizeof (FOURB), 1, fn);
  FREAD (&imagf, sizeof (FOURB), 1, fn);
  FREAD (&nlen, sizeof (FOURB), 1, fn);

  *name = (char *) GC_MALLOC (sizeof (char) * nlen);
  FREAD (*name, sizeof (char), nlen, fn);

  /* Now create matrix */

  /*
   * RLaB or Matlab Numeric Matrix.
   */

  //fprintf(stderr, "T = %i, P = %i\n", T, P);

  if (T == 0)
  {
    if (P == 0)     /* double */
    {
      if (imagf == 0)
      {
        m = (void *) mdr_Create (nrow, ncol);
        FREAD (MDRPTR(m), sizeof (double), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
      }
      else if (imagf == 1)
      {
        m = (void *) mdc_Create (nrow, ncol);
        for (i = 0; i < nrow * ncol; i++)
          FREAD (&MdcV0r((MDC *) m,i), sizeof (double), 1, fn);
        for (i = 0; i < nrow * ncol; i++)
          FREAD (&MdcV0i((MDC *) m,i), sizeof (double), 1, fn);
        *rtype = MATRIX_DENSE_COMPLEX;
      }
    }
    else if (P == 1)    /* float */
    {
      if (imagf == 0)
      {
        float *ftmp = (float *) GC_MAIOP ((nrow * ncol) * sizeof (float));
        FREAD (ftmp, sizeof (float), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
        m = (void *) mdr_Create (nrow, ncol);
        for (i = 0; i < nrow * ncol; i++)
          MdrV0 (m, i) = (double) ftmp[i];
        GC_FREE (ftmp);
      }
      else if (imagf == 1)
      {
        float *frtmp = (float *) GC_MAIOP ((nrow * ncol) * sizeof (float));
        float *fitmp = (float *) GC_MAIOP ((nrow * ncol) * sizeof (float));

        FREAD (frtmp, sizeof (float), 1, fn);
        FREAD (fitmp, sizeof (float), 1, fn);

        m = (void *) mdc_Create (nrow, ncol);
        *rtype = MATRIX_DENSE_COMPLEX;
        for (i = 0; i < nrow * ncol; i++)
        {
          MdcV0r ((MDC *) m, i) = (double) frtmp[i];
          MdcV0i ((MDC *) m, i) = (double) fitmp[i];
        }
        GC_FREE (frtmp);
        GC_FREE (fitmp);
      }
    }
    else if (P == 2)    /* 32 bit signed int */
    {
      if (imagf == 0)
      {
        // integer matrix
        m = (void *) mdi_Create (nrow, ncol);
        FREAD ( MDIPTR(m), sizeof (int), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
        //int *itmp = (int *) GC_MAIOP ((nrow * ncol) * sizeof (int));
        //FREAD (itmp, sizeof (int), nrow * ncol, fn);
        //*rtype = MATRIX_DENSE_REAL;
        //m = (void *) mdr_Create (nrow, ncol);
        //for (i = 0; i < nrow * ncol; i++)
        //  MdrV0 (m, i) = (double) itmp[i];
        //GC_FREE (itmp);
      }
      else if (imagf == 1)
      {
        int *irtmp = (int *) GC_MAIOP ((nrow * ncol) * sizeof (int));
        int *iitmp = (int *) GC_MAIOP ((nrow * ncol) * sizeof (int));
        FREAD (irtmp, sizeof (int), 1, fn);
        FREAD (iitmp, sizeof (int), 1, fn);

	m = (void *) mdc_Create (nrow, ncol);
	*rtype = MATRIX_DENSE_COMPLEX;
	for (i = 0; i < nrow * ncol; i++)
	{
    MdcV0r ((MDC *) m, i) = (double) irtmp[i];
    MdcV0i ((MDC *) m, i) = (double) iitmp[i];
	}
	GC_FREE (irtmp);
	GC_FREE (iitmp);
      }
    }
    else if (P == 3)		/* 16 bit signed int */
    {
      if (imagf == 0)
      {
        short int *itmp =
            (short int *) GC_MAIOP ((nrow * ncol) * sizeof (short int));
        FREAD (itmp, sizeof (short int), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
        m = (void *) mdi_Create (nrow, ncol);
        for (i = 0; i < nrow * ncol; i++)
          MdiV0 (m, i) = (int) itmp[i]; // integer
        GC_FREE (itmp);
      }
      else if (imagf == 1)
      {
	short int *irtmp =
	  (short int *) GC_MAIOP ((nrow * ncol) * sizeof (short int));
	short int *iitmp =
	  (short int *) GC_MAIOP ((nrow * ncol) * sizeof (short int));

	FREAD (irtmp, sizeof (short int), 1, fn);
	FREAD (iitmp, sizeof (short int), 1, fn);

	m = (void *) mdc_Create (nrow, ncol);
	*rtype = MATRIX_DENSE_COMPLEX;
	for (i = 0; i < nrow * ncol; i++)
	{
    MdcV0r ((MDC *) m, i) = (double) irtmp[i];
    MdcV0i ((MDC *) m, i) = (double) iitmp[i];
	}
	GC_FREE (irtmp);
	GC_FREE (iitmp);
      }
    }
    else if (P == 4)		/* 16 bit unsigned int */
    {
      if (imagf == 0)
      {
        short unsigned int *itmp = (short unsigned int *)
            GC_MAIOP ((nrow * ncol) * sizeof (short unsigned int));
        FREAD (itmp, sizeof (short unsigned int), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
        m = (void *) mdi_Create (nrow, ncol); // rlabplus
        for (i = 0; i < nrow * ncol; i++)
          MdiV0 (m, i) = (int) itmp[i];
        GC_FREE (itmp);
      }
      else if (imagf == 1)
      {
	short unsigned int *irtmp = (short unsigned int *)
	  GC_MAIOP ((nrow * ncol) * sizeof (short unsigned int));
	short unsigned int *iitmp = (short unsigned int *)
	  GC_MAIOP ((nrow * ncol) * sizeof (short unsigned int));

	FREAD (irtmp, sizeof (short unsigned int), 1, fn);
	FREAD (iitmp, sizeof (short unsigned int), 1, fn);

	m = (void *) mdc_Create (nrow, ncol);
	*rtype = MATRIX_DENSE_COMPLEX;
	for (i = 0; i < nrow * ncol; i++)
	{
    MdcV0r ((MDC *) m, i) = (double) irtmp[i];
    MdcV0i ((MDC *) m, i) = (double) iitmp[i];
	}
	GC_FREE (irtmp);
	GC_FREE (iitmp);
      }
    }
    else if (P == 5)		/* 8 bit unsigned int */
    {
      if (imagf == 0)
      {
        char *itmp = (char *) GC_MAIOP ((nrow * ncol) * sizeof (char));
        FREAD (itmp, sizeof (char), nrow * ncol, fn);
        *rtype = MATRIX_DENSE_REAL;
        m = (void *) mdi_Create (nrow, ncol); // rlabplus
        for (i = 0; i < nrow * ncol; i++)
          MdrV0 (m, i) = (int) itmp[i];
        GC_FREE (itmp);
      }
      else if (imagf == 1)
      {
	char *irtmp = (char *) GC_MAIOP ((nrow * ncol) * sizeof (char));
	char *iitmp = (char *) GC_MAIOP ((nrow * ncol) * sizeof (char));

	FREAD (irtmp, sizeof (char), 1, fn);
	FREAD (iitmp, sizeof (char), 1, fn);

	m = (void *) mdc_Create (nrow, ncol);
	*rtype = MATRIX_DENSE_COMPLEX;
	for (i = 0; i < nrow * ncol; i++)
	{
    MdcV0r ((MDC *) m, i) = (double) irtmp[i];
    MdcV0i ((MDC *) m, i) = (double) iitmp[i];
	}
	GC_FREE (irtmp);
	GC_FREE (iitmp);
      }
    }
  }

  /*
   * RLaB String Matrix.
   */

  else if (T == SM)
  {
    m = (void *) mds_Create (nrow, ncol);
    for (i = 0; i < nrow * ncol; i++)
    {
      FREAD (&nlen, sizeof (FOURB), 1, fn);
      MdsV0(m,i) = (char *) GC_MALLOC (sizeof (char) * nlen);
      FREAD ( MdsV0(m,i), sizeof (char), nlen, fn);
      *rtype = MATRIX_DENSE_STRING;
    }
  }

  /*
   * Unsupported type of Matlab Matrix.
   */

  else
  {
    rerror ("read: unsupported type of matlab matrix");
  }

  return (m);
}

/* **************************************************************
 * Read a sparse matrix.
 * ************************************************************** */

void *
msparse_ReadB (FILE * fn, int T, int P, char **name, int *rtype)
{
  FOURB nlen, imagf, nrow, ncol, nnz;
  void *m = 0;

  /* Read rest of header info */
  FREAD (&nrow, sizeof (FOURB), 1, fn);
  FREAD (&ncol, sizeof (FOURB), 1, fn);
  FREAD (&nnz, sizeof (FOURB), 1, fn);
  FREAD (&imagf, sizeof (FOURB), 1, fn);
  FREAD (&nlen, sizeof (FOURB), 1, fn);

  *name = (char *) GC_MALLOC (sizeof (char) * nlen);
  FREAD (*name, sizeof (char), nlen, fn);

  /* Now create matrix */
  if (!imagf)
  {
    /* Sparse-Real. */
    m = msr_Create (nrow, ncol);
    msr_Setup ((MSR *) m, nnz);
    *rtype = MATRIX_SPARSE_REAL;

    FREAD (((MSR *) (m))->ia, sizeof (int), nrow + 1, fn);
    FREAD (((MSR *) (m))->ja, sizeof (int), nnz, fn);
    FREAD (((MSR *) (m))->d, sizeof (double), nnz, fn);
  }
  else
  {
    /* Sparse-Complex. */
    /* Sparse-Real. */
    m = msc_Create (nrow, ncol);
    msc_Setup ((MSC *) m, nnz);
    *rtype = MATRIX_SPARSE_COMPLEX;

    FREAD (((MSC *) (m))->ia, sizeof (int), nrow + 1, fn);
    FREAD (((MSC *) (m))->ja, sizeof (int), nnz, fn);
    FREAD (((MSC *) (m))->c, sizeof (Complex), nnz, fn);
  }
  return (m);
}

/**************************
 * Write out Binary Trees *
 **************************/

static void
    btree_writeb_nodes (ListNode *, Rfile *);
static void
    btree_writeb_nodes_h5 (Rfile * rf, Btree * bt, char *name_prev);

/*
 * This subroutine writes the contents of a tree, but
 * does write information about the tree itself.
 */

//void btree_WriteB (Btree * btree, FILE * fn, char *name)
void
    btree_WriteB (Rfile *rf, char *name, Btree * btree)
{
  if (rf->filetype == RFILE_FILE)
  {
    FILE * fn = rf->fileds_f;
    // Check if the list is protected. Then do not write it to file.
    if(!btree->isconst)
    {
      FOURB num_nodes, type;
      FOURB nlen;

      type = word + BT;

      // Write out the header information
      fwrite (&type, sizeof (FOURB), 1, fn);

      // We must walk the tree first, to figure out
      // the number of valid nodes.
      num_nodes = btree_GetRealNumNodes (btree);
      fwrite (&num_nodes, sizeof (FOURB), 1, fn);

      // Tree Name
      nlen = strlen (name) + 1;
      fwrite (&nlen, sizeof (FOURB), 1, fn);
      fwrite (name, sizeof (char), nlen, fn);

      // Now write out contents
      btree_writeb_nodes (btree->root_node, rf);
    }
  }
  else  if (rf->filetype == RFILE_H5)
  {
    btree_writeb_nodes_h5 (rf, btree, name);
  }
}

static void
    btree_writeb_nodes_h5 (Rfile * rf, Btree *bt, char *name_prev)
{
  int i;
  void (*vfptr) ();

  // how many subnodes in the tree?
  int   nnode = btree_CountNodes (bt);
  if (!nnode)
    return;

  // get their names
  char **members = (char **) GC_MALLOC (nnode * sizeof (char *));
  btree_get_node_names (bt->root_node, members);

  // go over nodes, write the atomic_types out, and call itself for
  // the other btree's
  for (i=0; i< nnode; i++)
  {
    ListNode * node = btree_FindNode (bt, members[i]);
    Ent *e = var_ent (node);
    int type = ent_type (e);
    char * name = var_key (node);

    // name under which the variable is going to be saved is
    // "${name_prev}/${name}"
    char *obj_name = GC_malloc((strlen(name_prev)+strlen(name)+2)*sizeof(char));
    memcpy(obj_name, name_prev, strlen(name_prev));
    obj_name[strlen(name_prev)] = '/';
    memcpy(&obj_name[strlen(name_prev)+1], name, strlen(name)+1);

    // write the atomic types out, explore the btree's and ignore the rest
    if (    (type == MATRIX_DENSE_REAL)
        ||  (type == MATRIX_DENSE_COMPLEX)
        ||  (type == MATRIX_DENSE_STRING)
        ||  (type == MATRIX_SPARSE_REAL)
        ||  (type == MATRIX_SPARSE_COMPLEX)
       )
    {
      // prepare to write
      vfptr = (VVFPTR) writeb_method[type].op;
      if (vfptr == 0)
      {
        fprintf (stderr, "Entity type: %s\n", etd (e));
        rerror ("writeb() operation not supported");
      }

      // write it
      (*vfptr) (rf, obj_name, ent_data (e));

    }
    else if (type == BTREE)
    {
      // go further down for recursion
      btree_writeb_nodes_h5 (rf, ent_data(e), obj_name);
    }

    // clean-up
    if (obj_name)
      GC_free(obj_name);
    if (members[i])
      GC_free(members[i]);
  }

  // clean-up
  if (members)
    GC_free (members);
}

//static void btree_writeb_nodes (ListNode * node, FILE * fn)
static void
    btree_writeb_nodes (ListNode * node, Rfile * rf)
{
  char *name;
  int type;
  Ent *ent;
  void (*vfptr) ();

  if (node != 0)
  {
    btree_writeb_nodes (node->prev, rf);

    /*
     * Figure out how to write each node of the tree.
     * First we must check, and skip certain types
     * (like UNDEF, and FUNCTION).
     */

    ent = var_ent (node);
    type = ent_type (ent);

    if ((type == UNDEF) || (type == BLTIN) || (type == U_FUNCTION))
    {
      /* Skip on to the next one... */
    }
    else
    {
      /* Write the data (finally). */
      vfptr = (VVFPTR) writeb_method[type].op;

      if (vfptr == 0)
      {
        fprintf (stderr, "Entity type: %s\n", etd (ent));
        rerror ("writeb() operation not supported");
      }

      name = var_key (node);
      //(*vfptr) (ent_data (ent), fn, name);
      (*vfptr) (rf, name, ent_data (ent));
      if (!ferr_check ("writeb", rf->fileds_f))
        rerror ("writeb: error during write");
    }

    btree_writeb_nodes (node->next, rf);
  }
}

/**********************
 * Read a Binary Tree *
 **********************/

Btree *
btree_ReadB (FILE * fn, char **lname)
{
  int P, T, i, rtype;
  Ent *ent;
  FOURB num_nodes, nlen, type;
  Btree *newlist;
  Btree *btree;
  ListNode *lnode;
  char *name;
  void *m;

  FREAD (&num_nodes, sizeof (FOURB), 1, fn);

  FREAD (&nlen, sizeof (FOURB), 1, fn);
  *lname = (char *) GC_MALLOC (sizeof (char) * nlen);
  FREAD (*lname, sizeof (char), nlen, fn);

  /* Create a list and start reding object(s) from fn */
  newlist = btree_Create ();

  for (i = 0; i < num_nodes; i++)
  {
    FREAD (&type, sizeof (FOURB), 1, fn);

    T = type % 10;		/* Matrix Type. */
    P = ((type - T) % 100) / 10;	/* Matrix Precision. */

    if ( T == 0 || T == SM )
    {
      if ((m = matrix_ReadB (fn, T, P, &name, &rtype)) == 0)
        rerror ("readb: cannot read matrix");
      if ((lnode = btree_FindNode (newlist, name)) != 0)
      {
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        install (newlist, name, ent);
      }
    }
    else if (T == BT)
    {
      if ((btree = btree_ReadB (fn, &name)) == 0)
        rerror ("readb: cannot read list");
      if ((lnode = btree_FindNode (newlist, name)) != 0)
      {
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, BTREE);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, BTREE);
        install (newlist, name, ent);
      }
    }
    else
    {
      fprintf (stderr, "readb: unkown type: %s nodes: %i type: %i\n",
               name, (int) num_nodes, T);
      rerror ("readb: unknown type ?");
    }
  }
  return (newlist);
}

/***************************
 * Check for errors or EOF *
 ***************************/

static int
ferr_check (char *name, FILE * fn)
{
  if (ferror (fn))
  {
    warning_1 ("error occurred during I/O");
    return 0;
  }
  if (feof (fn))
  {
    warning_1 ("EOF occurred during I/O");
    return 0;
  }
  return 1;
}

/*
 * Return TRUE (1) if file-IO operations OK
 * Return FALSE (0) if EOF or file-IO error encountered.
 */

static int
fread_check (FILE * fn)
{
  if (ferror (fn))
    return 0;
  if (feof (fn))
    return 0;

  return 1;
}

/*************************************
 * Reverse the byte order of a word. *
 *************************************/

static void
reverse_word (void *word_ptr, int nbytes)
{
  char rword[32], *tmp;
  int i;

  memcpy (rword, word_ptr, nbytes);
  tmp = word_ptr;
  for (i = 0; i < nbytes; i++)
  {
    tmp[i] = rword[(nbytes - 1) - i];
  }

  return;
}

/*********************************************
 * A cover for fread() that swaps the bytes  *
 * to perform big/little endian translation. *
 * ONLY works for elemental data, i.e. ints  *
 * doubles, etc... not strucures.            *
 *********************************************/

size_t
fread_swap (void *ptr, size_t size, size_t nitems, FILE * stream)
{
  char *fptr, *tmp;
  size_t i;
  size_t stat;

  /* Point to beginning of array */
  tmp = (char *) ptr;

  for (i = 0; i < nitems; i++)
  {
    /*
     * Inc tmp as we go along reading and swapping
     * as we go.
     */

    fptr = tmp + i * size;
    stat = fread (fptr, size, 1, stream);

    if (stat != 1)		/* Early return */
      return (stat);

    reverse_word (fptr, size);
  }
  return (nitems);
}

/* **************************************************************
 * RLaB interface to fread()
 * var = fread ( FILE, nitems, type, swap_flag )
 *
 * var:    The return value, a matrix.
 *
 * FILE:   File descriptor string ("stdin", "stdout", etc...).
 * nitems: Number of items to read.
 * type:   "char", "short int", "int", "float", "double"
 * swap_flag: (optional) swap bytes if TRUE (1).
 * ************************************************************** */

#define FREADERR(ptr, size, n, fn) \
      if ((*freadptr) (ptr, size, n, fn) != n) \
        if (!fread_check (fn)) { rerror ("fread: error during read"); }

#include "mathl.h"

#define RLAB_FREAD_BUFFER_SIZE_ATOMIC 1024

Ent *
Fread (int nargs, Datum args[])
{
  double double_tmp;
  char *fnstring, *type;
  char *term_char=0;
  unsigned char have_term_char = 0;
  int count, nitems=-1, swapf, count_tc=0;
  int close_after_rw=0;
  Rfile *rf=0;
  FILE *fn=0;
  MDR *m=0;
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;

  /* Check n_args */
  if (nargs < 3)
    rerror ("fread: wrong number of input arguments");

  e1 = bltin_get_ent (args[0]);
  fnstring = class_char_pointer (e1);

  //
  // Try and open() the file
  //
  rf = rfile_find(fnstring);
  if (!rf)
  {
    rf = get_rfile_ds(fnstring, RFILE_FILE, "rb", 0);
    close_after_rw = 1;
  }
  else
  {
    // file exists. Is it regular file?
    if (rf->filetype != RFILE_FILE)
      rerror("fread: can read from regular files only!");
  }
  if (!rf)
    rerror("fread: Terrible internal error. Cannot open file for read");

  fn = rf->fileds_f;

  /*
   * Get the number of items/words we are supposed to read.
   */


  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) == MATRIX_DENSE_REAL)
  {
    double_tmp = class_double (e2);
    nitems = (int) class_double (e2);
    if (detect_inf_r (&double_tmp, 1))
      nitems = -1;
    else
    {
      if (nitems <= 0)
        rerror ("fread: nitems <= 0 invalid");
    }
  }

  /*
   * Get the type (int, float, etc...) of words to read.
   */

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_STRING)
    rerror ("fread: Cannot continue! Missing data type.");
  type = class_char_pointer (e3);

  /*
   * Check for byte-swapping request.
   */
  if (nargs == 4)
  {
    e4 = bltin_get_ent (args[3]);
    swapf = (int) class_double (e4);
  }
  else
  {
    swapf = 0;
  }

  if (nargs == 5)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type(e5) == MATRIX_DENSE_STRING)
    {
      term_char = class_char_pointer (e5);
      if (term_char)
        have_term_char = strlen(term_char);
    }
  }

  /*
   * Set the fread function pointer according to
   * byte-swapping needs.
   */

  if (swapf)
    freadptr = (FREADER) fread_swap;
  else
    freadptr = (FREADER) fread;

  if (!strcmp (type, "unsigned char") || !strcmp (type, "uchar") || !strcmp (type, "uint8_t"))
  {
    unsigned char char_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      count_tc=0;
      while (fread_check (fn))
      {
        if ((*freadptr)(&char_tmp, sizeof (unsigned char), 1, fn))
          MdiV0 (m, count++) = (int) char_tmp;
        else
          break;

        // attach termination characters to the output. do not gobble it
        if (have_term_char)
        {
          if (char_tmp==term_char[count_tc])
            count_tc++;
          if (count_tc==have_term_char)
            break;
        }

        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&char_tmp, sizeof (unsigned char), 1, fn))
          MdiV0 (m, count++) = (int) char_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "char") || !strcmp (type, "schar") || !strcmp (type, "int8_t"))
  {
    char char_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      count_tc=0;
      while (fread_check (fn))
      {
        if ((*freadptr)(&char_tmp, sizeof (char), 1, fn))
          MdiV0 (m, count++) = (int) char_tmp;
        else
          break;

        // attach termination characters to the output. do not gobble it
        if (have_term_char)
        {
          if (char_tmp==term_char[count_tc])
            count_tc++;
          if (count_tc==have_term_char)
            break;
        }
        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1,nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&char_tmp, sizeof (char), 1, fn))
          MdiV0 (m, count++) = (int) char_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "unsigned short int")  || !strcmp (type, "uint16_t") || !strcmp (type, "uint16"))
  {
    unsigned short int int_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      count_tc=0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&int_tmp, sizeof (unsigned short int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;

        // attach termination characters to the output. do not gobble it
        if (have_term_char)
        {
          if (int_tmp==term_char[count_tc])
            count_tc++;
          if (count_tc==have_term_char)
            break;
        }

        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&int_tmp, sizeof (unsigned short int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "short int")  || !strcmp (type, "int16_t") || !strcmp (type, "int16"))
  {
    short int int_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&int_tmp, sizeof (short int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;

        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&int_tmp, sizeof (short int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "unsigned int")  || !strcmp (type, "uint32_t") || !strcmp (type, "uint32"))
  {
    unsigned int int_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      count_tc=0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&int_tmp, sizeof (unsigned int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&int_tmp, sizeof (unsigned int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "int")  || !strcmp (type, "int32_t") || !strcmp (type, "int32"))
  {
    int int_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdi_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&int_tmp, sizeof (int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdi_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&int_tmp, sizeof (int), 1, fn))
          MdiV0 (m, count++) = (int) int_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "float") || !strcmp (type, "float32"))
  {
    float real_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdr_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&real_tmp, sizeof (float), 1, fn))
          MdrV0 (m, count++) = (double) real_tmp;
        else
          break;

        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdr_Create (1, nitems);
      count = 0;
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&real_tmp, sizeof (float), 1, fn))
          MdrV0 (m, count++) = (double) real_tmp;
        else
          break;
      }

      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else if (!strcmp (type, "double") || !strcmp (type, "float64"))
  {
    double real_tmp;
    if (nitems == -1)
    {
      // read entire content in 1-byte chunks
      m = mdr_Create (1, RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      count = 0;
      while (fread_check (fn))
      {
        if ((*freadptr) (&real_tmp, sizeof (double), 1, fn))
          MdrV0 (m, count++) = (double) real_tmp;
        else
          break;

        if (count > (m->ncol - 1))
          mdr_Extend (m, 1, m->ncol + RLAB_FREAD_BUFFER_SIZE_ATOMIC);
      }
      mdr_Extend (m, 1, count);
    }
    else
    {
      // read 'nitems' in 1-byte chunks
      m = mdr_Create (1, nitems);
      count = 0;
      while (fread_check (fn) && count<nitems)
      {
        if((*freadptr) (&real_tmp, sizeof (double), 1, fn))
          MdrV0 (m, count++) = (double) real_tmp;
        else
          break;
      }
      // if operation bombed before getting 'nitems', return what you got
      if (count < nitems)
      {
        if (count > 0)
          mdr_Extend (m, 1, count);
        else
        {
          mdr_Destroy(m);
          m = mdr_Create(0,0);
        }
      }
    }
  }
  else
  {
    rfile_Destroy(fnstring);
    rerror ("fread: Invalid TYPE argument. Closing the file!");
  }

  if (close_after_rw)
    rfile_Destroy(fnstring);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Assign_Rlab_MDR(m);
}

/* **************************************************************
 * RLaB interface to fwrite()
 * var = fwrite ( FILE, type, data)
 *
 * var:    The return value, a matrix, indicating the number of
 *         successfull writes.
 *
 * FILE:   File descriptor string ("stdin", "stdout", etc...).
 * type:   "char", "short int", "int", "float", "double"
 * data:   The data to write.
 * ************************************************************** */

#define FWRITE(ptr, size, fn) \
      if ((fwrite (ptr, size, 1, fn) != 1)) \
      if (!fread_check (fn)) { rerror ("fwrite: error during write"); }

Ent *
Fwrite (int nargs, Datum args[])
{
  char *fnstring, *type;
  int i, data_size = 0;
  int close_after_rw=0;
  Rfile *rf=0;
  FILE *fn=0;
  Ent *e1=0, *e2=0, *e3=0, *rent;

  int int_tmp;
  unsigned int uint_tmp;
  short int sint_tmp;
  unsigned char uchar_tmp;
  char char_tmp, *char_ptr;
  float float_tmp;
  double double_tmp;


  /* Check n_args */
  if (nargs < 3)
    rerror ("fwrite: insufficient input argument (needs 3)");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("fwrite: First argument 'filename' must be string");
  fnstring = class_char_pointer (e1);

  //
  // Try and open() the file
  //
  rf = rfile_find( fnstring );
  if (!rf)
  {
    // file does not exist: create it as a regular file
    rf = get_rfile_ds(fnstring, RFILE_FILE, "wb", 0);
    close_after_rw = 1;
  }
  else
  {
    // file exists. Is it regular file?
    if (rf->filetype != RFILE_FILE)
      rerror("fwrite: can write to regular files only!");
  }

  if (!rf)
    rerror("fwrite: Terrible internal error. Cannot open file for write");

  fn = rf->fileds_f;

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_STRING)
    rerror("fwrite: Second argument 'type' must be string");
  type = class_char_pointer (e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) == MATRIX_DENSE_REAL)
    data_size = MNR (ent_data (e3)) * MNC (ent_data (e3));
  else if (ent_type (e3) == MATRIX_DENSE_STRING)
  {
    MDS *ms = ent_data (e3);
    data_size = (ms->nrow) * (ms->ncol);
  }
  else
    rerror ("fwrite: cannot write this class");

  if (ent_type (e3) == MATRIX_DENSE_REAL)
  {
    if (!strcmp (type, "char"))
    {
      for (i = 0; i < data_size; i++)
      {
        char_tmp = (char) MdrV0 (ent_data (e3), i);
        FWRITE (&char_tmp, sizeof (char), fn);
      }
    }
    else if (!strcmp (type, "unsigned char"))
    {
      for (i = 0; i < data_size; i++)
      {
        uchar_tmp = (unsigned char) MdrV0 (ent_data (e3), i);
        FWRITE (&uchar_tmp, sizeof (unsigned char), fn);
      }
    }
    else if (!strcmp (type, "short int"))
    {
      for (i = 0; i < data_size; i++)
      {
        sint_tmp = (short int) MdrV0 (ent_data (e3), i);
        FWRITE (&sint_tmp, sizeof (short int), fn);
      }
    }
    else if (!strcmp (type, "unsigned int"))
    {
      for (i = 0; i < data_size; i++)
      {
        uint_tmp = (unsigned int) MdrV0 (ent_data (e3), i);
        FWRITE (&uint_tmp, sizeof (unsigned int), fn);
      }
    }
    else if (!strcmp (type, "int"))
    {
      for (i = 0; i < data_size; i++)
      {
        int_tmp = (int) MdrV0 (ent_data (e3), i);
        FWRITE (&int_tmp, sizeof (int), fn);
      }
    }
    else if (!strcmp (type, "float"))
    {
      for (i = 0; i < data_size; i++)
      {
        float_tmp = (float) MdrV0 (ent_data (e3), i);
        FWRITE (&float_tmp, sizeof (float), fn);
      }
    }
    else if (!strcmp (type, "double"))
    {
      for (i = 0; i < data_size; i++)
      {
        double_tmp = (double) MdrV0 (ent_data (e3), i);
        FWRITE (&double_tmp, sizeof (double), fn);
      }
    }
  }
  else if (ent_type (e3) == MATRIX_DENSE_STRING)
  {
    if (!strcmp (type, "char"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (char) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "unsigned char"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (unsigned char) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "short int"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (short int) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "unsigned int"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (unsigned int) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "int"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (int) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "float"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (float) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
    else if (!strcmp (type, "double"))
    {
      for (i = 0; i < data_size; i++)
      {
        int slen = 0;
        char_ptr = cpstr (MdsV0 (ent_data (e3), i));
        slen = strlen (char_ptr);
        FWRITE (char_ptr, sizeof (double) * slen, fn);
        GC_FREE (char_ptr);
      }
    }
  }
  else
  {
    rerror ("fwrite: invalid DATA argument");
  }

  if (close_after_rw)
    rfile_Destroy(fnstring);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (data_size);
  ent_SetType (rent, MATRIX_DENSE_REAL);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return (rent);
}

/* **************************************************************
 * RLaB interface to fseek.
 *
 * fseek ( FILENAME, OFFSET, ORIGIN )
 * FILENAME:   string
 * OFFSET:     byte (character) offset from position defined
 *             by ORIGIN.
 * ORIGIN:     "SEEK_SET"  (beginning (default))
 *             "SEEK_CUR"  (current position)
 *             "SEEK_END"  (end of file)
 * ************************************************************** */

Ent *
Fseek (int nargs, Datum args[])
{
  char *fnstring, *origin;
  int offset;
  int seek_val = 0;
  Rfile *rf=0;
  FILE *fn;
  Ent *e1=0, *e2=0, *e3=0, *rent;

  /* Check nargs */
  if (nargs < 2 || nargs > 3)
    rerror ("fseek: 2 or 3 arguments allowed");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("fseek: First argument 'filename' must be string");
  fnstring = class_char_pointer (e1);

  //
  // Try and open() the file
  //
  rf = rfile_find(fnstring);
  if (!rf)
    rerror("fseek: File does not exist. Did you open it first?");
  else
  {
    // file exists. Is it regular file?
    if (rf->filetype != RFILE_FILE)
      rerror("fseek: Requires RFILE to be a regular file!");
  }

  fn = rf->fileds_f;

  e2 = bltin_get_ent (args[1]);
  offset = (int) class_double (e2);

  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    origin = class_char_pointer (e3);

    if (!strcmp (origin, "SEEK_SET"))
      seek_val = SEEK_SET;
    else if (!strcmp (origin, "SEEK_CUR"))
      seek_val = SEEK_CUR;
    else if (!strcmp (origin, "SEEK_END"))
      seek_val = SEEK_END;
    else
    {
      rerror ("fseek: invalid ORIGIN argument");
    }
  }
  else
  {
    seek_val = SEEK_SET;
  }

  fseek (fn, offset, seek_val);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

/*
 * ASCII write
 * The interface writes out object(s), using the original
 * ASCII write method.
 *
 * Each object it completely responsible for writing/reading itself..
 */

Ent *
WriteASCII (int nargs, Datum args[])
{
  Ent *e1=0, *e=0, *rent;
  ListNode *lnode;
  char *string, *name;
  int i, type;
  Rfile *rf=0;
  int close_after_rw=0;
  FILE *fn;
  void (*vfptr) ();

  if (nargs < 2)
    rerror ("write_ascii: requires at least 2 arguments");

  e1 = bltin_get_ent (args[0]);
  string = class_char_pointer (e1);
  //
  // Try and open() the file
  //
  rf = rfile_find(string);
  if (!rf)
  {
    rf = get_rfile_ds(string, RFILE_FILE, "w", 0);
    close_after_rw = 1;
  }
  if (!rf)
    rerror("write_ascii: Terrible internal error. Cannot open file");

  fn = rf->fileds_f;

  /* Loop over the number of arguments, writing them out... */
  for (i = 1; i < nargs; i++)
  {
    if (args[i].type != VAR)
      rerror ("write: arguments must be variables");

    lnode = (ListNode *) args[i].u.ptr;
    e = convert_datum_to_ent (args[i]);
    name = var_key (lnode);
    type = ent_type (e);

    /* Write the data (finally). */
    vfptr = (VVFPTR) write_method[type].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e));
      rerror ("write() operation not supported");
    }

    (*vfptr) (ent_data (e), fn, name);
    if (!ferr_check ("write", fn))
      rerror ("write: error during write");
  }

  if (close_after_rw)
    rfile_Destroy (string);

  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

#define N_WRITE   5		/* # of columns to write pre line */

/*
 * Write out a REAL matrix
 */

static void
mdr_WriteASCII (MDR * m, FILE * fn, char *name)
{
  int i, j, k, nrow, ncol, npri, rem, start;

  nrow = MdrNR (m);
  ncol = MdrNC (m);
  npri = MdrNC (m) / N_WRITE;
  rem = MdrNC (m) % N_WRITE;

  fprintf (fn, "# matrix : %s no_of_rows: %d no_of_cols: %d",
	   name, MNR (m), MNC (m));

  fprintf (fn, " REAL\n");

  start = 1;
  for (i = 0; i < npri; i++)
  {
    fprintf (fn, "# matrix columns %d thru %d\n",
	     N_WRITE * i + 1, (N_WRITE - 1) + (N_WRITE * i) + 1);

    for (k = 1; k <= nrow; k++)	/* print all rows */
    {
      for (j = start; j <= N_WRITE + start - 1; j++)
      {
	fprintf (fn, " %-24.16g", Mdr1 (m, k, j));
      }
      fprintf (fn, "\n");
    }
    start += N_WRITE;		/* inc our col position */
  }

  /* Now come back and write out the last few columns */
  if (!rem)
    return;
  fprintf (fn, "# matrix columns %d thru %d\n", MNC (m) - rem + 1, MNC (m));

  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      fprintf (fn, " %-24.16g", Mdr1 (m, k, i));
    }
    fprintf (fn, "\n");
  }
  fflush (fn);
  return;
}

/*
 * Write out a COMPLEX matrix
 */

static void
mdc_WriteASCII (MDC * m, FILE * fn, char *name)
{
  int i, j, k, nrow, ncol, npri, rem, start;

  nrow = MdcNR (m);
  ncol = MdcNC (m);
  npri = MdcNC (m) / N_WRITE;
  rem = MdcNC (m) % N_WRITE;

  fprintf (fn, "# matrix : %s no_of_rows: %d no_of_cols: %d COMPLEX\n",
	   name, MNR (m), MNC (m));

  start = 1;
  for (i = 0; i < npri; i++)
  {
    fprintf (fn, "# matrix columns %d thru %d\n",
	     N_WRITE * i + 1, (N_WRITE - 1) + (N_WRITE * i) + 1);

    for (k = 1; k <= nrow; k++)	/* print all rows */
    {
      for (j = start; j <= N_WRITE + start - 1; j++)
      {
	fprintf (fn, " %-24.16g %-24.16g", Mdc1r (m, k, j), Mdc1i (m, k, j));
      }
      fprintf (fn, "\n");
    }
    start += N_WRITE;		/* inc our col position */
  }

  /* Now come back and write out the last few columns */
  if (!rem)
    return;
  fprintf (fn, "# matrix columns %d thru %d\n", MNC (m) - rem + 1, MNC (m));

  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      fprintf (fn, " %-24.16g %-24.16g", Mdc1r (m, k, i), Mdc1i (m, k, i));
    }
    fprintf (fn, "\n");
  }
  fflush (fn);
  return;
}

/*
 * Write out a STRING matrix.
 * Just write it out in a single column.
 * Note, that we _have_ to write out the string length
 * as the first field, since rlab strings are pretty
 * versatile.
 */

static void
mds_WriteASCII (MDS * m, FILE * fn, char *name)
{
  int i, nrow, ncol;

  nrow = MNR (m);
  ncol = MNC (m);

  fprintf (fn, "# matrix : %s no_of_rows: %d no_of_cols: %d",
	   name, MNR (m), MNC (m));

  fprintf (fn, " STRING\n");

  for (i = 0; i < nrow * ncol; i++)
  {
    fprintf (fn, "%i ",  (int) strlen (MdsV0 (m, i)));
    fprintf (fn, "%s\n", MdsV0 (m, i));
  }

  fflush (fn);
  return;
}

/*
 * Write out a Binary Tree (list) in ASCII format.
 * This function looks a lot like the WriteASCII function
 * since binary trees can contain any other sort of object.
 */

static void btree_writeascii_nodes (ListNode * lnode, FILE * fn);

static void
btree_WriteASCII (Btree * btree, FILE * fn, char *name)
{
  int num_nodes;

  /* Write out the header information */
  num_nodes = btree_GetRealNumNodes (btree);

  fprintf (fn, "# list : %s no_of_elements: %d \n", name, num_nodes);

  /* Now write out contents */
  btree_writeascii_nodes (btree->root_node, fn);
}

static void
btree_writeascii_nodes (ListNode * lnode, FILE * fn)
{
  Ent *e;
  char *name;
  int type;
  void (*vfptr) ();

  if (lnode != 0)
  {
    /* Go down the left branch. */
    btree_writeascii_nodes (lnode->prev, fn);

    /* Write out the node... */
    e = listNode_GetEnt (lnode);
    type = ent_type (e);
    name = var_key (lnode);
    vfptr = (VVFPTR) write_method[type].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e));
      rerror ("write() operation not supported");
    }

    (*vfptr) (ent_data (e), fn, name);
    if (!ferr_check ("write", fn))
      rerror ("write: error during write");

    /* Go down the right branch. */
    btree_writeascii_nodes (lnode->next, fn);
  }
}

/* **************************************************************
 * Read RLaB objects from a file. The default behavior is to read
 * ALL of the objects in the file, and install them on the GLOBAL
 * list. If the file contains a LIST, it is also put on the GLOBAL
 * list.
 *
 *
 * ************************************************************** */

MDR *mdr_ReadASCII (FILE * fn, int nrow, int ncol);
MDC *mdc_ReadASCII (FILE * fn, int nrow, int ncol);
MDS *mds_ReadASCII (FILE * fn, int nrow, int ncol);
Btree *btree_ReadASCII (FILE * fn, int nel);
MDS *legacy_string_Read (FILE * fn, int length);

Ent *
ReadASCII (int nargs, Datum args[])
{
  Datum arg2;
  Ent *e1=0, *e2=0, *ent, *rent;
  ListNode *lnode;
  Btree *btree, *symtab;
  char *string, name[256], type[20], dtype[20];
  int length, nel, nrow, ncol, retv, rtype=0;
  Rfile *rf=0;
  int close_after_rw=0;
  FILE *fn=0;
  void *m=0;

  /* Check n_args */
  if (nargs > 2 || nargs <= 0)
    rerror ("read: requires 1 or 2 arguments");

  /* get file name from argument list, always the 1st argument */
  e1 = bltin_get_ent (args[0]);
  string = class_char_pointer (e1);

  //
  // Try and open() the file
  //
  rf = rfile_find(string);
  if (!rf)
  {
    rf = get_rfile_ds(string, RFILE_FILE, "r", 0);
    close_after_rw = 1;
  }
  if (!rf)
    rerror("fread: Terrible internal error. Cannot open file for read");

  fn = rf->fileds_f;


  if (nargs == 1)
  {
    /* Use the global symbol table */
    symtab = get_symtab_ptr ();
  }
  else				/* nargs == 2 */
  {
    /* Use a user-supplied list */
    arg2 = args[1];

    /* Check to see if variable exists, etc... */
    if (arg2.type != VAR)
      rerror ("readb: 2nd argument must be a variable");

    /* If the specified variable is already a list, then just
     * create (or overwrite) members. If it is not a list, then
     * destroy it, and create a new list.
     */

    e2 = var_ent (arg2.u.ptr);

    if (ent_type (e2) == BTREE)
    {
      /* Create a new member. */
      symtab = ent_data (e2);
    }
    else
    {
      /* Destroy it, and create a list. */
      ent_Destroy (e2);

      e2 = ent_Create ();
      symtab = btree_Create ();
      ent_data (e2) = symtab;
      ent_SetType (e2, BTREE);
      ent_IncRef (e2);
      listNode_AttachEnt (arg2.u.ptr, e2);
    }
  }

  /*
   * Save current position in file, we will need to go back to
   * this point when we hand the read task off to the object
   * libraries
   */
  while (fscanf (fn, "# %s :", type) != EOF)
  {
    if (!strcmp ("scalar", type)) /* Legacy */
    {
      rerror ("scalar ASCII read not supported yet");
    }
    else if (!strcmp ("matrix", type))
    {
      /* Read the matrix particulars. */
      retv = fscanf (fn, "%s no_of_rows: %d no_of_cols: %d %s\n",
                     name, &nrow, &ncol, dtype);
      if (retv == 0 || retv == EOF)
      {
        fprintf (stderr, "read: ERROR in ASCII matrix format");
        fprintf (stderr, "in file %s\n", string);
        rerror ("read: invalid ASCII input");
      }

      if (!strcmp (dtype, "REAL"))
      {
        rtype = MATRIX_DENSE_REAL;
        if ((m = (void *) mdr_ReadASCII (fn, nrow, ncol)) == 0)
        {
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }
      else if (!strcmp (dtype, "COMPLEX"))
      {
        rtype = MATRIX_DENSE_COMPLEX;
        if ((m = (void *) mdc_ReadASCII (fn, nrow, ncol)) == 0)
        {
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }
      else if (!strcmp (dtype, "STRING"))
      {
        rtype = MATRIX_DENSE_STRING;
        if ((m = (void *) mds_ReadASCII (fn, nrow, ncol)) == 0)
        {
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }

      /*
      * Now, put the new matrix in the proper list.
      * Either global, or variable (defined above).
      */

      if ((lnode = btree_FindNode (symtab, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        install (symtab, (name), ent);
      }
    }
    else if (!strcmp ("string", type))  /* Legacy */
    {
      /* Get string particulatrs */
      fscanf (fn, "%s Length : %d\n", name, &length);
      m = (void *) legacy_string_Read (fn, length);
      rtype = MATRIX_DENSE_STRING;

      if ((lnode = btree_FindNode (symtab, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        install (symtab, (name), ent);
      }
    }
    else if (!strcmp ("list", type))
    {
      fscanf (fn, "%s no_of_elements: %d \n", name, &nel);
      btree = btree_ReadASCII (fn, nel);
      rtype = BTREE;

      if ((lnode = btree_FindNode (symtab, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, rtype);
        install (symtab, (name), ent);
      }
    }
    else
    {
      rfile_Destroy (string);
      fprintf (stderr, "read_ascii: do not know how to read");
      fprintf (stderr, " %s type of ASCII data\n", type);
      rerror ("read_ascii: unknown data type");
    }
  }

  if (close_after_rw)
    rfile_Destroy (string);

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

/*
 * Functions to read particular matrix data types.
 */

#define NLINE 256

MDR *
    mdr_ReadASCII (FILE * fn, int nrow, int ncol)
{
  int i, j, k, npri, rem, start;
  char jnk_str[56];
  MDR *new = mdr_Create (nrow, ncol);

  npri = ncol / N_WRITE;
  rem = ncol % N_WRITE;

  start = 1;
  for (i = 0; i < npri; i++)
  {
    fgets (jnk_str, NLINE, fn);
    /* fscanf(fn, "# matrix columns %*s %*s %*s\n"); */
    for (k = 1; k <= nrow; k++) /* read all rows */
    {
      for (j = start; j <= N_WRITE + start - 1; j++)
      {
        fscanf (fn, " %le", &Mdr1 (new, k, j));
      }
      fscanf (fn, "\n");
    }
    start += N_WRITE;   /* inc our col position */
  }

  /* Now come back and read the last few columns */
  if (!rem)
  {
    return (new);
  }

  fgets (jnk_str, NLINE, fn);
  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      fscanf (fn, " %le", &Mdr1 (new, k, i));
    }
    fscanf (fn, "\n");
  }
  return (new);
}

MDC *
    mdc_ReadASCII (FILE * fn, int nrow, int ncol)
{
  int i, j, k, npri, rem, start;
  char jnk_str[56];
  MDC *new = mdc_Create (nrow, ncol);

  npri = ncol / N_WRITE;
  rem = ncol % N_WRITE;

  start = 1;
  j = 0;
  for (i = 0; i < npri; i++)
  {
    fgets (jnk_str, NLINE, fn);
    for (k = 1; k <= nrow; k++) /* read all rows */
    {
      for (j = start; j <= N_WRITE + start - 1; j++)
      {
        fscanf (fn, " %le %le", &Mdc1r (new, k, j), &Mdc1i (new, k, j));
      }
      fscanf (fn, "\n");
    }
    start += N_WRITE;   /* inc our col position */
  }

  /* Now come back and read the last few columns */
  if (!rem)
  {
    return (new);
  }

  fgets (jnk_str, NLINE, fn);
  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      fscanf (fn, " %le %le", &Mdc1r (new, k, j), &Mdc1i (new, k, j));
    }
    fscanf (fn, "\n");
  }
  return (new);
}

MDS *
    mds_ReadASCII (FILE * fn, int nrow, int ncol)
{
  char *stmp;
  int i = 0;
  int j = 0;
  int slen = 0;
  int size = nrow * ncol;
  MDS *new = mds_Create (nrow, ncol);

  /* Now start reading the string. */
  for (j = 0; j < size; j++)
  {
    fscanf (fn, "%i ", &slen);
    stmp = (char *) GC_MALLOC ((size_t) sizeof (char) * (slen + 1));

    for (i = 0; i < slen; i++)
    {
      stmp[i] = (char) fgetc (fn);
    }

    /* Tack on the NULL. */
    stmp[slen] = '\0';

    /* Now, put the string into the matrix. */
    MdsV0 (new, j) = stmp;

    /* Make sure we get the newline. */
    if (slen > 0)
      fscanf (fn, "\n");
  }
  return (new);
}

Btree *
    btree_ReadASCII (FILE * fn, int nel)
{
  int i, length, nrow, ncol, retv, rtype;
  char dtype[20], *name, type[20], *string;
  Btree *btree, *newlist;
  void *m;
  ListNode *lnode;
  Ent *ent;

  /* Initialize to shut up compiler. */
  rtype = 0;
  name = 0;
  m = 0;

  /* Create a list and start reding object(s) from fn */
  newlist = btree_Create ();

  for (i = 0; i < nel; i++)
  {
    fscanf (fn, "# %s :", type);

    if (!strcmp ("scalar", type)) /* Legacy */
    {
      rerror ("scalar ASCII read not supported yet");
    }
    else if (!strcmp ("matrix", type))
    {
      /* Read the matrix particulars. */
      retv = fscanf (fn, "%s no_of_rows: %d no_of_cols: %d %s\n",
                     name, &nrow, &ncol, dtype);
      if (retv == 0 || retv == EOF)
      {
        rerror ("read_ascii: invalid ASCII input");
      }

      if (!strcmp (dtype, "REAL"))
      {
        rtype = MATRIX_DENSE_REAL;
        if ((m = (void *) mdr_ReadASCII (fn, nrow, ncol)) == 0)
        {
          string = get_file_ds_name (fn);
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }
      else if (!strcmp (dtype, "COMPLEX"))
      {
        rtype = MATRIX_DENSE_COMPLEX;
        if ((m = (void *) mdc_ReadASCII (fn, nrow, ncol)) == 0)
        {
          string = get_file_ds_name (fn);
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }
      else if (!strcmp (dtype, "STRING"))
      {
        rtype = MATRIX_DENSE_STRING;
        if ((m = (void *) mds_ReadASCII (fn, nrow, ncol)) == 0)
        {
          string = get_file_ds_name (fn);
          close_file_ds (string);
          rerror ("read: cannot read matrix");
        }
      }

      /*
      * Now, put the new matrix in the proper list.
      * Either global, or variable (defined above).
      */

      if ((lnode = btree_FindNode (newlist, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        install (newlist, (name), ent);
      }
    }
    else if (!strcmp ("string", type))  /* Legacy */
    {
      /* Get string particulatrs */
      fscanf (fn, "%s Length : %d\n", name, &length);
      m = (void *) legacy_string_Read (fn, length);
      rtype = MATRIX_DENSE_STRING;

      if ((lnode = btree_FindNode (newlist, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = m;
        ent_SetType (ent, rtype);
        install (newlist, (name), ent);
      }
    }
    else if (!strcmp ("list", type))
    {
      fscanf (fn, "%s no_of_elements: %d \n", name, &nel);
      btree = btree_ReadASCII (fn, nel);
      rtype = BTREE;

      if ((lnode = btree_FindNode (newlist, name)) != 0)
      {
        /* Symbol already exists. */
        ent_Destroy (var_ent (lnode));
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, rtype);
        ent_IncRef (ent);
        listNode_AttachEnt (lnode, ent);
      }
      else
      {
        /* Must create new symbol/variable. */
        ent = ent_Create ();
        ent_data (ent) = btree;
        ent_SetType (ent, rtype);
        install (newlist, (name), ent);
      }
    }
    else
    {
      string = get_file_ds_name (fn);
      close_file_ds (string);
      fprintf (stderr, "read_ascii: do not know how to read");
      fprintf (stderr, " %s type of ASCII data\n", type);
      rerror ("read_ascii: unknown data type");
    }
  }
  return (newlist);
}

MDS *
    legacy_string_Read (FILE * fn, int length)
{
  int i;
  char *string;
  MDS *new = mds_Create (1, 1);

  string = (char *) GC_MALLOC ((size_t) sizeof (char) * (length + 1));

  for (i = 0; i < length; i++)
  {
    string[i] = (char) fgetc (fn);
  }
  string[length] = '\0';

  MdsV0 (new, 0) = string;

  if (length > 0)
    fscanf (fn, "\n");

  return (new);
}

Ent *
Open (int nargs, Datum args[])
{
  ListNode *node;

  char *mode=0, *name=0;
  int buffsize = 0;
  int istatus = 0;
  Ent *e1=0, *e2=0, *e3=0, *rent;
  Rfile *rf = 0;
  double dto = 1; // second, default time-out for socket

  int rval = RLAB_STATUS_FAILURE;

  if (nargs < 1)
    rerror ("open: requires at least 1 argument");

  // we force user to use URLs for the file names:
  // 1. filename
  // 2. serial port (e.g., serial:///dev/ttyS0)
  // 3. URL         (e.g., http://www.suse.de)
  // 4. tcp socket  (e.g., tcp://127.0.0.1:5555)
  // 5. hierarhical data format ver. 5 (e.g., h5://filename)
  e1 = bltin_get_ent (args[0]);
  name = class_char_pointer (e1);
  if (!name)
    rerror("open: unknown filename");
  if (!strlen(name))
    rerror("open: empty string cannot be filename");

  // once we have a name,we can see if the url has been opened before
  rf = rfile_find(name);
  if (rf)
    goto Open_clean_up;

//   fprintf(stderr, "name = %s\n", name);
  if (    (!strncmp(name, "tcp://", 6))
      ||  (!strncmp(name, "udp://", 6))
     )
  {
    int imode=1;
    int iprot=1;

    // we can choose between the two protocols:
    if (!strncmp(name, "tcp://", 6))
      iprot = 1;
    else if (!strncmp(name, "udp://", 6))
      iprot = 0;

    // is there a second argument list with options
    if (nargs == 2)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == BTREE)
      {
        //
        // user provided options argument: use them to instantiate socket
        //
        // mode: (l)isten or (c)onnect
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SOCKET_MODE);
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            mode = class_char_pointer( var_ent (node) );
            if (mode)
              if (strlen(mode))
            {
              if (mode[0] == 'l' || mode[0] == 'L')
                imode = 0;
              else if (mode[0] == 'c' || mode[0] == 'C')
                imode = 1;
            }
          }
        }
        // timeout: how long do we sit at the socket and wait
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SOCKET_TIMEOUT);
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            dto = class_double( var_ent (node) );
            dto = dto < 0 ? -dto : dto;
          }
        }
      }
    }
    rf = get_rfile_ds (name, RFILE_SOCKET, iprot, imode, dto);
  }
	else
#ifdef HAVE_HDF5_SO
  // is it a hierarhical data format ver. 5?
  // allow for two formats:
  // (1) general URL: {protocol}://{filename}, with protocol=h5,hdf5; or,
  // (2) plain filename with extension ".h5"
  if (    (!strncmp(name, "h5://", 5))
      ||  (!strncmp(name, "hdf5://", 7))
      ||  (strstr (name, ".h5"))
     )
  {
    //
    // get the access property of the file: assume that the user knows what she is doing
    //
    if (nargs < 2)
      rerror("open: (h5) Requires file access mode");

    e2 = bltin_get_ent (args[1]);
    mode = class_char_pointer (e2);
    if (isvalidstring(mode)<1)
      rerror("open: (h5) Unknown file access mode");
    if (    mode[0]!='a' && mode[0]!='A'
        &&  mode[0]!='r' && mode[0]!='R'
        &&  mode[0]!='w' && mode[0]!='W'
       )
      rerror ("open: (h5) Unknown file access mode. Supported are 'a', 'w', and 'r'");

    // get/create the file
    rf = get_rfile_ds (name, RFILE_H5, mode);
  }
  else
#endif
#ifdef HAVE_LIBCURL
  // is it curl port?
  // if so then the second argument is a list of options
  if (    (!strncmp(name, "http://", 7))
      ||  (!strncmp(name, "https://", 8))
      ||  (!strncmp(name, "ftp://", 6))
     )
  {
    //
    // initialize CURL with 'name' if it already does not exist
    //
    rf = get_rfile_ds(name, RFILE_CURL);

    //
    // process the parameters sent to curl in a list:
    //
    if (nargs == 2)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == BTREE)
      {


        //
        // user provided options argument: use them to instantiate
        // the curl port
        //
        // CURLOPT_URL
        node = btree_FindNode (ent_data (e2), "CURLOPT_URL");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * url = class_char_pointer( var_ent (node) );
            if (strlen(url))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_URL, (url));
          }
        }
        else
          istatus |= curl_easy_setopt(rf->curl, CURLOPT_URL, (name));
        // CURLOPT_RLABPLUS_DEBUG
        node = btree_FindNode (ent_data (e2), "CURLOPT_RLABPLUS_DEBUG");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            int idebug = (int) class_double( var_ent (node) );
            if (idebug)
              CURLOPT_RLABPLUS_DEBUG = 1;
          }
        }
        // CURLOPT_PROXY
        node = btree_FindNode (ent_data (e2), "CURLOPT_PROXY");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * proxy = class_char_pointer( var_ent (node) );
            if (strlen(proxy))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXY, (proxy));
          }
        }
        // CURLOPT_PROXYPORT
        node = btree_FindNode (ent_data (e2), "CURLOPT_PROXYPORT");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            long iport = (long) class_double( var_ent (node) );
            istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYPORT, iport);
          }
        }
        // CURLOPT_PROXYTYPE
        node = btree_FindNode (ent_data (e2), "CURLOPT_PROXYTYPE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * proxytype = class_char_pointer( var_ent (node) );
            if (!strcmp(proxytype, "http"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_HTTP);
//          else if (!strcmp(proxytype, "http1"))
//            curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_HTTP_1_0);
            else if (!strcmp(proxytype, "socks4a"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS4A);
            else if (!strcmp(proxytype, "socks4"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS4);
            else if (!strcmp(proxytype, "socks5"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS5);
            else if (!strcmp(proxytype, "socks5h"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS5_HOSTNAME);
            else
            {
              fprintf(stderr, "Only following proxy types are supported by curl: n");
              fprintf(stderr, "'socks4a' , 'socks4' , 'socks5' , and 'socks5h' \n");
              rerror ("Wrong proxy type.\n");
            }
          }
        }
        // CURLOPT_HTTPPROXYTUNNEL
        node = btree_FindNode (ent_data (e2), "CURLOPT_HTTPPROXYTUNNEL");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int itunnel = (int) class_double( var_ent (node) );
            if (itunnel)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPPROXYTUNNEL, 1);
          }
        }
        // CURLOPT_INTERFACE
        node = btree_FindNode (ent_data (e2), "CURLOPT_INTERFACE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * eth0 = class_char_pointer( var_ent (node) );
            if (strlen(eth0))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_INTERFACE, (eth0));
          }
        }
        // CURLOPT_PORT
        node = btree_FindNode (ent_data (e2), "CURLOPT_PORT");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            long iport = (long) class_double( var_ent (node) );
            istatus |= curl_easy_setopt(rf->curl, CURLOPT_PORT, iport);
          }
        }
        // CURLOPT_FOLLOWLOCATION
        node = btree_FindNode (ent_data (e2), "CURLOPT_FOLLOWLOCATION");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int ifolowloc = (int) class_double( var_ent (node) );
            if (ifolowloc)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_FOLLOWLOCATION, 1L);
          }
        }
        // CURLOPT_UNRESTRICTED_AUTH
        node = btree_FindNode (ent_data (e2), "CURLOPT_UNRESTRICTED_AUTH");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int isauth = (int) class_double( var_ent (node) );
            if (isauth)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_UNRESTRICTED_AUTH, 1L);
          }
        }
        // CURLOPT_MAXREDIRS
        node = btree_FindNode (ent_data (e2), "CURLOPT_MAXREDIRS");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            long maxred = (long) class_double( var_ent (node) );
            if (maxred > 0)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_MAXREDIRS, maxred);
          }
        }
        // CURLOPT_TCP_NODELAY
        node = btree_FindNode (ent_data (e2), "CURLOPT_TCP_NODELAY");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int itunnel = (int) class_double( var_ent (node) );
            if (itunnel)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_TCP_NODELAY, 1L);
          }
        }
        // CURLOPT_USERAGENT
        node = btree_FindNode (ent_data (e2), "CURLOPT_USERAGENT");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * ua = class_char_pointer( var_ent (node) );
            if (strlen(ua))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_USERAGENT, (ua));
          }
        }
        // CURLOPT_REFERER
        node = btree_FindNode (ent_data (e2), "CURLOPT_REFERER");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * referer = class_char_pointer( var_ent (node) );
            if (strlen(referer))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_REFERER, (referer));
          }
        }
        // CURLOPT_COOKIE
        node = btree_FindNode (ent_data (e2), "CURLOPT_COOKIE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * cookie = class_char_pointer( var_ent (node) );
            if (strlen(cookie))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIE, (cookie));
          }
        }
        // CURLOPT_COOKIEFILE
        node = btree_FindNode (ent_data (e2), "CURLOPT_COOKIEFILE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * cookiefile = class_char_pointer( var_ent (node) );
            if (strlen(cookiefile))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIEFILE, (cookiefile));
          }
        }
        // CURLOPT_COOKIEJAR
        node = btree_FindNode (ent_data (e2), "CURLOPT_COOKIEJAR");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * cookiejar = class_char_pointer( var_ent (node) );
            if (strlen(cookiejar))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIEJAR, (cookiejar));
          }
        }
        // CURLOPT_COOKIESESSION
        node = btree_FindNode (ent_data (e2), "CURLOPT_COOKIESESSION");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int icookies = (int) class_double( var_ent (node) );
            if (icookies)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIESESSION, 1L);
          }
        }
        // CURLOPT_HTTPGET
        node = btree_FindNode (ent_data (e2), "CURLOPT_HTTPGET");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int ihttpget = (int) class_double( var_ent (node) );
            if (ihttpget)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPGET, 1L);
            else
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPGET, 0L);
          }
        }
        //
        // PROTOCOL OPTIONS
        //
        // CURLOPT_CRLF
        node = btree_FindNode (ent_data (e2), "CURLOPT_CRLF");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int itt = (int) class_double( var_ent (node) );
            if (itt)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_CRLF, 1L);
            else
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_CRLF, 0L);
          }
        }
        // CURLOPT_POSTQUOTE
        node = btree_FindNode (ent_data (e2), "CURLOPT_POSTQUOTE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            MDS *pq = ent_data( var_ent (node) );
            int i, n = pq->nrow * pq->ncol;
            struct curl_slist *pqlist=0;

            for (i=0; i<n; i++)
              curl_slist_append(pqlist, MdsV0(pq,i));
            istatus |= curl_easy_setopt(rf->curl, CURLOPT_POSTQUOTE, pqlist);

            // free the curl list
            curl_slist_free_all (pqlist);
          }
        }
        //
        // SSL OPTIONS
        //
        // CURLOPT_SSL_VERIFYPEER
        node = btree_FindNode (ent_data (e2), "CURLOPT_SSL_VERIFYPEER");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int ivp = (int) class_double( var_ent (node) );
            if (ivp)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYPEER, 1L);
            else
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYPEER, 0L);
          }
        }
        // CURLOPT_SSL_VERIFYHOST
        node = btree_FindNode (ent_data (e2), "CURLOPT_SSL_VERIFYHOST");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            int ivh = (int) class_double( var_ent (node) );
            if (ivh)
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYHOST, 1L);
            else
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYHOST, 0L);
          }
        }
        //
        // NAMES and PASSWORDS OPTIONS (Authentication)
        //
        // CURLOPT_USERPWD
        node = btree_FindNode (ent_data (e2), "CURLOPT_USERPWD");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * userpwd = class_char_pointer( var_ent (node) );
            if (strlen(userpwd))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_USERPWD, (userpwd));
          }
        }
        // CURLOPT_PROXYUSERPWD
        node = btree_FindNode (ent_data (e2), "CURLOPT_PROXYUSERPWD");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * proxyuserpwd = class_char_pointer( var_ent (node) );
            if (strlen(proxyuserpwd))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYUSERPWD, (proxyuserpwd));
          }
        }
        // CURLOPT_HTTPAUTH
        node = btree_FindNode (ent_data (e2), "CURLOPT_HTTPAUTH");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            char * httpauth = class_char_pointer( var_ent (node) );
            if (!strcmp(httpauth, "basic"))
              curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
            else if (!strcmp(httpauth, "digest"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_DIGEST);
//           else if (!strcmp(httpauth, "digest_ie"))
//             curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_DIGEST_IE);
            else if (!strcmp(httpauth, "gss"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_GSSNEGOTIATE);
            else if (!strcmp(httpauth, "ntlm"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_NTLM);
            else if (!strcmp(httpauth, "anysafe"))
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_ANYSAFE);
            else
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_ANY);
          }
        }
        // CURLOPT_VERBOSE
        node = btree_FindNode (ent_data (e2), "CURLOPT_VERBOSE");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            long iverb = (long) class_double( var_ent (node) );
            istatus |= curl_easy_setopt(rf->curl, CURLOPT_VERBOSE, iverb);
          }
        }
        // CURLOPT_STDERR
        node = btree_FindNode (ent_data (e2), "CURLOPT_STDERR");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            char * stdfile = class_char_pointer( var_ent (node) );
            if (strlen(stdfile))
            {
              if(rf->curl_stderr_fd)
                fclose(rf->curl_stderr_fd);
              rf->curl_stderr_fd = fopen(stdfile, "a");
              istatus |= curl_easy_setopt(rf->curl, CURLOPT_STDERR, rf->curl_stderr_fd);
            }
          }
        }
        //  CURLOPT_TIMEOUT
        node = btree_FindNode (ent_data (e2), "CURLOPT_TIMEOUT");
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_REAL)
          {
            long itout = (long) class_double( var_ent (node) );
            istatus |= curl_easy_setopt(rf->curl, CURLOPT_TIMEOUT, itout);
          }
        }
      }
    }

  }
  else
#endif
  if (    (!strncmp(name, "serial://", 9))
      ||  ( strstr (name, "/dev/tty"))
     )
  {
    // open a serial port
    //    open(filename, <<data_format;baud_rate;flow_control;raw;debug>>)
    //
    int i1;
    int ibrate = 600, idatb = 8, istop = 1, ihard = 0, iraw = 1;
    int iparity = 0, hupcl=0;
    char *eol=0;

    RLABPLUS_SERIAL_DEBUG = 0;

    if (nargs == 2)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == BTREE)
      {
        // data_format = "{5,6,7,8}{N,O,E}{1,2}"
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_DPS);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_STRING)
          {
            char * cdp = class_char_pointer (var_ent(node));
            if (!cdp)
              rerror("open: '" RLAB_NAME_SERIAL_DPS "' must be string");
            if (strlen(cdp)!=3)
              rerror("open: '" RLAB_NAME_SERIAL_DPS "' cannot be determined");
            // databits
            switch (cdp[0])
            {
              case '5':
                idatb = 5;
                break;
              case '6':
                idatb = 6;
                break;
              case '7':
                idatb = 7;
                break;
              case '8':
                idatb = 8;
                break;
            }
            // parity
            switch (cdp[1])
            {
              case 'n':
              case 'N':
                iparity = 0;
                break;
              case 'O':
              case 'o':
                iparity = 1;
                break;
              case 'e':
              case 'E':
                iparity = 2;
                break;
            }
            // stop bits
            switch  (cdp[2])
            {
              case '1':
                istop = 1;
                break;
              case '2':
                istop = 2;
                break;
            }
          }
        }

        // baud
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_BAUD);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_REAL)
          {
            i1 = (int) class_double (var_ent(node));
            if (i1 > 0)
              ibrate = i1;
          }
        }
        // hupcl
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_HUPCL);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_REAL)
          {
            i1 = (int) class_double (var_ent(node));
            if (i1 > 0)
              hupcl = 1;
          }
        }

        // flow_control = {x,h}
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_FC);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_STRING)
          {
            char * flow = class_char_pointer (var_ent(node));
            if (flow)
            {
              if (!strlen(flow))
                rerror("open: cannot determine ''" RLAB_NAME_SERIAL_FC "''");

              // databits for flow control
              switch (flow[0])
              {
                case 'n':
                case 'N':
                  // none?
                  ihard = 0;
                  break;

                case 'x':
                case 'X':
                  // software control? (xon|xoff)
                  ihard = 1;
                  break;

                case 'H':
                case 'h':
                  // hardware control
                  ihard = 2;
                  break;
              }
            }
          }
        }

        // raw or not
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_RAW);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_REAL)
          {
            i1 = (int) class_double (var_ent(node));
            if (i1 == 0)
              iraw = 0;
          }
        }

        // eol character
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_EOL);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_STRING)
            eol = class_char_pointer (var_ent(node));
        }

        // write debug information as ...
        node = btree_FindNode (ent_data (e2), RLAB_NAME_SERIAL_DEBUG);
        if (node != 0)
        {
          if (ent_type(var_ent(node))== MATRIX_DENSE_REAL)
          {
            i1 = (int) class_double (var_ent(node));
            if (i1 >= 1 && i1 <= 3)
              RLABPLUS_SERIAL_DEBUG = i1;
          }
          else if (ent_type(var_ent(node))== MATRIX_DENSE_STRING)
          {
            char * dbg = class_char_pointer (var_ent(node));
            if (dbg)
              if (strlen(dbg)>=1)
              {
                switch (dbg[0])
                {
                  case 'h':
                  case 'H':
                    // hex
                    RLABPLUS_SERIAL_DEBUG = 1;
                    break;
                  case 'c':
                  case 'C':
                    // char
                    RLABPLUS_SERIAL_DEBUG = 2;
                    break;
                  case 'i':
                  case 'I':
                    // int
                    RLABPLUS_SERIAL_DEBUG = 3;
                    break;
                }
              }
          }
        }
      }
    }

    // we have all the information
    rf = get_rfile_ds(name, RFILE_COMM, ibrate, idatb, iparity, istop, ihard, iraw, hupcl);
    if (rf)
    {
      if (rf->eol)
        GC_free(rf->eol);
      if (isvalidstring(eol)>0)
        rf->eol = cpstr(eol);
      else
        rf->eol = 0;
    }
  }
  else
  {
    //
    // open a file or a pipe
    //    open(filename, mode, options)
    //
    if (nargs < 2)
      rerror("open: Second argument 'file access mode' required for regular files");

    e2 = bltin_get_ent (args[1]);
    mode = class_char_pointer (e2);
    if (!mode)
      rerror("open: unknown file access mode");
    if (!strlen(mode))
      rerror("open: unknown file access mode");
    if (    mode[0]!='a' && mode[0]!='A'
            &&  mode[0]!='r' && mode[0]!='R'
            &&  mode[0]!='w' && mode[0]!='W'
       )
      rerror ("open: Unknown file access mode. Supported are 'a', 'w', and 'r'");

    char *eol=0, *csp = 0;
    MDS *fmt=0;

    // are there any options?
    if (nargs == 3)
    {
      e3 = bltin_get_ent (args[2]);
      if (ent_type (e3) == BTREE)
      {
        // buffer_size
        node = btree_FindNode (ent_data (e3), RLAB_NAME_FILE_BS);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            buffsize = (int) class_double (var_ent(node));
            if (buffsize < 0)
              buffsize = 0;
          }
        }

        // eol
        node = btree_FindNode (ent_data (e3), RLAB_NAME_FILE_EOL);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
            eol = class_char_pointer (var_ent(node));
        }

        // csp
        node = btree_FindNode (ent_data (e3), RLAB_NAME_FILE_CSP);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
            csp = class_char_pointer (var_ent(node));
        }

        // fmt
        node = btree_FindNode (ent_data (e3), RLAB_NAME_FILE_FMT);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
            fmt = ent_data(var_ent(node));
        }
      }
    }

    // now open it
    rf = get_rfile_ds (name, RFILE_FILE, mode, buffsize);

    if (rf)
    {
      //
      // open succeeded:
      //
      if (eol)
      {
        if (rf->eol)
          GC_free (rf->eol);
        rf->eol = cpstr(eol);
      }
      else
        rf->eol = cpstr(RLABPLUS_FILE_EOL_UNIX);

      if (csp)
      {
        if (rf->csp)
          GC_free (rf->csp);
        rf->csp = cpstr(csp);
      }
      else
        rf->csp = cpstr(RLABPLUS_FILE_CSP);

      if (fmt)
      {
        if (rf->fmt)
          mds_Destroy (rf->fmt);
        rf->fmt = mds_Copy( fmt );
      }
      else
      {
        rf->fmt = mds_CreateScalar( RLABPLUS_FILE_FMT );
      }
    }
  }

  // if we now have 'rf' then Open was a success
  if (rf)
    rval = RLAB_STATUS_SUCCESS;

Open_clean_up:

  // clean up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (rval);
  ent_type(rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// close rfile when needed
//
Ent *
Close (int nargs, Datum args[])
{
  char *fname=0;
  Ent *e1=0, *rent=0;
  Rfile *rf=0;

  int rval = RLAB_STATUS_FAILURE;

  if (nargs != 1)
    rerror ("close: 1 argument allowed");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror("close: requires string 'filename'");

  fname = class_char_pointer (e1);
  if (!fname)
    goto Close_get_out;
  if (fname[0] == '\0')
    goto Close_get_out;

  // now check the name and decide how to close it
//   fprintf(stdout,"trying to find %s\n", fname);
  rf = rfile_find (fname);
  if (rf)
  {
//     fprintf(stdout,"found %s\n", fname);
    if (close_file_ds (fname))
      rval = RLAB_STATUS_SUCCESS;
  }

Close_get_out:

  // clean-up
  ent_Clean (e1);

  rent = ent_Create ();
  ent_data(rent) = mdi_CreateScalar(rval);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

//
// list open files on the system: ignore the first
// file RFILE_NULL, and the last one (just an empty structure)
//
Ent *
ListOpenFiles (int nargs, Datum args[])
{

  Ent *rent;
  MDS *w=0;
  Rfile *rf;

  int i=0;

  if (nrfile>1)
  {
    w=mds_Create(nrfile-1,2);          // Number of open files
    rf = rfile_list;
    while (rf)
    {
      if (rf->name)
      {
        // get the name of the file
        Mds0(w, i, 0) = cpstr(rf->name);

        // figure out its type
        if (rf->filetype == RFILE_NULL)
          Mds0(w, i, 1) = cpstr("null");
        else if (rf->filetype == RFILE_FILE)
          Mds0(w, i, 1) = cpstr("file");
        else if (rf->filetype == RFILE_COMM)
          Mds0(w, i, 1) = cpstr("serial");
        else if (rf->filetype == RFILE_SOCKET)
          Mds0(w, i, 1) = cpstr("socket");
#ifdef HAVE_LIBUCRL
        else if (rf->filetype == RFILE_CURL)
          Mds0(w, i, 1) = cpstr("curl");
#endif
#ifdef HAVE_HDF5_SO
        else if (rf->filetype == RFILE_H5)
          Mds0(w, i, 1) = cpstr("HDF5");
#endif
      }
      i++;
      rf = rf->next;
    }
  }
  else
    w=mds_Create(0,0);

  rent = ent_Create ();
  ent_data(rent) = w;
  ent_type(rent) = MATRIX_DENSE_STRING;
  return rent;
}

#ifdef HAVE_HDF5_SO
//
// allow user to check if something is an hdf5 file without opening it
//
Ent *
ent_hdf5_isfile (int nargs, Datum args[])
{
  // rlab
  Ent *e1=0, *rent;
  char *name=0;
  MDS *s1=0;

  int ival=0;

  //
  // get the name of hdf5 file
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror ("h5isfile: 'file name' must be a single string");
  s1 = (MDS *) ent_data (e1);
  if ((s1->nrow) * (s1->ncol) != 1)
    rerror ("h5isfile: 'file name' must be a single string");
  name = MdsV0(s1,0);
  if (!name)
    rerror ("h5isfile: empty filename");
  if (!strlen(name))
    rerror ("h5isfile: empty filename");

  // print warning messages?
  if(!RLAB_HDF5_DEBUG)
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  ival = H5Fis_hdf5(name);
  if (ival <=0)
    ival = 0;

  // cleanup rlab
  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar(ival);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);

}
#endif

//
// remove nonprintable characters including
//  '\n' '\r' from both string ends
//
void chomp2_string(char *s)
{
  int i,j,ilen;

  ilen = isvalidstring(s);
  if (ilen>0)
  {
    for (i = ilen-1; i >= 0; i--)
    {
      if (((int) s[i]) <= 32)
        s[i] = '\0';
      else
        break;
    }
    ilen=isvalidstring(s);
    j=0;
    for (i=0; i < ilen; i++)
    {
      if (((int) s[i]) <= 32)
        j++;
      else
        break;
    }
    if (j>0)
    {
      for (i=j; i<=ilen; i++ )
        s[i-j] = s[i];
    }
  }
  return;
}

//
// remove '\n' '\r' from string end
//
void chomp_string(char *s)
{
  int i,ilen;

  ilen = isvalidstring(s);
  if (ilen>0)
  {
    for (i = ilen-1; i >= 0; i--)
    {
      if (  s[i] == 10 ||  s[i] == 13 )
        s[i] = 0;
      else
        break;
    }
  }
  return;
}


// Ent *
// ent_string_getstring (int nargs, Datum args[])
// {
//   Ent *e1=0, *e2=0;
//   MDS *s=0;
//   char *f=0;
//   char *rval;
//   FILE *istream=0;
//   int i, j, k;
//   MDR * nline = 0;
//
//   if (nargs > 3)
//     rerror (THIS_SOLVER ": " RLAB_ERROR_AT_MOST_THREE_ARG_REQUIRED);
//
//   // read from stdin
//   if (nargs == 0)
//   {
//     string_buff[0] = '\0';
//     fseek(stdin, 0, SEEK_END); // Flush standard input
//     rval = fgets (string_buff, MAX_STRING_BUFF, stdin);
//
//     // strip non-characters from end of the string
//     chomp_string(string_buff);
//
//     // return what you've got
//     return ent_Create_Rlab_String(string_buff);
//   }
//
//   // file name
//   e1 = bltin_get_ent (args[0]);
//   if (ent_type(e1) != MATRIX_DENSE_STRING)
//     rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
//   f = class_char_pointer (e1);
//   if (!f)
//     rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
//   if (!strlen(f))
//     rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
//
//   // indices of lines that are going to be red
//   if (nargs == 2)
//   {
//     // get it
//     e2 = bltin_get_ent (args[1]);
//     if (ent_type(e2) == MATRIX_DENSE_REAL)
//     {
//       nline = ent_data (e2);
//       if (!(MNR(nline)*MNC(nline)))
//         nline = 0;
//       else
//       {
//         // is 'nline' sorted ?
//         for (i=0; i<MNR(nline)*MNC(nline)-1; i++)
//         {
//           if (mdrV0(nline,i) >= MdrV0(nline,i+1))
//           {
//             nline = 0;
//             break;
//           }
//         }
//       }
//     }
//   }
//
//   if (*f == '|')
//   {
//     //
//     // pipe can be read only once: very slow if many lines are given
//     //
//     f++;
//     istream = popen (f, "r");
//     if (!istream)
//       goto nothing_to_do;
// //     rerror (THIS_SOLVER ": First argument 'filename' (pipe) cannot be open !\n");
//
//     s = mds_Create(0,0);
//     while (!feof (istream))
//     {
//       if (!fgets (string_buff, MAX_STRING_BUFF, istream))
//         break;
//
//       // strip non-characters from end of the string
//       chomp_string(string_buff);
//
//       if (!strlen(string_buff))
//         break;
//
//       if (strlen(string_buff))
//         s = mds_Stack(s, mds_CreateScalar(string_buff));
//     }
//     pclose (istream);
//   }
//   else
//   {
//     //
//     // we read in from file two times: first time we count the rows,
//     // second time we read the file in.
//     //
//     istream = fopen (f, "r");
//     if (!istream)
//       goto nothing_to_do;
//
//     if (!nline)
//     {
//       j = 0;
//       while (!feof (istream))
//       {
//         if (!fgets (string_buff, MAX_STRING_BUFF, istream))
//           break;
//         j++;
//       }
//       rewind (istream);
//       s = mds_Create (j, 1);
//     }
//     else
//       s = mds_Create (MNR(nline)*MNC(nline), 1);
//
//     if (nline)
//     {
//       j = 0;  // j counts the lines in the file
//       k = 0;  // k counts the entries in the array 'nline'
//       while (!feof (istream))
//       {
//         if (!fgets (string_buff, MAX_STRING_BUFF, istream))
//           break;
//
//         j++;
//         if (j != mdiV0(nline,k))
//           continue;
//
//         // strip non-characters from end of the string
//         chomp_string(string_buff);
//
//         MdsV0 (s, k) = cpstr (string_buff);
//         k++;
//         if (k == MNR(nline)*MNC(nline))
//           break;
//       }
//     }
//     else
//     {
//       j = 0;
//       while (!feof (istream) && j < MNR(s)*MNC(s))
//       {
//         if (!fgets (string_buff, MAX_STRING_BUFF, istream))
//           break;
//
//         // strip non-characters from end of the string
//         chomp_string(string_buff);
//
//         MdsV0 (s, j) = cpstr (string_buff);
//         j++;
//       }
//     }
//     fclose (istream);
//   }
//
// nothing_to_do:
//
//   ent_Clean (e1);
//   ent_Clean (e2);
//
//   return ent_Assign_Rlab_MDS (s);
// }

//
// reads: read input from keyboard or from file as string vector
//
extern char string_buff[MAX_STRING_BUFF];
extern int processed_pattern_starts_string(char * str, MDS * s);
extern int mds_find_first_left_pattern_in_string(char * str, MDS * processed_pattern, int * iptr_ex);
extern char * remove_pattern_from_string(MDS * processed_pattern, char * str);
#undef THIS_SOLVER
#define THIS_SOLVER "reads"
Ent *
ReadS (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int close_after_rw = 0;
  char *join_csp=0, *name=0;
  MDS *s=0, *comment=0, *comment_pattern=0, *note=0, *note_pattern=0, *lstrip=0, *lstrip_pattern=0;
  MDS *start=0, *start_pattern=0, *stop=0, *stop_pattern=0, *grep=0, *grep_pattern=0;
  MDR *userows=0, *set_userows=0;
  int i,j,min_line_len=1,iskiprows=0, join_rows=1, join_rows_counter=0;
  int start_pattern_notfound=1, stop_pattern_notfound=1;
//   FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  Rfile *rf=0;

  // Check n_args
  if (nargs > 2)
    rerror (THIS_SOLVER ": zero to two arguments required!");

  // read from stdin
  if (nargs == 0)
  {
    string_buff[0] = '\0';
    fseek(stdin, 0, SEEK_END); // Flush standard input
    fgets (string_buff, MAX_STRING_BUFF, stdin);

    // strip non-characters from end of the string
    chomp_string(string_buff);

    // return what you've got
    return ent_Create_Rlab_String(string_buff);
  }

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": first argument must be a single string");

  name = class_char_pointer(e1);
  if (isvalidstring(name)<1)
  {
    fprintf (stderr,THIS_SOLVER ": Cannot open for read: Filename cannot be empty !\n");
    goto reads_exit;
  }

  // Check if file exists in our list.
  // If the file does not exist pursue default behavior:
  //    User is trying to read an existing file on local file system
  //    However, after reading it, remove it from our list of open files.
  rf = rfile_find(name);

  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after reading
    // just assume it is a pipe or a regular file
    rf = get_rfile_ds (name, RFILE_FILE, "r", 0);
    close_after_rw = 1;
  }
  if (!rf)
  {
    fprintf (stderr, "%s, cannot open for read\n", name);
    goto reads_exit;
  }
  
  // Try and open() the file
  if (!rf->fileds_f)
  {
    fprintf (stderr, "%s, cannot open for read\n", name);
    goto reads_exit;
  }
  
  // indices of lines that are going to be red
  if (nargs == 2)
  {
    // get it
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_REAL)
    {
      userows = ent_data(e2);
      if (SIZE(userows)>0)
        set_userows = mdr_VectorSet(userows);
      if (SIZE(set_userows)<1)
      {
        mdr_Destroy(set_userows);
        set_userows=0;
        userows=0;
      }
    }
    else if (ent_type (e2) == BTREE)
    {
      ListNode *
          node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_MIN_LINE_LEN);
      if (node)
      {
        min_line_len = class_int (var_ent (node));
        if (min_line_len < 1 )
          min_line_len = 1;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_JOINROWS);
      if (node)
      {
        join_rows = class_int (var_ent (node));
        if (join_rows < 1 )
          join_rows = 1;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_JOINCSP);
      if (node)
      {
        join_csp = (char *) class_char_pointer (var_ent (node));
        if (isvalidstring(join_csp)<1)
          join_csp = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_SKIPROWS);
      if (node)
      {
        iskiprows = class_int (var_ent (node));
        if (iskiprows <= 0 )
          iskiprows = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_USEROWS);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          userows = ent_data(var_ent (node));
          if (SIZE(userows)>0)
            set_userows = mdr_VectorSet(userows);
          if (SIZE(set_userows)<1)
          {
            mdr_Destroy(set_userows);
            set_userows=0;
            userows=0;
          }
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_GREP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          grep = ent_data(var_ent (node));
          if (SIZE(grep)<1)
            grep=0;

          // process strings to create patterns
          if (grep)
          {
            grep_pattern = mds_Create(1,SIZE(grep));
            j=0;
            for (i=0; i<SIZE(grep);i++)
            {
              if (isvalidstring(MdsV0(grep,i))>0)
                MdsV0(grep_pattern,j++) = process_rlab_string_pattern(MdsV0(grep,i));
            }
            if (j)
              mds_Extend(grep_pattern,1,j);
            else
            {
              mds_Destroy(grep_pattern);
              grep_pattern = 0;
            }
          }
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_COMMENT);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          comment = ent_data(var_ent (node));
          if (SIZE(comment)<1)
            comment=0;

          // process strings to create patterns
          if (comment)
          {
            comment_pattern = mds_Create(1,SIZE(comment));
            j=0;
            for (i=0; i<SIZE(comment);i++)
            {
              if (isvalidstring(MdsV0(comment,i))>0)
                MdsV0(comment_pattern,j++) = process_rlab_string_pattern(MdsV0(comment,i));
            }
            if (j)
              mds_Extend(comment_pattern,1,j);
            else
            {
              mds_Destroy(comment_pattern);
              comment_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_NOTE);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          note = ent_data(var_ent (node));
          if (SIZE(note)<1)
            note=0;

          // process strings to create patterns
          if (note)
          {
            note_pattern = mds_Create(1,SIZE(note));
            j=0;
            for (i=0; i<SIZE(note);i++)
            {
              if (isvalidstring(MdsV0(note,i))>0)
                MdsV0(note_pattern,j++) = process_rlab_string_pattern(MdsV0(note,i));
            }
            if (j)
              mds_Extend(note_pattern,1,j);
            else
            {
              mds_Destroy(note_pattern);
              note_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_LSTRIP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          lstrip = ent_data(var_ent (node));
          if (SIZE(lstrip)<1)
            lstrip=0;

          // process strings to create patterns
          if (lstrip)
          {
            lstrip_pattern = mds_Create(1,SIZE(lstrip));
            j=0;
            for (i=0; i<SIZE(lstrip);i++)
            {
              if (isvalidstring(MdsV0(lstrip,i))>0)
                MdsV0(lstrip_pattern,j++) = process_rlab_string_pattern(MdsV0(lstrip,i));
            }
            if (j)
              mds_Extend(lstrip_pattern,1,j);
            else
            {
              mds_Destroy(lstrip_pattern);
              lstrip_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_START);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          start = ent_data(var_ent (node));
          if (SIZE(start)<1)
            start=0;

          // process strings to create patterns
          if (start)
          {
            start_pattern = mds_Create(1,SIZE(start));
            j=0;
            for (i=0; i<SIZE(start);i++)
            {
              if (isvalidstring(MdsV0(start,i))>0)
                MdsV0(start_pattern,j++) = process_rlab_string_pattern(MdsV0(start,i));
            }
            if (j)
              mds_Extend(start_pattern,1,j);
            else
            {
              mds_Destroy(start_pattern);
              start_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_STOP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          stop = ent_data(var_ent (node));
          if (SIZE(stop)<1)
            stop=0;

          // process strings to create patterns
          if (stop)
          {
            stop_pattern = mds_Create(1,SIZE(stop));
            j=0;
            for (i=0; i<SIZE(stop);i++)
            {
              if (isvalidstring(MdsV0(stop,i))>0)
                MdsV0(stop_pattern,j++) = process_rlab_string_pattern(MdsV0(stop,i));
            }
            if (j)
              mds_Extend(stop_pattern,1,j);
            else
            {
              mds_Destroy(stop_pattern);
              stop_pattern = 0;
            }
          }
        }
      }
    }
  }

  int irow=0, idx_set_userows=0, len, len_clean, ifoundcomment=-1;
  char *joined_str=0, c;

  if (set_userows)
    s = mds_Create(SIZE(set_userows),1);
  else
    s = mds_Create(20,1);

  //
  // read file/pipe in chunks of 20 rows
  //
  irow=0;
  j=-1;
  while (!feof (rf->fileds_f))
  {
    string_buff[0] = '\0';

    if (!fgets (string_buff, MAX_STRING_BUFF, rf->fileds_f))
      break;
    j++;

//     fprintf(stderr, "string_buff = %s\n", string_buff);

    char * str = string_buff, *str_clean=0;
    len = isvalidstring(string_buff);

    // do we skip it?
    if (j<iskiprows)
      continue;

    if (len < min_line_len)
      continue;

    // do we grep it?
    if (grep_pattern)
    {
      if (mds_find_first_left_pattern_in_string(str, grep_pattern, NULL) < 0)
        continue;
    }

    // did user provide list of rows they wants processed?
    if (set_userows)
    {
      if (j==0)
      {
        while ((j>(mdiV0(set_userows,idx_set_userows)-1)) && (idx_set_userows < SIZE(set_userows)) )
          idx_set_userows++;
      }

      if (idx_set_userows == SIZE(set_userows))
        break;

      if (j!=(mdiV0(set_userows,idx_set_userows)-1))
        continue;

      idx_set_userows++;
    }

    // if user provided start pattern: then we look for it before we can
    // process the lines
    if (start_pattern && start_pattern_notfound==1)
    {
      if (processed_pattern_starts_string(str, start_pattern) < 1)
        continue;

      start_pattern_notfound=0;
      if (stop_pattern)
        stop_pattern_notfound=1;
      continue;
    }
    
    // if user provided stop pattern: then we stop processing the lines
    if (stop_pattern && stop_pattern_notfound==1)
    {
      if (processed_pattern_starts_string(str, stop_pattern) > 0)
      {
        if (!start_pattern)
          break;

        stop_pattern_notfound=0;
        start_pattern_notfound=1;
      }
    }

    if (!stop_pattern_notfound)
      continue;

    // strip non-characters from end of the string
    chomp_string(str);

    // did we implement comments?
    //    if so, temporarily replace the first character of comment with '\n'
    //    to prevent match processing past it
    if (comment_pattern)
    {
      ifoundcomment = mds_find_first_left_pattern_in_string(str,comment_pattern,NULL);
      if (ifoundcomment >= 0)
      {
          // we have found comment
        if (ifoundcomment < min_line_len)
        {
            //  what is left after the comment is removed is not long enough
            //  to warrant further processing
          continue;
        }

          // we have modified string; update its length
        c = str[ifoundcomment];
        str[ifoundcomment] = '\0';
        len = isvalidstring(str);
      }
    }

    // at this point we have a row
    if (join_rows > 1)
    {
      join_rows_counter++;
      if (join_rows_counter > 1)
        string_concat(&joined_str, join_csp);
      string_concat(&joined_str, str);
      if (join_rows_counter < join_rows)
        continue;

        // fix 'str' before replacing it with another pointer
      if (ifoundcomment > -1)
      {
        str[ifoundcomment] = c;
        ifoundcomment = -1;
      }
      str = joined_str;
      join_rows_counter=0;
    }

    //
    // processing depends whether the notes are provided
    //
    if (note_pattern)
    {
        // notes are removed from string before the strings are processed
      str_clean = remove_pattern_from_string(note_pattern, str);
      len_clean = isvalidstring(str_clean);
      if (len_clean < min_line_len)
      {
        GC_FREE (str_clean);
        continue;
      }

        // fix 'str' before replacing it with another pointer
      if (ifoundcomment > -1)
      {
        str[ifoundcomment] = c;
        ifoundcomment = -1;
      }

      str = str_clean;
      len = len_clean;
    }

    // finally: apply lstrip to remove anything superficial in the
    // string prior to processing
    int iblank = 0;
    if (lstrip_pattern)
    {
      iblank = processed_pattern_starts_string(str, lstrip_pattern);
    }

    irow++;
    if (irow > MNR(s))
      mds_Extend(s, MNR(s) + 20, 1);

    if (isvalidstring(str+iblank)>0)
      MdsV1(s,irow) = cpstr(str+iblank);
    else
      MdsV1(s,irow) = cpstr("");

    // we implemented comments
    //   revert string to its previous value but only if notes and join_lines have not affected
    //   its position
    if (ifoundcomment > -1)
    {
      str[ifoundcomment] = c;
      ifoundcomment = -1;
    }
    
    if (str_clean)
    {
      GC_FREE (str_clean);
      str_clean = 0;
    }
    
    if (joined_str)
    {
      GC_FREE(joined_str);
      joined_str = 0;
    }

  }
  if (irow>0)
    mds_Extend(s, irow, 1);
  else
  {
    mds_Destroy(s);
    s = 0;
  }

reads_exit:

  if (close_after_rw)
    rfile_Destroy(name);
  if (note_pattern)
    mds_Destroy(note_pattern);
  if (lstrip_pattern)
    mds_Destroy(lstrip_pattern);
  if(stop_pattern)
    mds_Destroy(stop_pattern);
  if(start_pattern)
    mds_Destroy(start_pattern);
  if (set_userows)
    mdr_Destroy(set_userows);

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDS (s);
}
