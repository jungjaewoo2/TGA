/* main.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992-1997  Ian R. Searle
   Copyright (C) 2004-2013  Marijan Kostrun

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

#include "version.h"
#include "rlab.h"
#include "mem.h"
#include "code.h"
#include "util.h"
#include "rfile.h"
#include "rfileio.h"
#include "bltin.h"
#include "class.h"
#include "bltin1.h"
#include "bltin2.h"
#include "bltin3.h"
#include "bltin4.h"
#include "bltin5.h"
#include "rfft.h"
#include "mdc.h"
#include "mathl.h"

// #include <termcap.h>

#include <stdio.h>
#include <sys/types.h>

#ifdef __STDC__
#include <stdlib.h>
#include <argp.h>
#else
extern char *getenv ();
#endif

#ifdef HAVE_READLINE
extern int read_history (const char *);
extern void stifle_history (int);
extern int write_history (const char *);
#endif

#undef  THIS_FILE
#define THIS_FILE "main.c"
#include "rlab_macros.h"
// #ifdef unix
// static char PATH_DELIM[] = "/";
// #endif

// #ifdef OS2
// static char PATH_DELIM[] = "\\";
// #endif

// #ifdef DOS
// static char PATH_DELIM[] = "\\";
// #endif

#ifdef WIN32
// static char PATH_DELIM[] = "\\";
char *DEFAULT_RC0 = "rlab.rc";
char *DEFAULT_HELP = "doc\\help";
char *DEFAULT_LIB = "rlib";
char *DEFAULT_PAGER = "more";
char *DEFAULT_HELP_PAGER = "more";
char *DEFAULT_SEARCH_PATH = "toolbox";
#endif /* WIN32 */

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#else
#if defined (WIN32)
#include <io.h>
#else	/* if defined (unix) */
#include <sys/dir.h>
#endif
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

/* Scanner declarations */
extern void set_rlab_input (int type);
extern int yyparse ();
extern FILE *yyin;
extern int lineno;
extern int flush_line;
extern char *curr_file_name;

/* print.c */
extern int close_file_ds (char *fn);

/* code.c */
extern void init_frame (void);
extern void init_stack (void);
extern void init_flstack (void);

/* getline.c: so that rlab script can be called with command line arguments */
char **rlab_argv=0;
int    rlab_argc=0;
int    rlab_argv_idx[50];
int    rlab_argv_count=0;

void init_misc (void);
void init_var (void);
void init_file_list (void);
void init_static_tree (void);
int run_program (char *input);
void run (Program * p);
void run_no_clean (Program * p);
void run_debug (Program * p);
void warning_1 (char *s);
void print_greeting (void);

void init_environment (void);
int  is_class_eval=0;


extern void diss_assemble (Inst * p, int progoff);

char *progname;

static char DEF_RLAB2_RC0[] = "RLAB2_RC0";

static char *rlab_rc0;		/* RLAB2_RC0, DEFAULT_RC0 */
static char *help_dir;		/* RLAB2_HELP_DIR, DEFAULT_HELP */
static char *lib_dir;		/* RLAB2_LIB_DIR, DEFAULT_LIB */
static char *pager;		/* RLAB2_PAGER, DEFAULT_PAGER */
static char *help_pager;	/* RLAB2_HELP_PAGER */
static char *search_path;	/* RLAB2_PATH, DEFAULT_SEARCH_PATH */
char *rl_histfile;		/* Readline history stuff. */
int rl_histsize;

static int print_machine = 0;   /* If TRUE, print contents of machine queue */
int use_readline = 1;           /* If FALSE DO NOT use GNU readline */
static int use_rc0 = 1;         /* If TRUE run rlab_rc0 on start-up */
static int use_pager = 1;       /* if TRUE use default pager (more) */
static int load_lib = 1;        /* if TRUE load libraries */
int line_nos = 1;               /* If TRUE use line #s in op-codes */
static int message = 1;         /* If TRUE print the start-up message */
int symtab_debug = 0;           /* If TRUE print out symbol-table debugging info on err */

/* If TRUE we will go interactive after start-up */
static int force_interactive = 0;
/* If TRUE we are running interactively */
int interactive = 1;

//
// argp command line interface to rlab
//
const char *argp_program_version = "rlab2-" version_string;
const char *argp_program_bug_address = "<mkostrun@gmail.com>";
static char doc[] = "rlab2-" version_string " -- matrix oriented, interactive programming environment";
static char args_doc[] = "[FILES] [- ARGS]";
static struct argp_option options[] = {
  {"debug",       'd',          0,      0,  "Print out disassembled internal programs" },
  {"interactive", 'i',          0,      0,  "If executing scripts stay interactive upon their completion" },
  {"nolibs",      'l',          0,      0,  "Do not load default libraries at the start" },
  {"nomsg",       'm',          0,      0,  "Do not print greeting" },
  {"norc0",       'q',          0,      0,  "Do not load initialization file" },
  {"exec",        'e', "\"CMDS\"",      0,  "Execute commands prior to entering rlab shell or loading rlab scripts" },
  { 0 }
};
/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  int   message;
  int   print_machine;
  int   force_interactive;
  int   load_lib;
  int   use_rc0;
  int   count;
  char *cmd;  // Argument for -e
};
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
  {
    case 'd':
      arguments->print_machine = 1;
      arguments->count++;
      break;
    case 'i':
      arguments->force_interactive=1;
      arguments->message = 0;
      arguments->count++;
      break;
    case 'l':
      arguments->load_lib = 0;
      arguments->count++;
      break;
    case 'm':
      arguments->message = 0;
      arguments->count++;
      break;
    case 'q':
      arguments->use_rc0 = 0;
      arguments->count++;
      break;
    case 'e':
      arguments->cmd = arg;
      arguments->count += 2;
      break;
    default:
      break;
  }
  return 0;
}
/*
   The ARGP structure itself.
 */
static struct argp argp = {options, parse_opt, args_doc, doc};


/*
static char usage_string[] = "rlab -icVdhlmnpq [-e \"cmds\"] [file(s)] [-]\n\
    -V        Print version information and exit\n\
    -d        Print out disassembled internal programs\n\
    -h        Print out this message\n\
    -i        if executing scripts stay interactive upon their completion\n\
    -l        Do not load default libraries\n\
    -m        Do not print greeting\n\
    -q        Do not load initialization file\n\
    -e \"cmds\" Execute commands in provided string first - no white spaces allowed";
 */

int run_string (char *evalue);

/* **************************************************************
 * Dummy MAIN__ to make linker/f2c happy.
 * ************************************************************** */

int MAIN__ (void)
{
  return (0);
}

/* **************************************************************
 * Another main trick to help the garbage collector.
 * ************************************************************** */

#ifdef HAVE_GC
#include "private/gc_priv.h"
int real_main (int argc, char *argv[]);
int main (int argc, char *argv[])
{
  int dummy;

  GC_stackbottom = (ptr_t) (&dummy);
  return (real_main (argc, argv));
}
#endif /* HAVE_GC */

/* **************************************************************
 * main, RLaB
 * ************************************************************** */

#ifdef HAVE_GC
int real_main (int argc, char *argv[])
#else
int main (int argc, char *argv[])
#endif /* HAVE_GC */
{
  int c, r;
  char *evalue=0;

  /* Process command line args */
  progname = argv[0];
  set_progname (progname);

  //
  // argp processing: init defaults
  // 
  struct arguments arguments;
  arguments.force_interactive = force_interactive;
  arguments.message = message;
  arguments.print_machine = print_machine;
  arguments.load_lib = load_lib;
  arguments.use_rc0 = use_rc0;
  arguments.cmd = 0;  // Argument for -e
  arguments.count = 0;

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  message = arguments.message;
  print_machine = arguments.print_machine;
  if (arguments.cmd)
    evalue = arguments.cmd;
  force_interactive  = arguments.force_interactive;
  load_lib = arguments.load_lib;
  use_rc0 = arguments.use_rc0;

//   fprintf(stderr, "argument.count = %i\n", arguments.count);

  /* Perform initialization */
  init_misc ();
  init_environment ();
  init_var ();
  init_symbol_table ();
  init_file_list ();
  init_static_tree ();

  class_init ();
  class_bltin1_init ();
  class_bltin2_init ();
  class_bltin3_init ();
  class_bltin4_init ();
  class_bltin5_init ();
  class_fft_init ();
  class_io_init ();

  if (use_rc0)
  {
    if ((r = run_program (rlab_rc0)) == 0)
    {
      fprintf (stderr, "\nCould not open RLaB init script\n");
      fprintf (stderr, "try setting environment variable RLAB2_RC0 ");
      fprintf (stderr, "and re-run RLaB\n\n\n");
    }
    else if (r < 0)
    {
      fprintf (stderr, "ERROR in rlab_rc0 file\n");
      return (0);
    }
  }

  /* Process rfiles in library */
  if (load_lib)
  {
    rfile_dir (lib_dir);
  }

  // did we finish the list of argp parameters with "--" ?
  c = arguments.count+1;
  if (argv[c])
    if (!strcmp(argv[c],"--"))
      c++;

  // after procession the start parameters using getopt() assume the names of rlab scripts
  // follow, terminated '-'
  while (c < argc)
  {
//     fprintf(stderr, "post: argv[%i] = %s\n", c, argv[c]);
    if (argv[c][0] == '-')
    {
      if (isvalidstring(argv[c])>1)
        c--;
      break;
    }
    else
    {
      rlab_argv_idx[rlab_argv_count++] = c;
      interactive = 0;
    }
    c++;
  }
//   fprintf(stderr, "post: rlab_argv_count = %i, interactive=%i, force_interactive=%i\n",rlab_argv_count, interactive, force_interactive);

  // what follows are the array of strings that is accessible to rlab through argv() command
  if ((++c) < argc)
  {
    rlab_argc = argc - c;
    rlab_argv = &argv[c];
  }

  // first, execute statements in 'evalue'
  if (evalue)
  {
    run_string (evalue);
  }

  if (force_interactive)
    interactive=1;

  // now run the scripts if provided
  if (rlab_argv_count>0)
  {
    for (c=0; c<rlab_argv_count; c++)
      run_program (argv[rlab_argv_idx[c]]);
  }

  /* Finally, go interactive */
  if (interactive)
  {
    if (message)
      print_greeting ();
    run_program (0);

#ifdef HAVE_READLINE
    if (rl_histfile)
    {
      stifle_history (rl_histsize);
      write_history (rl_histfile);
    }
#endif /* HAVE_READLINE */
  }

  return (0);
}

/* **************************************************************
 * Misc initialization.
 * ************************************************************** */

#include "fpe.h"

void
init_misc ()
{
  /*
   * Set the garbage collector for
   * incremental/generational operations.
   */

#ifdef HAVE_GC
  /* GC_enable_incremental (); */
#endif

  /* Set the input function pointer */
  set_rlab_input (0);

  /* Initialize the interpreter stacks */
  init_frame ();
  init_stack ();
  init_flstack ();

  /* Try and setup for fpe exceptions. */
  setup_fpe_handling ();

  /* Create Inf and NaNs. */
  init_inf_nan ();
}

 /* **************************************************************
  * Get all RLaB environment/default variables.
  * ************************************************************** */
#if defined(unix)
void
init_environment ()
{
  if ((rlab_rc0 = getenv (DEF_RLAB2_RC0)) == 0)
    rlab_rc0 = DEFAULT_RC0;

  if ((help_dir = getenv ("RLAB2_HELP_DIR")) == 0)
    help_dir = DEFAULT_HELP;

  if ((lib_dir = getenv ("RLAB2_LIB_DIR")) == 0)
    lib_dir = DEFAULT_LIB;

  if (use_pager)
  {
    if ((pager = getenv ("RLAB2_PAGER")) == 0)
    {
      if ((pager = getenv ("PAGER")) == 0)
	pager = DEFAULT_PAGER;
    }
    if ((help_pager = getenv ("RLAB2_HELP_PAGER")) == 0)
    {
      help_pager = DEFAULT_HELP_PAGER;
    }
  }
  else
  {
    pager = cpstr ("cat");
    help_pager = cpstr ("cat");
  }

  if ((search_path = getenv ("RLAB2_PATH")) == 0)
    search_path = DEFAULT_SEARCH_PATH;
#ifdef HAVE_READLINE
  {
    int size;

    /*
    * Get the readline history variables...
    * The defaults are $HOME/.rlab_history
    * and 128 histsize.
    */

    if (use_readline)
    {
      char *home, *tmp;
      if ((rl_histfile = getenv ("RLAB_HISTFILE")) == 0)
      {
        if ((home = getenv ("HOME")) == 0)
        {
          rl_histfile = 0;
        }
        else
        {
          size = strlen (home);
          rl_histfile = (char *) GC_MALLOC (sizeof (char) * (size + 16));
          strcpy (rl_histfile, home);
          strcat (rl_histfile, "/.rlab_history");
        }
      }

      if ((tmp = getenv ("RLAB_HISTSIZE")) == 0)
      {
        rl_histsize = 128;
      }
      else
      {
        rl_histsize = atoi (tmp);
      }
    }

    /* Now, read the history. */
    if (rl_histfile)
      read_history (rl_histfile);
  }
#endif /* HAVE_READLINE */

}

#elif defined(WIN32)
extern char *calculate_home ();

void init_environment ()
{
  int size = 0;
  char *rpath = calculate_home ();

  /* Piece together needed locations. */
  if ((rlab_rc0 = getenv (DEF_RLAB2_RC0)) == 0)
  {
    size = strlen (rpath) + strlen (DEFAULT_RC0) + 2; /* 1 for NULL and 1 for \ */
    rlab_rc0 = (char *) GC_MALLOC (sizeof (char) * size);
    strcpy (rlab_rc0, rpath);
    strcat (rlab_rc0, "\\");
    strcat (rlab_rc0, DEFAULT_RC0);
  }

  if ((help_dir = getenv ("RLAB2_HELP_DIR")) == 0)
  {
    size = strlen (rpath) + strlen (DEFAULT_HELP) + 2;  /* 1 for NULL and 1 for \ */
    help_dir = (char *) GC_MALLOC (sizeof (char) * size);
    strcpy (help_dir, rpath);
    strcat (help_dir, "\\");
    strcat (help_dir, DEFAULT_HELP);
  }

  if ((lib_dir = getenv ("RLAB2_LIB_DIR")) == 0)
  {
    size = strlen (rpath) + strlen (DEFAULT_LIB) + 2; /* 1 for NULL and 1 for \ */
    lib_dir = (char *) GC_MALLOC (sizeof (char) * size);
    strcpy (lib_dir, rpath);
    strcat (lib_dir, "\\");
    strcat (lib_dir, DEFAULT_LIB);
  }

  /*
   * Set the rfile search path.
   */
  if ((search_path = getenv ("RLAB2_PATH")) == 0)
  {
    size = 2 +			/* .;  */
      strlen (rpath) + 1 + strlen (DEFAULT_LIB) + 1 +
      strlen (rpath) + 1 + strlen ("toolbox") + 1 +
      strlen (rpath) + 1 + strlen ("controls-toolbox") + 1 +
      strlen (rpath) + 1 + strlen ("examples") + 1;
    search_path = (char *) GC_MALLOC (sizeof (char) * size);

    strcpy (search_path, ".;");
    strcat (search_path, rpath);
    strcat (search_path, "\\");
    strcat (search_path, DEFAULT_LIB);
    strcat (search_path, ";");
    strcat (search_path, rpath);
    strcat (search_path, "\\toolbox;");

    strcat (search_path, rpath);
    strcat (search_path, "\\controls-toolbox;");

    strcat (search_path, rpath);
    strcat (search_path, "\\examples");

    search_path[size] = '\0';
  }
  else
  {
    search_path = DEFAULT_SEARCH_PATH;
  }

  /*
   * Set the help and default pagers...
   */
  if (use_pager)
  {
    if ((pager = getenv ("RLAB2_PAGER")) == 0)
    {
      if ((pager = getenv ("PAGER")) == 0)
	pager = DEFAULT_PAGER;
    }
    if ((help_pager = getenv ("RLAB2_HELP_PAGER")) == 0)
    {
      help_pager = DEFAULT_HELP_PAGER;
    }
  }
  else
  {
    pager = cpstr ("cat");
    help_pager = cpstr ("cat");
  }
}

#endif /* unix */


void
init_var ()
{
  /* Initialize miscellaneous variables in code.c */
  set_print_machine (print_machine);
  set_line_nos (line_nos);
  set_util_line_nos (line_nos);
  set_use_pager (use_pager);
  set_code_pager (pager);

  /* Initialize miscellaneous variables in rfile.c */
  set_search_path (search_path);
  set_help_dir (help_dir);
  set_lib_dir (lib_dir);
  set_pager (pager);
  set_help_pager (help_pager);
}

/* **************************************************************
 * Run the parser/machine.
 * ************************************************************** */

extern Program *program_Get (void);
extern Inst *get_program_counter (void);
extern void set_program_counter (Inst * prgm);
extern int reset_frame_ptr (void);
extern int reset_stack_ptr (void);
extern Datum *get_stackp (void);
extern void set_stackp (Datum * new_stackp);

#ifdef WIN32
#include <float.h>
#endif

#undef  THIS_SOLVER
#define THIS_SOLVER "run_program"
int run_program (char *input)
{
  int retval;
  Inst *oldpc;
  Program *old_program, *program;

  if (input == 0)
  {
    if (!new_file ("stdin"))
      return (0);
  }
  else
  {
    close_file_ds (input);  /* In case it was open by accident */
    if (!new_file (input))
      return (0);
  }

  program = program_Create (500); /* Create a new program array */
  old_program = program_Get (); /* Save the old program array */
  oldpc = get_program_counter (); /* Save the old counter */

  /* Point the parser at the new program array */
  program_Set (program);

  /* Run the parser/machine */
  while (1)
  {
    if (!setjmp (*jmp_inc_buff ()))
    {
      /* Normal operation */
      signal (SIGFPE, fpecatch);
      signal (SIGINT, intcatch_wait);
#ifdef HAVE_PIPE
#if defined linux
      signal (SIGPIPE, SIG_IGN);
#else
#if defined (WIN32)
      ;
#else
      signal (SIGPIPE, pipecatch);
#endif /* WIN32 */
#endif /* linux */
#endif /* HAVE_PIPE */

//       printf(THIS_FILE ": " THIS_SOLVER ": before run (program)\n");
      run (program);
//       printf(THIS_FILE ": " THIS_SOLVER ": after run (program)\n");

      /* Decrement the jmp buffer counter */
      dec_buff ();
      retval = 1;
      break;
    }
    else
    {
      /* An error (longjmp) has occurred */
#ifdef WIN32
      _fpreset ();
#endif
      reset_frame_ptr ();
      reset_stack_ptr ();
      retval = -1;
      if (input == 0)
        continue;
      else
      {
        /*
         * We error'ed, keep longjmp'ing until we hit
         * input == 0, or we run out of jmps.
         */

        new_file (0);
        flush_line = 0;
        program_Destroy (program);
        program_Set (old_program);
        set_program_counter (oldpc);

        if (get_ijmp () > 0)
          longjmp (*jmp_dec_buff (), 1);
        else
          return (retval);
      }
    }
  }

  /* Reset the old program array, etc, ... */
  program_Destroy (program);
  program_Set (old_program);
  set_program_counter (oldpc);

  close_file_ds (input);
  return (retval);
}

int run_program_debug (char *input)
{
  int retval;
  Inst *oldpc;
  Program *old_program, *program;

  if (input == 0)
  {
    if (!new_file ("stdin"))
      return (0);
  }
  else
  {
    close_file_ds (input);  /* In case it was open by accident */
    if (!new_file (input))
      return (0);
  }

  program = program_Create (500);	/* Create a new program array */
  old_program = program_Get ();	/* Save the old program array */
  oldpc = get_program_counter ();	/* Save the old counter */

  /* Point the parser at the new program array */
  program_Set (program);

  /* Run the parser/machine */
  while (1)
  {
    if (!setjmp (*jmp_inc_buff ()))
    {
      /* Normal operation */
      signal (SIGFPE, fpecatch);
      signal (SIGINT, intcatch_wait);
#ifdef HAVE_PIPE
#if defined (linux)
      signal (SIGPIPE, SIG_IGN);
#else
#if defined (WIN32)
      ;
#else
      signal (SIGPIPE, pipecatch);
#endif /* WIN32 */
#endif /* linux */
#endif /* HAVE_PIPE */

      run_debug (program);

      /* Decrement the jmp buffer counter */
      dec_buff ();
      retval = 1;
      break;
    }
    else
    {
      /* An error (longjmp) has occurred */
      reset_frame_ptr ();
      reset_stack_ptr ();
      retval = -1;
      if (input == 0)
        continue;
      else
      {
       /*
        * We error'ed, keep longjmp'ing until we hit
        * input == 0, or we run out of jmps.
        */
        new_file (0);
        flush_line = 0;
        program_Destroy (program);
        program_Set (old_program);
        set_program_counter (oldpc);

        if (get_ijmp () > 0)
          longjmp (*jmp_dec_buff (), 1);
        else
          return (retval);
      }
    }
  }

  /* Reset the old program array, etc, ... */
  program_Destroy (program);
  program_Set (old_program);
  set_program_counter (oldpc);

  close_file_ds (input);
  return (retval);
}

 /*
  * Run the parser / machine on a character string.
  */

char *eval_string;
extern void set_rlab_input (int type);
extern char *line_contents;
int do_eval = 0;

#undef  THIS_SOLVER
#define THIS_SOLVER  "run_program_eval"
int run_program_eval()
{
  int retval;
  Datum *old_stackp;
  Inst *oldpc;
  Program *old_program, *program;

  /*
   * Set up new program space, and save current state
   * so we can get back.
   */

  program = program_Create (500);	/* Create a new program array */
  old_program = program_Get ();	/* Save the old program array */
  oldpc = get_program_counter ();	/* Save the old counter */
  old_stackp = get_stackp ();	/* Save current stackp */

  /* Point the parser at the new program array */
  program_Set (program);

  /* Run the parser/machine */
  while (1)
  {
    if (!setjmp (*jmp_inc_buff ()))
    {
      /* Normal operation */
      signal (SIGFPE, fpecatch);
      signal (SIGINT, intcatch_wait);
#ifdef HAVE_PIPE
#if defined (linux)
      signal (SIGPIPE, SIG_IGN);
#else
#if defined (WIN32)
      ;
#else
      signal (SIGPIPE, pipecatch);
#endif /* WIN32 */
#endif /* linux */
#endif /* HAVE_PIPE */
      initcode ();
      yyparse ();
      if (print_machine)
      {
        diss_assemble (program->prog, get_progoff ());
      }
      execute (program->prog);

      /* Decrement the jmp buffer counter */
      dec_buff ();
      retval = 1;		/* Successfull return */
      break;
    }
    else
    {
      /*
       * An error (longjmp) has occurred.
       * Do not reset the frame and stack pointers
       * since we will no longjmp back to the prompt.
       * Instead we will return to the caller.
       */

      retval = 0;
      program_Destroy (program);
      program_Set (old_program);
      set_program_counter (oldpc);
      set_stackp (old_stackp);

      return (retval);
    }
  }

  /* Reset the old program array, etc, ... */
  program_Destroy (program);
  program_Set (old_program);
  set_program_counter (oldpc);

  return (retval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "run_program_eval_class"
int run_program_eval_class ()
{
  int retval=1;

  /**
    * Set up new program space, and save current state
    * so we can get back.
    */
  Program *old_program=program_Get ();    /* Save the old program array */
  Program *program=program_Create (500);  /* Create a new program array */;
  Inst *oldpc = get_program_counter ();   /* Save the old counter */;
  Datum *old_stackp = get_stackp ();      /* Save current stackp */

  /* Point the parser at the new program array */
  program_Set (program);
  for (initcode (); yyparse (); initcode ())
  {
    if (print_machine)
    {
      diss_assemble (program->prog, get_progoff ());
    }
    execute (program->prog);

    if (setjmp (*jmp_inc_buff ()))
    {
      retval = 0;
      program_Destroy (program);
      program_Set (old_program);
      set_program_counter (oldpc);
      set_stackp (old_stackp);
      return (retval);
    }

    /* Decrement the jmp buffer counter */
    dec_buff ();
  }

  /* Reset the old program array, etc, ... */
  program_Destroy (program);
  program_Set (old_program);
  set_program_counter (oldpc);

  return (retval);
}


extern int get_class_scope( void );
extern Btree * get_class_publ_symtab( void );

#undef  THIS_SOLVER
#define THIS_SOLVER  "Eval"
Ent * Eval (int nargs, Datum args[])
{
  int retv, class_scope = get_class_scope();
  char *script=0;
  Ent *e1=0, *rent=0;

  is_class_eval = 0;
  eval_string = 0;

  if (nargs != 1)
    rerror ("eval: 1 argument allowed");

  e1 = convert_datum_to_ent (args[0]);
  //e1 = bltin_get_ent (args[0]);

  if (ent_type(e1) != MATRIX_DENSE_STRING)
  {
    if (ent_type(e1)==UNDEF)
    {
      ent_Clean (e1);
      rerror ("eval: Internal Error: Cannot evaluate UNDEF!");
    }

    // return input if not a string
    rent = ent_Copy(e1);
    ent_Clean (e1);
    return (rent);
  }
  script = class_char_pointer (e1);
  int len=isvalidstring(script);
  if (len < 1)
    goto _exit_eval;

  //
  // try to evaluate string
  //
  eval_string = script;
  do_eval++;
  set_rlab_input (1);
  if (class_scope)
  {
    // string is class-defining
    is_class_eval = 1;
    retv = run_program_eval_class ();
  }
  else
  {
    // string is standard rlab expression
    retv = run_program_eval ();
  }

  if (!retv)
  {
    do_eval = 0;
    set_rlab_input (0);
    goto _exit_eval;
  }

  if (--do_eval == 0)
    set_rlab_input (0);

  if (class_scope)
  {
    rent = ent_Create();
    ent_type(rent) = BTREE;
    ent_data(rent) = get_class_publ_symtab();
    ent_Ref (rent) = 1;
  }
  else
  {
    rent = get_eval_ret ();
    if (rent)
    {
      if (ent_type(rent)==UNDEF)
      {
        ent_Destroy(rent);
        rent=0;
      }
    }
  }

_exit_eval:

  ent_Clean (e1);

  if (!rent)
  {
    rent = ent_Create();
  }

  return rent;
}

// Ent *
// Eval_old (int nargs, Datum args[])
// {
//   int retv;
//   Ent *e=0;
//   Ent *eval_ret=0, *ret_ent=0;
// 
//   if (nargs != 1)
//     rerror ("eval: 1 argument allowed");
// 
//   e = convert_datum_to_ent (args[0]);
//   eval_string = class_char_pointer (e);
// 
//   do_eval++;
//   set_rlab_input (1);
//   retv = run_program_eval ();
// 
//   if (retv == 0)
//   {
//     /* Error during eval. */
//     do_eval = 0;
//     set_rlab_input (0);
// 
//     ret_ent = ent_Create ();
//     ent_double (ret_ent) = 0.0;
//     ent_SetType (ret_ent, DOUBLE);
//     return (ret_ent);
//   }
// 
//   if (--do_eval == 0)
//     set_rlab_input (0);
// 
//   /* Now return eval_ret */
//   eval_ret = get_eval_ret ();
// 
//   ent_Clean (e);
//   return (eval_ret);
// }

int run_string (char *evalue)
{
  Ent *eval_ret=0;
  int retv;

  if (isvalidstring(evalue)<1)
    return 0;

  eval_string = evalue;

  do_eval++;
  set_rlab_input (1);
  retv = run_program_eval ();

  if (retv == 0)
  {
    /* Error during eval. */
    do_eval = 0;
    set_rlab_input (0);
    return 1;
  }
  if (--do_eval == 0)
    set_rlab_input (0);

  /* Now return eval_ret */
  eval_ret = get_eval_ret ();

  ent_Clean (eval_ret);
  return (retv);
}

void run (Program *program)
{
  for (initcode (); yyparse (); initcode ())
  {
    if (print_machine)
    {
      diss_assemble (program->prog, get_progoff ());
    }
    execute (program->prog);
  }   
}

void run_no_clean (Program *program)
{
  for (initcode (); yyparse (); initcode ())
  {
    if (print_machine)
    {
      diss_assemble (program->prog, get_progoff ());
    }
    execute (program->prog);
  }
}

void run_debug (Program *program)
{
  for (initcode (); yyparse (); initcode ())
  {
    execute_debug (program->prog);
  }
}

 /* **************************************************************
  * Print a greeting to RLaB users.
  * ************************************************************** */

#if RLAB_VERSION_NUM == 3
static char
    *rlab_greeting_text[] = {
        "_",
        "",
        "Welcome to RLaB3+rlabplus Rel. " version_string " for Linux",
        "RLaB3 and rlabplus(C) 2004-" BUILD_YEAR " Marijan Kostrun, RLaB and RLaB2 (C) 1992-2001 Ian Searle",
        "Please check http://rlabplus.sourceforge.net for the latest news, & c",
        "RLaB3 comes with ABSOLUTELY NO WARRANTY: for details type `help warranty'",
        "New users can type `help INTRO' to get started",
        "This is free software, and you are welcome to redistribute it under",
        "certain conditions; type `help conditions' for details",
        "_",
        "",
        NULL
    };
#else
static char
    *rlab_greeting_text[] = {
        "_",
        "",
        "Welcome to RLaB2+rlabplus Rel. " version_string " for Linux",
        "rlabplus(C) 2004-" BUILD_YEAR " Marijan Kostrun, RLaB and RLaB2 (C) 1992-2001 Ian Searle",
        "Please check http://rlabplus.sourceforge.net for the latest news, & c",
        "RLaB2 comes with ABSOLUTELY NO WARRANTY: for details type `help warranty'",
        "New users can type `help INTRO' to get started",
        "This is free software, and you are welcome to redistribute it under",
        "certain conditions; type `help conditions' for details",
        "_",
        "",
        NULL
    };
#endif


static char *rlab_newline =
#if defined (WIN32)
    "\r\n"
#else
    "\n"
#endif
    ;

#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
void print_greeting ()
{
#if defined (linux)
  struct winsize sz;
  ioctl(0, TIOCGWINSZ, &sz);
  int cols = sz.ws_col;
#else
  int cols = 80;
#endif
  int i=0;
  int j;
  while (rlab_greeting_text[i])
  {
    if (!rlab_greeting_text[i])
    { break; }
    if (strlen(rlab_greeting_text[i]) > 1)
    {
      // position rlab_greeting_text in the center of the console
      if (cols > strlen(rlab_greeting_text[i]))
        for (j=0; j<(int)(0.5 * (cols-strlen(rlab_greeting_text[i]))); j++)
        { printf(" "); }
      printf( "%s", rlab_greeting_text[i]);
    }
    else if (strlen(rlab_greeting_text[i]) == 1)
    {
      if (cols > strlen(rlab_greeting_text[i]))
        for (j=0; j<cols; j+=strlen(rlab_greeting_text[i]))
        { printf("%s", rlab_greeting_text[i]); }
    }
    printf( "%s", rlab_newline);
    i++;
  }
}

 /* **************************************************************
  * Cover for pclose() so we don't get messed up by pclose() that
  * cause a SIGPIPE, when closing an already broken pipe.
  * ************************************************************** */

#ifdef HAVE_PIPE
void rpclose (FILE * fp)
{

#if !defined (WIN32)
  signal (SIGPIPE, SIG_IGN);
#endif
  pclose (fp);
#if defined (linux)
  signal (SIGPIPE, SIG_IGN);
#else
# if defined (WIN32)
  ;
# else
  signal (SIGPIPE, pipecatch);
# endif /* WIN32 */
#endif /* linux */
}

#else

void rpclose (FILE *fp)
{
  fclose (fp);
}

#endif /* HAVE_PIPE */

