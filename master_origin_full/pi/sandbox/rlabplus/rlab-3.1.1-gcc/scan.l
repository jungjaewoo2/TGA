/* scan.l */

/*  This file is a part of RLaB ("Our"-LaB)
    Copyright (C) 1992, 1993, 1994  Ian R. Searle
    Copyright (C) 2015-2016 M. Kostrun

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
 ***********************************************************************/

%{
#include "rlab.h"
#include "code.h"
#include "y.tab.h"
#include "list.h"
#include "util.h"
#include "mem.h"
#include "bltin.h"
#include "print.h"
#include "rfileio.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>

#ifndef FOPEN_MAX
#define FOPEN_MAX 20
#endif

#ifdef titan
extern double strtod ();
#endif

#undef  THIS_FILE
#define THIS_FILE "scan.l"

extern int flush_line;                 /* rlab.y */
static char dummy_line[] = "DUMMY \n"; /* Dummy string */
int lineno = 0;                        /* Keep track of line numbers */
char *curr_file_name = 0;              /* Keep track of the current file */

char *line_contents = 0;               /* For getting at current line */
extern int use_readline;               /* From main.c */
extern void intcatch_wait (int);       /* main.c */

char *rm_escape (char *s);
int   count_eol (char *s);

extern int read_script;
char newclass_def[65536] = {'\0'};
int newclass_count_curly_braces=0;


extern int Rlab_Debug;

/*
 * The two ways RLaB can get its input for the scanner.
 */

typedef void (*RINPUT) ();
RINPUT rinput;   /* Point this at the correct input function. */

void rlab_yy_input     (char *buf, int *result, int *max_size);
void rlab_string_input (char *buf, int *result, int *max_size);
void set_rlab_input (int type);

#undef YY_INPUT
#define YY_INPUT(buf, result, max_size) \
        (*rinput) (buf, &(result), &(max_size))

#undef yywrap
int yywrap (void);
int new_file (char *file_name);

/* Prompt info */
extern int prompt;               /* From rlab.y */
static char prompt1[] = ">> ";    /* Use when prompt == 0, by Marijan Kostrun, just like M*tlab */
static char prompt2[] = "  ";    /* Use when prompt > 0 */

static char *endp;

%}

DIGIT  [0-9]
FN     ([a-zA-Z0-9][a-zA-Z0-9\_\-\.]*)

%x  FLUSH
%x  HELP_COND
%x  RFILE_COND
%x  REQ_COND
%x  SYS_COND
%x  RLABSCRIPT

%%

  /* Always check flush_line */
  if(flush_line)
  {
    BEGIN FLUSH;
  }
  if(read_script)
  {
    BEGIN RLABSCRIPT;
  }

<FLUSH>.*[\r]*\n   {
                /* Eat up the rest of the input line */
                flush_line = 0;
                BEGIN 0; lineno++;
                return '\n';
              }

  /* SCRIPT command */
<RLABSCRIPT>.*  {
                  newclass_def[0] = '\0';
                  newclass_count_curly_braces = 1;
                  int len = 0, k;
                  unsigned char c;
                  while (1)
                  {
                    c = input();
                    if (c == '\n')
                    {
                      ++lineno;
                      if (len)
                      {
                        k = len;
                        while (k--)
                        {
                          if ((newclass_def[k] == ' ')||(newclass_def[k] == '\t'))
                          {
                            len--;
                            continue;
                          }
                          break;
                        }
                      }
                      else
                        continue;
                    }
                    newclass_def[len++] = c;
                    if (c == '{')
                    {
                      newclass_count_curly_braces++;
                    }
                    else if (c == '}')
                    {
                      newclass_count_curly_braces--;
                      if (!newclass_count_curly_braces)
                      {
                        newclass_def[--len] = '\0'; /* ignore closing right curly brace */
                        while(newclass_def[--len]=='\n'){;}
                          newclass_def[++len] = '\0';
                        yylval.name = cpstr(newclass_def);
                        BEGIN(INITIAL);
                        read_script = 0;
                        return R_SCRIPT;
                      }
                    }
                  }
                }


 /* skip blanks and tabs */
[ \t]+  { ; }

 /* skip some control characters */
\013    { /* \f  return */ ; }
\014    { /* \f  form-feed */ ; }
\015    { /* \f  return    */ ; }

  /* what is thy comment? 2016-09-22 */
  /* Recognize comments C-style. No memory leaks. Don't get angry if luser forgets to terminate comment */
"/*"    {
          unsigned char c1=0, c2=input();
          while(1)
          {
            if (c2==0xff)
            {
              printf ("WARNING: File ended with unterminated comment in Line %i\n", lineno);
              break;
            }
            if (c2=='/' && c1=='*')
              break;
            if (c2=='\n')
              ++lineno;
            c1 = c2;
            c2 = input();
          }
          BEGIN(INITIAL);
          if(c2 == 0xff)
          {
            unput(c2);
          }
          return '\n';
        }

  /* Recognize comments: C++, Matlab, shell. */
  /* "//".*\n   { lineno++; } */
  /* "#".*\n    { lineno++; } */
  /* "%".*\n    { lineno++; } */
  /* what is thy comment? 2016-09-22: No memory leaks for files ending with comment line */
"//" |
"#"  |
"%"     {
          unsigned char c;
          while (1)
          {
            c = input();
            if (c=='\n')
            {
              ++lineno;
              break;
            }
            if (c==0xff)
            {
              break;
            }
          }
          /* Maybe you want to place back EOF? It is processed at the end of the scanner */
          BEGIN(INITIAL);
          if(c == 0xff)
          {
            unput(c);
          }
          return '\n';
        }


        /* recognize keywords */
"if"            return IF;
"else"          return ELSE;
"do"            return DO;
"while"         return WHILE;
"until"         return UNTIL;
"for"           return FOR;
"in"            return IN;
"break"         return BREAK;
"continue"      return CONTINUE;
"switch"        return SWITCH;
"case"          return SWITCH_CASE;
"default"       return SWITCH_DEFAULT;
"function"      return FUNCTION;
"return"        return RETURN;
"local"         return LOCAL_DEC;
"global"        return GLOBAL_DEC;
"static"        return FSTATIC_DEC;
"classdef"      return CLASSDEF;
"private"       return CLASS_PRIV;
"public"        return CLASS_PUBL;
"script"        return RSCRIPT_DEC;
"then"          return LOOPTHEN;
"quit"          return QUIT;
"EOF"           return ENDFILE;
"\$self"        return SELF;
"\$\$"          return GST;

  /* Special keywords/commands */
  /* HELP cmd*/
[ \t]*help  { BEGIN HELP_COND; }

<HELP_COND>([ \t]*)/\n  { BEGIN 0; lineno++; return HELP; }

<HELP_COND>{FN}/([ \t]*\n)  { yylval.name = cpstr(yytext); BEGIN 0; lineno++; return HELP_NAME; }

<HELP_COND>[ \t]*  { ; }

  /* rfile filename1 filename2 ... */
[ \t]*"rfile"  { BEGIN RFILE_COND; return RFILE; }

  /* Handle blanks and tabs. */
<RFILE_COND>[ \t]*   { ; }

  /* Multiple file names. */
<RFILE_COND>{FN}   { yylval.name = cpstr (yytext); return FNAME; }

  /* Handle coments on the same line. */
<RFILE_COND>"#".*  { ; }
<RFILE_COND>"//".* { ; }
<RFILE_COND>"%".*  { ; }

  /* End of special rfile handling. */
<RFILE_COND>[\r]*\n   { BEGIN 0; lineno++; return '\n'; }
<RFILE_COND>";"  { BEGIN 0; return '\n'; }

  /* The require command. */
[ \t]*"require"  { BEGIN REQ_COND; return REQUIRE; }

  /* Handle blanks and tabs. */
<REQ_COND>[ \t]* { ; }

  /* Multiple file names. */
<REQ_COND>{FN}   { yylval.name = cpstr (yytext); return FNAME; }

  /* Handle coments on the same line. */
<REQ_COND>"#".*  { ; }
<REQ_COND>"//".* { ; }
<REQ_COND>"%".*  { ; }

  /* End of special rfile handling. */
<REQ_COND>[\r]*\n     { BEGIN 0; lineno++; return '\n'; }

  /* The system (\) command. */
  /* We must do this here because of the LDIV operator. */
^"\\"           { BEGIN SYS_COND; }
<SYS_COND>.*    { yylval.name = cpstr(yytext); return SYS_CMD; }
<SYS_COND>[\r]*\n    { BEGIN 0; lineno++; return '\n'; }

  /* continuation characters */
"..."[ \t]*[\r]*\n   { lineno++; }

({DIGIT}+)/("...") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                   }

  /* read in bit numbers: kmk IV-2015*/
("0"[0-1]+(b|B))  |
("1"[0-1]+(b|B))  {
                    register int _n = strtod (yytext, (char **) &endp);
                    register int _d = 1;
                    register int _r = _n % 10;
                    while (_n)
                    {
                      _d  = 2 * _d;
                      _n  = _n / 10;
                      _r += (_n % 10) * _d;
                    }
                    yylval.d_val = _r;
                    return INTEGER;
                  }

  /* read in integers as ending with l/L or in hex format - prepended 0x/0X - kmk III-2015*/
({DIGIT}+(l|L))   |
"0x"[a-fA-F0-9]*  |
"0X"[a-fA-F0-9]*  {
                    yylval.d_val = strtod (yytext, (char **) &endp);
                    return INTEGER;
                  }

  /* read in real numbers */
{DIGIT}+"."?         |
{DIGIT}*"."{DIGIT}+  |
{DIGIT}+[eE][+-]?{DIGIT}+  {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return NUMBER; }
{DIGIT}*"."{DIGIT}+[eE][+-]?{DIGIT}+  {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return NUMBER; }
{DIGIT}+"."{DIGIT}*[eE][+-]?{DIGIT}+  {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return NUMBER; }

    /* read complex constants */
({DIGIT}+"."?)+(i|j)        |
({DIGIT}*"."{DIGIT}+)+(i|j) |
({DIGIT}+[eE][+-]?{DIGIT}+)+(i|j)  {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return iNUMBER; }
({DIGIT}*"."{DIGIT}+[eE][+-]?{DIGIT}+)+(i|j) {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return iNUMBER; }
({DIGIT}+"."{DIGIT}*[eE][+-]?{DIGIT}+)+(i|j) {
                       errno = 0;
                       yylval.d_val = strtod (yytext, (char **) &endp);
                       if (errno)
                         warning_1 ("decimal over/under flow problem");
                       return iNUMBER; }

 /* read in variables, function names, whatever */
[_a-zA-Z][a-zA-Z0-9\_]* { yylval.name = cpstr(yytext); return(NAME); }

 /* get character strings enclosed in " " \"([^\\\"\n]|\\.)*\" */
\"([^\\\"]|\\.)*\"  {
    yylval.name = cpstr_strip(yytext);
    rm_escape(yylval.name);
    lineno+=count_eol(yylval.name);
    return(R_STRING);
                    }

 /* read in semicolons, through away newline if it follows a ';'
    or a '?'. We have to do this because the parser interprets
    '\n', ';', '?' as end-of-statement. If we sent ;\n to the
    parser it would interpret it a one valid statement and one
    null statement. */

";"       { return SEMIC; }
"?"       { return Q_MARK; }
[\r]*\n   { lineno++; return '\n'; }

  /* read in special operators */
"+="  { return ADDTO; }
"-="  { return SUBFROM; }
"*="  { return EL_MUL_BY; }
"/="  { return EL_DIV_BY; }
"++"  { return INC; }
"--"  { return DEC; }
">="  { return GE; }
">"   { return GT; }
"<="  { return LE; }
"<"   { return LT; }
"!="  { return NE; }
"!"   { return NOT; }
"&&"  { return AND; }
"||"  { return OR; }
"=="  { return EQ; }
"'"   { return TRANSPOSE; }
"<<"  { return LEFT_LIST; }
">>"  { return RIGHT_LIST; }
"</"  { return L_OLIST; }
"/>"  { return R_OLIST; }
".*"  { return EL_MUL_OP; }
"./"  { return EL_DIV_OP; }
"\\"  { return LDIV_OP; }
".\\" { return EL_LDIV_OP; }
".^"  { return EL_POW_OP; }
".'"  { return EL_TRANSPOSE; }

 /* Handle special conditions: such as 3./A */
({DIGIT}+)/(".*") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                  }
({DIGIT}+)/("./") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                  }
({DIGIT}+)/(".\\") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                   }
({DIGIT}+)/(".^") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                  }
({DIGIT}+)/(".'") {
                      sscanf(yytext, "%lf", &(yylval.d_val));
                      return NUMBER;
                  }

.       { return yytext[0]; }

<<EOF>> { return REOF; }

%%

/* **************************************************************
 * Process STRINGS before giving them to yyparse(). This piece of
 * code obtained from MAWK (GNU copylefted) by Mike Brennan.
 * However I have modified it so any bugs are my fault.
 * ************************************************************** */

#define  ET_END     9

static struct
{
  char in;
  char out;
} escape_test[ET_END+1] = {
  { 'n' , '\n' },
  { 't' , '\t' },
  { 'f' , '\f' },
  { 'b' , '\b' },
  { 'r' , '\r' },
  { 'a' , '\07' },
  { 'v' , '\013' },
  { '\\', '\\' },
  { '\"', '\"' },
  { 0 , 0 }
};


/* count '\n' in the string if it is multiline */
int count_eol(char *s)
{
  register char *p = s;
  int count=0;
  while (*p)
  {
    if (*p == '\n'){ ++count; }
    *p++;
  }

  return count;
}

/* process the escape characters in a string, in place . */
char * rm_escape(char *s)
{
  register char *p, *q ;
  int i ;

  q = p = s ;

  while ( *p )
  {
    if ( *p == '\\' )
    {
      escape_test[ET_END].in = * ++p ;
      i = 0 ;
      while ( escape_test[i].in != *p )  i++ ;
      if ( i != ET_END )
      {
        p++ ; *q++ = escape_test[i].out ;
      }
      else if ( *p == 0 )
        *q++ = '\\' ;
      else
	    {
	      *q++ = '\\' ; *q++ = *p++ ;
	    }
    }
    else
    {
      *q++ = *p++ ;
    }
  }
  *q = 0 ;
  return s;
}

/* **************************************************************
 * Set the input function pointer from an external source file.
 * ************************************************************** */
void set_rlab_input (int type)
{
  if (type == 0)
  {
    /* Read from stdin */
    rinput = (RINPUT) rlab_yy_input;
  }
  else
  {
    /* Read from a character string */
    rinput = (RINPUT) rlab_string_input;
  }
}

#ifdef HAVE_READLINE
char *do_gets (char *line, char *prompt);
extern char *readline (char *prompt);
extern void add_history (char *string);
#endif

#ifdef HAVE_INPUT_EDIT
char  *do_gets (char *line, char *prompt);
extern char *getline (char *prompt);
extern void gl_histadd (char *string);
#endif

/* **************************************************************
 * Get input for flex scanner
 * ************************************************************** */
void rlab_yy_input(char *buf, int *result, int *max_size)
{
  char *prmpt;
  int write_diary;

  write_diary = get_write_diary ();

  signal (SIGINT, SIG_IGN);

  if (flush_line)
  {
    strcpy (buf, dummy_line);
    *result = strlen (buf);
    line_contents = buf;
    signal (SIGINT, intcatch_wait);
    return;
  }

  if(prompt)
    prmpt = prompt2;
  else
    prmpt = prompt1;

#if defined(HAVE_READLINE)  || defined(HAVE_INPUT_EDIT)
  if((yyin == stdin) && use_readline)
  {
    /* Interactive, using readline command editing */
    if(do_gets(buf, prmpt) == 0)
      *result = YY_NULL;
    else
      *result = strlen(buf);
    lineno = 0;
  }
  else if((yyin == stdin) && !use_readline)
  {
    /* Interactive, no readline */
    fprintf(stderr, prmpt);

    if((fgets(buf, *max_size, yyin)) == 0)
      *result = YY_NULL;
    else
      *result = strlen(buf);
    lineno = 0;
  }
#else
  if(yyin == stdin)
  {
    /* interactive, put out a prompt */
    fprintf(stderr, prmpt);
    if((fgets(buf, *max_size, yyin)) == 0)
      *result = YY_NULL;
    else
      *result = strlen(buf);
    lineno = 0;
  }
#endif  /* HAVE_READLINE || HAVE_INPUT_EDIT */
  else
  {
    if((fgets(buf, *max_size, yyin)) == NULL)
      *result  = YY_NULL;
    else
      *result = strlen(buf);
  }

  if(write_diary)
  {
    FILE *diary_file_ptr = get_diary_file_ptr ();
    fprintf(diary_file_ptr, "%s", buf);
  }
  else if (Rlab_Debug)
  {
    fprintf(stderr, "%i: %s", lineno, buf);
  }

  /* Set line_contents so parser can get it in case of parser error */
  line_contents = buf;

  signal(SIGINT, intcatch_wait);
}

extern char *eval_string;    /* bltin2.c */

/*
 * Give the scanner input from a character string.
 * eval_string is a global variable.
 * Note that we have to append a '\n' to each string.
 */

void rlab_string_input (char * buf, int * result, int * max_size)
{
  size_t eval_length;

  if (flush_line)
  {
    strcpy (buf, dummy_line);
    *result = strlen (buf);
    line_contents = buf;
    return;
  }

  if (eval_string == 0)
  {
    /* The string has already been handed to the scanner. */
    *result = YY_NULL;
    return;
  }

  eval_length = strlen (eval_string);
  if (eval_length < (*max_size - 1))
  {
    /* Copy the string ONCE */
    strcpy (buf, eval_string);
    strcat (buf, "\n");
    line_contents = buf;
    *result = strlen (buf);
    eval_string = 0;
  }
  else
  {
    /*
     * We are going to have to hand the scanner portions
     * of the eval_string.
     */
  }
}

#ifdef HAVE_READLINE
/* **************************************************************
 * Use GNU Readline library to provide command line editing.
 * Read a string, and return a pointer to it.
 * Returns NULL on EOF.
 * ************************************************************** */
char *do_gets(char * buf, char * prompt)
{
  static char *line_read = 0;
  char *tmp;
  int len;

  /*
   * If the buffer pool has already been allocated,
   * return the memory to the free pool.
   */

  if(line_read != 0)
  {
    free(line_read);
    line_read = 0;
  }

  /* Get a line from the user. */
  line_read = readline(prompt);

  /* If the line has any text in it, save it on the history. */
  if(line_read && *line_read)
    add_history(line_read);

  /* We must put the newline back on the end of line_read */
  if(line_read != 0)
  {
    len = strlen(line_read);
    tmp = GC_MALLOC(sizeof(char)*(len+2));
    if (tmp == 0)
      rerror ("out of memory");
    strcpy(tmp, line_read);
    tmp[len] = '\n';
    tmp[len+1] = '\0';

    /* Both this and the above statments need fixing */
    strcpy(buf, tmp);
    GC_FREE(tmp);
  }
  else
    buf = 0;

  return(buf);
}

#endif  /* HAVE_READLINE */
#ifdef HAVE_INPUT_EDIT

/* **************************************************************
 * Use input-edit library to provide command line editing.
 * Read a string, and return a pointer to it.
 * Returns NULL on EOF.
 * ************************************************************** */
char * do_gets(char *buf, char *prompt)
{
  char *tmp;

  tmp = getline(prompt);
  gl_histadd(tmp);
  if (*tmp != 0)
  {
    strcpy(buf, tmp);
  }
  else
  {
    buf = 0;
  }
  return buf;
}

#endif  /* HAVE_INPUT_EDIT */

static int lnstack[ FOPEN_MAX ];
static char *fnstack[ FOPEN_MAX ];
static FILE *fstack[ FOPEN_MAX ];
static int bptr = 0;

/* **************************************************************
 * Replacement yywrap().
 * Yywrap is called by scanner, when it thinks it is at end-of-input.
 * If yywrap returns 0, then it is assumed that the functions has
 * gone ahead and setup yyin to point to another input file. If
 * yywrap returns 1, then the scanner terminates, returning 0 to
 * it's caller.
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "yywrap"
int yywrap(void)
{
  if(!new_file(0))
  {
    return(1);
  }

  return(0);
}

/* **************************************************************
 * Tell flex to read from another source. If called with a 0 ptr,
 * decrement the buffer. If called with a string, then try to open
 * the file and increment the buffer.
 * Return 1: if buffer is successfuly inc'ed, or dec'ed
 * Return 0: if unsuccessful
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "new_file"
int new_file(char *file_name)
{
  static int once = 0;
  FILE *tmp=0;

  if(file_name != 0)
  {
    /* Check for stack over-run */
    if(bptr > FOPEN_MAX - 1)
    {
      fprintf(stderr, "ERROR: too many input files\n");
      return(0);
    }
    if(!strcmp("stdin", file_name))
    {
      yyin = stdin;
      fstack[bptr] = yyin;
    }
    else
    {
      /* Open new file, and point flex(1) at it */
      if((tmp = get_file_ds(file_name, "r", 0)) == 0)
      {
        fprintf(stderr, "ERROR: %s, cannot open for read\n", file_name);
        return(0);
      }
      else
      {
        yyin = tmp;
        fstack[bptr] = yyin;
      }
    }
    if(once)                 /* Why did I do this ? */
      yyrestart(yyin);
    else
      once = 1;

    /* Point cur_file_name to current input file */
    if (fnstack[bptr])
      GC_FREE(fnstack[bptr]);
    fnstack[bptr] = cpstr(file_name);
    if (curr_file_name)
      GC_FREE(curr_file_name);
    curr_file_name = cpstr(file_name);

    /* Save curent line number on stack */
    lnstack[bptr] = lineno;

    /* Now set lineno for new file */
    lineno = 1;

    bptr++;
    return(1);
  }
  else
  {
    /* Decrement to current file */
    bptr--;

    /* Close the old file */
    if (fstack[ bptr ] != stdin)
    {
      close_file_ds( fnstack[ bptr ] );
      GC_FREE(fnstack[bptr]);
      fnstack[bptr]=0;
    }

    /* Dec bptr again to point to the next file */
    bptr--;
    if( bptr >= 0 )
    {
      if (curr_file_name)
        GC_FREE(curr_file_name);
      curr_file_name = cpstr(fnstack[bptr]);
      yyin = fstack[bptr];
      yyrestart(yyin);
      /* lineno = lnstack[bptr+1] + 1; */
      lineno = lnstack[bptr+1];

      bptr++;    /* Inc the ptr for next call to new_file() */
    }
    else     /* We're out of files */
    {
      if( fstack[0] == stdin )
      {
        if (curr_file_name)
          GC_FREE(curr_file_name);
        curr_file_name = cpstr("stdin");
        yyin = stdin;
        yyrestart(yyin);
        /* lineno = lnstack[0] + 1; */
        lineno = lnstack[0];
        bptr = 1;
      }
      else
      {
        bptr = 0;
        /* lineno = lnstack[0] + 1; */
        lineno = lnstack[0];
      }
    }
    return(0);
  }
}

/*
 * Set the zero index FILEstack ptr to stdin.
 */

void set_input_stdin()
{
  yyin = stdin;
  fstack[0] = yyin;
}

/*
 * Call yyrestart in case of longjmp
 */
void rlab_yyrestart()
{
  yyrestart(yyin);
}

/* **************************************************************
 * Free up flex's input buffers.
 * ************************************************************** */
void scanner_cleanup()
{
  int i;

  for(i = 0; i < FOPEN_MAX; i++)
  {
    if (fnstack[i])
      GC_FREE(fnstack[i]);
    fnstack[i]=0;
  }
  if (curr_file_name)
    GC_FREE(curr_file_name);
  curr_file_name = 0;
}
