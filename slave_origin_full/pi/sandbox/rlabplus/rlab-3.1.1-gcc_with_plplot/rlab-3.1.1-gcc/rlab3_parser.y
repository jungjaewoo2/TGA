/* rlab_parser.y */

/*  This file is a part of RLaB ("Our"-LaB) and rlab3
    Copyright (C) 1995  Ian R. Searle
    Copyright (C) 2001-2016 Marijan Kostrun

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
#include "ent.h"
#include "mem.h"
#include "rfile.h"
#include "list.h"
#include "code.h"
#include "symbol.h"
#include "util.h"
#include "function.h"

#include <setjmp.h>
#include "rlab_solver_parameters_names.h"
#undef  THIS_FILE
#define THIS_FILE "rlab_parser.y"
#include "rlab_macros.h"

/* Trick yyparse() so we can pass yychar into yyerror() */
#define yyerror(s)  ryyerror(s, yychar)
void ryyerror (char *s, int yychar);

/* Lookup variables in symbol table */
static Var * name_lookup      (List *l1, List *l2, List *l3, List *l4, char *name, int scope);
static Var * name_lookup_func (List *l1, List *l2, List *l3, List *l4, char *name, int scope);

void static_tree_DestroyNodeByKey (char *key);

extern Datum * get_stackp (void);
extern Datum pop(void);
extern void  push(Datum);
extern char *line_contents;   /* scan.l */
extern int   lineno;            /* scan.l */
extern char *curr_file_name;  /* scan.l */
extern int   line_nos;          /* main.c */

char  priv_class_name[12]={'\0'};
char  *classdef_filename=0;

extern int do_eval;           /* bltin2.c, main.c */

void fstatic_var_push (char *file_name, char *name);
ListNode * static_lookup (char *file_name, char *name);
static void resolve_continue_tags (int bop, int eop, int tloop);
static void resolve_return_tags (int bop, int eop);

static void resolve_break_tags (int bop, int eop, int tloop);
static void tag_brk_cont_case (List **list_ptr, int bop);
static void resolve_case_break_tags (List **blist, List **clist, List **elist,
                                     int begin_op, int end_op, int dflt_end);

static int for_then=0;

List *static_tree;
int flush_line = 0;           /* Tells yylex() when to flush rest of line */
int prompt = 0;               /* prompt=0: put out prompt-1, p != 0: prompt-2 */

// script support
int read_script = 0;

//
// switch, for-then, while-then, do-while/until all use 
//  break and continue tags, some also use case-end tags
// 
static List *blist = 0;       // The break list
static List *clist = 0;       // The continue list.
static List *rlist = 0;       // The return list.
static Datum sw_switch_datum;
int sw_progoff_default;
extern Datum op_pull_datum;
#define switch_clist_offset  (100000L)
static List *sw_blist = 0;    /* The break list for switch commands */
static List *sw_clist = 0;    /* The case list for switch commands */
static List *sw_elist = 0;    /* The case end list for switch commands */
static int looping = 0;       /* for keeping track of BREAK, CONTINUE */
static int switching = 0;     /* for keeping track of BREAK inside switch() */

// functions after rlab3
List *arg_sym_tab=0;            /* For tracking argument lists */
List *lsymtab=0;                /* For tracking local() statement lists */
List *gsymtab=0;                /* For tracking global() statement lists */
List *ssymtab=0;                /* For tracking static() statement lists inside functions */

// classdef: definition of class member creation function
int   class_scope=0;        // 0-normal operation, 1-inside classdef(){..}
List *class_args_symtab=0;  // class member instantiation variables
List *class_glob_symtab=0;  // For tracking global() variables inside classdef
Btree *class_publ_symtab=0; // For tracking class-member public() variables
Btree *class_stat_symtab=0; // For tracking class-member static() variables
void   btree_class_publ_var_push (Btree *class_publ_symtab, char *name);
void   btree_class_stat_var_push (Btree *class_stat_symtab, char *name);

int scope = GLOBAL;           /* Determines current scope */
static int psave = 0;         /* Save program pointer (offset) */
static int lsave = 0;         /* Save line number */
int loff = 100;               /* The line # offset, to avoid op-code collisions */
int i;

static int dlocal = 0;
static int dglobal = 0;

%}    /* Declarations */

%start program

%union {
  List       *list;   /* ptr to a linked list */
  ListNode   *ent;    /* pointer to an RLaB entity */
  double      d_val;  /* double numbers */
  int         n_int;  /* number of items in a list */
  char       *name;   /* char ptr to name string passed from scanner */
  struct _nn {
    int off;
    int count;
  } nn;               /* For machine offset, and count */
  Var *vinfo;         /* Global, arg, and local variable info */
}

%token <d_val> INTEGER NUMBER iNUMBER REOF ENDFILE
%token <name>  NAME R_STRING HELP_NAME FNAME UNDEFINED SYS_CMD R_SCRIPT CLASS_PUBL CLASS_PRIV
%token <n_int> DO UNTIL WHILE IF ELSE QUIT FOR IN BREAK CONTINUE LOOPTHEN
%token <n_int> SWITCH SWITCH_CASE SWITCH_DEFAULT
%token <n_int> FUNCTION RETURN SEMIC Q_MARK LOCAL_DEC GLOBAL_DEC SELF
%token <n_int> HELP RFILE REQUIRE GST FSTATIC_DEC RSCRIPT_DEC CLASSDEF JNK

%type  <n_int>  program stmt c_stmt stmts line reof cstmts cline NL
%type  <n_int>  texpr_or_empty rfile_stmt req_stmt
%type  <n_int>  expr vec_expr mat_expr mat_list texpr vexpr value
%type  <n_int>  cond do while until if for for_spec end function class switch_label
%type  <n_int>  if_stmt_open if_stmt while_stmt for_stmt return_stmt do_stmt switch_stmt
%type  <n_int>  break continue else then
%type  <n_int>  switch_case_list switch_case_list1 switch_case switch_default_case switch_case_open
%type  <n_int>  vid for_vec_expr assign list_expr list_el olist 
%type  <n_int>  separator newline newlines list_member nl
%type  <list>   opt_arg_name_list class_inst_var_list
%type  <n_int>  class_private_stmt opt_class_private_name_list
%type  <n_int>  class_public_stmt opt_class_public_name_list
%type  <n_int>  fstatic_stmt opt_fstatic_name_list 
%type  <n_int>  gl_stmt opt_local_name_list
%type  <vinfo>  var
%type  <ent>    function_ent class_ent
%type  <nn>     opt_arg_list list_list texpr_list olist_list
%type  <name>   rnames rscript_stmts

/* Operator Precedence and Associativity */
%nonassoc JNK
%nonassoc DO
%nonassoc INTEGER NUMBER iNUMBER NAME R_STRING WHILE IF QUIT
%nonassoc SWITCH SWITCH_CASE SWITCH_DEFAULT
%nonassoc FOR BREAK CONTINUE FUNCTION RETURN SELF GST FSTATIC_DEC 
%nonassoc ELSE 
%right '='
%left  ','              /* concatenation operator */
%left  ':'              /* vector creation operator */
%left  OR               /* || */
%left  AND              /* && */
%left  EQ NE            /* == != */
%left  GT GE LT LE      /* > >= < <= */
%left  ADDTO SUBFROM EL_MUL_BY EL_DIV_BY
%left  '+' '-'
%left  '*' '/' EL_MUL_OP EL_DIV_OP EL_LDIV_OP LDIV_OP
%nonassoc UNARY_MINUS UNARY_PLUS
%right '^' EL_POW_OP
%nonassoc TRANSPOSE EL_TRANSPOSE
%nonassoc NOT INC DEC
%nonassoc LEFT_LIST RIGHT_LIST L_OLIST R_OLIST
%left '[' ']' '.'

%%      /* Grammar */

program
: /* empty */       { $$ = 0; }
| program REOF      { code(OP_ENDOFFILE); return(0); }
| program line      { code(STOP); return(1); }
| program HELP      { code(OP_HELP); code(STOP); return (1); }
| program HELP_NAME { $$ = 0; code(OP_HELP_NAME); codep ($2); code(STOP); return (1); }
| program SYS_CMD   { code(OP_SYS_CMD); codep($2); code(STOP); return (1); }
| program error
;

reof
: /* empty */ { $$ = 0; }
| REOF        { prompt--; code(OP_ENDOFFILE); return(0); }
;

 /*
  * Line
  */
line
: NL          { $$ = get_progoff (); }
| ENDFILE     { code(OP_JMP_ENDOFFILE); }
| stmts NL    { $$ = $1; }
| stmts REOF  {
                $$ = $1;
                printf("WARNING: File %s ended with statement and not an empty line or comment: Memory leak likely!\n",
                      curr_file_name);
                code(OP_ENDOFFILE);
              }
;
 /*
  * Statements
  */

stmts
: stmt
| stmts stmt
;

stmt
: ENDFILE { code(OP_JMP_ENDOFFILE); }
| texpr separator
{
  if (do_eval && !class_scope)
  {
    /* Only for eval() */
    code (OP_SAVE_EVAL);
  }
  else
  {
    /* Normal operation */
    switch ($2)
    {
      case 0:
        code (OP_POP_CLEAN);
        break;

      case 1:
        code (OP_PRINT);
        break;

      default:
        code (OP_PRINT_ASSIGN);
    }
    if (line_nos)
    {
      code (OP_LINE_NO);
      code (lineno+loff);
    }
  }
}
| if_stmt_open
| if_stmt
| do_stmt
| while_stmt
| for_stmt
| switch_stmt
| break
| continue
| return_stmt
| c_stmt
| gl_stmt
| fstatic_stmt
| class_public_stmt
| class_private_stmt
| rfile_stmt
| req_stmt
| QUIT                 { code(OP_QUIT); }
;


separator
: /* empty */ { $$ = 2; /*asign to 'ans' and print*/}
| SEMIC       { $$ = 0; }
| Q_MARK      { $$ = 1; /*print*/}
;

rfile_stmt
: RFILE         { code(OP_RFILE); }
| RFILE rnames  { code(OP_RFILE_NAME); codep($2); }
;

req_stmt
: REQUIRE rnames   { code(OP_REQ_NAME); codep($2); }
;

rnames
: FNAME           { $$ = $1; }
| rnames FNAME    { $$ = strappend ($1, $2); GC_FREE ($1); GC_FREE ($2); }
;

left_paren: '('
;

right_paren: ')'
;

left_brace: '{'
;

right_brace: '}'
;

else
: ELSE
  {
    $$ = code (OP_JMP);
    code (STOP);
  }
;

then
: LOOPTHEN ':'
  {
    if (!looping)
    {
      yyerror("then not allowed outside for-loop");
      longjmp( *jmp_dec_buff (), 1 );
    }
    if (for_then!=-1 && for_then!=-2)
    {
      yyerror("only one 'then:' allowed inside for-then or while-then loop");
      longjmp( *jmp_dec_buff (), 1 );
    }
    if (for_then == -1)
    {
      // for-then loop
      code (OP_FOR_THEN);
      for_then = get_progoff ();
      code (STOP); /* placeholder */
      code (OP_FOR_LOOP);
    }
    else if (for_then == -2)
    {
      // while-then loop
      code(OP_JMP);
      for_then = get_progoff ();
      code(STOP);                           /* Placeholder */
    }
  }

cline
: NL    { $$ = get_progoff (); }
| then  
| stmts  %prec JNK
; 

cstmts
: cline
| cstmts cline
;


c_stmt: left_brace cstmts right_brace  { $$ = $2; }
;

 /* kmk gratefuly acknowledges the following sources:
  *   "Resolving the General Dangling Else/If-Else Ambiguity", http://www.parsifalsoft.com/ifelse.html
  *   "if elsif else statement parsing", http://stackoverflow.com/questions/9164907/if-elsif-else-statement-parsing
  */
if_stmt
: if cond newline c_stmt newline end else newline if_stmt_open end
{
  $$ = $1;
  code_sp ($2 + 1, $9 - ($2 + 1));  // if 'false' jump to next if statement
  code_sp ($7 + 1, $10 - ($7 + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
| if cond newline c_stmt newline end else newline if_stmt end
{
  $$ = $1;
  code_sp ($2 + 1, $9 - ($2 + 1));  // if 'false' jump to next if statement
  code_sp ($7 + 1, $10 - ($7 + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
| if cond newline c_stmt newline end else newline c_stmt end
{
  $$ = $1;
  code_sp ($2 + 1, $9 - ($2 + 1));  // if 'false' jump to after 'else' statement
  code_sp ($7 + 1, $10 - ($7 + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
;

if_stmt_open
: if cond newline c_stmt newline reof end %prec IF
{
  /* if 'cond' $2 is not true, jump to 'end' $7 */
  $$ = $1;
  code_sp ($2 + 1, $7 - ($2 + 1));
  prompt--;
}
;


/*
 * switch command - start, kmk X-2016
 */
switch_stmt
: switch_label left_paren expr right_paren
  {
    /* remove this from stack */
    code(OP_PULL_DATUM);
    sw_switch_datum = op_pull_datum; 
  }
  newline '{' newline switch_case_list newline '}' end
  {
    $$ = $9;
    switching--;
    prompt--;
    int sw_end_jump = (sw_progoff_default > 0 ? sw_progoff_default : $12);
    resolve_case_break_tags (&sw_blist, &sw_clist, &sw_elist, $$, $12, sw_end_jump);
  }
;

switch_label
: SWITCH
  {
    sw_progoff_default=-1;
    prompt++;
    switching++;
    $$ = get_progoff ();
  }
;

switch_case_list
: /* empty */           { $$=0; }
| switch_default_case   { $$=$1; }
| switch_case_list1     { $$=$1; }
| switch_case_list1 switch_default_case { $$=$1; }
;

switch_case_list1
: switch_case                   { $$=$1; }
| switch_case_list1 switch_case { $$=$1; }
;

switch_case_open
: SWITCH_CASE
  {
    tag_brk_cont_case (&sw_clist, get_progoff()); /* location of prev OP_IFSJMP */
  }
  expr ':' 
  { /* value is already on the stack */
    code(OP_PUSH_DATUM);
    code(OP_EQ);                            /* are they equal */
    tag_brk_cont_case (&sw_elist, get_progoff());       /* location of next case statement OP_IFSJMP */
    code(OP_IFSJMP);                        /* where to jump if not true */
    code(STOP);                             /* placeholder for next case command */
  }
;

switch_case
: switch_case_open newline cstmts newline break /* full case statement */
| switch_case_open newline cstmts newline       /* case statement without break at the end */
  { /* do a fake break switch that points to next case statement */
    $$ = get_progoff();
    code(OP_JMP);
    code(STOP);                           /* Placeholder */
    int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
    tag_brk_cont_case (&sw_blist, -n*switch_clist_offset-$$);   /* for fake brake switch */
  }
| switch_case_open newline                /* empty case statement */
  { /* do a fake break switch that points to next case statement */
    $$ = get_progoff();
    code(OP_JMP);
    code(STOP);                           /* Placeholder */
    int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
    tag_brk_cont_case (&sw_blist, -n*switch_clist_offset-$$);   /* for fake brake switch */
  }
;

switch_default_case
: SWITCH_DEFAULT ':' newline
  {
    sw_progoff_default = get_progoff();
  }
  cstmts
;

/*
 * switch command - end
 */

do
: DO
{
  $$ = get_progoff ();
  looping ++;
}
;

do_stmt
: do newline c_stmt newline while cond separator end
  {
    $$ = $1;
    code_sp ($6 + 1, $8 - $6 + 1 );
    code (OP_JMP);
    code ($1 - $8 - 1);
    resolve_break_tags ($$, $8 + 2, 0);
    resolve_continue_tags ($$, $8 + 2, 0);
    prompt--;
    looping-=2; /* key word 'while' does looping++ too */
  }
| do newline c_stmt newline until cond separator end
  {
    $$ = $1;
    code_sp ($6 + 1, $1 - $8 - 1);
    code (OP_JMP);
    code ($2 + 1);
    resolve_break_tags ($$, $8 + 2, 0);
    resolve_continue_tags ($$, $8 + 2, 0);
    prompt--;
    looping--;
  }
;

  /* while and while-then construct, kmk X-2016 */
while_stmt
: while cond newlines c_stmt end
  {
    $$ = $1;
    if (for_then > 0)
    {
      code_sp ($2 + 1, for_then - $2);
      code_sp (for_then, $1 - for_then);
      resolve_continue_tags ($$, for_then - 1, 0);
      resolve_break_tags ($$, $5, 0);
    }
    else
    {
      code_sp ($2 + 1, $5 - $2 + 1);
      code (OP_JMP);
      code ($1 - $5 - 1);
      resolve_continue_tags ($$, $5 + 2, 0);
      resolve_break_tags ($$, $5 + 2, 0);
    }
    looping--;
    prompt--;
    for_then=0;
  }
;


  /* for and for-then construct, kmk X-2016 */
for_stmt
: for for_spec newlines c_stmt end
  {
    $$ = $1;
    if (for_then > 0)
    {
      code_sp(for_then, $5 - for_then);
    }
    code_sp ($2 + 1, $5 - ($2+1));
    code (OP_FOR_LOOP);
    code (OP_FOR_LOOP_DONE);
    resolve_break_tags ($$, $5, 1);
    if (for_then > 0)
    {
      resolve_continue_tags ($$, for_then-1, 1);
    }
    else
    {
      resolve_continue_tags ($$, $5, 1);
    }
    resolve_return_tags ($$, $5);
    looping--;
    prompt--;
    for_then=0;
  }
;

for
: FOR
  {
    $$ = get_progoff ();
    looping++;
    prompt++;
    for_then=-1;
  }
;

for_spec
: left_paren vid IN for_vec_expr right_paren
  {
    $$ = code (OP_FOR_LOOP_I);
    code (STOP);
  }
;

for_vec_expr: texpr
             ;

return_stmt
: RETURN texpr separator
  {
    int tmp;
    $$ = $2;
    if(scope == GLOBAL)
    {
      yyerror("return not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    code (OP_FUNCTION_RETURN);
    tmp = code (STOP);
    tag_brk_cont_case(&rlist, tmp);
  }
;

cond
: left_paren texpr right_paren
  {
    $$ = code(OP_IFSJMP);
    code (STOP);               /* place holder for jump val */
  }
;

while
: WHILE
  {
    for_then=-2;
    $$ = get_progoff ();
    looping ++;
    prompt++;
  }
;

until
: UNTIL
  {
    $$ = get_progoff ();
    looping ++;
    prompt++;
  }
;

if
: IF
  {
    $$ = get_progoff ();
    prompt++;
  }
;

break
: BREAK separator
  {
    if (switching)
    {
      int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
      $$ = code (OP_JMP);
      code(STOP);               /* Placeholder */
      tag_brk_cont_case (&sw_blist, n*switch_clist_offset + $$);
    }
    else if(looping)
    {
      if (for_then > 0)
        if (get_progoff() > for_then )
        {
          yyerror("break-stmt cannot be after then-stmt in the same loop");
          longjmp( *jmp_dec_buff (), 1 );
        }
      $$ = code (OP_JMP);
      code(STOP);               /* Placeholder */
      tag_brk_cont_case (&blist, $$);
    }
    else
    {
      yyerror("break-stmt not allowed outside loop or switch environments");
      longjmp( *jmp_dec_buff (), 1 );
    }
  }
;

continue: CONTINUE separator
  {
    if(!looping)
    {
      yyerror("continue-stmt not allowed outside loop");
      longjmp( *jmp_dec_buff (), 1 );
    }
    if (for_then > 0)
      if (get_progoff() > for_then )
      {
        yyerror("continue-stmt cannot be after then-stmt in the same loop");
        longjmp( *jmp_dec_buff (), 1 );
      }
    $$ = code(OP_JMP);
    code(STOP);            /* Placeholder */
    tag_brk_cont_case (&clist, $$);
  }
;

end
: /* empty */    { $$ = get_progoff (); }
;

 /*
  * Expressions
  */

 /* Top level expression */
texpr
: expr %prec '='
| vexpr
| assign
;

expr
: '(' texpr ')'     { $$ = $2; }
;

 /* Vector expression */
vexpr
: expr ':' expr  { code(OP_VECTOR_CREATE); code(2); }
| expr ':' expr ':' expr { code(OP_VECTOR_CREATE); code(3); }
;

 /* Assignments */
assign
: vid '=' texpr                                   { code(OP_ASSIGN); }
 /* matrix assign */
| expr '[' vec_expr SEMIC vec_expr ']' '=' texpr  { code(OP_MATRIX_ASSIGN); code(1); }
| expr '[' vec_expr SEMIC ']' '=' texpr           { code(OP_MATRIX_ASSIGN); code(2); }
| expr '[' SEMIC vec_expr ']' '=' texpr           { code(OP_MATRIX_ASSIGN); code(3); }
| expr '[' vec_expr ']' '=' texpr                 { code(OP_MATRIX_VEC_ASSIGN); }
 /* list assign */
| expr '.' '[' expr ']' '=' texpr                 { code(OP_LIST_ASSIGN); code(1); }
| expr '.' NAME '=' texpr                         { code(OP_LIST_ASSIGN); code(2); codep($3); }
| olist '=' texpr                                 { code(OP_OLIST_ASSIGN); }
;

 /* 
  * Matrix expression 
  */
mat_expr
: '[' mat_list ']'  { $$ = $2; }
;

expr
: mat_expr
;

mat_list
: /* empty matrix */       { $$ = code(OP_EMPTY_MATRIX_CREATE); }
| vec_expr                 { code(OP_MATRIX_CREATE); }
| mat_list SEMIC vec_expr  { code(OP_MATRIX_APPEND); }
| mat_list SEMIC nl vec_expr  { code(OP_MATRIX_APPEND); }
;

 /* Vector concatenation */
 vec_expr: expr
         | vexpr
         | vec_expr ',' vec_expr  { code(OP_VEC_APPEND); }
         ;

 /* sub-matrix */
expr
: expr '[' vec_expr SEMIC vec_expr ']'  {  code(OP_MATRIX_SUB); code(1); }
| expr '[' vec_expr SEMIC ']'           {  code(OP_MATRIX_SUB); code(2); }
| expr '[' SEMIC vec_expr ']'           {  code(OP_MATRIX_SUB); code(3); }
| expr '[' SEMIC ']' { /* no-op */ }
| expr '[' ']'       { /* no-op */ }
| expr '[' vec_expr ']'                 {  code(OP_MATRIX_VEC_SUB); }
| expr '[' ':' ']'                      {  code(OP_MATRIX_COL); }
;

expr: function
;

expr: class
;

 /*
  * List Expression
  */
expr
: list_expr
| list_member
;

 /* Experimental */
list_el
: NAME '=' texpr  { $$ = $3; code(OP_LIST_EL_CREATE); codep($1); }
| texpr
;

list_el_term
: SEMIC { /* default action causes type conflict */ }
| NL {}
;

left_list
: LEFT_LIST
| LEFT_LIST SEMIC
| LEFT_LIST NL
;

right_list
: RIGHT_LIST
| list_el_term RIGHT_LIST
;

list_expr
: left_list list_list right_list
  {
    $$ = $2.off;
    code(OP_LIST_CREATE);
    code($2.count);
  }
;

list_member
: expr '.' '[' expr ']' { code(OP_LIST_MEMB); code(1); }
| expr '.' NAME         { code(OP_LIST_MEMB); code(2); codep($3); }
;

list_list
: /* empty */        { $$.count = 0; $$.off = get_progoff (); }
| list_el            { $$.count = 1; $$.off = $1; }
| list_list list_el_term list_el { $$.count = $1.count + 1; $$.off = $1.off; }
;

 /*
  * Open-List 
  */

olist
: L_OLIST olist_list R_OLIST
  {
    $$ = $2.off;
    code(OP_OLIST);
    code($2.count);
  }
;

 olist_list
: vid
  {
    $$.count = 1;
    $$.off = $1;
  }
| olist_list SEMIC vid
  {
    $$.count = $1.count + 1;
    $$.off = $1.off;
  }
;

 /*
  * Casting/Coercion expresion.
  */

 /*  expr: '(' NAME expr ')' { } */

 /* 
  * General Expressions
  */
expr
: value { $$=$1; }
| vid
| vid INC   { code(OP_INC); }
| vid DEC   { code(OP_DEC); }
| NAME '(' opt_arg_list ')'
  {
    Var *ret = name_lookup_func (arg_sym_tab, lsymtab, gsymtab, ssymtab, $1, scope);
    if (ret->type == GLOBAL)
    {
      code(OP_PUSH_VAR);
      codep(ret->var);
    }
    else if(ret->type == STATIC_VAR)
    {
      code(OP_PUSH_STATIC_VAR);
      code(ret->offset);
    }

    else if(ret->type == LOCAL_VAR)
    {
      code(OP_PUSH_LOCAL_VAR);
      code(ret->offset);
    }
    else if(ret->type == ARG_VAR)
    {
      code(OP_PUSH_ARG);
      code(ret->offset);
    }
    GC_FREE(ret);
    code(OP_FUNCTION_CALL);
    code($3.count);   /* number of args on stack */
    $$ = $3.off;
  }
| expr '(' opt_arg_list ')'
  {
    /* MDEs are here to stay! MK IX-2016 */
    code(OP_FUNCTION_CALL_2);
    code($3.count);   /* number of args on stack */
    $$ = $3.off;
  }
| self '(' opt_arg_list ')'
  {
    code(OP_FUNCTION_CALL_SELF);
    code($3.count);   /* number of args on stack */
    $$ = $3.off;
  }
| list_member '(' opt_arg_list ')'
  {
    code(OP_FUNCTION_CALL_1);
    code($3.count);   /* number of args on stack */
  }
| expr '+' expr     { code(OP_ADD); }
| expr '-' expr     { code(OP_SUB); }
| expr '*' expr     { code(OP_MUL); }
| expr EL_MUL_OP expr   { code(OP_EL_MUL); }
| expr '/' expr         { code(OP_DIV); }
| expr LDIV_OP expr     { code(OP_LDIV); }
| expr EL_DIV_OP expr   { code(OP_EL_DIV); }
| expr EL_LDIV_OP expr  { code(OP_EL_LDIV); }
| expr '^' expr         { code(OP_POWER); }
| expr EL_POW_OP expr   { code(OP_EL_POWER); }
| expr TRANSPOSE    { code(OP_TRANSPOSE); }
| expr EL_TRANSPOSE { code(OP_EL_TRANSPOSE); }
| expr GT expr      { code(OP_GT);  }
| expr GE expr      { code(OP_GE);  }
| expr LT expr      { code(OP_LT);  }
| expr LE expr      { code(OP_LE);  }
| expr EQ expr      { code(OP_EQ);  }
| expr NE expr      { code(OP_NE);  }
| expr AND expr     { code(OP_AND); }
| expr OR expr      { code(OP_OR);  }
| expr ADDTO expr     { code(OP_ADDTO);     }
| expr SUBFROM expr   { code(OP_SUBFROM);   }
| expr EL_MUL_BY expr { code(OP_EL_MUL_BY); }
| expr EL_DIV_BY expr { code(OP_EL_DIV_BY); }
| NOT expr          { $$ = $2; code(OP_NOT); }
| '-' expr %prec UNARY_MINUS  { $$ = $2;   code(OP_NEGATE); }
| '+' expr %prec UNARY_PLUS   { $$ = $2; }
;

value
: NUMBER    { $$ = code(OP_PUSH_CONSTANT); coded($1); }
| INTEGER   { $$ = code(OP_PUSH_INTEGER); coded($1); }
| iNUMBER   { $$ = code(OP_PUSH_iCONSTANT); coded($1); }
| R_STRING  { $$ = code(OP_PUSH_STRING); codep($1); }
| rscript_stmts  { $$ = code(OP_PUSH_STRING); codep($1); }
;

 /*
  * Variables
  */
vid
: var
  {
    if($1->type == GLOBAL)
    {
      $$ = code(OP_PUSH_VAR);
      codep($1->var);
    }
    else if($1->type == LOCAL_VAR)
    {
      $$ = code(OP_PUSH_LOCAL_VAR);
      code($1->offset);
    }
    else if($1->type == STATIC_VAR)
    {
      $$ = code(OP_PUSH_STATIC_VAR);
      code($1->offset);
    }
    else if($1->type == ARG_VAR)
    {
      $$ = code(OP_PUSH_ARG);
      code($1->offset);
    }
    GC_FREE($1);
  }
;

var
: NAME
  {
    $$ = name_lookup (arg_sym_tab, lsymtab, gsymtab, ssymtab, $1, scope);
  }
| GST   { $$ = gst (); }
;

nl
: NL    { $$ = 0;}
| nl NL          { $$ = 0; }
;

newline
: /* empty */    { $$ = 0; }
| NL             { $$ = 0; }
;

newlines
: /* empty */    { $$ = 0; }
| NL             { $$ = 0; }
| newlines NL    { $$ = 0; }
;

NL
: '\n'  { /* empty */ }
;


 /*
  * Functions
  */
function
: function_ent
  {
    $$ = code(OP_PUSH_VAR); codep($1);
  }
;

function_ent
: FUNCTION
  {
    prompt++;
  }
  left_paren opt_arg_name_list right_paren newline '{' newlines
  {
    if(scope == LOCAL)
    {
      if(arg_sym_tab != 0)
        list_Destroy(arg_sym_tab);
      if(lsymtab != 0)
        list_Destroy(lsymtab);
      if(gsymtab != 0)
        list_Destroy(gsymtab);
      if(ssymtab != 0)
        list_Destroy(ssymtab);
      yyerror("function decl not allowed inside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    scope = LOCAL;
    psave = get_progoff ();
    lsave = lineno+loff;
    arg_sym_tab = $4;
    lsymtab = list_Create ();
    gsymtab = list_Create ();
    ssymtab = 0;
    function_setup1 (lsave, curr_file_name);
  }
  cstmts '}'
  {
    code (OP_DEF_FUNC_RET);
    code (STOP);
    $$ = function_setup2 (arg_sym_tab, lsymtab, gsymtab, ssymtab, get_progoff ());
    arg_sym_tab = 0;
    lsymtab = 0;
    gsymtab = 0;
    ssymtab = 0;
    if (rlist)
      list_DestroyAllNodes (rlist);
    scope = GLOBAL;
    prompt--;
  }
;

self: SELF {
  if(scope == GLOBAL)
  {
    yyerror("$self not allowed outside function");
    longjmp( *jmp_dec_buff (), 1 );
  }
}
;

 /* For each argument in a User-Function create a special 
    type of variable */
opt_arg_name_list
: { $$ = list_Create(); } /* empty () */
| NAME { $$ = arg_var_push(0, $1); }
| opt_arg_name_list ',' NAME { arg_var_push($1, $3); }
;

 /* List of local variables for a function */
gl_stmt
: LOCAL_DEC
  {
    if (scope == GLOBAL)
    {
      yyerror("local decl not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    dlocal = 1;
  }
  '(' opt_local_name_list ')' separator
  {
    $$ = get_progoff ();
    dlocal = 0;
  }
| GLOBAL_DEC
  {
    if (scope == GLOBAL)
    {
      yyerror("global decl not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    dglobal = 1;
  }
  '(' opt_local_name_list ')' separator
  {
    $$ = get_progoff ();
    dglobal = 0;
  }
;

opt_local_name_list
: NAME
  {
    $$ = get_progoff ();
    if (dlocal)
      local_var_push(lsymtab, $1);
    else
      global_var_push(gsymtab, $1);
  }
| opt_local_name_list ',' NAME
  {
    if (dlocal)
      local_var_push(lsymtab, $3);
    else
      global_var_push(gsymtab, $3);
  }
;


 /* As we accumulate args on the stack, keep count */
opt_arg_list
: /* empty */ { $$.count = 0; $$.off = get_progoff (); }
| texpr       { $$.count = 1; $$.off = $1; }
| texpr_list
  {
    $$.count = $1.count;
    $$.off = $1.off;
  }
;

texpr_list
: texpr_or_empty ',' texpr_or_empty  { $$.count = 2; $$.off = $1; }
| texpr_list ',' texpr_or_empty      { $$.count = $1.count + 1; $$.off = $1.off; }
;

texpr_or_empty
: /* empty */   { code(OP_PUSH_UNDEF); }
| texpr
;

 /* List of file- or function-static variables */
fstatic_stmt
: FSTATIC_DEC '(' opt_fstatic_name_list ')' separator { $$ = get_progoff (); }
;

opt_fstatic_name_list
: NAME
  {
    if(scope == GLOBAL)
    {
      if (class_scope>0)
      {
        btree_class_stat_var_push (class_stat_symtab, $1);
      }
      else
      {
        fstatic_var_push(curr_file_name, $1);
      }
    }
    else if (scope == LOCAL)
    {
      var_push(&ssymtab, $1);
    }
    $$ = 0;
  }
| opt_fstatic_name_list ',' NAME
  {
    if(scope == GLOBAL)
    {
      if (class_scope>0)
      {
        btree_class_stat_var_push (class_stat_symtab, $3);
      }
      else
      {
        fstatic_var_push(curr_file_name, $3);
      }
    }
    else if (scope == LOCAL)
    {
      var_push(&ssymtab, $3);
    }
    $$ = 0;
  }
;

 /*
  * Scripts - I
  */
rscript_dec
: RSCRIPT_DEC newline  { read_script=1; }
;

rscript_stmts
: rscript_dec R_SCRIPT  { read_script=0; $$ = $2; }
;

 /*
  * Scripts - II : CLASSDEF
  */
class_private_stmt
: CLASS_PRIV '(' opt_class_private_name_list ')' separator newlines
  { $$ = get_progoff (); }
;

opt_class_private_name_list
: NAME
  {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      fstatic_var_push(priv_class_name, $1);
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    $$ = 0;
  }
| opt_class_private_name_list ',' NAME
  {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      fstatic_var_push(priv_class_name, $3);
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    $$ = 0;
  }
;

class_public_stmt
: CLASS_PUBL '(' opt_class_public_name_list ')' separator newlines
  {
    $$ = get_progoff ();
  }
;

opt_class_public_name_list
: NAME
  {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      btree_class_publ_var_push (class_publ_symtab, $1);
    }
    else
    {
      yyerror("public decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    $$ = 0;
  }
| opt_class_public_name_list ',' NAME
  {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      btree_class_publ_var_push (class_publ_symtab, $3);
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    $$ = 0;
  }
;


class_inst_var_list
: /* empty () */                { $$ = list_Create(); }
| NAME                          { $$ = arg_var_push(0, $1); }
| class_inst_var_list ',' NAME  { arg_var_push($1, $3); }
;

class
: class_ent
  {
    $$ = code(OP_PUSH_VAR); codep($1);
  }
;

class_ent
: CLASSDEF
  {
    if(scope != GLOBAL)
    {
      if(class_args_symtab != 0)
        list_Destroy(class_args_symtab);
      if(class_glob_symtab != 0)
        list_Destroy(class_glob_symtab);
      if(class_publ_symtab != 0)
        btree_Destroy(class_publ_symtab);
      yyerror("classdef not allowed inside function declaration!\n");
      longjmp( *jmp_dec_buff (), 1 );
    }
    prompt++;
  }
  left_paren class_inst_var_list right_paren
  {
    read_script=1;
    scope = GLOBAL;
    class_args_symtab = $4;
    psave = get_progoff ();
    lsave = lineno+loff;
  }
  newlines R_SCRIPT
  {
    read_script=0;
    $$ = classdef_setup ($4, curr_file_name, $8, lsave);
    class_args_symtab=0;
    if (rlist)
      list_DestroyAllNodes (rlist);
    scope = GLOBAL;
    prompt--;
  }
;


%%

#include "rlab_parser_f.c"
