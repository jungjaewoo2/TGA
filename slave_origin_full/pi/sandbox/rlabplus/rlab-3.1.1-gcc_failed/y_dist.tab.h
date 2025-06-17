/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INTEGER = 258,
     NUMBER = 259,
     iNUMBER = 260,
     REOF = 261,
     ENDFILE = 262,
     NAME = 263,
     R_STRING = 264,
     HELP_NAME = 265,
     FNAME = 266,
     UNDEFINED = 267,
     SYS_CMD = 268,
     R_SCRIPT = 269,
     CLASS_PUBL = 270,
     CLASS_PRIV = 271,
     DO = 272,
     UNTIL = 273,
     WHILE = 274,
     IF = 275,
     ELSE = 276,
     QUIT = 277,
     FOR = 278,
     IN = 279,
     BREAK = 280,
     CONTINUE = 281,
     LOOPTHEN = 282,
     SWITCH = 283,
     SWITCH_CASE = 284,
     SWITCH_DEFAULT = 285,
     FUNCTION = 286,
     RETURN = 287,
     SEMIC = 288,
     Q_MARK = 289,
     LOCAL_DEC = 290,
     GLOBAL_DEC = 291,
     SELF = 292,
     HELP = 293,
     RFILE = 294,
     REQUIRE = 295,
     GST = 296,
     FSTATIC_DEC = 297,
     RSCRIPT_DEC = 298,
     CLASSDEF = 299,
     JNK = 300,
     OR = 301,
     AND = 302,
     NE = 303,
     EQ = 304,
     LE = 305,
     LT = 306,
     GE = 307,
     GT = 308,
     EL_DIV_BY = 309,
     EL_MUL_BY = 310,
     SUBFROM = 311,
     ADDTO = 312,
     LDIV_OP = 313,
     EL_LDIV_OP = 314,
     EL_DIV_OP = 315,
     EL_MUL_OP = 316,
     UNARY_PLUS = 317,
     UNARY_MINUS = 318,
     EL_POW_OP = 319,
     EL_TRANSPOSE = 320,
     TRANSPOSE = 321,
     DEC = 322,
     INC = 323,
     NOT = 324,
     R_OLIST = 325,
     L_OLIST = 326,
     RIGHT_LIST = 327,
     LEFT_LIST = 328
   };
#endif
/* Tokens.  */
#define INTEGER 258
#define NUMBER 259
#define iNUMBER 260
#define REOF 261
#define ENDFILE 262
#define NAME 263
#define R_STRING 264
#define HELP_NAME 265
#define FNAME 266
#define UNDEFINED 267
#define SYS_CMD 268
#define R_SCRIPT 269
#define CLASS_PUBL 270
#define CLASS_PRIV 271
#define DO 272
#define UNTIL 273
#define WHILE 274
#define IF 275
#define ELSE 276
#define QUIT 277
#define FOR 278
#define IN 279
#define BREAK 280
#define CONTINUE 281
#define LOOPTHEN 282
#define SWITCH 283
#define SWITCH_CASE 284
#define SWITCH_DEFAULT 285
#define FUNCTION 286
#define RETURN 287
#define SEMIC 288
#define Q_MARK 289
#define LOCAL_DEC 290
#define GLOBAL_DEC 291
#define SELF 292
#define HELP 293
#define RFILE 294
#define REQUIRE 295
#define GST 296
#define FSTATIC_DEC 297
#define RSCRIPT_DEC 298
#define CLASSDEF 299
#define JNK 300
#define OR 301
#define AND 302
#define NE 303
#define EQ 304
#define LE 305
#define LT 306
#define GE 307
#define GT 308
#define EL_DIV_BY 309
#define EL_MUL_BY 310
#define SUBFROM 311
#define ADDTO 312
#define LDIV_OP 313
#define EL_LDIV_OP 314
#define EL_DIV_OP 315
#define EL_MUL_OP 316
#define UNARY_PLUS 317
#define UNARY_MINUS 318
#define EL_POW_OP 319
#define EL_TRANSPOSE 320
#define TRANSPOSE 321
#define DEC 322
#define INC 323
#define NOT 324
#define R_OLIST 325
#define L_OLIST 326
#define RIGHT_LIST 327
#define LEFT_LIST 328



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{


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



} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
