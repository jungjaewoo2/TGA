/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison implementation for Yacc-like parsers in C
   
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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */


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




# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
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

/* Copy the second part of user declarations.  */



#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1866

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  90
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  95
/* YYNRULES -- Number of rules.  */
#define YYNRULES  233
/* YYNRULES -- Number of states.  */
#define YYNSTATES  419

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   328

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      89,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      85,    86,    63,    61,    47,    62,    84,    64,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    48,     2,
       2,    46,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    82,     2,    83,    71,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    87,     2,    88,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    65,    66,    67,    68,    69,    70,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    16,    19,    22,
      23,    25,    27,    29,    32,    35,    37,    40,    42,    45,
      47,    49,    51,    53,    55,    57,    59,    61,    63,    65,
      67,    69,    71,    73,    75,    77,    79,    80,    82,    84,
      86,    89,    92,    94,    97,    99,   101,   103,   105,   107,
     110,   112,   114,   116,   118,   121,   125,   136,   147,   158,
     166,   167,   180,   182,   183,   185,   187,   190,   192,   195,
     196,   201,   207,   212,   215,   216,   222,   224,   233,   242,
     248,   254,   256,   262,   264,   268,   272,   274,   276,   278,
     281,   284,   285,   287,   289,   291,   295,   299,   305,   309,
     318,   326,   334,   341,   349,   355,   359,   363,   365,   366,
     368,   372,   377,   379,   381,   385,   392,   398,   404,   409,
     413,   418,   423,   425,   427,   429,   431,   435,   437,   439,
     441,   443,   446,   449,   451,   454,   458,   464,   468,   469,
     471,   475,   479,   481,   485,   487,   489,   492,   495,   500,
     505,   510,   515,   519,   523,   527,   531,   535,   539,   543,
     547,   551,   555,   558,   561,   565,   569,   573,   577,   581,
     585,   589,   593,   597,   601,   605,   609,   612,   615,   618,
     620,   622,   624,   626,   628,   630,   632,   634,   636,   639,
     640,   642,   643,   645,   648,   650,   652,   653,   654,   666,
     668,   669,   671,   675,   676,   683,   684,   691,   693,   697,
     698,   700,   702,   706,   710,   711,   713,   719,   721,   725,
     728,   731,   738,   740,   744,   751,   753,   757,   758,   760,
     764,   766,   767,   768
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
      91,     0,    -1,    -1,    91,     6,    -1,    91,    93,    -1,
      91,    38,    -1,    91,    10,    -1,    91,    13,    -1,    91,
       1,    -1,    -1,     6,    -1,   158,    -1,     7,    -1,    94,
     158,    -1,    94,     6,    -1,    95,    -1,    94,    95,    -1,
       7,    -1,   136,    96,    -1,   110,    -1,   109,    -1,   122,
      -1,   123,    -1,   124,    -1,   111,    -1,   133,    -1,   134,
      -1,   128,    -1,   108,    -1,   165,    -1,   172,    -1,   178,
      -1,   176,    -1,    97,    -1,    98,    -1,    22,    -1,    -1,
      33,    -1,    34,    -1,    39,    -1,    39,    99,    -1,    40,
      99,    -1,    11,    -1,    99,    11,    -1,    85,    -1,    86,
      -1,    87,    -1,    88,    -1,    21,    -1,    27,    48,    -1,
     158,    -1,   105,    -1,    94,    -1,   106,    -1,   107,   106,
      -1,   102,   107,   103,    -1,   132,   129,   156,   108,   156,
     135,   104,   156,   110,   135,    -1,   132,   129,   156,   108,
     156,   135,   104,   156,   109,   135,    -1,   132,   129,   156,
     108,   156,   135,   104,   156,   108,   135,    -1,   132,   129,
     156,   108,   156,    92,   135,    -1,    -1,   113,   100,   137,
     101,   112,   156,    87,   156,   114,   156,    88,   135,    -1,
      28,    -1,    -1,   119,    -1,   115,    -1,   115,   119,    -1,
     118,    -1,   115,   118,    -1,    -1,    29,   117,   137,    48,
      -1,   116,   156,   107,   156,   133,    -1,   116,   156,   107,
     156,    -1,   116,   156,    -1,    -1,    30,    48,   156,   120,
     107,    -1,    17,    -1,   121,   156,   108,   156,   130,   129,
      96,   135,    -1,   121,   156,   108,   156,   131,   129,    96,
     135,    -1,   130,   129,   157,   108,   135,    -1,   125,   126,
     157,   108,   135,    -1,    23,    -1,   100,   153,    24,   127,
     101,    -1,   136,    -1,    32,   136,    96,    -1,   100,   136,
     101,    -1,    19,    -1,    18,    -1,    20,    -1,    25,    96,
      -1,    26,    96,    -1,    -1,   137,    -1,   138,    -1,   139,
      -1,    85,   136,    86,    -1,   137,    48,   137,    -1,   137,
      48,   137,    48,   137,    -1,   153,    46,   136,    -1,   137,
      82,   142,    33,   142,    83,    46,   136,    -1,   137,    82,
     142,    33,    83,    46,   136,    -1,   137,    82,    33,   142,
      83,    46,   136,    -1,   137,    82,   142,    83,    46,   136,
      -1,   137,    84,    82,   137,    83,    46,   136,    -1,   137,
      84,     8,    46,   136,    -1,   150,    46,   136,    -1,    82,
     141,    83,    -1,   140,    -1,    -1,   142,    -1,   141,    33,
     142,    -1,   141,    33,   155,   142,    -1,   137,    -1,   138,
      -1,   142,    47,   142,    -1,   137,    82,   142,    33,   142,
      83,    -1,   137,    82,   142,    33,    83,    -1,   137,    82,
      33,   142,    83,    -1,   137,    82,    33,    83,    -1,   137,
      82,    83,    -1,   137,    82,   142,    83,    -1,   137,    82,
      48,    83,    -1,   159,    -1,   181,    -1,   147,    -1,   148,
      -1,     8,    46,   136,    -1,   136,    -1,    33,    -1,   158,
      -1,    81,    -1,    81,    33,    -1,    81,   158,    -1,    80,
      -1,   144,    80,    -1,   145,   149,   146,    -1,   137,    84,
      82,   137,    83,    -1,   137,    84,     8,    -1,    -1,   143,
      -1,   149,   144,   143,    -1,    79,   151,    78,    -1,   153,
      -1,   151,    33,   153,    -1,   152,    -1,   153,    -1,   153,
      76,    -1,   153,    75,    -1,     8,    85,   169,    86,    -1,
     137,    85,   169,    86,    -1,   163,    85,   169,    86,    -1,
     148,    85,   169,    86,    -1,   137,    61,   137,    -1,   137,
      62,   137,    -1,   137,    63,   137,    -1,   137,    68,   137,
      -1,   137,    64,   137,    -1,   137,    65,   137,    -1,   137,
      67,   137,    -1,   137,    66,   137,    -1,   137,    71,   137,
      -1,   137,    72,   137,    -1,   137,    74,    -1,   137,    73,
      -1,   137,    56,   137,    -1,   137,    55,   137,    -1,   137,
      54,   137,    -1,   137,    53,   137,    -1,   137,    52,   137,
      -1,   137,    51,   137,    -1,   137,    50,   137,    -1,   137,
      49,   137,    -1,   137,    60,   137,    -1,   137,    59,   137,
      -1,   137,    58,   137,    -1,   137,    57,   137,    -1,    77,
     137,    -1,    62,   137,    -1,    61,   137,    -1,     4,    -1,
       3,    -1,     5,    -1,     9,    -1,   175,    -1,   154,    -1,
       8,    -1,    41,    -1,   158,    -1,   155,   158,    -1,    -1,
     158,    -1,    -1,   158,    -1,   157,   158,    -1,    89,    -1,
     160,    -1,    -1,    -1,    31,   161,   100,   164,   101,   156,
      87,   157,   162,   107,    88,    -1,    37,    -1,    -1,     8,
      -1,   164,    47,     8,    -1,    -1,    35,   166,    85,   168,
      86,    96,    -1,    -1,    36,   167,    85,   168,    86,    96,
      -1,     8,    -1,   168,    47,     8,    -1,    -1,   136,    -1,
     170,    -1,   171,    47,   171,    -1,   170,    47,   171,    -1,
      -1,   136,    -1,    42,    85,   173,    86,    96,    -1,     8,
      -1,   173,    47,     8,    -1,    43,   156,    -1,   174,    14,
      -1,    16,    85,   177,    86,    96,   157,    -1,     8,    -1,
     177,    47,     8,    -1,    15,    85,   179,    86,    96,   157,
      -1,     8,    -1,   179,    47,     8,    -1,    -1,     8,    -1,
     180,    47,     8,    -1,   182,    -1,    -1,    -1,    44,   183,
     100,   180,   101,   184,   157,    14,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   194,   194,   195,   196,   197,   198,   199,   200,   204,
     205,   212,   213,   214,   215,   227,   228,   232,   233,   263,
     264,   265,   266,   267,   268,   269,   270,   271,   272,   273,
     274,   275,   276,   277,   278,   279,   284,   285,   286,   290,
     291,   295,   299,   300,   303,   306,   309,   312,   316,   324,
     354,   355,   356,   360,   361,   365,   373,   380,   387,   397,
     412,   411,   428,   438,   439,   440,   441,   445,   446,   451,
     450,   465,   466,   474,   486,   485,   497,   505,   516,   531,
     558,   585,   595,   602,   606,   622,   630,   640,   649,   657,
     686,   706,   715,   716,   717,   721,   726,   727,   732,   734,
     735,   736,   737,   739,   740,   741,   748,   752,   756,   757,
     758,   759,   763,   764,   765,   770,   771,   772,   773,   774,
     775,   776,   779,   782,   789,   790,   795,   796,   800,   801,
     805,   806,   807,   811,   812,   816,   825,   826,   830,   831,
     832,   840,   849,   854,   871,   872,   873,   874,   875,   904,
     911,   917,   922,   923,   924,   925,   926,   927,   928,   929,
     930,   931,   932,   933,   934,   935,   936,   937,   938,   939,
     940,   941,   942,   943,   944,   945,   946,   947,   948,   952,
     953,   954,   955,   956,   963,   990,   994,   998,   999,  1003,
    1004,  1008,  1009,  1010,  1014,  1022,  1030,  1034,  1029,  1073,
    1085,  1086,  1087,  1093,  1092,  1107,  1106,  1123,  1131,  1143,
    1144,  1145,  1153,  1154,  1158,  1159,  1164,  1168,  1187,  1212,
    1216,  1223,  1228,  1241,  1257,  1264,  1277,  1294,  1295,  1296,
    1300,  1308,  1323,  1307
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INTEGER", "NUMBER", "iNUMBER", "REOF",
  "ENDFILE", "NAME", "R_STRING", "HELP_NAME", "FNAME", "UNDEFINED",
  "SYS_CMD", "R_SCRIPT", "CLASS_PUBL", "CLASS_PRIV", "DO", "UNTIL",
  "WHILE", "IF", "ELSE", "QUIT", "FOR", "IN", "BREAK", "CONTINUE",
  "LOOPTHEN", "SWITCH", "SWITCH_CASE", "SWITCH_DEFAULT", "FUNCTION",
  "RETURN", "SEMIC", "Q_MARK", "LOCAL_DEC", "GLOBAL_DEC", "SELF", "HELP",
  "RFILE", "REQUIRE", "GST", "FSTATIC_DEC", "RSCRIPT_DEC", "CLASSDEF",
  "JNK", "'='", "','", "':'", "OR", "AND", "NE", "EQ", "LE", "LT", "GE",
  "GT", "EL_DIV_BY", "EL_MUL_BY", "SUBFROM", "ADDTO", "'+'", "'-'", "'*'",
  "'/'", "LDIV_OP", "EL_LDIV_OP", "EL_DIV_OP", "EL_MUL_OP", "UNARY_PLUS",
  "UNARY_MINUS", "'^'", "EL_POW_OP", "EL_TRANSPOSE", "TRANSPOSE", "DEC",
  "INC", "NOT", "R_OLIST", "L_OLIST", "RIGHT_LIST", "LEFT_LIST", "'['",
  "']'", "'.'", "'('", "')'", "'{'", "'}'", "'\\n'", "$accept", "program",
  "reof", "line", "stmts", "stmt", "separator", "rfile_stmt", "req_stmt",
  "rnames", "left_paren", "right_paren", "left_brace", "right_brace",
  "else", "then", "cline", "cstmts", "c_stmt", "if_stmt", "if_stmt_open",
  "switch_stmt", "$@1", "switch_label", "switch_case_list",
  "switch_case_list1", "switch_case_open", "$@2", "switch_case",
  "switch_default_case", "$@3", "do", "do_stmt", "while_stmt", "for_stmt",
  "for", "for_spec", "for_vec_expr", "return_stmt", "cond", "while",
  "until", "if", "break", "continue", "end", "texpr", "expr", "vexpr",
  "assign", "mat_expr", "mat_list", "vec_expr", "list_el", "list_el_term",
  "left_list", "right_list", "list_expr", "list_member", "list_list",
  "olist", "olist_list", "value", "vid", "var", "nl", "newline",
  "newlines", "NL", "function", "function_ent", "$@4", "$@5", "self",
  "opt_arg_name_list", "gl_stmt", "$@6", "$@7", "opt_local_name_list",
  "opt_arg_list", "texpr_list", "texpr_or_empty", "fstatic_stmt",
  "opt_fstatic_name_list", "rscript_dec", "rscript_stmts",
  "class_private_stmt", "opt_class_private_name_list", "class_public_stmt",
  "opt_class_public_name_list", "class_inst_var_list", "class",
  "class_ent", "$@8", "$@9", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,    61,    44,    58,   301,
     302,   303,   304,   305,   306,   307,   308,   309,   310,   311,
     312,    43,    45,    42,    47,   313,   314,   315,   316,   317,
     318,    94,   319,   320,   321,   322,   323,   324,   325,   326,
     327,   328,    91,    93,    46,    40,    41,   123,   125,    10
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    90,    91,    91,    91,    91,    91,    91,    91,    92,
      92,    93,    93,    93,    93,    94,    94,    95,    95,    95,
      95,    95,    95,    95,    95,    95,    95,    95,    95,    95,
      95,    95,    95,    95,    95,    95,    96,    96,    96,    97,
      97,    98,    99,    99,   100,   101,   102,   103,   104,   105,
     106,   106,   106,   107,   107,   108,   109,   109,   109,   110,
     112,   111,   113,   114,   114,   114,   114,   115,   115,   117,
     116,   118,   118,   118,   120,   119,   121,   122,   122,   123,
     124,   125,   126,   127,   128,   129,   130,   131,   132,   133,
     134,   135,   136,   136,   136,   137,   138,   138,   139,   139,
     139,   139,   139,   139,   139,   139,   140,   137,   141,   141,
     141,   141,   142,   142,   142,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   143,   143,   144,   144,
     145,   145,   145,   146,   146,   147,   148,   148,   149,   149,
     149,   150,   151,   151,   137,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,   137,   152,
     152,   152,   152,   152,   153,   154,   154,   155,   155,   156,
     156,   157,   157,   157,   158,   159,   161,   162,   160,   163,
     164,   164,   164,   166,   165,   167,   165,   168,   168,   169,
     169,   169,   170,   170,   171,   171,   172,   173,   173,   174,
     175,   176,   177,   177,   178,   179,   179,   180,   180,   180,
     181,   183,   184,   182
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     2,     2,     2,     0,
       1,     1,     1,     2,     2,     1,     2,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     0,     1,     1,     1,
       2,     2,     1,     2,     1,     1,     1,     1,     1,     2,
       1,     1,     1,     1,     2,     3,    10,    10,    10,     7,
       0,    12,     1,     0,     1,     1,     2,     1,     2,     0,
       4,     5,     4,     2,     0,     5,     1,     8,     8,     5,
       5,     1,     5,     1,     3,     3,     1,     1,     1,     2,
       2,     0,     1,     1,     1,     3,     3,     5,     3,     8,
       7,     7,     6,     7,     5,     3,     3,     1,     0,     1,
       3,     4,     1,     1,     3,     6,     5,     5,     4,     3,
       4,     4,     1,     1,     1,     1,     3,     1,     1,     1,
       1,     2,     2,     1,     2,     3,     5,     3,     0,     1,
       3,     3,     1,     3,     1,     1,     2,     2,     4,     4,
       4,     4,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     2,     2,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     2,     2,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     2,     0,
       1,     0,     1,     2,     1,     1,     0,     0,    11,     1,
       0,     1,     3,     0,     6,     0,     6,     1,     3,     0,
       1,     1,     3,     3,     0,     1,     5,     1,     3,     2,
       2,     6,     1,     3,     6,     1,     3,     0,     1,     3,
       1,     0,     0,     8
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     8,   180,   179,   181,     3,    12,   185,
     182,     6,     7,     0,     0,    76,    86,    88,    35,    81,
      36,    36,    62,   196,     0,   203,   205,   199,     5,    39,
       0,   186,     0,   189,   231,     0,     0,     0,     0,   130,
     108,     0,    46,   194,     4,     0,    15,    33,    34,     0,
      28,    20,    19,    24,     0,   189,    21,    22,    23,     0,
      27,     0,     0,    25,    26,    36,    92,    93,    94,   107,
     138,   124,   125,     0,   144,   145,   184,    11,   122,   195,
       0,    29,    30,     0,   183,    32,    31,   123,   230,   209,
       0,     0,    37,    38,    89,    90,     0,    36,     0,     0,
      42,    40,    41,     0,   219,   190,     0,   178,   145,   177,
     176,   185,     0,   142,   131,   132,   112,   113,     0,   109,
       0,    14,    17,    16,    13,     0,    52,    51,    53,     0,
      50,    44,     0,     0,     0,   191,     0,   191,   189,    18,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   163,   162,     0,     0,   209,   185,   127,
     139,     0,   209,     0,     0,   147,   146,   209,   220,   210,
       0,   211,     0,   225,     0,   222,     0,   200,    84,     0,
       0,    43,   217,     0,   227,     0,     0,     0,   141,     0,
     106,     0,    95,    49,    47,    55,    54,     0,   189,     0,
       0,   192,     0,     0,     0,    96,   171,   170,   169,   168,
     167,   166,   165,   164,   175,   174,   173,   172,   152,   153,
     154,   156,   157,   159,   158,   155,   160,   161,     0,     0,
     119,     0,   137,     0,     0,     0,   128,   133,     0,   135,
     129,     0,   105,    98,     0,   148,   214,   214,     0,    36,
       0,    36,   201,     0,   207,     0,     0,     0,    36,   228,
       0,     0,     0,   137,     0,   143,   110,     0,   187,   114,
      45,    60,     0,     0,    91,   193,    85,    91,   189,     0,
     118,     0,   121,     0,   120,     0,     0,   149,   126,   134,
     140,   151,   150,   215,   213,   212,   226,   191,   223,   191,
       0,   189,     0,    36,    36,   218,   216,     0,   232,     0,
       0,   120,     0,   111,   188,   189,    87,     0,     0,     0,
      83,    80,    79,     9,    97,   117,   116,     0,     0,   104,
     136,   224,   221,   202,     0,   208,   204,   206,   229,   191,
     117,   116,     0,   136,     0,    36,    36,    82,    10,    91,
       0,     0,     0,   115,   102,     0,   191,     0,   115,   189,
      91,    91,    59,    48,   189,   101,   100,     0,   103,   197,
     233,    63,    77,    78,     0,    99,     0,    69,     0,   189,
      65,   189,    67,    64,    91,    91,    91,     0,     0,   189,
       0,    68,    66,    73,    58,    57,    56,   198,     0,    74,
      91,   189,    70,     0,    61,    72,    50,    75,    71
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,   359,    44,   126,    46,    94,    47,    48,   101,
     136,   281,    49,   205,   374,   127,   128,   129,    50,    51,
      52,    53,   325,    54,   389,   390,   391,   398,   392,   393,
     413,    55,    56,    57,    58,    59,   135,   329,    60,   137,
      61,   328,    62,    63,    64,   331,    65,    66,    67,    68,
      69,   118,   119,   170,   248,    70,   249,    71,    72,   171,
      73,   112,    74,   108,    76,   277,   104,   210,   105,    78,
      79,    96,   386,    80,   263,    81,    98,    99,   265,   180,
     181,   182,    82,   193,    83,    84,    85,   186,    86,   184,
     270,    87,    88,   106,   349
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -285
static const yytype_int16 yypact[] =
{
    -285,   530,  -285,  -285,  -285,  -285,  -285,  -285,  -285,   -55,
    -285,  -285,  -285,   -47,   -42,  -285,  -285,  -285,  -285,  -285,
      74,    74,  -285,  -285,  1198,  -285,  -285,  -285,  -285,    39,
      39,  -285,    10,    14,  -285,  1413,  1413,  1413,    98,   -22,
    1413,  1198,  -285,  -285,  -285,   791,  -285,  -285,  -285,   878,
    -285,  -285,  -285,  -285,    37,    14,  -285,  -285,  -285,    37,
    -285,    37,    37,  -285,  -285,    74,  1451,  -285,  -285,  -285,
    1241,  -285,    38,    52,  -285,   -21,  -285,  -285,  -285,  -285,
      56,  -285,  -285,   102,  -285,  -285,  -285,  -285,  -285,  1137,
     116,   134,  -285,  -285,  -285,  -285,    37,    74,    65,    67,
    -285,   142,   142,   150,  -285,  -285,    37,    75,   -34,    75,
      20,  -285,   -17,  -285,  -285,  -285,  1489,  -285,    -6,   115,
      78,  -285,  -285,  -285,  -285,   113,   996,  -285,  -285,   617,
    -285,  -285,  1413,    80,    98,    14,  1198,    14,    14,  -285,
    1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,
    1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,  1413,
    1413,  1413,  1413,  -285,  -285,  1043,     7,  1137,   -33,  -285,
    -285,    -9,  1137,  1198,  1198,  -285,  -285,  1137,  -285,   118,
      82,   122,   123,  -285,   -27,  -285,    -2,   166,  -285,   169,
     169,  -285,  -285,     4,   174,  1090,     9,    98,  -285,   925,
    -285,  1413,  -285,  -285,  -285,  -285,  -285,  1526,    14,   151,
     -54,  -285,   100,   -54,    80,  1565,  1713,  1748,  1781,  1781,
     398,   398,   398,   398,   127,   127,   127,   127,   167,   167,
      75,    75,    75,    75,    75,    75,    75,    75,  1284,   119,
    -285,   -15,   158,  1413,   111,  1198,  -285,  -285,  1151,  -285,
    -285,   120,  -285,  -285,   124,  -285,  1198,  1198,   197,    74,
     200,    74,  -285,    13,  -285,    29,    31,   205,    74,  -285,
      45,  1284,   -11,  -285,  1413,  -285,   115,   925,  -285,  -285,
    -285,  -285,    95,  1198,  -285,  -285,  -285,  -285,    14,  1413,
    -285,    26,  -285,  1327,   168,  1198,   443,  -285,  -285,  -285,
    -285,  -285,  -285,  -285,  -285,  -285,  -285,    14,  -285,    14,
     207,    14,   208,    74,    74,  -285,  -285,   209,  -285,    47,
    1370,  -285,  1640,   115,  -285,    14,  -285,    37,    37,   100,
    -285,  -285,  -285,    25,  1677,   173,   176,    49,  1198,  -285,
     177,    14,    14,  -285,   138,  -285,  -285,  -285,  -285,    14,
    -285,  -285,    53,  -285,   139,    74,    74,  -285,  -285,  -285,
     206,  1198,  1198,   182,  -285,  1198,    14,    -8,  -285,    14,
    -285,  -285,  -285,  -285,    14,  -285,  -285,  1198,  -285,    14,
    -285,    89,  -285,  -285,   -13,  -285,   878,  -285,   181,    14,
      89,    14,  -285,  -285,  -285,  -285,  -285,   704,  1413,    14,
     148,  -285,  -285,   878,  -285,  -285,  -285,  -285,  1603,  -285,
    -285,   878,  -285,   878,  -285,   217,  -285,   878,  -285
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -285,  -285,  -285,  -285,   244,   -41,     0,  -285,  -285,   216,
     -40,  -200,  -285,  -285,  -285,  -285,  -128,  -259,  -131,  -130,
    -126,  -285,  -285,  -285,  -285,  -285,  -285,  -285,  -142,  -125,
    -285,  -285,  -285,  -285,  -285,  -285,  -285,  -285,  -285,   -57,
     -20,  -285,  -285,  -155,  -285,  -284,   114,   286,   -14,  -285,
    -285,  -285,  -137,    16,  -285,  -285,  -285,  -285,  -285,  -285,
    -285,  -285,  -285,    -1,  -285,  -285,   -45,  -129,     8,  -285,
    -285,  -285,  -285,  -285,  -285,  -285,  -285,  -285,    77,  -138,
    -285,  -136,  -285,  -285,  -285,  -285,  -285,  -285,  -285,  -285,
    -285,  -285,  -285,  -285,  -285
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -216
static const yytype_int16 yytable[] =
{
      75,   206,   208,   332,   123,   138,   380,    17,   213,    77,
     133,   114,   286,   245,   132,   242,   197,   273,   293,   134,
     258,    95,   320,    75,   246,   174,   117,   199,   241,   244,
      89,   358,   201,    42,   251,    43,   201,   113,    90,   254,
      75,   175,   176,    91,    75,   260,   -91,   115,    75,   360,
     100,   267,    89,   124,   175,   176,   187,   130,   272,   259,
     310,   198,   276,   311,   279,   139,   194,    43,   294,    75,
     318,   247,   321,   201,    42,   372,   312,   200,   312,   284,
      43,    43,   287,   288,   261,   123,   382,   383,    75,   243,
     268,   274,   317,   214,   201,   103,   201,   188,   173,   280,
     201,   291,   195,    43,   196,   167,   111,    92,    93,   335,
     404,   405,   406,   326,    16,   313,   178,   314,   387,   388,
     304,   305,   131,   172,   183,    75,   414,   397,    75,   357,
     350,   280,   363,   209,   319,    75,   368,   130,    97,    31,
     323,   177,   185,   211,   411,   211,   161,   162,   163,   164,
     189,   117,   190,   191,   417,   120,   337,   195,   192,   196,
     167,   203,   201,   282,   202,  -215,    75,    42,   255,   256,
     257,    75,    75,    75,   262,   283,    75,   264,   341,   250,
     342,   117,   269,   352,   169,   117,   280,   117,   153,   154,
     155,   156,   157,   158,   159,   160,   275,   297,   161,   162,
     163,   164,   292,   179,   295,   306,   301,   278,   308,   195,
     302,   196,   167,   315,   338,   343,   345,   348,   285,   361,
     367,   285,   362,   365,   117,   366,   369,   373,   377,   399,
     155,   156,   157,   158,   159,   160,   410,   379,   161,   162,
     163,   164,    20,   333,    75,    45,   102,    75,   401,   195,
     212,   196,   167,   394,   395,    75,    75,   117,   396,   307,
     418,   309,   327,   117,   300,   402,   344,   266,   316,   206,
     355,   356,     0,     0,     0,     0,     0,     0,     0,   117,
     354,   179,    75,   206,     0,   324,   179,   252,   253,   206,
       0,   179,     0,     0,    75,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   117,     0,     0,     0,
       0,     0,     0,   346,   347,   211,     0,   211,     0,     0,
       0,   107,   109,   110,   381,     0,   116,     0,     0,   384,
       0,     0,     0,     0,     0,     0,     0,    75,     0,     0,
       0,     0,     0,     0,   400,     0,   403,     0,     0,   285,
     285,     0,     0,     0,   409,   370,   371,   211,     0,   298,
      75,    75,   169,     0,    75,     0,   415,     0,     0,     0,
     303,   303,     0,     0,   211,   285,    75,     0,     0,     0,
       0,     0,     0,     0,     0,    75,     0,   285,     0,     0,
       0,     0,     0,     0,   130,     0,    75,   330,     0,     0,
       0,     0,    75,     0,     0,   130,     0,     0,     0,   339,
      75,   130,    75,     0,     0,     0,    75,     0,   207,   416,
       0,   130,     0,     0,     0,   130,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,     0,
       0,   116,   364,     0,     0,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,     0,     0,   161,
     162,   163,   164,     0,     0,   375,   376,     0,     0,   378,
     195,   116,   196,   167,     0,   116,     0,   116,     0,     0,
       0,   385,   141,   142,   143,   144,   145,   146,   147,   148,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,     0,     0,   161,   162,   163,   164,     0,     0,
       0,     0,     0,     0,   116,   195,   340,   196,   167,   296,
       2,     3,     0,     4,     5,     6,     7,     8,     9,    10,
      11,     0,     0,    12,     0,    13,    14,    15,     0,    16,
      17,     0,    18,    19,     0,    20,    21,   116,    22,     0,
     322,    23,    24,   116,     0,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,   334,     0,     0,     0,   116,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    35,    36,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   116,    37,     0,    38,
       0,    39,    40,     0,     0,    41,     0,    42,     0,    43,
       4,     5,     6,     0,   122,     9,    10,     0,     0,     0,
       0,     0,    13,    14,    15,     0,    16,    17,     0,    18,
      19,     0,    20,    21,   125,    22,     0,     0,    23,    24,
       0,     0,    25,    26,    27,     0,    29,    30,    31,    32,
      33,    34,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    35,    36,
       0,     0,     0,     0,   408,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    37,     0,    38,     0,    39,    40,
       0,     0,    41,     0,    42,   204,    43,     4,     5,     6,
       0,   122,     9,    10,     0,     0,     0,     0,     0,    13,
      14,    15,     0,    16,    17,     0,    18,    19,     0,    20,
      21,   125,    22,     0,     0,    23,    24,     0,     0,    25,
      26,    27,     0,    29,    30,    31,    32,    33,    34,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    35,    36,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    37,     0,    38,     0,    39,    40,     0,     0,    41,
       0,    42,   407,    43,     4,     5,     6,   121,   122,     9,
      10,     0,     0,     0,     0,     0,    13,    14,    15,     0,
      16,    17,     0,    18,    19,     0,    20,    21,     0,    22,
       0,     0,    23,    24,     0,     0,    25,    26,    27,     0,
      29,    30,    31,    32,    33,    34,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    35,    36,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    37,     0,
      38,     0,    39,    40,     0,     0,    41,     0,    42,     0,
      43,     4,     5,     6,     0,   122,     9,    10,     0,     0,
       0,     0,     0,    13,    14,    15,     0,    16,    17,     0,
      18,    19,     0,    20,    21,   125,    22,     0,     0,    23,
      24,     0,     0,    25,    26,    27,     0,    29,    30,    31,
      32,    33,    34,     0,     0,     0,     0,     0,     4,     5,
       6,     0,     0,     9,    10,     0,     0,     0,     0,    35,
      36,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    37,    23,    38,     0,    39,
      40,     0,    27,    41,     0,    42,    31,    43,    33,    34,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    35,    36,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     4,
       5,     6,    37,   122,     9,    10,    39,    40,     0,     0,
      41,    13,    14,    15,    43,    16,    17,     0,    18,    19,
       0,    20,    21,     0,    22,     0,     0,    23,    24,     0,
       0,    25,    26,    27,     0,    29,    30,    31,    32,    33,
      34,     0,     0,     0,     0,     0,     4,     5,     6,     0,
       0,     9,    10,     0,     0,     0,     0,    35,    36,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,    23,    38,   238,    39,    40,     0,
      27,    41,     0,    42,    31,     0,    33,    34,     0,     0,
       0,   239,     0,     4,     5,     6,     0,     0,     9,    10,
       0,     0,     0,     0,    35,    36,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      37,    23,     0,   271,    39,    40,   240,    27,    41,     0,
       0,    31,     0,    33,    34,     0,     0,     0,   239,     0,
       4,     5,     6,     0,     0,     9,    10,     0,     0,     0,
       0,    35,    36,     0,     4,     5,     6,     0,     0,   168,
      10,     0,     0,     0,     0,     0,     0,    37,    23,     0,
       0,    39,    40,   240,    27,    41,     0,     0,    31,     0,
      33,    34,    23,     0,  -214,     0,     0,     0,    27,     0,
       0,     0,    31,     0,    33,    34,     0,     0,    35,    36,
       0,     4,     5,     6,     0,     0,     9,    10,     0,     0,
       0,     0,    35,    36,    37,     0,    38,     0,    39,    40,
       0,     0,    41,     0,     0,     0,     0,     0,    37,    23,
      38,   299,    39,    40,     0,    27,    41,     0,     0,    31,
       0,    33,    34,     0,     4,     5,     6,     0,     0,   168,
      10,     0,     0,     0,     0,     0,     0,     0,     0,    35,
      36,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    23,     0,     0,    37,     0,    38,    27,    39,
      40,     0,    31,    41,    33,    34,     0,     4,     5,     6,
       0,     0,     9,    10,     0,     0,     0,     0,     0,     0,
       0,     0,    35,    36,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    23,     0,     0,    37,     0,
      38,    27,    39,    40,     0,    31,    41,    33,    34,     0,
       4,     5,     6,     0,     0,     9,    10,     0,     0,     0,
       0,     0,     0,     0,     0,    35,    36,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    23,     0,
       0,    37,     0,     0,    27,    39,    40,   290,    31,    41,
      33,    34,     0,     4,     5,     6,     0,     0,     9,    10,
       0,     0,     0,     0,     0,     0,     0,     0,    35,    36,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    23,     0,     0,    37,     0,     0,    27,    39,    40,
     336,    31,    41,    33,    34,     0,     4,     5,     6,     0,
       0,     9,    10,     0,     0,     0,     0,     0,     0,     0,
       0,    35,    36,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    23,     0,     0,    37,     0,     0,
      27,    39,    40,   351,    31,    41,    33,    34,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    35,    36,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      37,     0,     0,     0,    39,    40,     0,     0,    41,   140,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
       0,     0,   161,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,   165,     0,   166,   167,   140,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,     0,     0,
     161,   162,   163,   164,     0,     0,     0,     0,     0,     0,
       0,   195,     0,   196,   167,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,     0,     0,   161,   162,   163,
     164,     0,     0,     0,     0,     0,     0,     0,   195,     0,
     196,   167,   280,   289,   141,   142,   143,   144,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,     0,     0,   161,   162,   163,   164,
       0,     0,     0,     0,     0,     0,     0,   195,     0,   196,
     167,   412,   141,   142,   143,   144,   145,   146,   147,   148,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,     0,     0,   161,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,   195,     0,   196,   167,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,     0,
       0,   161,   162,   163,   164,     0,     0,     0,     0,     0,
       0,     0,   195,   353,   196,   167,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,     0,     0,   161,   162,
     163,   164,     0,     0,     0,     0,     0,     0,     0,   195,
       0,   196,   167,   142,   143,   144,   145,   146,   147,   148,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,     0,     0,   161,   162,   163,   164,     0,     0,
       0,     0,     0,     0,     0,   195,     0,   196,   167,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,     0,     0,   161,
     162,   163,   164,     0,     0,     0,     0,     0,     0,     0,
     195,     0,   196,   167,   145,   146,   147,   148,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
       0,     0,   161,   162,   163,   164,     0,     0,     0,     0,
       0,     0,     0,   195,     0,   196,   167
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-285)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
       1,   129,   133,   287,    45,    62,    14,    20,   137,     1,
      55,    33,   212,    46,    54,     8,    33,     8,    33,    59,
      47,    21,    33,    24,    33,    46,    40,    33,   165,   167,
      85,     6,    47,    87,   172,    89,    47,    38,    85,   177,
      41,    75,    76,    85,    45,    47,    21,    39,    49,   333,
      11,    47,    85,    45,    75,    76,    96,    49,   195,    86,
      47,    78,   199,   263,   201,    65,   106,    89,    83,    70,
     270,    80,    83,    47,    87,   359,    47,    83,    47,   210,
      89,    89,   213,   214,    86,   126,   370,   371,    89,    82,
      86,    82,    47,   138,    47,    85,    47,    97,    46,    86,
      47,   238,    82,    89,    84,    85,     8,    33,    34,    83,
     394,   395,   396,    18,    19,    86,    14,    86,    29,    30,
     256,   257,    85,    85,     8,   126,   410,   386,   129,   329,
      83,    86,    83,   134,   271,   136,    83,   129,    24,    41,
     277,    85,     8,   135,   403,   137,    71,    72,    73,    74,
      85,   165,    85,    11,   413,    41,   293,    82,     8,    84,
      85,    48,    47,   208,    86,    47,   167,    87,    86,    47,
      47,   172,   173,   174,     8,    24,   177,     8,   307,   171,
     309,   195,     8,   320,    70,   199,    86,   201,    61,    62,
      63,    64,    65,    66,    67,    68,   197,    86,    71,    72,
      73,    74,    83,    89,    46,     8,    86,   199,     8,    82,
      86,    84,    85,     8,    46,     8,     8,     8,   210,    46,
     349,   213,    46,    46,   238,    87,    87,    21,    46,    48,
      63,    64,    65,    66,    67,    68,    88,   366,    71,    72,
      73,    74,    25,   288,   245,     1,    30,   248,   390,    82,
     136,    84,    85,   384,   384,   256,   257,   271,   384,   259,
     415,   261,   282,   277,   248,   390,   311,   190,   268,   397,
     327,   328,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   293,
     325,   167,   283,   411,    -1,   277,   172,   173,   174,   417,
      -1,   177,    -1,    -1,   295,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   320,    -1,    -1,    -1,
      -1,    -1,    -1,   313,   314,   307,    -1,   309,    -1,    -1,
      -1,    35,    36,    37,   369,    -1,    40,    -1,    -1,   374,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   338,    -1,    -1,
      -1,    -1,    -1,    -1,   389,    -1,   391,    -1,    -1,   341,
     342,    -1,    -1,    -1,   399,   355,   356,   349,    -1,   245,
     361,   362,   248,    -1,   365,    -1,   411,    -1,    -1,    -1,
     256,   257,    -1,    -1,   366,   367,   377,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   386,    -1,   379,    -1,    -1,
      -1,    -1,    -1,    -1,   386,    -1,   397,   283,    -1,    -1,
      -1,    -1,   403,    -1,    -1,   397,    -1,    -1,    -1,   295,
     411,   403,   413,    -1,    -1,    -1,   417,    -1,   132,   411,
      -1,   413,    -1,    -1,    -1,   417,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,   162,    -1,
      -1,   165,   338,    -1,    -1,    57,    58,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    -1,    71,
      72,    73,    74,    -1,    -1,   361,   362,    -1,    -1,   365,
      82,   195,    84,    85,    -1,   199,    -1,   201,    -1,    -1,
      -1,   377,    49,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    -1,    -1,    71,    72,    73,    74,    -1,    -1,
      -1,    -1,    -1,    -1,   238,    82,    83,    84,    85,   243,
       0,     1,    -1,     3,     4,     5,     6,     7,     8,     9,
      10,    -1,    -1,    13,    -1,    15,    16,    17,    -1,    19,
      20,    -1,    22,    23,    -1,    25,    26,   271,    28,    -1,
     274,    31,    32,   277,    -1,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,   289,    -1,    -1,    -1,   293,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    61,    62,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   320,    77,    -1,    79,
      -1,    81,    82,    -1,    -1,    85,    -1,    87,    -1,    89,
       3,     4,     5,    -1,     7,     8,     9,    -1,    -1,    -1,
      -1,    -1,    15,    16,    17,    -1,    19,    20,    -1,    22,
      23,    -1,    25,    26,    27,    28,    -1,    -1,    31,    32,
      -1,    -1,    35,    36,    37,    -1,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    61,    62,
      -1,    -1,    -1,    -1,   398,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    77,    -1,    79,    -1,    81,    82,
      -1,    -1,    85,    -1,    87,    88,    89,     3,     4,     5,
      -1,     7,     8,     9,    -1,    -1,    -1,    -1,    -1,    15,
      16,    17,    -1,    19,    20,    -1,    22,    23,    -1,    25,
      26,    27,    28,    -1,    -1,    31,    32,    -1,    -1,    35,
      36,    37,    -1,    39,    40,    41,    42,    43,    44,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    61,    62,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    77,    -1,    79,    -1,    81,    82,    -1,    -1,    85,
      -1,    87,    88,    89,     3,     4,     5,     6,     7,     8,
       9,    -1,    -1,    -1,    -1,    -1,    15,    16,    17,    -1,
      19,    20,    -1,    22,    23,    -1,    25,    26,    -1,    28,
      -1,    -1,    31,    32,    -1,    -1,    35,    36,    37,    -1,
      39,    40,    41,    42,    43,    44,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    61,    62,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    77,    -1,
      79,    -1,    81,    82,    -1,    -1,    85,    -1,    87,    -1,
      89,     3,     4,     5,    -1,     7,     8,     9,    -1,    -1,
      -1,    -1,    -1,    15,    16,    17,    -1,    19,    20,    -1,
      22,    23,    -1,    25,    26,    27,    28,    -1,    -1,    31,
      32,    -1,    -1,    35,    36,    37,    -1,    39,    40,    41,
      42,    43,    44,    -1,    -1,    -1,    -1,    -1,     3,     4,
       5,    -1,    -1,     8,     9,    -1,    -1,    -1,    -1,    61,
      62,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    77,    31,    79,    -1,    81,
      82,    -1,    37,    85,    -1,    87,    41,    89,    43,    44,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    61,    62,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,     3,
       4,     5,    77,     7,     8,     9,    81,    82,    -1,    -1,
      85,    15,    16,    17,    89,    19,    20,    -1,    22,    23,
      -1,    25,    26,    -1,    28,    -1,    -1,    31,    32,    -1,
      -1,    35,    36,    37,    -1,    39,    40,    41,    42,    43,
      44,    -1,    -1,    -1,    -1,    -1,     3,     4,     5,    -1,
      -1,     8,     9,    -1,    -1,    -1,    -1,    61,    62,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    77,    31,    79,    33,    81,    82,    -1,
      37,    85,    -1,    87,    41,    -1,    43,    44,    -1,    -1,
      -1,    48,    -1,     3,     4,     5,    -1,    -1,     8,     9,
      -1,    -1,    -1,    -1,    61,    62,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      77,    31,    -1,    33,    81,    82,    83,    37,    85,    -1,
      -1,    41,    -1,    43,    44,    -1,    -1,    -1,    48,    -1,
       3,     4,     5,    -1,    -1,     8,     9,    -1,    -1,    -1,
      -1,    61,    62,    -1,     3,     4,     5,    -1,    -1,     8,
       9,    -1,    -1,    -1,    -1,    -1,    -1,    77,    31,    -1,
      -1,    81,    82,    83,    37,    85,    -1,    -1,    41,    -1,
      43,    44,    31,    -1,    47,    -1,    -1,    -1,    37,    -1,
      -1,    -1,    41,    -1,    43,    44,    -1,    -1,    61,    62,
      -1,     3,     4,     5,    -1,    -1,     8,     9,    -1,    -1,
      -1,    -1,    61,    62,    77,    -1,    79,    -1,    81,    82,
      -1,    -1,    85,    -1,    -1,    -1,    -1,    -1,    77,    31,
      79,    80,    81,    82,    -1,    37,    85,    -1,    -1,    41,
      -1,    43,    44,    -1,     3,     4,     5,    -1,    -1,     8,
       9,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    61,
      62,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    31,    -1,    -1,    77,    -1,    79,    37,    81,
      82,    -1,    41,    85,    43,    44,    -1,     3,     4,     5,
      -1,    -1,     8,     9,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    61,    62,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    31,    -1,    -1,    77,    -1,
      79,    37,    81,    82,    -1,    41,    85,    43,    44,    -1,
       3,     4,     5,    -1,    -1,     8,     9,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    61,    62,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    31,    -1,
      -1,    77,    -1,    -1,    37,    81,    82,    83,    41,    85,
      43,    44,    -1,     3,     4,     5,    -1,    -1,     8,     9,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    61,    62,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    31,    -1,    -1,    77,    -1,    -1,    37,    81,    82,
      83,    41,    85,    43,    44,    -1,     3,     4,     5,    -1,
      -1,     8,     9,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    61,    62,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    31,    -1,    -1,    77,    -1,    -1,
      37,    81,    82,    83,    41,    85,    43,    44,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    61,    62,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      77,    -1,    -1,    -1,    81,    82,    -1,    -1,    85,    48,
      49,    50,    51,    52,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      -1,    -1,    71,    72,    73,    74,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    82,    -1,    84,    85,    48,    49,    50,
      51,    52,    53,    54,    55,    56,    57,    58,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    -1,    -1,
      71,    72,    73,    74,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    82,    -1,    84,    85,    49,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    -1,    -1,    71,    72,    73,
      74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    82,    -1,
      84,    85,    86,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    -1,    -1,    71,    72,    73,    74,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    82,    -1,    84,
      85,    48,    49,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    -1,    -1,    71,    72,    73,    74,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    82,    -1,    84,    85,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    -1,
      -1,    71,    72,    73,    74,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    82,    83,    84,    85,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    64,    65,    66,    67,    68,    -1,    -1,    71,    72,
      73,    74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    82,
      -1,    84,    85,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    -1,    -1,    71,    72,    73,    74,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    82,    -1,    84,    85,    51,
      52,    53,    54,    55,    56,    57,    58,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    -1,    71,
      72,    73,    74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      82,    -1,    84,    85,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      -1,    -1,    71,    72,    73,    74,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    82,    -1,    84,    85
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    91,     0,     1,     3,     4,     5,     6,     7,     8,
       9,    10,    13,    15,    16,    17,    19,    20,    22,    23,
      25,    26,    28,    31,    32,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    61,    62,    77,    79,    81,
      82,    85,    87,    89,    93,    94,    95,    97,    98,   102,
     108,   109,   110,   111,   113,   121,   122,   123,   124,   125,
     128,   130,   132,   133,   134,   136,   137,   138,   139,   140,
     145,   147,   148,   150,   152,   153,   154,   158,   159,   160,
     163,   165,   172,   174,   175,   176,   178,   181,   182,    85,
      85,    85,    33,    34,    96,    96,   161,   136,   166,   167,
      11,    99,    99,    85,   156,   158,   183,   137,   153,   137,
     137,     8,   151,   153,    33,   158,   137,   138,   141,   142,
     136,     6,     7,    95,   158,    27,    94,   105,   106,   107,
     158,    85,   100,   156,   100,   126,   100,   129,   129,    96,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    71,    72,    73,    74,    82,    84,    85,     8,   136,
     143,   149,    85,    46,    46,    75,    76,    85,    14,   136,
     169,   170,   171,     8,   179,     8,   177,   100,    96,    85,
      85,    11,     8,   173,   100,    82,    84,    33,    78,    33,
      83,    47,    86,    48,    88,   103,   106,   137,   108,   153,
     157,   158,   136,   157,   156,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,    33,    48,
      83,   142,     8,    82,   169,    46,    33,    80,   144,   146,
     158,   169,   136,   136,   169,    86,    47,    47,    47,    86,
      47,    86,     8,   164,     8,   168,   168,    47,    86,     8,
     180,    33,   142,     8,    82,   153,   142,   155,   158,   142,
      86,   101,   156,    24,   108,   158,   101,   108,   108,    48,
      83,   142,    83,    33,    83,    46,   137,    86,   136,    80,
     143,    86,    86,   136,   171,   171,     8,    96,     8,    96,
      47,   101,    47,    86,    86,     8,    96,    47,   101,   142,
      33,    83,   137,   142,   158,   112,    18,   130,   131,   127,
     136,   135,   135,   156,   137,    83,    83,   142,    46,   136,
      83,   157,   157,     8,   156,     8,    96,    96,     8,   184,
      83,    83,   142,    83,   156,   129,   129,   101,     6,    92,
     135,    46,    46,    83,   136,    46,    87,   157,    83,    87,
      96,    96,   135,    21,   104,   136,   136,    46,   136,   157,
      14,   156,   135,   135,   156,   136,   162,    29,    30,   114,
     115,   116,   118,   119,   108,   109,   110,   107,   117,    48,
     156,   118,   119,   156,   135,   135,   135,    88,   137,   156,
      88,   107,    48,   120,   135,   156,   158,   107,   133
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

    { (yyval.n_int) = 0; }
    break;

  case 3:

    { code(OP_ENDOFFILE); return(0); }
    break;

  case 4:

    { code(STOP); return(1); }
    break;

  case 5:

    { code(OP_HELP); code(STOP); return (1); }
    break;

  case 6:

    { (yyval.n_int) = 0; code(OP_HELP_NAME); codep ((yyvsp[(2) - (2)].name)); code(STOP); return (1); }
    break;

  case 7:

    { code(OP_SYS_CMD); codep((yyvsp[(2) - (2)].name)); code(STOP); return (1); }
    break;

  case 9:

    { (yyval.n_int) = 0; }
    break;

  case 10:

    { prompt--; code(OP_ENDOFFILE); return(0); }
    break;

  case 11:

    { (yyval.n_int) = get_progoff (); }
    break;

  case 12:

    { code(OP_JMP_ENDOFFILE); }
    break;

  case 13:

    { (yyval.n_int) = (yyvsp[(1) - (2)].n_int); }
    break;

  case 14:

    {
                (yyval.n_int) = (yyvsp[(1) - (2)].n_int);
                printf("WARNING: File %s ended with statement and not an empty line or comment: Memory leak likely!\n",
                      curr_file_name);
                code(OP_ENDOFFILE);
              }
    break;

  case 17:

    { code(OP_JMP_ENDOFFILE); }
    break;

  case 18:

    {
  if (do_eval && !class_scope)
  {
    /* Only for eval() */
    code (OP_SAVE_EVAL);
  }
  else
  {
    /* Normal operation */
    switch ((yyvsp[(2) - (2)].n_int))
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
    break;

  case 35:

    { code(OP_QUIT); }
    break;

  case 36:

    { (yyval.n_int) = 2; /*asign to 'ans' and print*/}
    break;

  case 37:

    { (yyval.n_int) = 0; }
    break;

  case 38:

    { (yyval.n_int) = 1; /*print*/}
    break;

  case 39:

    { code(OP_RFILE); }
    break;

  case 40:

    { code(OP_RFILE_NAME); codep((yyvsp[(2) - (2)].name)); }
    break;

  case 41:

    { code(OP_REQ_NAME); codep((yyvsp[(2) - (2)].name)); }
    break;

  case 42:

    { (yyval.name) = (yyvsp[(1) - (1)].name); }
    break;

  case 43:

    { (yyval.name) = strappend ((yyvsp[(1) - (2)].name), (yyvsp[(2) - (2)].name)); GC_FREE ((yyvsp[(1) - (2)].name)); GC_FREE ((yyvsp[(2) - (2)].name)); }
    break;

  case 48:

    {
    (yyval.n_int) = code (OP_JMP);
    code (STOP);
  }
    break;

  case 49:

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
    break;

  case 50:

    { (yyval.n_int) = get_progoff (); }
    break;

  case 55:

    { (yyval.n_int) = (yyvsp[(2) - (3)].n_int); }
    break;

  case 56:

    {
  (yyval.n_int) = (yyvsp[(1) - (10)].n_int);
  code_sp ((yyvsp[(2) - (10)].n_int) + 1, (yyvsp[(9) - (10)].n_int) - ((yyvsp[(2) - (10)].n_int) + 1));  // if 'false' jump to next if statement
  code_sp ((yyvsp[(7) - (10)].n_int) + 1, (yyvsp[(10) - (10)].n_int) - ((yyvsp[(7) - (10)].n_int) + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
    break;

  case 57:

    {
  (yyval.n_int) = (yyvsp[(1) - (10)].n_int);
  code_sp ((yyvsp[(2) - (10)].n_int) + 1, (yyvsp[(9) - (10)].n_int) - ((yyvsp[(2) - (10)].n_int) + 1));  // if 'false' jump to next if statement
  code_sp ((yyvsp[(7) - (10)].n_int) + 1, (yyvsp[(10) - (10)].n_int) - ((yyvsp[(7) - (10)].n_int) + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
    break;

  case 58:

    {
  (yyval.n_int) = (yyvsp[(1) - (10)].n_int);
  code_sp ((yyvsp[(2) - (10)].n_int) + 1, (yyvsp[(9) - (10)].n_int) - ((yyvsp[(2) - (10)].n_int) + 1));  // if 'false' jump to after 'else' statement
  code_sp ((yyvsp[(7) - (10)].n_int) + 1, (yyvsp[(10) - (10)].n_int) - ((yyvsp[(7) - (10)].n_int) + 1)); // 'else' jumps to the end of compound statement
  prompt--;
}
    break;

  case 59:

    {
  /* if 'cond' $2 is not true, jump to 'end' $7 */
  (yyval.n_int) = (yyvsp[(1) - (7)].n_int);
  code_sp ((yyvsp[(2) - (7)].n_int) + 1, (yyvsp[(7) - (7)].n_int) - ((yyvsp[(2) - (7)].n_int) + 1));
  prompt--;
}
    break;

  case 60:

    {
    /* remove this from stack */
    code(OP_PULL_DATUM);
    sw_switch_datum = op_pull_datum; 
  }
    break;

  case 61:

    {
    (yyval.n_int) = (yyvsp[(9) - (12)].n_int);
    switching--;
    prompt--;
    int sw_end_jump = (sw_progoff_default > 0 ? sw_progoff_default : (yyvsp[(12) - (12)].n_int));
    resolve_case_break_tags (&sw_blist, &sw_clist, &sw_elist, (yyval.n_int), (yyvsp[(12) - (12)].n_int), sw_end_jump);
  }
    break;

  case 62:

    {
    sw_progoff_default=-1;
    prompt++;
    switching++;
    (yyval.n_int) = get_progoff ();
  }
    break;

  case 63:

    { (yyval.n_int)=0; }
    break;

  case 64:

    { (yyval.n_int)=(yyvsp[(1) - (1)].n_int); }
    break;

  case 65:

    { (yyval.n_int)=(yyvsp[(1) - (1)].n_int); }
    break;

  case 66:

    { (yyval.n_int)=(yyvsp[(1) - (2)].n_int); }
    break;

  case 67:

    { (yyval.n_int)=(yyvsp[(1) - (1)].n_int); }
    break;

  case 68:

    { (yyval.n_int)=(yyvsp[(1) - (2)].n_int); }
    break;

  case 69:

    {
    tag_brk_cont_case (&sw_clist, get_progoff()); /* location of prev OP_IFSJMP */
  }
    break;

  case 70:

    { /* value is already on the stack */
    code(OP_PUSH_DATUM);
    code(OP_EQ);                            /* are they equal */
    tag_brk_cont_case (&sw_elist, get_progoff());       /* location of next case statement OP_IFSJMP */
    code(OP_IFSJMP);                        /* where to jump if not true */
    code(STOP);                             /* placeholder for next case command */
  }
    break;

  case 72:

    { /* do a fake break switch that points to next case statement */
    (yyval.n_int) = get_progoff();
    code(OP_JMP);
    code(STOP);                           /* Placeholder */
    int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
    tag_brk_cont_case (&sw_blist, -n*switch_clist_offset-(yyval.n_int));   /* for fake brake switch */
  }
    break;

  case 73:

    { /* do a fake break switch that points to next case statement */
    (yyval.n_int) = get_progoff();
    code(OP_JMP);
    code(STOP);                           /* Placeholder */
    int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
    tag_brk_cont_case (&sw_blist, -n*switch_clist_offset-(yyval.n_int));   /* for fake brake switch */
  }
    break;

  case 74:

    {
    sw_progoff_default = get_progoff();
  }
    break;

  case 76:

    {
  (yyval.n_int) = get_progoff ();
  looping ++;
}
    break;

  case 77:

    {
    (yyval.n_int) = (yyvsp[(1) - (8)].n_int);
    code_sp ((yyvsp[(6) - (8)].n_int) + 1, (yyvsp[(8) - (8)].n_int) - (yyvsp[(6) - (8)].n_int) + 1 );
    code (OP_JMP);
    code ((yyvsp[(1) - (8)].n_int) - (yyvsp[(8) - (8)].n_int) - 1);
    resolve_break_tags ((yyval.n_int), (yyvsp[(8) - (8)].n_int) + 2, 0);
    resolve_continue_tags ((yyval.n_int), (yyvsp[(8) - (8)].n_int) + 2, 0);
    prompt--;
    looping-=2; /* key word 'while' does looping++ too */
  }
    break;

  case 78:

    {
    (yyval.n_int) = (yyvsp[(1) - (8)].n_int);
    code_sp ((yyvsp[(6) - (8)].n_int) + 1, (yyvsp[(1) - (8)].n_int) - (yyvsp[(8) - (8)].n_int) - 1);
    code (OP_JMP);
    code ((yyvsp[(2) - (8)].n_int) + 1);
    resolve_break_tags ((yyval.n_int), (yyvsp[(8) - (8)].n_int) + 2, 0);
    resolve_continue_tags ((yyval.n_int), (yyvsp[(8) - (8)].n_int) + 2, 0);
    prompt--;
    looping--;
  }
    break;

  case 79:

    {
    (yyval.n_int) = (yyvsp[(1) - (5)].n_int);
    if (for_then > 0)
    {
      code_sp ((yyvsp[(2) - (5)].n_int) + 1, for_then - (yyvsp[(2) - (5)].n_int));
      code_sp (for_then, (yyvsp[(1) - (5)].n_int) - for_then);
      resolve_continue_tags ((yyval.n_int), for_then - 1, 0);
      resolve_break_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int), 0);
    }
    else
    {
      code_sp ((yyvsp[(2) - (5)].n_int) + 1, (yyvsp[(5) - (5)].n_int) - (yyvsp[(2) - (5)].n_int) + 1);
      code (OP_JMP);
      code ((yyvsp[(1) - (5)].n_int) - (yyvsp[(5) - (5)].n_int) - 1);
      resolve_continue_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int) + 2, 0);
      resolve_break_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int) + 2, 0);
    }
    looping--;
    prompt--;
    for_then=0;
  }
    break;

  case 80:

    {
    (yyval.n_int) = (yyvsp[(1) - (5)].n_int);
    if (for_then > 0)
    {
      code_sp(for_then, (yyvsp[(5) - (5)].n_int) - for_then);
    }
    code_sp ((yyvsp[(2) - (5)].n_int) + 1, (yyvsp[(5) - (5)].n_int) - ((yyvsp[(2) - (5)].n_int)+1));
    code (OP_FOR_LOOP);
    code (OP_FOR_LOOP_DONE);
    resolve_break_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int), 1);
    if (for_then > 0)
    {
      resolve_continue_tags ((yyval.n_int), for_then-1, 1);
    }
    else
    {
      resolve_continue_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int), 1);
    }
    resolve_return_tags ((yyval.n_int), (yyvsp[(5) - (5)].n_int));
    looping--;
    prompt--;
    for_then=0;
  }
    break;

  case 81:

    {
    (yyval.n_int) = get_progoff ();
    looping++;
    prompt++;
    for_then=-1;
  }
    break;

  case 82:

    {
    (yyval.n_int) = code (OP_FOR_LOOP_I);
    code (STOP);
  }
    break;

  case 84:

    {
    int tmp;
    (yyval.n_int) = (yyvsp[(2) - (3)].n_int);
    if(scope == GLOBAL)
    {
      yyerror("return not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    code (OP_FUNCTION_RETURN);
    tmp = code (STOP);
    tag_brk_cont_case(&rlist, tmp);
  }
    break;

  case 85:

    {
    (yyval.n_int) = code(OP_IFSJMP);
    code (STOP);               /* place holder for jump val */
  }
    break;

  case 86:

    {
    for_then=-2;
    (yyval.n_int) = get_progoff ();
    looping ++;
    prompt++;
  }
    break;

  case 87:

    {
    (yyval.n_int) = get_progoff ();
    looping ++;
    prompt++;
  }
    break;

  case 88:

    {
    (yyval.n_int) = get_progoff ();
    prompt++;
  }
    break;

  case 89:

    {
    if (switching)
    {
      int n = list_GetNumNodes(sw_clist);   /* node with which break statements is associated */
      (yyval.n_int) = code (OP_JMP);
      code(STOP);               /* Placeholder */
      tag_brk_cont_case (&sw_blist, n*switch_clist_offset + (yyval.n_int));
    }
    else if(looping)
    {
      if (for_then > 0)
        if (get_progoff() > for_then )
        {
          yyerror("break-stmt cannot be after then-stmt in the same loop");
          longjmp( *jmp_dec_buff (), 1 );
        }
      (yyval.n_int) = code (OP_JMP);
      code(STOP);               /* Placeholder */
      tag_brk_cont_case (&blist, (yyval.n_int));
    }
    else
    {
      yyerror("break-stmt not allowed outside loop or switch environments");
      longjmp( *jmp_dec_buff (), 1 );
    }
  }
    break;

  case 90:

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
    (yyval.n_int) = code(OP_JMP);
    code(STOP);            /* Placeholder */
    tag_brk_cont_case (&clist, (yyval.n_int));
  }
    break;

  case 91:

    { (yyval.n_int) = get_progoff (); }
    break;

  case 95:

    { (yyval.n_int) = (yyvsp[(2) - (3)].n_int); }
    break;

  case 96:

    { code(OP_VECTOR_CREATE); code(2); }
    break;

  case 97:

    { code(OP_VECTOR_CREATE); code(3); }
    break;

  case 98:

    { code(OP_ASSIGN); }
    break;

  case 99:

    { code(OP_MATRIX_ASSIGN); code(1); }
    break;

  case 100:

    { code(OP_MATRIX_ASSIGN); code(2); }
    break;

  case 101:

    { code(OP_MATRIX_ASSIGN); code(3); }
    break;

  case 102:

    { code(OP_MATRIX_VEC_ASSIGN); }
    break;

  case 103:

    { code(OP_LIST_ASSIGN); code(1); }
    break;

  case 104:

    { code(OP_LIST_ASSIGN); code(2); codep((yyvsp[(3) - (5)].name)); }
    break;

  case 105:

    { code(OP_OLIST_ASSIGN); }
    break;

  case 106:

    { (yyval.n_int) = (yyvsp[(2) - (3)].n_int); }
    break;

  case 108:

    { (yyval.n_int) = code(OP_EMPTY_MATRIX_CREATE); }
    break;

  case 109:

    { code(OP_MATRIX_CREATE); }
    break;

  case 110:

    { code(OP_MATRIX_APPEND); }
    break;

  case 111:

    { code(OP_MATRIX_APPEND); }
    break;

  case 114:

    { code(OP_VEC_APPEND); }
    break;

  case 115:

    {  code(OP_MATRIX_SUB); code(1); }
    break;

  case 116:

    {  code(OP_MATRIX_SUB); code(2); }
    break;

  case 117:

    {  code(OP_MATRIX_SUB); code(3); }
    break;

  case 118:

    { /* no-op */ }
    break;

  case 119:

    { /* no-op */ }
    break;

  case 120:

    {  code(OP_MATRIX_VEC_SUB); }
    break;

  case 121:

    {  code(OP_MATRIX_COL); }
    break;

  case 126:

    { (yyval.n_int) = (yyvsp[(3) - (3)].n_int); code(OP_LIST_EL_CREATE); codep((yyvsp[(1) - (3)].name)); }
    break;

  case 128:

    { /* default action causes type conflict */ }
    break;

  case 129:

    {}
    break;

  case 135:

    {
    (yyval.n_int) = (yyvsp[(2) - (3)].nn).off;
    code(OP_LIST_CREATE);
    code((yyvsp[(2) - (3)].nn).count);
  }
    break;

  case 136:

    { code(OP_LIST_MEMB); code(1); }
    break;

  case 137:

    { code(OP_LIST_MEMB); code(2); codep((yyvsp[(3) - (3)].name)); }
    break;

  case 138:

    { (yyval.nn).count = 0; (yyval.nn).off = get_progoff (); }
    break;

  case 139:

    { (yyval.nn).count = 1; (yyval.nn).off = (yyvsp[(1) - (1)].n_int); }
    break;

  case 140:

    { (yyval.nn).count = (yyvsp[(1) - (3)].nn).count + 1; (yyval.nn).off = (yyvsp[(1) - (3)].nn).off; }
    break;

  case 141:

    {
    (yyval.n_int) = (yyvsp[(2) - (3)].nn).off;
    code(OP_OLIST);
    code((yyvsp[(2) - (3)].nn).count);
  }
    break;

  case 142:

    {
    (yyval.nn).count = 1;
    (yyval.nn).off = (yyvsp[(1) - (1)].n_int);
  }
    break;

  case 143:

    {
    (yyval.nn).count = (yyvsp[(1) - (3)].nn).count + 1;
    (yyval.nn).off = (yyvsp[(1) - (3)].nn).off;
  }
    break;

  case 144:

    { (yyval.n_int)=(yyvsp[(1) - (1)].n_int); }
    break;

  case 146:

    { code(OP_INC); }
    break;

  case 147:

    { code(OP_DEC); }
    break;

  case 148:

    {
    Var *ret = name_lookup_func (arg_sym_tab, lsymtab, gsymtab, ssymtab, (yyvsp[(1) - (4)].name), scope);
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
    code((yyvsp[(3) - (4)].nn).count);   /* number of args on stack */
    (yyval.n_int) = (yyvsp[(3) - (4)].nn).off;
  }
    break;

  case 149:

    {
    /* MDEs are here to stay! MK IX-2016 */
    code(OP_FUNCTION_CALL_2);
    code((yyvsp[(3) - (4)].nn).count);   /* number of args on stack */
    (yyval.n_int) = (yyvsp[(3) - (4)].nn).off;
  }
    break;

  case 150:

    {
    code(OP_FUNCTION_CALL_SELF);
    code((yyvsp[(3) - (4)].nn).count);   /* number of args on stack */
    (yyval.n_int) = (yyvsp[(3) - (4)].nn).off;
  }
    break;

  case 151:

    {
    code(OP_FUNCTION_CALL_1);
    code((yyvsp[(3) - (4)].nn).count);   /* number of args on stack */
  }
    break;

  case 152:

    { code(OP_ADD); }
    break;

  case 153:

    { code(OP_SUB); }
    break;

  case 154:

    { code(OP_MUL); }
    break;

  case 155:

    { code(OP_EL_MUL); }
    break;

  case 156:

    { code(OP_DIV); }
    break;

  case 157:

    { code(OP_LDIV); }
    break;

  case 158:

    { code(OP_EL_DIV); }
    break;

  case 159:

    { code(OP_EL_LDIV); }
    break;

  case 160:

    { code(OP_POWER); }
    break;

  case 161:

    { code(OP_EL_POWER); }
    break;

  case 162:

    { code(OP_TRANSPOSE); }
    break;

  case 163:

    { code(OP_EL_TRANSPOSE); }
    break;

  case 164:

    { code(OP_GT);  }
    break;

  case 165:

    { code(OP_GE);  }
    break;

  case 166:

    { code(OP_LT);  }
    break;

  case 167:

    { code(OP_LE);  }
    break;

  case 168:

    { code(OP_EQ);  }
    break;

  case 169:

    { code(OP_NE);  }
    break;

  case 170:

    { code(OP_AND); }
    break;

  case 171:

    { code(OP_OR);  }
    break;

  case 172:

    { code(OP_ADDTO);     }
    break;

  case 173:

    { code(OP_SUBFROM);   }
    break;

  case 174:

    { code(OP_EL_MUL_BY); }
    break;

  case 175:

    { code(OP_EL_DIV_BY); }
    break;

  case 176:

    { (yyval.n_int) = (yyvsp[(2) - (2)].n_int); code(OP_NOT); }
    break;

  case 177:

    { (yyval.n_int) = (yyvsp[(2) - (2)].n_int);   code(OP_NEGATE); }
    break;

  case 178:

    { (yyval.n_int) = (yyvsp[(2) - (2)].n_int); }
    break;

  case 179:

    { (yyval.n_int) = code(OP_PUSH_CONSTANT); coded((yyvsp[(1) - (1)].d_val)); }
    break;

  case 180:

    { (yyval.n_int) = code(OP_PUSH_INTEGER); coded((yyvsp[(1) - (1)].d_val)); }
    break;

  case 181:

    { (yyval.n_int) = code(OP_PUSH_iCONSTANT); coded((yyvsp[(1) - (1)].d_val)); }
    break;

  case 182:

    { (yyval.n_int) = code(OP_PUSH_STRING); codep((yyvsp[(1) - (1)].name)); }
    break;

  case 183:

    { (yyval.n_int) = code(OP_PUSH_STRING); codep((yyvsp[(1) - (1)].name)); }
    break;

  case 184:

    {
    if((yyvsp[(1) - (1)].vinfo)->type == GLOBAL)
    {
      (yyval.n_int) = code(OP_PUSH_VAR);
      codep((yyvsp[(1) - (1)].vinfo)->var);
    }
    else if((yyvsp[(1) - (1)].vinfo)->type == LOCAL_VAR)
    {
      (yyval.n_int) = code(OP_PUSH_LOCAL_VAR);
      code((yyvsp[(1) - (1)].vinfo)->offset);
    }
    else if((yyvsp[(1) - (1)].vinfo)->type == STATIC_VAR)
    {
      (yyval.n_int) = code(OP_PUSH_STATIC_VAR);
      code((yyvsp[(1) - (1)].vinfo)->offset);
    }
    else if((yyvsp[(1) - (1)].vinfo)->type == ARG_VAR)
    {
      (yyval.n_int) = code(OP_PUSH_ARG);
      code((yyvsp[(1) - (1)].vinfo)->offset);
    }
    GC_FREE((yyvsp[(1) - (1)].vinfo));
  }
    break;

  case 185:

    {
    (yyval.vinfo) = name_lookup (arg_sym_tab, lsymtab, gsymtab, ssymtab, (yyvsp[(1) - (1)].name), scope);
  }
    break;

  case 186:

    { (yyval.vinfo) = gst (); }
    break;

  case 187:

    { (yyval.n_int) = 0;}
    break;

  case 188:

    { (yyval.n_int) = 0; }
    break;

  case 189:

    { (yyval.n_int) = 0; }
    break;

  case 190:

    { (yyval.n_int) = 0; }
    break;

  case 191:

    { (yyval.n_int) = 0; }
    break;

  case 192:

    { (yyval.n_int) = 0; }
    break;

  case 193:

    { (yyval.n_int) = 0; }
    break;

  case 194:

    { /* empty */ }
    break;

  case 195:

    {
    (yyval.n_int) = code(OP_PUSH_VAR); codep((yyvsp[(1) - (1)].ent));
  }
    break;

  case 196:

    {
    prompt++;
  }
    break;

  case 197:

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
    arg_sym_tab = (yyvsp[(4) - (8)].list);
    lsymtab = list_Create ();
    gsymtab = list_Create ();
    ssymtab = 0;
    function_setup1 (lsave, curr_file_name);
  }
    break;

  case 198:

    {
    code (OP_DEF_FUNC_RET);
    code (STOP);
    (yyval.ent) = function_setup2 (arg_sym_tab, lsymtab, gsymtab, ssymtab, get_progoff ());
    arg_sym_tab = 0;
    lsymtab = 0;
    gsymtab = 0;
    ssymtab = 0;
    if (rlist)
      list_DestroyAllNodes (rlist);
    scope = GLOBAL;
    prompt--;
  }
    break;

  case 199:

    {
  if(scope == GLOBAL)
  {
    yyerror("$self not allowed outside function");
    longjmp( *jmp_dec_buff (), 1 );
  }
}
    break;

  case 200:

    { (yyval.list) = list_Create(); }
    break;

  case 201:

    { (yyval.list) = arg_var_push(0, (yyvsp[(1) - (1)].name)); }
    break;

  case 202:

    { arg_var_push((yyvsp[(1) - (3)].list), (yyvsp[(3) - (3)].name)); }
    break;

  case 203:

    {
    if (scope == GLOBAL)
    {
      yyerror("local decl not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    dlocal = 1;
  }
    break;

  case 204:

    {
    (yyval.n_int) = get_progoff ();
    dlocal = 0;
  }
    break;

  case 205:

    {
    if (scope == GLOBAL)
    {
      yyerror("global decl not allowed outside function");
      longjmp( *jmp_dec_buff (), 1 );
    }
    dglobal = 1;
  }
    break;

  case 206:

    {
    (yyval.n_int) = get_progoff ();
    dglobal = 0;
  }
    break;

  case 207:

    {
    (yyval.n_int) = get_progoff ();
    if (dlocal)
      local_var_push(lsymtab, (yyvsp[(1) - (1)].name));
    else
      global_var_push(gsymtab, (yyvsp[(1) - (1)].name));
  }
    break;

  case 208:

    {
    if (dlocal)
      local_var_push(lsymtab, (yyvsp[(3) - (3)].name));
    else
      global_var_push(gsymtab, (yyvsp[(3) - (3)].name));
  }
    break;

  case 209:

    { (yyval.nn).count = 0; (yyval.nn).off = get_progoff (); }
    break;

  case 210:

    { (yyval.nn).count = 1; (yyval.nn).off = (yyvsp[(1) - (1)].n_int); }
    break;

  case 211:

    {
    (yyval.nn).count = (yyvsp[(1) - (1)].nn).count;
    (yyval.nn).off = (yyvsp[(1) - (1)].nn).off;
  }
    break;

  case 212:

    { (yyval.nn).count = 2; (yyval.nn).off = (yyvsp[(1) - (3)].n_int); }
    break;

  case 213:

    { (yyval.nn).count = (yyvsp[(1) - (3)].nn).count + 1; (yyval.nn).off = (yyvsp[(1) - (3)].nn).off; }
    break;

  case 214:

    { code(OP_PUSH_UNDEF); }
    break;

  case 216:

    { (yyval.n_int) = get_progoff (); }
    break;

  case 217:

    {
    if(scope == GLOBAL)
    {
      if (class_scope>0)
      {
        btree_class_stat_var_push (class_stat_symtab, (yyvsp[(1) - (1)].name));
      }
      else
      {
        fstatic_var_push(curr_file_name, (yyvsp[(1) - (1)].name));
      }
    }
    else if (scope == LOCAL)
    {
      var_push(&ssymtab, (yyvsp[(1) - (1)].name));
    }
    (yyval.n_int) = 0;
  }
    break;

  case 218:

    {
    if(scope == GLOBAL)
    {
      if (class_scope>0)
      {
        btree_class_stat_var_push (class_stat_symtab, (yyvsp[(3) - (3)].name));
      }
      else
      {
        fstatic_var_push(curr_file_name, (yyvsp[(3) - (3)].name));
      }
    }
    else if (scope == LOCAL)
    {
      var_push(&ssymtab, (yyvsp[(3) - (3)].name));
    }
    (yyval.n_int) = 0;
  }
    break;

  case 219:

    { read_script=1; }
    break;

  case 220:

    { read_script=0; (yyval.name) = (yyvsp[(2) - (2)].name); }
    break;

  case 221:

    { (yyval.n_int) = get_progoff (); }
    break;

  case 222:

    {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      fstatic_var_push(priv_class_name, (yyvsp[(1) - (1)].name));
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    (yyval.n_int) = 0;
  }
    break;

  case 223:

    {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      fstatic_var_push(priv_class_name, (yyvsp[(3) - (3)].name));
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    (yyval.n_int) = 0;
  }
    break;

  case 224:

    {
    (yyval.n_int) = get_progoff ();
  }
    break;

  case 225:

    {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      btree_class_publ_var_push (class_publ_symtab, (yyvsp[(1) - (1)].name));
    }
    else
    {
      yyerror("public decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    (yyval.n_int) = 0;
  }
    break;

  case 226:

    {
    if ((scope == GLOBAL) && (class_scope>0))
    {
      btree_class_publ_var_push (class_publ_symtab, (yyvsp[(3) - (3)].name));
    }
    else
    {
      yyerror("private decl allowed only inside classdef");
      longjmp( *jmp_dec_buff (), 1 );
    }
    (yyval.n_int) = 0;
  }
    break;

  case 227:

    { (yyval.list) = list_Create(); }
    break;

  case 228:

    { (yyval.list) = arg_var_push(0, (yyvsp[(1) - (1)].name)); }
    break;

  case 229:

    { arg_var_push((yyvsp[(1) - (3)].list), (yyvsp[(3) - (3)].name)); }
    break;

  case 230:

    {
    (yyval.n_int) = code(OP_PUSH_VAR); codep((yyvsp[(1) - (1)].ent));
  }
    break;

  case 231:

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
    break;

  case 232:

    {
    read_script=1;
    scope = GLOBAL;
    class_args_symtab = (yyvsp[(4) - (5)].list);
    psave = get_progoff ();
    lsave = lineno+loff;
  }
    break;

  case 233:

    {
    read_script=0;
    (yyval.ent) = classdef_setup ((yyvsp[(4) - (8)].list), curr_file_name, (yyvsp[(8) - (8)].name), lsave);
    class_args_symtab=0;
    if (rlist)
      list_DestroyAllNodes (rlist);
    scope = GLOBAL;
    prompt--;
  }
    break;



      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}





#include "rlab_parser_f.c"
