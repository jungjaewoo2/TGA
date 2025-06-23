#define INTEGER 257
#define NUMBER 258
#define iNUMBER 259
#define REOF 260
#define ENDFILE 261
#define NAME 262
#define R_STRING 263
#define HELP_NAME 264
#define FNAME 265
#define UNDEFINED 266
#define SYS_CMD 267
#define R_SCRIPT 268
#define CLASS_PUBL 269
#define CLASS_PRIV 270
#define DO 271
#define UNTIL 272
#define WHILE 273
#define IF 274
#define ELSE 275
#define QUIT 276
#define FOR 277
#define IN 278
#define BREAK 279
#define CONTINUE 280
#define LOOPTHEN 281
#define SWITCH 282
#define SWITCH_CASE 283
#define SWITCH_DEFAULT 284
#define FUNCTION 285
#define RETURN 286
#define SEMIC 287
#define Q_MARK 288
#define LOCAL_DEC 289
#define GLOBAL_DEC 290
#define SELF 291
#define HELP 292
#define RFILE 293
#define REQUIRE 294
#define GST 295
#define FSTATIC_DEC 296
#define RSCRIPT_DEC 297
#define CLASSDEF 298
#define JNK 299
#define OR 300
#define AND 301
#define EQ 302
#define NE 303
#define GT 304
#define GE 305
#define LT 306
#define LE 307
#define ADDTO 308
#define SUBFROM 309
#define EL_MUL_BY 310
#define EL_DIV_BY 311
#define EL_MUL_OP 312
#define EL_DIV_OP 313
#define EL_LDIV_OP 314
#define LDIV_OP 315
#define UNARY_MINUS 316
#define UNARY_PLUS 317
#define EL_POW_OP 318
#define TRANSPOSE 319
#define EL_TRANSPOSE 320
#define NOT 321
#define INC 322
#define DEC 323
#define LEFT_LIST 324
#define RIGHT_LIST 325
#define L_OLIST 326
#define R_OLIST 327
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union {
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
#endif /* !YYSTYPE_IS_DECLARED */
extern YYSTYPE yylval;
