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
#define INTEGER 257
#define	NUMBER	258
#define	iNUMBER	259
#define	NAME	260
#define	R_STRING	261
#define	HELP_NAME	262
#define	FNAME	263
#define	UNDEFINED	264
#define	SYS_CMD	265
#define	WHILE	266
#define	IF	267
#define	ELSE	268
#define	QUIT	269
#define	FOR	270
#define	IN	271
#define	BREAK	272
#define	CONTINUE	273
#define	FUNCTION	274
#define	RETURN	275
#define	SEMIC	276
#define	Q_MARK	277
#define	LOCAL_DEC	278
#define	GLOBAL_DEC	279
#define	SELF	280
#define	HELP	281
#define	RFILE	282
#define	REQUIRE	283
#define	GST	284
#define	FSTATIC_DEC	285
#define	JNK	286
#define	OR	287
#define	AND	288
#define	EQ	289
#define	NE	290
#define	GT	291
#define	GE	292
#define	LT	293
#define	LE	294
#define	EL_ADD_OP	295
#define	EL_SUB_OP	296
#define	EL_MUL_OP	297
#define	EL_DIV_OP	298
#define	EL_LDIV_OP	299
#define	LDIV_OP	300
#define	UNARY_MINUS	301
#define	UNARY_PLUS	302
#define	EL_POW_OP	303
#define	TRANSPOSE	304
#define	EL_TRANSPOSE	305
#define	NOT	306
#define	INC	307
#define	DEC	308
#define	LEFT_LIST	309
#define	RIGHT_LIST	310
#define	L_OLIST	311
#define	R_OLIST	312
#define CLASS 313

extern YYSTYPE yylval;
