/* **************************************************************
 * Lookup a NAME in the symbol table, and return a ptr to the 
 * entity. This function pays attention to the scope of the 
 * variable usage, and looks for the NAME in the appropriate
 * place.
 * ************************************************************** */

char * get_current_file_name(void)
{
  return curr_file_name;
}

char * set_priv_class_name(char * class_ptr_name)
{
  sprintf(priv_class_name, "%s", class_ptr_name);
  return priv_class_name;
}

char * set_classdef_filename(char * ptr)
{
  classdef_filename = ptr;
  return classdef_filename;
}

//
//
//
Btree * get_class_stat_symtab()
{
  if (! class_stat_symtab)
    class_stat_symtab = btree_Create();
  return class_stat_symtab;
}

Btree * just_get_class_stat_symtab()
{
  return class_stat_symtab;
}

void null_class_stat_symtab()
{
  class_stat_symtab = 0;
  return ;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "btree_class_stat_var_push"
void btree_class_stat_var_push (Btree *class_stat_symtab, char *name)
{
  if (! class_stat_symtab)
    class_stat_symtab = (Btree *) get_class_stat_symtab();

  ListNode *lnode= btree_FindNode (class_stat_symtab, name);
  if (!lnode)
  {
    Ent *ent = ent_Create();
    ent_type(ent) = UNDEF;
    install (class_stat_symtab, name, ent);
  }
}

//
//
//
Btree * get_class_publ_symtab()
{
  if (! class_publ_symtab)
    class_publ_symtab = btree_Create();
  return class_publ_symtab;
}
void null_class_publ_symtab()
{
  class_publ_symtab = 0;
  return;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "btree_class_publ_var_push"
void btree_class_publ_var_push (Btree *class_publ_symtab, char *name)
{
  if (! class_publ_symtab)
    class_publ_symtab = (Btree *) get_class_publ_symtab();

  ListNode *lnode= btree_FindNode (class_publ_symtab, name);
  if (!lnode)
  {
    Ent *ent = ent_Create();
    ent_type(ent) = UNDEF;
    install (class_publ_symtab, name, ent);
  }
}

int inc_class_scope()
{
  ++class_scope;
  return class_scope;
}

int dec_class_scope()
{
  if (class_scope>0)
    --class_scope;
  return class_scope;
}

int get_class_scope()
{
  return class_scope;
}



extern int   get_function_scope (void);       /* code.c */
extern List *get_function_arglist (void);     /* code.c */
extern List *get_function_locallist (void);   /* code.c */
extern List *get_function_globallist (void) ; /* code.c */
extern List *get_function_statlist (void) ; /* code.c */

#undef  THIS_SOLVER
#define THIS_SOLVER "name_lookup"
static Var * name_lookup(List *arg_list, List *local_list, List *global_list,
                         List *stat_list, char *name, int scope)
{
  ListNode *var=0;

  Var *retval = (Var *) GC_MALLOC(sizeof(Var));
  if (retval == 0)
    rerror (THIS_FILE ": " THIS_SOLVER ": " RLAB_ERROR_OUT_OF_MEMORY);

  if (do_eval && !scope)
  {
    /* Adjust things so eval string will get proper scope, etc... */
    scope = get_function_scope ();
    if (scope == LOCAL)
    {
      arg_list = get_function_arglist ();
      local_list = get_function_locallist ();
      global_list = get_function_globallist ();
      stat_list = get_function_statlist ();
    }
  }

  if (scope == LOCAL)
  {
    if((var = lookup(arg_list, name)) != 0)
    {
      /* ARG */
      retval->type = ARG_VAR;
      retval->offset = lvar_GetOffset(var_Data(var));
      retval->name = name;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(class_stat_symtab, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(class_publ_symtab, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope>0) && (var = static_lookup(classdef_filename, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope==0) && (var = static_lookup(curr_file_name, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ( (class_scope>0) && ((var = static_lookup(priv_class_name, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((var = lookup_private(stat_list, name)) != 0)
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": found %s in 'stat_list'\n", name);
      GC_FREE(name);
      retval->type = STATIC_VAR;
      retval->offset = (int) (var->scope);
    }
    else if ((var = lookup_private(local_list, name)) != 0)
    {
      GC_FREE(name);
      retval->type = LOCAL_VAR;
      retval->offset = (int) lvar_GetOffset(var_Data(var));
    }
    else if ((var = lookup_private(global_list, name)) != 0)
    {
      if ((var = lookup (0, name)) != 0)
      {
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
      else
      {
        var = install(0, name, UNDEF);
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
    }
    else if (class_scope>0)
    {
      /* inside classdef all undefined are private */
      fstatic_var_push(priv_class_name, name);
      var = static_lookup(priv_class_name, name);
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else
    {
      /* It is not in the global-scope, create a local. */
      retval->type = LOCAL_VAR;
      local_var_push (local_list, name);
      retval->offset = list_GetNumNodes (local_list);
    }
  }
  else  /* scope == GLOBAL */
  {
    if ( (class_scope>0) && ((var = static_lookup(priv_class_name, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(class_stat_symtab, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(class_publ_symtab, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope>0) && (var = static_lookup(classdef_filename, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope==0) && (var = static_lookup(curr_file_name, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if((var = lookup(0, name)) != 0)
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if (class_scope>0)
    {
      /* inside classdef all undefined are private */
      fstatic_var_push(priv_class_name, name);
      var = static_lookup(priv_class_name, name);
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else
    {
      /* we could not find it, so create it */
      var = install(0, name, UNDEF);
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
  }

  return(retval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "name_lookup_func"
static Var *
name_lookup_func(List *arg_list, List *local_list, List *global_list, List* stat_list, 
                 char *name, int scope)
{
  ListNode *var=0;
  Var *retval = (Var *) GC_MALLOC(sizeof(Var));
  if (retval == 0)
    rerror (THIS_FILE ": " THIS_SOLVER ": " RLAB_ERROR_OUT_OF_MEMORY);

  if (do_eval && !scope)
  {
    /* Adjust things so eval string will get proper scope, etc... */
    scope = get_function_scope ();
    if (scope == LOCAL)
    {
      arg_list = get_function_arglist ();
      local_list = get_function_locallist ();
      global_list = get_function_globallist ();
      stat_list = get_function_statlist ();
    }
  }

  if (scope == LOCAL)
  {
    /* 2nd check function's arg list. */
    if((var = lookup(arg_list, name)) != 0)
    {
      /* ARG */
      retval->type = ARG_VAR;
      retval->offset = lvar_GetOffset(var_Data(var));
      retval->name = name;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(get_class_stat_symtab(), name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((var = lookup_private(stat_list, name)) != 0)
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": found %s in 'stat_list'\n", name);
      GC_FREE(name);
      retval->type = STATIC_VAR;
      retval->offset = (int) (var->scope);
    }
    else if ((var = lookup_private(local_list, name)) != 0)
    { /* 1st check function's local list. */
      GC_FREE(name);
      retval->type = LOCAL_VAR;
      retval->offset = (int) lvar_GetOffset(var_Data(var));
    }
    /* 3rd check function's global list. */
    else if(global_list && ((var = lookup(global_list, name)) != 0))
    { /* GLOBAL_VAR */
      if ((var = lookup (0, name)) != 0)
      { /* Really get it from the global symbol table. */
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
      else
      { /* Create it on the global symbol table. */
        var = install(0, name, UNDEF);
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
    }
    else if ((class_scope>0) && (var = static_lookup(classdef_filename, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope==0) && (var = static_lookup(curr_file_name, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else
    { /* we did not find it, look on global symbol table */
      if((var = lookup(0, name)) != 0)
      { /* Found it. */
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
      else
      { /* Create it on the global symbol table. */
        var = install(0, name, UNDEF);
        GC_FREE(name);
        retval->type = GLOBAL;
        retval->var = var;
      }
    }
  }
  else
  { /* scope == GLOBAL */
    if ( (class_scope>0) && ((var = static_lookup(priv_class_name, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(get_class_stat_symtab(), name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if( (class_scope>0) && ((var = btree_FindNode(class_publ_symtab, name)) != 0) )
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope>0) && (var = static_lookup(classdef_filename, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if ((class_scope==0) && (var = static_lookup(curr_file_name, name)))
    {
      //printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in 'curr_file_name'\n", name);
      //if (curr_file_name)
      //  printf(THIS_FILE ": " THIS_SOLVER ": IN: Lookup for %s in %s\n", name, curr_file_name);
      GC_FREE (name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else if((var = lookup(0, name)) != 0)
    {
      GC_FREE(name);
      retval->type = GLOBAL;
      retval->var = var;
    }
    else
    { 
      /* we could not find it, so create it */
      var = install(0, name, UNDEF);
      retval->type = GLOBAL;
      retval->var = var;
    }
  }

  return(retval);
}

/* **************************************************************
 * File Static Functions.
 * ************************************************************** */

List * get_static_tree (void)
{
  return (static_tree);
}

void init_static_tree (void)
{
  static_tree = list_Create ();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "static_tree_DestroyNodeByKey"
void static_tree_DestroyNodeByKey (char *key)
{
  if (!static_tree)
    return;
  if (!key)
    return;

  ListNode *bad_node=list_GetNodeByKey (static_tree, key);
  if (bad_node)
  {
    list_DestroyNodeByAddr(static_tree, bad_node);
  }
}

void destroy_fstatic_tree (void)
{
  list_Destroy (static_tree);
}

/*
 * Push a variable onto a file-static tree.
 * Don't return anything, this function is only called when
 * a fstatic() is encountered.
 */
#undef  THIS_SOLVER
#define THIS_SOLVER "fstatic_var_push"
void fstatic_var_push (char *file_name, char *name)
{
  ListNode *file_tree=list_GetNodeByKey (static_tree, file_name), *lnode=0;
  Ent *ent;

  if (file_tree)
  {
    /* Look for name */
    if (name)
    {
      if (!(lnode = btree_FindNode (var_Data (file_tree), name)))
      {
        /* Did not find the file-name, create it */
        install (var_Data (file_tree), name, UNDEF);
      }
    }
  }
  else
  {
    /* Create file_tree, and add name */
    ent = ent_Create ();
    ent_data (ent) = btree_Create ();
    ent_SetType (ent, BTREE);
    file_tree = list_Install (static_tree, file_name, ent);
    if (name)
      lnode = install (var_Data (file_tree), name, UNDEF);
  }
}

#undef  THIS_SOLVER
#define THIS_SOLVER "static_lookup"
ListNode * static_lookup (char *file_name, char *name)
{
  ListNode *file_tree, *lnode;

  if ((file_tree = list_GetNodeByKey (static_tree, file_name)))
  {
    if ((lnode = btree_FindNode (var_Data (file_tree), name)))
      return (lnode);
  }
  return (0);
}

/* **************************************************************
 * Called for YACC syntax error. Write parser error messages to
 * stderr (and diary file if necc.). Also force machine and scanner
 * to reset.
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "ryyerror"
void ryyerror(char *s, int yychar)
{
  char  *loc, *diag_str;
  int i, j;

  /* Print out error message */
  if( strcmp("stdin", curr_file_name) )
    fprintf(stderr, "%s on line: %d, file: %s\n", s, lineno, curr_file_name);
  else
    fprintf(stderr, "%s\n", s );

  /* Now try and print out offending tokens */
  if(strrchr(line_contents, '\n') == 0)
    fprintf(stderr, "%s\n", line_contents);
  else
    fprintf(stderr, "%s", line_contents);

  /* find next token/char on line */
  if((loc = strrchr(line_contents, yychar)) == 0)
    i = strlen(line_contents);
  else
    i = loc - line_contents;

  /* Create a string of blanks with a '^' where the 
     parse error is. Check value of loc, error may be at
     beginning of string (loc = 0) */
  if(i == 0)
  {
    diag_str = (char *) GC_MALLOC(2*sizeof(char));
    if (diag_str == 0)
      rerror ("out of memory");
    diag_str[0] = '^';
    diag_str[1] = '\0';
  }
  else
  {
    diag_str = (char *) GC_MALLOC((i+1)*sizeof(char));
    if (diag_str == 0)
      rerror ("out of memory");
    for(j = 0; j < i-1; j++)
      diag_str[j] = ' ';
    diag_str[j] = '^';
    diag_str[j+1] = '\0';
  }

  fprintf(stderr, "%s\n", diag_str);
  fflush(stderr);

  GC_FREE(diag_str);

  /* Reset the loop detection flag */
  looping = 0;

  /* Reset the switch detection flag */
  switching = 0;

  /*
   * Flush the input line since the rest is probably
   * Uninteligble.
   */
  flush_line = 1;

  /* Reset the scanner, in case we error'ed during a load() */
  /* new_file(0); */

  /* Most likely there are some bad op-codes in the 
     machine... reset */
  initcode();

  /* Reset the prompt, in case we error'ed in a loop */
  prompt = 0;

  /* 
   * Handle some clean-up if we were in the midst of 
   * parsing a function, then longjmp back to the prompt.
   */

  if (scope == LOCAL)
  {
    scope = GLOBAL;
    if(arg_sym_tab != 0) 
      list_Destroy(arg_sym_tab);
    if(lsymtab != 0)
      list_Destroy(lsymtab);
    if(gsymtab != 0)
      list_Destroy(gsymtab);
    if(ssymtab != 0)
      list_Destroy(ssymtab);

    function_setup_clear_nf();
  }

  /* Clean up the break and continue lists. */
  if (blist)
    list_DestroyAllNodes (blist);
  if (clist)
    list_DestroyAllNodes (clist);
  if (rlist)
    list_DestroyAllNodes (rlist);

  if (sw_blist)
    list_DestroyAllNodes (sw_blist);
  if (sw_clist)
    list_DestroyAllNodes (sw_clist);
  if (sw_elist)
    list_DestroyAllNodes (sw_elist);

  longjmp( *jmp_dec_buff (), 1 );
}

/* **************************************************************
 * Functions for handling break and continue statements.
 * ************************************************************** */

/*
 * Store the offset of each break statement so that
 * we can set the jump value later on when we have
 * finished the while or for statement.
 */
static void tag_brk_cont_case (List **list, int break_op)
{
  ListNode *lnode;

  if (!*list)
    *list = list_Create ();

  lnode = listNode_Create ();

  /* Use scope member as a shortcut. */
  lnode->scope = break_op;

  /* Add the listNode. */
  list_PushNode (*list, lnode);
}

/*
 * Now that we have finished a for or while loop
 * set the jump values for each break statement
 * within the loop or switch statement
 */
static void resolve_break_tags (int begin_op, int end_op, int tloop)
{
  if (!blist)
    return;

  ListNode *lnode=0, *pnode=0;

  lnode = list_GetLastNode (blist);

  /* 
   * See if there are any breaks in the current loop.
   */
  while (lnode)
  {
    pnode = listNode_GetPrevNode (lnode);
    if (lnode->scope > begin_op && lnode->scope < end_op)
    {
      code_sp (lnode->scope + 1, end_op - lnode->scope - 1 + tloop);     /* for */
      list_DetachNodeByAddr (blist, lnode);
      listNode_Destroy (lnode);
    }
    lnode = pnode;
  }

  // kilt the list at the end of the loop or switch statement
  //list_Destroy(*list);
  //*list=0;
}

#include "list.h"
static void resolve_case_break_tags (List **blist, List **clist, List **elist, int begin_op,
                                     int end_op, int dflt_end)
{
  if (!*blist)
    return;
  if (!*clist)
    return;
  if (!*elist)
    return;

  int n,j,n_prev;

  ListNode *l_bnode=0, *p_bnode=0;
  ListNode *l_cnode=0, *p_cnode=0;
  ListNode *l_enode=0, *p_enode=0;

  /*
   * do the breaks first: break node scope keeps track of location
   * from which the jump is performed and the cnode index to which
   * it belongs:
   *  this is because construct may cause multiple break; statement
   *  to be associated with the same cnode, e.g., break; as a part
   *  of the if statement
   *      scope = switch_clist_offset * (cnode index) + (true jump)
   */
  l_bnode = list_GetLastNode (*blist);
  n = l_bnode->scope / switch_clist_offset;     // cnode index
  if (n<0) n = -n;

  l_enode = list_GetLastNode (*elist);
  if (l_enode)
    p_enode = listNode_GetPrevNode (l_enode);
  else
    p_enode=0;

  while (l_bnode)
  {
    p_bnode = listNode_GetPrevNode (l_bnode);

    n_prev = n;
    n = l_bnode->scope / switch_clist_offset;     // cnode index
    j = l_bnode->scope - n * switch_clist_offset; // true jump (+real / -fake)
    if (n<0) n = -n;

    if (n != n_prev)
    {
      // jump to next enode only if n has changed
      if (l_enode)
        p_enode = listNode_GetPrevNode (l_enode);
      else
        p_enode=0;
    }

    if (j > 0)
    {
      /* true break statements: go to the end of switch environment */
      code_sp (j + 1, end_op - j - 1);
    }
    else
    {
      /* fake break statements: go to the next elist entry */
      if (p_enode)
      {
        code_sp (-j + 1, p_enode->scope + j + 1);
      }
      else
      {
        code_sp (-j + 1, dflt_end + j - 1);
      }
    }
    list_DetachNodeByAddr (*blist, l_bnode);
    listNode_Destroy (l_bnode);
    l_bnode = p_bnode;

    if (p_enode)
      l_enode = p_enode;
    else
      l_enode = 0;
   }

  /* 
   * Chain the case statements:
   *
   */
  // this is first case statement
  l_cnode = list_GetLastNode (*clist);
  l_enode = list_GetLastNode (*elist);
  while (l_cnode)
  {
    p_cnode = listNode_GetPrevNode (l_cnode);
    p_enode = listNode_GetPrevNode (l_enode);

    // chain the cases
    if (p_cnode)
    {
      //    end of case C1 -> beginning of case C2
      code_sp (l_enode->scope + 1, p_cnode->scope - l_enode->scope-1);
    }
    else
    {
      // after last case, jump is either to the end or to the default statement
      code_sp (l_enode->scope + 1, dflt_end - l_enode->scope-1);
    }

    list_DetachNodeByAddr (*clist, l_cnode);
    listNode_Destroy (l_cnode);
    l_cnode = p_cnode;

    list_DetachNodeByAddr (*elist, l_enode);
    listNode_Destroy (l_enode);
    l_enode = p_enode;
  }

  // kilt the list at the end of the loop or switch statement
  list_Destroy(*clist);
  *clist=0;
  list_Destroy(*elist);
  *elist=0;
  list_Destroy(*blist);
  *blist=0;
}


static void resolve_continue_tags (int begin_op, int end_op, int tloop)
{
  ListNode *lnode, *pnode;

  if (!clist)
    return;

  lnode = list_GetLastNode (clist);

  /* 
   * See if there are any continues in the current loop.
   */

  while (lnode)
  {
    pnode = listNode_GetPrevNode (lnode);
    if (lnode->scope > begin_op && lnode->scope < end_op)
    {
      if (tloop)
        code_sp (lnode->scope + 1, end_op - lnode->scope - 1); /* for */
      else
        code_sp (lnode->scope + 1, begin_op - lnode->scope - 1); /* while */

      list_DetachNodeByAddr (clist, lnode);
      listNode_Destroy (lnode);
    }
    lnode = pnode;
  }
}

static void resolve_return_tags (int bop, int eop)
{
  ListNode *lnode, *pnode;

  if (!rlist)
    return;

  lnode = list_GetLastNode (rlist);

  /* 
   * See if there are any returns in the current loop.
   */

  while (lnode)
  {
    pnode = listNode_GetPrevNode (lnode);
    if (lnode->scope > bop && lnode->scope < eop)
    {
      code_sp(lnode->scope, OP_FOR_LOOP_DONE_RET);
      list_DetachNodeByAddr (rlist, lnode);
      listNode_Destroy (lnode);
    }
    lnode = pnode;
  }
}
