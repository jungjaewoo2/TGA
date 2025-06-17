//-------------------------------------------------------------------//

//  Synopsis:   Show the variables in a list.

//  Syntax: who ( )
//              who ( LIST )

//  Description:

//  Who returns a STRING matrix of all the variables currently in the
//  global symbol table. Additionally who() will return the variable
//  list of an RLaB LIST when a LIST is provided as an argument. 
//
//                who ( $$ )
//  is the same as
//
//                who ()
//
//  Since `$$' is the special symbol that represents the global symbol
//  table.  
//
//  See Also: what

//-------------------------------------------------------------------//

static (NCOL)
NCOL = 5;

who = function ( list )
{
  if (!exist (list))
  {
#
# Setup the list variable...
#
    j = 1;
    dv = [""];
    for (i in members ($$))
    {
    //printf("%s %i\n", i, isprot($$.[i]));
      ctmp = class ($$.[i]);
      if (ctmp == "num" || ...
          ctmp == "string" || ...
          ctmp == "list")
      {
        if( isprot($$.[i]) ) { continue; }
        dv[j] = i;
        j = j + 1;
      }
    }
    } else {
#
# use list variable provided by user
#
          j = 1;
    dv = [""];
    for (i in members (list))
    {
      //printf("%s %i\n", i, isprot($$.[i]));
      ctmp = class (list.[i]);
      if (ctmp == "num" || ...
          ctmp == "string" || ...
          ctmp == "list")
      {
        if( isprot(list.[i]) ) { continue; }
        dv[j] = i;
        j = j + 1;
      }
    }
  }

#
# Now, make the matrix more presentable...
#
  rem = NCOL - mod (dv.n, NCOL);
  for (i in 1:rem)
  {
    dv = [dv, ""];
  }

  nrow = dv.n / NCOL;
  dv = reshape (dv, nrow, NCOL);
  return dv;
};
