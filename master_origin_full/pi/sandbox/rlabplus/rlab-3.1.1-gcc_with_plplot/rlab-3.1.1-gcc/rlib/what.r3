//-------------------------------------------------------------------//

//  Synopsis:   Show what functions are active in the workspace.

//  Syntax:	what ( )
//              what ( LIST )

//  Description:
//
//  Show what functions are in a list. What returns a string matrix
//  all of the functions in a list. If no arguments are supplied (most
//  common usage) then the names of the functions in the global
//  workspace are returned.

//  See Also: who

//-------------------------------------------------------------------//

static (NCOL)
NCOL = 5;

what = function ( list )
{

  #
  # Setup the list variable...
  #

  if (!exist (list))
  {
    j = 1;
    dv = [""];
    for (i in members ($$))
    {
      ctmp = class ($$.[i]);
      if (ctmp == "function")
      {
        if (strsplt (i)[1] != "_")
        {
          dv[j] = i;
          j = j + 1;
        }
      }
    }
  } else {
    #
    # Get the data-object names for list...
    #
    j = 1;
    dv = [""];
    for (i in members (list))
    {
      ctmp = class (list.[i]);
      if (ctmp == "function")
      {
        if (strsplt (i)[1] != "_")
        {
          dv[j] = i;
          j = j + 1;
        }
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
