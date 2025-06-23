//-------------------------------------------------------------------//

//  Synopsis:   Print a description of the variables in a list.

//  Syntax:  whos ( L )

//  Description:

//  The whos function prints a tabular listing (to stdout) of the
//  contents of the list, or symbol-table L. Information on functions
//  is not output. If no argument is provided whos prints out
//  information from the global-symbol table.

//  See Also: sizeof, who, what

//-------------------------------------------------------------------//

static (whos_1, Btotal);

whos = function ( LIST )
{
  Btotal = 0;
  if (!exist (LIST))
  {
    whos_1 ();
  } else {
    whos_1 (LIST);
  }
  printf ("Total MBytes = %f\n", Btotal/1.e6);
};

whos_1 = function ( LIST )
{
  if(!exist(LIST))
  {
    printf ("\tName            Class\tType\tSize\t\tNBytes\n");
    for (i in members ($$))
    {
      nbytes = sizeof ($$.[i]);
      Btotal = Btotal + nbytes;
      if (class ($$.[i]) == "function")
      { continue; }
      if (class ($$.[i]) == "list")
      {
        if (isprot($$.[i])) { continue; }
        m = size ($$.[i]);
        printf ("\t%-15s", i);
        printf ("\t%s\t%s\t%i\t\t%i\n", ...
            class ($$.[i]), type ($$.[i]), m[1], nbytes);
        } else {
          m = size ($$.[i]);
        printf ("\t%-15s", i);
        printf ("\t%s\t%s\t%i\t%i\t%i\n", ...
            class ($$.[i]), $$.[i].type, m[1], m[2], nbytes);
      }
    }
  } else {
    printf ("\tName            Class\tType\tSize\t\tNBytes\n");
    for (i in members (LIST))
    {
      nbytes = sizeof (LIST.[i]);
      Btotal = Btotal + nbytes;
      if (class (LIST.[i]) == "function")
      { continue }
      if (class (LIST.[i]) == "list")
      {
        m = size (LIST.[i]);
        printf ("\t%-15s", i);
        printf ("\t%s\t%s\t%i\t\t%i\n", ...
              class (LIST.[i]), type (LIST.[i]), m[1], nbytes);
      } else {
        m = size (LIST.[i]);
        printf ("\t%-15s", i);
        printf ("\t%s\t%s\t%i\t%i\t%i\n", ...
              class (LIST.[i]), type (LIST.[i]), m[1], m[2], nbytes);
      }
    }
  }
};
