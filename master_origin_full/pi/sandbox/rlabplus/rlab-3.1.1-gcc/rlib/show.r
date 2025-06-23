//-------------------------------------------------------------------//

//  Synopsis:   Show information about a variable.

//  Syntax:     show ( A )

//  Description:

// 	Show takes a single argument and returns a brief description
// 	of the argument, it's class, name, and any other pertinent
// 	characteristics.
 
// 	Note that all the information provided by show() can be
// 	obtained by other means. class(), type(), and size().
// 	Additionally, direct member reference will provide the same
// 	information ss show.

// See Also: class, size, type

//-------------------------------------------------------------------//

show = function ( A )
{
  if (!exist (A))
  {
    return ("\tUNDEFINED");
  }

  if (class (A) != "list")
  {
    for (i in members (A))
    {
      if (class (A.[i]) == "string")
      {
	printf ("\t%-20s:\t%s\n", i, A.[i]);
      else if (class (A.[i]) == "num") {
	if (A.[i].n == 1)
	{
	  printf ("\t%-20s:\t%g\n", i, A.[i]);
	else
	  printf ("\t%-20s:\t%s, %s, %s, %ix%i\n", i, A.[i].class, ...
	           A.[i].type, A.[i].storage, A.[i].nr, A.[i].nc);
      } } }
    }
  else
    for (i in members (A))
    {
      printf ("\t%-20s:%s\t%s\n", i, class (A.[i]), type(A.[i]));
    }
  }
};
