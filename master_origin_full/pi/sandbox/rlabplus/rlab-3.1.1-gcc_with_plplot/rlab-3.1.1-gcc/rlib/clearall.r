//-------------------------------------------------------------------

//  Synopsis:	Clear all the variables and user-functions from the workspace.

//  Syntax:	clearall()

//  Description:

//  The function clearall, clears all user data objects from the workspace:
//  Scalars, strings, matrices and lists. User-functions have to be cleared
//  using clear() function.


//  See Also: clear
//-------------------------------------------------------------------


clearall = function ( )
{
  for (i in members ($$))
  {
    if (class ($$.[i]) != "function" && class ($$.[i]) != "list")
    {
      if (i != "pi" && i != "eps" && i != "_rlab_search_path" && i != "_rlab_config")
      { clear ($$.[i]); }
    }
    if(class ($$.[i]) == "list")
    {
      if(isprot($$.[i]))
      { continue; }
      // element is a list. check whether there are functions in its member list.
      // assume single depth.
//       if (i == "_rlab_search_path" || i == "_rlab_config")
//       { continue; }
//       if (i == "mks" || i == "const")
//       { continue; }
      for(j in members($$.[i]))
      {
        if (class($$.[i].[j])!="function")
        {
          clear ($$.[i].[j]);
        else if(type($$.[i].[j])=="user")
        {
          clear ($$.[i].[j]);
        }}
      }
      if(members($$.[i])==[])
      {clear ($$.[i]);}
    }
  }
};

