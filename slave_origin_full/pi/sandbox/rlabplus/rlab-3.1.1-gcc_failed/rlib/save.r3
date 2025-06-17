//-------------------------------------------------------------------//

//  Synopsis:   Save the workspace contents.

//  Syntax:	save ( )
//		save ( FILE )

//  Description:

//  The save function writes the contents of all the workspace
//  variables to a file. The default file, if none is specified is
//  "SAVE".
//

//-------------------------------------------------------------------//

static(VERBOSE);
VERBOSE=0;

save = function ( FILE, iverb )
{
  if (!exist (FILE))
  {
    FILE = "rlab.save";
  }
  if (exist (iverb))
  {
    VERBOSE = 1;
  } else {
    VERBOSE = 0;
  }


  open (FILE, "wb");

  if (VERBOSE)
  { printf("Saving global workspace to file %s\n", FILE); }
  for (i in members ($$))
  {
    if(!exist($$.[i]))
    { continue; }
    if (isprot($$.[i]))
    { continue; }
    if (class ($$.[i]) != "function")
    {
      if (VERBOSE)
      {
        printf("  %s", i);
        if (exist($$.[i].type))
        {
          printf("\t%s",$$.[i].type);
        } else {
          printf("\t%s","list");
        }
        printf("\n");
      }
      writeb (FILE, $$.[i]);
    }
  }

  close (FILE);

  if (VERBOSE)
  { printf("Saving Done!\n\n"); }

};

