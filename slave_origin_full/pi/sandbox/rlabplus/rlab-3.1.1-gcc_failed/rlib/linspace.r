//-------------------------------------------------------------------
//
//  Synopsis:   Generate a linearly spaced vector of size 'n'.
//
//  Syntax:	linspace ( s1 , s2 )
//		linspace ( s1 , s2 , n )
//
//-------------------------------------------------------------------


linspace = function ( x1 , x2, n )
{
  if(!exist(x1)||!exist(x2))
  {
    printf("linspace: Creates a uniformly spaced mesh between the endpoints.\n");
    printf("linspace: Format:\n");
    printf("linspace:   y=linspace(x1,x2/,n/),\n");
    printf("linspace: where 'x1' and 'x2' are the endpoints and 'n' is the size of\n");
    printf("linspace: the mesh so that  y[1]=x1  and  y[n]=x2 .\n");
    error("requires two or three arguments");
  }

  if (!exist (n))
  { n = 100; }

  if (n < 2)
  { error ("linspace: N must be >= 2"); }

  return [x1+[0:n-2]'*(x2-x1)/(n-1);x2];
};
