//-------------------------------------------------------------------
//
//  Synopsis:   Generates a logarithmically spaced vector in range
//              10^x1..10^x2 of size 'n'.
//
//  Syntax:	logspace ( x1 , x2 )
//		logspace ( x1 , x2 , n )
//
//-------------------------------------------------------------------

logspace = function(x1, x2, n)
{
  if(!exist(x1)||!exist(x2)){
    printf("logspace: Creates a logarithmically spaced mesh between the endpoints.\n");
    printf("logspace: Format:\n");
    printf("logspace:   y=logspace(x1,x2/,n/),\n");
    printf("logspace: where '10^x1' and '10^x2' are the endpoints and 'n'\n");
    printf("logspace:is the size of mesh.\n");
    error("requires three arguments");
  }
  if(!exist (n)) { n = 50; }

  if(n<2)
  { error ("logspace: third argument must be >= 2"); }

  if(!isreal(x1))
  { error("logspace: first argument must be real"); }

  if(!isreal(x2))
  { error("logspace: second argument must be real"); }

  xval = x1+(x2-x1)/(n-1)*[0:n-1]';

  return 10.^(xval);
};
