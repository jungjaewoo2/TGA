//
// generating words from english language using aspell
//
static(INIT, WHICH, THIS_LIBRARY, ASPELL);

if (exist(INIT))
{ EOF }

THIS_LIBRARY = "libwords";
if (!exist(WHICH))
{
  WHICH=reads("|which which");
  if (strlen(WHICH)<strlen("which"))
  {
    clear(WHICH);
    printf("Function 'which' is not available: library %s cannot be initialized\n", ...
        THIS_LIBRARY);
    EOF
  }
}

// check for 'aspell'
ASPELL = reads("|" + WHICH + " aspell 2>/dev/null");
if (isempty(ASPELL))
{ clear(ASPELL); }
if (strlen(ASPELL)<strlen("aspell"))
{ clear(ASPELL);}

//
// generate matrix of random characters 
//
randchar = function(nr, nc, l, ctype)
{
  if (!exist(ctype))
  { ctype="alphanumeric"; }
  if (!exist( l ))
  { l  = 7; }
  if (l < 1)
  { l  = 7; }
  if (ctype=="alphanumeric")
  { ascint=[48:57,65:90,97:122]; }
  if (ctype=="printable")
  { ascint=[33:126,161:255]; }
  if (ctype=="alphabet" || ctype=="text" || ctype=="letters")
  { ascint=[65:90,97:122]; }
  if (ctype=="numeric" || ctype=="digits")
  { ascint=[48:57]; }

  noasc=length(ascint);
  rndc=blank(nr,nc);
  for(i in 1:nr)
  {
    for (j in 1:nc)
    {
      rndc[i;j] = char( ascint[ 1 + int(noasc*uniform(1,l)) ] );
    }
  }
  return rndc;
};

//
// pick random word from default aspell dictionary
//
randomword = function(count)
{
  global(randchar);

  if (!exist(count))
  { count = 1; }

  rval = blank(1,count);

  j = 0;
  while (j < count)
  {
    // do these inspire 'aspell' to find a word matching it?
    if (!exist(ASPELL))
    {
      len = 3 + floor(12 * uniform());
      rval[j] = tolower(randchar(1,1,len,"alphabet"));
    }

    // generate 8 random characters
    a = tolower(randchar(1,1,8,"alphabet"));

    cmd = "| echo '"+ a +"' | " + ASPELL + " -a 2>/dev/null";
    r = sum(reads(cmd));
    r = strsplt(r, a);
    n = strsplt(r[2], ":");

      // if n[1]==0 then try generate another 8-random character word
    idx = strtod(n[1],<<lstrip="'BLANK">>);
    if (length(idx)==1)
    { continue; }

      // aspell is inspired: n[2] is the word list of length of idx[1]
    j++;
    w = strsplt(n[2], ",");
    chomp(w);
    rval[j] = shuffle(w,1);
  }

  //remove apostrophe from words
  rval = gsub("'", rval).string;
  return tolower(rval);
};

//
// use for pattern:
//  {w,n}'some fixed text'{w,n}
// where 'w' is for random word, and 'n' is for random number
// if the fixed text does not contain letters w or n then apostrophes are not needed
randomlabel = function(pattern, options)
{
  if (!exist(pattern))
  { return blank(); }
  if (type(pattern)!="string")
  { return blank(); }

  // default values
  fix_n_digits = 0;
  leading_zeros = 0;
  max_n_digits = 4;
  n_digits = 0;
  fmt = "%.0f";

  // user can modify them
  if (class(options)=="list")
  {
    if (exist(options.n_digits))
    { max_n_digits = n_digits = options.n_digits; }
    if (exist(options.max_n_digits))
    { max_n_digits = options.max_n_digits; }

    if (exist(options.leading_zeros))
    {
      if (options.leading_zeros)
      {
        fmt = "%0"+num2str(floor(max_n_digits),"%.0f")+".0f";
        leading_zeros = 1;
      }
    }
  }

  global(randomword);

  p = strsplt(shuffle(pattern,1));
  rval = "";
  toggle = 0;
  for (p1 in p)
  {
    switch(p1)
    {
      case "'":
        toggle = !toggle;
        break;

      case "w":
        if (!toggle)
        {
          rval = rval + randomword();
          break;
        }
      case "n":
        if (!toggle)
        {
          if (!fix_n_digits)
          { n_digits = 1 + floor(max_n_digits * uniform()); }

          if (leading_zeros)
          { rval = rval + num2str(floor(10^(n_digits-1)*uniform()),fmt); }
          if (!leading_zeros)
          { rval = rval + num2str((10^n_digits - 10^(n_digits-1))*uniform()+10^(n_digits-1),fmt); }
          break;
        }
      default:
        rval = rval + p1;
        break;
    }
  }
  return rval;
};


INIT = 1;


