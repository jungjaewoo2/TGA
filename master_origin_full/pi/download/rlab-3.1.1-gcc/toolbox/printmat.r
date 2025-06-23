//--------------------------------------------------------------------------------
// printmat

// Synopsis: Print a matrix with title and labels.

// Syntax: printmat(A, Aname, ROWLAB, COLLAB);

// Description:

// This routine prints the matrix A using the title contained in the
// string Aname and the row labels contained in ROWLAB and the
// column labels contained in COLLAB. Note: ROWLAB and COLLAB are
// vectors of strings, such as
//
//         rowlab[1]="alpha";
//         rowlab[2]="beta";
//         rowlab[3]="gamma";
//
//  Note: Aname, ROWLAB, and COLLAB are optional variables.
//
//--------------------------------------------------------------------------------

require num2str

printmat = function (a, mname, rowlab, collab, fn)
{
  // Set defaults, error check.

  if (!exist (a)) { error ("printmat: must supply a matrix"); }
  if (!exist (mname)) { mname = ""; }
  if (!exist (rowlab)) { rl = []; else rl = rowlab; }
  if (!exist (collab)) { cl = []; else cl = collab; }
  if (!exist (fn))     { fn = "stdout"; }
  
  if (!isempty (rl)) {
    if (rl.n != a.nr) { error ("printmat: ROWLAB.n != A.nr"); }
  }
  if (!isempty (cl)) {
    if (cl.n != a.nc) { error ("printmat: COLLAB.n != A.nc"); }
  }

  nrows=a.nr;
  ncols=a.nc;
  
  // Create row and column labels if necessary
  if (rl.n == 0)
  {
    for (i in 1:nrows) 
    {
      sprintf (tmp, "%3i", i);
      rl[i]="--"+ tmp +" --> ";
    }
  }
  if (cl.n == 0)
  {   
    for (i in 1:ncols) { cl[i]="----"+num2str(i)+"---- "; }
  }
  
  col_per_scrn=5;
  len=12;
  
  if ((nrows==0)||(ncols==0)) 
  { 
    if (length (mname)) 
    {
      fprintf(fn," \n%s = \n \n",mname);
      return 0;
    }
    fprintf(fn," \n%s \n","     [] \n");
    return 0;
  }

  // Print matrix name
  col=1;
  n = min([col_per_scrn-1,ncols-1]);
  if (length (mname)) 
  {
    fprintf(fn,"\n%s = \n \n",mname);
  }

  // Print column labels
  s="";
  icol=0;
  while (col <= ncols) 
  {
    icol=icol+1;
    s="            ";
    for (j in 0:n) 
    {
      ishift=13-length(cl[j+col]);
      for (k in 1:ishift) 
      {
	s=s+" ";
      }
      s=s+cl[j+col];
    }
    fprintf(fn,"%s\n",s);
    
    // Print Row Labels
    for (i in 1:nrows) 
    {
      s=""+rl[i];
      ishift=12-length(rl[i]);
      for (k in 1:ishift) 
      {
	s=s+" ";
      }
      for (j in 0:n) 
      {
	element = a[i;col+j];
	if (element == 0) {
	  s=s+"           0";
	else if (element >= 1.0e+06) {
	  sdum="";
	  sprintf(sdum," %12.5e",element);
	  s=s+sdum;
	else if (element <= -1.0e+05) {
	  sdum="";
	  sprintf(sdum," %12.5e",element);
	  s=s+sdum;
	else if (abs(element) < 0.0001) {
	  sdum="";
	  sprintf(sdum," %12.5e",element);
	  s=s+sdum;
	else
	  sdum="";
	  sprintf(sdum," %12.5f",element);
	  s=s+sdum;
        } } } }
      }
      fprintf(fn,"%s\n",s);
    }
    col = col+col_per_scrn;
    fprintf(fn,"%s"," \n");
    if ((ncols-col) < n) 
    {
      n=ncols-col;
    }
  }
  
  return 0;  
};
