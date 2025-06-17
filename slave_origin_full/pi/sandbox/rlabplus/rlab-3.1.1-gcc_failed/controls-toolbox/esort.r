//------------------------------------------------------------------------------
//
// esort
//
// Syntax: a=esort(p)
//
// This routine sorts the complex continuous eigenvalues in descending
// order. It sorts based on the real part. The unstable eigenvalues (positive
// real part) will appear first.
//
// The routine passes back the sorted eigenvalues and the cooresponding
// indices in a list:
//
//    a.s   = sorted eigenvalues
//    a.ndx = index
//
//------------------------------------------------------------------------------

esort = function(p)
{
   t = sort(-real(p));
   return << s=p[t.idx]; ndx=t.idx >>;
};

