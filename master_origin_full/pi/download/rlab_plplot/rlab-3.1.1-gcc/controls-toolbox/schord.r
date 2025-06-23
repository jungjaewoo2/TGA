//----------------------------------------------------------------------
//
// schord
//
// Syntax: </Qo; To/> = schord(Qi, Ti, index)
//
//	Ordered schur decomposition.
//	</Qo; To/> = schord(Qi, Ti, index)  Given the square (possibly 
//	complex) upper-triangular matrix Ti and orthogonal matrix Qi
//	(as output, for example, from the function SCHUR),
//	SCHORD finds an orthogonal matrix Q so that the eigenvalues
//	appearing on the diagonal of Ti are ordered according to the
//	increasing values of the array index where the i-th element
//	of index corresponds to the eigenvalue appearing as the
//	element  Ti(i,i).
//
//	The input argument list may be [Qi, Ti, index] or [Ti, index].
//	The output list may be [Qo, To] or [To].
//
//	The orthogonal matrix Q is accumulated in Qo as  Qo = Qi*Q.
//     	If Qi is not given it is set to the identity matrix.
//
//          *** WARNING: SCHORD will not reorder REAL Schur forms.
//
//----------------------------------------------------------------------
schord = function(Qi, Ti, index)
{
  n = Ti.nr;
  if (Ti.nr != Ti.nc){ error("Nonsquare Ti"); }
  if (nargs == 2)
  {
    index = Ti; 
    To = Qi;
  else
    To = Ti;
  }
  
  if (nargs > 2) {
     Qo = Qi; 
  else 
     Qo = eye(n,n);
  }
  //
  for (i in 1:(n-1))
  {
    // -- find following element with smaller value of index --
    move = 0; 
    k = i;
    for (j in (i+1):n) {
	if (index[j] < index[k]) {
	   k = j; 
	   move = 1; 
	}
    }

    // -- propagate eigenvalue up the diagonal from k-th to i-th entry --
    if (move)
    {
        for (l in k:(i+1):-1)
        {
	    l1 = l-1; 
	    l2 = l;
	    t = givens(To[l1;l1]-To[l2;l2], To[l1;l2]);
	    t = [t[2;];t[1;]];
	    To[;l1:l2] = To[;l1:l2]*t; 
	    To[l1:l2;] = t'*To[l1:l2;];
	    Qo[;l1:l2] = Qo[;l1:l2]*t;
	    ix = index[l1]; 
	    index[l1] = index[l2]; 
	    index[l2] = ix;
        }
    }
  }
  return << Qo=Qo; To=To>>;
};

