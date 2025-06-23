//-------------------------------------------------------------------//

//  Synopsis:   Visualize the sparsity structure.

//  Syntax:     spy ( S, mark )

//  Description:

//  Visualize the sparsITY structure.
//  spy(S) plots the sparsity pattern of any matrix S.
//  spy(S,mark) uses the specified marker (1-8).
//
//-------------------------------------------------------------------//

spy = function( S, mark )     
{
  if (max(S.nr,S.nc) <= 1) {
     error("spy: cannot plot a vector or a scalar");
  }
  if (S.type == "complex") {
    S = real (S);
  }
  if (S.storage == "dense") {
     S3 = spconvert(sparse(S)); 
  } else { 
     S3 = spconvert(S);
  }
  S3[;1] = (S.nr+1)*ones(S3.nr,1) - S3[;1];
  plimits(1,S.nc,1,S.nr);
  plstyle("point");
  if (exist(mark)) {
     mark = mod(mark,8);
     if(mark==0) {mark=8;}
     pvec = [mark:8,1:mark-1];
     plpoint(pvec);
  }
  pltitle("Matrix Non-zero Entries");
  plwid(4); plegend();
  plot( [ S3[;2], S3[;1] ] );
  // reset
  pltitle();
  plimits();
  plstyle();
  plpoint();
  plimits();
  plwid();
};
