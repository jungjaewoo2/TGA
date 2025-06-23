//---------------------------------------------------------------------------

//  Synopsis:   Compute the vector dot product.

//  Syntax:	dot ( A , B )

//  Description:

//  Compute the dot product of two vectors, A and B.
//---------------------------------------------------------------------------

dot = function ( A , B )
{
  if (min (size (A) != 1) || min (size (B) != 1)) 
  {
    error ("dot: A and B must be vectors");
  }
  return (A[:]' * B[:])
};
