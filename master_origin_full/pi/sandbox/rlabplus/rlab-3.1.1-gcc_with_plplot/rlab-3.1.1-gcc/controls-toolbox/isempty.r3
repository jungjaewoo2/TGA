//------------------------------------------------------------------
//
// isempty(A)
//
// Syntax: isempty(A)
//
//         Returns 1 if A is empty.
//         Returns 0 if the matrix has any elements.
//
//         isempty will work on list-objects as well as matrices.
//
// Original Author: Jeff Layton 
//------------------------------------------------------------------

isempty = function(A)
{
  return !all (size (A));
};

