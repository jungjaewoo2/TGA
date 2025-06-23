//-------------------------------------------------------------------//

// Synopsis:    Solve the equation A*x = b using the Cholesky factorization

// Syntax:      cholsub ( U, b )
//
// Description:
//
//        The cholsub function computes the solution to the set of
//        linear equations described by:
//
//                A * x = b
//
//        The 1st argument to cholsub (U) is the result from
//        `chol(A)'. The second argument to cholsub is the matrix b. 
//        b can contain multiple right hand sides.
//        Matrix A must be a real symmetric positive definite, or
///       complex Hermitian positive definite matrix.
//        U is the Cholesky decomposition of A, i.e. A=U'*U.
//
//        Cholsub returns a matrix x which contains the solution(s) to
//        the aforementioned equations.
//
//
//        Example:  > A = [2,-1;-1,2]
//                   A =
//                          2         -1
//                         -1          2
//                  > b = [1,2]'
//                   b =
//                          1
//                          2
//                  > x = cholsub(chol(A), b)
//                   x =
//                       1.33
//                       1.67
//
// See Also: chol, factor, backsub, lu, solve
//
// T.-S. Yang (yang@isec.com)  5/3/96

//-------------------------------------------------------------------//

cholsub = function(U, b)
{
    x = b;
    n = x.nr;
    if (U.nr != U.nc) 
    {
       error("cholsub: the first argument must be square");
    }
    if (U.nr != b.nr)
    {
       error("cholsub: arguments with different number of rows")
    }
    
    // forward elimination
    x[1;] = x[1;]./U[1;1];
    for (i in 2:n)
    {
        x[i;] = (x[i;] - U[1:i-1;i]'*x[1:i-1;])./U[i;i];
    }
    
    // back substitution
    x[n;] = x[n;]./U[n;n];
    for (i in n-1:1:-1)
    {
       x[i;] = (x[i;] - U[i;i+1:n]*x[i+1:n;])./U[i;i];
    }
    return x;
};
