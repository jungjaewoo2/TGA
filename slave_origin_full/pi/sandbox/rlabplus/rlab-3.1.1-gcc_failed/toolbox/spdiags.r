//-------------------------------------------------------------------//

// Synopsis: Sparse matrix formed from diagonals.
//           spdiags, which generalizes the function "diag", deals with 
//           three matrices, in various combinations, as both input and 
//           output.

// Syntax:   </ B ; d /> = spdiags(A)
//           B = spdiags(A,d)
//           A = spdiags(B,d,A)
//           A = spdiags(B,d,m,n)

// Description:

//    </ B ; d /> = spdiags(A) extracts all nonzero diagonals from the 
//    m-by-n matrix A.  B is a min(m,n)-by-p matrix whose columns are 
//    the p nonzero diagonals of A.  d is a vector of length p whose 
//    integer components specify the diagonals in A.
// 
//    B = spdiags(A,d) extracts the diagonals specified by d.
//    A = spdiags(B,d,A) replaces the diagonals of A specified by d 
//        with the columns of B.  The output is sparse.
//    A = spdiags(B,d,m,n) creates an m-by-n sparse matrix from the
//        columns of B and places them along the diagonals specified 
//        by d.
// 
//    Roughly, A, B and d are related by
//        for (k in 1:p) {
//            B[;k] = diag(A,d[k]);
//        }
//    Some elements of B, corresponding to positions "outside" of A, 
//    are not actually used.  They are not referenced when B is input 
//    and are set to zero when B is output.
// 
//    Example: These commands generate a sparse tridiagonal 
//             representation of the classic second difference operator 
//             on n points.

//             e = ones(n,1);
//             A = spdiags([e,-2*e,e], -1:1, n, n)
// 
//    See also diag()
//
//-------------------------------------------------------------------//

spdiags = function(a1, a2, a3, a4)
{
  if (nargs == 0) 
  {
     error("spdiags: needs at least 1 argument");
  else if (nargs == 1)
  {
        // [B,d] = spdiags(A)
        m = a1.nr;
        n = a1.nc;
        p = min(m,n);
        d = [];
        B[p;1] = 0;
        i = 0;
        for (k in -(m-1):(n-1))
        {   
            b = full(diag(a1,k));
            if (any(b))
            {
               i = i + 1;
               d[i] = k;
               if (k >= 0)
               {
                  for (j in 1:b.n)
                  {
                      B[j+k;i] = b[j];
                  }
               else
                  for (j in 1:b.n)
                  {
                      B[j;i] = b[j];
                  }
               }
            }
        }
        if (d==[]) { B=[]; }
        return <<B;d[:]>>;
  else if (nargs == 2)
  {
        // B = spdiags(A,d)
        B[min(a1.nr,a1.nc);1] = 0;
        for (i in 1:a2.n)
        {
            k = a2[i];
            b = full(diag(a1,k));
            if (k >= 0)
            {
               for (j in 1:b.n)
               {
                   B[j+k;i] = b[j];
               }
            else
               for (j in 1:b.n)
               {
                   B[j;i] = b[j];
               }
            }
        }
        return B;  
  else if (nargs == 3)
  {
        // A = spdiags(B,d,A)
        m = a3.nr;
        n = a3.nc;
        if (a3.type != a1.type)
        {
           error("spdiags: data type of 1st and 3rd argument must be the same");
        }
        // costly!
        #for (i in 1:a2.n)
        #{
        #    k = a2[i];
        #    if (k >= 0)
        #    {
        #       for (j in 1:a1.nr)
        #       {
        #          if (j>m || (j+k)>n) { break;}
        #          a3[j;j+k] = a1[j+k;i];
        #       }
        #    else
        #       for (j in 1:a1.nr)
        #       {
        #          if (j>m || (j-k)>n) { break;}
        #          a3[j-k;j] = a1[j;i];
        #       }               
        #    }
        #}
        // cheaper
        if (a3.type == "real") 
        {
           tmp = [m,n,0];
           for (i in 1:a2.n)
           {
             k = a2[i];
             if (k >= 0)
             {
               for (j in 1:a1.nr)
               {
                  if (j>m || (j+k)>n) { break;}
                  tmp = [tmp;
                         j,j+k,-a3[j;j+k];
                         j,j+k,a1[j+k;i]];
               }
             else
               for (j in 1:a1.nr)
               {
                  if (j>m || (j-k)>n) { break;}
                  tmp = [tmp;
                         j-k,j,-a3[j-k;j];
                         j-k,j,a1[j;i]];
               }               
             }
           }           
        else
           // complex
           tmp = [m,n,0,0];        
           for (i in 1:a2.n)
           {
             k = a2[i];
             if (k >= 0)
             {
               for (j in 1:a1.nr)
               {
                  if (j>m || (j+k)>n) { break;}
                  tmp = [tmp;
                         j,j+k,-real(a3[j;j+k]),-imag(a3[j;j+k]);
                         j,j+k, real(a1[j+k;i]), imag(a1[j+k;i])];
               }
             else
               for (j in 1:a1.nr)
               {
                  if (j>m || (j-k)>n) { break;}
                  tmp = [tmp;
                         j-k,j,-real(a3[j-k;j]),-imag(a3[j-k;j]);
                         j-k,j, real(a1[j;i]),   imag(a1[j;i])   ];
               }               
             }
           }
        }     
        return (a3+spconvert(tmp));
  else if (nargs == 4)
  {
        // A = spdiags(B,d,m,n)        
        m = a3;
        n = a4;
        //if (a1.nc != a2.n)
        //{
        //   error("spdiags: incorrect sizes of B and d");
        //}
        if (a1.nr > min(m,n))
        {
           error("spdiags: row dimension of B too large");
        }
        if (a1.type == "complex") 
        {
           sa = [m,n,0,0];
           for (i in 1:a2.n)
           {
               k = a2[i];
               if (k >= 0)
               {
                  for (j in 1:a1.nr)
                  {
                     if (j>m || (j+k)>n) { break;}
                     sa = [sa;[j,j+k,real(a1[j+k;i]),imag(a1[j+k;i])]];
                  }
               else
                  for (j in 1:a1.nr)
                  {
                     if (j>m || (j-k)>n) { break;}
                     sa = [sa;[j-k,j,real(a1[j;i]),imag(a1[j;i])]];
                  }               
               }
           }
        else
           sa = [m,n,0];
           for (i in 1:a2.n)
           {
               k = a2[i];
               if (k >= 0)
               {
                  for (j in 1:a1.nr)
                  {
                     if (j>m || (j+k)>n) { break;}
                     sa = [sa;[j,j+k,a1[j+k;i]]];
                  }
               else
                  for (j in 1:a1.nr)
                  {
                     if (j>m || (j-k)>n) { break;}
                     sa = [sa;[j-k,j,a1[j;i]]];
                  }               
               }
           }
        }
        
        A = spconvert(sa);
        return A;
  else
  {
     error("spdiags: Too many arguments (max 4)");
  }}}}}}
};
