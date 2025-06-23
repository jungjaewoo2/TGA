//
// test for determinant of sparse matrix
//


NITER=1000;

for(i in 1:NITER){

    x = rand(3,3);
    b = rand(3,1);
    det(x)
    sx = sparse( x );
    det(sx)

    y = rand(3,3) + 1i*rand(3,3);
    det(y)

    sy = sparse( y );
    sx=solve(sy,b);
    x =solve(y,b);
    det(sy)
    }
x
sx
