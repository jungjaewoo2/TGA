//
//
// 
filt = function (x)
{
    n = size(x);
    if (n[1] > 1) {
        x = x';
    }
    len = length(x);
    b=[1,-1]; a=[1,-0.99];
    nfilt = 2; nfact = 3;
    zi = max(size(a), size(b));
    x1 = x[4:2:-1];
    x2 = x[(len-1):(len-nfact):-1];
    y = [(2*x[1]-x1),x[;],(2*x[len]-x2)];
    lst = filter(b,a,y,zi*y[1]);
    y = lst.y;
    y = y[length(y):1:-1];
    lst = filter(b,a,y,zi*y[1]);
    y = lst.y;
    y = y[length(y):1:-1];
    y[1:nfact,(len+nfact+(1:nfact))] = [];
    return y;
};
a = filt([23:54:3])


    