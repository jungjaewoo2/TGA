// pj: chaos
// playing with financial markets

//
// get the data
//
dir   = "../data/";
fname = "gspc_dohlcv.csv";
sx = reads( dir + fname );
// x[1:;1:6] :
//  x[;1] -> date as a number in YYYYMMDD format
//  x[;2] -> open
//  x[;3] -> high
//  x[;4] -> low
//  x[;5] -> close
x  = strtod(sx);

// invert
x = x[x.nr:1:-1;];
// change in closing price (day - previous day)
dc = x[2:x.nr;2] - x[1:(x.nr-1);2] ;
// daily volatility: (hi - lo)
vl = x[;3] - x[;4];

if(!exist(y)) { y = <<>>; }

NMAX = 6000;

//
// average mutual information of change in the closing price
//
n  = [1:NMAX];
an = ami(dc, n);
y.[1] = [n', an'];
av = ami( vl, n);
y.[2] = [n', av'];

//
// autocorrelation function
//
n  = [1:NMAX];
mdc = mean(dc);
an = xcorr(dc-mdc, dc-mdc, n, "c");
y.[3] = [n', an'];
av = xcorr(vl, vl, n, "c");
y.[4] = [n', av'];


//
// false nearest neighbours
//
terminal = stderr();
embd = [1:20];
y.[5]=[embd', falsenn(vl, embd, 1, terminal)'];
y.[6]=[embd', falsenn(dc, embd, 1, terminal)'];


rfile module_xmgr
