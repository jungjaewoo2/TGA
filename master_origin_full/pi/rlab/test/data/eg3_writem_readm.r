//
//
//
fn = "./test_writem.csv";

nr = int(100 * uniform()) + 1;
nc = int(100 * uniform()) + 1;

x = rand(nr, nc);
x = ifelse(x>0.9, nan(), x);
x = ifelse(x<0.1, inf(), x);
x = ifelse(x<0.2, -inf(), x);

//
// write to file
//
opts_writem = <<>>;
opts_writem.format = "%.3f";
opts_writem.csp = ",";
opts_writem.eol = "\r\n";
opts_writem.nan="NaN";
opts_writem.inf_pos="inf";
opts_writem.inf_neg="-inf";
writem(fn, x, opts_writem);


//
// read from the same file
//
opts_readm = <<>>;
opts_readm.csp = ",";
y = readm(fn, opts_readm);

idx_num = find(finite(x));
printf("This should be unity -> %g, Is it?\n", ...
  sum(abs(floor(x[idx_num], <<bin=1e-3>>) - y[idx_num])<=1.1e-3) / length(idx_num));

