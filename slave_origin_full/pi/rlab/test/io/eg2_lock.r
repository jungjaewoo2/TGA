//
//
//

rng(1,"uniform",[-1,1]);

x = 10 * rand(10,2);

opts = <<>>;
opts.eol = "\r\n";
opts.csp = ",";

opts.format = ["%.6f", "%8.2f"];
fn = "test1.txt";
lock(fn);
writem(fn, x, opts);

