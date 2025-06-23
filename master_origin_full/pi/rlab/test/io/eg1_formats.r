//
//
//

rng(1,"uniform",[-1,1]);

x = 10 * rand(10,2);

opts = <<>>;
opts.eol = "\r\n";
opts.csp = ",";

opts.format = ["%.6f", "%8.2f"];
fn = "./test1.txt";
writem(fn, x, opts);

opts.format = ["%12s", "%5s"];
writem("./test2.txt", text(x), opts);

opts.format = ["%.6f", "%.12f"];
writem("./test3.txt", x + 1i*rand(x), opts);

