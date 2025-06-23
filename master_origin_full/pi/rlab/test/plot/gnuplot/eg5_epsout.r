//
// file: main5.r
//

gnuwins(1, "eps");

gnucmd ("set output 'gnuplot5.eps'");
gnutitle  ("Making of a strawberry-rhubarb pie");
gnuxlabel ("Parameter {/Symbol strawberry} (strawberry)");
gnuylabel ("Parameter {/Symbol rhubarb} (rhubarb)");
gnulimits (0,1.5,0,1.5);
gnuytics  (0.1, 2);
gnuxtics  (0.1, 2);
gnucmd("unset grid; set grid xtics ytics;");
gnuformat (["using ($2):($1) with lines", ...
    "using ($2):($1) with lines", ...
    "with lines", ...
    "with points", ...
    "using 1:($1+$2<1.3 ? $2 : 1/0) with points "]);
gnulegend ( [ ...
    "Gnuplot-loaded file \"tryme1.dat\"", ...
    "Gnuplot-loaded file \"tryme2.dat\"", ...
    "Gnuplot-evaluated expression y = 1.3-x", ...
    "RLaB-generated random points", ...
    "Gnuplot-filtered RLaB-generated random points" ...
] );
data = <<>>;
data.[1] = "tryme1.dat";  // file containing data to be plotted
data.[2] = "tryme2.dat";  // the same
data.[3] = "1.3-x";      // gnuplot expression we want plotted

rng(1, "uniform", [0,1.5]);
rng(2, "uniform", [0,1.11]);

rng(1);
x = rand (1000,1);
rng(2);
y = rand (1000,1);

t = [x, y];
data.[4] = t;  // some rlab generated data
data.[5] = t;  // some rlab generated data to be filtered by gnuplot

//gnuplot(data, "gnuplot5.eps", "postscript enh color");
gnuplot(data);

gnuclose(1);


