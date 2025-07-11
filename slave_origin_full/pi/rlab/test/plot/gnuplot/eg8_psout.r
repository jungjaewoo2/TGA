//
// main8.r: should produce a postscript file that differs only in few lines
// from a postscript file produced by execution of main8.gnu within gnuplot
//

gnuwins(1, "set term postscript enhanced color dashed");
cmds = ["set output \"gnuplot8_colorindex_rlab.ps\"", ...
    "set size 0.5,0.5", ...
    "set noborder", ...
    "set nokey", ...
    "set style line  1 linetype  1 linewidth 8", ...
    "set style line  2 linetype  2 linewidth 8", ...
    "set style line  3 linetype  3 linewidth 8", ...
    "set style line  4 linetype  4 linewidth 8", ...
    "set style line  5 linetype  5 linewidth 8", ...
    "set style line  6 linetype  6 linewidth 8", ...
    "set style line  7 linetype  7 linewidth 8", ...
    "set style line  8 linetype  8 linewidth 8", ...
    "set style line  9 linetype  9 linewidth 8", ...
    "set style line 10 linetype 10 linewidth 8", ...
    "unset xtics", ...
    "set ytics nomirror 1", ...
    "set yrange [ -1.5 : 10.5 ]", ...
    "plot  1 w l ls  1,  2 w l ls  2,  3 w l ls  3, 4 w l ls  4,\\", ...
    "5 w l ls  5,  6 w l ls  6,  7 w l ls  7, 8 w l ls  8,\\", ...
    "9 w l ls  9, 10 w l ls 10, -1 w line -1, 0 with line 0" ...
    ];
gnucmd (cmds);

gnuclose(1);

