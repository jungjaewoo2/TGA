//
//
//
gnuwins(2);
fn = "./data.txt";

data = readm(fn,3); // skip first three rows

xy = data[;2,3];

NR = 12;
NC = 8;

x = xy[;1];
y = xy[;2];

x_min = 1/NR * floor(NR * min(x)) - 0.5/NR;
x_max = 1/NR * ceil (NR * max(x)) + 0.5/NR;
y_min = 1/NC * floor(NC * min(y)) - 0.5/NC;
y_max = 1/NC * ceil (NC * max(y)) + 0.5/NC;

hxy = hist2(xy, [x_min:x_max:1/NR],[y_min:y_max:1/NC]);
1
gnuwin (1);
gnucmd ("reset");
gnulimits(x_min,x_max,y_min,y_max);
gnucmd("set pm3d map");
gnuxlabel ("X");
gnuylabel ("Y");
hist_norm(0);
if (hist_norm())
{
  gnutext   ("Frequency",[0.86,0.87],1,"screen");
else
  gnutext   ("Count of "+num2str(length(x)),[0.86,0.87],1,"screen");
}
gnusplot(hxy);
hist_norm(1);
2
d = 0.5;
c = cluster.pao(xy, d);
3
c_data = <<>>;
for (i in range(c.size))
{
  _n = find(c.feature == i);
  if (isempty(_n))
  { continue; }
  c_data.[i] = xy[_n;];
}
4
// 2-D histogram: (r0,c0)
// scatter plot:
gnuwin (2);
gnucmd ("reset");
gnulimits(x_min,x_max,y_min,y_max);
gnulegend (["Clusters : ",""] + members(c_data) + num2str(c.size'," (%3.0f)"));
gnuxlabel ("X");
gnuylabel ("Y");
gnuformat ([...
    "with points ps 2 pt 1 lc rgb 'red'", ...
    "with points ps 2 pt 1 lc rgb 'green'", ...
    "with points ps 2 pt 1 lc rgb 'blue'", ...
    "with points ps 2 pt 1 lc rgb 'orange'", ...
    "with points ps 2 pt 1 lc rgb 'brown'", ...
    "with points ps 2 pt 1 lc rgb 'pink'", ...
    "with points ps 2 pt 1 lc rgb 'magenta'", ...
    "with points ps 2 pt 1 lc rgb 'grey'", ...
    "with points ps 2 pt 1 lc rgb 'black'", ...
[]]);
gnuplot ( c_data, "cluster_" + text(d,"d=%04.2f") + ".eps" );

5

