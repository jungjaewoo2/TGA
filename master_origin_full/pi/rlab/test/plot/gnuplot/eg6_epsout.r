//
// file: main6.r
//

x = [0:2:1/16]';
y = [0:1:1/16]';
z = zeros(length(x),length(y));
for (i in 1:length(x))
{
  for (j in 1:length(y))
  { z[i;j] = x[i].^2 + 2*x[i] + y[j].^2; }
}

gnuwins(1, "eps");
gnucmd ("set output 'gnuplot6.eps'");
gnucmd ("set size 1.4,1");
gnuxlabel ("x");
gnuylabel ("y");
gnuzlabel ("{/Symbol F}(x,y)");
gnulimits (0,1,0,1);
gnuytics  (1/4, 5);
gnuxtics  (1/4, 5);
gnuztics  (1/2, 5);
gnulegend ( ["x^2+2*x+y^2 calculated inside rlab", ...
             "sin(pi*x)*sin(pi*y) by gnuplot" ...
             ] );
gnucmd ("set view 53,13" );
gnucmd ("set grid" );
gnuformat ("with lines");
data = <<>>;
  // first data set: from rlab
  data.[1].x = x;
  data.[1].y = y;
  data.[1].z = z;
  // second data set: a function calculated by gnuplot
  data.[2] = "sin(pi*x)*sin(pi*y)";
gnusplot(data);

gnuclose(1);


if (exist(cell))
{
  z_levels = [0.1,0.2,0.4,0.9];
  gnuwins(1);
  a=conrec(z_levels,z',x,y);
  gnulegend(["Contour levels with z = ", ""] + text(z_levels,"%.1f"));
  gnuformat(["with points","with points","with points","with points"]);
  gnuplot(<<a1=a[1];a2=a[2];a3=a[3];a4=a[4]>>);
}

