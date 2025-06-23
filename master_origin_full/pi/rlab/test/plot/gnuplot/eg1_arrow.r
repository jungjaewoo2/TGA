//
//
//
gnuwins (1);

old_debug = gnudefault("debug");
gnudefault("debug",1);



gnuwin (1);
gnulimits (-2,2,-3,4);
gnuxtics  (1,5);
gnuytics  (1,5);

_as = "set style arrow " + text(range(gnudefault().style_arrow)') ...
  + " " + gnudefault().style_arrow';
gnucmd (_as);
_rs = "as "+ text(range(gnudefault().style_arrow)');

cs = members(gnudefault("color"));
nc = length(cs);


gnuarrow([0,0],[-1,-1],_rs[1] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[-1, 0],_rs[2] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[-1, 1],_rs[3] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[ 0,-1],_rs[4] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[ 1,-1],_rs[5] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[ 1, 0],_rs[6] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[ 1, 1],_rs[7] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuarrow([0,0],[ 0, 1],_rs[8] + " lc rgb '" + cs[ int(uniform() * nc) + 1 ]+ "'");
gnuplot (,"arrows.eps");

gnudefault("debug",old_debug);

