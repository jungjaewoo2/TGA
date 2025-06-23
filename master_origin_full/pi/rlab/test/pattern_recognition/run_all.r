//
//
//

// fns = ls("./eg*.r");
fns = [...
  "eg1_cluster.r", ...
    "eg2_cluster.r", ...
        "eg3_cluster.r", ...
            "eg4_cluster.r", ...
                []];


_JMAX = (int(1000*uniform()) + 1);
for (_jj in 1:_JMAX)
{
  spinner();
  //fn = shuffle(fns,1);
  fn = fns[ mod(_jj-1, length(fns))  + 1 ];
  //fn = fns[ mod(_jj, 4)  + 1 ];
  printf("%g/%g: executing script %s\n", _jj, _JMAX, fn);
  NITER = int(10*uniform()) + 5;
  load(fn);
  printf("%g/%g: Done\n\n", _jj, _JMAX);
}


