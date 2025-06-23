//
//
//

//fns = ls("./eg_eig*.r");
//_NN = int(1000*uniform()) + 1;
_NN = 1000;
fns =["eg_dge3.r"];
for (_jj in 1:_NN)
{
  spinner();
  fn = shuffle(fns,1);
  printf("%g: executing script %s\n", _jj, fn);
  //NITER = int(10*uniform()) + 5;
  NITER = 10;
  load(fn);
  printf("%g: Done\n\n", _jj);
}


