//
//
//
if (!exist(NITER))
{ NITER = 100;}

for (i in 1:NITER)
{
  1
  printf("%g: ", i);
  2
  x = randchar(1,1000,2);
  3
  y = randchar(1,1000,2);
  4
  x = set(x);
  y = set(y);
  5
  x = [x,x] + "";
  y = [y,y] + "";
  6
  x = set(x + "");
  y = set(y + "");
  7
  t = intersect(x,y);
//   stop()
  8
  x = union(x,y);
  9
  x = complement(x,y);
  printf("Done\n");
}
