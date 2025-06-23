//
// hash functions in rlab
//

NITER = 10000;

for (i in 1:NITER)
{
  x = randchar(1,1,8);
  printf("for string='%s', md5 = %s\n", x, hash("md5", x));
}
