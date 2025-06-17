//
// hash functions in rlab
//

NITER = 1000;

if(!isfile("test.txt"))
{
  stop("for this example to work there has to exist file 'test.txt'");
}

for (i in 1:NITER)
{
  // read the content of the file 'test.txt' into
  // a string variable 'x'
  x = reads("test.txt");

  // hash it with 'md5'
  y = hash("md5", x);
}

printf("x = \n");
x
printf("md5 = %s\n", y);

// now use external program 'md5sum' to find
// the md5 hash of the file that was read in
// into 'x'
printf("\nSame using the external program 'md5sum' on file 'test.txt':\n");
system("md5sum test.txt");

