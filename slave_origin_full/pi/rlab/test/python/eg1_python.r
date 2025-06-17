//
//
//
py.init();


// scalars
x1 = rand();
py.var(,"x", x1);
colors("red");
py.cmd("print 'x=',x");
colors();
y1 = py.var(,"x");
d1 = max(max(abs(x1-y1)));
if (d1 != 0)
{
  error("test failed!");
else
  printf("Test passed for real scalars\n");
}

// vectors
x2 = rand(1,6);
py.var(,"x",x2);
colors("red");
py.cmd("print 'x=',x");
colors();
y2 = py.var(,"x");
d2 = max(max(abs(x2-y2)));
if (d2 != 0)
{
  error("test failed!");
else
  printf("Test passed for real vectors\n");
}

// matrices
x3 = rand(4,6);
py.var(,"x",x3);
colors("red");
py.cmd("print 'x=',x");
colors();
y3 = py.var(,"x");
d3 = max(max(abs(x3-y3)));
if (d3 != 0)
{
  error("test failed!");
else
  printf("Test passed for real matrices\n");
}

printf("\n\nExpressions:\n");
//
py.cmd([...
  "x=1.0", "y=2.0" ...
]);
py.cmd("print '(python): x+y=', x+y");
c=py.eval("x+y");
printf("(rlab)  : x+y= %g\n", c);





