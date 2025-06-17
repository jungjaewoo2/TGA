//
// file: eg_md1.r
//
// test of a matrix dense integer data type

NITER = 1000;

for(i in 1:NITER)
{
  spinner();

  n = int(10*rand() + 1)+0;

  if(i == NITER){ n = 5;}

  x = int  ( (rand(n,n)-0.5) * 15 );
  x = round( int( (rand(n,n)-0.5) * 15 ) );
  x = abs  ( int( (rand(n,n)-0.5) * 15 ) );
}

printf("[x,x] =\n");
[x,x]

printf("[x;x] =\n");
[x;x]

printf("2*x =\n");
x = int(2).*x

printf("change x[;1] in x =\n");
x[;1] = int( [1;2;3;4;5] )

printf("change x[1;] in x =\n");
x[1;] = int( [1,2,3,4,5] )

zero = zeros(x);
izero = int( zero );

printf("for-loop in x[;1] =\n");
for( i in x[;1])
{
  printf("(%s) i = %i\n", i.type, i);
}

printf("x==int(0) && x==0 : ");
y1 = [x==izero];
y2 = [x==zero];
printf("%g\n", all(all(y1==y2)) );

printf("x!=int(0) && x!=0 : ");
y1 = [x!=izero];
y2 = [x!=zero];
printf("%g\n", all(all(y1==y2)) );

printf("x>int(0) && x>0 : ");
y1 = [ x > izero ];
y2 = [ x > zero ];
printf("%g\n", all(all(y1==y2)) );

printf("x>=int(0) && x>=0 : ");
y1 = [ x >= izero ];
y2 = [ x >= zero ];
printf("%g\n", all(all(y1==y2)) );

printf("x<int(0) && x<0 : ");
y1 = [ x < izero ];
y2 = [ x < zero ];
printf("%g\n", all(all(y1==y2)) );

printf("x<=int(0) && x<=0 :");
y1 = [ x <= izero ];
y2 = [ x <= zero ];
printf("%g\n", all(all(y1==y2)) );

printf("!x\n");
!x

printf("x\n");
x

printf("max(x[;1]) = %g\n", max(x[;1]));
printf("max(x[1;]) = %g\n", max(x[1;]));
max(x);
max(x,x);
max(x,x+0);

printf("min(x[;1]) = %g\n", min(x[;1]));
printf("min(x[1;]) = %g\n", min(x[1;]));
min(x);
min(x,x);
min(x,x+0);

printf("all zeros should follow:\n");
printf("sum:\n");
sum(x)-sum(x+0)
printf("cumsum:\n");
cumsum(x)-cumsum(x+0)
printf("prod:\n");
prod(x) - prod(x+0)
printf("cumprod:\n");
cumprod(x) - cumprod(x+0)
printf("isinf:\n");
isinf(x)
printf("isnan: \n");
isnan(x)
printf("finite: \n");
finite(x)

printf("sin : %g\n", max(max(sin(x)-sin(x+0))) );
printf("cos : %g\n", max(max(cos(x)-cos(x+0))) );
printf("tan : %g\n", max(max(tan(x)-tan(x+0))) );
printf("asin: %g\n", max(max(asin(x)-asin(x+0))) );
//printf("acos: %g\n", max(max(acos(x)-acos(x+0))) );
printf("atan: %g\n", max(max(atan(x)-atan(x+0))) );

diag(x);
diag(x+0);

save();

clearall();
read("rlab.save");
whos();

writem("test.txt",x);
close("test.txt")


