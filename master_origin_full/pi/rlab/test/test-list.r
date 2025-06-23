a.b = <<>>;
a.b.c = [1,2;3,4];
a.b.d = 13;

b = a;

printf ("a refc = %i\t\t b refc = %i\n", entinfo (a).refc, entinfo (b).refc);
printf ("a addr = %s\t b addr = %s\n", entinfo (a).addr, entinfo (b).addr);

# Now, change b.b.d

b.b.d = 42;
printf ("b.b.d = %g\n", b.b.d);
printf ("a.b.d = %g\n", a.b.d);

printf ("a refc = %i\t\t b refc = %i\n", entinfo (a).refc, entinfo (b).refc);
printf ("a addr = %s\t b addr = %s\n", entinfo (a).addr, entinfo (b).addr);

# Now change b.b.c

b.b.c[1;1] = 100;
printf ("a.b.c = \n");
writem ("stdout", a.b.c);
printf ("b.b.c = \n");
writem ("stdout", b.b.c);

printf ("a refc = %i\t\t b refc = %i\n", entinfo (a).refc, entinfo (b).refc);
printf ("a addr = %s\t b addr = %s\n", entinfo (a).addr, entinfo (b).addr);

#
# Get a little fancier...
# Test function argument passing...
#

x = function ( A )
{
  A.b.d = 1.e6;
  printf ("A.b.d = %g\n", A.b.d);
};

printf ("a.b.d = %g\n", a.b.d);
x(a);
printf ("a.b.d = %g\n", a.b.d);

y = function ( A )
{
  b = A;
  b.x = <<>>;
  b.x.y = rand(4,4);
  return b;
};

z = y(a);
if (members(a).n == members(z).n) { error (); }

#
# Now check out clear() usage...
#

clear (b.b.c);
if (members(a).n != members(b).n) { error (); }

