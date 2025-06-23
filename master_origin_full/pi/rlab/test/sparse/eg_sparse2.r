#
# Test the Rlab sparse matrix factorization functions.
# Assuming, this version of Rlab has been compiled to use one of the
# sparse matrix factorization libraries, like SuperLU.
#

rfile speye sprand
rfile randsvd bandred

n = 50;                              # Size of problems.
X = 4;                               # Hueristic.


sa = randsvd (n,10);
sa = sparse (bandred (sa,8));
b = rand (n, 2);

zsa = randsvd (n,10) + randsvd (n,10)*1j;
zsa = sparse (bandred (zsa, 8));
zb = rand(n,2) + rand(n,2)*1j;

#
# Test the "standard" solve...
#


x1 = solve (sa,b);
max (abs (sa*x1-b));

if (max (max (abs (sa*x1-b))/(X*mnorm(sa,"1")*sa.nr)) > eps){
    error ("superlu problem");
else printf ("\tpassed sparse-real solve() test ...\n");}

z1 = solve (zsa,zb);
max (abs (zsa*z1-zb));

if (max (max (abs (zsa*z1-zb))/(X*mnorm(zsa,"1")*zsa.nr)) > eps){
    error ("superlu problem");
else printf ("\tpassed sparse-complex solve test ...\n");}

#
# Test factor() / backsub() pair (equiv to solve)...
#

fa = factor (sa);
x2 = backsub (fa,b);
max (abs (sa*x2-b));
if (max (max (abs (sa*x2-b))/(X*mnorm(sa,"1")*sa.nr)) > eps)
{
  error ("superlu factor/backsub problem");
  else
  printf ("\tpassed sparse-real factor/backsub test ...\n");
}

fa = spfactor (sa);
x2 = backsub (fa,b);
max (abs (sa*x2-b));
if (max (max (abs (sa*x2-b))/(X*mnorm(sa,"1")*sa.nr)) > eps)
{
  error ("superlu spfactor/backsub 1 problem");
  else
  printf ("\tpassed sparse-real spfactor/backsub test 1...\n");
}

fa = spfactor (sa);
x2 = backsub (fa,b);
max (abs (sa*x2-b));
if (max (max (abs (sa*x2-b))/(X*mnorm(sa,"1")*sa.nr)) > eps)
{
  error ("superlu spfactor/backsub 2 problem");
  else
  printf ("\tpassed sparse-real spfactor/backsub test 2...\n");
}

fz = factor (zsa);
zx2 = backsub (fz,zb);
max (abs (zsa*zx2-zb));
if (max (max (abs (zsa*zx2-zb))/(X*mnorm(zsa,"1")*zsa.nr)) > eps)
{
  error ("superlu factor/backsub problem");
  else
  printf ("\tpassed sparse-complex factor/backsub test ...\n");
}

fz = spfactor (zsa);
zx2 = backsub (fz,zb);
max (abs (zsa*zx2-zb));
if (max (max (abs (zsa*zx2-zb))/(X*mnorm(zsa,"1")*zsa.nr)) > eps)
{
  error ("superlu spfactor/backsub 1 problem");
  else
  printf ("\tpassed sparse-complex spfactor/backsub test 1...\n");
}

fz = spfactor (zsa);
zx2 = backsub (fz,zb);
max (abs (zsa*zx2-zb));
if (max (max (abs (zsa*zx2-zb))/(X*mnorm(zsa,"1")*zsa.nr)) > eps)
{
  error ("superlu spfactor/backsub 2 problem");
  else
  printf ("\tpassed sparse-complex spfactor/backsub test 2...\n");
}

