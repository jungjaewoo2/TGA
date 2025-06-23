// printf("Execute this script a few times. Eventually it will cause\n");
// printf("segmentation fault within the RLaB2: The reason of crash is\n");
// printf("that the interpreter cannot handle double error in function\n");
// printf("definition: first error occurs in evaluation of what is\n");
// printf("going to be returned, here 'a+y', where 'y' is undefined;\n");
// printf("second, there is nothing to return.\n");
printf("Interpreter has problem how to determine line number under multiple simultaneous error conditions.\n");
printf("To fix this: write functions which do not evaluate final expression on\n");
printf("'return' line. Rather, assign it to a local variable which is then returned.\n");
printf("Here that works as following:\n\n");
printf("\tBAD:\n\n");
printf("\tthefunc = function(V,p){\n");
printf("\t    a = rand(3,3);\n");
printf("\t    return a+y;\n");
printf("\t    };\n\n\n");
printf("\tGOOD:\n\n");
printf("\tthefunc = function(V,p){\n");
printf("\t    a = rand(3,3);\n");
printf("\t    x = a+y;\n");
printf("\t    return x;\n");
printf("\t    };\n\n");

if (!exist(wtf))
{ wtf="b"; }
wtf=prompt("What to do (g-ood, b-ad)", wtf, ["g", "b"]);

thefunc = function(V,p)
{
  global(wtf);

  a=rand(3,3);
  if (wtf=="g")
  {
    printf ("Error created before 'return' statement, so its line number can be determined\n");
    x = a+y;
    return x;
  }
  printf ("Error created during 'return' statement!\n");
  printf ("The line number at which error occurs cannot be determined!\n");

  return (a+y);
};

a = 0.01
b = 0.05;
p = [a,b];


c  = thefunc(loV,p);
