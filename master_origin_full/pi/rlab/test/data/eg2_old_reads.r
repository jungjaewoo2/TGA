//
//
//

fn = "./mtrain09+12+15.txt";

printf("Open file %s in editor, and press 'enter' to continue", fn);
pause("");

printf("We read the content of the file two lines at a time and\n");
printf("stop when we read an empty string matrix [].\n");
printf("We process the data and put each class into its own matrix:\n");
open(fn, "r");

data = <<>>;

while(1)
{
  x = reads(fn,[1,2]);
  if (isempty(x))
  { break; }

//   x

  // label
  l = strsplt(x[2], " ");
  l = tolower(l[2]);
  // data
  d = strsplt(x[1], " ");
  d = strtod(d[2:len(d)]);

  if (!exist(data.[l]))
  {
    data.[l] = d;
  else
    data.[l] = [data.[l]; d];
  }
}
close(fn);

"data = "
data
