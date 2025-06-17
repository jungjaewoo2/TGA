//
// file: eq1_http.r
//
// download stock quotes from yahoo.com in two different ways



url1="http://finance.yahoo.com/q/hp?s=PFE&a=00&b=4&c=1982&d=00&e=10&f=2010&g=d";
url2="http://ichart.finance.yahoo.com/table.csv?s=PFE&a=00&b=4&c=1982&d=00&e=10&f=2010&g=d&ignore=.csv";

open(url1);
x1 = readm(url1);

open(url2);
x2 = readm(url2);

x1
size(x1)

colors("red");
printf("One could try to process the previous html dump,\n");
// pause()
colors()

printf("Bur one could also try to get the information from the site directly in csv format:\n");
// pause()
x2
size(x2)
