//
// file: eg1_http.r
//
//  downloading a http page

http="http://ichart.finance.yahoo.com/table.csv?s=PFE&d=0&e=12&f=2010&g=d&a=0&b=4&c=1982&ignore=.csv";

// put your proxy here:
opt = <<>>;
// opt.CURLOPT_PROXY     = "your.proxy.here";
// opt.CURLOPT_PROXYPORT = PROXYPORT;
// opt.CURLOPT_PROXYTYPE = "http";
open(http, opt);
x = readm(http);
readm(http, "./yahoo_quote.csv");
close (http);

size(x)
y=strsplt(x, char(10))';

size(y)

y[1]
y[2]
lastr(y)

