//
// file: eg3_http.r
//
//  downloading page
//  directly, or through proxy

ftp="http://ftp.gnu.org/CRYPTO.README";

if (1)
{
  open(ftp);
  readm(ftp, "crypto.readme");
  close (ftp);
else

  // if you have http proxy available this might work
  curl_opt = <<>>;
  curl_opt.CURLOPT_PROXY     = "your_http_proxy_here";
  curl_opt.CURLOPT_PROXYPORT = port_of_your_http_proxy_here;
  curl_opt.CURLOPT_PROXYTYPE = "http";

  open(ftp, curl_opt);
  readm(ftp, "crypto.readme");
  close (ftp);
}


