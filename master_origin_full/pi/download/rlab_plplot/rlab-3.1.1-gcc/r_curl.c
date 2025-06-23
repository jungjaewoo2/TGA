//
//
//
#ifdef HAVE_LIBCURL


//
// process configuration options provided to curl, e.g., through
//  open("http://....", url_opts)
//
static int rfile_curl_config (Rfile * rf, char *url, Ent * btree_params)
{
  char *proxy=0;
  size_t iport;
  int istatus=-4, itunnel;
  ListNode * node=0;

  if (!rf)
    return (-3);

  if (!rf->curl)
    return (-2);

  if (isvalidstring(url)<1)
    return (-1);

  // just open the URL
  istatus = curl_easy_setopt(rf->curl, CURLOPT_URL, url);

  if (!btree_params)
    return istatus;

  if (ent_type (btree_params) != BTREE)
    return istatus;

  // CURLOPT_RLABPLUS_DEBUG
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_RLABPLUS_DEBUG",CURLOPT_RLABPLUS_DEBUG);

  // CURLOPT_PROXY
  proxy = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_PROXY",proxy,1,NULL);
  if (proxy)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXY, proxy);

  // CURLOPT_PROXYPORT
  iport=0;
  RLABCODE_PROCESS_BTREE_ENTRY_D(btree_params,node,"CURLOPT_PROXYPORT",iport,class_int,>,0);
  if (iport > 0)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYPORT, iport);

  // CURLOPT_PROXYTYPE
  char * proxytype=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_PROXYTYPE",proxytype,1,NULL);
  if (proxytype)
  {
    if (!strcmp(proxytype, "http"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_HTTP);
    else if (!strcmp(proxytype, "socks4a"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS4A);
    else if (!strcmp(proxytype, "socks4"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS4);
    else if (!strcmp(proxytype, "socks5"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS5);
    else if (!strcmp(proxytype, "socks5h"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS5_HOSTNAME);
    else
    {
      printf("Warning: curl supports 'http', 'socks4a' , 'socks4' , 'socks5' , and 'socks5h' as proxy types!\n");
    }
  }

  // CURLOPT_HTTPPROXYTUNNEL
  itunnel=0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_HTTPPROXYTUNNEL",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPPROXYTUNNEL, 1L);

  // CURLOPT_INTERFACE
  char *ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_INTERFACE",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_INTERFACE, ethN);

  // CURLOPT_PORT
  iport = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_D(btree_params,node,"CURLOPT_PORT",iport,class_int,>,0);
  if (iport > 0)
  {
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_PORT, iport);
  }

  // CURLOPT_FOLLOWLOCATION
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_FOLLOWLOCATION",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_FOLLOWLOCATION, 1L);

  // CURLOPT_UNRESTRICTED_AUTH
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_UNRESTRICTED_AUTH",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_UNRESTRICTED_AUTH, 1L);

  // CURLOPT_MAXREDIRS
  iport = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_D(btree_params,node,"CURLOPT_MAXREDIRS",iport,class_int,>,0);
  if (iport)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_MAXREDIRS, iport);

  // CURLOPT_TCP_NODELAY
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_TCP_NODELAY",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_TCP_NODELAY, 1L);

  // CURLOPT_USERAGENT
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_USERAGENT",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_USERAGENT, ethN);

  // CURLOPT_REFERER
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_REFERER",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_REFERER, ethN);

  // CURLOPT_COOKIE
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_COOKIE",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIE, ethN);

  // CURLOPT_COOKIEFILE
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_COOKIEFILE",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIEFILE, ethN);

  // CURLOPT_COOKIEJAR
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_COOKIEJAR",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIEJAR, ethN);

  // CURLOPT_COOKIESESSION
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_COOKIESESSION",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_COOKIESESSION, 1L);

  // CURLOPT_CRLF
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_CRLF",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_CRLF, 1L);
  else
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_CRLF, 0L);

  //
  // PROTOCOL OPTIONS
  //
  // CURLOPT_POSTQUOTE
  MDS *pq=0;
  RLABCODE_PROCESS_BTREE_ENTRY_MD(btree_params,node,"CURLOPT_POSTQUOTE",pq,MDS,SIZE,1,NULL);
  if (pq)
  {
    int i;
    struct curl_slist *pqlist=0;
    for (i=0; i<SIZE(pq); i++)
      curl_slist_append(pqlist, MdsV0(pq,i));

    istatus |= curl_easy_setopt(rf->curl, CURLOPT_POSTQUOTE, pqlist);

    // free the curl list
    curl_slist_free_all (pqlist);
  }

  // CURLOPT_HTTPGET
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_HTTPGET",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPGET, 1L);
  else
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPGET, 0L);

  //
  // SSL OPTIONS
  //

  // CURLOPT_SSL_VERIFYPEER
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_SSL_VERIFYPEER",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYPEER, 1L);
  else
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYPEER, 0L);

  // CURLOPT_SSL_VERIFYHOST
  itunnel = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_BOOL(btree_params,node,"CURLOPT_SSL_VERIFYPEER",itunnel);
  if (itunnel)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYHOST, 1L);
  else
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_SSL_VERIFYHOST, 0L);

  //
  // NAMES and PASSWORDS OPTIONS (Authentication)
  //

  // CURLOPT_USERPWD
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_USERPWD",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_USERPWD, ethN);

  // CURLOPT_PROXYUSERPWD
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_PROXYUSERPWD",ethN,1,NULL);
  if (ethN)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_PROXYUSERPWD, ethN);

  // CURLOPT_HTTPAUTH
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_HTTPAUTH",ethN,1,NULL);
  if (ethN)
  {
    if (!strcmp(ethN, "basic"))
      curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
    else if (!strcmp(ethN, "digest"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_DIGEST);
    else if (!strcmp(ethN, "gss"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_GSSNEGOTIATE);
    else if (!strcmp(ethN, "ntlm"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_NTLM);
    else if (!strcmp(ethN, "anysafe"))
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_ANYSAFE);
    else
      istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPAUTH, CURLAUTH_ANY);
  }

  // CURLOPT_VERBOSE
  iport = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_D(btree_params,node,"CURLOPT_VERBOSE",iport,class_int,>,0);
  if (iport > 0)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_VERBOSE, iport);
  else
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_VERBOSE, 0);

  // CURLOPT_STDERR
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_STDERR",ethN,1,NULL);
  if (ethN)
  {
    if(rf->curl_stderr_fd)
      fclose(rf->curl_stderr_fd);
    rf->curl_stderr_fd = fopen(ethN, "a");
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_STDERR, rf->curl_stderr_fd);
  }

  //  CURLOPT_TIMEOUT
  iport = 0;
  RLABCODE_PROCESS_BTREE_ENTRY_D(btree_params,node,"CURLOPT_TIMEOUT",iport,class_int,>,0);
  if (iport>0)
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_TIMEOUT, iport);

  // writem: CURLOPT_READDATA - from file
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_READDATA",ethN,1,NULL);
  if (ethN)
  {
    if(rf->curl_read_fd)
      fclose(rf->curl_read_fd);
    rf->curl_read_fd = fopen(ethN, "rb");
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_READFUNCTION, NULL);
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_READDATA, rf->curl_read_fd);
  }

  // readm: CURLOPT_WRITEDATA - to file
  ethN=0;
  RLABCODE_PROCESS_BTREE_ENTRY_S(btree_params,node,"CURLOPT_WRITEDATA",ethN,1,NULL);
  if (ethN)
  {
    if(rf->curl_write_fd)
      fclose(rf->curl_write_fd);
    rf->curl_write_fd = fopen(ethN, "wb");
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_WRITEFUNCTION, NULL);
    istatus |= curl_easy_setopt(rf->curl, CURLOPT_WRITEDATA, rf->curl_write_fd);
  }

  // process HEADER directive
  MDS *header=0;
  RLABCODE_PROCESS_BTREE_ENTRY_MD(btree_params,node,"CURLOPT_HTTPHEADER",header,MDS,SIZE,1,NULL);
  if (header)
  {
    int i;
    struct curl_slist *entry=0;
    for (i=0; i<SIZE(header); i++)
      curl_slist_append(entry, MdsV0(header,i));

    istatus |= curl_easy_setopt(rf->curl, CURLOPT_HTTPHEADER, entry);

    // free the curl list
    curl_slist_free_all (entry);
  }

  return istatus;
}


#undef  THIS_SOLVER
#define THIS_SOLVER "curl.post"
Ent * ent_curl_post (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  char *post_field=0;
  char *fname=0;
  Rfile *rf=0;
  MDS *s=0;
  int i;

  // Check n_args
  if (nargs != 2)
    rerror (THIS_SOLVER ": two arguments required!");

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": first argument must be a single string");
  fname = class_char_pointer( e1);
  if (isvalidstring(fname)<1)
    rerror (THIS_SOLVER ": empty filename");

  // Check if file exists in our list.
  // If the file does not exist pursue default behavior:
  //    User is trying to read an existing file on local file system
  //    However, after reading it, remove it from our list of open files.
  rf = rfile_find(fname);
  if (!rf)
    rerror (THIS_SOLVER ": cannot open file for reading");

  if (rf->filetype != RFILE_CURL)
    rerror (THIS_SOLVER ": curl cannot handle provided URL\n");

  //
  // now check the second argument
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": second argument must be a string matrix");
  s = ent_data (e2);
  if (SIZE(s)<1)
    rerror (THIS_SOLVER ": second argument must be a string matrix");

  curl_easy_setopt(rf->curl, CURLOPT_POST, 1L);
  for (i=0; i<SIZE(s); i++)
  {
    post_field = MdsV0(s,i);
    if (isvalidstring(post_field)<1)
      continue;
    curl_easy_setopt(rf->curl, CURLOPT_POSTFIELDS, post_field);
  }
  CURLcode res = curl_easy_perform (rf->curl);

  // rlab stuff:
  ent_Clean (e1);
  ent_Clean (e2);

  if(res != CURLE_OK)
  {
    fprintf(stderr, THIS_SOLVER " failed with message %s\n", curl_easy_strerror(res));
    return ent_Create_Rlab_Failure();
  }
  return ent_Create_Rlab_Success();
}


#undef  THIS_SOLVER
#define THIS_SOLVER "curl.get.redirect_url"
Ent * ent_curl_get_redirect_url (int nargs, Datum args[])
{
  Ent *e1=0;
  char *location=0;
  char *fname=0;
  Rfile *rf=0;

  // Check n_args
  if (nargs != 1)
    rerror (THIS_SOLVER ": one argument required!");

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": first argument must be a single string");
  fname = class_char_pointer( e1);
  if (isvalidstring(fname)<1)
    rerror (THIS_SOLVER ": empty filename");

  // Check if file exists in our list.
  // If the file does not exist pursue default behavior:
  //    User is trying to read an existing file on local file system
  //    However, after reading it, remove it from our list of open files.
  rf = rfile_find(fname);
  if (!rf)
    rerror (THIS_SOLVER ": cannot open file for reading");

  if (rf->filetype != RFILE_CURL)
    rerror (THIS_SOLVER ": curl cannot handle provided URL\n");

  CURLcode res = curl_easy_getinfo(rf->curl, CURLINFO_REDIRECT_URL, &location);

  if(res != CURLE_OK)
  {
    location = 0;
    fprintf(stderr, THIS_SOLVER " failed with message %s\n", curl_easy_strerror(res));
  }

  // rlab stuff:
  ent_Clean (e1);
  return ent_Create_Rlab_String(location);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "curl.get.response_code"
Ent * ent_curl_get_response_code (int nargs, Datum args[])
{
  Ent *e1=0;
  size_t response_code;
  char *fname=0;
  Rfile *rf=0;

  // Check n_args
  if (nargs != 1)
    rerror (THIS_SOLVER ": one argument required!");

  //
  // get the name of the stream
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": first argument must be a single string");
  fname = class_char_pointer( e1);
  if (isvalidstring(fname)<1)
    rerror (THIS_SOLVER ": empty filename");

  // Check if file exists in our list.
  // If the file does not exist pursue default behavior:
  //    User is trying to read an existing file on local file system
  //    However, after reading it, remove it from our list of open files.
  rf = rfile_find(fname);
  if (!rf)
    rerror (THIS_SOLVER ": cannot open file for reading");

  if (rf->filetype != RFILE_CURL)
    rerror (THIS_SOLVER ": curl cannot handle provided URL\n");

  CURLcode res = curl_easy_getinfo(rf->curl, CURLINFO_RESPONSE_CODE, &response_code);

  if(res != CURLE_OK)
  {
    response_code = RLAB_STATUS_FAILURE;
    fprintf(stderr, THIS_SOLVER " failed with message %s\n", curl_easy_strerror(res));
  }
  else
    response_code = RLAB_STATUS_SUCCESS;

  // rlab stuff:
  ent_Clean (e1);
  return ent_Create_Rlab_Double(response_code);
}


#endif