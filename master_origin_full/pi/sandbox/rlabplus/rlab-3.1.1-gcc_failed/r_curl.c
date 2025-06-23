//
//
//
#ifdef HAVE_LIBCURL


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
  long response_code;
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
    response_code = -1;
    fprintf(stderr, THIS_SOLVER " failed with message %s\n", curl_easy_strerror(res));
  }

  // rlab stuff:
  ent_Clean (e1);
  return ent_Create_Rlab_Double(response_code);
}


#endif