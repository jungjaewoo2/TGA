#ifdef HAVE_LIBCURL
extern Ent * ent_curl_get_response_code (int nargs, Datum args[]);
extern Ent * ent_curl_get_redirect_url  (int nargs, Datum args[]);
extern Ent * ent_curl_post (int nargs, Datum args[]);

//
// curl 
//
Bltin rlab_curl_get_bltin[] = {
  {BLTIN, "response_code", ent_curl_get_response_code},
  {BLTIN, "redirect_url", ent_curl_get_redirect_url},
  {0, 0, 0},
};

Bltin rlab_curl_bltin[] = {
  {BLTIN, "post", ent_curl_post},
  {0, 0, 0},
};

#endif