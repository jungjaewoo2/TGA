#define RLAB_DATE_IDX_YEAR 0
#define RLAB_DATE_IDX_MONT 1
#define RLAB_DATE_IDX_DAYS 2
#define RLAB_DATE_IDX_HOUR 3
#define RLAB_DATE_IDX_MINS 4
#define RLAB_DATE_IDX_SECS 5
#define RLAB_DATE_IDX_IDST 6
#define RLAB_DATE_IDX_GMTO 7

//
// time: what is the current time in 1-time format?
//
// Ent *
// ent_time (int nargs, Datum args[])
// {
//   return ent_Create_Rlab_Double( (double) time(NULL) );
// }

//
// localtime: any time -> 8-time
//
#undef THIS_SOLVER
#define THIS_SOLVER "clock"
Ent *
ent_Clock (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, nx=1;

  MDR *x=0, *w=0;
  time_t tt1;
  struct tm *time_str, time_str2;
  struct timezone tzp;
  struct timeval tv;

  if (nargs != 1 && nargs!=0)
  {
    fprintf(stdout, THIS_SOLVER ": Converts time from epoch seconds to 8-column format.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   tm = " THIS_SOLVER "(/t/)\n");
    rerror (RLAB_ERROR_ARG1_MDR_VECTOR);
  }

  if (nargs == 1)
  {
    // get time
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

    x  = ent_data (e1);
    if (MNC(x)!=1 && MNC(x)!=6 && MNC(x)!=8)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
    nx = MNR(x);

    if (MNC(x) != 8)
      gettimeofday (&tv, &tzp);
  }

  w  = mdi_Create (nx,8);

  for (i=0; i<nx; i++)
  {
    if (x)
    {
      if (MNC(x)==1)
        tt1 = (time_t) mdiV0(x, i);
      else if (MNC(x)==8 || MNC(x)==6)
      {
        time_str2.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
        time_str2.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
        time_str2.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
        time_str2.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
        time_str2.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
        time_str2.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
        time_str2.tm_wday  =  -1;
        time_str2.tm_yday  =  -1;
        if (MNC(x) == 8)
        {
          time_str2.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
          time_str2.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
        }
        else
        {
          time_str2.tm_isdst =  -1;  // i don't know
          time_str2.tm_gmtoff = 60 * tzp.tz_minuteswest;
        }
        tt1 = mktime( &time_str2 );
      }
    }
    else
    {
      tt1 = time (NULL);
    }

    time_str = localtime(&tt1);

    Mdi0(w,i,RLAB_DATE_IDX_YEAR) = time_str->tm_year + 1900;
    Mdi0(w,i,RLAB_DATE_IDX_MONT) = time_str->tm_mon + 1;
    Mdi0(w,i,RLAB_DATE_IDX_DAYS) = time_str->tm_mday;
    Mdi0(w,i,RLAB_DATE_IDX_HOUR) = time_str->tm_hour;
    Mdi0(w,i,RLAB_DATE_IDX_MINS) = time_str->tm_min;
    Mdi0(w,i,RLAB_DATE_IDX_SECS) = time_str->tm_sec;
    Mdi0(w,i,RLAB_DATE_IDX_IDST) = time_str->tm_isdst;
    Mdi0(w,i,RLAB_DATE_IDX_GMTO) = time_str->tm_gmtoff / 3600;
  }

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR (w);
}

//
// mktime: 6,8-time to 1-time
//
#undef THIS_SOLVER
#define THIS_SOLVER "seconds"
Ent *
ent_Seconds (int nargs, Datum args[])
{
  Ent *e1=0;

  MDR *x=0, *r=0;
  int nx=1, i;

  struct tm time_str;
  struct timezone tzp;
  struct timeval tv;

  if (nargs != 0 && nargs != 1)
  {
    fprintf(stdout, THIS_SOLVER ": Converts time from the 8-entry format to second count.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   t = mktime(tm)\n");
    fprintf(stdout, THIS_SOLVER ": where 'tm' is a 10-column matrix containing the time specification,\n");
    fprintf(stdout, THIS_SOLVER ": while 't' is the time as a number of seconds since unix epoch.\n");
    rerror ("one or no arguments required");
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": improper first argument 'time'\n");
    x  = ent_data (e1);
    if (MNC(x) != 8 && MNC(x) != 6 && MNC(x)!=1)
      rerror (THIS_SOLVER ": first argument 'time' does not have 10 columns\n");

    nx = MNR(x);
  }

  r  = mdr_Create(nx, 1);

  if (x)
  {
    // if x is single column, assume seconds, and just return its copy
    if (MNC(x) ==1)
    {
      ent_IncRef(e1);
//       return ent_Assign_Rlab_MDR (mdr_Copy(x));
      return e1;
    }

    if (MNC(x) != 8)
      gettimeofday (&tv, &tzp);
  }

  for (i=0; i<nx; i++)
  {
    if (x)
    {
      time_str.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
      time_str.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
      time_str.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
      time_str.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
      time_str.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
      time_str.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
      time_str.tm_wday  =  -1;
      time_str.tm_yday  =  -1;
      if (MNC(x) == 8)
      {
        time_str.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
        time_str.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
      }
      else
      {
        time_str.tm_isdst =  -1;  // i don't know
        time_str.tm_gmtoff = 60 * tzp.tz_minuteswest;
      }

      MdrV0 (r, i) = (double) mktime( &time_str );
    }
    else
      MdrV0 (r, i) = (double) time(NULL);
  }

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR (r);
}



Ent *
ent_gmtime (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *w=0, *x=0;
  int i,nx;

  time_t tt1;

  struct tm *time_str, time_str2;
  struct timezone tzp;
  struct timeval tv;

  if (nargs!=0 && nargs!=1)
  {
    fprintf(stdout, "gmtime: Returns current UTC time in 8-col format.\n");
    fprintf(stdout, "localtime:   tm = gmtime()\n");
    rerror ("no arguments required");
  }
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("gmtime: improper first argument 'time'\n");
    x  = ent_data (e1);
    if (MNC(x)!=8 && MNC(x)!=6)
      rerror ("gmtime: first argument 'time' must have 8 or 6 columns!\n");
  }

  if (!x)
  {
    w = mdi_Create (1, 8);

    time(&tt1);
    time_str = gmtime ( &tt1 );

    MdiV0(w,RLAB_DATE_IDX_SECS) = time_str->tm_sec;
    MdiV0(w,RLAB_DATE_IDX_MINS) = time_str->tm_min;
    MdiV0(w,RLAB_DATE_IDX_HOUR) = time_str->tm_hour;
    MdiV0(w,RLAB_DATE_IDX_DAYS) = time_str->tm_mday;
    MdiV0(w,RLAB_DATE_IDX_MONT) = time_str->tm_mon + 1;
    MdiV0(w,RLAB_DATE_IDX_YEAR) = time_str->tm_year + 1900;
    MdiV0(w,RLAB_DATE_IDX_IDST) = time_str->tm_isdst;
    MdiV0(w,RLAB_DATE_IDX_GMTO) = 0;
 }
  else
  {
    nx = MNR(x);
    w  = mdi_Create (nx, 8);

    if (MNC(x) != 8)
    {
      // user did not provide DST or TZ information: read it from the system
      // and assume it applies to dates in question
      gettimeofday (&tv, &tzp);
    }

    for (i=0; i<nx; i++)
    {
      if (MNC(x)==6 || MNC(x)==8)
      {
        time_str2.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
        time_str2.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
        time_str2.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
        time_str2.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
        time_str2.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
        time_str2.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
        time_str2.tm_wday  =  -1;
        time_str2.tm_yday  =  -1;
        if (MNC(x) == 8)
        {
          time_str2.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
          time_str2.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
        }
        else
        {
          time_str2.tm_isdst =  -1;  // i don't know
          time_str2.tm_gmtoff = 60 * tzp.tz_minuteswest;
        }

        tt1 = mktime( &time_str2 );
      }
      else
        tt1 = mdiV0(x,i);

      time_str = gmtime(&tt1);

      Mdi0(w,i,RLAB_DATE_IDX_SECS) = time_str->tm_sec;
      Mdi0(w,i,RLAB_DATE_IDX_MINS) = time_str->tm_min;
      Mdi0(w,i,RLAB_DATE_IDX_HOUR) = time_str->tm_hour;
      Mdi0(w,i,RLAB_DATE_IDX_DAYS) = time_str->tm_mday;
      Mdi0(w,i,RLAB_DATE_IDX_MONT) = time_str->tm_mon + 1;
      Mdi0(w,i,RLAB_DATE_IDX_YEAR) = time_str->tm_year + 1900;
      Mdi0(w,i,RLAB_DATE_IDX_IDST) = time_str->tm_isdst;
      Mdi0(w,i,RLAB_DATE_IDX_GMTO) = 0;
    }
  }

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR (w);
}


//
//
//
#undef THIS_SOLVER
#define THIS_SOLVER "time2dstr"
Ent *
ent_strftime (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;

  char *dummyfmt = "%c (%z/%Z)";
  char s1[100], *fmt=0;

  MDR *x=0;
  MDS *s=0;
  int nx, i;

  struct tm time_str, *time_str2;
  struct timezone tzp;
  struct timeval tv;
  time_t tt1;

  if (nargs != 1 && nargs != 2)
  {
    fprintf(stdout, THIS_SOLVER ": Converts time to a string according to a provided format.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   s = strftime(tm/, fmt/)\n");
    fprintf(stdout, THIS_SOLVER ": where 'tm' is a 6 or 8 column matrix containing the time specification,\n");
    fprintf(stdout, THIS_SOLVER ": while 'fmt' is the time format according to C-specifications.\n");
    rerror ("one or two arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": improper first argument 'time'\n");
  x  = ent_data (e1);
  if (MNC(x)!=8 && MNC(x)!=6 && MNC(x)!=1)
    rerror (THIS_SOLVER ": first argument 'tm' does not have 1,6 or 8 columns\n");

  if (nargs == 1)
  {
    fmt = dummyfmt;
  }
  else if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_STRING)
      rerror (THIS_SOLVER ": improper first argument 'fmt'\n");
    fmt = class_char_pointer(e2);
    if (!fmt)
      fmt = dummyfmt;
  }

  nx = MNR(x);
  s  = mds_Create(nx, 1);

  if (MNC(x)==8 || MNC(x)==6)
  {
    if (MNC(x) != 8)
    {
      // user did not provide DST or TZ information: read it from the system
      // and assume it applies to dates in question
      gettimeofday (&tv, &tzp);
    }

    for (i=0; i<nx; i++)
    {
      time_str.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
      time_str.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
      time_str.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
      time_str.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
      time_str.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
      time_str.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
      time_str.tm_wday  =  -1;
      time_str.tm_yday  =  -1;
      if (MNC(x) == 8)
      {
        time_str.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
        time_str.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
      }
      else
      {
        time_str.tm_isdst =  -1;  // i don't know
        time_str.tm_gmtoff = 60 * tzp.tz_minuteswest;
      }

      // we rewrite time structure, otherwise strftime segfaults!
      tt1 = mktime(&time_str);
      time_str2 = localtime(&tt1);
      strftime (s1, sizeof (s1), fmt, time_str2);
      MdsV0 (s, i) = cpstr (s1);
    }
  }
  else
  {
    for (i=0; i<nx; i++)
    {
      tt1 =   (time_t) mdiV0(x,i);

      time_str2 = localtime(&tt1);
      strftime (s1, sizeof (s1), fmt, time_str2);
      MdsV0 (s, i) = cpstr (s1);
    }
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDS (s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "etime"
Ent *
ent_Etime (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;

  MDR *x1=0, *x2=0, *w=0;
  int nx1, nx2, nx, i, i1=0, i2=0;

  struct tm time_str1, time_str2;
  time_t tt1, tt2;
  struct timezone tzp;
  struct timeval tv;


  if (nargs != 2)
  {
    fprintf(stdout, THIS_SOLVER ": Finds the difference between two times in seconds.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   s = " THIS_SOLVER "(time1, time2)\n");
    fprintf(stdout, THIS_SOLVER ": where 'time1' and 'time2' are in a 9-column matrix containing the\n");
    fprintf(stdout, THIS_SOLVER ": time specification.\n");
    rerror ("two arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": improper first argument 'time1'\n");
  x1  = ent_data (e1);
  if (MNC(x1)!=8 && MNC(x1)!=6 && MNC(x1)!=1)
    rerror (THIS_SOLVER ": first argument 'time1' must have 1,6, or 8 columns\n");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": second argument 'time2' must have 1,6, or 8 columns\n");
  x2  = ent_data (e2);
  if (MNC(x2)!=8 && MNC(x2)!=6 && MNC(x2)!=1)
    rerror (THIS_SOLVER ": second argument 'time2' must have 1,6, or 8 columns\n");

  nx1 = MNR(x1);
  nx2 = MNR(x2);
  nx  = MAX(nx1,nx2);
  w   = mdr_Create(nx, 1);

  if (MNC(x1)==6 || MNC(x2)==6)
  {
    // user did not provide DST or TZ information: read it from the system
    // and assume it applies to dates in question
    gettimeofday (&tv, &tzp);
  }

  for (i=0; i<nx; i++)
  {
    // figure out the first time
    if (MNC(x1)==6 || MNC(x1)==8)
    {
      i1 = MIN(i,nx1-1);
      time_str1.tm_sec   =  mdi0(x1,i1,RLAB_DATE_IDX_SECS);
      time_str1.tm_min   =  mdi0(x1,i1,RLAB_DATE_IDX_MINS);
      time_str1.tm_hour  =  mdi0(x1,i1,RLAB_DATE_IDX_HOUR);
      time_str1.tm_mday  =  mdi0(x1,i1,RLAB_DATE_IDX_DAYS);
      time_str1.tm_mon   =  mdi0(x1,i1,RLAB_DATE_IDX_MONT) - 1;
      time_str1.tm_year  =  mdi0(x1,i1,RLAB_DATE_IDX_YEAR) - 1900;
      time_str1.tm_wday  =  -1;
      time_str1.tm_yday  =  -1;
      if (MNC(x1) == 8)
      {
        time_str1.tm_isdst =  mdi0(x1,i1,RLAB_DATE_IDX_IDST);
        time_str1.tm_gmtoff = 3600 * mdi0(x1,i1,RLAB_DATE_IDX_GMTO);
      }
      else
      {
        time_str1.tm_isdst =  -1;  // i don't know
        time_str1.tm_gmtoff = 60 * tzp.tz_minuteswest;
      }
      tt1 = mktime(&time_str1);
    }
    else
      tt1 = (time_t) mdiV0(x1,i1);

    // figure out the second time
    if (MNC(x2)==6 || MNC(x2)==8)
    {
      i2 = MIN(i,nx2-1);
      time_str2.tm_sec   =  mdi0(x2,i2,RLAB_DATE_IDX_SECS);
      time_str2.tm_min   =  mdi0(x2,i2,RLAB_DATE_IDX_MINS);
      time_str2.tm_hour  =  mdi0(x2,i2,RLAB_DATE_IDX_HOUR);
      time_str2.tm_mday  =  mdi0(x2,i2,RLAB_DATE_IDX_DAYS);
      time_str2.tm_mon   =  mdi0(x2,i2,RLAB_DATE_IDX_MONT) - 1;
      time_str2.tm_year  =  mdi0(x2,i2,RLAB_DATE_IDX_YEAR) - 1900;
      time_str2.tm_wday  =  -1;
      time_str2.tm_yday  =  -1;
      if (MNC(x2) == 8)
      {
        time_str2.tm_isdst =  mdi0(x2,i2,RLAB_DATE_IDX_IDST);
        time_str2.tm_gmtoff = 3600 * mdi0(x2,i2,RLAB_DATE_IDX_GMTO);
      }
      else
      {
        time_str2.tm_isdst =  -1;  // i don't know
        time_str2.tm_gmtoff = 60 * tzp.tz_minuteswest;
      }
      tt2 = mktime(&time_str2);
    }
    else
      tt2 = (time_t) mdiV0(x2,i2);

    MdrV0 (w, i) = difftime(tt1, tt2);
  }

  return ent_Assign_Rlab_MDR (w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "dayofweek"
Ent *
ent_dayofweek (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *w=0, *x=0;
  int i,nx=1;

  time_t tt1;

  struct tm *time_str, time_str2;
  struct timezone tzp;
  struct timeval tv;

  if (nargs!=0 && nargs!=1)
  {
    fprintf(stdout, THIS_SOLVER ": Finds day of the week of the date matrix\n");
    fprintf(stdout, THIS_SOLVER ":   d = dayofweek(t)\n");
    rerror ("no arguments required");
  }
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": improper first argument 'time'\n");
    x  = ent_data (e1);
    if (MNC(x)!=8 && MNC(x)!=6 && MNC(x)!=1)
      rerror (THIS_SOLVER ": first argument 'time' must have 1, 6 or 8 columns!\n");

    nx = MNR(x);
  }

  w = mdr_Create (nx, 1);
  if (!x)
  {
    tt1 = time (NULL);
    time_str = localtime ( &tt1 );
    MdrV0(w,0) = (time_str->tm_wday > 0 ? time_str->tm_wday : 7) - 1
    + (time_str->tm_hour + time_str->tm_min/60.0 + time_str->tm_sec/2600.0)/24.0;
  }
  else
  {
    if (MNC(x) == 6)
    {
      // user did not provide DST or TZ information: read it from the system
      // and assume it applies to dates in question
      gettimeofday (&tv, &tzp);
    }

    for (i=0; i<nx; i++)
    {
      if (MNC(x)==6 || MNC(x)==8)
      {
        time_str2.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
        time_str2.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
        time_str2.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
        time_str2.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
        time_str2.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
        time_str2.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
        time_str2.tm_wday  =  -1;
        time_str2.tm_yday  =  -1;
        if (MNC(x) == 8)
        {
          time_str2.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
          time_str2.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
        }
        else
        {
          time_str2.tm_isdst =  -1;  // i don't know
          time_str2.tm_gmtoff = 60 * tzp.tz_minuteswest;
        }

        tt1 = mktime( &time_str2 );
      }
      else
        tt1 = mdiV0(x,i);

      time_str = localtime ( &tt1 );
      MdrV0(w,i) = (time_str->tm_wday > 0 ? time_str->tm_wday : 7) - 1
                    + (time_str->tm_hour + time_str->tm_min/60.0 + time_str->tm_sec/2600.0)/24.0;
    }
  }

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR (w);
}


#undef THIS_SOLVER
#define THIS_SOLVER "dayofyear"
Ent *
ent_dayofyear (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *w=0, *x=0;
  int i,nx=1;

  time_t tt1;

  struct tm *time_str, time_str2;
  struct timezone tzp;
  struct timeval tv;

  if (nargs!=0 && nargs!=1)
  {
    fprintf(stdout, THIS_SOLVER ": Finds day of the week of the date matrix\n");
    fprintf(stdout, THIS_SOLVER ":   d = dayofyear(t)\n");
    rerror ("no arguments required");
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": improper first argument 'time'\n");
    x  = ent_data (e1);
    if (MNC(x)!=8 && MNC(x)!=6)
      rerror (THIS_SOLVER ": first argument 'time' must have 8 or 6 columns!\n");

    nx = MNR(x);
  }

  w = mdr_Create (nx, 1);
  if (!x)
  {
    tt1 = time (NULL);
    time_str = localtime ( &tt1 );
    MdrV0(w,0) = time_str->tm_yday + (time_str->tm_hour + time_str->tm_min/60.0 + time_str->tm_sec/2600.0)/24.0;
  }
  else
  {
    if (MNC(x) != 8)
    {
      // user did not provide DST or TZ information: read it from the system
      // and assume it applies to dates in question
      gettimeofday (&tv, &tzp);
    }

    for (i=0; i<nx; i++)
    {
      if (MNC(x)==6 || MNC(x)==8)
      {
        time_str2.tm_sec   =  mdi0(x,i,RLAB_DATE_IDX_SECS);
        time_str2.tm_min   =  mdi0(x,i,RLAB_DATE_IDX_MINS);
        time_str2.tm_hour  =  mdi0(x,i,RLAB_DATE_IDX_HOUR);
        time_str2.tm_mday  =  mdi0(x,i,RLAB_DATE_IDX_DAYS);
        time_str2.tm_mon   =  mdi0(x,i,RLAB_DATE_IDX_MONT) - 1;
        time_str2.tm_year  =  mdi0(x,i,RLAB_DATE_IDX_YEAR) - 1900;
        time_str2.tm_wday  =  -1;
        time_str2.tm_yday  =  -1;
        if (MNC(x) == 8)
        {
          time_str2.tm_isdst =  mdi0(x,i,RLAB_DATE_IDX_IDST);
          time_str2.tm_gmtoff = 3600 * mdi0(x,i,RLAB_DATE_IDX_GMTO);
        }
        else
        {
          time_str2.tm_isdst =  -1;  // i don't know
          time_str2.tm_gmtoff = 60 * tzp.tz_minuteswest;
        }

        tt1 = mktime( &time_str2 );
      }
      else
        tt1 = mdiV0(x,i);

      time_str = localtime ( &tt1 );
      MdrV0(w,i) = time_str->tm_yday + (time_str->tm_hour + time_str->tm_min/60.0 + time_str->tm_sec/2600.0)/24.0;
    }
  }

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR (w);
}


//
//
//
#define _XOPEN_SOURCE
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./clibs/gnulib/gnulib.h"

extern char *strptime(const char *s, const char *format, struct tm *tm);

#undef THIS_SOLVER
#define THIS_SOLVER "dstr2time"
Ent *
ent_dstr2time (int nargs, Datum args[])
{
  Ent *e1=0,*e2=0;

  MDR *w=0;
  MDS *s=0;
  int ns, i;
  char *f=0;
  char *fmt=0;

  struct tm time_str;

  if (nargs != 1 && nargs != 2)
  {
    fprintf(stdout, THIS_SOLVER ": Converts POSIX-compliant date string to 8-column format.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   t = " THIS_SOLVER "(s/,fmt/)\n");
    fprintf(stdout, THIS_SOLVER ": where 's' is POSIX date.\n");
    rerror ("one or two arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_VECTOR);
    goto end_of_ent_dstr2time;
  }
  s = ent_data (e1);
  ns = MNR(s) * MNC(s);
  if (!ns)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_VECTOR);
    goto end_of_ent_dstr2time;
  }

  if (nargs>1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      f = class_char_pointer(e2);
      if (f)
        if (strlen(f))
          fmt = f;
    }
  }

  w  = mdi_Create(ns, 8);
  mdr_Zero(w);

  if (fmt)
  {
    for (i=0; i<ns; i++)
    {
      if(strptime(MdsV0(s,i), fmt, &time_str))
      {
        mktime(&time_str);
        Mdi0(w,i,RLAB_DATE_IDX_SECS) = time_str.tm_sec;
        Mdi0(w,i,RLAB_DATE_IDX_MINS) = time_str.tm_min;
        Mdi0(w,i,RLAB_DATE_IDX_HOUR) = time_str.tm_hour;
        Mdi0(w,i,RLAB_DATE_IDX_DAYS) = time_str.tm_mday;
        Mdi0(w,i,RLAB_DATE_IDX_MONT) = time_str.tm_mon + 1;
        Mdi0(w,i,RLAB_DATE_IDX_YEAR) = time_str.tm_year + 1900;
        Mdi0(w,i,RLAB_DATE_IDX_IDST) = time_str.tm_isdst;
        Mdi0(w,i,RLAB_DATE_IDX_GMTO) = time_str.tm_gmtoff/3600;
      }
    }
  }
  else
  {
    // no format: use gnu's parse_datetime
    for (i=0; i<ns; i++)
    {
      if(gnu_parse_datetime(MdsV0(s,i), &time_str))
      {
        mktime(&time_str);
        Mdi0(w,i,RLAB_DATE_IDX_SECS) = time_str.tm_sec;
        Mdi0(w,i,RLAB_DATE_IDX_MINS) = time_str.tm_min;
        Mdi0(w,i,RLAB_DATE_IDX_HOUR) = time_str.tm_hour;
        Mdi0(w,i,RLAB_DATE_IDX_DAYS) = time_str.tm_mday;
        Mdi0(w,i,RLAB_DATE_IDX_MONT) = time_str.tm_mon + 1;
        Mdi0(w,i,RLAB_DATE_IDX_YEAR) = time_str.tm_year + 1900;
        Mdi0(w,i,RLAB_DATE_IDX_IDST) = time_str.tm_isdst;
        Mdi0(w,i,RLAB_DATE_IDX_GMTO) = time_str.tm_gmtoff/3600;
      }
    }
  }

end_of_ent_dstr2time:

  ent_Clean(e1);
  ent_Clean(e2);
  return ent_Assign_Rlab_MDR (w);
}

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#undef THIS_SOLVER
#define THIS_SOLVER "stat"
Ent *
ent_Stat (int nargs, Datum args[])
{
  Ent *e1=0;
  char *str=0;
  Btree *bw = btree_Create ();
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  struct stat fileLStat, fileStat;

  if (nargs != 1)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": access file information on file system.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf(rlab_stderr, THIS_SOLVER ":   t = " THIS_SOLVER "(fn)\n");
    fprintf(rlab_stderr, THIS_SOLVER ": where 'fn' is the filename.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    goto _stat_exit;
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    goto _stat_exit;
  }

  str = class_char_pointer(e1);
  if (isvalidstring(str)<1)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    goto _stat_exit;
  }

  if(lstat(str,&fileLStat) >= 0 && stat(str,&fileStat)>=0)
  {
    if (S_ISLNK(fileLStat.st_mode))
    {
      install (bw, "link", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISDIR(fileStat.st_mode))
    {
      install (bw, "dir", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISBLK(fileStat.st_mode))
    {
      install (bw, "block", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISREG(fileStat.st_mode))
    {
      install (bw, "file", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISFIFO(fileStat.st_mode))
    {
      install (bw, "fifo", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISSOCK(fileStat.st_mode))
    {
      install (bw, "socket", ent_Create_Rlab_Double(1.0));
    }
    if (S_ISCHR(fileStat.st_mode))
    {
      install (bw, "char", ent_Create_Rlab_Double(1.0));
    }
    install (bw, "size", ent_Create_Rlab_Double((double) fileStat.st_size));
    install (bw, "nlinks", ent_Create_Rlab_Double((double) fileStat.st_nlink));
    install (bw, "inode", ent_Create_Rlab_Double((double) fileStat.st_ino));
    install (bw, "mode", ent_Create_Rlab_Double((double) (fileStat.st_mode & (S_IRWXU | S_IRWXG | S_IRWXO))));
    install (bw, "atime", ent_Create_Rlab_Double((double) fileStat.st_atime));
    install (bw, "mtime", ent_Create_Rlab_Double((double) fileStat.st_mtime));
    install (bw, "ctime", ent_Create_Rlab_Double((double) fileStat.st_ctime));
    install (bw, "dev", ent_Create_Rlab_Double((double) fileStat.st_dev));
    install (bw, "uid", ent_Create_Rlab_Double((double) fileStat.st_uid));
    install (bw, "gid", ent_Create_Rlab_Double((double) fileStat.st_gid));
    if ((fileStat.st_mode & S_ISVTX)>0)
      install (bw, "sticky_bit", ent_Create_Rlab_Double((double) 1.0));
    if ((fileStat.st_mode & S_ISUID)>0)
      install (bw, "set_uid", ent_Create_Rlab_Double((double) 1.0));
    if ((fileStat.st_mode & S_ISGID)>0)
      install (bw, "set_gid", ent_Create_Rlab_Double((double) 1.0));
  }

_stat_exit:

  ent_Clean(e1);
  return ent_Assign_Rlab_BTREE(bw);
}


//
// getpid
//
#include <sys/types.h>
#include <unistd.h>
#undef THIS_SOLVER
#define THIS_SOLVER "getpid"
Ent *
ent_Getpid (int nargs, Datum args[])
{
  return ent_Create_Rlab_Double ((double) getpid());
}

//
// fnctl
//
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#undef THIS_SOLVER
#define THIS_SOLVER "fnctl_lock"
Ent *
ent_Fnctl (int nargs, Datum args[])
{
  Ent *e1=0;

  MDR *w=0;
  char *s=0;
  int fd;

  struct flock lock;
  memset (&lock, 0, sizeof(struct flock));

  if (nargs != 1)
  {
    fprintf(stdout, THIS_SOLVER ": access file information on file system.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   t = " THIS_SOLVER "(lockfn)\n");
    fprintf(stdout, THIS_SOLVER ": where 'lockfn' is the filename used for locking.\n");
    fprintf(stdout, THIS_SOLVER ": one argument required\n");
    goto fnctl_exit;
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_STRING)
  {
    s = class_char_pointer(e1);
    fd = open (s, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR);
    if (fd == -1)
      fprintf(stdout, THIS_SOLVER ": Internal error %s\n", strerror(errno));
    else
    {
      lock.l_type = F_WRLCK;
      fcntl (fd, F_SETLKW, &lock);
    }
    w = mdr_CreateScalar(fd);
  }
  else if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    fd = (int) class_double(e1);
    lock.l_type = F_UNLCK;
    fcntl (fd, F_SETLKW, &lock);
    close (fd);
    w = mdr_CreateScalar(-2);
  }

fnctl_exit:

  ent_Clean(e1);
  return ent_Assign_Rlab_MDR (w);
}

