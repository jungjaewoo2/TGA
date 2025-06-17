// file: libtime.r
// functions for manipulations with time and dates
// rlabplus (c) 2000-2013 marijan kostrun

// many functions based on matlab collection timeutil by
//   Author:      Peter John Acklam
//   Time-stamp:  2002-05-24 13:30:06 +0200
//   E-mail:      pjacklam@online.no
//   URL:         http://home.online.no/~pjacklam/matlab/software/util/timeutil

static(...
    RLAB_ERROR_FIRSTARG_1_6_8_COL, ...
    RLAB_ERROR_FIRSTARG_DATESTRING,...
        RLAB_ERROR_FIRSTARG_3_COL,...
            RLAB_ERROR_FIRSTARG_1_COL, ...
                RLAB_ERROR_SECONDARG_1_COL ...
      );
RLAB_ERROR_FIRSTARG_1_6_8_COL = "First argument must be 1-col, 6-col or 8-col real matrix!";
RLAB_ERROR_FIRSTARG_DATESTRING = "First argument must be date string!";
RLAB_ERROR_FIRSTARG_3_COL = "First argument must be 3-column real matrix!";
RLAB_ERROR_FIRSTARG_1_COL = "First argument must be 1-column real vector!";
RLAB_ERROR_SECONDARG_1_COL= "Second argument must be 1-column real vector!";

static(_days_in_month, _sumdays_in_month);
_days_in_month =    [31, 28, 31,  30,  31,  30,  31,  31,  30,  31,  30,  31];
_sumdays_in_month = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];

static(TIMESTAMP_DEFAULT_FORMAT,TIMESTAMP_FORMAT);
TIMESTAMP_DEFAULT_FORMAT = "%Y%m%d-%H%M%S";

timestamp = function( x )
{
  if (!exist(x))
  {
    if (!exist(TIMESTAMP_FORMAT))
    { TIMESTAMP_FORMAT = TIMESTAMP_DEFAULT_FORMAT; }
  } else {
    if (class(x)!="string")
    {
      TIMESTAMP_FORMAT = TIMESTAMP_DEFAULT_FORMAT;
    } else {
      TIMESTAMP_FORMAT = x[1];
    }
  }
  return time2dstr(seconds(),TIMESTAMP_FORMAT);
};

timezone = function ( tz )
{
  if (exist(tz))
  {
    if (class(tz)=="string")
    {
      if (strlen(tz))
      {
        putenv("TZ=" + toupper(tz));
      }
    }
  }

  tz = reads("| /bin/date '+%Z' 2>/dev/null");
  return tz;
};

//
//
//
datesinrange = function(t, r)
{
  THIS_FUNCTION="datesinrange";
  if (class(t)!="num"|| type(t)=="complex")
  {
    t = seconds();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = seconds(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }
  if (class(r)!="num"|| type(r)=="complex")
  {
    r = 0;
  } else {
    if (r.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_SECONDARG_1_COL);
      return [];
    }
  }

  if (any(r != 0))
  {
    t = t + 86400 .* r[:];
  }

  return clock(t);
};

dayofmonth = function( t )
{
  THIS_FUNCTION = "dayofmonth";

  if (!exist(t))
  {
    tr = clock();
  } else {
    if (t.nc!=1 && t.nc!=6 && t.nc!=8)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }

    tr = clock(t);
  }

  return (tr[;3]-1 + (tr[;4]+tr[;5]./60+tr[;6]./3600) ./ 24);
};

isleapyear = function( y )
{
  THIS_FUNCTION = "isleapyear";
  if (class(y)!="num"|| type(y)=="complex")
  {
    y = clock()[1];
  } else {
    if (y.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
      return [];
    }
  }

  rval = ( !mod(y, 4) && mod(y, 100) ) || !mod(y, 400);

  return rval;
};

monthofyear = function( t )
{
  THIS_FUNCTION = "monthofyear";

  if (!exist(t))
  {
    tr = clock();
  } else {

    if (t.nc!=1 && t.nc!=6 && t.nc!=8)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }

    tr = clock(t);
  }

  return tr[;2];
};

daysinmonth = function(m, y)
{
  THIS_FUNCTION = "daysinmonth";
  if ( (class(m)!="num") || (type(m)=="complex") )
  {
    m = clock()[2];
  } else {
    if (m.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
      return [];
    }
  }
  if ((class(y)!="num") || (type(y)=="complex"))
  {
    y = clock()[1];
  } else {
    if (y.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_SECONDARG_1_COL);
      return [];
    }
  }

  rval = _days_in_month[m];

  // Add leap day as necessary.
  rval = rval + ( (m == 2) && isleapyear(y) );

  return rval;
};

daysinyear = function ( y )
{
  THIS_FUNCTION = "daysinyear";
  if (class(y)!="num"|| type(y)=="complex")
  {
    y = clock()[1];
  } else {
    if (y.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
      return [];
    }
  }

  rval = 365 + isleapyear( y );
  return rval;
};

weeksinyear = function ( y )
{
  THIS_FUNCTION = "weeksinyear";
  if (class(y)!="num"|| type(y)=="complex")
  {
    y = clock()[1];
  } else {
    if (y.nc!=1)
    {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
      return [];
    }
  }

  rval = (365 + isleapyear(y) ...
    + dayofweek( [y,1,4,0,0,0] ) - dayofweek( [y+1,1,4,0,0,0]))/7;

  return rval;
};

date2jd = function( t )
{
  //   JD = DATE2JD( t ) returns the Julian
  //   day number of the given date (Gregorian calendar) plus a fractional part
  //   depending on the time of day.
  //
  //   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
  //   MINUTE or SECOND will be replaced by zeros.
  //
  //   If no date is specified, the current date and time is used.
  //
  //   Start of the JD (Julian day) count is from 0 at 12 noon 1 January -4712
  //   (4713 BCE), Julian proleptic calendar.  Note that this day count conforms
  //   with the astronomical convention starting the day at noon, in contrast
  //   with the civil practice where the day starts with midnight.
  //
  //   Astronomers have used the Julian period to assign a unique number to
  //   every day since 1 January 4713 BCE.  This is the so-called Julian Day
  //   (JD).  JD 0 designates the 24 hours from noon UTC on 1 January 4713 BCE
  //   (Julian proleptic calendar) to noon UTC on 2 January 4713 BCE.

  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-05-24 13:30:06 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  THIS_FUNCTION = "date2jd";

  if (class(t)!="num")
  {
    t = clock();
  }

  if (t.nc==1 || t.nc==6 || t.nc==8)
  {
    t = clock(t);
  }
  else
  {
    printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
    return [];
  }

  // The following algorithm is a modified version of the one found in the
  // Calendar FAQ.
  a = floor((14 - t[;2])/12);
  y = t[;1] + 4800 - a;
  m = t[;2] + 12 * a - 3;

  // For a date in the Gregorian calendar:
  rval = t[;3] + floor((153*m + 2)/5) ...
    + y*365 + floor(y/4) - floor(y/100) + floor(y/400) - 32045 ...
    + ( t[;6] + 60*t[;5] + 3600*(t[;4] - 12) ) ./ 86400;

  return rval;
};

weekofyear = function( t )
{
  //   The week number is an integer between 1 and 53, inclusive.
  //
  //   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
  //   MINUTE or SECOND will be replaced by zeros.  If no date is specified,
  //   the current date and time is used.
  //
  //   This function is ISO 8601 compliant:  The first week of a given year is
  //   the first week which has more days in the given year than in the
  //   previous year.
  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)

  // convert it to 6-col format
  THIS_FUNCTION = "weekofyear";
  if (class(t)!="num"|| type(t)=="complex")
  {
    t = clock();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = clock(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }
  y = t[;1];

  // ISO 8601 states that a week that lies partly in one year and partly
  // in another is assigned a number in the year in which most of its
  // days lie.  Consequence: Week 1 of any year is the week that
  // contains 4 January.  Hence, the first week started on the following
  // day of the year
  //
  //    yd1 =  4: week 1 started on 4 January of current year
  //    yd1 =  3: week 1 started on 3 January of current year
  //    yd1 =  2: week 1 started on 2 January of current year
  //    yd1 =  1: week 1 started on 1 January of current year
  //    yd1 =  0: week 1 started on 31 December of previous year
  //    yd1 = -1: week 1 started on 30 December of previous year
  //    yd1 = -2: week 1 started on 29 December of previous year
  yd1 = 5 - floor(dayofweek( [y,1,4,0,0,0] ));

  // Get the ordinal day number and calculate the "raw" week number.
  yd2 = floor(dayofyear( t ));
  w = 1 + floor( (yd2 - yd1 + 1)/7 );

  // Take care of the case when the week number is 0.
  i = find(w == 0);
  if (!isempty(i))
  {
    y[i] = y[i] - 1;
    w[i] = weeksinyear(y[i]);
  }

  // Take care of the case when the week number is 53.
  i = find(w == 53);
  if (!isempty(i))
  {
    j = find( weeksinyear(y[i]) == 52);
    if(any(j[:]))
    {
      y[i[j]] = y[i[j]] + 1;
      w[i[j]] = 1;
    }
  }

  w = w + (dayofweek(t) - 1) / 7;
  return w;
};

date2mjd = function ( t )
{
  THIS_FUNCTION = "date2mjd";
  if (class(t)!="num"|| type(t)=="complex")
  {
    t = clock();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = clock(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }

  rval = date2jd( t ) - 2400000.5;
  return rval;
};

yearnum = function ( t )
{
  THIS_FUNCTION = "yearnum";
  if (class(t)!="num"|| type(t)=="complex")
  {
    t = clock();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = clock(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }

  year = t[;1] + (dayofyear(t) - 1) ...
    ./ daysinyear(t[;1]);

  return year;
};

//
// find the day 'd' in the week determined by time 't' and week 'k' from it
//  if k==0 and 't' is not given find the day 'd' in the current week
//  if k<0 then find the day 'd' in the k-th week before this one
//    and opposite for k>0
//  if 't' is provided then use that time as the point from which the day
//  'd' in the same week is found, and then go 'k' weeks earlier (k<0) or later (k>0)
otherday = function(t, d, k)
{
  THIS_FUNCTION = "otherday";

  if (class(t)!="num"|| type(t)=="complex")
  {
    t = seconds();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = seconds(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }
  if (class(d)!="num" || type(d)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_SECONDARG_1_COL);
    return [];
  }
  if (d.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_SECONDARG_1_COL);
    return [];
  }

  // k=0:  is how many weeks ago, or in the future
  if(!exist(k))
  { k = 0; }

  // d is the number of the day
  d = mod(abs(d),7);
  d = ifelse(d>0,d,7);

  // find local day
  ts = seconds(t);
  ld = dayofweek(ts);

  // find day 'd' from
  //dd = ld - d + 7 .* (ld < d);
  dd = ld - d;
  ts = ts - 24 * 3600 * (dd - k*7);

  return clock(ts);
};

earlierday = function(d, t)
{
  THIS_FUNCTION = "earlierday";

  // find first appearance of the day 'd' before time 't'
  if (class(d)!="num" || type(d)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (d.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  if (!exist(t))
  {
    t = clock();
  } else {
    t = clock(t);
  }

  // d is the number of the day
  d = mod(abs(d),7);
  d = ifelse(d,d,7);

  // find local day
  ts = seconds(t);
  ld = dayofweek(t);

  // find day 'd' from
  dd = ld - d + 7 .* (ld < d);
  dd = dd + (!dd)*7;
  ts = ts - 24 * 3600 * dd;

  return clock(ts);
};

laterday = function(d, t)
{
  THIS_FUNCTION = "laterday";

  // find first appearance of the day 'd' after time 't'
  if (class(d)!="num" || type(d)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (d.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  if (!exist(t))
  {
    t = clock();
  } else {
    t = clock(t);
  }

  // d is the number of the day
  d = mod(abs(d),7);
  d = ifelse(d,d,7);

  // find local day
  ts = seconds(t);
  ld = dayofweek(t);

  // find day 'd' from
  dd = d - ld + 7 .* (ld > d);
  ts = ts + 24 * 3600 * dd;
  return clock(ts);
};

days2hms = function(days)
{
  THIS_FUNCTION = "days2hms";

  //DAYS2HMS Convert days into hours, minutes, and seconds.
  //
  //   [HOUR, MINUTE, SECOND] = DAYS2HMS(DAYS) converts the number of days to
  //   hours, minutes, and seconds.
  //
  //   The following holds (to within rounding precision):
  //
  //     DAYS = HOUR / 24 + MINUTE / (24 * 60) + SECOND / (24 * 60 * 60)
  //          = (HOUR + (MINUTE + SECOND / 60) / 60) / 24
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:52:02 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(days)!="num" || type(days)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (days.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  s = 86400 * days;
  h = floor(s/3600);
  s = s - 3600 * h;
  m = floor(s/60);
  s = s - 60 * m;

  rval = [h, m, s];
  return rval;
};

hms2hrs = function ( hms )
{
  THIS_FUNCTION = "hms2hrs";

  //HMS2DAYS Convert hours, minutes, and seconds to hours.
  //
  //   DAYS = HMS2DAYS(HOUR, MINUTE, SECOND) converts the number of hours,
  //   minutes, and seconds to a number of days.
  //
  //   The following holds (to within rounding precision):
  //
  //     HRS = HOUR + MINUTE / (60) + SECOND / (60 * 60)
  //          = (HOUR + (MINUTE + SECOND / 60) / 60)
  //   Author:      Peter John Acklam
  //   Time-stamp:  2004-09-22 08:45:33 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(hms)!="num" || type(hms)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }
  if (hms.nc !=3)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }

  hrs = hms[;1] + (hms[;2] + hms[;3] ./ 60) ./ 60;
  return hrs;
};

hms2days = function ( hms )
{
  THIS_FUNCTION = "hms2days";

  //HMS2DAYS Convert hours, minutes, and seconds to days.
  //
  //   DAYS = HMS2DAYS(HOUR, MINUTE, SECOND) converts the number of hours,
  //   minutes, and seconds to a number of days.
  //
  //   The following holds (to within rounding precision):
  //
  //     DAYS = HOUR / 24 + MINUTE / (24 * 60) + SECOND / (24 * 60 * 60)
  //          = (HOUR + (MINUTE + SECOND / 60) / 60) / 24
  //   Author:      Peter John Acklam
  //   Time-stamp:  2004-09-22 08:45:33 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(hms)!="num" || type(hms)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }
  if (hms.nc !=3)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }

  days = (hms[;1] + (hms[;2] + hms[;3] ./ 60) ./ 60) ./ 24;
  return days;
};

hms2sec = function ( hms )
{
  THIS_FUNCTION = "hms2sec";

  //HMS2SEC  Convert from hours, minutes and seconds to seconds.
  //
  //   HMS2SEC(HOUR, MINUTE, SECOND) converts the number of hours, minutes and
  //   seconds to seconds.

  //   Author:      Peter John Acklam
  //   Time-stamp:  2003-06-22 19:22:49 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(hms)!="num" || type(hms)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }
  if (hms.nc !=3)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }

  rval = hms[;3] + 60 * hms[;2] + 3600 * hms[;1];
  return rval;
};

sec2hms = function ( s )
{
  THIS_FUNCTION = "sec2hms";

  //SEC2HMS  Convert seconds to hours, minutes and seconds.
  //
  //   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
  //   SEC into hours, minutes and seconds.

  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:50:09 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(s)!="num" || type(s)=="complex")
  {
    printf("%s: Warning: %s\n",THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (s.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  h = floor(s/3600);
  s = s - 3600 * h;
  m = floor(s/60);
  s = s - 60*m;

  return [h,m,s];
};

jd2date = function ( jd )
{
  THIS_FUNCTION = "jd2date";

  //JD2DATE Gregorian calendar date from modified Julian day number.
  //
  //   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2DATE(JD) returns the
  //   Gregorian calendar date (year, month, day, hour, minute, and second)
  //   corresponding to the Julian day number JD.
  //
  //   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
  //   (4713 BCE), Julian proleptic calendar.  Note that this day count conforms
  //   with the astronomical convention starting the day at noon, in contrast
  //   with the civil practice where the day starts with midnight.
  //
  //   Astronomers have used the Julian period to assign a unique number to
  //   every day since 1 January 4713 BCE.  This is the so-called Julian Day
  //   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BCE
  //   (Julian calendar) to noon UTC on 2 January 4713 BCE.
  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-05-24 15:24:45 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  // Adding 0.5 to JD and taking FLOOR ensures that the date is correct.
  // Here are some sample values:
  //
  //  MJD     Date       Time
  //  -1.00 = 1858-11-16 00:00 (not 1858-11-15 24:00!)
  //  -0.75 = 1858-11-16 06:00
  //  -0.50 = 1858-11-16 12:00
  //  -0.25 = 1858-11-16 18:00
  //   0.00 = 1858-11-17 00:00 (not 1858-11-16 24:00!)
  //  +0.25 = 1858-11-17 06:00
  //  +0.50 = 1858-11-17 12:00
  //  +0.75 = 1858-11-17 18:00
  //  +1.00 = 1858-11-18 00:00 (not 1858-11-17 24:00!)

  if (class(jd)!="num" || type(jd)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (jd.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  // integer part
  ijd = floor(jd + 0.5);
  if (jd!=ijd)
  {
    fjd = jd - ijd + 0.5;
    hms = days2hms(fjd);
  } else {
    hms = [0,0,0];
  }

  // The following algorithm is from the Calendar FAQ.

  a = ijd + 32044;
  b = floor((4 * a + 3) / 146097);
  c = a - floor((b * 146097) / 4);

  d = floor((4 * c + 3) / 1461);
  e = c - floor((1461 * d) / 4);
  m = floor((5 * e + 2) / 153);

  day = e - floor((153 * m + 2) / 5) + 1;
  mon = m + 3 - 12 * floor(m / 10);
  yea = b * 100 + d - 4800 + floor(m / 10);

  rval = [yea, mon, day, hms];

  return rval;
};

jd2jdate = function( jday )
{
  THIS_FUNCTION = "jd2jdate";

  //JD2JDATE Julian day number to Julian calendar date.
  //
  //   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2JDATE(JDAY) returns the
  //   Julian calendar date (year, month, day, hour, minute, and second)
  //   corresponding to the Julian day number JDAY.
  //
  //   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
  //   (4713 BCE), Julian proleptic calendar.  Note that this day count conforms
  //   with the astronomical convention starting the day at noon, in contrast
  //   with the civil practice where the day starts with midnight.
  //
  //   Astronomers have used the Julian period to assign a unique number to
  //   every day since 1 January 4713 BCE.  This is the so-called Julian Day
  //   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BCE
  //   (Julian calendar) to noon UTC on 2 January 4713 BCE.

  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)

  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-05-24 15:24:45 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(jday)!="num" || type(jday)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (jday.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  ijday = floor(jday);                 // integer part
  fjday = jday - ijday;                // fraction part

  second = 86400 * fjday;
  hour   = floor(second / 3600);         // get number of hours
  second = second - 3600 * hour;         // remove the hours
  minute = floor(second / 60);           // get number of minutes
  second = second - 60 * minute;         // remove the minutes
  hour   = hour + 12;                  // Julian days start at noon

  // The following algorithm is from the Calendar FAQ.
  b = 0;
  c = ijday + 32082;

  d = floor((4 * c + 3) / 1461);
  e = c - floor((1461 * d) / 4);
  m = floor((5 * e + 2) / 153);

  day   = e - floor((153 * m + 2) / 5) + 1;
  month = m + 3 - 12 * floor(m / 10);
  year  = b * 100 + d - 4800 + floor(m / 10);

  rval = [year, month, day, hour, minute, second];
  return rval;
};

jd2mjd = function (jd)
{
  //JD2MJD Modified Julian day number from Julian day number.
  //
  //   [MJD] = JD2MJD(JD) returns the modified Julian day number corresponding
  //   to the given Julian day number.
  //
  //   See also MJD2JD, DATE2MJD, DATE2JD.
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:51:32 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  rval = jd - 2400000.5;
  return rval;
};

jdate2jd = function( t )
{
  THIS_FUNCTION = "jdate2jd";

  //JDATE2JD Julian day number from Julian date.
  //
  //   JDATE2JD(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the Julian day
  //   number of the given date (Julian calendar) plus a fractional part
  //   depending on the time of day.
  //
  //   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
  //   MINUTE or SECOND will be replaced by zeros.
  //
  //   If no date is specified, the current date and time is used.
  //
  //   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
  //   (4713 BCE), Julian proleptic calendar.  Note that this day count conforms
  //   with the astronomical convention starting the day at noon, in contrast
  //   with the civil practice where the day starts with midnight.
  //
  //   Astronomers have used the Julian period to assign a unique number to
  //   every day since 1 January 4713 BCE.  This is the so-called Julian Day
  //   (JD).  JD 0 designates the 24 hours from noon UTC on 1 January 4713 BCE
  //   (Julian calendar) to noon UTC on 2 January 4713 BCE.
  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-05-24 13:31:03 +0200
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(t)!="num" || type(t)=="complex")
  {
    t = clock();
  } else {
    if (t.nc==1 || t.nc==6 || t.nc==8)
    {
      t = clock(t);
    } else {
      printf("%s: %s\n",THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_6_8_COL);
      return [];
    }
  }

  // The following algorithm is from the Calendar FAQ.  The one in the
  // Calendar FAQ is correct back to, and including, -4800-03-01 Julian
  // proleptic calendar.  The algorithm below is correct for all dates in
  // the Julian proleptic calendar.
  a = floor((14 - t[;2])/12);
  y = t[;1] + 4800 - a;
  m = t[;2] + 12*a - 3;

  jd = t[;3] + floor((153*m + 2)/5) + y*365 + floor(y/4) - 32083 ...
    + ( t[;6] + 60*t[;5] + 3600*(t[;4] - 12) )/86400;

  return jd;
};

mjd2jd = function (mjd)
{
  THIS_FUNCTION = "mjd2jd";

  //MJD2JD Modified Julian day number from Julian day number.
  //
  //   JD = MJD2JD(MJD) returns the Julian day number corresponding to the
  //   given modified Julian day number.
  //
  //   See also JD2MJD, DATE2JD, DATE2MJD.
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:50:25 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(mjd)!="num" || type(mjd)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (mjd.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  rval = mjd + 2400000.5;
  return rval;
};

mjd2date = function( mjd )
{
  THIS_FUNCTION="mjd2date";

  //MJD2DATE Gregorian calendar date from Julian day number.
  //
  //   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = MJD2DATE(MJD) returns the
  //   Gregorian calendar date (year, month, day, hour, minute, and second)
  //   corresponding to the Julian day number JDAY.
  //
  //   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
  //   (4713 BCE), Julian proleptic calendar.  Note that this day count conforms
  //   with the astronomical convention starting the day at noon, in contrast
  //   with the civil practice where the day starts with midnight.
  //
  //   Astronomers have used the Julian period to assign a unique number to
  //   every day since 1 January 4713 BCE.  This is the so-called Julian Day
  //   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BCE
  //   (Julian calendar) to noon UTC on 2 January 4713 BCE.

  //   Sources:  - http://tycho.usno.navy.mil/mjd.html
  //             - The Calendar FAQ (http://www.faqs.org)

  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:50:30 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  // We could have got everything by just using
  //
  //   jd = mjd2jd(mjd);
  //   [year, month, day, hour, minute, second] = jd2date(jd);
  //
  // but we lose precision in the fraction part when MJD is converted to JD
  // because of the large offset (2400000.5) between JD and MJD.

  if (class(mjd)!="num" || type(mjd)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION, RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }
  if (mjd.nc !=1)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_1_COL);
    return [];
  }

  jd = mjd2jd(mjd);
  ymd = jd2date(jd);

  fmjd = mjd - floor(mjd);
  hms  = days2hms(fmjd);

  rval = [ymd, hms];
  return rval;
};

isdate = function( ymd )
{
  THIS_FUNCTION="isdate";

  //ISDATE True for valid dates (Gregorian calendar).
  //
  //   ISDATE(YEAR, MONTH, DAY) returns 1 if input is a valid year-month-date
  //   triple and 0 otherwise.  Gregorian calendar is assumed.
  //
  //   See also ISJDATE.
  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:51:49 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam
  //
  //   I fixed bug here so putting something like
  //    isdate([15,15,15])
  //  returns 0 rather then causing an error

  if (class(ymd)!="num" || type(ymd)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }
  if (ymd.nc !=3)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }

  t =   (ymd[;1] == round(ymd[;1])) && (ymd[;2] == round(ymd[;2])) ...
     && (ymd[;3] == round(ymd[;3])) && (1 <= ymd[;2]) && (ymd[;2] <= 12) && (1 <= ymd[;3]);

  // proceed only for valid dates so far
  _j = find(t);
  if(isempty(_j))
  { return t; }

  is_february = (ymd[_j;2] == 2);
  is_leapyear = isleapyear(ymd[_j;1]);
  days_in_month = _days_in_month[ ymd[_j;2] ] + ( is_february && is_leapyear );

  t[_j] = t[_j] && ( ymd[_j;3] <= days_in_month );
  return t;
};

isjdate = function( ymd )
{
  THIS_FUNCTION="isjdate";

  //ISJDATE True for valid dates (Julian calendar).
  //
  //   ISDATE(YEAR, MONTH, DAY) returns 1 if input is a valid year-month-date
  //   triple and 0 otherwise.  Julian calendar is assumed.
  //
  //   See also ISDATE.

  //   Author:      Peter John Acklam
  //   Time-stamp:  2002-03-03 12:51:47 +0100
  //   E-mail:      pjacklam@online.no
  //   URL:         http://home.online.no/~pjacklam

  if (class(ymd)!="num" || type(ymd)=="complex")
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }
  if (ymd.nc !=3)
  {
    printf("%s: Warning: %s\n", THIS_FUNCTION,RLAB_ERROR_FIRSTARG_3_COL);
    return [];
  }

  t =   (ymd[;1] == round(ymd[;1])) && (ymd[;2] == round(ymd[;2])) ...
     && (ymd[;3] == round(ymd[;3])) && (1 <= ymd[;2]) && (ymd[;2] <= 12) && (1 <= ymd[;3]);

  // proceed only for valid dates so far
  _j = find(t);
  if(isempty(_j))
  { return t; }

  is_february = (ymd[_j;2] == 2);
  is_leapyear = !mod(ymd[;1], 4);
  days_in_month = _days_in_month[ ymd[_j;2] ] + ( is_february && is_leapyear );

  t[_j] = t[_j] && ( ymd[_j;3] <= days_in_month );
  return t;
};


