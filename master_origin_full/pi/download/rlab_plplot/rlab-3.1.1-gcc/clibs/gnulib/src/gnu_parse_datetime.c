//
//
//
#include <config.h>
#include <time.h>
#include <timespec.h>
#include <stdbool.h>

bool
parse_datetime (struct timespec *result, char const *p,
                struct timespec const *now);

int gnu_parse_datetime(char *s, struct tm * ts)
{
  struct timespec result;
  struct tm * time_str;
  time_t tt1;

  if (parse_datetime(&result, s, (struct timespec const *)NULL))
  {
    tt1 = (time_t) result.tv_sec;
    time_str = localtime ( &tt1 );

    ts->tm_sec = time_str->tm_sec;
    ts->tm_min = time_str->tm_min;
    ts->tm_hour = time_str->tm_hour;
    ts->tm_mday = time_str->tm_mday;
    ts->tm_mon = time_str->tm_mon;
    ts->tm_year = time_str->tm_year;
    ts->tm_wday = time_str->tm_wday;
    ts->tm_yday = time_str->tm_yday;
    ts->tm_isdst = time_str->tm_isdst;
    ts->tm_gmtoff = time_str->tm_gmtoff;

    return 1;
  }

  return 0;
}
