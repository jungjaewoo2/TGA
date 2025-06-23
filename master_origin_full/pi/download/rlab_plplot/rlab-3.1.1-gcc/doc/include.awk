# include - replace #include "F" by contents of file F
#           replacie #includecode "F" by contents of file F, with line-numbers,
#           and certain edits for SGML idiosyncracies.

/^#include / {
    gsub(/"/, "", $2)
    while (getline x <$2 > 0)
        print x
    next
}
/^#includecode / {
    gsub(/"/, "", $2);
    lineno = 1;
    while (getline x <$2 > 0)
    {
      gsub("&","\&ero;",x);
      printf("%3i: %s\n", lineno++, x);
    }
    next
}
{ print }
