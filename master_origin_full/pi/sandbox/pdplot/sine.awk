awk '
   BEGIN {
     print "sine wave"
     print "xscale 1 value"
     print "yscale 1 time"
     for (i=0; i<=3; i+=.02) {
        print i, sin(2*3.14159*i)
     }
   }
'
