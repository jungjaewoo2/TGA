#
# A simple AWK script to change some names...
#

{ gsub("pclose", "plclose"); 
  gsub("pend", "plend"); 
  gsub("pstart", "plstart"); 
  gsub("ptitle", "pltitle"); 
  gsub("pwin", "plwin"); 
  gsub("showpwin", "showplwin"); 
  print; }
